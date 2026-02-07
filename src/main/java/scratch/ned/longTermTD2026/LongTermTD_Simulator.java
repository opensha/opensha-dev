package scratch.ned.longTermTD2026;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;

import com.google.common.io.Files;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.calc.recurInterval.EqkProbDistCalc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.AperiodicityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.AperiodicityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.HistoricalOpenInterval;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.RenewalModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepUtils;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.UCERF3_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.WG02_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;
import org.opensha.sha.faultSurface.RuptureSurface;

import scratch.UCERF3.erf.utils.ProbModelsPlottingUtils;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc.RenewalModelType;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;

public class LongTermTD_Simulator {
	
	final static double MILLISEC_PER_YEAR = TimeDepUtils.MILLISEC_PER_YEAR;

	
	public static void simulateEvents(TimeDepFaultSystemSolutionERF erf, String inputTimeSinceLastMillisFileName, String outputTimesinceLastMillisFileName, 
			double numYears, File resultsDir, long randomSeed, boolean verbose, boolean makePlots, double timeStepYrs) {

		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		UCERF3_ProbabilityModel u3ProbModel = null;
		WG02_ProbabilityModel wgProbModel = null;
		boolean isPoisson = false;
		
		boolean aveRecurIntervals = true; 		// default value applied by WG02
		boolean aveNormTimeSinceLast = true;
		AperiodicityModel aperModel = null;
		if (probModel instanceof UCERF3_ProbabilityModel) {			// U3 calculation type
			u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			aveRecurIntervals = u3ProbModel.getAveragingTypeChoice().isAveRI();
			aveNormTimeSinceLast = u3ProbModel.getAveragingTypeChoice().isAveNTS();
			aperModel = u3ProbModel.getAperiodicityModel();
		} else if (probModel instanceof WG02_ProbabilityModel) {
			wgProbModel = (WG02_ProbabilityModel)probModel;
			aperModel = wgProbModel.getAperiodicityModel();
		}
		else if (probModel instanceof FSS_ProbabilityModel.Poisson) {
			isPoisson = true;
		}
		else
			throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());
		
//		// can also do things via the parameter list
//		ParameterList probModelParams = probModel.getAdjustableParameters();
//		if (probModelParams.containsParameter(RenewalModels.PARAM_NAME))
//			probModelParams.setValue(RenewalModels.PARAM_NAME, RenewalModels.BPT);		

		FaultSystemSolution fltSysSolution=erf.getSolution();
		// this has zeros where events were filtered our by the ERF (mags too low); 
		double[] longTermRateOfFltSysRup = erf.getLongTermRateOfFltSysRupInERF(); // this can be different than what's in the FSS
		
		FaultSystemRupSet fltSysRupSet = fltSysSolution.getRupSet();
		int numRupsInFaultSystem = fltSysRupSet.getNumRuptures();
		int numSections = fltSysRupSet.getNumSections();

		
		double[] longTermPartRateForSectArray = probModel.getSectLongTermPartRates(); // this is a duplicate
		double[] sectionArea = fltSysRupSet.getAreaForAllSections();
		long[] dateOfLastForSect = probModel.getSectDOLE(); // this is a duplicate
		
		
		if(resultsDir != null)
			if(!resultsDir.exists()) 
				resultsDir.mkdir();
		
		File plotsDir = null;
		File otherPlotsDir = null;
		if(makePlots) {
			if(resultsDir == null)
				throw new RuntimeException ("Cannot save plots because resultsDir == null");
			plotsDir = new File(resultsDir, "plots");
			if(!plotsDir.exists()) plotsDir.mkdir();
			
			otherPlotsDir = new File(plotsDir, "otherPlots");
			if(!otherPlotsDir.exists()) otherPlotsDir.mkdir();
		}
		

		RandomDataGenerator randomDataSampler = new RandomDataGenerator();
		randomDataSampler.reSeed(randomSeed);  // for reproducibility; this is for the next event time
		Random random = new Random(randomSeed);  // this is for the nth rup sampler
		
//		ProbabilityModelOptions probType = (ProbabilityModelOptions)erf.getParameter(ProbabilityModelParam.NAME).getValue();

		String tempString1 = "\t"+erf.getAdjustableParameterList().getParameterListMetadataString();
		String tempString2 = "\n\t"+erf.getAdjustableParameterList().getParameter(
				TimeDepFaultSystemSolutionERF.PROB_MODEL_PARAM_NAME).getValue().toString();
		String erfParamMetadataString = tempString1.replace(";", "\n\t")+tempString2.replace(";", "\n\t");
		String infoString="";
		
		if(resultsDir != null) {
			infoString = "Information for simulation results in this directory ("+resultsDir+")\n\n"+
					"inputTimeSinceLastMillisFileName = "+inputTimeSinceLastMillisFileName+"\n\n"+
					"outputTimesinceLastMillisFileName = "+outputTimesinceLastMillisFileName+"\n\n"+
					"simulation duration = "+numYears+" (years)\n\n"+
					"randomSeed = "+randomSeed+"\n\n"+
					"ERF Parameters:\n\n";
			infoString += erfParamMetadataString;
			
//			System.out.println("Full ERF metadata string:\t"+erf.getAdjustableParameterList().getParameterListMetadataString());
//			System.out.println("Prob model metadata string:\t"+erf.getAdjustableParameterList().getParameter(
//					TimeDepFaultSystemSolutionERF.PROB_MODEL_PARAM_NAME).getValue().toString());

			
//			for(Parameter param:erf.getAdjustableParameterList()) {
//				String valueString = "null";
//				if(param.getValue() != null) valueString=param.getValue().toString();
//				infoString += "\t"+param.getName()+" = "+valueString+"\n";//+param.getValue().toString()+"\n";
//			}
//			infoString += "\tTimespan duration (ignored) = "+erf.getTimeSpan().getDuration()+"\n";
//			if(probType != ProbabilityModelOptions.POISSON)
//				infoString += "\tTimespan start year = "+erf.getTimeSpan().getStartTimeYear()+"\n\n";
//			else
//				infoString += "\tTimespan start year = 1970 (default for Poisson simulations)\n\n";			
		}
		
		if(verbose) System.out.println(infoString);
				
		// LABELING AND FILENAME STUFF
//		String typeCalcForU3_Probs;
//		if(aveRecurIntervals)
//			typeCalcForU3_Probs = "aveRI";
//		else
//			typeCalcForU3_Probs = "aveRate";
//		if(aveNormTimeSinceLast)
//			typeCalcForU3_Probs += "_aveNTS";
//		else
//			typeCalcForU3_Probs += "_aveTS";
//		
//		String probTypeString;
//		if (probType == ProbabilityModelOptions.POISSON) {
//			probTypeString= "Pois";
//		}
//		else if(probType == ProbabilityModelOptions.U3_BPT) {
//			probTypeString= "U3BPT";
//		}
//		else if(probType == ProbabilityModelOptions.WG02_BPT) {
//			probTypeString= "WG02BPT";
//		}
//		else
//			throw new RuntimeException("Probability type unrecognized");
//
//		String aperString = "aper";
//		if(probType != ProbabilityModelOptions.POISSON) {
//			boolean first = true;
//			for(double aperVal:aperValues) {
//				if(first)
//					aperString+=aperVal;
//				else
//					aperString+=","+aperVal;
//				first=false;
//			}
//			aperString = aperString.replace(".", "pt");			
//		}
//		
//		if(verbose) System.out.println("\naperString: "+aperString+"\n");
		
		int durationKyrs = (int) Math.round(numYears/1000);
		
//		String fullProbTypeString = probTypeString;
//		if(probType == ProbabilityModelOptions.U3_BPT)
//			fullProbTypeString += " ("+aperString+", "+typeCalcForU3_Probs+")";
//		else if(probType == ProbabilityModelOptions.WG02_BPT)
//			fullProbTypeString += " ("+aperString+")";

		// INTIALIZE THINGS:
		double[] srcGainsAtMaxRateTimeArray=null, srcGainsAtMaxAveGainTimeArray=null;
		double maxTotalRate=-1;
		double maxAveGainAtTime=-100;
		double maxTotalRateYr=-1, maxAveGainAtTimeYr=1;

		// set original start time and total duration
		long origStartTimeMillis = 0;
		if(!isPoisson)
			origStartTimeMillis = erf.getTimeSpan().getStartTimeInMillis();
		double origStartYear = ((double)origStartTimeMillis)/MILLISEC_PER_YEAR+1970.0;
		if(verbose) {
			System.out.println("orig start time: "+origStartTimeMillis+ " millis ("+origStartYear+" yrs)");
			System.out.println("numYears: "+numYears);
		}
				
		TreeMap<Long, Integer> nthRupAtEpochMap = new TreeMap<Long, Integer>();
		ArrayList<Double> totExpRateAtEventTimeList = new ArrayList<Double>();
		ArrayList<Double> normRI_ForEventList = new ArrayList<Double>();
//		ArrayList<Double> mag_ForEventList = new ArrayList<Double>();
		ArrayList<Integer> fltSysRupIndexForEventList = new ArrayList<Integer>();

		
    	ArbDiscrEmpiricalDistFunc_3D normRI_AlongStrike = new ArbDiscrEmpiricalDistFunc_3D(0.05d,0.95d,10);
		double[] obsSectRateArray = new double[numSections];
		double[] obsSectSlipRateArray = new double[numSections];
		double[] obsSectRateArrayMlt7pt3 = new double[numSections];
		double[] obsSectRateArrayMgt7pt3 = new double[numSections];

		double[] obsRupRateArray = new double[erf.getTotNumRups()];
		double[] aveRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time
		double[] minRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time
		double[] maxRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time

		// initialize tracking of sections that had one or more ruptures in simulation
		boolean[] sectionRupturedDuringSim = new boolean[numSections];  
		for(int i=0;i<numSections;i++)
			sectionRupturedDuringSim[i] = false;
		
		// this is for writing out simulated events that occur
		FileWriter eventFileWriter=null;
		if(resultsDir != null) {
			eventFileWriter=null;
			try {
				eventFileWriter = new FileWriter(resultsDir+"/sampledEventsData.txt");
				eventFileWriter.write("nthRupIndex\tfssRupIndex\tyear\tepoch\tnormRupRI\trupMag\trupArea\n");
			} catch (IOException e1) {
				e1.printStackTrace();
			}			
		}

		// temporarily set the forecast as Poisson, to get long-term rates, & no background 
		IncludeBackgroundOption includeBackground = (IncludeBackgroundOption)erf.getParameter(IncludeBackgroundParam.NAME).getValue();
//		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
		erf.updateForecast();

		// fill in totalRate, longTermRateOfNthRups, magOfNthRups, and longTermSlipRateForSectArray
		double totalRate=0;
		IntegerPDF_FunctionSampler nthRupRandomSampler = new IntegerPDF_FunctionSampler(erf.getTotNumRups());
		double[] longTermRateOfNthRups = new double[erf.getTotNumRups()];	// this will include any aftershock reductions
		if(obsRupRateArray.length != longTermRateOfNthRups.length)
			throw new RuntimeException("obsRupRateArray.length="+ obsRupRateArray.length+"\nlongTermRateOfNthRups.length="+longTermRateOfNthRups.length);
		double[] magOfNthRups = new double[erf.getTotNumRups()];
		double[] longTermSlipRateForSectArray = new double[numSections];
		ArrayList<Double> aperValuesList = new ArrayList<Double>();
		for(int nthRup=0; nthRup<erf.getTotNumRups(); nthRup++) {
			ProbEqkRupture rup = erf.getNthRupture(nthRup);
			double rate = rup.getMeanAnnualRate(erf.getTimeSpan().getDuration());
			longTermRateOfNthRups[nthRup] = rate;
			totalRate += rate;
			magOfNthRups[nthRup] = rup.getMag();
			nthRupRandomSampler.set(nthRup, rate);
			int fltSysIndex = erf.getFltSysRupIndexForNthRup(nthRup);
			// aperiodicities
			if(!isPoisson) {
				double aper =aperModel.getRuptureAperiodicity(fltSysIndex);
				if(!aperValuesList.contains(aper))
					aperValuesList.add(aper);
			}
			// slip rates
			List<Integer> sectIndices = fltSysRupSet.getSectionsIndicesForRup(fltSysIndex);
			double mag = fltSysRupSet.getMagForRup(erf.getFltSysRupIndexForNthRup(nthRup));
			double slips[] = getAveSlipOnSectionsForRup(erf, nthRup, mag, sectIndices.size());
			for(int s=0;s<sectIndices.size();s++) {
				int sectID = sectIndices.get(s);
				longTermSlipRateForSectArray[sectID] += rate*slips[s];
			}					
		}
		int numAperValues = aperValuesList.size();

		if(verbose) System.out.println("totalRate long term = "+totalRate);
		if(resultsDir != null) infoString += "\n\nTotal long-term rate (per year) = "+totalRate+"\n\n";
		
		// this is for storing section normalized RIs
		ArrayList<Double> normalizedSectRecurIntervals = new ArrayList<Double>();

//		// This does not make sense
//		ArrayList<ArrayList<Double>> normalizedSectRecurIntervalsMagDepList = new ArrayList<ArrayList<Double>>();
//		for(int i=0;i<aperValuesList.size();i++) {
//			normalizedSectRecurIntervalsMagDepList.add(new ArrayList<Double>());
//		}
		

		
		double totalLongTermRate = totalRate;
		double simDuration = 1/totalLongTermRate;  // used to compute next prob gain
		
		// make the target MFD - 
		SummedMagFreqDist targetMFD=null;
		double origTotMoRate=Double.NaN;
		targetMFD = ERF_Calculator.getTotalMFD_ForERF(erf, 5.05, 8.95, 40, true);
		origTotMoRate = ERF_Calculator.getTotalMomentRateInRegion(erf, null);
		System.out.println("originalTotalMomentRate: "+origTotMoRate);
		targetMFD.setName("Target MFD");
		String tempString = "total rate = "+(float)targetMFD.getTotalIncrRate();
		tempString += "\ntotal rate >= 6.7 = "+(float)targetMFD.getCumRate(6.75);
		tempString += "\ntotal MoRate = "+(float)origTotMoRate;
		targetMFD.setInfo(tempString);			

		// reset ERF to original state
		erf.setCustomProbabilityModel(probModel);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(includeBackground);
		erf.updateForecast();

		// MFD & MoRate for simulation
		SummedMagFreqDist simMFD = new SummedMagFreqDist(5.05,8.95,40);
		double obsMoRate = 0;
		
		// Make local sectIndexArrayForSrcList for faster simulations
		ArrayList<int[]> sectIndexArrayForSrcList = new ArrayList<int[]>();
		for(int s=0; s<erf.getNumFaultSystemSources();s++) {
			List<Integer> indexList = fltSysRupSet.getSectionsIndicesForRup(erf.getFltSysRupIndexForSource(s));
			int[] indexArray = new int[indexList.size()];
			for(int i=0;i<indexList.size();i++)
				indexArray[i] = indexList.get(i);
			sectIndexArrayForSrcList.add(indexArray);
		}
		
		// set the ave cond recurrence intervals
		double[] aveCondRecurIntervalForFltSysRups = TimeDepUtils.computeAveCondRecurIntervalForFltSysRups(fltSysRupSet, 
														longTermPartRateForSectArray, aveRecurIntervals);
		// print minimum and maximum conditional rate of rupture
		if(verbose || resultsDir != null) {
			double minCondRI=Double.MAX_VALUE,maxCondRI=0;
			for(double ri: aveCondRecurIntervalForFltSysRups) {
				if(!Double.isInfinite(ri)) {
					if(ri < minCondRI) minCondRI = ri;
					if(ri > maxCondRI) maxCondRI = ri;
				}
			}
			if(verbose) System.out.println("minCondRI="+minCondRI+"\nmaxCondRI="+maxCondRI);
			if(resultsDir != null) infoString +="minCondRI="+minCondRI+"\nmaxCondRI="+maxCondRI+"\n\n";
		}
		
		// set simulation time
		double currentYear=origStartYear;
		long currentTimeMillis = origStartTimeMillis;

		// read section date of last file if not null
		if(inputTimeSinceLastMillisFileName != null && !isPoisson) {
			readSectTimeSinceLastEventFromFile(inputTimeSinceLastMillisFileName, currentTimeMillis, dateOfLastForSect);
			for(int s=0; s<dateOfLastForSect.length;s++)
				probModel.setSectDOLE(s, dateOfLastForSect[s]);
		}

		// remove date of last event if Poisson
		if(isPoisson)
			for(int s=0; s<dateOfLastForSect.length;s++)
				dateOfLastForSect[s]= Long.MIN_VALUE;

		// this is to track progress
		int percDoneThresh=0;
		int percDoneIncrement=5;

		long startRunTime = System.currentTimeMillis();	
		if(verbose) System.out.println("Starting simulation loop");
		
		int numSectThatRuptured = 0;
		int numSimulatedRups=0;
		
		while (currentYear<numYears+origStartYear) {
			
			// write progress
			if(verbose) {
				int percDone = (int)Math.round(100*(currentYear-origStartYear)/numYears);
				if(percDone >= percDoneThresh) {
					double timeInMin = ((double)(System.currentTimeMillis()-startRunTime)/(1000.0*60.0));
					int numGoodDateOfLast=0;
					for(long dateOfLast:dateOfLastForSect) {
						if(dateOfLast != Long.MIN_VALUE)
							numGoodDateOfLast+=1;					
					}
					int percentGood = (int)Math.round((100.0*(double)numGoodDateOfLast/(double)dateOfLastForSect.length));
					System.out.println("\n"+percDoneThresh+"% done in "+(float)timeInMin+" minutes"+";  totalRate="+(float)totalRate+"; yr="+(float)currentYear+";  % sect with date of last = "+percentGood);	
					System.out.println("\tFraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+"\n");
			
					percDoneThresh += percDoneIncrement;
				}				
			}
			
			// update gains and sampler if not Poisson
			if(!isPoisson) {
				double[] probGainForFaultSystemSource = new double[erf.getNumFaultSystemSources()];
				double aveGainAtTime = 0;
				// first the gains
				for(int s=0;s<erf.getNumFaultSystemSources();s++) {
					int fltSysRupIndex = erf.getFltSysRupIndexForSource(s);
					probGainForFaultSystemSource[s] = probModel.getProbabilityGain(fltSysRupIndex, currentTimeMillis, simDuration);
				}

//				if(probType == ProbabilityModelOptions.U3_BPT) {
//					for(int s=0;s<erf.getNumFaultSystemSources();s++) {
//						int fltSysRupIndex = erf.getFltSysRupIndexForSource(s);
//						probGainForFaultSystemSource[s] = getU3_ProbGainForRup(fltSysRupIndex, erf.getAperiodicityForSource(s), 0.0, false, aveRecurIntervals, aveNormTimeSinceLast, currentTimeMillis, simDuration);
//					}
//				}
//				else if(probType == ProbabilityModelOptions.WG02_BPT) {
//					sectionGainArray=null; // set this null so it gets updated
//					for(int s=0;s<erf.getNumFaultSystemSources();s++) {
//						int fltSysRupIndex = erf.getFltSysRupIndexForSource(s);
//						probGainForFaultSystemSource[s] = getWG02_ProbGainForRup(fltSysRupIndex, erf.getAperiodicityForSource(s), false, currentTimeMillis, simDuration);
//					}
//				}		
				// now update totalRate and ruptureSampler (for all rups since start time changed)
				for(int n=0; n<erf.getTotNumRupsFromFaultSystem();n++) {
					double probGain = probGainForFaultSystemSource[erf.getSrcIndexForNthRup(n)];
					aveGainAtTime += probGain;
					double newRate = longTermRateOfNthRups[n] * probGain;	// applied as a rate gain
					nthRupRandomSampler.set(n, newRate);
					aveRupProbGainArray[n] += probGain;
					if(minRupProbGainArray[n]>probGain)
						minRupProbGainArray[n] = probGain;
					if(maxRupProbGainArray[n]<probGain)
						maxRupProbGainArray[n] = probGain;
				}
				totalRate = nthRupRandomSampler.getSumOfY_vals();
				aveGainAtTime /= erf.getNumFaultSystemSources();
				
				if(maxTotalRate<totalRate) {
					srcGainsAtMaxRateTimeArray = probGainForFaultSystemSource;
					maxTotalRate=totalRate;
					maxTotalRateYr = currentYear;
				}
				if(maxAveGainAtTime<aveGainAtTime) {
					srcGainsAtMaxAveGainTimeArray = probGainForFaultSystemSource;
					maxAveGainAtTime=aveGainAtTime;
					maxAveGainAtTimeYr = currentYear;

				}

			}
			
			
			// Try sampling an event time
			double timeToNextInYrs;
			long eventTimeMillis;
			if(Double.isNaN(timeStepYrs)) {  // no time step; sample time of next event
				timeToNextInYrs = randomDataSampler.nextExponential(1.0/totalRate);
				eventTimeMillis = currentTimeMillis + (long)(timeToNextInYrs*MILLISEC_PER_YEAR);
			}
			else {
				long numEvents = randomDataSampler.nextPoisson(totalRate*timeStepYrs);
				if(numEvents > 0) { // get a random number or events for next time step;  set timeStepYrs low enough that more than 1 event is rare
					timeToNextInYrs = timeStepYrs*randomDataSampler.nextUniform(0, 1);
					eventTimeMillis = currentTimeMillis + (long)(timeToNextInYrs*MILLISEC_PER_YEAR);
				}
				else { // got nothing
					currentYear += timeStepYrs;
					currentTimeMillis = currentTimeMillis+(long)(timeStepYrs*MILLISEC_PER_YEAR);
					continue; // skip the rest
				}
			}

			// sample a rupture
			int nthRup = nthRupRandomSampler.getRandomInt(random);
			int srcIndex = erf.getSrcIndexForNthRup(nthRup);
			int fltSystRupIndex = erf.getFltSysRupIndexForSource(srcIndex);
			double rupMag = magOfNthRups[nthRup];
			double rupArea = fltSysRupSet.getAreaForRup(fltSystRupIndex);
//			mag_ForEventList.add(rupMag);
			fltSysRupIndexForEventList.add(fltSystRupIndex);

			nthRupAtEpochMap.put(eventTimeMillis,nthRup);
			totExpRateAtEventTimeList.add(totalRate); // assumed constant since last event
			obsRupRateArray[nthRup] += 1;
			numSimulatedRups+=1;
			simMFD.addResampledMagRate(rupMag, 1.0, true);
			obsMoRate += MagUtils.magToMoment(rupMag);

			// compute and save the normalized rup recurrence interval if all sections had date of last
			double aveNormRI;
			List<Integer> sectIndicesForRup = fltSysRupSet.getSectionsIndicesForRup(fltSystRupIndex);
			if(aveNormTimeSinceLast) {	// average time since last
				aveNormRI = getAveNormTimeSinceLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect, longTermPartRateForSectArray);
			}
			else {
				long aveDateOfLastMillis = getAveDateOfLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect);
				if(aveDateOfLastMillis != Long.MIN_VALUE) {
					double timeSinceLast = (eventTimeMillis-aveDateOfLastMillis)/MILLISEC_PER_YEAR;
					aveNormRI = timeSinceLast/aveCondRecurIntervalForFltSysRups[fltSystRupIndex];
				}
				else
					aveNormRI = Double.NaN;
			}		
			normRI_ForEventList.add(aveNormRI);
			

			// write event info out UPDATE THIS
			try {
				if(resultsDir != null) eventFileWriter.write(nthRup+"\t"+fltSystRupIndex+"\t"+
						(currentYear+timeToNextInYrs)+"\t"+eventTimeMillis+"\t"+aveNormRI+
						"\t"+rupMag+"\t"+rupArea+"\n");
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			// compute things that are function of subsection
			int[] sectID_Array = sectIndexArrayForSrcList.get(erf.getSrcIndexForFltSysRup(fltSystRupIndex));
			int numSectInRup=sectID_Array.length;
			double mag = fltSysRupSet.getMagForRup(erf.getFltSysRupIndexForNthRup(nthRup));
			double slips[] = getAveSlipOnSectionsForRup(erf, nthRup, mag, numSectInRup);

			if (makePlots) {  // only do this if plots being made
				HistogramFunction sumRI_AlongHist = new HistogramFunction(normRI_AlongStrike.getMinX(), normRI_AlongStrike.getMaxX(), normRI_AlongStrike.getNumX());
				HistogramFunction numRI_AlongHist = new HistogramFunction(normRI_AlongStrike.getMinX(), normRI_AlongStrike.getMaxX(), normRI_AlongStrike.getNumX());
				int ithSectInRup=0;
				for(int sect : sectID_Array) {
					obsSectSlipRateArray[sect] += slips[ithSectInRup];
					long timeOfLastMillis = dateOfLastForSect[sect];
					if(timeOfLastMillis != Long.MIN_VALUE) {
						double normYrsSinceLast = ((eventTimeMillis-timeOfLastMillis)/MILLISEC_PER_YEAR)*longTermPartRateForSectArray[sect];
						normalizedSectRecurIntervals.add(normYrsSinceLast);
						//						if(sect == testSectionIndex)
						//							normalizedSectRecurIntervalsForTestSect.add(normYrsSinceLast);;
//						if(numAperValues>0)
//							normalizedSectRecurIntervalsMagDepList.get(magDepAperiodicity.getAperIndexForRupMag(rupMag)).add(normYrsSinceLast);

						double normDistAlong = ((double)ithSectInRup+0.5)/(double)numSectInRup;
						sumRI_AlongHist.add(normDistAlong, normYrsSinceLast);
						numRI_AlongHist.add(normDistAlong, 1.0);
					}
					ithSectInRup += 1;
				}
				// now put above averages in normRI_AlongStrike
				if(numSectInRup>10) {
					for(int i =0;i<sumRI_AlongHist.size();i++) {
						double num = numRI_AlongHist.getY(i);
						if(num > 0) {
							normRI_AlongStrike.set(sumRI_AlongHist.getX(i), sumRI_AlongHist.getY(i)/num, 1.0);
						}
					}				
				}
			}


			// reset last event time and increment simulated/obs rate on sections
			for(int sect:sectIndexArrayForSrcList.get(srcIndex)) {
				dateOfLastForSect[sect] = eventTimeMillis;
				probModel.setSectDOLE(sect, eventTimeMillis);
				obsSectRateArray[sect] += 1.0; // add the event

				if(sectionRupturedDuringSim[sect]==false) // first time ruptured
					numSectThatRuptured += 1;
				sectionRupturedDuringSim[sect] = true;

				if(rupMag<7.3)
					obsSectRateArrayMlt7pt3[sect] += 1;
				else
					obsSectRateArrayMgt7pt3[sect] += 1;
			}

			// increment time
			currentYear += timeToNextInYrs;
			currentTimeMillis = eventTimeMillis;
		}
		
		double simLoopTimeInMin = ((double)(System.currentTimeMillis()-startRunTime)/(1000.0*60.0));
		if(verbose) System.out.println("Simulation loop took "+simLoopTimeInMin+" min\n");
		if(resultsDir != null) infoString +="Simulation loop took "+(float)simLoopTimeInMin+" min\n\n";
		
		numSectThatRuptured=0;
		for(int s=0;s<sectionRupturedDuringSim.length;s++)
			if(sectionRupturedDuringSim[s])
				numSectThatRuptured += 1;
			else {
				if(verbose)
					System.out.println("Section "+s+" never ruptured; targetRate="+(float)longTermPartRateForSectArray[s]+"; name = "+fltSysRupSet.getFaultSectionData(s).getName());
			}

		if (verbose) System.out.println("Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+"\n");
		if(resultsDir != null) {
			infoString += "Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+
					" ("+(numSections-numSectThatRuptured)+" sections didn't rupture)\n\n";
			infoString += "maxTotalRate="+ maxTotalRate+" at year="+maxTotalRateYr+
					"\nmaxAveGainAtTime="+maxAveGainAtTime+"at year="+maxAveGainAtTimeYr+"\n\n";
			try {
				eventFileWriter.close();
			} catch (IOException e2) {
				// TODO Auto-generated catch block
				e2.printStackTrace();
			}
		}
		
		// check that no time to next were zero millisec
		if(nthRupAtEpochMap.size() != totExpRateAtEventTimeList.size())
			throw new RuntimeException("nthRupAtEpochMap.size() != totExpRateAtEventTime.size()");

		// write section date of last file if not null
		if(outputTimesinceLastMillisFileName != null && resultsDir != null)
			writeSectTimeSinceLastEventToFile(resultsDir, outputTimesinceLastMillisFileName, currentTimeMillis, dateOfLastForSect);
		
		if(resultsDir != null) {
			simMFD.scale(1.0/numYears);
			simMFD.setName("Simulated MFD");
			obsMoRate /= numYears;
			double obsTotRate = simMFD.getTotalIncrRate();
			double rateRatio = obsTotRate/targetMFD.getTotalIncrRate();
			String infoString2 = "total rate = "+(float)obsTotRate+" (ratio="+(float)rateRatio+")";
			double obsTotRateAbove6pt7 = simMFD.getCumRate(6.75);
			double rateAbove6pt7_Ratio = obsTotRateAbove6pt7/targetMFD.getCumRate(6.75);
			infoString2 += "\ntotal rate >= 6.7 = "+(float)obsTotRateAbove6pt7+" (ratio="+(float)rateAbove6pt7_Ratio+")";
			double moRateRatio = obsMoRate/origTotMoRate;
			infoString2 += "\ntotal MoRate = "+(float)obsMoRate+" (ratio="+(float)moRateRatio+")";
			simMFD.setInfo(infoString2);
			
			infoString += "\n\nSimulationStats:\n";
			infoString += "totRate\tratio\ttotRateM>=6.7\tratio\ttotMoRate\tratio\n";
			infoString += (float)obsTotRate+"\t"+(float)rateRatio+"\t"+(float)obsTotRateAbove6pt7+"\t"+(float)rateAbove6pt7_Ratio+"\t"+(float)obsMoRate+"\t"+(float)moRateRatio;			
		}
		
		
		// ******************* Plot Making and finish infoString **************************
		if(makePlots) {
			
			// MFDs
			ProbModelsPlottingUtils.writeMFD_ComprisonPlot(targetMFD, simMFD, plotsDir);
			
			
			// make normalized rup recurrence interval plots
			ArrayList<Double> normalizedRupRecurIntervals = new ArrayList<Double>();
//			ArrayList<ArrayList<Double>> normalizedRupRecurIntervalsMagDepList = new ArrayList<ArrayList<Double>>();
			HashMap<Double,ArrayList<Double>> normalizedRupRecurIntervalsAperMap = new HashMap<Double,ArrayList<Double>>();

			for(double aper:aperValuesList) {
				normalizedRupRecurIntervalsAperMap.put(aper, new ArrayList<Double>());
			}
			for(int e=0;e<normRI_ForEventList.size();e++) {
				double normRI = normRI_ForEventList.get(e);
				if(Double.isNaN(normRI))
					continue;
				normalizedRupRecurIntervals.add(normRI);
				if(numAperValues>0) {
//					double rupMag = mag_ForEventList.get(e);
//					double aper = aperModel.getRuptureAperiodicity(rupSetIndex);
					double aper = aperModel.getRuptureAperiodicity(fltSysRupIndexForEventList.get(e));
					normalizedRupRecurIntervalsAperMap.get(aper).add(normRI);

				}
			}
			double aper=Double.NaN;
			if(numAperValues==1)
				aper=aperValuesList.get(0);	// only one value, so include for comparison
			infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedRupRecurIntervals, aper, 
					plotsDir, "Normalized Rupture RIs", "normalizedRupRecurIntervals");		
			// now mag-dep:
			if(numAperValues >1) {
				for(double aperVal : aperValuesList) {
					String label = "aper="+aperVal;
					infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedRupRecurIntervalsAperMap.get(aperVal), aperVal, 
							plotsDir, "Norm Rup RIs; "+label, "normRupRecurIntsForTargetAper"+aperVal);
				}
			}
			
			
			// make normalized section recurrence interval plots
			aper=Double.NaN;
			infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedSectRecurIntervals, aper, 
					plotsDir, "Normalized Section RIs", "normalizedSectRecurIntervals");		
			// now mag-dep:
//			if(numAperValues >1) {
//				for(int i=0;i<numAperValues;i++) {
//					String label = magDepAperiodicity.getMagDepAperInfoString(i);
//					infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedSectRecurIntervalsMagDepList.get(i), aperValues[i], 
//							otherPlotsDir, "Norm Sect RIs; "+label, "normSectRecurIntsForMagRange"+i);
//				}
//			}
			
			// write infoString
			if(resultsDir != null) {
				FileWriter info_fr;
				try {
					info_fr = new FileWriter(new File(resultsDir,"INFO.txt"));
					info_fr.write(infoString);
					info_fr.close();
				} catch (IOException e1) {
					e1.printStackTrace();
				}			
			}
			
			// plot long-term rate versus time
			ProbModelsPlottingUtils.writeSimExpRateVsTime(nthRupAtEpochMap,totExpRateAtEventTimeList,
					totalLongTermRate, plotsDir);

			
			// plot simulated versus imposed rup rates
			for(int i=0;i<obsRupRateArray.length;i++) {
				obsRupRateArray[i] = obsRupRateArray[i]/numYears;
			}
			ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermRateOfNthRups, obsRupRateArray, plotsDir, "simVsImposedRupRates", 
					"", "Imposed Rup Rate (/yr)", "Simulated Rup Rate (/yr)",Double.NaN,Double.NaN);


			// plot observed versus imposed section slip rates
			for(int i=0;i<obsSectSlipRateArray.length;i++) {
				obsSectSlipRateArray[i] = obsSectSlipRateArray[i]/numYears;
			}
			ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermSlipRateForSectArray, obsSectSlipRateArray, plotsDir, "simVsImposedSectionSlipRates", 
					"", "Imposed Slip Rate (mm/yr)", "Simulated Slip Rate (mm/yr)", Double.NaN,Double.NaN);				

			
			// plot observed versus imposed section rates
			double[] numObsOnSectionArray = new double[obsSectRateArray.length];
			for(int i=0;i<obsSectRateArray.length;i++) {
				numObsOnSectionArray[i] = obsSectRateArray[i];
				obsSectRateArray[i] = obsSectRateArray[i]/numYears;
				obsSectRateArrayMlt7pt3[i] = obsSectRateArrayMlt7pt3[i]/numYears;
				obsSectRateArrayMgt7pt3[i] = obsSectRateArrayMgt7pt3[i]/numYears;
			}
			ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermPartRateForSectArray, obsSectRateArray, plotsDir, "simVsImposedSectionPartRates", 
					"", "Imposed Sect Part Rate (/yr)", "Simulated Sect Part Rate (/yr)", Double.NaN,Double.NaN);
			
			
			// write map of observed vs imposed section rates
			double[] sectPartRateRatioArray = new double[obsSectRateArray.length];
			double[] sectPartRateRatioSigmaArray = new double[obsSectRateArray.length];
			for(int i=0;i<obsSectRateArray.length;i++) {
				sectPartRateRatioArray[i] = obsSectRateArray[i]/longTermPartRateForSectArray[i];
				// Compute sigma as Poisson hi to low confidence bound ratio divided by 4.0.
				// use this for now until PoissonRateFromNinT_Calc is moved out of scratch 
				double alpha = numObsOnSectionArray[i]+1;  // also called shape parameter
				double beta = 1d/numYears; 
				GammaDistribution gd = new GammaDistribution(alpha,beta);
				sectPartRateRatioSigmaArray[i] = (gd.inverseCumulativeProbability(0.975)/gd.inverseCumulativeProbability(0.025))/4.0;;
//				sectPartRateRatioSigmaArray[i] = 0d;
//				double low95bound = PoissonRateFromNinT_Calc.getRateForCumulativeProb(numObsOnSectionArray[i], numYears, 0.025);
//				double hi95bound = PoissonRateFromNinT_Calc.getRateForCumulativeProb(numObsOnSectionArray[i], numYears, 0.975);
//				sectPartRateRatioSigmaArray[i] = (hi95bound/low95bound)4.0;
			}
			ProbModelsPlottingUtils.writeMapOfSimOverTargetPartRates (sectPartRateRatioArray, sectPartRateRatioSigmaArray, 
					fltSysRupSet.getFaultSectionDataList(), plotsDir);
			
			// Turn max source gains into section gains and make map if not poisson
			if(!isPoisson) {
				double[] sectGainsAtMaxRateTimeArray = new double[fltSysRupSet.getNumSections()];
				double[] sectGainsAtMaxAveGainTimeArray = new double[fltSysRupSet.getNumSections()];
				double[] sectLongTermRateArray = new double[fltSysRupSet.getNumSections()];
				for(int n=0; n<erf.getTotNumRupsFromFaultSystem();n++) {
					int[] sectID_Array = sectIndexArrayForSrcList.get(erf.getSrcIndexForNthRup(n));
					for(int s:sectID_Array) {
						// max rate time
						sectGainsAtMaxRateTimeArray[s] += longTermRateOfNthRups[n] * srcGainsAtMaxRateTimeArray[erf.getSrcIndexForNthRup(n)];
						// max ave gain time
						sectGainsAtMaxAveGainTimeArray[s] += longTermRateOfNthRups[n]*srcGainsAtMaxAveGainTimeArray[erf.getSrcIndexForNthRup(n)];
						sectLongTermRateArray[s] += longTermRateOfNthRups[n];
					}
				}
				// turn them into gains
				for(int s=0;s<sectLongTermRateArray.length;s++) {
					sectGainsAtMaxRateTimeArray[s] /= sectLongTermRateArray[s];
					sectGainsAtMaxAveGainTimeArray[s] /= sectLongTermRateArray[s];
				}
				ProbModelsPlottingUtils.writeMapOfSectionProbGains(sectGainsAtMaxRateTimeArray, fltSysRupSet.getFaultSectionDataList(), plotsDir, "mapSectGainsAtMaxRateTime");
				ProbModelsPlottingUtils.writeMapOfSectionProbGains(sectGainsAtMaxAveGainTimeArray, fltSysRupSet.getFaultSectionDataList(), plotsDir, "mapSectGainsAtMaxAveGainTime");
			}

			
			// write section rates with names
			FileWriter sectRates_fw;
			try {
				sectRates_fw = new FileWriter(new File(resultsDir,"/obsVsImposedSectionPartRates.txt"));
				sectRates_fw.write("sectID\timposedRate\tsimulatedRate\tsimOverImpRateRatio\tsectName\n");
				for(int i=0;i<fltSysRupSet.getNumSections();i++) {
					FaultSection fltData = fltSysRupSet.getFaultSectionData(i);
					sectRates_fw.write(fltData.getSectionId()+"\t"+longTermPartRateForSectArray[i]+"\t"+obsSectRateArray[i]+
							"\t"+sectPartRateRatioArray[i]+"\t"+fltData.getName()+"\n");
				}
				sectRates_fw.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			
			
			// OTHER PLOTS BELOW (in "otherPlots" dir)

			
			// plot ave norm RI along strike
			ProbModelsPlottingUtils.writeNormRI_SlongStrikePlot(normRI_AlongStrike, otherPlotsDir);

			// make aveGainVsMagFunction, weighted by rupture rate.  note that these are for all 
			//ruptures (not just the one that occurred at each event time)
			for(int i=0;i<aveRupProbGainArray.length;i++) aveRupProbGainArray[i] /= numSimulatedRups;
			ProbModelsPlottingUtils.writeAveGainVsMagHistPlot(magOfNthRups, aveRupProbGainArray, longTermRateOfNthRups, otherPlotsDir);

			// plot ave, min, and max gain of each rupture vs mag (averaged over simulation regardless 
			//of whether the event occurred).  Not sure how useful this is. It's a big file
			ProbModelsPlottingUtils.writeRupGainMeanMinMaxVsMagPlot(magOfNthRups, aveRupProbGainArray, minRupProbGainArray, maxRupProbGainArray, otherPlotsDir);
	
			// write out these gains and other stuff (big file)
			FileWriter gain_fr;
			try {
				gain_fr = new FileWriter(new File(otherPlotsDir,"aveRupGainData.txt"));
				gain_fr.write("nthRupIndex\taveRupGain\tminRupGain\tmaxRupGain\trupMag\trupLongTermRate\trupCondRI\trupName\n");
				for(int i=0;i<aveRupProbGainArray.length;i++) {
					gain_fr.write(i+"\t"+aveRupProbGainArray[i]+"\t"+minRupProbGainArray[i]+"\t"+maxRupProbGainArray[i]+"\t"+magOfNthRups[i]
							+"\t"+longTermRateOfNthRups[i]+"\t"+aveCondRecurIntervalForFltSysRups[erf.getFltSysRupIndexForNthRup(i)]+"\t"+
							erf.getSource(erf.getSrcIndexForNthRup(i)).getName()+"\n");
				}
				gain_fr.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}	
			
			
			// this makes mag dependent section participation stuff; not sure it's still useful
			ArrayList<String> outStringList = new ArrayList<String>();
			int numSect=fltSysRupSet.getNumSections();
			double[] targetSectRateArrayMlt7pt3 = new double[numSect];
			double[] targetSectRateArrayMgt7pt3 = new double[numSect];
			for(int s=0;s<numSect;s++) {
				double partRateMlow=0;
				double partRateMhigh=0;
				for (int r : fltSysRupSet.getRupturesForSection(s)) {
					double mag = fltSysRupSet.getMagForRup(r);
					if(mag<7.3)
						partRateMlow += fltSysSolution.getRateForRup(r);
					else
						partRateMhigh = fltSysSolution.getRateForRup(r);
				}
				targetSectRateArrayMlt7pt3[s]=partRateMlow;
				targetSectRateArrayMgt7pt3[s]=partRateMhigh;
				outStringList.add(s+"\t"+obsSectRateArray[s]+"\t"+longTermPartRateForSectArray[s]+"\t"+
						(obsSectRateArray[s]/longTermPartRateForSectArray[s])+"\t"+
						targetSectRateArrayMlt7pt3[s]+"\t"+
						obsSectRateArrayMlt7pt3[s]+"\t"+
						obsSectRateArrayMlt7pt3[s]/targetSectRateArrayMlt7pt3[s]+"\t"+
						targetSectRateArrayMgt7pt3[s]+"\t"+
						obsSectRateArrayMgt7pt3[s]+"\t"+
						obsSectRateArrayMgt7pt3[s]/targetSectRateArrayMgt7pt3[s]+"\t"+
						fltSysRupSet.getFaultSectionData(s).getName()+"\n");
			}
			File dataFile = new File(resultsDir,"magDepSectRates_magBelowAndAbove7pt3.txt");
			try {
				FileWriter fileWriter = new FileWriter(dataFile);
				fileWriter.write("secID\tsimRate\ttargetRate\tratio\ttargetLowMrate\tsimLowMrate\tlowRation\ttargetHighMrate\tsimHighMrate\thighRatio\tsectName\n");
				for(String line:outStringList) {
					fileWriter.write(line);
				}
				fileWriter.close();
			}catch(Exception e) {
				e.printStackTrace();
			}
			ProbModelsPlottingUtils.writeSimOverImposedVsImposedSecPartRates_Plot(otherPlotsDir, targetSectRateArrayMlt7pt3, 
					obsSectRateArrayMlt7pt3, 10d, numYears, "Mlt7pt3");
			ProbModelsPlottingUtils.writeSimOverImposedVsImposedSecPartRates_Plot(otherPlotsDir, targetSectRateArrayMgt7pt3, 
					obsSectRateArrayMgt7pt3, 10d, numYears, "Mgt7pt3");
		}
		
		
		if(verbose) System.out.println("INFO STRING:\n\n"+infoString);
	}

	
	private static void readSectTimeSinceLastEventFromFile(String fileName, long currentTimeMillis, long[] dateOfLastForSect) {
		
		try {
			List<String> fileLines = Files.readLines(new File(fileName), Charset.defaultCharset());
			int s=0;
			int numBad=0;
			for (String line:fileLines) {
				String[] st = StringUtils.split(line,"\t");
				int sectIndex = Integer.valueOf(st[0]);
				long timeSince = Long.valueOf(st[1]);
				if(timeSince != Long.MIN_VALUE) {
					dateOfLastForSect[s] = currentTimeMillis-timeSince;
//					dateOfLastForSect[s] = Long.MIN_VALUE;
				}
				else {
					dateOfLastForSect[s] = Long.MIN_VALUE;
					numBad +=1;
				}
				if(s != sectIndex)
					throw new RuntimeException("bad index");
				s+=1;
			}
			int percBad = (int)Math.round(100.0*(double)numBad/(double)dateOfLastForSect.length);
			System.out.println(numBad+" sections out of "+dateOfLastForSect.length+" had no date of last event in input file ("+percBad+"%)");
		} catch (Exception e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
	}

	
	private static void writeSectTimeSinceLastEventToFile(File outputDir, String fileName, long currentTimeMillis, long[] dateOfLastForSect) {		
		if(!outputDir.exists())
			outputDir.mkdir();
		File dataFile = new File(outputDir,fileName);
		try {
			FileWriter fileWriter = new FileWriter(dataFile);
			int numBad=0;
			for(int i=0; i<dateOfLastForSect.length;i++) {
				// time since last millis
				if(dateOfLastForSect[i] != Long.MIN_VALUE) {
					long timeSince = currentTimeMillis-dateOfLastForSect[i];	// ti
					if(timeSince < 0) {
						if(timeSince > -TimeDepUtils.MILLISEC_PER_YEAR) {
							System.out.println("Converting slightly negative time since last ("+timeSince+") to zero");
							timeSince=0;
						}
						else {
							throw new RuntimeException("bad time since last");
						}
					}
					fileWriter.write(i+"\t"+timeSince+"\n");					
				}
				else {
					fileWriter.write(i+"\t"+Long.MIN_VALUE+"\n");
					numBad+=1;
				}
			}
			fileWriter.close();
			int percBad = (int)Math.round(100.0*(double)numBad/(double)dateOfLastForSect.length);
			System.out.println(numBad+" sections out of "+dateOfLastForSect.length+" had no date of last event in output file ("+percBad+"%)");
			
		}catch(Exception e) {
			e.printStackTrace();
		}
	}

	
	private static double[] getAveSlipOnSectionsForRup(TimeDepFaultSystemSolutionERF erf, int nthRup, double rupMag, int numSectInRup) {
		double slips[];
		FaultSystemRupSet fltSysRupSet = erf.getSolution().getRupSet();
		if(fltSysRupSet instanceof InversionFaultSystemRupSet) {
			slips = ((InversionFaultSystemRupSet) fltSysRupSet).getSlipOnSectionsForRup(erf.getFltSysRupIndexForNthRup(nthRup));
		}
		else {	// apply ave to all sections
			//					double mag = fltSysRupSet.getMagForRup(erf.getFltSysRupIndexForNthRup(nthRup));
			double area = fltSysRupSet.getAreaForRup(erf.getFltSysRupIndexForNthRup(nthRup));
			if(area==0)
				throw new RuntimeException("area=0");
			double aveSlip = FaultMomentCalc.getSlip(area, MagUtils.magToMoment(rupMag));
			slips = new double[numSectInRup];
			for(int i=0;i<slips.length;i++)
				slips[i]=aveSlip;
		}
		return slips;
	}

	
	
	/**
	 * This method returns average last-event date (epoch milliseconds) averaged over all fault sections
	 * (and weighted by section area). 
	 * 
	 *  Long.MIN_VALUE is returned if any of the fault sections lack a date of last event 
	 *  
	 * @param sectIndicesForRup
	 * @param presentTimeMillis if non null, only fault sections with date of last events at or before the present
	 * time in milliseconds will be considered.
	 */
	public static long getAveDateOfLastEventWhereKnown(List<Integer> sectIndicesForRup, Long presentTimeMillis, double[] sectionArea, long[] dateOfLastForSect) {
//		System.out.println("getAveDateOfLastEventWhereKnown");
		double totRupArea=0;
		boolean allSectionsHadDateOfLast = true;
		double sumDateOfLast = 0;
		for(int s:sectIndicesForRup) {
			long dateOfLast = dateOfLastForSect[s];
			double area = sectionArea[s];
			totRupArea+=area;
			if(dateOfLast != Long.MIN_VALUE && (presentTimeMillis == null || dateOfLast <= presentTimeMillis)) {
//System.out.println("dateOfLast="+dateOfLast+"; presentTimeMillis="+presentTimeMillis);
				sumDateOfLast += (double)dateOfLast*area;
			}
			else {
				allSectionsHadDateOfLast = false;
			}
		}
		if(allSectionsHadDateOfLast)
			return Math.round(sumDateOfLast/totRupArea);  // epoch millis
		else
			return Long.MIN_VALUE;
	}
	
	
	/**
	 * This method returns normalized time since last event (timeSince/meanRecurInt)
	 *  averaged over fault sections (and weighted by section area). 
	 * 
	 *  Double.NaN is returned if any one of the fault sections lack a date of last event 
	 *  
	 * @param sectIndicesForRup
	 * @param presentTimeMillis - present time in epoch milliseconds
	 */
	public static double getAveNormTimeSinceLastEventWhereKnown(List<Integer> sectIndicesForRup, long presentTimeMillis, 
			double[] sectionArea, long[] dateOfLastForSect, double[] longTermPartRateForSectArray) {
		
		double totRupArea=0;
		boolean allSectionsHadDateOfLast = true;
		double sumNormTimeSinceLast = 0;
		for(int s : sectIndicesForRup) {
			long dateOfLast = dateOfLastForSect[s];
			double area = sectionArea[s];
			totRupArea+=area;
// System.out.println("dateOfLast="+dateOfLast+"; presentTimeMillis="+presentTimeMillis);
			if(dateOfLast != Long.MIN_VALUE && dateOfLast <= presentTimeMillis) {
				sumNormTimeSinceLast += area*((double)(presentTimeMillis-dateOfLast)/MILLISEC_PER_YEAR)*longTermPartRateForSectArray[s];
			}
			else {
				allSectionsHadDateOfLast = false;
			}
		}
		if(allSectionsHadDateOfLast)
			return sumNormTimeSinceLast/totRupArea; 
		else {
			return Double.NaN;
		}
	}
	
	
	
	public static void simulateEventsFast(TimeDepFaultSystemSolutionERF erf, String inputTimeSinceLastMillisFileName, String outputTimesinceLastMillisFileName, 
			double numYears, File resultsDir, long randomSeed, boolean verbose, boolean makePlots, double timeStepYrs) {

		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		boolean isPoisson = false;
		
		boolean aveRecurIntervals=true;
		boolean aveNormTimeSinceLast=true;
		if (probModel instanceof UCERF3_ProbabilityModel) {			// U3 calculation type
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			aveRecurIntervals = u3ProbModel.getAveragingTypeChoice().isAveRI();
			aveNormTimeSinceLast = u3ProbModel.getAveragingTypeChoice().isAveNTS();
			u3ProbModel.setSaveDebugInfo(false); // safer in simulation mode
		} else if (probModel instanceof WG02_ProbabilityModel) {
//			WG02_ProbabilityModel wgProbModel = (WG02_ProbabilityModel)probModel; // not needed
			aveRecurIntervals = true; 		// default value applied by WG02
			aveNormTimeSinceLast = true;

		}
		else if (probModel instanceof FSS_ProbabilityModel.Poisson) {
			isPoisson = true;
			if(verbose) System.out.println("isPoisson = "+isPoisson);
		}
		else
			throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());
		
		FaultSystemSolution fltSysSolution=erf.getSolution();
		FaultSystemRupSet fltSysRupSet = fltSysSolution.getRupSet();
		int numSections = fltSysRupSet.getNumSections();
		double[] longTermPartRateForSectArray = probModel.getSectLongTermPartRates(); // this is a duplicate
		double[] sectionArea = fltSysRupSet.getAreaForAllSections();
		long[] dateOfLastForSect = probModel.getSectDOLE(); // this is a duplicate
		
		if(resultsDir != null) {
			if(!resultsDir.exists()) 
				resultsDir.mkdir();
		}
		else throw new RuntimeException("resultsDir is null'");
		
		RandomDataGenerator randomDataSampler = new RandomDataGenerator();
		randomDataSampler.reSeed(randomSeed);  // for reproducibility; this is for the next event time
		Random random = new Random(randomSeed);  // this is for the nth rup sampler
		
		// INTIALIZE THINGS:
		double[] srcGainsAtMaxRateTimeArray=null, srcGainsAtMaxAveGainTimeArray=null;
		double maxTotalRateAtTime=-1;
		double maxAveGainAtTime=-100;
		double maxTotalRateAtTimeYr=-1, maxAveGainAtTimeYr=1;

		// set original start time and total duration
		long origStartTimeMillis = 0; // 1970 for Poisson simulations
		if(!isPoisson)
			origStartTimeMillis = erf.getTimeSpan().getStartTimeInMillis();
		double origStartYear = ((double)origStartTimeMillis)/MILLISEC_PER_YEAR+1970.0;
		
		String infoString = "Information for simulation results in this directory ("+resultsDir+")\n\n"+
				"inputTimeSinceLastMillisFileName = "+inputTimeSinceLastMillisFileName+"\n\n"+
				"outputTimesinceLastMillisFileName = "+outputTimesinceLastMillisFileName+"\n\n"+
				"simulation duration = "+numYears+" (years)\n\n"+
				"simulation start time: "+origStartTimeMillis+ " millis ("+origStartYear+" yrs)\n\n"+
				"randomSeed = "+randomSeed+"\n\n"+
				"ERF Parameters:\n\n";
		String tempString1 = "\t"+erf.getAdjustableParameterList().getParameterListMetadataString();
		String tempString2 = "\n\t"+erf.getAdjustableParameterList().getParameter(
				TimeDepFaultSystemSolutionERF.PROB_MODEL_PARAM_NAME).getValue().toString();
		String erfParamMetadataString = tempString1.replace(";", "\n\t")+tempString2.replace(";", "\n\t");
		infoString += erfParamMetadataString;
				
		// initialize tracking of sections that had one or more ruptures in simulation
		boolean[] sectionRupturedDuringSim = new boolean[numSections];  
		for(int i=0;i<numSections;i++)
			sectionRupturedDuringSim[i] = false;
		
		// this is for writing out simulated events that occur
		FileWriter eventFileWriter=null;
		eventFileWriter=null;
		try {
			eventFileWriter = new FileWriter(resultsDir+"/sampledEventsData.txt");
			eventFileWriter.write("nthRupIndex\tfssRupIndex\tyear\tepoch\tnormRupRI\trupMag\trupArea\ttotalRateAtThisTime\n");
		} catch (IOException e1) {
			e1.printStackTrace();
		}			

		// temporarily set the forecast as Poisson, to get long-term rates, & no background 
		IncludeBackgroundOption includeBackgroundOption = (IncludeBackgroundOption)erf.getParameter(IncludeBackgroundParam.NAME).getValue();
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
		erf.updateForecast();

		// fill in totalRate, longTermRateOfNthRups, magOfNthRups, and longTermSlipRateForSectArray
		double totalLongTermRate=0, longTermMoRate=0;
		IntegerPDF_FunctionSampler nthRupRandomSampler = new IntegerPDF_FunctionSampler(erf.getTotNumRups());
		double[] longTermRateOfNthRups = new double[erf.getTotNumRups()];	// this will include any aftershock reductions
		double[] magOfNthRups = new double[erf.getTotNumRups()];
		for(int nthRup=0; nthRup<erf.getTotNumRups(); nthRup++) {
			ProbEqkRupture rup = erf.getNthRupture(nthRup);
			double rate = rup.getMeanAnnualRate(erf.getTimeSpan().getDuration());
			longTermRateOfNthRups[nthRup] = rate;
			totalLongTermRate += rate;
			longTermMoRate += rate*MagUtils.magToMoment(rup.getMag());
			magOfNthRups[nthRup] = rup.getMag();
			nthRupRandomSampler.set(nthRup, rate);
		}
		
		
		double simDuration = 1/totalLongTermRate;  // used to compute next prob gain
		
		infoString += "\n\nTotal long-term rate (per year) = "+totalLongTermRate+"\n\n";

		// reset ERF to original state
		erf.setCustomProbabilityModel(probModel);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(includeBackgroundOption);
		erf.updateForecast();
		
		// Make local sectIndexArrayForSrcList for faster simulations
		ArrayList<int[]> sectIndexArrayForSrcList = new ArrayList<int[]>();
		for(int s=0; s<erf.getNumFaultSystemSources();s++) {
			List<Integer> indexList = fltSysRupSet.getSectionsIndicesForRup(erf.getFltSysRupIndexForSource(s));
			int[] indexArray = new int[indexList.size()];
			for(int i=0;i<indexList.size();i++)
				indexArray[i] = indexList.get(i);
			sectIndexArrayForSrcList.add(indexArray);
		}
		
		// set the ave cond recurrence intervals
		double[] aveCondRecurIntervalForFltSysRups = TimeDepUtils.computeAveCondRecurIntervalForFltSysRups(fltSysRupSet, 
				longTermPartRateForSectArray, aveRecurIntervals);
		double minCondRI=Double.MAX_VALUE,maxCondRI=0;
		for(double ri: aveCondRecurIntervalForFltSysRups) {
			if(!Double.isInfinite(ri)) {
				if(ri < minCondRI) minCondRI = ri;
				if(ri > maxCondRI) maxCondRI = ri;
			}
		}
		
		infoString +="minCondRI="+minCondRI+"\nmaxCondRI="+maxCondRI+"\n\n";

		//write out infoString if verbose
		if(verbose) System.out.println(infoString);
		
		// set simulation time
		double currentYear=origStartYear;
		long currentTimeMillis = origStartTimeMillis;

		// read section date of last file if not null
		if(inputTimeSinceLastMillisFileName != null && !isPoisson) {
			readSectTimeSinceLastEventFromFile(inputTimeSinceLastMillisFileName, currentTimeMillis, dateOfLastForSect);
			for(int s=0; s<dateOfLastForSect.length;s++)
				probModel.setSectDOLE(s, dateOfLastForSect[s]);
		}

		// remove date of last event if Poisson - this not neccessary?
		if(isPoisson)
			for(int s=0; s<dateOfLastForSect.length;s++)
				dateOfLastForSect[s]= Long.MIN_VALUE;

		// this is to track progress
		int percDoneThresh=0;
		int percDoneIncrement=5;

		long startRunTime = System.currentTimeMillis();	
		if(verbose) System.out.println("Starting simulation loop");
		
		int numSectThatRuptured = 0;
		int numSimulatedRups=0;
		double totalRateAtTime=totalLongTermRate; // this for Poisson (won't get updated in loop)
		double simMoRate = 0;

		
		while (currentYear<numYears+origStartYear) {
			
			// write progress
			if(verbose) {
				int percDone = (int)Math.round(100*(currentYear-origStartYear)/numYears);
				if(percDone >= percDoneThresh) {
					double timeInMin = ((double)(System.currentTimeMillis()-startRunTime)/(1000.0*60.0));
					int numGoodDateOfLast=0;
					for(long dateOfLast:dateOfLastForSect) {
						if(dateOfLast != Long.MIN_VALUE)
							numGoodDateOfLast+=1;					
					}
					int percentGood = (int)Math.round((100.0*(double)numGoodDateOfLast/(double)dateOfLastForSect.length));
					System.out.println("\n"+percDoneThresh+"% done in "+(float)timeInMin+" minutes"+"; yr="+(float)currentYear+";  % sect with date of last = "+percentGood);	
					System.out.println("\tFraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+"\n");
			
					percDoneThresh += percDoneIncrement;
				}				
			}
			
			// update gains and sampler if not Poisson
			if(!isPoisson) {
				double[] probGainForFaultSystemSource = new double[erf.getNumFaultSystemSources()];
				double aveGainAtTime = 0;
				// first the gains
				for(int s=0;s<erf.getNumFaultSystemSources();s++) {
					int fltSysRupIndex = erf.getFltSysRupIndexForSource(s);
					probGainForFaultSystemSource[s] = probModel.getProbabilityGain(fltSysRupIndex, currentTimeMillis, simDuration);
				}

				// now update totalRate and ruptureSampler (for all rups since start time changed)
				for(int n=0; n<erf.getTotNumRupsFromFaultSystem();n++) {
					double probGain = probGainForFaultSystemSource[erf.getSrcIndexForNthRup(n)];
					aveGainAtTime += probGain;
					double newRate = longTermRateOfNthRups[n] * probGain;	// applied as a rate gain
					nthRupRandomSampler.set(n, newRate);
				}
				
				totalRateAtTime = nthRupRandomSampler.getSumOfY_vals();
				aveGainAtTime /= erf.getNumFaultSystemSources();
				if(maxTotalRateAtTime<totalRateAtTime) {
					srcGainsAtMaxRateTimeArray = probGainForFaultSystemSource;
					maxTotalRateAtTime=totalRateAtTime;
					maxTotalRateAtTimeYr = currentYear;
				}
				if(maxAveGainAtTime<aveGainAtTime) {
					srcGainsAtMaxAveGainTimeArray = probGainForFaultSystemSource;
					maxAveGainAtTime=aveGainAtTime;
					maxAveGainAtTimeYr = currentYear;
				}
			}
			
			// Try sampling an event time
			double timeToNextInYrs;
			long eventTimeMillis;
			if(Double.isNaN(timeStepYrs)) {  // no time step; sample time of next event
				timeToNextInYrs = randomDataSampler.nextExponential(1.0/totalRateAtTime);
				eventTimeMillis = currentTimeMillis + (long)(timeToNextInYrs*MILLISEC_PER_YEAR);
			}
			else {
				long numEvents = randomDataSampler.nextPoisson(totalRateAtTime*timeStepYrs);
				if(numEvents > 0) { // get a random number or events for next time step;  set timeStepYrs low enough that more than 1 event is rare
					timeToNextInYrs = timeStepYrs*randomDataSampler.nextUniform(0, 1);
					eventTimeMillis = currentTimeMillis + (long)(timeToNextInYrs*MILLISEC_PER_YEAR);
				}
				else { // got nothing
					currentYear += timeStepYrs;
					currentTimeMillis = currentTimeMillis+(long)(timeStepYrs*MILLISEC_PER_YEAR);
					continue; // skip the rest
				}
			}

			// sample a rupture
			int nthRup = nthRupRandomSampler.getRandomInt(random);
			int srcIndex = erf.getSrcIndexForNthRup(nthRup);
			int fltSystRupIndex = erf.getFltSysRupIndexForSource(srcIndex);
			double rupMag = magOfNthRups[nthRup];
			double rupArea = fltSysRupSet.getAreaForRup(fltSystRupIndex);
			numSimulatedRups+=1;
			simMoRate += MagUtils.magToMoment(rupMag);

			// compute and save the normalized rup recurrence interval if all sections had date of last
			double aveNormRI;
			List<Integer> sectIndicesForRup = fltSysRupSet.getSectionsIndicesForRup(fltSystRupIndex);
			if(aveNormTimeSinceLast) {	// average time since last
				aveNormRI = getAveNormTimeSinceLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect, longTermPartRateForSectArray);
			}
			else {
				long aveDateOfLastMillis = getAveDateOfLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect);
				if(aveDateOfLastMillis != Long.MIN_VALUE) {
					double timeSinceLast = (eventTimeMillis-aveDateOfLastMillis)/MILLISEC_PER_YEAR;
					aveNormRI = timeSinceLast/aveCondRecurIntervalForFltSysRups[fltSystRupIndex];
				}
				else
					aveNormRI = Double.NaN;
			}		

			// write event info out
			try {
				eventFileWriter.write(nthRup+"\t"+fltSystRupIndex+"\t"+
						(currentYear+timeToNextInYrs)+"\t"+eventTimeMillis+"\t"+aveNormRI+
						"\t"+rupMag+"\t"+rupArea+"\t"+totalRateAtTime+"\n");
			} catch (IOException e1) {
				e1.printStackTrace();
			}

			// reset last event time and increment simulated/obs rate on sections
			for(int sect:sectIndexArrayForSrcList.get(srcIndex)) {
				dateOfLastForSect[sect] = eventTimeMillis;
				probModel.setSectDOLE(sect, eventTimeMillis);

				if(sectionRupturedDuringSim[sect]==false) // first time ruptured
					numSectThatRuptured += 1;
				sectionRupturedDuringSim[sect] = true;
			}

			// increment time
			currentYear += timeToNextInYrs;
			currentTimeMillis = eventTimeMillis;
		}

		// close fileWriter
		try {
			eventFileWriter.close();
		} catch (IOException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		
		double simLoopTimeInMin = ((double)(System.currentTimeMillis()-startRunTime)/(1000.0*60.0));
		if(verbose) System.out.println("Simulation loop took "+simLoopTimeInMin+" min\n");
		infoString +="Simulation loop took "+(float)simLoopTimeInMin+" min\n\n";
		
		numSectThatRuptured=0;
		for(int s=0;s<sectionRupturedDuringSim.length;s++)
			if(sectionRupturedDuringSim[s])
				numSectThatRuptured += 1;

		if (verbose) System.out.println("Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+"\n");
		infoString += "Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+
				" ("+(numSections-numSectThatRuptured)+" sections didn't rupture)\n\n";
		infoString += "maxTotalRate = "+ maxTotalRateAtTime+" at year="+maxTotalRateAtTimeYr+
				"\nmaxAveGainAtTime = "+maxAveGainAtTime+" at year="+maxAveGainAtTimeYr+"\n\n";
		
		// write section date of last file if not null
		if(outputTimesinceLastMillisFileName != null && resultsDir != null)
			writeSectTimeSinceLastEventToFile(resultsDir, outputTimesinceLastMillisFileName, currentTimeMillis, dateOfLastForSect);
		
		simMoRate /= numYears;
		double obsTotRate = numSimulatedRups/numYears;
		double rateRatio = obsTotRate/totalLongTermRate;
		double moRateRatio = simMoRate/longTermMoRate;
		infoString += "\n\nSimulationStats:\n";
		infoString += "totRate\tratio\ttotMoRate\tratio\n";
		infoString += (float)obsTotRate+"\t"+(float)rateRatio+"\t"+(float)simMoRate+"\t"+(float)moRateRatio;			

		FileWriter info_fr;
		try {
			info_fr = new FileWriter(new File(resultsDir,"INFO.txt"));
			info_fr.write(infoString);
			info_fr.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}		
		
		//Write out srcGainsAtMaxRateTimeArray & srcGainsAtMaxAveGainTimeArray
		if(!isPoisson) {
			File dataFile = new File(resultsDir,"srcGainsAtMaxRateAnAveGainTime.txt");
			try {
				FileWriter maxGainsWriter = new FileWriter(dataFile);
				maxGainsWriter.write("fltSysIndex"+"\t"+"srcGainAtMaxRateTime"+"\t"+"srcGainAtMaxAveGainTime"+"\n");					
				for(int i=0; i<srcGainsAtMaxRateTimeArray.length;i++) {
					maxGainsWriter.write(i+"\t"+srcGainsAtMaxRateTimeArray[i]+"\t"+srcGainsAtMaxAveGainTimeArray[i]+"\n");					
				}
			maxGainsWriter.close();
			}catch(Exception e) {
				e.printStackTrace();
			}			
		}
	}

	/**
	 * Remove unused arguments
	 * @param erf
	 * @param inputTimeSinceLastMillisFileName
	 * @param outputTimesinceLastMillisFileName
	 * @param numYears
	 * @param resultsDir
	 * @param randomSeed
	 * @param verbose
	 * @param makePlots
	 * @param timeStepYrs
	 */
	public static void generateSimulationPlots(TimeDepFaultSystemSolutionERF erf, String inputTimeSinceLastMillisFileName, 
			double numYears, File resultsDir, boolean verbose) {

		if(!resultsDir.exists()) 
			throw new RuntimeException(resultsDir+" does not exist");
		
		File plotsDir = new File(resultsDir, "plots");
		if(!plotsDir.exists()) plotsDir.mkdir();
			
		File otherPlotsDir = new File(plotsDir, "otherPlots");
		if(!otherPlotsDir.exists()) otherPlotsDir.mkdir();

		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		UCERF3_ProbabilityModel u3ProbModel = null;
		WG02_ProbabilityModel wgProbModel = null;
		boolean isPoisson = false;
		
		boolean aveRecurIntervals = true; 		// default value applied by WG02
		boolean aveNormTimeSinceLast = true;
		AperiodicityModel aperModel = null;
		if (probModel instanceof UCERF3_ProbabilityModel) {			// U3 calculation type
			u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			aveRecurIntervals = u3ProbModel.getAveragingTypeChoice().isAveRI();
			aveNormTimeSinceLast = u3ProbModel.getAveragingTypeChoice().isAveNTS();
			aperModel = u3ProbModel.getAperiodicityModel();
		} else if (probModel instanceof WG02_ProbabilityModel) {
			wgProbModel = (WG02_ProbabilityModel)probModel;
			aperModel = wgProbModel.getAperiodicityModel();
		}
		else if (probModel instanceof FSS_ProbabilityModel.Poisson) {
			isPoisson = true;
		}
		else
			throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());
		
//		// can also do things via the parameter list
//		ParameterList probModelParams = probModel.getAdjustableParameters();
//		if (probModelParams.containsParameter(RenewalModels.PARAM_NAME))
//			probModelParams.setValue(RenewalModels.PARAM_NAME, RenewalModels.BPT);		

		FaultSystemSolution fltSysSolution=erf.getSolution();
		// this has zeros where events were filtered our by the ERF (mags too low); 
		double[] longTermRateOfFltSysRup = erf.getLongTermRateOfFltSysRupInERF(); // this can be different than what's in the FSS
		
		FaultSystemRupSet fltSysRupSet = fltSysSolution.getRupSet();
		int numRupsInFaultSystem = fltSysRupSet.getNumRuptures();
		int numSections = fltSysRupSet.getNumSections();
		
		double[] longTermPartRateForSectArray = probModel.getSectLongTermPartRates(); // this is a duplicate
		double[] sectionArea = fltSysRupSet.getAreaForAllSections();
		long[] dateOfLastForSect = probModel.getSectDOLE(); // this is a duplicate
		
		// set original start time and total duration
		long origStartTimeMillis = 0; // 1970 for Poisson simulations
		if(!isPoisson)
			origStartTimeMillis = erf.getTimeSpan().getStartTimeInMillis();
		double origStartYear = ((double)origStartTimeMillis)/MILLISEC_PER_YEAR+1970.0;

		String infoString = "Information for simulation results in this directory ("+resultsDir+")\n\n"+
				"inputTimeSinceLastMillisFileName = "+inputTimeSinceLastMillisFileName+"\n\n"+
				"simulation duration = "+numYears+" (years)\n\n"+
				"simulation start time: "+origStartTimeMillis+ " millis ("+origStartYear+" yrs)\n\n"+
				"ERF Parameters:\n\n";
		String tempString1 = "\t"+erf.getAdjustableParameterList().getParameterListMetadataString();
		String tempString2 = "\n\t"+erf.getAdjustableParameterList().getParameter(
				TimeDepFaultSystemSolutionERF.PROB_MODEL_PARAM_NAME).getValue().toString();
		String erfParamMetadataString = tempString1.replace(";", "\n\t")+tempString2.replace(";", "\n\t");
		infoString += erfParamMetadataString;
		
		double maxTotalRate=-1;
				
		TreeMap<Long, Integer> nthRupAtEpochMap = new TreeMap<Long, Integer>();
		ArrayList<Double> totExpRateAtEventTimeList = new ArrayList<Double>();
		ArrayList<Double> normRI_ForEventList = new ArrayList<Double>();
		ArrayList<Double> mag_ForEventList = new ArrayList<Double>();
		ArrayList<Integer> fltSysRupIndexForEventList = new ArrayList<Integer>();

    	ArbDiscrEmpiricalDistFunc_3D normRI_AlongStrike = new ArbDiscrEmpiricalDistFunc_3D(0.05d,0.95d,10);
		double[] obsSectRateArray = new double[numSections];
		double[] obsSectSlipRateArray = new double[numSections];
		double[] obsSectRateArrayMlt7pt3 = new double[numSections];
		double[] obsSectRateArrayMgt7pt3 = new double[numSections];

		double[] obsRupRateArray = new double[erf.getTotNumRups()];
//		double[] aveRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time
//		double[] minRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time
//		double[] maxRupProbGainArray = new double[erf.getTotNumRups()];	// averages the prob gains at each event time

		// initialize tracking of sections that had one or more ruptures in simulation
		boolean[] sectionRupturedDuringSim = new boolean[numSections];  
		for(int i=0;i<numSections;i++)
			sectionRupturedDuringSim[i] = false;

		// temporarily set the forecast as Poisson, to get long-term rates, & no background 
		IncludeBackgroundOption includeBackground = (IncludeBackgroundOption)erf.getParameter(IncludeBackgroundParam.NAME).getValue();
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
		erf.updateForecast();

		// fill in longTermRateOfNthRups, magOfNthRups, and longTermSlipRateForSectArray
		double totalLongTermRate=0;
		IntegerPDF_FunctionSampler nthRupRandomSampler = new IntegerPDF_FunctionSampler(erf.getTotNumRups());
		double[] longTermRateOfNthRups = new double[erf.getTotNumRups()];	// this will include any aftershock reductions
		if(obsRupRateArray.length != longTermRateOfNthRups.length)
			throw new RuntimeException("obsRupRateArray.length="+ obsRupRateArray.length+"\nlongTermRateOfNthRups.length="+longTermRateOfNthRups.length);
		double[] magOfNthRups = new double[erf.getTotNumRups()];
		double[] longTermSlipRateForSectArray = new double[numSections];
		ArrayList<Double> aperValuesList = new ArrayList<Double>();
		for(int nthRup=0; nthRup<erf.getTotNumRups(); nthRup++) {
			ProbEqkRupture rup = erf.getNthRupture(nthRup);
			double rate = rup.getMeanAnnualRate(erf.getTimeSpan().getDuration());
			longTermRateOfNthRups[nthRup] = rate;
			totalLongTermRate += rate;
			magOfNthRups[nthRup] = rup.getMag();
			nthRupRandomSampler.set(nthRup, rate);
			int fltSysIndex = erf.getFltSysRupIndexForNthRup(nthRup);
			// aperiodicities list
			if(!isPoisson) {
				double aper =aperModel.getRuptureAperiodicity(fltSysIndex);
				if(!aperValuesList.contains(aper))
					aperValuesList.add(aper);
			}
			// slip rates
			List<Integer> sectIndices = fltSysRupSet.getSectionsIndicesForRup(fltSysIndex);
			double mag = fltSysRupSet.getMagForRup(erf.getFltSysRupIndexForNthRup(nthRup));
			double slips[] = getAveSlipOnSectionsForRup(erf, nthRup, mag, sectIndices.size());
			for(int s=0;s<sectIndices.size();s++) {
				int sectID = sectIndices.get(s);
				longTermSlipRateForSectArray[sectID] += rate*slips[s];
			}					
		}
		int numAperValues = aperValuesList.size();

		infoString += "\n\nTotal long-term rate (per year) = "+totalLongTermRate+"\n\n";
		
		// this is for storing section normalized RIs
		ArrayList<Double> normalizedSectRecurIntervals = new ArrayList<Double>();

		// Norm RIs at paleo sites
		PaleoseismicConstraintData paleoDataMod = erf.getSolution().getRupSet().getModule(PaleoseismicConstraintData.class);
		boolean processPaleoSiteRIs = false;
		HashMap<Integer,ArrayList<Double>> paleoSitesNormRI_List_Map = null;
		HashMap<Integer,String> paleoSitesNameMap = null;
		if(paleoDataMod != null) {
			processPaleoSiteRIs = true;
			paleoSitesNormRI_List_Map = new HashMap<Integer,ArrayList<Double>>();
			paleoSitesNameMap = new HashMap<Integer,String>();
			for(SectMappedUncertainDataConstraint constr : paleoDataMod.getPaleoRateConstraints()) {
				int paleoSectID = constr.sectionIndex;
				paleoSitesNormRI_List_Map.put(paleoSectID, new ArrayList<Double>());
				paleoSitesNameMap.put(paleoSectID, constr.getName());

//				// this was used to make sure the paleo sites were properly mapped
//				RuptureSurface sectSurf =erf.getSolution().getRupSet().getFaultSectionData(paleoSectID).getFaultSurface(1.0);
//				double distTest = sectSurf.getDistanceJB(constr.dataLocation);
//				if(distTest > 3) { // 3 km
//					// check if any other sections are closer
//					double minDist = Double.POSITIVE_INFINITY;
//					for(FaultSection fltSectData : fltSysRupSet.getFaultSectionDataList()) {
//						double dist = fltSectData.getFaultSurface(1.0).getDistanceJB(constr.dataLocation);
//						if(minDist>dist)
//							minDist=dist;
//					}
//					System.out.println((float)distTest+"\t"+constr.getName()+";  minDist to all sections: "+(float)minDist);
//					// throw new RuntimeException("Paleo site not near fault section: dist="+(float)distTest+" for "+constr.getName());
//					System.exit(0);
//				}
//				System.out.println((float)distTest+"\t"+constr.getName());
			}
		}
		
		double simDuration = 1/totalLongTermRate;  // used to compute next prob gain
		
		// make the target MFD - 
		SummedMagFreqDist targetMFD=null;
		double origTotMoRate=Double.NaN;
		targetMFD = ERF_Calculator.getTotalMFD_ForERF(erf, 5.05, 8.95, 40, true);
		origTotMoRate = ERF_Calculator.getTotalMomentRateInRegion(erf, null);
		System.out.println("originalTotalMomentRate: "+origTotMoRate);
		targetMFD.setName("Target MFD");
		String tempString = "total rate = "+(float)targetMFD.getTotalIncrRate();
		tempString += "\ntotal rate >= 6.7 = "+(float)targetMFD.getCumRate(6.75);
		tempString += "\ntotal MoRate = "+(float)origTotMoRate;
		targetMFD.setInfo(tempString);			

		// reset ERF to original state
		erf.setCustomProbabilityModel(probModel);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(includeBackground);
		erf.updateForecast();

		// MFD & MoRate for simulation
		SummedMagFreqDist simMFD = new SummedMagFreqDist(5.05,8.95,40);
		double simMoRate = 0;
		
		// Make local sectIndexArrayForSrcList for faster simulations
		ArrayList<int[]> sectIndexArrayForSrcList = new ArrayList<int[]>();
		for(int s=0; s<erf.getNumFaultSystemSources();s++) {
			List<Integer> indexList = fltSysRupSet.getSectionsIndicesForRup(erf.getFltSysRupIndexForSource(s));
			int[] indexArray = new int[indexList.size()];
			for(int i=0;i<indexList.size();i++)
				indexArray[i] = indexList.get(i);
			sectIndexArrayForSrcList.add(indexArray);
		}
		
		// set the ave cond recurrence intervals
		double[] aveCondRecurIntervalForFltSysRups = TimeDepUtils.computeAveCondRecurIntervalForFltSysRups(fltSysRupSet, 
														longTermPartRateForSectArray, aveRecurIntervals);
		// print minimum and maximum conditional rate of rupture
		double minCondRI=Double.MAX_VALUE,maxCondRI=0;
		for(double ri: aveCondRecurIntervalForFltSysRups) {
			if(!Double.isInfinite(ri)) {
				if(ri < minCondRI) minCondRI = ri;
				if(ri > maxCondRI) maxCondRI = ri;
			}
		}
		infoString +="minCondRI="+minCondRI+"\nmaxCondRI="+maxCondRI+"\n\n";
		
		// set simulation time
		double currentYear=origStartYear;
		long currentTimeMillis = origStartTimeMillis;

		// read section date of last file if not null
		if(inputTimeSinceLastMillisFileName != null && !isPoisson) {
			readSectTimeSinceLastEventFromFile(inputTimeSinceLastMillisFileName, currentTimeMillis, dateOfLastForSect);
			for(int s=0; s<dateOfLastForSect.length;s++)
				probModel.setSectDOLE(s, dateOfLastForSect[s]);
		}

		// remove date of last event if Poisson
		if(isPoisson)
			for(int s=0; s<dateOfLastForSect.length;s++)
				dateOfLastForSect[s]= Long.MIN_VALUE;

		long startRunTime = System.currentTimeMillis();	
		if(verbose) System.out.println("Starting simulation loop");
		
		int numSectThatRuptured = 0;
		int numSimulatedRups=0;
		
		
		// this is for reading simulated events data
		File dataFile = new File(resultsDir+"/sampledEventsData.txt");
		BufferedReader reader=null;
		boolean firstLine = true; // for skipping header
		try {
			reader = new BufferedReader(scratch.UCERF3.utils.UCERF3_DataUtils.getReader(dataFile.toURL()));

			String line;
			while ((line = reader.readLine()) != null) {	
				
				if(firstLine) {
					firstLine = false;
					continue;
				}

				//parse line from file:
				String[] st = StringUtils.split(line,"\t");

				int nthRup = Integer.valueOf(st[0]);
				int fltSystRupIndex = Integer.valueOf(st[1]);
				double eventYear = Double.valueOf(st[2]);
				long eventTimeMillis = Long.valueOf(st[3]);
				double aveNormRI = Double.valueOf(st[4]);
				double rupMag = Double.valueOf(st[5]); 
				double rupArea = Double.valueOf(st[6]); 
				double totalRateAtTime = Double.valueOf(st[7]); 

				int srcIndex = erf.getSrcIndexForNthRup(nthRup);


				//			nthRup+"\t"+fltSystRupIndex+"\t"+
				//			(currentYear+timeToNextInYrs)+"\t"+eventTimeMillis+"\t"+aveNormRI+
				//			"\t"+rupMag+"\t"+rupArea+"\t"+totalRateAtTime+"\n"

				//			double timeToNextInYrs;

				mag_ForEventList.add(rupMag);
				fltSysRupIndexForEventList.add(fltSystRupIndex);
				nthRupAtEpochMap.put(eventTimeMillis,nthRup);
				totExpRateAtEventTimeList.add(totalRateAtTime); // assumed constant since last event
				normRI_ForEventList.add(aveNormRI);
				obsRupRateArray[nthRup] += 1;
				numSimulatedRups+=1;
				simMFD.addResampledMagRate(rupMag, 1.0, true);
				simMoRate += MagUtils.magToMoment(rupMag);

				// USE THIS TO CHECK FOR CONSISTENT ERF SETTINGS
				double aveNormRI_test;
				List<Integer> sectIndicesForRup = fltSysRupSet.getSectionsIndicesForRup(fltSystRupIndex);
				if(aveNormTimeSinceLast) {	// average time since last
					aveNormRI_test = getAveNormTimeSinceLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect, longTermPartRateForSectArray);
				}
				else {
					long aveDateOfLastMillis = getAveDateOfLastEventWhereKnown(sectIndicesForRup, eventTimeMillis, sectionArea, dateOfLastForSect);
					if(aveDateOfLastMillis != Long.MIN_VALUE) {
						double timeSinceLast = (eventTimeMillis-aveDateOfLastMillis)/MILLISEC_PER_YEAR;
						aveNormRI_test = timeSinceLast/aveCondRecurIntervalForFltSysRups[fltSystRupIndex];
					}
					else
						aveNormRI_test = Double.NaN;
				}		
				if(!Double.isNaN(aveNormRI_test) || !Double.isNaN(aveNormRI))
					if((float)aveNormRI_test != (float)aveNormRI)
						throw new RuntimeException("aveNormRI_test != aveNormRI; "+aveNormRI_test+"  vs  "+aveNormRI);

				// compute things that are function of subsection
				int[] sectID_Array = sectIndexArrayForSrcList.get(erf.getSrcIndexForFltSysRup(fltSystRupIndex));
				int numSectInRup=sectID_Array.length;
				double slips[] = getAveSlipOnSectionsForRup(erf, nthRup, rupMag, numSectInRup);

				HistogramFunction sumRI_AlongHist = new HistogramFunction(normRI_AlongStrike.getMinX(), normRI_AlongStrike.getMaxX(), normRI_AlongStrike.getNumX());
				HistogramFunction numRI_AlongHist = new HistogramFunction(normRI_AlongStrike.getMinX(), normRI_AlongStrike.getMaxX(), normRI_AlongStrike.getNumX());
				int ithSectInRup=0;
				for(int sect : sectID_Array) {
					obsSectSlipRateArray[sect] += slips[ithSectInRup];
					long timeOfLastMillis = dateOfLastForSect[sect];
					if(timeOfLastMillis != Long.MIN_VALUE) {
						double normYrsSinceLast = ((eventTimeMillis-timeOfLastMillis)/MILLISEC_PER_YEAR)*longTermPartRateForSectArray[sect];
						normalizedSectRecurIntervals.add(normYrsSinceLast);
						double normDistAlong = ((double)ithSectInRup+0.5)/(double)numSectInRup;
						sumRI_AlongHist.add(normDistAlong, normYrsSinceLast);
						numRI_AlongHist.add(normDistAlong, 1.0);
						// add to list if paleo site
						if(paleoSitesNormRI_List_Map.keySet().contains(sect)) {
							paleoSitesNormRI_List_Map.get(sect).add(normYrsSinceLast);
						}
					}
					ithSectInRup += 1;
				}
				// now put above averages in normRI_AlongStrike
				if(numSectInRup>10) {
					for(int i =0;i<sumRI_AlongHist.size();i++) {
						double num = numRI_AlongHist.getY(i);
						if(num > 0) {
							normRI_AlongStrike.set(sumRI_AlongHist.getX(i), sumRI_AlongHist.getY(i)/num, 1.0);
						}
					}				
				}


				// reset last event time and increment simulated/obs rate on sections
				for(int sect:sectIndexArrayForSrcList.get(srcIndex)) {
					dateOfLastForSect[sect] = eventTimeMillis;
					//				probModel.setSectDOLE(sect, eventTimeMillis);
					obsSectRateArray[sect] += 1.0; // add the event

					if(sectionRupturedDuringSim[sect]==false) // first time ruptured
						numSectThatRuptured += 1;
					sectionRupturedDuringSim[sect] = true;

					if(rupMag<7.3)
						obsSectRateArrayMlt7pt3[sect] += 1;
					else
						obsSectRateArrayMgt7pt3[sect] += 1;
				}

				// increment time
				currentYear = eventYear;
				currentTimeMillis = eventTimeMillis;
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		double simLoopTimeInMin = ((double)(System.currentTimeMillis()-startRunTime)/(1000.0*60.0));
		if(verbose) System.out.println("Loop took "+simLoopTimeInMin+" min\n");
//		if(resultsDir != null) infoString +="Simulation loop took "+(float)simLoopTimeInMin+" min\n\n";
		
		numSectThatRuptured=0;
		// write out sections that didn't rupture
		try {
			FileWriter sectsNotRupturedWriter = new FileWriter(new File(resultsDir,"sectionsThatDidNotRupture.txt"));
			for(int s=0;s<sectionRupturedDuringSim.length;s++)
				if(sectionRupturedDuringSim[s])
					numSectThatRuptured += 1;
				else
					sectsNotRupturedWriter.write("Section "+s+" never ruptured; targetRate="+(float)longTermPartRateForSectArray[s]+"; name = "+fltSysRupSet.getFaultSectionData(s).getName()+"\n");					
			sectsNotRupturedWriter.close();
		}catch(Exception e) {
			e.printStackTrace();
		}			

		if (verbose) System.out.println("Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+"\n");

		infoString += "Final fraction of sections that ruptured: "+ (double)numSectThatRuptured/(double)numSections+
				" ("+(numSections-numSectThatRuptured)+" sections didn't rupture)\n\n";
	
		// check that no time to next were zero millisec
//		if(nthRupAtEpochMap.size() != totExpRateAtEventTimeList.size())
//			throw new RuntimeException("nthRupAtEpochMap.size() != totExpRateAtEventTime.size()");

		simMFD.scale(1.0/numYears);
		simMFD.setName("Simulated MFD");
		simMoRate /= numYears;
		double obsTotRate = simMFD.getTotalIncrRate();
		double rateRatio = obsTotRate/targetMFD.getTotalIncrRate();
		String infoString2 = "total rate = "+(float)obsTotRate+" (ratio="+(float)rateRatio+")";
		double obsTotRateAbove6pt7 = simMFD.getCumRate(6.75);
		double rateAbove6pt7_Ratio = obsTotRateAbove6pt7/targetMFD.getCumRate(6.75);
		infoString2 += "\ntotal rate >= 6.7 = "+(float)obsTotRateAbove6pt7+" (ratio="+(float)rateAbove6pt7_Ratio+")";
		double moRateRatio = simMoRate/origTotMoRate;
		infoString2 += "\ntotal MoRate = "+(float)simMoRate+" (ratio="+(float)moRateRatio+")";
		simMFD.setInfo(infoString2);

		infoString += "\n\nSimulationStats:\n";
		infoString += "totRate\tratio\ttotRateM>=6.7\tratio\ttotMoRate\tratio\n";
		infoString += (float)obsTotRate+"\t"+(float)rateRatio+"\t"+(float)obsTotRateAbove6pt7+"\t"+(float)rateAbove6pt7_Ratio+"\t"+(float)simMoRate+"\t"+(float)moRateRatio;			

		double[] srcGainsAtMaxRateTimeArray = null;
		double[] srcGainsAtMaxAveGainTimeArray = null;
		//Read srcGainsAtMaxRateTimeArray & srcGainsAtMaxAveGainTimeArray
		if(!isPoisson) {
			srcGainsAtMaxRateTimeArray = new double[erf.getNumFaultSystemSources()];
			srcGainsAtMaxAveGainTimeArray = new double[erf.getNumFaultSystemSources()];
			File srcGainsFile = new File(resultsDir,"srcGainsAtMaxRateAnAveGainTime.txt");
			if(srcGainsFile.exists())
				try {
					List<String> fileLines = Files.readLines(srcGainsFile, Charset.defaultCharset());
					int s=0;
					firstLine = true;
					for (String line:fileLines) {
						if(firstLine) {
							firstLine=false;
							continue;
						}
						String[] st = StringUtils.split(line,"\t");
						int sectIndex = Integer.valueOf(st[0]);
						if(s != sectIndex)
							throw new RuntimeException("bad index");
						srcGainsAtMaxRateTimeArray[s] = Double.valueOf(st[1]);
						srcGainsAtMaxAveGainTimeArray[s] = Double.valueOf(st[2]);
						s+=1;
					}
				} catch (Exception e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
		}

		
		// ******************* Plot Making **************************

		// MFDs
		ProbModelsPlottingUtils.writeMFD_ComprisonPlot(targetMFD, simMFD, plotsDir);


		// make normalized rup recurrence interval plots
		ArrayList<Double> normalizedRupRecurIntervals = new ArrayList<Double>();
		//			ArrayList<ArrayList<Double>> normalizedRupRecurIntervalsMagDepList = new ArrayList<ArrayList<Double>>();
		HashMap<Double,ArrayList<Double>> normalizedRupRecurIntervalsAperMap = new HashMap<Double,ArrayList<Double>>();

		for(double aper:aperValuesList) {
			normalizedRupRecurIntervalsAperMap.put(aper, new ArrayList<Double>());
		}
		// normalize RIs for subduction zone events
		ArrayList<Double> subductionNormRupRIs = new ArrayList<Double>();
		RupSetTectonicRegimes tectonicRegimes = fltSysRupSet.getModule(RupSetTectonicRegimes.class);
		if(tectonicRegimes==null)
			throw new RuntimeException("RupSetTectonicRegimes cannot be null");

		for(int e=0;e<normRI_ForEventList.size();e++) {
			double normRI = normRI_ForEventList.get(e);
			if(Double.isNaN(normRI))
				continue;
			normalizedRupRecurIntervals.add(normRI);
			if(tectonicRegimes.get(fltSysRupIndexForEventList.get(e)) == TectonicRegionType.SUBDUCTION_INTERFACE)
				subductionNormRupRIs.add(normRI);
			if(numAperValues>0) {
				//					double rupMag = mag_ForEventList.get(e);
				//					double aper = aperModel.getRuptureAperiodicity(rupSetIndex);
				double aper = aperModel.getRuptureAperiodicity(fltSysRupIndexForEventList.get(e));
				normalizedRupRecurIntervalsAperMap.get(aper).add(normRI);

			}
		}
		double aper=Double.NaN;
		if(numAperValues==1)
			aper=aperValuesList.get(0);	// only one value, so include for comparison
		infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedRupRecurIntervals, aper, 
				plotsDir, "Normalized Rupture RIs", "normalizedRupRecurIntervals");		

		// for subduction zones
		infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(subductionNormRupRIs, aper, 
				plotsDir, "Subduction Norm Rup RIs", "normalizedRupRecurIntervalsForSubductionZones");		
		
		// now mag-dep:
		if(numAperValues >1) {
			for(double aperVal : aperValuesList) {
				String label = "aper="+aperVal;
				infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedRupRecurIntervalsAperMap.get(aperVal), aperVal, 
						plotsDir, "Norm Rup RIs; "+label, "normRupRecurIntsForTargetAper"+aperVal);
			}
		}


		// make normalized section recurrence interval plots
		aper=Double.NaN;
		infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedSectRecurIntervals, aper, 
				plotsDir, "Normalized Section RIs", "normalizedSectRecurIntervals");		
		// now mag-dep:
		//			if(numAperValues >1) {
		//				for(int i=0;i<numAperValues;i++) {
		//					String label = magDepAperiodicity.getMagDepAperInfoString(i);
		//					infoString += ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(normalizedSectRecurIntervalsMagDepList.get(i), aperValues[i], 
		//							otherPlotsDir, "Norm Sect RIs; "+label, "normSectRecurIntsForMagRange"+i);
		//				}
		//			}
		
		// make normalized section hazard rate  plots
		ProbModelsPlottingUtils.writeNormalizedDistHazardRatePlotWithFits(normalizedSectRecurIntervals, 
				plotsDir, "Normalized Section Hazard Rate", "normalizedSectHazardRate");
		
		// paleo sites norm RIs
//		paleoPlotsDir
		File paleoPlotsDir = new File(plotsDir, "paleoSitesPlots");
		if(!paleoPlotsDir.exists()) paleoPlotsDir.mkdir();
		for(int paleoSiteID : paleoSitesNormRI_List_Map.keySet()) {
			String paleoFileName = paleoSitesNameMap.get(paleoSiteID).replace(" ", "")+"_normRI";
			ProbModelsPlottingUtils.writeNormalizedDistPlotWithFits(paleoSitesNormRI_List_Map.get(paleoSiteID), aper, 
					paleoPlotsDir, paleoSitesNameMap.get(paleoSiteID), paleoFileName);	
			paleoFileName = paleoSitesNameMap.get(paleoSiteID).replace(" ", "")+"_hazRate";
			ProbModelsPlottingUtils.writeNormalizedDistHazardRatePlotWithFits(paleoSitesNormRI_List_Map.get(paleoSiteID), 
					paleoPlotsDir, paleoSitesNameMap.get(paleoSiteID), paleoFileName);
		}

		// write infoString
		FileWriter info_fr;
		try {
			info_fr = new FileWriter(new File(resultsDir,"INFO_ForPlots.txt"));
			info_fr.write(infoString);
			info_fr.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}			

		// plot long-term rate versus time
		ProbModelsPlottingUtils.writeSimExpRateVsTime(nthRupAtEpochMap,totExpRateAtEventTimeList,
				totalLongTermRate, plotsDir);


		// plot simulated versus imposed rup rates
		for(int i=0;i<obsRupRateArray.length;i++) {
			obsRupRateArray[i] = obsRupRateArray[i]/numYears;
		}
		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermRateOfNthRups, obsRupRateArray, plotsDir, "simVsImposedRupRates", 
				"", "Imposed Rup Rate (/yr)", "Simulated Rup Rate (/yr)",Double.NaN,Double.NaN);


		// plot observed versus imposed section slip rates
		for(int i=0;i<obsSectSlipRateArray.length;i++) {
			obsSectSlipRateArray[i] = obsSectSlipRateArray[i]/numYears;
		}
		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermSlipRateForSectArray, obsSectSlipRateArray, plotsDir, "simVsImposedSectionSlipRates", 
				"", "Imposed Slip Rate (mm/yr)", "Simulated Slip Rate (mm/yr)", Double.NaN,Double.NaN);				


		// plot observed versus imposed section rates
		double[] numObsOnSectionArray = new double[obsSectRateArray.length];
		for(int i=0;i<obsSectRateArray.length;i++) {
			numObsOnSectionArray[i] = obsSectRateArray[i];
			obsSectRateArray[i] = obsSectRateArray[i]/numYears;
			obsSectRateArrayMlt7pt3[i] = obsSectRateArrayMlt7pt3[i]/numYears;
			obsSectRateArrayMgt7pt3[i] = obsSectRateArrayMgt7pt3[i]/numYears;
		}
		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(longTermPartRateForSectArray, obsSectRateArray, plotsDir, "simVsImposedSectionPartRates", 
				"", "Imposed Sect Part Rate (/yr)", "Simulated Sect Part Rate (/yr)", Double.NaN,Double.NaN);


		// write map of observed vs imposed section rates
		double[] sectPartRateRatioArray = new double[obsSectRateArray.length];
		double[] sectPartRateRatioSigmaArray = new double[obsSectRateArray.length];
		for(int i=0;i<obsSectRateArray.length;i++) {
			sectPartRateRatioArray[i] = obsSectRateArray[i]/longTermPartRateForSectArray[i];
			// Compute sigma as Poisson hi to low confidence bound ratio divided by 4.0.
			// use this for now until PoissonRateFromNinT_Calc is moved out of scratch 
			double alpha = numObsOnSectionArray[i]+1;  // also called shape parameter
			double beta = 1d/numYears; 
			GammaDistribution gd = new GammaDistribution(alpha,beta);
			sectPartRateRatioSigmaArray[i] = (gd.inverseCumulativeProbability(0.975)/gd.inverseCumulativeProbability(0.025))/4.0;;
			//				sectPartRateRatioSigmaArray[i] = 0d;
			//				double low95bound = PoissonRateFromNinT_Calc.getRateForCumulativeProb(numObsOnSectionArray[i], numYears, 0.025);
			//				double hi95bound = PoissonRateFromNinT_Calc.getRateForCumulativeProb(numObsOnSectionArray[i], numYears, 0.975);
			//				sectPartRateRatioSigmaArray[i] = (hi95bound/low95bound)4.0;
		}
		ProbModelsPlottingUtils.writeMapOfSimOverTargetPartRates (sectPartRateRatioArray, sectPartRateRatioSigmaArray, 
				fltSysRupSet.getFaultSectionDataList(), plotsDir);

		// Turn max source gains into section gains and make map if not poisson
		if(!isPoisson) {
			double[] sectGainsAtMaxRateTimeArray = new double[fltSysRupSet.getNumSections()];
			double[] sectGainsAtMaxAveGainTimeArray = new double[fltSysRupSet.getNumSections()];
			double[] sectLongTermRateArray = new double[fltSysRupSet.getNumSections()];
			for(int n=0; n<erf.getTotNumRupsFromFaultSystem();n++) {
				int[] sectID_Array = sectIndexArrayForSrcList.get(erf.getSrcIndexForNthRup(n));
				for(int s:sectID_Array) {
					// max rate time
					sectGainsAtMaxRateTimeArray[s] += longTermRateOfNthRups[n] * srcGainsAtMaxRateTimeArray[erf.getSrcIndexForNthRup(n)];
					// max ave gain time
					sectGainsAtMaxAveGainTimeArray[s] += longTermRateOfNthRups[n]*srcGainsAtMaxAveGainTimeArray[erf.getSrcIndexForNthRup(n)];
					sectLongTermRateArray[s] += longTermRateOfNthRups[n];
				}
			}
			// turn them into gains
			for(int s=0;s<sectLongTermRateArray.length;s++) {
				sectGainsAtMaxRateTimeArray[s] /= sectLongTermRateArray[s];
				sectGainsAtMaxAveGainTimeArray[s] /= sectLongTermRateArray[s];
			}
			ProbModelsPlottingUtils.writeMapOfSectionProbGains(sectGainsAtMaxRateTimeArray, fltSysRupSet.getFaultSectionDataList(), plotsDir, "mapSectGainsAtMaxRateTime");
			ProbModelsPlottingUtils.writeMapOfSectionProbGains(sectGainsAtMaxAveGainTimeArray, fltSysRupSet.getFaultSectionDataList(), plotsDir, "mapSectGainsAtMaxAveGainTime");
		}


		// write section rates with names
		FileWriter sectRates_fw;
		try {
			sectRates_fw = new FileWriter(new File(resultsDir,"/obsVsImposedSectionPartRates.txt"));
			sectRates_fw.write("sectID\timposedRate\tsimulatedRate\tsimOverImpRateRatio\tsectName\n");
			for(int i=0;i<fltSysRupSet.getNumSections();i++) {
				FaultSection fltData = fltSysRupSet.getFaultSectionData(i);
				sectRates_fw.write(fltData.getSectionId()+"\t"+longTermPartRateForSectArray[i]+"\t"+obsSectRateArray[i]+
						"\t"+sectPartRateRatioArray[i]+"\t"+fltData.getName()+"\n");
			}
			sectRates_fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		// OTHER PLOTS BELOW (in "otherPlots" dir)


		// plot ave norm RI along strike
		ProbModelsPlottingUtils.writeNormRI_SlongStrikePlot(normRI_AlongStrike, otherPlotsDir);

		//			// make aveGainVsMagFunction, weighted by rupture rate.  note that these are for all 
		//			//ruptures (not just the one that occurred at each event time)
		//			for(int i=0;i<aveRupProbGainArray.length;i++) aveRupProbGainArray[i] /= numSimulatedRups;
		//			ProbModelsPlottingUtils.writeAveGainVsMagHistPlot(magOfNthRups, aveRupProbGainArray, longTermRateOfNthRups, otherPlotsDir);
		//
		//			// plot ave, min, and max gain of each rupture vs mag (averaged over simulation regardless 
		//			//of whether the event occurred).  Not sure how useful this is. It's a big file
		//			ProbModelsPlottingUtils.writeRupGainMeanMinMaxVsMagPlot(magOfNthRups, aveRupProbGainArray, minRupProbGainArray, maxRupProbGainArray, otherPlotsDir);
		//	
		//			// write out these gains and other stuff (big file)
		//			FileWriter gain_fr;
		//			try {
		//				gain_fr = new FileWriter(new File(otherPlotsDir,"aveRupGainData.txt"));
		//				gain_fr.write("nthRupIndex\taveRupGain\tminRupGain\tmaxRupGain\trupMag\trupLongTermRate\trupCondRI\trupName\n");
		//				for(int i=0;i<aveRupProbGainArray.length;i++) {
		//					gain_fr.write(i+"\t"+aveRupProbGainArray[i]+"\t"+minRupProbGainArray[i]+"\t"+maxRupProbGainArray[i]+"\t"+magOfNthRups[i]
		//							+"\t"+longTermRateOfNthRups[i]+"\t"+aveCondRecurIntervalForFltSysRups[erf.getFltSysRupIndexForNthRup(i)]+"\t"+
		//							erf.getSource(erf.getSrcIndexForNthRup(i)).getName()+"\n");
		//				}
		//				gain_fr.close();
		//			} catch (IOException e1) {
		//				e1.printStackTrace();
		//			}	


		// this makes mag dependent section participation stuff; not sure it's still useful
		ArrayList<String> outStringList = new ArrayList<String>();
		int numSect=fltSysRupSet.getNumSections();
		double[] targetSectRateArrayMlt7pt3 = new double[numSect];
		double[] targetSectRateArrayMgt7pt3 = new double[numSect];
		for(int s=0;s<numSect;s++) {
			double partRateMlow=0;
			double partRateMhigh=0;
			for (int r : fltSysRupSet.getRupturesForSection(s)) {
				double mag = fltSysRupSet.getMagForRup(r);
				if(mag<7.3)
					partRateMlow += fltSysSolution.getRateForRup(r);
				else
					partRateMhigh = fltSysSolution.getRateForRup(r);
			}
			targetSectRateArrayMlt7pt3[s]=partRateMlow;
			targetSectRateArrayMgt7pt3[s]=partRateMhigh;
			outStringList.add(s+"\t"+obsSectRateArray[s]+"\t"+longTermPartRateForSectArray[s]+"\t"+
					(obsSectRateArray[s]/longTermPartRateForSectArray[s])+"\t"+
					targetSectRateArrayMlt7pt3[s]+"\t"+
					obsSectRateArrayMlt7pt3[s]+"\t"+
					obsSectRateArrayMlt7pt3[s]/targetSectRateArrayMlt7pt3[s]+"\t"+
					targetSectRateArrayMgt7pt3[s]+"\t"+
					obsSectRateArrayMgt7pt3[s]+"\t"+
					obsSectRateArrayMgt7pt3[s]/targetSectRateArrayMgt7pt3[s]+"\t"+
					fltSysRupSet.getFaultSectionData(s).getName()+"\n");
		}
		File dataFile2 = new File(resultsDir,"magDepSectRates_magBelowAndAbove7pt3.txt");
		try {
			FileWriter fileWriter = new FileWriter(dataFile2);
			fileWriter.write("secID\tsimRate\ttargetRate\tratio\ttargetLowMrate\tsimLowMrate\tlowRation\ttargetHighMrate\tsimHighMrate\thighRatio\tsectName\n");
			for(String line2:outStringList) {
				fileWriter.write(line2);
			}
			fileWriter.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		ProbModelsPlottingUtils.writeSimOverImposedVsImposedSecPartRates_Plot(otherPlotsDir, targetSectRateArrayMlt7pt3, 
				obsSectRateArrayMlt7pt3, 10d, numYears, "Mlt7pt3");
		ProbModelsPlottingUtils.writeSimOverImposedVsImposedSecPartRates_Plot(otherPlotsDir, targetSectRateArrayMgt7pt3, 
				obsSectRateArrayMgt7pt3, 10d, numYears, "Mgt7pt3");


		if(verbose) System.out.println("INFO STRING:\n\n"+infoString);
	}




	public static void main(String[] args) {

	}

}
