package scratch.ned.longTermTD2026;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.Range;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.earthquake.calc.recurInterval.BPT_DistCalc;
import org.opensha.sha.earthquake.calc.recurInterval.EqkProbDistCalc;
import org.opensha.sha.earthquake.calc.recurInterval.LognormalDistCalc;
import org.opensha.sha.earthquake.calc.recurInterval.WeibullDistCalc;
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
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen.DataToInclude;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.utils.ProbModelsPlottingUtils;
import scratch.ned.longTermTD2026.WeibullFit.WeibullParams;
import scratch.ned.nshm23.FSS_Fetcher2023;

public class LongTermTD_2026_Analyses {
	
	
	private static TimeDepFaultSystemSolutionERF getFullPrefUS26_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/fullPrefUS_FSS.zip";
		FaultSystemSolution sol = FSS_Fetcher2023.getPreferredFull_FSS(full_FSS_fileName);	
		
		TimeDepFaultSystemSolutionERF erf = new TimeDepFaultSystemSolutionERF();
		erf.setSolution(sol);
		
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.UCERF3_METHOD);
		
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);

		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		if (probModel instanceof UCERF3_ProbabilityModel) {
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			// setting by enum is prefferred
			u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_MIDDLE);
			// this is for simulation mode
			u3ProbModel.setSaveDebugInfo(true);
		} 
		
		ParameterList modelParams = probModel.getAdjustableParameters();
		if (modelParams.containsParameter(RenewalModels.PARAM_NAME))
			modelParams.setValue(RenewalModels.PARAM_NAME, RenewalModels.BPT);		
		
		erf.getTimeSpan().setStartTime(2025);	// this shouldn't matter
		erf.getTimeSpan().setDuration(50);; 	// this shouldn't matter
		erf.updateForecast();

		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);


		return erf;
	}
	
	
	private static void generateReportPage(String outputDir, String fss_fileNameWithPath, int startYear, int duration, 
			int histOpenIntYear, FSS_ProbabilityModels probModChoice, AperiodicityModels aperModelChoice, RenewalModels 
			renewalModelChoice, BPTAveragingTypeOptions averagingChoice, PaleoMappingAlgorithm paleoMapping, 
			DataToInclude paleoDataToInclude, String referenceDir, String titleString, String infoString,
			boolean includeUCERF3_Comp) {
		
		FaultSystemSolution sol = FSS_Fetcher2023.getPreferredFull_FSS(fss_fileNameWithPath);		
		TimeDepFaultSystemSolutionERF erf = new TimeDepFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);		
		erf.setProbabilityModelChoice(probModChoice);
		
		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		if (probModel instanceof UCERF3_ProbabilityModel) {
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			u3ProbModel.setAperiodicityModelChoice(aperModelChoice);
			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(histOpenIntYear, true));
			u3ProbModel.setSaveDebugInfo(true);
			u3ProbModel.setRenewalModelChoice(renewalModelChoice);
			u3ProbModel.setAveragingTypeChoice(averagingChoice);
			
//			// REMOVE FOLLOWING METHOD IF JAMIE'S TEST IS NO GOOD
//			double[] aveCondRecurIntervalForFltSysRupsAlt = TimeDepUtils.testJamieAveCondRecurIntervalForFltSysRups(erf.getSolution());
//			u3ProbModel.tempSetAveCondRecurIntervalForFltSysRups(aveCondRecurIntervalForFltSysRupsAlt);
		} 
		else if (probModel instanceof WG02_ProbabilityModel) {
			((WG02_ProbabilityModel)probModel).setAperiodicityModelChoice(aperModelChoice); // this throws exception if non applicable choice?
			// hard coded for BPT & no historic open interval
		}
		else if (probModel instanceof FSS_ProbabilityModel.Poisson) {
			// do nothing
		}
		else
			throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());
		erf.getTimeSpan().setStartTime(startYear);
		erf.getTimeSpan().setDuration(duration);
		// not sure this is needed here
		erf.updateForecast();
		
		TimeDepFaultSystemSolutionERF erfU3 = null;
		if(includeUCERF3_Comp) {
			erfU3 = new TimeDepFaultSystemSolutionERF();
			erfU3.setSolution(getUCERF3_BranchAveFaultSysSol());
			erfU3.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);		
			erfU3.setProbabilityModelChoice(probModChoice);
			
			FSS_ProbabilityModel probModelU3 = erfU3.getProbabilityModel();
			if (probModelU3 instanceof UCERF3_ProbabilityModel) {
				UCERF3_ProbabilityModel u3pm = (UCERF3_ProbabilityModel)probModelU3;
				u3pm.setAperiodicityModelChoice(aperModelChoice);
				u3pm.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(histOpenIntYear, true));
				u3pm.setSaveDebugInfo(true);
				u3pm.setRenewalModelChoice(renewalModelChoice);
				u3pm.setAveragingTypeChoice(averagingChoice);
			} 
			else if (probModelU3 instanceof WG02_ProbabilityModel) {
				((WG02_ProbabilityModel)probModelU3).setAperiodicityModelChoice(aperModelChoice); // this throws exception if non applicable choice?
				// hard coded for BPT & no historic open interval
			}
			else if (probModelU3 instanceof FSS_ProbabilityModel.Poisson) {
				// do nothing
			}
			else
				throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());	

			erfU3.getTimeSpan().setStartTime(startYear);
			erfU3.getTimeSpan().setDuration(duration);
			// not sure this is needed here
			erfU3.updateForecast();
		}

		try {
			TimeDependentReportPageGen.generatePage(new File(outputDir), erf, paleoMapping, 
					paleoDataToInclude, new File(referenceDir), titleString, infoString, erfU3);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// list source index involving section index
//		int targetSectID = 7898;  // Prince William sound
//		int targetSectID = 6920;  // Susitna Glacier Subsection
//		int targetSectID = 5594;  // New Madrid - SSCn (New Madrid west)
//		for(int s=0; s<erf.getNumFaultSystemSources(); s++) {
//			int fssRupID = erf.getFltSysRupIndexForSource(s);
//			List<Integer> sectIDList = erf.getSolution().getRupSet().getSectionsIndicesForRup(fssRupID);
//			if(sectIDList.contains(targetSectID)) {
//				System.out.println(s+"\tsource utilizes sect "+targetSectID+"\t"+sectIDList.size()+"\t"+erf.getSource(s).getName()+"\t"+sectIDList);
//			}
//		}
	}
	
	
	private static FaultSystemSolution getUCERF3_BranchAveFaultSysSol() {
		// get UCERF3 branch-ave solution
		File storeDir = MeanUCERF3.getStoreDir();
		File solFile = MeanUCERF3.checkDownload(new File(storeDir, "cached_FM3_1_dep100.0_depMean_rakeMean.zip")).join();
		FaultSystemSolution sol=null;
		try {
			sol = FaultSystemSolution.load(solFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sol;
	}
	
	
	private static void generateReportPageUCERF3(String outputDir, int startYear, int duration, 
			int histOpenIntYear, FSS_ProbabilityModels probModChoice, AperiodicityModels aperModelChoice, RenewalModels 
			renewalModelChoice, BPTAveragingTypeOptions averagingChoice, String titleString, String infoString) {
		
		FaultSystemSolution sol = getUCERF3_BranchAveFaultSysSol();
		
//		// get UCERF3 branch-ave solution
//		File storeDir = MeanUCERF3.getStoreDir();
//		File solFile = MeanUCERF3.checkDownload(new File(storeDir, "cached_FM3_1_dep100.0_depMean_rakeMean.zip")).join();
//		FaultSystemSolution sol=null;
//		try {
//			sol = FaultSystemSolution.load(solFile);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		TimeDepFaultSystemSolutionERF erf = new TimeDepFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);		
		erf.setProbabilityModelChoice(probModChoice);
		
		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		if (probModel instanceof UCERF3_ProbabilityModel) {
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			u3ProbModel.setAperiodicityModelChoice(aperModelChoice);
			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(histOpenIntYear, true));
			u3ProbModel.setSaveDebugInfo(true);
			u3ProbModel.setRenewalModelChoice(renewalModelChoice);
			u3ProbModel.setAveragingTypeChoice(averagingChoice);	
		} 
		else
			throw new RuntimeException("Unsupported type of FSS_ProbabilityModel: "+probModel.getName());
		
		erf.getTimeSpan().setStartTime(startYear);
		erf.getTimeSpan().setDuration(duration);

		// not sure this is needed here
		erf.updateForecast();

		try {
			TimeDependentReportPageGen.generatePageUCERF3(new File(outputDir), erf, 
					titleString, infoString);

		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	
	


	
	
	private static void generateDOLE_ReportPages() {

		long currentTimeEpoch = System.currentTimeMillis();
		String dateString = new java.text.SimpleDateFormat("MM_dd_yyyy").format(new java.util.Date (currentTimeEpoch)); // Epoch in seconds, remove '*1000' for milliseconds
		String full_FSS_fileName = "/Users/field/nshm-haz_data/fullPrefUS_FSS.zip";

		FaultSystemSolution sol = FSS_Fetcher2023.getPreferredFull_FSS(full_FSS_fileName);		
		
		File tdMainDir = new File("/Users/field/markdown/nshm23_time_dependence_"+dateString);

		try {
			
//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
//			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "onlyHistoricRupDOLE"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.HIST_RUPS_ONLY);

////			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "allDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.ALL_DATA);
//
////			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
//			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "allDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.ALL_DATA);
//
////			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
//			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "forDebugging_onlyPaleoDOLE_nearestSubsect"), sol, PaleoMappingAlgorithm.CLOSEST_SECT, DataToInclude.PALEO_ONLY);
//
////			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
//			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "onlyPaleoDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.PALEO_ONLY);
//
////			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
//			TimeDependentReportPageGen.generateOldPage(new File(tdMainDir, "onlyPaleoDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.PALEO_ONLY);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}
	
	public static void generatePreliminaryResults() {
		
		String rootDir = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/PreliminaryResults/";
		String referenceDir = rootDir+"ReferenceModel/";
		String fss_fileNameWithPath = rootDir+"fullPrefUS_FSS.zip";
		String outputDir, titleString, infoString;


		int startYear = 2026;
		int duration = 30; 
		int histOpenIntYear = 1875;
		FSS_ProbabilityModels probModChoice = FSS_ProbabilityModels.UCERF3_METHOD;
		AperiodicityModels aperModelChoice = AperiodicityModels.NSHM26_MIDDLE;
		RenewalModels renewalModelChoice = RenewalModels.BPT;
		BPTAveragingTypeOptions averagingChoice = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
		PaleoMappingAlgorithm paleoMapping = PaleoMappingAlgorithm.NEIGHBORING_SECTS;
		DataToInclude paleoDataToInclude = DataToInclude.HIST_RUPS_ONLY;
		boolean includeUCERF3_Comp = false;
		
		// Reference calculation
		outputDir = referenceDir;
		titleString = "Reference 2026 Time-Dependent Model";
		infoString = "This is the Reference (preferred) model. ";
		includeUCERF3_Comp=true;
		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
				 histOpenIntYear,  probModChoice,  aperModelChoice,  
				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
				referenceDir, titleString, infoString, includeUCERF3_Comp);
	
		
		// UCERF3 results
//		outputDir = rootDir+"ReferenceUCERF3/";
//		titleString = "Reference UCERF3";
//		infoString = "Reference parameters applied to UCERF3 branch-average fault system solution, and using UCERF3 DOLE data. ";
//		generateReportPageUCERF3(outputDir, startYear, duration, 
//				histOpenIntYear, probModChoice, aperModelChoice,  
//				renewalModelChoice, averagingChoice, titleString, infoString, includeUCERF3_Comp);	
		
//		// BPTAveragingTypeOptions.AVE_RI_AVE_TIME_SINCE option
//		outputDir = rootDir+"ReferencelModel_AVE_RI_AVE_TIME_SINCE/";
//		titleString = "Reference, But With BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE";
//		infoString = "This is the Reference (preferred) model, but with AVE_RI_AVE_NORM_TIME_SINCE chosen for BPTAveragingTypeOptions. ";
//		averagingChoice = BPTAveragingTypeOptions.AVE_RI_AVE_TIME_SINCE;
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);

	
		
//		outputDir = rootDir+"ReferencelModel_LowCOV/";
//		titleString = "Reference, But With Low COV/Aperiodicity Branch";
//		infoString = "This is the Reference (preferred) model, but with the low COV/aperiodicity branch. ";
//		aperModelChoice = AperiodicityModels.NSHM26_LOW;
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);

		
//		outputDir = rootDir+"ReferencelModel_HighCOV/";
//		titleString = "Reference, But With High COV/Aperiodicity Branch";
//		infoString = "This is the Reference (preferred) model, but with the high COV/aperiodicity branch. ";
//		aperModelChoice = AperiodicityModels.NSHM26_HIGH;
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);
		
		
//		outputDir = rootDir+"ReferencelModel_PlusPaleoDOLE/";
//		titleString = "Reference, But With Paleo DOLE Included";
//		infoString = "This is the Reference (preferred) model, but with the paleoseismic date of last event (DOLE) data included. ";
//		paleoDataToInclude = DataToInclude.ALL_DATA;
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);

//		outputDir = rootDir+"ReferencelModel_OnlyPaleoDOLE/";
//		titleString = "Reference, But With Only Paleo DOLE";
//		infoString = "This is the Reference (preferred) model, but with only paleoseismic date of last event (DOLE) data included. ";
//		paleoDataToInclude = DataToInclude.PALEO_ONLY;
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);

		
//		outputDir = rootDir+"LognormalModel/";
//		titleString = "Reference, But With Lognormal Renewal Model";
//		renewalModelChoice = RenewalModels.LOGNORMAL;
//		infoString = "This is the Reference model, but with the renewal model switched to Lognormal. ";
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);


//		outputDir = rootDir+"NoOpenIntervalModel/";
//		histOpenIntYear=startYear;
//		titleString = "Reference, But With No Open Interval";
//		infoString = "This is the Reference model, but where the open interval year is "+histOpenIntYear+".";
//		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);
		

//		outputDir = rootDir+"WeibullModel/";
//		titleString = "Reference, But With Weibull Renewal Model";
//		renewalModelChoice = RenewalModels.WEIBULL;
//		infoString = "This is the Reference model, but with the renewal model switched to Weibull. ";
//		generateReportPage(outputDir, fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);
		

//		outputDir = rootDir+"BPT_SlipRateDepAperiodicity/";
//		aperModelChoice = AperiodicityModels.NSHM26_SLIPRATE_TEST;
//		titleString = "Reference, But With Slip-rate Dependent Aperiodicity";
//		infoString = "This is the Reference model, but with slip-rate dependent aperiodicity (COV). ";
//		generateReportPage(outputDir, fss_fileNameWithPath,  startYear,  duration, 
//				 histOpenIntYear,  probModChoice,  aperModelChoice,  
//				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
//				referenceDir, titleString, infoString);


	}

	
	private static void setDOLE_asFractionOfRI(TimeDepFaultSystemSolutionERF erf, double fractRI) {
		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		double[] longTermPartRateForSectArray = probModel.getSectLongTermPartRates(); // this is a duplicate
		FaultSystemRupSet fltSysRupSet = erf.getSolution().getRupSet();
		long origStartTimeMillis = erf.getTimeSpan().getStartTimeInMillis();
		for(int s=0; s<fltSysRupSet.getNumSections();s++) {
			long doleMillis = origStartTimeMillis-(long)(fractRI*(1.0/longTermPartRateForSectArray[s])*TimeDepUtils.MILLISEC_PER_YEAR);
			fltSysRupSet.getFaultSectionData(s).setDateOfLastEvent(doleMillis);
			probModel.setSectDOLE(s, doleMillis);
		}
	}
	
	
	public static void generateRenewalModelPlots(boolean extrapolate) {
		
		String rootDir = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/RenewalModelTestPlots";
		String dirName = "Precision"+EqkProbDistCalc.getPrecision();  // change the value in this class if desired.
		if(extrapolate)
			dirName += "_extrapolated";
		File parentDir = new File(rootDir,dirName);
		if(!parentDir.exists()) 
			parentDir.mkdir();
//		File outputDir = new File(parentDir,"BlaBla"); 

		double mean=1.0; 
		double deltaX = 5e-5; 
		int numPoints= 200001;
		
		EqkProbDistCalc[] renewalModelArray = {new BPT_DistCalc(), new LognormalDistCalc(), new WeibullDistCalc()};
		double[] aperArray = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}; 
		double[] durationArray = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.3};

//		EqkProbDistCalc[] renewalModelArray = {new WeibullDistCalc()};
//		double[] aperArray = {0.8}; 
//		double[] durationArray = {1e-2};

		for(EqkProbDistCalc renewalModel:renewalModelArray) {
			
			ArrayList<XY_DataSet> pdfFuncs = new  ArrayList<XY_DataSet>();
			ArrayList<XY_DataSet> surviveFuncs = new  ArrayList<XY_DataSet>();
			ArrayList<XY_DataSet> hazardFuncs = new  ArrayList<XY_DataSet>();
			ArrayList<ArrayList<XY_DataSet>> condProbFuncsList = new  ArrayList<ArrayList<XY_DataSet>>();
			ArrayList<ArrayList<XY_DataSet>> probVsOpenIntervalFuncsList = new  ArrayList<ArrayList<XY_DataSet>>();
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			String xAxisLabel = "";
			String yAxisLabel = "";
			Range xAxisRange = new Range(1e-2,10);
			Range yAxisRange = null;
			boolean logX = false;
			boolean logY = false;
			double widthInches= 5;
			double heightInches = 4;
			File fileNamePrefix = null; 
			boolean popupWindow = false;

//			String plotName = "PDF";
			for(double aper:aperArray) {
				renewalModel.setAllParameters(mean, aper, deltaX, numPoints);
				pdfFuncs.add(renewalModel.getPDF());
				surviveFuncs.add(renewalModel.getSurvivorFunc());
				hazardFuncs.add(renewalModel.getHazFunc());
			}
			
			for(double dur:durationArray) {
				ArrayList<XY_DataSet> condProbFuncs = new  ArrayList<XY_DataSet>();
				ArrayList<XY_DataSet> probVsOpenIntervalFuncs = new  ArrayList<XY_DataSet>();
				for(double aper:aperArray) {
					renewalModel.setAllParameters(mean, aper, deltaX, numPoints);
					condProbFuncs.add(renewalModel.getCondProbGainFunc(dur,extrapolate));
					probVsOpenIntervalFuncs.add(renewalModel.getCondProbForUnknownTimeSinceLastEventFunc(dur, 10));
				}
				condProbFuncsList.add(condProbFuncs);
				probVsOpenIntervalFuncsList.add(probVsOpenIntervalFuncs);
			}
			

			String fileName = renewalModel.getName()+"_PDF";
			ProbModelsPlottingUtils.writeAndOrPlotFuncs(
					pdfFuncs, plotChars, "PDF - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
					yAxisRange, true, true, widthInches, heightInches, new File(parentDir,fileName), popupWindow);

			fileName = renewalModel.getName()+"_Survival";
			ProbModelsPlottingUtils.writeAndOrPlotFuncs(
					surviveFuncs, plotChars, "Survival - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
					yAxisRange, true, true, widthInches, heightInches, new File(parentDir,fileName), popupWindow);

			fileName = renewalModel.getName()+"_Hazard";
			ProbModelsPlottingUtils.writeAndOrPlotFuncs(
					hazardFuncs, plotChars, "Hazard - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
					yAxisRange, logX, logY, widthInches, heightInches, new File(parentDir,fileName), popupWindow);
			fileName = renewalModel.getName()+"_Hazard_Log";
			ProbModelsPlottingUtils.writeAndOrPlotFuncs(
					hazardFuncs, plotChars, "Hazard - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
					yAxisRange, true, true, widthInches, heightInches, new File(parentDir,fileName), popupWindow);

			for(int i=0;i< durationArray.length;i++) {
				fileName = renewalModel.getName()+"_CondProbGain_"+durationArray[i];
				ProbModelsPlottingUtils.writeAndOrPlotFuncs(
						condProbFuncsList.get(i), plotChars, "CondProbGain - "+durationArray[i]+" - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
						yAxisRange, logX, logY, widthInches, heightInches, new File(parentDir,fileName), popupWindow);				
		
				fileName = renewalModel.getName()+"_CondProbGain_"+durationArray[i]+"_Log";
				ProbModelsPlottingUtils.writeAndOrPlotFuncs(
						condProbFuncsList.get(i), plotChars, "CondProbGain - "+durationArray[i]+" - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
						yAxisRange, true, true, widthInches, heightInches, new File(parentDir,fileName), popupWindow);				

//				fileName = renewalModel.getName()+"_ProbVsOpenInt_"+durationArray[i];
//				ProbModelsPlottingUtils.writeAndOrPlotFuncs(
//						probVsOpenIntervalFuncsList.get(i), plotChars, "ProbVsOpenInt - "+durationArray[i]+" - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
//						yAxisRange, logX, logY, widthInches, heightInches, new File(parentDir,fileName), popupWindow);				
		
				fileName = renewalModel.getName()+"_ProbVsOpenInt_"+durationArray[i]+"_Log";
				ProbModelsPlottingUtils.writeAndOrPlotFuncs(
						probVsOpenIntervalFuncsList.get(i), plotChars, "ProbVsOpenInt - "+durationArray[i]+" - "+renewalModel.getName(), xAxisLabel, yAxisLabel, xAxisRange,
						yAxisRange, true, true, widthInches, heightInches, new File(parentDir,fileName), popupWindow);				

			}

		}

	}
	
	
	public static void testU3_NSHM23_CA_Fault_Mappings() {
		
		
		// this is to get the "Linkto2014" IDs from Alex's geojson files
		HashMap<Integer,Integer> u3ID_for_nshm23iD = new HashMap<Integer,Integer>();
		// these are parent sections, so their ID is the subsection's parent ID
		List<? extends FaultSection> sects=null;
		try {
			sects = NSHM23_FaultModels.WUS_FM_v3.getFaultSections();
		} catch (IOException e) {
			e.printStackTrace();
		}
		for (FaultSection sect : sects) {
			int id2014 = ((GeoJSONFaultSection)sect).getProperties().getInt("Linkto2014", -1);
			if (id2014 >= 0) {
				int nshm23_ID = sect.getSectionId(); // because these are already parent IDs
				u3ID_for_nshm23iD.put(nshm23_ID, id2014);
	if(sect.getSectionName().contains("San Andreas"))
		System.out.println(nshm23_ID+"\t"+id2014+"\t"+sect.getSectionName());
			}
		}

//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		FaultSystemSolution fss = erf.getSolution();
//		FaultSystemRupSet rupSet = fss.getRupSet();
		
		FaultSystemSolution wusSol=null;
		try {
			wusSol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println();
		HashMap<Integer,String> id_name_map = new HashMap<Integer,String>();
		for (FaultSection fs: wusSol.getRupSet().getFaultSectionDataList()) {
			if(!id_name_map.keySet().contains(fs.getParentSectionId())){
				id_name_map.put(fs.getParentSectionId(),fs.getParentSectionName());
			}
		}
		for(int id:id_name_map.keySet())
			if(id_name_map.get(id).contains("San Andreas"))
				System.out.println(id+"\t"+id_name_map.get(id));
		
		
		// see main in TD_ERF_Example below to see how to get the full ERF in the new framework
		FaultSystemSolution sol=null;
		try {
			sol = scratch.kevin.tdProbModelPlayground.TD_ERF_Example.fetchU3_BA();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// OLD WAY
//		String fileName="/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";
//		FaultSystemSolution sol=null;
//		try {
//			sol = FaultSystemSolution.load(new File(fileName));
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		FaultSystemSolutionERF erfU3 = new FaultSystemSolutionERF(fileName);
//		erfU3.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
//		erfU3.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
//		erfU3.getParameter(MagDependentAperiodicityParam.NAME).setValue(MagDependentAperiodicityOptions.MID_VALUES);
//		BPTAveragingTypeOptions aveType = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
//		erfU3.setParameter(BPTAveragingTypeParam.NAME, aveType);
//		// erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
//		erfU3.updateForecast();	

		System.out.println();
		HashMap<Integer,String> id_name_mapU3 = new HashMap<Integer,String>();
		for (FaultSection fs: sol.getRupSet().getFaultSectionDataList()) {
			if(!id_name_mapU3.keySet().contains(fs.getParentSectionId())){
				id_name_mapU3.put(fs.getParentSectionId(),fs.getParentSectionName());
			}
		}
		for(int id:id_name_mapU3.keySet())
			if(id_name_mapU3.get(id).contains("San Andreas"))
				System.out.println(id+"\t"+id_name_mapU3.get(id));

		
		String goodString = "";
		String badString = "";
		int numGood=0;
		int numBad=0;
		
		for(int fltID:u3ID_for_nshm23iD.keySet()) {
			int u3Id = u3ID_for_nshm23iD.get(fltID);
				if(id_name_map.get(fltID).equals(id_name_mapU3.get(u3Id))) {  // names the same
					goodString += fltID+"\t"+id_name_map.get(fltID)+"\n";
					numGood+=1;
//					System.out.println(fltID+"\t"+id_name_map.get(fltID));
				}
				else {
					badString += fltID+"\t"+id_name_map.get(fltID)+"\tDIFF:  "+id_name_mapU3.get(u3Id)+"\n";
					numBad+=1;
//					System.out.println(fltID+"\t"+id_name_map.get(fltID)+"\tDIFF:  "+id_name_mapU3.get(fltID));
				}
		}
//		System.out.println("\n\ngoodString (size="+numGood+"):\n"+goodString);
//		System.out.println("\n\nbadString (size="+numBad+"):\n"+badString);

		
//		for(int fltID:id_name_map.keySet()) {
//			if(id_name_mapU3.keySet().contains(fltID)) { // same ID
//				if(id_name_map.get(fltID).equals(id_name_mapU3.get(fltID))) {  // names the same
//					goodString += fltID+"\t"+id_name_map.get(fltID)+"\n";
//					numGood+=1;
////					System.out.println(fltID+"\t"+id_name_map.get(fltID));
//				}
//				else {
//					badString += fltID+"\t"+id_name_map.get(fltID)+"\tDIFF:  "+id_name_mapU3.get(fltID)+"\n";
//					numBad+=1;
////					System.out.println(fltID+"\t"+id_name_map.get(fltID)+"\tDIFF:  "+id_name_mapU3.get(fltID));
//				}
//			}
//		}
//		System.out.println("\n\ngoodString (size="+numGood+"):\n"+goodString);
//		System.out.println("\n\nbadString (size="+numBad+"):\n"+badString);

//		for(int fltID:id_name_map.keySet()) {
//			String name = id_name_map.get(fltID);
//			if(id_name_mapU3.containsValue(name)) {
//				for(int idU3 : id_name_mapU3.keySet())
//					if(id_name_mapU3.get(idU3).equals(name))
//						System.out.println(id_name_map.get(fltID)+"\t"+fltID+"\t"+idU3);
//			}	
//		}

	}
	
	/**
	 * U3 Parent ID is the key and associated NSHM23 Parent ID is the value.  Mapping is found 
	 * by brute-force comparison of subsection end points.  Ambiguities arose when only a horizontal
	 * distance comparison was used.  This can also be modified to provide subsection mappings.
	 * @return HashMap<Integer,Integer>
	 */
	public static HashMap<Integer,Integer> getNSHM23_For_U3_ParentFaultID_Mapping() {
		
		HashMap<Integer,Integer> wusFromU3_ID_Map = new HashMap<Integer,Integer>();

		FaultSystemSolution wusSol=null;
		try {
			wusSol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/results_WUS_FM_v3_branch_averaged_gridded_simplified.zip"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		FaultSystemSolution u3_sol=null;
		try {
			u3_sol = scratch.kevin.tdProbModelPlayground.TD_ERF_Example.fetchU3_BA();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		List<? extends FaultSection> wusFltSectList = wusSol.getRupSet().getFaultSectionDataList();
		List<? extends FaultSection> u3_FltSectList = u3_sol.getRupSet().getFaultSectionDataList();
		String errorLog = "";
		for(int i=0;i<u3_FltSectList.size();i++) {
			FaultSection u3_Sect = u3_FltSectList.get(i);
			int u3_ID = u3_Sect.getParentSectionId();
			String u3_Name = u3_Sect.getParentSectionName();
			for(int j=0;j<wusFltSectList.size();j++) {
				FaultSection wusSect = wusFltSectList.get(j);
				int wusID = wusSect.getParentSectionId();
				String wusName = wusSect.getParentSectionName();
				
				Location firstLoc1 = u3_Sect.getFaultTrace().getFirst();
				Location lastLoc1 = u3_Sect.getFaultTrace().getLast();
				Location firstLoc2 = wusSect.getFaultTrace().getFirst();
				Location lastLoc2 = wusSect.getFaultTrace().getLast();
				if(LocationUtils.linearDistance(firstLoc1, firstLoc2)<1.0) {
					if(LocationUtils.linearDistance(lastLoc1, lastLoc2)<1.0) {
//						System.out.println(u3_ID+"\t"+wusID+"\t"+u3_Name+"\t"+wusName);
						if(wusFromU3_ID_Map.keySet().contains(u3_ID)) {
							if(wusFromU3_ID_Map.get(u3_ID) != wusID) {
								throw new RuntimeException(u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID);
//								errorLog+=u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID+"\t"+u3_Name+"\t"+wusName+"\n";
							}
						}
						else {
							wusFromU3_ID_Map.put(u3_ID, wusID);
						}
					}
				}
				else if(LocationUtils.linearDistance(firstLoc1, lastLoc2)<1.0) {
					if(LocationUtils.linearDistance(lastLoc1, firstLoc2)<1.0) {
//						System.out.println(u3_ID+"\t"+wusID+"\t"+u3_Name+"\t"+wusName);
						if(wusFromU3_ID_Map.keySet().contains(u3_ID)) {
							if(wusFromU3_ID_Map.get(u3_ID) != wusID) {
									throw new RuntimeException(u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID);
//								errorLog+=u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID+"\t"+u3_Name+"\t"+wusName+"\n";
							}
						}
						else {
							wusFromU3_ID_Map.put(u3_ID, wusID);
						}
					}		
				}
			}
		}
//		for(int key:wusFromU3_ID_Map.keySet())
//			System.out.println(key+"\t"+wusFromU3_ID_Map.get(key));
//		System.out.println("Num matches = "+wusFromU3_ID_Map.size());
//		System.err.println(errorLog);
			
		return wusFromU3_ID_Map;
		
	}
	
	
	
	public static HashMap<Integer,Integer> old_getNSHM23_For_U3_SubsectionID_Mapping() {
		
		HashMap<Integer,Integer> wusFromU3_ID_Map = new HashMap<Integer,Integer>();

		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
		
		FaultSystemSolution u3_sol=null;
		try {
			u3_sol = scratch.kevin.tdProbModelPlayground.TD_ERF_Example.fetchU3_BA();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return TimeDependentReportPageGen.getSubsectionID_Mapping(u3_sol,erf.getSolution());
		
//		List<? extends FaultSection> wusFltSectList = wusSol.getRupSet().getFaultSectionDataList();
//		List<? extends FaultSection> u3_FltSectList = u3_sol.getRupSet().getFaultSectionDataList();
//		String errorLog = "";
//		for(int i=0;i<u3_FltSectList.size();i++) {
//			FaultSection u3_Sect = u3_FltSectList.get(i);
//			int u3_ID = u3_Sect.getSectionId();
//			String u3_Name = u3_Sect.getSectionName();
//			for(int j=0;j<wusFltSectList.size();j++) {
//				FaultSection wusSect = wusFltSectList.get(j);
//				int wusID = wusSect.getSectionId();
//				String wusName = wusSect.getSectionName();
//				
//				Location firstLoc1 = u3_Sect.getFaultTrace().getFirst();
//				Location lastLoc1 = u3_Sect.getFaultTrace().getLast();
//				Location firstLoc2 = wusSect.getFaultTrace().getFirst();
//				Location lastLoc2 = wusSect.getFaultTrace().getLast();
//				if(LocationUtils.linearDistance(firstLoc1, firstLoc2)<1.0) {
//					if(LocationUtils.linearDistance(lastLoc1, lastLoc2)<1.0) {
////						System.out.println(u3_ID+"\t"+wusID+"\t"+u3_Name+"\t"+wusName);
//						if(wusFromU3_ID_Map.keySet().contains(u3_ID)) {
//							if(wusFromU3_ID_Map.get(u3_ID) != wusID) {
////								throw new RuntimeException(u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID);
//								errorLog+=u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID+"\t"+u3_Name+"\t"+wusName+"\n";
//							}
//						}
//						else {
//							wusFromU3_ID_Map.put(u3_ID, wusID);
//						}
//					}
//				}
//				else if(LocationUtils.linearDistance(firstLoc1, lastLoc2)<1.0) {
//					if(LocationUtils.linearDistance(lastLoc1, firstLoc2)<1.0) {
////						System.out.println(u3_ID+"\t"+wusID+"\t"+u3_Name+"\t"+wusName);
//						if(wusFromU3_ID_Map.keySet().contains(u3_ID)) {
//							if(wusFromU3_ID_Map.get(u3_ID) != wusID) {
////									throw new RuntimeException(u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID);
//								errorLog+=u3_ID+" previously associated with "+wusFromU3_ID_Map.get(u3_ID)+" and now with "+wusID+"\t"+u3_Name+"\t"+wusName+"\n";
//							}
//						}
//						else {
//							wusFromU3_ID_Map.put(u3_ID, wusID);
//						}
//					}		
//				}
//			}
//		}
////		for(int key:wusFromU3_ID_Map.keySet())
////			System.out.println(key+"\t"+wusFromU3_ID_Map.get(key));
//		System.out.println("Num matches = "+wusFromU3_ID_Map.size()+" out of "+u3_sol.getRupSet().getNumSections()+" in UCERF3");
//		System.err.println(errorLog);
//			
//		return wusFromU3_ID_Map;
		
	}

	
	public static void weibullSamplingOskinTests() {
		
		// Test the sampling
		WeibullDistCalc wCalc = new WeibullDistCalc();
		wCalc.setAll(1.0, 0.5, 0.01, 1000);
		double[] samples = wCalc.getRandomSamples(1000000);
		HistogramFunction hFunc = new HistogramFunction(wCalc.getPDF().getMinX(), wCalc.getPDF().size(), wCalc.getPDF().getDelta());
		for(double val:samples)
			hFunc.add(val, 1.0);
		hFunc.normalizeToPDF();
		
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
		funcList.add(hFunc);
		funcList.add(wCalc.getPDF());
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String xAxisLabel = "RI (yrs)";
		String yAxisLabel = "Density";
		Range xAxisRange = new Range(0, 5);
		Range yAxisRange = new Range(0, 5);
		boolean logX = false;
		boolean logY = false;
		double widthInches = 7; // inches
		double heightInches = 6; // inches
		boolean popupWindow = true;
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcList,plotChars,"Test",xAxisLabel,yAxisLabel,xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, null, popupWindow);

		// the following is a maximum likelihood fit
		WeibullParams fitParams = WeibullFit.fit(samples);
        System.out.println("COV: " + wCalc.getCOV_fromShapeParam(fitParams.shape));
        System.out.println("Mean: " + wCalc.getMean_fromShapeAndScaleParams(fitParams.shape, fitParams.scale));

        
		HistogramFunction meanDist = new HistogramFunction(0.005, 1000, 0.01);
		HistogramFunction covDist = new HistogramFunction(0.005, 150, 0.01);
		HistogramFunction meanDistFromMoments = new HistogramFunction(0.005, 1000, 0.01);
		HistogramFunction covDistFromMoments = new HistogramFunction(0.005, 150, 0.01);
		int numSampleserWindow = 10;
		int numWindows = 10000;
		for(int i=0;i<numWindows;i++) {
			samples = wCalc.getRandomSamples(numSampleserWindow);
			fitParams = WeibullFit.fit(wCalc.getRandomSamples(numSampleserWindow));
			double cov = wCalc.getCOV_fromShapeParam(fitParams.shape);
			double mean = wCalc.getMean_fromShapeAndScaleParams(fitParams.shape, fitParams.scale);
			covDist.add(cov, 1.0);
			meanDist.add(mean, 1.0);
			
			DescriptiveStatistics stats = new DescriptiveStatistics(samples);
			double mean2 = stats.getMean();
			meanDistFromMoments.add(mean2, 1.0);
			covDistFromMoments.add(stats.getStandardDeviation()/mean2, 1.0);
		}
		covDist.normalizeToPDF();
		covDist.setInfo("COV mean = "+covDist.computeMean()+"\nCOV StdDev = "+covDist.computeStdDev());
		meanDist.normalizeToPDF();
		meanDist.setInfo("Mean mean = "+meanDist.computeMean()+"\nMean StdDev = "+meanDist.computeStdDev());
	
		covDistFromMoments.normalizeToPDF();
		covDistFromMoments.setInfo("COV mean (from moments) = "+covDistFromMoments.computeMean()+"\nCOV StdDev = "+covDistFromMoments.computeStdDev());
		meanDistFromMoments.normalizeToPDF();
		meanDistFromMoments.setInfo("Mean mean (from moments) = "+meanDistFromMoments.computeMean()+"\nMean StdDev = "+meanDistFromMoments.computeStdDev());

		ArrayList<XY_DataSet> funcListMean = new ArrayList<XY_DataSet>();
		funcListMean.add(meanDist);
		xAxisLabel = "Mean";
		yAxisLabel = "Density";
		xAxisRange = new Range(0, 5);
		yAxisRange = new Range(0, 5);
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcListMean,plotChars,"Mean Distribution from Max Likihood (N="+numSampleserWindow+")",xAxisLabel,yAxisLabel,xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, null, popupWindow);

		ArrayList<XY_DataSet> funcListMeanFromMoments = new ArrayList<XY_DataSet>();
		funcListMeanFromMoments.add(meanDistFromMoments);
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcListMeanFromMoments,plotChars,"Mean Distribution from Moments (N="+numSampleserWindow+")",xAxisLabel,yAxisLabel,xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, null, popupWindow);
		
		HistogramFunction cdf = covDist.getCumulativeDistFunctionWithHalfBinOffset();
		String str = ""+
				"\nProb(0.05-0.15) = "+(float)(cdf.getInterpolatedY(0.15)-cdf.getInterpolatedY(0.05))+
				"\nProb(0.15-0.25) = "+(float)(cdf.getInterpolatedY(0.25)-cdf.getInterpolatedY(0.15))+
				"\nProb(0.25-0.35) = "+(float)(cdf.getInterpolatedY(0.35)-cdf.getInterpolatedY(0.25))+
				"\nProb(0.35-0.45) = "+(float)(cdf.getInterpolatedY(0.45)-cdf.getInterpolatedY(0.35))+
				"\nProb(0.45-0.55) = "+(float)(cdf.getInterpolatedY(0.55)-cdf.getInterpolatedY(0.45))+
				"\nProb(0.55-0.65) = "+(float)(cdf.getInterpolatedY(0.65)-cdf.getInterpolatedY(0.55))+
				"\nProb(0.65-0.75) = "+(float)(cdf.getInterpolatedY(0.75)-cdf.getInterpolatedY(0.65))+
				"\nProb(0.75-0.85) = "+(float)(cdf.getInterpolatedY(0.85)-cdf.getInterpolatedY(0.75))+
				"\nProb(0.85-0.95) = "+(float)(cdf.getInterpolatedY(0.95)-cdf.getInterpolatedY(0.85))+
				"\nProb(0.95-1.05) = "+(float)(cdf.getInterpolatedY(1.05)-cdf.getInterpolatedY(0.95));
		covDist.setInfo(covDist.getInfo()+"\n"+str);
		ArrayList<XY_DataSet> funcListCOV = new ArrayList<XY_DataSet>();
		funcListCOV.add(covDist);
		xAxisLabel = "COV";
		yAxisLabel = "Density";
		xAxisRange = new Range(0, 3);
		yAxisRange = new Range(0, 5);
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcListCOV,plotChars,"COV Distribution from Max Likihood (N="+numSampleserWindow+")",xAxisLabel,yAxisLabel,xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, null, popupWindow);

		cdf = covDistFromMoments.getCumulativeDistFunctionWithHalfBinOffset();
		str = ""+
				"\nProb(0.05-0.15) = "+(float)(cdf.getInterpolatedY(0.15)-cdf.getInterpolatedY(0.05))+
				"\nProb(0.15-0.25) = "+(float)(cdf.getInterpolatedY(0.25)-cdf.getInterpolatedY(0.15))+
				"\nProb(0.25-0.35) = "+(float)(cdf.getInterpolatedY(0.35)-cdf.getInterpolatedY(0.25))+
				"\nProb(0.35-0.45) = "+(float)(cdf.getInterpolatedY(0.45)-cdf.getInterpolatedY(0.35))+
				"\nProb(0.45-0.55) = "+(float)(cdf.getInterpolatedY(0.55)-cdf.getInterpolatedY(0.45))+
				"\nProb(0.55-0.65) = "+(float)(cdf.getInterpolatedY(0.65)-cdf.getInterpolatedY(0.55))+
				"\nProb(0.65-0.75) = "+(float)(cdf.getInterpolatedY(0.75)-cdf.getInterpolatedY(0.65))+
				"\nProb(0.75-0.85) = "+(float)(cdf.getInterpolatedY(0.85)-cdf.getInterpolatedY(0.75))+
				"\nProb(0.85-0.95) = "+(float)(cdf.getInterpolatedY(0.95)-cdf.getInterpolatedY(0.85))+
				"\nProb(0.95-1.05) = "+(float)(cdf.getInterpolatedY(1.05)-cdf.getInterpolatedY(0.95));
		covDistFromMoments.setInfo(covDistFromMoments.getInfo()+"\n"+str);
		ArrayList<XY_DataSet> funcListCOVFromMoments = new ArrayList<XY_DataSet>();
		funcListCOVFromMoments.add(covDistFromMoments);
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcListCOVFromMoments,plotChars,"COV Distribution from Moments (N="+numSampleserWindow+")",xAxisLabel,yAxisLabel,xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, null, popupWindow);
	
	}
	
	/**
	 * This provides the id and name of all parent sections that contain a given string
	 * @param string
	 * @param erf
	 */
	public static void listParentSectionsThatContainStringInName(String string, TimeDepFaultSystemSolutionERF erf) {
		ArrayList<Integer> alreadyFound = new ArrayList<Integer>();
		System.out.println("Parent sections that have "+string+" in the name:");
		for(FaultSection fs:erf.getSolution().getRupSet().getFaultSectionDataList()) {
			if(fs.getParentSectionName() != null && fs.getParentSectionName().contains(string)) {
				int id=fs.getParentSectionId();
				if(!alreadyFound.contains(id)) {
					System.out.println(id+"\t"+fs.getParentSectionName());
					alreadyFound.add(id);
				}
			}	
		}
	}

	

	public static void main(String[] args) {
		
//		weibullSamplingOskinTests();
		
//		getNSHM23_For_U3_ParentFaultID_Mapping();
//		old_getNSHM23_For_U3_SubsectionID_Mapping();
//		System.exit(0);
	
//		// temp junk
//		TimeDepFaultSystemSolutionERF erf2 = getFullPrefUS26_ERF();
//		HashMap<Integer,String> nameFromID_Map = new HashMap<Integer,String>();
//		for (FaultSection sect:erf2.getSolution().getRupSet().getFaultSectionDataList()) {
//			if(!nameFromID_Map.keySet().contains(sect.getParentSectionId())) {
//				String name = sect.getParentSectionName();
//				nameFromID_Map.put(sect.getParentSectionId(), name);
////				System.out.println(sect.getParentSectionId()+"\t"+name);
//				if(name == null) {
//					System.out.println(sect.getParentSectionId()+"\t"+sect.getSectionId()+"\t"+sect.getSectionName());
//				}
//			}
//		}
//		System.exit(0);
		
//		testU3_NSHM23_CA_Fault_Mappings(); System.exit(0);		

//		listParentSectionsThatContainStringInName("Peninsula",getFullPrefUS26_ERF());
		
		generatePreliminaryResults();
		System.exit(0);
		
		
//		// write min and max ave slip rate for each rupture
//		String full_FSS_fileName = "/Users/field/nshm-haz_data/fullPrefUS_FSS.zip";
//		FaultSystemSolution sol = FSS_Fetcher2023.getPreferredFull_FSS(full_FSS_fileName);	
//		double min=1000, max=0;
//		FaultSystemRupSet rupSet = sol.getRupSet();
//		for(int i=0;i<rupSet.getNumRuptures();i++) {
//			double sr = rupSet.getAveSlipRateForRup(i);
//			if(min>sr) min=sr;
//			if(max<sr) max=sr;
//		}
//		System.out.println("MIN RUP SLIP_RATE = "+min*1000);
//		System.out.println("MAX RUP SLIP_RATE = "+max*1000);



//		// RENEWAL MODEL TEST PLOTS
//		// Hard-code-adjust the parameter EqkProbDistCalc.NUMERICAL_PRECISION 
//		// to generate different directory/results to see the effect on curves
//		// (this was done to determine the default value)
//		generateRenewalModelPlots(true);
//		generateRenewalModelPlots(false);
		
//		// This tests expm1 and log1p
//		for(double i=-16;i<=0;i++) {
//			double expt = Math.pow(10, i);
//			double prob = 1-Math.exp(-expt);
//			double expt2 = -Math.log(1-prob);
//			System.out.println(i+"\t"+expt+"\t"+prob+"\t"+expt2+"\t"+(expt2/expt));
//		}
//		for(double i=-16;i<=0;i++) {
//			double expt = Math.pow(10, i);
//			double prob = -Math.expm1(-expt);
//			double expt2 = -Math.log1p(-prob);
//			System.out.println(i+"\t"+expt+"\t"+prob+"\t"+expt2+"\t"+(expt2/expt));
//		}
//		System.exit(0);
		
//		// second term does not add anything
//		double v1=0.9999999999997521;
//		double v2=5e-5*1.1e-12; // = 5.5e-17
//		System.out.println(v2);
//		System.err.println(v1+v2);
//		System.exit(0);
		
		
		String rootDir = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/Simulations/";
		
			
			
		// Delete or Keep this method?
//		generateDOLE_ReportPages();
		
		// BPT Single aperiodicity simulations
////		double[] aperArray = {0.1,0.3,0.5,0.7,0.9};
//		double[] aperArray = {0.2,0.4,0.6,0.8,1.0};
//		for(double aper:aperArray) {
//			TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//			u3ProbModel.setRenewalModelChoice(RenewalModels.BPT);
//			u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.SINGLE_VALUED);
//			((AperiodicityModel.SingleValued)u3ProbModel.getAperiodicityModel()).setValue(aper);
//			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
//			File parentDir = new File(rootDir+"bptSimulationsUS26_SingleAper/");
//			if(!parentDir.exists()) 
//				parentDir.mkdir();
//			int numYrs = 50000;
//			String aperString = Double.toString(aper).replace("0.", "_pt");
//			for(int i=1; i<2;i++) {
//				File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"+aperString); 
//				long seed = 984087634+i*1000;
//				String inputFile=null;
//				if(i==1)
//					inputFile = rootDir+"bptSimulationsUS26/Run1_50000yrs/outputTimesinceLast.txt";
//				else
//					inputFile = rootDir+"bptSimulationsUS26_SingleAper/Run"+(i-1)+"_"+numYrs+"yrs"+aperString+"/outputTimesinceLast.txt"; // make simulations consecutive
//
//				LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//						seed, true, true, Double.NaN);
//				LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//			}
//		}
		
		
//		// BPT SIMULTATIONS Slip-rate dependent aperiodicity 
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_SLIPRATE_TEST);
//		u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
//		File parentDir = new File(rootDir+"BPT_SimulationsUS26_SlipRateDepAperTest/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"BPT_SimulationsUS26_SlipRateDepAperTest/Run"+(i-1)+"_"+numYrs+"yrs"+"/outputTimesinceLast.txt"; // make simulations consecutive
//
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
		// Weibull SIMULTATIONS Slip-rate dependent aperiodicity 
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		u3ProbModel.setRenewalModelChoice(RenewalModels.WEIBULL);
//		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_SLIPRATE_TEST);
//		u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
//		File parentDir = new File(rootDir+"WeibullSimulationsUS26_SlipRateDepAperTest/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				setDOLE_asFractionOfRI(erf, 0.75);
//			else
//				inputFile = rootDir+"WeibullSimulationsUS26_SlipRateDepAperTest/Run"+(i-1)+"_"+numYrs+"yrs"+"/outputTimesinceLast.txt"; // make simulations consecutive
//
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			if(i==1)
//				setDOLE_asFractionOfRI(erf, 0.75); // this has to be redone
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}
		

//		// Weibull SIMULTATIONS fixed/single aperiodicity
		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
		u3ProbModel.setRenewalModelChoice(RenewalModels.WEIBULL);
		//	double[] aperArray = {0.1,0.3,0.5,0.7,0.9};
		double[] aperArray = {0.4};
		for(double aper:aperArray) {

			u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.SINGLE_VALUED);
			((AperiodicityModel.SingleValued)u3ProbModel.getAperiodicityModel()).setValue(aper);
			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
			File parentDir = new File(rootDir+"weibullSimulationsUS26_SingleAper/");
			if(!parentDir.exists()) 
				parentDir.mkdir();
			int numYrs = 50000;
			String aperString = Double.toString(aper).replace("0.", "_pt");
			for(int i=2; i<3;i++) {
				File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"+aperString); 
				long seed = 984087634+i*1000;
				String inputFile=null;
				if(i==1)
					setDOLE_asFractionOfRI(erf, 0.75);
				else
					inputFile = rootDir+"weibullSimulationsUS26_SingleAper/Run"+(i-1)+"_"+numYrs+"yrs"+aperString+"/outputTimesinceLast.txt"; // make simulations consecutive

				LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
						seed, true, true, Double.NaN);
				if(i==1)
					setDOLE_asFractionOfRI(erf, 0.75); // this has to be redone
				LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
			}
		}

		
//		// FULL US2026 Weibull SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		String aperSuffix ="";
//		// Alt Aper?
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_LOW);
////		aperSuffix = "_aperLOW";
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_HIGH);
////		aperSuffix = "_aperHIGH";
//		u3ProbModel.setRenewalModelChoice(RenewalModels.WEIBULL);		
//		File parentDir = new File(rootDir+"weibullSimulationsUS26"+aperSuffix+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				setDOLE_asFractionOfRI(erf, 0.75);
////				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"weibullSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//
//			setDOLE_asFractionOfRI(erf, 0.75); // need to do again
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
//		// FULL US2026 Lognormal SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		String aperSuffix ="";
//		// Alt Aper?
//		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_LOW);
//		aperSuffix = "_aperLOW";
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_HIGH);
////		aperSuffix = "_aperHIGH";
//		u3ProbModel.setRenewalModelChoice(RenewalModels.LOGNORMAL);		
//		File parentDir = new File(rootDir+"lognormalSimulationsUS26"+aperSuffix+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"lognormalSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			setDOLE_asFractionOfRI(erf, 0.75); // needed this for LOW APER (absolutely all sections need DOLE)
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
//		// FULL US2026 BPT SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		String aperSuffix ="";
//		// Alt Aper?
////		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_LOW);
////		aperSuffix = "_aperLOW";
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_HIGH);
////		aperSuffix = "_aperHIGH";
//		
//		File parentDir = new File(rootDir+"bptSimulationsUS26"+aperSuffix+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"bptSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
////			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
////					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}
//		// do plots for multiple simulations
//		LongTermTD_Simulator.sectPlotsForMultSimulations(rootDir+"bptSimulationsUS26/", "Run", "_"+numYrs+"yrs", 20, erf.getSolution(), numYrs);
//		LongTermTD_Simulator.rupPlotsForMultSimulations(rootDir+"bptSimulationsUS26/", "Run", "_"+numYrs+"yrs", 20, erf, numYrs);
	
		
		// POISSON SIMULATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
//		File parentDir = new File(rootDir+"poissonSimulationsUS26/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=76; i<81;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			LongTermTD_Simulator.simulateEvents(erf, null,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, null, numYrs, outputDir, true);	
//		}
//		LongTermTD_Simulator.sectPlotsForMultSimulations(rootDir+"poissonSimulationsUS26/", "Run", "_"+numYrs+"yrs", 60, erf.getSolution(), numYrs);
//		LongTermTD_Simulator.rupPlotsForMultSimulations(rootDir+"poissonSimulationsUS26/", "Run", "_"+numYrs+"yrs", 60, erf, (double)numYrs);


		
//		// TEST OF BIAS CORR - FULL US2026 BPT SIMULTATIONS - 
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		String aperSuffix ="";
//		// Alt Aper?
////		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_LOW);
////		aperSuffix = "_aperLOW";
////		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.NSHM26_HIGH);
////		aperSuffix = "_aperHIGH";
//		
//		File parentDir = new File(rootDir+"bptSimulationsUS26_BiasCorrTest"+aperSuffix+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=20; i<21;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"bptSimulationsUS26_BiasCorrTest/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}
//		LongTermTD_Simulator.sectPlotsForMultSimulations(rootDir+"bptSimulationsUS26_BiasCorrTest/", "Run", "_"+numYrs+"yrs", 10, erf.getSolution(), numYrs);
//		LongTermTD_Simulator.rupPlotsForMultSimulations(rootDir+"bptSimulationsUS26_BiasCorrTest/", "Run", "_"+numYrs+"yrs", 10, erf, numYrs);

		
		// FULL US2026 BPT SIMULTATIONS WITH BPTAveragingTypeOptions.AVE_RI_AVE_TIME_SINCE
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		BPTAveragingTypeOptions choice = BPTAveragingTypeOptions.AVE_RI_AVE_TIME_SINCE;
//		u3ProbModel.setAveragingTypeChoice(choice);
//
//		File parentDir = new File(rootDir+"bptSimulationsUS26_"+choice.getCompactLabel()+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"bptSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
////			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
////					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
		
//		// JAMIE TEST RUN, US2026 BPT SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		
//		// need to change averaging type
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		BPTAveragingTypeOptions choice = BPTAveragingTypeOptions.AVE_RI_AVE_TIME_SINCE;
//		u3ProbModel.setAveragingTypeChoice(choice);
//
//		// override the aveCondRecurIntervalForFltSysRups
//		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
//		double[] aveCondRecurIntervalForFltSysRupsAlt = TimeDepUtils.testJamieAveCondRecurIntervalForFltSysRups(erf.getSolution());
//		// REMOVE FOLLOWING METHOD IF THIS TEST IS NO GOOD
//		((UCERF3_ProbabilityModel)probModel).tempSetAveCondRecurIntervalForFltSysRups(aveCondRecurIntervalForFltSysRupsAlt);
//
//		String aperSuffix ="";
//		File parentDir = new File(rootDir+"Jamie_Test_bptSimulationsUS26"+aperSuffix+"/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<2;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				inputFile = rootDir+"bptSimulationsUS26/Run1"+"_"+numYrs+"yrs/outputTimesinceLast.txt"; 
//			else
//				inputFile = rootDir+"bptSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
	}
}
