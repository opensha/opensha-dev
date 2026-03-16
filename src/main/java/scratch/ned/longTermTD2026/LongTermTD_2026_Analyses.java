package scratch.ned.longTermTD2026;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
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
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen.DataToInclude;

import scratch.UCERF3.erf.utils.ProbModelsPlottingUtils;
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

		return erf;
	}
	
	
//	private static void testReportPages(TimeDepFaultSystemSolutionERF erf) {
//		long currentTimeEpoch = System.currentTimeMillis();
//		String dateString = new java.text.SimpleDateFormat("MM_dd_yyyy").format(new java.util.Date (currentTimeEpoch)); // Epoch in seconds, remove '*1000' for milliseconds
//		String full_FSS_fileName = "/Users/field/nshm-haz_data/fullPrefUS_FSS.zip";
//
//		FaultSystemSolution sol = FSS_Fetcher2023.getPreferredFull_FSS(full_FSS_fileName);		
//		
//		File tdMainDir = new File("/Users/field/markdown/nshm23_time_dependence_"+dateString);
//
//			try {
//				TimeDependentReportPageGen.generatePage(new File(tdMainDir, "allDOLE_fullParent"), erf, 
//						PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.ALL_DATA, new File(tdMainDir, "allDOLE_fullParent"));
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//
//	}
//	
//	
	private static void generateReportPage(String outputDir, String fss_fileNameWithPath, int startYear, int duration, 
			int histOpenIntYear, FSS_ProbabilityModels probModChoice, AperiodicityModels aperModelChoice, RenewalModels 
			renewalModelChoice, BPTAveragingTypeOptions averagingChoice, PaleoMappingAlgorithm paleoMapping, 
			DataToInclude paleoDataToInclude, String referenceDir, String titleString, String infoString) {
		
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

		try {
			TimeDependentReportPageGen.generatePage(new File(outputDir), erf, 
					paleoMapping, paleoDataToInclude, new File(referenceDir), titleString, infoString);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// list source index involving section index
//		int targetSectID = 7898;  // Prince William sound
//		int targetSectID = 6920;  // Susitna Glacier Subsection
		int targetSectID = 5594;  // New Madrid - SSCn (New Madrid west)
		for(int s=0; s<erf.getNumFaultSystemSources(); s++) {
			int fssRupID = erf.getFltSysRupIndexForSource(s);
			List<Integer> sectIDList = erf.getSolution().getRupSet().getSectionsIndicesForRup(fssRupID);
			if(sectIDList.contains(targetSectID)) {
				System.out.println(s+"\tsource utilizes sect "+targetSectID+"\t"+sectIDList.size()+"\t"+erf.getSource(s).getName()+"\t"+sectIDList);
			}
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
		
//		outputDir = referenceDir;
//		titleString = "Reference 2026 Time-Dependent Model";
//		infoString = "This is the Reference (preferred) model. ";
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


		outputDir = rootDir+"NoOpenIntervalModel/";
		histOpenIntYear=startYear;
		titleString = "Reference, But With No Open Interval";
		infoString = "This is the Reference model, but where the open interval year is "+histOpenIntYear+".";
		generateReportPage(outputDir,  fss_fileNameWithPath,  startYear,  duration, 
				 histOpenIntYear,  probModChoice,  aperModelChoice,  
				renewalModelChoice,  averagingChoice, paleoMapping,  paleoDataToInclude,  
				referenceDir, titleString, infoString);

//		outputDir = rootDir+"WeibullModel/";
//		referenceDir = rootDir+"NoOpenIntervalModel/";
//		titleString = "Reference, But With Weibull Renewal Model";
//		renewalModelChoice = RenewalModels.WEIBULL;
//		infoString = "This is the Reference model, but with the renewal model switched to Weibull and no open interval. ";
//		histOpenIntYear=startYear;
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


	public static void main(String[] args) {

//		// RENEWAL MODEL TEST PLOTS
//		// Hard-code-adjust the parameter EqkProbDistCalc.NUMERICAL_PRECISION 
//		// to generate different directory/results to see the effect on curves
//		// (this was done to determine the default value)
//		generateRenewalModelPlots(true);
//		generateRenewalModelPlots(false);
		
		// This tests expm1 and log1p
		for(double i=-16;i<=0;i++) {
			double expt = Math.pow(10, i);
			double prob = 1-Math.exp(-expt);
			double expt2 = -Math.log(1-prob);
			System.out.println(i+"\t"+expt+"\t"+prob+"\t"+expt2+"\t"+(expt2/expt));
		}
		for(double i=-16;i<=0;i++) {
			double expt = Math.pow(10, i);
			double prob = -Math.expm1(-expt);
			double expt2 = -Math.log1p(-prob);
			System.out.println(i+"\t"+expt+"\t"+prob+"\t"+expt2+"\t"+(expt2/expt));
		}
		System.exit(0);
		
//		// second term does not add anything
//		double v1=0.9999999999997521;
//		double v2=5e-5*1.1e-12; // = 5.5e-17
//		System.out.println(v2);
//		System.err.println(v1+v2);
//		System.exit(0);
		
//		generatePreliminaryResults();
		
		String rootDir = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/";
		
			
			
		// Delete or Keep this method?
//		generateDOLE_ReportPages();
		
		// BPT Single aperiodicity simulations
//		double[] aperArray = {0.1,0.3,0.5,0.7,0.9};
		double[] aperArray = {0.2,0.4,0.6,0.8,1.0};
		for(double aper:aperArray) {
			TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
			u3ProbModel.setRenewalModelChoice(RenewalModels.BPT);
			u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.SINGLE_VALUED);
			((AperiodicityModel.SingleValued)u3ProbModel.getAperiodicityModel()).setValue(aper);
			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
			File parentDir = new File(rootDir+"bptSimulationsUS26_SingleAper/");
			if(!parentDir.exists()) 
				parentDir.mkdir();
			int numYrs = 50000;
			String aperString = Double.toString(aper).replace("0.", "_pt");
			for(int i=1; i<2;i++) {
				File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"+aperString); 
				long seed = 984087634+i*1000;
				String inputFile=null;
				if(i==1)
					inputFile = rootDir+"bptSimulationsUS26/Run1"+"_"+numYrs+"yrs/outputTimesinceLast.txt";
				else
					inputFile = rootDir+"bptSimulationsUS26_SingleAper/Run"+(i-1)+"_"+numYrs+"yrs/outputTimesinceLast.txt"; // make simulations consecutive

				LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
						seed, true, true, Double.NaN);
				if(i==1)
					setDOLE_asFractionOfRI(erf, 0.6); // this has to be redone
				LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
			}
		}
		
		
//		// FULL US2026 Weibull SIMULTATIONS fixed/single aperiodicity
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		u3ProbModel.setRenewalModelChoice(RenewalModels.WEIBULL);
//		double aper = 0.2;
//		u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.SINGLE_VALUED);
//		((AperiodicityModel.SingleValued)u3ProbModel.getAperiodicityModel()).setValue(aper);
//		u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(erf.getTimeSpan().getStartTimeYear(), true));
//		File parentDir = new File(rootDir+"weibullSimulationsUS26/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		String aperString = Double.toString(aper).replace("0.", "_pt");
//		for(int i=2; i<3;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"+aperString); 
//			long seed = 984087634+i*1000;
//			String inputFile=null;
//			if(i==1)
//				setDOLE_asFractionOfRI(erf, 0.75);
//////				inputFile = rootDir+"poissonSimulationsUS26/Run1_1000000yrs/outputTimesinceLast.txt";
////			inputFile = rootDir+"bptSimulationsUS26/Run1_50000yrs/outputTimesinceLast.txt";
//			else
//				inputFile = rootDir+"weibullSimulationsUS26/Run"+(i-1)+"_"+numYrs+"yrs"+aperString+"/outputTimesinceLast.txt"; // make simulations consecutive
//
//			LongTermTD_Simulator.simulateEvents(erf, inputFile,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			if(i==1)
//				setDOLE_asFractionOfRI(erf, 0.75); // this has to be redone
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
//		// FULL US2026 Lognormal SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)erf.getProbabilityModel();
//		u3ProbModel.setRenewalModelChoice(RenewalModels.LOGNORMAL);		
//		File parentDir = new File(rootDir+"lognormalSimulationsUS26/");
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
//			LongTermTD_Simulator.generateSimulationPlots(erf, inputFile, numYrs, outputDir, true);	
//		}

		
//		// FULL US2026 BPT SIMULTATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		File parentDir = new File(rootDir+"bptSimulationsUS26/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=2; i<21;i++) {
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

//		// POISSON SIMULATIONS
//		TimeDepFaultSystemSolutionERF erf = getFullPrefUS26_ERF();
//		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
//		File parentDir = new File(rootDir+"poissonSimulationsUS26/");
//		if(!parentDir.exists()) 
//			parentDir.mkdir();
//		int numYrs = 50000;
//		for(int i=1; i<21;i++) {
//			File outputDir = new File(parentDir,"Run"+i+"_"+numYrs+"yrs"); 
//			long seed = 984087634+i*1000;
//			LongTermTD_Simulator.simulateEvents(erf, null,"outputTimesinceLast.txt", numYrs, outputDir, 
//					seed, true, true, Double.NaN);
//			LongTermTD_Simulator.generateSimulationPlots(erf, null, numYrs, outputDir, true);	
//		}

	}

}
