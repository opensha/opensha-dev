package scratch.ned.longTermTD2026;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.HistoricalRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoDOLE_Data;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.PaleoMappingAlgorithm;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.TimeDependentReportPageGen.DataToInclude;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbModelsPlottingUtils;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.ned.nshm23.CONUS_TD_ERF_Demo;

public class LongTermTD2026_Analyses {
	
	/**
	 * This verifies that nth rup indices do not get messed up by turning on and off background seismicity
	 */
	private static void test_nthRupState() {
		String fileName="/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fileName);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
		erf.getParameter(MagDependentAperiodicityParam.NAME).setValue(MagDependentAperiodicityOptions.MID_VALUES);
		erf.setParameter(BPTAveragingTypeParam.NAME, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
		erf.updateForecast();	
		
		int numRup = erf.getTotNumRups();
		int[] srcIndex = new int[numRup];
		int[] fssRupIndex = new int[numRup];
		int[]  rupIndexInSource = new int[numRup];
		
		for(int nthRup=0;nthRup<numRup;nthRup++) {
			srcIndex[nthRup] = erf.getFltSysRupIndexForNthRup(nthRup);
			rupIndexInSource[nthRup] = erf.getRupIndexInSourceForNthRup(nthRup);
			fssRupIndex[nthRup] = erf.getFltSysRupIndexForNthRup(nthRup);
		}
		
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.INCLUDE);
		erf.updateForecast();	
		for(int nthRup=0;nthRup<numRup;nthRup++) {
			if(srcIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(rupIndexInSource[nthRup] != erf.getRupIndexInSourceForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(fssRupIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
		}

		
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.updateForecast();	
		for(int nthRup=0;nthRup<numRup;nthRup++) {
			if(srcIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(rupIndexInSource[nthRup] != erf.getRupIndexInSourceForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(fssRupIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
		}
		
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.INCLUDE);
		erf.updateForecast();	
		for(int nthRup=0;nthRup<numRup;nthRup++) {
			if(srcIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(rupIndexInSource[nthRup] != erf.getRupIndexInSourceForNthRup(nthRup))
				throw new RuntimeException("Problem");
			if(fssRupIndex[nthRup] != erf.getFltSysRupIndexForNthRup(nthRup))
				throw new RuntimeException("Problem");
		}
		System.out.println("success!");


	}

	
	private static void makeTestTD_CalculationFiles() {
		
		long currentTimeEpoch = System.currentTimeMillis();
		String dateString = new java.text.SimpleDateFormat("MM_dd_yyyy").format(new java.util.Date (currentTimeEpoch)); // Epoch in seconds, remove '*1000' for milliseconds
		File outputDir1 = new File("/Users/field/markdown/nshm23_time_dependence_"+dateString);
		if(!outputDir1.exists()) outputDir1.mkdir();
		System.out.println(!outputDir1.exists()+"\t"+outputDir1);

		File outputDir2 = new File(outputDir1,"TD_TestFiles");
		if(!outputDir2.exists()) outputDir2.mkdir();
		System.out.println(!outputDir2.exists()+"\t"+outputDir2);

		String fssFileName = "results_WUS_FM_v3_branch_averaged_gridded_simplified.zip";
		int startYear = 2025;
		int openIntYear = 1875;
		double openInt = startYear-openIntYear;
//		double openInt = 0;
		double forecastDurationYears = 50;
		ProbabilityModelOptions probModel = ProbabilityModelOptions.U3_BPT;
		MagDependentAperiodicityOptions aperModel = MagDependentAperiodicityOptions.MID_VALUES;
		
		String stringForReadmeFile =
				"This directory contains a verification files for long-term time-dependent "+
				"calculations for the following settings:\n\n"+
				"Only Historic Rupture DOLE data (no paleoseismic)\n" +
				"Fault System Solution filename = "+fssFileName+"\n"+
				"Start Year = "+startYear+"\n" +
				"Forecast Duration = "+forecastDurationYears+" years\n"+
				"Open Interval Year = "+openIntYear+"\n" +
				"Open Interval = "+openInt+"\n" +
				"Probability Model = "+probModel+"\n" +
				"Aperiodicity Model:\n\n";
				
				for(int i=0;i<aperModel.getAperValuesArray().length;i++) {
					if(i==0)
						stringForReadmeFile += "\t"+aperModel.getAperValuesArray()[i]+"\t for M ≤ "+aperModel.getAperMagBoundariesArray()[i]+"\n";
					else if (i == aperModel.getAperValuesArray().length-1)
						stringForReadmeFile += "\t"+aperModel.getAperValuesArray()[i]+"\t for M > "+aperModel.getAperMagBoundariesArray()[i-1]+"\n";
					else
						stringForReadmeFile += "\t"+aperModel.getAperValuesArray()[i]+"\t for "+aperModel.getAperMagBoundariesArray()[i-1]+" ≤ M > "+aperModel.getAperMagBoundariesArray()[i]+"\n";
				}		
				stringForReadmeFile += "\n This directory was generated by running the OpenSHA java method \n\n\tscratch.ned.longTermTD2026.longTermTD2026_Analyses.makeTestTD_CalculationFiles() \n\non "+dateString;
		
		// get solution
		FaultSystemSolution sol = null;
		try {
			sol = FaultSystemSolution.load(new File("/Users/field/nshm-haz_data/"+fssFileName));
			// set historical rup DOLE
			FaultSystemRupSet rupSet = sol.getRupSet();
			List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
			List<PaleoDOLE_Data> paleoDataList = new ArrayList<PaleoDOLE_Data>();
			List<HistoricalRupture> histRupDataList = DOLE_SubsectionMapper.loadHistRups();			
			DOLE_SubsectionMapper.mapDOLE(subSects, histRupDataList, paleoDataList, PaleoMappingAlgorithm.NEIGHBORING_SECTS, false);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
		erf.setParameter(ProbabilityModelParam.NAME, probModel);
		erf.setParameter(MagDependentAperiodicityParam.NAME, aperModel);
		erf.setParameter(HistoricOpenIntervalParam.NAME, openInt); 
		erf.getTimeSpan().setStartTime(startYear);
		erf.getTimeSpan().setDuration(forecastDurationYears);
		erf.updateForecast();
		
List<Integer> testSectID_List = sol.getRupSet().getSectionsIndicesForRup(174250);	
for(int s:testSectID_List) {
	long dateLastEpoch = sol.getRupSet().getFaultSectionData(s).getDateOfLastEvent();
	double dateLastYrs = (double)dateLastEpoch/ProbabilityModelsCalc.MILLISEC_PER_YEAR + 1970.0;
	if(dateLastEpoch == Long.MIN_VALUE)
		dateLastYrs = Double.NaN;
	System.out.println(s+"\t"+dateLastYrs);
}
System.exit(0);
		
		ProbabilityModelsCalc probCalc = new ProbabilityModelsCalc(erf);
		
		String headerString = "";
		headerString += "srcIndex,";	// added here
		headerString += "longTermRate,";	// added here
		headerString += "gainTest,";	// added here
		headerString += "numRup,";	// added here
		// from ProbabilityModelsCalc
		headerString += "fltSysRupIndex,";
		headerString += "probGain,";
		headerString += "condProb,";
		headerString += "rupMag,";
		headerString += "aveCondRecurInterval,";
		headerString += "aveCondRecurIntervalWhereUnknown,";		
		headerString += "aveTimeSinceLastWhereKnownYears,";
		headerString += "aveNormTimeSinceLastEventWhereKnown,";
		headerString += "totRupArea,";
		headerString += "totRupAreaWithDateOfLast,";
		headerString += "fractRupAreaWithDateOfLast,";
		headerString += "aperValue[getAperIndexForRupMag(rupMag)],";
		headerString += "numSubsectForRup,";


//		double minRatio = Double.MAX_VALUE;
//		double maxRatio = -Double.MAX_VALUE;

		FileWriter fw;
		FileWriter fw_readme;
		try {
			fw = new FileWriter(new File(outputDir2+"/rupTestFile.csv"));
			fw_readme = new FileWriter(new File(outputDir2+"/INFO.txt"));
			fw_readme.write(stringForReadmeFile);
			fw.write(headerString+"\n"); 
//			System.out.println(headerString);
			for(int s=0;s<erf.getNumFaultSystemSources();s++) {
				if(erf.getSource(s).isSourcePoissonian())
					throw new RuntimeException("source is poissonian: "+s+"\t"+erf.getSource(s).getName());
				int fsrIndex = erf.getFltSysRupIndexForSource(s);
				double probGainTest = probCalc.getU3_ProbGainForRup(fsrIndex, openInt, false, true, true, 
						erf.getTimeSpan().getStartTimeInMillis(), forecastDurationYears);
				double testRate=0;  // this should equal fss rate time gain
				for(int r=0;r<erf.getSource(s).getNumRuptures();r++) {
					double prob = erf.getSource(s).getRupture(r).getProbability();
// this not the right way because U3 TD does not use Poisson calc
//					double rate;
//					if(prob<1e-8)
//						rate = prob/forecastDurationYears;	// avoiding numerical issues
//					else
//						rate = erf.getSource(s).getRupture(r).getMeanAnnualRate(forecastDurationYears);
					testRate+=prob/forecastDurationYears;
				}
				testRate /= probGainTest;
				if((float)testRate != (float)sol.getRateForRup(fsrIndex))
					throw new RuntimeException("long-term rate discrepancy for srcID="+s);
//				double rate = sol.getRateForRup(fsrIndex);
//				double ratio = testRate/rate;
//				if(minRatio>ratio) minRatio=ratio;
//				if(maxRatio<ratio) maxRatio=ratio;	
				
				String line = s+","+testRate+","+probGainTest+","+erf.getSource(s).getNumRuptures()+
						","+probCalc.u3_ProgGainForRupInfoString+"\n";
				fw.write(line); 
				
			}
			fw.close();
			fw_readme.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		probCalc.writeCurrentSectDataToCSV_File(outputDir2, "sectTestFile.csv");
		System.out.println("DONE (tests checked out)");
//		System.out.println("minRatio="+(float)minRatio);
//		System.out.println("maxRatio="+(float)maxRatio);
	}
	
	/**
	 * this returns a CEUS ERF with background seismicity excluded but all other 
	 * parameters set as default; erf.updateForecast()is not called
	 * @return
	 */
	private static FaultSystemSolutionERF getCEUS_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/ceus_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getCEUS_FSS(full_FSS_fileName);		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		return erf;
	}
	
	

	/**
	 * this returns a Alaska (AK) ERF with background seismicity excluded but all other 
	 * parameters set as default; erf.updateForecast()is not called
	 * @return
	 */
	private static FaultSystemSolutionERF getAK_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/alaska_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getAK_FSS(full_FSS_fileName);		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		return erf;
	}
	
	/**
	 * this returns a Cascadia ERF with background seismicity excluded but all other 
	 * parameters set as default; erf.updateForecast()is not called
	 * @return
	 */
	private static FaultSystemSolutionERF getCascadia_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/cascadia_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getCascadia_FSS(full_FSS_fileName);		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		return erf;
	}

	/**
	 * this returns a WUS w/ Cascadia ERF with background seismicity excluded but all other 
	 * parameters set as default; erf.updateForecast()is not called
	 * @return
	 */
	private static FaultSystemSolutionERF getWUS_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/wusWithCascadia_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getWUS_FSS(full_FSS_fileName);		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		return erf;
	}

	
	/**
	 * this returns the ERF with background seismicity excluded but all other 
	 * parameters set as default; erf.updateForecast()is not called
	 * @return
	 */
	private static FaultSystemSolutionERF getUS26_ERF() {
		String full_FSS_fileName = "/Users/field/nshm-haz_data/full_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getFull_FSS(full_FSS_fileName);		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		return erf;
	}
	
	private static void generateDOLE_ReportPages() {

		long currentTimeEpoch = System.currentTimeMillis();
		String dateString = new java.text.SimpleDateFormat("MM_dd_yyyy").format(new java.util.Date (currentTimeEpoch)); // Epoch in seconds, remove '*1000' for milliseconds
		String full_FSS_fileName = "/Users/field/nshm-haz_data/full_FSS_test.zip";
		FaultSystemSolution sol = CONUS_TD_ERF_Demo.getFull_FSS(full_FSS_fileName);		
		
		File tdMainDir = new File("/Users/field/markdown/nshm23_time_dependence_"+dateString);

		try {
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "allDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.ALL_DATA);

			// recreating solution to avoid propagating previous DOLE settings
//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "allDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.ALL_DATA);

//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "forDebugging_onlyPaleoDOLE_nearestSubsect"), sol, PaleoMappingAlgorithm.CLOSEST_SECT, DataToInclude.PALEO_ONLY);
			
//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyHistoricRupDOLE"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.HIST_RUPS_ONLY);

//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyPaleoDOLE_fullParent"), sol, PaleoMappingAlgorithm.FULL_PARENT, DataToInclude.PALEO_ONLY);

//			sol = FaultSystemSolution.load(new File(full_FSS_fileName));
			TimeDependentReportPageGen.generatePage(new File(tdMainDir, "onlyPaleoDOLE_neighbors"), sol, PaleoMappingAlgorithm.NEIGHBORING_SECTS, DataToInclude.PALEO_ONLY);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}

	/**
	 * this verifies that the new method produces the same result as the old one
	 */
	private static void testOldVsNewSimulationMethod() {
		
		String fileName="/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";
		String timeSinceLastFileName = "/Users/field/FilesFromOldComputerTransfer/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/erSimulations/timeSinceLastForSimulation.txt"; 
 // timeSinceLastFileName=null;
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fileName);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
		erf.getParameter(MagDependentAperiodicityParam.NAME).setValue(MagDependentAperiodicityOptions.MID_VALUES);
		BPTAveragingTypeOptions aveType = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
		erf.setParameter(BPTAveragingTypeParam.NAME, aveType);
 // erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();	
		
		ProbabilityModelsCalc testCalc = new ProbabilityModelsCalc(erf);
		long seed = 1234567l;
//		long seed = 4562l;
		long startTime = System.currentTimeMillis();
		double numYrs=1000; // Following took ~20 hrs;  =200000

		// need to move the following resultant dir by hand
		testCalc.testER_Simulation(timeSinceLastFileName, null, numYrs, "OldSimMethodResult", seed);
		
		File outputDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/OldVsNewSimMethodTest/NewSimMethodResult");
//		File outputDir=null;
		boolean makePlots=true;
		testCalc.simulateEvents(timeSinceLastFileName, null, numYrs, outputDir, seed, true, makePlots);
		
		double runtimeMin = (double)(System.currentTimeMillis()- startTime)/60000d;
		System.out.println("runtime (min) = "+(float)runtimeMin);

	}
	
	
	private static void bptSimulations(FaultSystemSolutionERF erf, double numYrs, File parentDir, 
			String dirName, long seed, String timeSinceLastFileName, MagDependentAperiodicityOptions aper) {
		
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
		erf.getParameter(MagDependentAperiodicityParam.NAME).setValue(aper);
		BPTAveragingTypeOptions aveType = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
		erf.setParameter(BPTAveragingTypeParam.NAME, aveType);
		erf.updateForecast();	
		ProbabilityModelsCalc testCalc = new ProbabilityModelsCalc(erf);
		long startTime = System.currentTimeMillis();
		if(!parentDir.exists()) 
			parentDir.mkdir();
		File outputDir = new File(parentDir,dirName);
		boolean makePlots=true;
		testCalc.simulateEvents(timeSinceLastFileName, "outputTimesinceLast.txt", numYrs, outputDir, seed, true, makePlots);
		double runtimeMin = (double)(System.currentTimeMillis()- startTime)/60000d;
		System.out.println("runtime (min) = "+(float)runtimeMin);
	}

	
	
	private static void poissonSimulations(FaultSystemSolutionERF erf, File parentDir, 
			String dirName, long seed, int numYrs) {
		
		String timeSinceLastFileName = null;
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		erf.updateForecast();	
		ProbabilityModelsCalc testCalc = new ProbabilityModelsCalc(erf);
		long startTime = System.currentTimeMillis();
		if(!parentDir.exists()) 
			parentDir.mkdir();
		File outputDir = new File(parentDir,dirName);
		boolean makePlots=true;
		testCalc.simulateEvents(timeSinceLastFileName, "outputTimesinceLast.txt", numYrs, outputDir, seed, true, makePlots);
		double runtimeMin = (double)(System.currentTimeMillis()- startTime)/60000d;
		System.out.println("runtime (min) = "+(float)runtimeMin);
	}


	
	
	/**
	 * This looks at stats over the 10 runs
	 */
	private static void bptSimulationsU3_PartRateStats() {
		
		File outputDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsU3/StatsPlotsFor10_Runs");
		if(!outputDir.exists())
			outputDir.mkdir();
		
		// loop over the ten runs
		int numSections = 2606;
		int numSimulations = 10;
		ArrayList<double[]> simRateArrayList = new ArrayList<double[]>();
		ArrayList<double[]> ratioArrayList = new ArrayList<double[]>();
		double[] targetRateArray = new double[numSections];
		String[] sectNameArray = new String[numSections];
		double totMeanRatio = 0;
		for(int i=1;i<=numSimulations;i++) {
			double[] simRateArray = new double[numSections];
			double[] ratioArray = new double[numSections];
			try {
				File dataFile = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsU3/Run"+i+"_100k/obsVsImposedSectionPartRates.txt");
				BufferedReader reader = new BufferedReader(scratch.UCERF3.utils.UCERF3_DataUtils.getReader(dataFile.toURL()));

//		        Path filePath = Paths.get("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsU3/Run"+i+"_100k/obsVsImposedSectionPartRates.txt"); 
//		        Charset asciiCharset = Charset.forName("US-ASCII");
//		        BufferedReader reader = Files.newBufferedReader(filePath, asciiCharset);

		        reader.readLine(); // skip header
				int s=0;
				String line;
				while ((line = reader.readLine()) != null) {
					String[] st = StringUtils.split(line,"\t");
					int sectIndex = Integer.valueOf(st[0]);
					targetRateArray[s] = Double.valueOf(st[1]);
					simRateArray[s] = Double.valueOf(st[2]);
					double ratio = Double.valueOf(st[3]);
					ratioArray[s] = ratio;
					totMeanRatio += ratio;
					sectNameArray[s] = st[4];
					if(s != sectIndex)
						throw new RuntimeException("bad index");
					s+=1;
				}
				simRateArrayList.add(simRateArray);
				ratioArrayList.add(ratioArray);
				
				reader.close();
			} catch (Exception e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		totMeanRatio /= (numSections*numSimulations);

		double[] simRateMeanArray = new double[numSections];
		double[] simRateStdevArray = new double[numSections];
		double[] simRateStdomArray = new double[numSections];
		double[] ratioMeanArray = new double[numSections];
		double[] ratioStdevArray = new double[numSections];
		double[] ratioStdomArray = new double[numSections];
		double[] ratioMeanFirst5 = new double[numSections];
		double[] ratioMeanLast5 = new double[numSections];
//		double[] simNormError = new double[numSections];
		double[] simFractStdom = new double[numSections];

		double minMeanRatio = Double.MAX_VALUE;
		double maxMeanRatio = -Double.MAX_VALUE;
		for(int s=0;s<numSections;s++) {
			DescriptiveStatistics simRateStats = new DescriptiveStatistics();
			DescriptiveStatistics ratioStats = new DescriptiveStatistics();
			for(int i=0;i<numSimulations;i++) {
				simRateStats.addValue(simRateArrayList.get(i)[s]);
				ratioStats.addValue(ratioArrayList.get(i)[s]);
				if(i<=4)
					ratioMeanFirst5[s]+=ratioArrayList.get(i)[s]/(5.0*totMeanRatio);
				else
					ratioMeanLast5[s]+=ratioArrayList.get(i)[s]/(5.0*totMeanRatio);
			}
			simRateMeanArray[s] = simRateStats.getMean();
			simRateStdevArray[s] = simRateStats.getStandardDeviation();
			simRateStdomArray[s] = simRateStdevArray[s]/Math.sqrt(numSimulations);
			ratioMeanArray[s] = ratioStats.getMean();
			if(minMeanRatio>ratioMeanArray[s]) minMeanRatio=ratioMeanArray[s];
			if(maxMeanRatio<ratioMeanArray[s]) maxMeanRatio=ratioMeanArray[s];
			ratioStdevArray[s] = ratioStats.getStandardDeviation();
			ratioStdomArray[s] = ratioStats.getStandardDeviation()/Math.sqrt(numSimulations);
//			simNormError[s] = (simRateMeanArray[s]/1.088-targetRateArray[s])/simRateStdomArray[s];
			simFractStdom[s] = simRateStdomArray[s]/simRateMeanArray[s];
		}
		
		// compute min rate
		double minRate = Double.MAX_VALUE;
		double maxFractError = -Double.MAX_VALUE;
		int indexForMaxFractError=-1;
		for(int r=0; r<simRateMeanArray.length;r++) {
			if(minRate>simRateMeanArray[r]) minRate=simRateMeanArray[r];
			if(maxFractError<simFractStdom[r]) {
				maxFractError=simFractStdom[r];
				indexForMaxFractError=r;
			}
		}
		System.out.println("minRate="+minRate+"\tN="+(float)minRate*1e6);
		System.out.println("maxFractStdom="+maxFractError+"\ttargetRate="+targetRateArray[indexForMaxFractError]+"\nValues for maxFractStdom:");
		for(int i=0;i<numSimulations;i++) {
			System.out.println(simRateArrayList.get(i)[indexForMaxFractError]);
		}
		
		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(ratioMeanFirst5, ratioMeanLast5, outputDir, "RatioSimFirst5_vs_RatioSimLast5_SectionPartRates", 
				"", "RatioSimFirst5 Sect Part Rate (/yr)", "RatioSimLast5 Sect Part Rate (/yr)", 0.5,2.0);
		
		PearsonsCorrelation pCorr = new PearsonsCorrelation();
		System.out.println("pCorr="+pCorr.correlation(ratioMeanFirst5, ratioMeanLast5));


		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(targetRateArray, simRateMeanArray, outputDir, "aveSimVsImposedSectionPartRates", 
				"", "Imposed Sect Part Rate (/yr)", "Ave Simulated Sect Part Rate (/yr)", Double.NaN,Double.NaN);

//		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(targetArray, ratioMeanArray, outputDir, "aveSimVsImposedRatioSectionPartRates", 
//				"", "Imposed Sect Part Rate (/yr)", "AveSimulated/Imposed Sect Part Rate (/yr)", Double.NaN,Double.NaN);
//		
		DefaultXY_DataSet ratioFunc = new DefaultXY_DataSet(targetRateArray,ratioMeanArray);
		ratioFunc.setInfo("minMeanRatio="+(float)minMeanRatio+"\nmaxMeanRatio="+(float)maxMeanRatio);
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		funcs.add(ratioFunc);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 0.5f, Color.RED));

		Range xAxisRange = new Range(1e-5,4e-2);
		Range yAxisRange = new Range(0,2);
		boolean logX = true;
		boolean logY = false;
		double widthInches = 7.0; // inches
		double heightInches = 6.0; // inches
		boolean popupWindow = false;
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcs,plotChars,"","Imposed Sect Part Rate (/yr)","AveSimulated/Imposed Sect Part Rate (/yr)",xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, new File(outputDir,"aveSimVsImposedRatioSectionPartRates"), popupWindow);

		
		DefaultXY_DataSet normErrorFunc = new DefaultXY_DataSet(targetRateArray,simFractStdom);
		ArrayList<XY_DataSet> funcs2 = new ArrayList<XY_DataSet>();
		funcs2.add(normErrorFunc);
		xAxisRange = new Range(1e-5,4e-2);
		yAxisRange = null;
		logX = true;
		logY = false;
		widthInches = 7.0; // inches
		heightInches = 6.0; // inches
		popupWindow = false;
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcs2,plotChars,"","Imposed Sect Part Rate (/yr)","simFractStdom",xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, new File(outputDir,"simFractStdomVsImposedRatioSectionPartRates"), popupWindow);

		
		DefaultXY_DataSet ratioVsNormStdomFunc = new DefaultXY_DataSet(simFractStdom,ratioMeanArray);
		ArrayList<XY_DataSet> funcs3 = new ArrayList<XY_DataSet>();
		funcs3.add(ratioVsNormStdomFunc);
		xAxisRange = new Range(0,0.005);
		yAxisRange = new Range(0.8,1.2);
		logX = false;
		logY = false;
		widthInches = 7.0; // inches
		heightInches = 6.0; // inches
		popupWindow = false;
		ProbModelsPlottingUtils.writeAndOrPlotFuncs(funcs3,plotChars,"","simFractStdom","AveSimulated/Imposed Sect Part Rate (/yr)",xAxisRange,yAxisRange,
				logX,logY,widthInches,heightInches, new File(outputDir,"simRateRatioVsSimFractStdom"), popupWindow);

		// look only at data with very small fractStdom
		int num =0;
		double fractStdomForSubset=0.005;
		for(int s=0;s<simFractStdom.length;s++)
			if(simFractStdom[s]<fractStdomForSubset) num+=1;
		double[] subsetObsRate = new double[num];
		double[] subsetTargetRate = new double[num];
		double[] subsetRatioMeanFirst5 = new double[num];
		double[] subsetRatioMeanLast5 = new double[num];

		num=0;
		double aveRatio=0;
		double wtAveRatio=0;
		double totWt=0;
		for(int s=0;s<simFractStdom.length;s++) {
			if(simFractStdom[s]<fractStdomForSubset) {
				subsetObsRate[num] = simRateMeanArray[s];
				subsetTargetRate[num] = targetRateArray[s];
				aveRatio+=subsetObsRate[num]/subsetTargetRate[num];
				subsetRatioMeanFirst5[num] = ratioMeanFirst5[s];
				subsetRatioMeanLast5[num] = ratioMeanLast5[s];
				num+=1;
			}
			wtAveRatio += (simRateMeanArray[s]/targetRateArray[s])/(simFractStdom[s]*simFractStdom[s]);
			totWt += 1/(simFractStdom[s]*simFractStdom[s]);
		}
		aveRatio /= subsetObsRate.length;
		wtAveRatio /= totWt;
		System.out.println("aveRatio for select data: "+aveRatio);
		System.out.println("wtAveRatio: "+ wtAveRatio);
		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(subsetTargetRate, subsetObsRate, outputDir, "selectAveSimVsImposedSectionPartRates", 
				"", "Select Imposed Part Rate (/yr)", "Select Ave Sim Part Rate (/yr)", Double.NaN,Double.NaN);

		ProbModelsPlottingUtils.writeSimVsImposedRateScatterPlot(subsetRatioMeanFirst5, subsetRatioMeanLast5, outputDir, "selectRatioSimFirst5_vs_RatioSimLast5_SectionPartRates", 
				"", "Select RatioSimFirst5 Sect Part Rate (/yr)", "Select RatioSimLast5 Sect Part Rate (/yr)", 0.5,2.0);

		
		// make map
		FaultSystemSolution sol=null;
		try {
			sol = FaultSystemSolution.load(new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
			FaultSystemRupSet rupSet = sol.getRupSet();
			List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
			Color[] sectColorArray = new Color[subSects.size()];
			GeographicMapMaker mapMaker = new GeographicMapMaker(subSects);
			mapMaker.setWriteGeoJSON(true);
			mapMaker.clearSectScalars();
			
			List<Color> sectColorList = new ArrayList<>();
			for (int i=0;i<ratioMeanArray.length;i++) {
				Color color = ProbModelsPlottingUtils.getRatioMapColor(ratioMeanArray[i]/totMeanRatio, simFractStdom[i]);
				sectColorList.add(color);
			}
			mapMaker.plotSectColors(sectColorList, null, null);
			mapMaker.setSectNaNChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 0, 255)));
			mapMaker.plot(outputDir, "sectPartRatioMap", " ");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		FileWriter sectRates_fw;
		try {
			sectRates_fw = new FileWriter(new File(outputDir,"/AveSectionPartRatesData.csv"));
			sectRates_fw.write("sectID,simRate,targetRate,ratio,simFractStdom,sectName\n");
			for(int i=0;i<targetRateArray.length;i++) {
				String sectName = sectNameArray[i].replace(","," ");
				sectRates_fw.write(i+","+simRateMeanArray[i]+","+targetRateArray[i]+","+ratioMeanArray[i]+","+simFractStdom[i]+","+sectName+"\n");
			}
			sectRates_fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}


		System.out.println("totMeanRatio="+totMeanRatio+"\nminMeanSectRatio="+minMeanRatio+"\nmaxMeanSectRatio="+maxMeanRatio);


	}
	
	
	public static void junkTest() {
		int numSections = 2287;
		double[] simRateArray = new double[numSections];
		double[] targetRateArray = new double[numSections];
		try {
			File dataFile = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsAK/Run1_100k/magDepSectRates_magBelowAndAbove7pt3.txt");
			BufferedReader reader = new BufferedReader(scratch.UCERF3.utils.UCERF3_DataUtils.getReader(dataFile.toURL()));
			reader.readLine(); // skip header
			int s=0;
			String line;
			while ((line = reader.readLine()) != null) {
				String[] st = StringUtils.split(line,"\t");
				int sectIndex = Integer.valueOf(st[0]);
				if(s != sectIndex)
					throw new RuntimeException("bad index");
				targetRateArray[s] = Double.valueOf(st[7]);
				simRateArray[s] = Double.valueOf(st[8]);
				s+=1;
			}
			reader.close();
		} catch (Exception e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		File dir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsAK/Run1_100k");
		ProbModelsPlottingUtils.writeSimOverImposedVsImposedSecPartRates_Plot(dir, targetRateArray, 
				simRateArray, 10d, 100000d, "Mgt7pt3");


	}



	public static void main(String[] args) {
		
		String rootDir = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/";
		
		
		// WUS w/ Cascadia (Middle branch) Poisson simulations
		File parentDir = new File(rootDir+"bptSimulationsWUS_withCascadia/");
		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.ALL_PT5_VALUES;
		long seed = 984087634;
		int numYrs = 50000;
//		bptSimulations(getWUS_ERF(),numYrs,parentDir, "Run1_aper0pt5", seed, rootDir+"poissonSimulationsWUS_withCascadia/Run1/outputTimesinceLast.txt", aper);
//		bptSimulations(getWUS_ERF(),numYrs,parentDir, "Run2_aper0pt5", seed, rootDir+"bptSimulationsWUS_withCascadia/Run1_aper0pt5/outputTimesinceLast.txt", aper);
		bptSimulations(getWUS_ERF(),numYrs,parentDir, "Run3_aper0pt5", seed, rootDir+"bptSimulationsWUS_withCascadia/Run2_aper0pt5/outputTimesinceLast.txt", aper);
		
		// WUS w/ Cascadia (Middle branch) Poisson simulations
//		File parentDir = new File(rootDir+"poissonSimulationsWUS_withCascadia/");
//		long seed = 984087634;
//		int numYrs = 10000000;
//		for(int i=1;i<=1;i++) {
//			poissonSimulations(getWUS_ERF(), parentDir, "Run"+i, seed+=i*7836271, numYrs);
//			seed += 7836271;
//		}

		

//		// Cascadia (Middle branch) Poisson simulations - THIS DOESN'T WORK BECUASE FULL RUPTURE CAUSES ZERO RATE; IMPLEMENT WITH WUS
//		File parentDir = new File(rootDir+"bptSimulationsCascadia/");
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.ALL_PT5_VALUES;
//		long seed = 984087634;
//		int numYrs = 1000000;
//		bptSimulations(getCascadia_ERF(),numYrs,parentDir, "Run1", seed, rootDir+"poissonSimulationsCascadia/Run1/outputTimesinceLast.txt", aper);

//		// Cascadia (Middle branch) Poisson simulations
//		File parentDir = new File(rootDir+"poissonSimulationsCascadia/");
//		long seed = 984087634;
//		int numYrs = 10000000;
//		for(int i=1;i<=1;i++) {
//			poissonSimulations(getCascadia_ERF(), parentDir, "Run"+i, seed+=i*7836271, numYrs);
//			seed += 7836271;
//		}

		
		//	AK BPT Simulations
//		long seed = 836529;
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.MID_VALUES;
//		double numYrs=100000;
//		File parentDir = new File(rootDir+"bptSimulationsAK/");
//		bptSimulations(getAK_ERF(),numYrs,parentDir, "Run1_100k", seed, rootDir+"poissonSimulationsAK/Run1/outputTimesinceLast.txt", aper);
//		bptSimulations(getAK_ERF(),numYrs,parentDir, "Run2_100k", seed, rootDir+"poissonSimulationsAK/Run2/outputTimesinceLast.txt", aper);
//		bptSimulations(getAK_ERF(),numYrs,parentDir, "Run3_100k", seed, rootDir+"bptSimulationsAK/Run2_100k/outputTimesinceLast.txt",aper);
//		bptSimulations(getAK_ERF(),numYrs,parentDir, "Run4_100k", seed, rootDir+"bptSimulationsAK/Run3_100k/outputTimesinceLast.txt",aper);
		
//		// AK Poisson simulations
//		File parentDir = new File(rootDir+"poissonSimulationsAK/");
//		long seed = 984087634-7836271;
//		int numYrs = 2000000;
//		for(int i=3;i<=3;i++) {
//			poissonSimulations(getAK_ERF(), parentDir, "Run"+i, seed+=i*7836271, numYrs);
//			seed += 7836271;
//		}

		
//		//	CEUS BPT Simulations
//		long seed = 836529;
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.ALL_PT5_VALUES;
//		double numYrs=10000000;
//		File parentDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsCEUS/");
//		bptSimulations(getCEUS_ERF(),numYrs,parentDir, "Run1_10000k_aper0pt5", seed, rootDir+"poissonSimulationsCEUS/Run1/outputTimesinceLast.txt", aper);

//		//	CEUS BPT Simulations
//		long seed = 836529;
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.MID_VALUES;
//		double numYrs=10000000;
//		File parentDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsCEUS/");
//		bptSimulations(getCEUS_ERF(),numYrs,parentDir, "Run1_10000k", seed, rootDir+"poissonSimulationsCEUS/Run1/outputTimesinceLast.txt", aper);

		// CEUS Poisson simulations
//		File parentDir = new File(rootDir+"poissonSimulationsCEUS/");
//		long seed = 984087634;
//		int numYrs = 50000000;
////		int numYrs = 500000;
//		for(int i=1;i<=1;i++) {
//			poissonSimulations(getCEUS_ERF(), parentDir, "Run"+i, seed, numYrs);
//			seed += 7836271;
//		}

		
		//	U3 BPT Simulations - THESE TAKE TOO LONG
//		long seed = 836529;
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.MID_VALUES;
//		FaultSystemSolutionERF erf = getUS26_ERF();
//		double numYrs=1000;
//		File parentDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsUS26/");
//		bptSimulations(erf,numYrs,parentDir, "Run1_1k", seed, rootDir+"poissonSimulationsUS26/Run1/outputTimesinceLast.txt", aper);

		
		// US26 Poisson simulations
//		File parentDir = new File(rootDir+"poissonSimulationsUS26/");
//		long seed = 984087634;
//		int numYrs = 1000000;
//		for(int i=1;i<=10;i++) {
//			poissonSimulations(getUS26_ERF(), parentDir, "Run"+i, seed, numYrs);
//			seed += 7836271;
//		}

		
//		bptSimulationsU3_PartRateStats();	
		
		// STILL NEED TO UPDATE THIS ONE:
//		try {
//		testCalc.testER_NextXyrSimulation(new File("TestRunNextXyrSim"), timeSinceLastFileName, 100, true, seed, 50.0);
//	} catch (IOException e) {
//		e.printStackTrace();
//	}
		

		
		// U3 BPT fixed aper simulations
//		long seed = 836529;
//		String timeSinceFile = rootDir+"poissonSimulations/Run1/outputTimesinceLast.txt";
//		bptSimulationsU3("Aper0pt5", seed, timeSinceFile, MagDependentAperiodicityOptions.ALL_PT5_VALUES);
//		bptSimulationsU3("Aper0pt2", seed, timeSinceFile, MagDependentAperiodicityOptions.ALL_PT2_VALUES);
//		bptSimulationsU3("Aper0pt8", seed, timeSinceFile, MagDependentAperiodicityOptions.ALL_PT8_VALUES);

		
		//	U3 BPT Simulations
//		long seed = 836529;
//		MagDependentAperiodicityOptions aper = MagDependentAperiodicityOptions.MID_VALUES;
//		String fileName="/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";
//		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fileName);
//		double numYrs=100000;
//		File parentDir = new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/ERF_Coordination/LongTermTD_2026/Analysis/bptSimulationsU3/");
//		bptSimulations(erf,numYrs,parentDir, "Run1_100k", seed, rootDir+"poissonSimulationsU3/Run1/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run2_100k", seed, rootDir+"poissonSimulationsU3/Run2/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run3_100k", seed, rootDir+"poissonSimulationsU3/Run3/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run4_100k", seed, rootDir+"poissonSimulationsU3/Run4/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run5_100k", seed, rootDir+"poissonSimulationsU3/Run5/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run6_100k", seed, rootDir+"poissonSimulationsU3/Run6/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run7_100k", seed, rootDir+"poissonSimulationsU3/Run7/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run8_100k", seed, rootDir+"poissonSimulationsU3/Run8/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run9_100k", seed, rootDir+"poissonSimulationsU3/Run9/outputTimesinceLast.txt", aper);
//		bptSimulations(erf,numYrs,parentDir, "Run10_100k", seed, rootDir+"poissonSimulationsU3/Run10/outputTimesinceLast.txt", aper);

		
//		// U3 Poisson Simulations
//		long seed = 984087634;
//		int numYrs = 1000000;
//		File parentDir = new File(rootDir+"poissonSimulationsU3/");
//		for(int i=1;i<=10;i++) {
//			String fileName="/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/UCERF3-TI/Figures/Fig11_FaultClusterFig/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";
//			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fileName);
//			seed += 7836271;
//			poissonSimulations(erf, parentDir, "Run"+i, seed, numYrs);
//		}
		
		
//		testOldVsNewSimulationMethod();
		
//		test_nthRupState();
		
//		makeTestTD_CalculationFiles();
		
//		generateDOLE_ReportPages();


	}

}
