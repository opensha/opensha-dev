package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.AbstractDiscretizedFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.calc.recurInterval.LognormalDistCalc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.Declustering.GardnerKnopoffDeclustering;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.math.Quantiles;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.utils.GardnerKnopoffAftershockFilter;
import scratch.kevin.ucerf3.eal.LossCOV_Model;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class U3ETAS_LossSimulationAnalysis {
	
	static final boolean D = true; // debug flag
	
	final static String dirName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Results";

	
//	static String fssFileName = "/Users/field/workspace/git/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
//	static String fssFileName = "/Users/field/FilesFromOldComputerTransfer/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
	static String fssFileName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";

	//	static String catalogsFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/Data/U3ETAS_SimulationData/results_m5_preserve_chain_FullTD_totRateScaleFactor1pt14.bin"; 
	//  static String catalogsFileName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnHazardAnalysisPaper/OtherStuff/Data/U3ETAS_SimulationData/results_m5_preserve_chain_FullTD_totRateScaleFactor1pt14.bin";
	static String catalogsFileName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/2024_06_06-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/results_m5_preserve_chain.bin";
	static String catalogsFileName2 = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/2024_08_23-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/results_m5_preserve_chain.bin";
	
	// mean loss data file
	static String lossDataCSV_Filename = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/average_nth_rup_losses_remapped_m2p5.csv";

	final static double[] CEA_PROB_ARRAY = {1e-3, 2e-3, 2.5e-3, 1d/350d, 4e-3, 5e-3, 1e-2, 2e-2, 4e-2, 0.1, 0.2};
	
//	Random random;
	JDKRandomGenerator randomGen;
	
	HashMap<Integer,Double> meanLossForNthRup, lossCOVForNthRup, magForNthRup;
	HashMap<Integer,Integer> fssIndexForNthRup, gridCellIndexForNthRup;
	FaultSystemSolutionERF_ETAS erf=null;

	static double totalPorfolioValue =  483e9/1e3; // total portfolio value in units of $1000
	static double catalogDuration=500; // hard coded
	static double catalogStartYear=2012;
	static long millisPerYr = (long)(1000*60*60*24*365.25);
	long startTimeMillis, endTimeMillis, fullDurationMillis;
	
	ArrayList<ObsEqkRupList> catalogList;
	ArrayList<ObsEqkRupList> catalogDeclusteredList;
	ArrayList<ObsEqkRupList> catalogRandomizedList;
	static double minMag = 5.05;
	static double deltaMag = 0.1;
	static int numMag = 40;
	
	// following are units of 1000 dollars
//	static double lossCurveLogMin = 0; // 1,000
//	static double lossCurveLogMax = 9;  // 1 trillion
	static double lossCurveLogMin = 3; // 1,000
	static double lossCurveLogMax = 12;  // 1 trillion
	static int lossCurveNum = 91; //40;
	static double lossCurveDelta = (lossCurveLogMax-lossCurveLogMin)/(double)(lossCurveNum-1);
	ArbitrarilyDiscretizedFunc curveXvals = getBlankLossExceedCurve(0d);
	
	LossCOV_Model covModel = LossCOV_Model.PORTER_POWER_LAW_2020_09_01_DOLLARS;
	double[] randomLossForEventID;
	int totalNumEvents;
	
	public U3ETAS_LossSimulationAnalysis(int seed) {
		
		// Note the grid cell stated for a catalog rupture does not necessarily equal that 
		// implied by the hypocenter (known leakage problem due to randomness added in U3ETAS)
		
		this.randomGen = new JDKRandomGenerator(seed);

		
		File outputDir = new File(dirName);
		if(!outputDir.exists()) outputDir.mkdir();
		
		TimeSpan tsp = new TimeSpan(TimeSpan.YEARS,TimeSpan.YEARS);
		tsp.setStartTime((int)catalogStartYear);
		tsp.setDuration((int)catalogDuration);
		startTimeMillis = tsp.getStartTimeInMillis();
		endTimeMillis = tsp.getEndTimeCalendar().getTimeInMillis();
		fullDurationMillis = endTimeMillis-startTimeMillis;
		
		erf = getERF();
		
		// Load catalogs & loss data
		CSVFile<String> lossDataCSV_File=null;
		File[] catFileArray = {new File(catalogsFileName),new File(catalogsFileName2)};
//		File[] catFileArray = {new File(catalogsFileName)};
		try {
			catalogList = loadCatalogs(catFileArray, 5.0, false);
			lossDataCSV_File = CSVFile.readFile(new File(lossDataCSV_Filename), true);
		} catch (IOException | DocumentException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		if (D) {
			System.out.println("num Catalogs = "+catalogList.size());
			System.out.println("Loss Data:");
			System.out.println("\tNumCols: "+lossDataCSV_File.getNumCols());
			System.out.println("\tNumRows: "+lossDataCSV_File.getNumRows());			

		}
		
		// set HashMaps from loss data
		meanLossForNthRup = new HashMap<Integer,Double>();
		lossCOVForNthRup = new HashMap<Integer,Double>();
		magForNthRup = new HashMap<Integer,Double>();
		fssIndexForNthRup = new HashMap<Integer,Integer>();
		gridCellIndexForNthRup = new HashMap<Integer,Integer>();
		for(int r=1;r<lossDataCSV_File.getNumRows();r++) {
				int nthIndex = lossDataCSV_File.getInt(r, 0);
				fssIndexForNthRup.put(nthIndex, lossDataCSV_File.getInt(r, 1));
				gridCellIndexForNthRup.put(nthIndex, lossDataCSV_File.getInt(r, 2));
				magForNthRup.put(nthIndex, lossDataCSV_File.getDouble(r, 3));
				meanLossForNthRup.put(nthIndex, lossDataCSV_File.getDouble(r, 4)*1000d);
				lossCOVForNthRup.put(nthIndex, lossDataCSV_File.getDouble(r, 5));
		}
		
		// check min and max COVs
		double minCOVfile=Double.MAX_VALUE, minCOVcalc=Double.MAX_VALUE;
		double maxCOVfile=0, maxCOVcalc=0;
		for(int n: meanLossForNthRup.keySet()) {
			if(minCOVfile>lossCOVForNthRup.get(n))  minCOVfile=lossCOVForNthRup.get(n);
			if(maxCOVfile<lossCOVForNthRup.get(n))  maxCOVfile=lossCOVForNthRup.get(n);
			double covCalc = covModel.getCOV(meanLossForNthRup.get(n));
			if(minCOVcalc>covCalc)  minCOVcalc=covCalc;
			if(maxCOVcalc<covCalc)  maxCOVcalc=covCalc;
		}
		System.out.println("minCOVfile="+minCOVfile);
		System.out.println("maxCOVfile="+maxCOVfile);
		System.out.println("minCOVcalc="+minCOVcalc);
		System.out.println("maxCOVcalc="+maxCOVcalc);
		
	
		// check mags & indices against catalogs
		double meanCatLengthYears=0;
		for (ObsEqkRupList catalog : catalogList) {
			meanCatLengthYears += ((double)(catalog.get(catalog.size()-1).getOriginTime()-catalog.get(0).getOriginTime()))/(double)millisPerYr;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				int nth = rup.getNthERF_Index();
				double magDiff=Math.abs(rup.getMag()/magForNthRup.get(nth) - 1.0);
				if(magDiff > 1e-4)
					throw new RuntimeException("ERROR: mags differ for nth rup: "+nth+"; cat mag = "+rup.getMag()+"; loss-file mag = "+ magForNthRup.get(nth));		
				if(rup.getFSSIndex() != fssIndexForNthRup.get(nth))
					throw new RuntimeException("ERROR: fss indices differ for nth rup: "+nth+"; cat fssID = "+rup.getFSSIndex()+"; loss-file fssID = "+ fssIndexForNthRup.get(nth));		
				if(rup.getGridNodeIndex() != gridCellIndexForNthRup.get(nth))
					throw new RuntimeException("ERROR: grid cell indices differ for nth rup: "+nth+"; cat cellID = "+rup.getGridNodeIndex()+"; loss-file cellID = "+gridCellIndexForNthRup.get(nth));
			}
		}		
		
		if(D) System.out.println("ave catalog duration (yrs): "+(float)(meanCatLengthYears/catalogList.size()));

		catalogDeclusteredList = getGK_DeclusteredCatalog(catalogList);
		
		getRandomizedCatalogs();
		
		setRandomLossForEvents();
		
		
	}
	
	public void plotMFDs() {
		// MFDs from ERFs:
		IncrementalMagFreqDist catalogMFD = makeCatalogMFD(catalogList);
		catalogMFD.setName("catalogMFD");
		
		IncrementalMagFreqDist catalogGKdeclMFD = makeCatalogMFD(catalogDeclusteredList);
		catalogGKdeclMFD.setName("catalogGKdeclMFD");
		
		IncrementalMagFreqDist catalogGKfiltMFD = makeCatalogMFD(getU3_GK_FilteredCatalog(catalogList));
		catalogGKfiltMFD.setName("catalogGKfiltMFD");
		
		IncrementalMagFreqDist catalogSpontaneustMFD = makeCatalogMFD(getSpontaneousEventsCatalog(catalogList));
		catalogSpontaneustMFD.setName("catalogSpontaneustMFD");
		
		if(erf==null)erf = getERF();
		
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
		erf.updateForecast();
		SummedMagFreqDist erfMFD_TD = org.opensha.sha.earthquake.calc.ERF_Calculator.getTotalMFD_ForERF(erf, minMag,minMag+(numMag-1)*deltaMag,numMag, true);
		erfMFD_TD.setName("erfMFD_TD");
		
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();
		SummedMagFreqDist erfMFD_TimeInd = org.opensha.sha.earthquake.calc.ERF_Calculator.getTotalMFD_ForERF(erf, minMag,minMag+(numMag-1)*deltaMag,numMag, true);
		erfMFD_TimeInd.setName("erfMFD_TimeInd");
		
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		mfdList.add(catalogMFD);
		mfdList.add(catalogGKdeclMFD);
		mfdList.add(catalogGKfiltMFD);
		mfdList.add(catalogSpontaneustMFD);
//		mfdList.add(erfMFD_TD);
//		mfdList.add(erfMFD_TimeInd);

		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.ORANGE));
		File outputDir = new File(dirName+ "/MFD_Figs");
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/MFD_Figs/incrMFD";
		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, "", "Magnitude", "Incremental Rate (per year)", 
				new Range(5,8.5), new Range(1e-5,2), false, true, fileNamePrefix, true);

		ArrayList<XY_DataSet> mfdCumList = new ArrayList<XY_DataSet>();
		mfdCumList.add(catalogMFD.getCumRateDistWithOffset());
		mfdCumList.add(catalogGKdeclMFD.getCumRateDistWithOffset());
		mfdCumList.add(catalogGKfiltMFD.getCumRateDistWithOffset());
		mfdCumList.add(catalogSpontaneustMFD.getCumRateDistWithOffset());
//		mfdCumList.add(erfMFD_TD.getCumRateDistWithOffset());
//		mfdCumList.add(erfMFD_TimeInd.getCumRateDistWithOffset());
		
		fileNamePrefix = dirName+"/MFD_Figs/cumMFD";
		PlottingUtils.writeAndOrPlotFuncs(mfdCumList, plotChars, "", "Magnitude", "Cumulative Rate (per year)", 
				new Range(5,8.5), new Range(1e-5,10), false, true, fileNamePrefix, true);
		
		IncrementalMagFreqDist catalogMFD_Ratio = makeMFD();
		catalogMFD_Ratio.setName("catalogMFD_Ratio");
		IncrementalMagFreqDist mfdGKdeclMFD_Ratio = makeMFD();
		mfdGKdeclMFD_Ratio.setName("mfdGKdeclMFD_Ratio");
		IncrementalMagFreqDist mfdGKfiltMFD_Ratio = makeMFD();
		mfdGKfiltMFD_Ratio.setName("mfdGKfiltMFD_Ratio");
		IncrementalMagFreqDist mfdSpontaneousMFD_Ratio = makeMFD();
		mfdSpontaneousMFD_Ratio.setName("mfdSpontaneousMFD_Ratio");
		for(int i=0;i<catalogMFD.size();i++) {
			catalogMFD_Ratio.set(i,catalogMFD.getY(i)/catalogMFD.getY(i));
			mfdGKdeclMFD_Ratio.set(i,catalogGKdeclMFD.getY(i)/catalogMFD.getY(i));
			mfdGKfiltMFD_Ratio.set(i,catalogGKfiltMFD.getY(i)/catalogMFD.getY(i));
			mfdSpontaneousMFD_Ratio.set(i,catalogSpontaneustMFD.getY(i)/catalogMFD.getY(i));
		}
		ArrayList<XY_DataSet> ratioList = new ArrayList<XY_DataSet>();
		ratioList.add(catalogMFD_Ratio);
		ratioList.add(mfdGKdeclMFD_Ratio);
		ratioList.add(mfdGKfiltMFD_Ratio);
		ratioList.add(mfdSpontaneousMFD_Ratio);
		fileNamePrefix = dirName+"/MFD_Figs/mfdRatios";
		PlottingUtils.writeAndOrPlotFuncs(ratioList, plotChars, "", "Magnitude", "Incremental Ratio", 
				new Range(5,8.5), new Range(0,1), false, false, fileNamePrefix, true);
	}
	
	/**
	 * Units: Millions of dollars
	 */
	public void writeAAL_ValuesBillions() {
		System.out.println("\nAverage Annual losses (AAL; billion $):");
		
		double aal_Catalog = getAveAnnaulLossFromCatalogList(catalogList)/1e9d;
		System.out.println("\nCatalog Stats (based on mean loss per rupture):");
		double[] stats = getLossStatsFromCatalogList(catalogList, false, true);
		System.out.println("\tAAL from catalogs: "+(float)(stats[0]/1e9d)+"\t(AAL=EventRate*MeanLoss)");
		System.out.println("\t\tEventRate = : "+(float)(stats[1]));
		System.out.println("\t\tMeanLoss = : "+(float)(stats[2]/1e9d));
		System.out.println("\t\tMedianLoss = : "+(float)(stats[3]/1e9d));
		System.out.println("\t\tLossCOV = : "+(float)(stats[4]));
		System.out.println("\t\tMinNonZeroLoss = : "+(float)(stats[5]/1e9d));
		System.out.println("\t\tMaxLoss = : "+(float)(stats[6]/1e9d));
		System.out.println("\t\tfractZeroLossRups = : "+(float)(stats[7]));
		System.out.println("\t\tAAL from faults = : "+(float)(stats[8]/1e9d));
		System.out.println("\t\tAAL from gridded seismicity = : "+(float)(stats[9]/1e9d));
		System.out.println("\t\tAAL for spontaneous seismicity = : "+(float)(stats[10]/1e9d));

		// NOT sure why this is so painfully slow; main results are equivalent to above
//		System.out.println("\nCatalog Stats (based on random loss from COV):");
//		stats = getLossStatsFromCatalogList(catalogList, true, true);
//		System.out.println("\tAAL from catalogs: "+(float)(stats[0]/1e9d)+"\t(AAL=EventRate*MeanLoss)");
//		System.out.println("\t\tEventRate = : "+(float)(stats[1]));
//		System.out.println("\t\tMeanLoss = : "+(float)(stats[2]/1e9d));
//		System.out.println("\t\tMedianLoss = : "+(float)(stats[3]/1e9d));
//		System.out.println("\t\tLossCOV = : "+(float)(stats[4]));
//		System.out.println("\t\tMinNonZeroLoss = : "+(float)(stats[5]/1e9d));
//		System.out.println("\t\tMaxLoss = : "+(float)(stats[6]/1e9d));
//		System.out.println("\t\tfractZeroLossRups = : "+(float)(stats[7]));
//		System.out.println("\t\tAAL from faults = : "+(float)(stats[8]/1e9d));
//		System.out.println("\t\tAAL from gridded seismicity = : "+(float)(stats[9]/1e9d));


		double aal_TD = getAveAnnualLossFromERF(true)/1e9d;  // Time Dependent
		if(D) System.out.println("\tAAL from ERF-TD = "+(float)aal_TD); 
		
		double aal_TimeInd = getAveAnnualLossFromERF(false)/1e9d;  // Time Independent
		if(D) System.out.println("\tAAL from ERF-TimeInd = "+(float)aal_TimeInd); 
		
		double aal_GKdeclCatalog = getAveAnnaulLossFromCatalogList(catalogDeclusteredList)/1e9d;
		System.out.println("\tAAL from GK declustered catalogs: "+(float)aal_GKdeclCatalog+"\n");
		
		double aal_GKfilteredCatalog = getAveAnnaulLossFromCatalogList(getU3_GK_FilteredCatalog(catalogList))/1e9d;
		System.out.println("\tAAL from GK U3 filtered catalogs: "+(float)aal_GKfilteredCatalog+"\n");


	}
	
	
	/**
	 * This plots the distribution of rupture losses (mean and 68% and 95% fractiles) as a function 
	 * of magnitude (not including rup rates).
	 * 
	 * @param catList
	 * @param mkPlot - not yet working
	 */
	public void plotRupLossVsMagStats(ArrayList<ObsEqkRupList> catList, boolean randomFromCOV, boolean mkPlot) {
		
		IncrementalMagFreqDist mfd = makeMFD();
		ArbDiscrEmpiricalDistFunc_3D arbFunc3D = new ArbDiscrEmpiricalDistFunc_3D(mfd.getMinX(),mfd.getMaxX(),mfd.size());
		for (ObsEqkRupList catalog : catList) {
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double rupLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(randomFromCOV) {
					if(rupLoss==0) {
						arbFunc3D.set(mfd.getClosestXIndex(rup.getMag()), rupLoss, 1.0);
					}
					else {
						double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(rupLoss).sample(); // get random sample
						arbFunc3D.set(mfd.getClosestXIndex(rup.getMag()), randLoss, 1.0);
					}
				}
				else
					arbFunc3D.set(mfd.getClosestXIndex(rup.getMag()), rupLoss, 1.0);
			}
		}				
		
//		ArbDiscrEmpiricalDistFunc[] array = arbFunc3D.getArbDiscrEmpDistFuncArray();
//		for(int i=0;i<array.length;i++) {
//			System.out.println(array[i].size());
//		}

		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
		funcList.add(arbFunc3D.getMeanCurve());
		funcList.add(arbFunc3D.getInterpolatedFractileCurve(0.5));
		funcList.add(arbFunc3D.getInterpolatedFractileCurve(0.025));
		funcList.add(arbFunc3D.getInterpolatedFractileCurve(0.975));
		funcList.add(arbFunc3D.getInterpolatedFractileCurve(0.16));
		funcList.add(arbFunc3D.getInterpolatedFractileCurve(0.84));
		funcList.get(0).setName("Mean");
		funcList.get(1).setName("Median");
		funcList.get(2).setName("2.5% fractile");
		funcList.get(3).setName("97.5% fractile");
		funcList.get(4).setName("16% fractile");
		funcList.get(5).setName("84% fractile");
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		String fileNamePrefix = dirName+"/LossForAllRupStats";
		if(randomFromCOV)
			fileNamePrefix += "RandCOV";
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "", "Magnitude", "Loss ($)",
				new Range(5d,8.5), new Range(1e1,1e11), false, true, 6d, 5d, fileNamePrefix, true);
	}
	
	/**
	 * This computes various stats from catalogs assuming all have duration = this.catalogDuration
	 * Output:
	 * stats[0] = AAL (average annual loss); equals EventRate*MeanLoss
	 * stats[1] = Event Rate (per year)
	 * stats[2] = Mean Loss
	 * stats[3] = Median Loss
	 * stats[4] = Loss COV
	 * stats[5] = minLoss; // minimum non-zero loss
	 * stats[6] = maxLoss;
	 * stats[7] = fractZeroLosses; // fraction of ruptures that are zero-loss 
	 * stats[8] = AAL from faults 
	 * stats[9] = AAL from gridded seismicity 
	 * @param catList
	 * @return stats
	 */
	public double[] getLossStatsFromCatalogList(ArrayList<ObsEqkRupList> catList, boolean randFromCOV, boolean mkPlot) {
		HistogramFunction hist = getBlankLossExceedCurveLogX();
		double[] stats = new double[11];
		ArbDiscrEmpiricalDistFunc distFunc = new ArbDiscrEmpiricalDistFunc();
		double aal=0,aal_faults=0,aal_gridded=0, aal_spontaneous=0;
		int n=0;
int cat =0;
		int numZeroLosses = 0;
		double minLoss = Double.MAX_VALUE, maxLoss=0;
		for (ObsEqkRupList catalog : catList) {
//cat +=1;
//if(randFromCOV) System.out.println("cat "+cat);
			for (ObsEqkRupture obsRup : catalog) {
				n+=1;
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double rupLoss;
				if(randFromCOV)
					rupLoss = randomLossForEventID[rup.getID()];
				else
					rupLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				
				if(rupLoss==0.0)
					numZeroLosses+=1;			
				else
					if(minLoss>rupLoss) minLoss=rupLoss;
				
				if(maxLoss<rupLoss) 
					maxLoss=rupLoss;
				aal += rupLoss;
				if(rup.getFSSIndex()>0)
					aal_faults +=rupLoss;
				else
					aal_gridded +=rupLoss;
				
				if(rup.getGeneration()==0)
					aal_spontaneous+=rupLoss;

				distFunc.set(rupLoss, 1.0);
				
				if(mkPlot) {
					double log10Loss = Math.log10(rupLoss);
					if(log10Loss<lossCurveLogMin)
						hist.add(0,1.0); // put in first bin
					else
						hist.add(log10Loss, 1.0);
				}
			}
		}		
		stats[0] = aal/((double)catList.size()*catalogDuration); 	// AAL
		stats[1] = n/((double)catList.size()*catalogDuration);		// Event Rate
		stats[2] = distFunc.getMean();
		stats[3] = distFunc.getMedian();
		stats[4] = distFunc.getCOV();
		stats[5] = minLoss;
		stats[6] = maxLoss;
		stats[7] = numZeroLosses/(double)n;
		stats[8] = aal_faults/((double)catList.size()*catalogDuration);
		stats[9] = aal_gridded/((double)catList.size()*catalogDuration);
		stats[10] = aal_spontaneous/((double)catList.size()*catalogDuration);
		
		if(Math.abs(stats[1]*stats[2]/stats[0] - 1.0) > 1e-6)
			throw new RuntimeException("AAL != EventRate*MeanLoss");
		
		if(mkPlot) {
			hist.scale(1.0/(catList.size()*catalogDuration));
//			// this should about equal the aal; it does
//			double aalTest =0;
//			for(int i=1;i<hist.size();i++)
//				aalTest += Math.pow(10, hist.getX(i))*hist.getY(i);
//			System.out.println("aalTest="+aalTest/1e9); // convert to billions
			
			ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
			funcList.add(hist);
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLUE));
			if(randFromCOV)
				PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Catalog Loss Rate Dist; Random From COV", 
					"Log10(Loss;thousand $)", "Bin Rate (per year)", null, null, false, false, null, true);
			else
				PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Catalog Loss Rate Dist; Mean Loss", 
						"Log10(Loss;thousand $)", "Bin Rate (per year)", null, null, false, false, null, true);

		}
		return stats;
	}

	
	/**
	 * This plots all points on the distribution (using ArbDiscrEmpiricalDistFunc)
	 * @param randFromCOV
	 * @param mkPlot
	 */
	public void plotCumAggrYearlyLossDistributions(boolean randFromCOV, boolean mkPlot) {
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(getCatalogs(), 1.0);
		ArbDiscrEmpiricalDistFunc distFuncTD = new ArbDiscrEmpiricalDistFunc();
		int counter = -1;
		double writeAt = 0.1;
		for (ObsEqkRupList catalog : subCatList) {
			counter+=1;
			double percDone = (double)counter/(double)subCatList.size();
			if(percDone > writeAt) {
				System.out.println("percDone: "+(float)percDone);
				writeAt+=0.1;
			}
			double totLoss=0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				if(randFromCOV)
					totLoss += randomLossForEventID[rup.getID()];
				else
					totLoss += meanLossForNthRup.get(rup.getNthERF_Index());
			}
			distFuncTD.set(totLoss, 1.0);	
		}			
		DiscretizedFunc cumDistFuncTD = distFuncTD.getCumDist();
		cumDistFuncTD.scale(1.0/distFuncTD.calcSumOfY_Vals());
		cumDistFuncTD.setName("cumDistFuncTD");
//		for(int i=0;i<cumDistFuncTD.size();i++)
//			cumDistFuncTD.set(cumDistFuncTD.getX(i), 1.0-cumDistFuncTD.getY(i));
		
		subCatList = getSubcatalogList(getRandomizedCatalogs(), 1.0);
		ArbDiscrEmpiricalDistFunc distFuncTI = new ArbDiscrEmpiricalDistFunc();
		counter = -1;
		writeAt = 0.1;
		for (ObsEqkRupList catalog : subCatList) {
			counter+=1;
			double percDone = (double)counter/(double)subCatList.size();
			if(percDone > writeAt) {
				System.out.println("percDone: "+(float)percDone);
				writeAt+=0.1;
			}
			double totLoss=0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				if(randFromCOV)
					totLoss += randomLossForEventID[rup.getID()];
				else
					totLoss += meanLossForNthRup.get(rup.getNthERF_Index());
			}
			distFuncTI.set(totLoss, 1.0);	
		}			
		DiscretizedFunc cumDistFuncTI = distFuncTI.getCumDist();
		cumDistFuncTI.scale(1.0/distFuncTI.calcSumOfY_Vals());
		cumDistFuncTI.setName("cumDistFuncTI");
//		for(int i=0;i<cumDistFuncTI.size();i++)
//			cumDistFuncTI.set(cumDistFuncTI.getX(i), 1.0-cumDistFuncTI.getY(i));

		double[] vals = {10e6,5e9};
		System.out.print("X\tY_TD\tY_TI");
		for(double val:vals)
			System.out.println(val+"\t"+cumDistFuncTD.getInterpolatedY_inLogXLogYDomain(val)+"\t"+cumDistFuncTI.getInterpolatedY_inLogXLogYDomain(val));

		
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
		funcList.add(cumDistFuncTD);
		funcList.add(cumDistFuncTI);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		if(randFromCOV)
			PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Cumulative Loss Dist; Random From COV", 
					"Loss($)", "Cumulative Dist", null, null, false, false, null, true);
			else
				PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Cumulative Loss Dist; Mean Loss", 
						"Loss ($)", "Cumulative Dist", null, null, true, false, null, true);
	}

	
	
	/**
	 * AAL in Billions per year
	 * This assumes catalogs have duration = this.catalogDuration
	 * @param catList
	 * @return
	 */
	public double getAveAnnaulLossFromCatalogList(ArrayList<ObsEqkRupList> catList) {

//		// commented out stuff was looking at empirical distributions
//		ArbDiscrEmpiricalDistFunc distFunc = new ArbDiscrEmpiricalDistFunc();
//
//		int numPts = 10000000;
//		double maxX = Math.pow(10,lossCurveLogMax);
//		double delta = maxX/(double)numPts;
//		double minX = 0.5*delta;
//		HistogramFunction histFunc = new HistogramFunction(minX, numPts, delta);

		double aal =0;
		int n=0;
		for (ObsEqkRupList catalog : catList) {
			for (ObsEqkRupture obsRup : catalog) {
				n+=1;
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double rupLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				aal += rupLoss;
//				distFunc.set(rupLoss, 1.0);
//				histFunc.add(rupLoss, 1.0);
			}
		}		
		aal /= (double)catList.size()*catalogDuration;
//		double meanRateEvents = n/((double)catList.size()*catalogDuration);
//		System.out.println("HERE: mean = "+distFunc.getMean()+"; median = "+
//				distFunc.getMedian()+"; mode = "+distFunc.getMostCentralMode()+"; COV = "+distFunc.getCOV()
//				+"; AAL = "+distFunc.getMean()*meanRateEvents);
//		histFunc.normalizeBySumOfY_Vals();
////		double median = histFunc.getCumulativeDistFunctionWithHalfBinOffset().getFirstInterpolatedX(0.5);
////		double median = histFunc.getCumulativeDistFunctionWithHalfBinOffset().getFirstInterpolatedX_inLogYDomain(0.5);
//		double median = histFunc.getCumulativeDistFunctionWithHalfBinOffset().getFirstInterpolatedX_inLogXLogYDomain(0.5);
//		System.out.println("HERE2: mean = "+histFunc.computeMean()+"; median = "+median+
//		"; mode = "+histFunc.getMode()+"; AAL = "+histFunc.computeMean()*meanRateEvents);
//		System.out.println("HERE2: first X bin: "+(histFunc.getX(0)-histFunc.getDelta()/2d)+" to "+(histFunc.getX(0)+histFunc.getDelta()/2d));
		return aal;
	}
	
	
	
	public double getAveAnnualLossFromERF(boolean timeDep) {
		double aal =0;
		if(erf==null)
			erf = getERF();
		
		if(timeDep)
			erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_BPT);
		else
			erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();

//		if(D)System.out.println("Prob Model: "+erf.getParameter(ProbabilityModelParam.NAME).getValue());
		
		double duration = erf.getTimeSpan().getDuration();
		double maxMagRejected = 0;
		for(int n=0;n<erf.getTotNumRups();n++) {
			ProbEqkRupture rup = erf.getNthRupture(n);
			if(meanLossForNthRup.containsKey(n)) {
				aal += rup.getMeanAnnualRate(duration)*meanLossForNthRup.get(n);
			}
			else {
				if(maxMagRejected<rup.getMag())
					maxMagRejected=rup.getMag();
			}
		}
		
//		if(D) System.out.println("\tAAL from ERF = "+(float)aal); // +";   maxMagRejected="+(float)maxMagRejected);

		return aal;
	}
	
	public ArbitrarilyDiscretizedFunc getLossExceedProbCurveFromERF() {
		ArbitrarilyDiscretizedFunc lossExceedProbCurve = getBlankLossExceedCurve(1d);
		if(erf==null)
			erf = getERF();
		erf.getTimeSpan().setDuration(1d);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();
		double duration = erf.getTimeSpan().getDuration();
		for(int n=0;n<erf.getTotNumRups();n++) {
			if(meanLossForNthRup.containsKey(n)) {
				double rupProb = erf.getNthRupture(n).getProbability();
				double aveLoss = meanLossForNthRup.get(n);
				if(aveLoss == 0.0)
					continue;
				DiscretizedFunc condCurve = covModel.calcLossExceedanceProbs(lossExceedProbCurve, aveLoss);
				for(int i=0;i<condCurve.size();i++) {
					double prevVal=lossExceedProbCurve.getY(i);
					double newVal = prevVal*(1.0-condCurve.getY(i)*rupProb);
					lossExceedProbCurve.set(i,newVal);
				}
			}
		}
		for(int i=0;i<lossExceedProbCurve.size();i++) {
			double newVal = 1.0-lossExceedProbCurve.getY(i);
			lossExceedProbCurve.set(i,newVal);
		}
		return lossExceedProbCurve;
	}
	
	
	/**
	 * This tests the covModel changes made in Aug 2024 using
	 * the loss exceed curve from a Poisson ERF.
	 * @return
	 */
	public ArrayList<XY_DataSet> testCOV_ModelChangesInAug2024() {
		LossCOV_Model covModel1 = LossCOV_Model.PORTER_POWER_LAW_2020_09_01_fixed;
		ArbitrarilyDiscretizedFunc lossExceedProbCurve = getBlankLossExceedCurve(1d);  // 
		LossCOV_Model covModel2 = LossCOV_Model.PORTER_POWER_LAW_2020_09_01;
		ArbitrarilyDiscretizedFunc lossExceedProbCurve2 = getBlankLossExceedCurve(1d); // previous
		ArbitrarilyDiscretizedFunc lossExceedProbCurve3 = getBlankLossExceedCurve(1d); // old/wrong distribution
		
//		double l = 1e7;
//		System.out.println("HERE1: "+covModel1.getCOV(l));
//		System.out.println("HERE2: "+covModel2.getCOV(l));
//		System.out.println("HERE3: "+covModel3.getCOV(l*1000/totalPorfolioValue));
		
		if(erf==null)
			erf = getERF();
		erf.getTimeSpan().setDuration(1d);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();
		double duration = erf.getTimeSpan().getDuration();
		for(int n=0;n<erf.getTotNumRups();n++) {
			if(meanLossForNthRup.containsKey(n)) {
				double rupProb = erf.getNthRupture(n).getProbability();
				double aveLoss = meanLossForNthRup.get(n);
				if(aveLoss == 0.0)
					continue;
				DiscretizedFunc condCurve = covModel1.calcLossExceedanceProbs(lossExceedProbCurve, aveLoss);
				DiscretizedFunc condCurve2 = covModel2.calcLossExceedanceProbs(lossExceedProbCurve2, aveLoss);
				
				// old/wrong distibution:
				DiscretizedFunc condCurve3 = new ArbitrarilyDiscretizedFunc();
				// this is the old/wrong distribution
				LogNormalDistribution dist = new LogNormalDistribution(Math.log(aveLoss), covModel1.getCOV(aveLoss));
				for (int i=0; i<lossExceedProbCurve3.size(); i++) {
					double x = lossExceedProbCurve3.getX(i);
					double y = 1d - dist.cumulativeProbability(x);
					condCurve3.set(x, y);
				}

				
				for(int i=0;i<condCurve.size();i++) {
					double prevVal=lossExceedProbCurve.getY(i);
					double newVal = prevVal*(1.0-condCurve.getY(i)*rupProb);
					lossExceedProbCurve.set(i,newVal);
					double prevVal2=lossExceedProbCurve2.getY(i);
					double newVal2 = prevVal2*(1.0-condCurve2.getY(i)*rupProb);
					lossExceedProbCurve2.set(i,newVal2);
					double prevVal3=lossExceedProbCurve3.getY(i);
					double newVal3 = prevVal3*(1.0-condCurve3.getY(i)*rupProb);
					lossExceedProbCurve3.set(i,newVal3);
				}
			}
		}
		for(int i=0;i<lossExceedProbCurve.size();i++) {
			double newVal = 1.0-lossExceedProbCurve.getY(i);
			lossExceedProbCurve.set(i,newVal);
			double newVal2 = 1.0-lossExceedProbCurve2.getY(i);
			lossExceedProbCurve2.set(i,newVal2);
			double newVal3 = 1.0-lossExceedProbCurve3.getY(i);
			lossExceedProbCurve3.set(i,newVal3);
		}
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
		funcList.add(lossExceedProbCurve);
		funcList.add(lossExceedProbCurve2);
		funcList.add(lossExceedProbCurve3);
		return funcList;
	}


	
	
	public ArbitrarilyDiscretizedFunc getLossRateCurveFromERF() {
		ArbitrarilyDiscretizedFunc lossRateCurve = getBlankLossExceedCurve(0d);
		if(erf==null)
			erf = getERF();
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();
		double duration = erf.getTimeSpan().getDuration();
		for(int n=0;n<erf.getTotNumRups();n++) {
			if(meanLossForNthRup.containsKey(n)) {
				double rupRate = erf.getNthRupture(n).getMeanAnnualRate(duration);
				double loss = meanLossForNthRup.get(n);
				if(loss == 0.0) {
					continue;
				}
				DiscretizedFunc condCurve = covModel.calcLossExceedanceProbs(lossRateCurve, loss);
				for(int i=0;i<condCurve.size();i++) {
					lossRateCurve.set(i,lossRateCurve.getY(i)+condCurve.getY(i)*rupRate);
				}
			}
		}
		return lossRateCurve;
	}
	
	
	
	
	public ArbitrarilyDiscretizedFunc getLossRateCurveFromERF_zeroCOV() {
		ArbitrarilyDiscretizedFunc lossRateCurve = getBlankLossExceedCurve(0d);
		if(erf==null)
			erf = getERF();
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		erf.updateForecast();
		double duration = erf.getTimeSpan().getDuration();
		for(int n=0;n<erf.getTotNumRups();n++) {
			if(meanLossForNthRup.containsKey(n)) {
				double loss = meanLossForNthRup.get(n);
				if(loss == 0.0)
					continue;
				double rupRate = erf.getNthRupture(n).getMeanAnnualRate(duration);
				for(int i=0;i<lossRateCurve.size();i++) {
					if(lossRateCurve.getX(i)<loss) {
						lossRateCurve.set(i,lossRateCurve.getY(i)+rupRate);
					}
					else
						break;
				}
			}
		}
		return lossRateCurve;
	}
	
	public ArbitrarilyDiscretizedFunc getLossRateCurveFromCatalogs(ArrayList<ObsEqkRupList> catList) {
		if(D) System.out.println("Working on getLossRateCurveFromCatalogs(); num catalogs = "+catList.size());
		ArbitrarilyDiscretizedFunc lossRateCurve = getBlankLossExceedCurve(0d);
		for (ObsEqkRupList catalog : catList) {
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss == 0.0)
					continue;
				DiscretizedFunc condCurve = covModel.calcLossExceedanceProbs(lossRateCurve, aveLoss);
				for(int i=0;i<condCurve.size();i++) {
					lossRateCurve.set(i,lossRateCurve.getY(i)+condCurve.getY(i));
				}
			}
		}		
		lossRateCurve.scale(1.0/((double)catList.size()*catalogDuration));
		return lossRateCurve;
	}
	
	
	public ArbitrarilyDiscretizedFunc getLossRateCurveFromCatalogsZeroCOV(ArrayList<ObsEqkRupList> catList) {
		ArbitrarilyDiscretizedFunc lossRateCurve = getBlankLossExceedCurve(0d);
		for (ObsEqkRupList catalog : catList) {
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss == 0.0)
					continue;
				for(int i=0;i<lossRateCurve.size();i++) {
					if(lossRateCurve.getX(i)<aveLoss)
						lossRateCurve.set(i,lossRateCurve.getY(i)+1d);
					else
						break;
				}
			}
		}		
		lossRateCurve.scale(1.0/((double)catList.size()*catalogDuration));
		return lossRateCurve;
	}
	
	
	public ArbitrarilyDiscretizedFunc getLossRateCurveFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList) {
		ArbitrarilyDiscretizedFunc lossRateCurve = getBlankLossExceedCurve(0d);
		for (ObsEqkRupList catalog : catList) {
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double loss = randomLossForEventID[rup.getID()];
				if(loss == 0.0)
					continue;
				for(int i=0;i<lossRateCurve.size();i++) {
					if(lossRateCurve.getX(i)<loss)
						lossRateCurve.set(i,lossRateCurve.getY(i)+1d);
					else
						break;
				}
			}
		}		
		lossRateCurve.scale(1.0/((double)catList.size()*catalogDuration));
		return lossRateCurve;
	}



	
	public UncertainArbDiscFunc getLossExceedProbCurveFromCatalogs(ArrayList<ObsEqkRupList> catList, double duration) {
		
		if(D) System.out.println("Working on getLossExceedProbCurveFromCatalogs(); num catalogs = "+catList.size());

		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		for (ObsEqkRupList catalog : subCatList) {
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(1d);
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss == 0.0)
					continue;
				DiscretizedFunc condCurve = covModel.calcLossExceedanceProbs(lossProbCurveForSubCat, aveLoss);
				for(int i=0;i<condCurve.size();i++) {
					double prevVal=lossProbCurveForSubCat.getY(i);
					double newVal = prevVal*(1.0-condCurve.getY(i));
					lossProbCurveForSubCat.set(i,newVal);
				}
			}
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				lossProbCurveForSubCat.set(i,1.0-lossProbCurveForSubCat.getY(i));
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}		
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		return uncertFunc;
	}
	
	

	
	
	
	public ArbitrarilyDiscretizedFunc getLossExceedProbCurveFromCatalogsZeroCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		for (ObsEqkRupList catalog : subCatList) {
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(0d);
			double maxLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(maxLoss<aveLoss) maxLoss=aveLoss;
			}
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < maxLoss)
					lossProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}		
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		return lossProbCurve;
	}
	
	
	public UncertainArbDiscFunc getLossExceedProbCurveFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
//		// TEMP
//		double lossExceedThresh = 4.0e4;
//		int[] numExceedForCatArray = new int[subCatList.size()];

		// TEMP
		int c=-1;
		for (ObsEqkRupList catalog : subCatList) {
			c+=1;

			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(0d);
			double maxLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxLoss<randLoss) maxLoss=randLoss;

				// TEMP
//				if(randLoss>=lossExceedThresh)
//					numExceedForCatArray[c] += 1;
				
			}
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < maxLoss)
					lossProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
			
//			// TEMP
//			if(numExceedForCatArray[c]>250) {
//				for (ObsEqkRupture obsRup : subCatList.get(c-1)) { // previous catalog
//					ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
//					double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
//					if(aveLoss==0)
//						continue;
//					System.out.println(rup.getMag()+"\t"+rup.getGeneration()+"\tPrevious catalog"); //  +"\t"+rup.getOriginTimeCal());
//				}
//				for (ObsEqkRupture obsRup : catalog) { // this catalog
//					ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
//					double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
//					if(aveLoss==0)
//						continue;
//					System.out.println(rup.getMag()+"\t"+rup.getGeneration()); //  +"\t"+rup.getOriginTimeCal());
//				}
//				System.exit(0);
//			}

			
		}		
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		
		// TEMP
//		int maxNum = 0;
//		for(int num:numExceedForCatArray)
//			if(num>maxNum) maxNum=num;
////		System.out.println("maxNum: "+maxNum);
//		HistogramFunction hist = new HistogramFunction(0d,(double)maxNum,maxNum+1);
//		for(int num:numExceedForCatArray)
//			hist.add(num, 1d);
//	//	System.out.println("\nhist.calcSumOfY_Vals():\n"+hist.calcSumOfY_Vals());
//
//		hist.normalizeBySumOfY_Vals();
//		System.out.println("\nnumExceedPDF:\n"+hist);
		
		return uncertFunc;
	}
	
	
	/**
	 * This returns a histogram representing the fraction of catalogs that had 0, 1, 2, 3... exceedances of the given
	 * threshold loss value.
	 * @param catList
	 * @param duration
	 * @param lossExceedThresh
	 * @return
	 */
	public HistogramFunction getNumLossExceedDistFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration, double lossExceedThresh) {

		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		int[] numExceedForCatArray = new int[subCatList.size()];
		int c=-1;
		for (ObsEqkRupList catalog : subCatList) {
			c+=1;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(randLoss>=lossExceedThresh)
					numExceedForCatArray[c] += 1;
			}
			
//			// TEMP
//			if(numExceedForCatArray[c]>250) {
//				for (ObsEqkRupture obsRup : subCatList.get(c-1)) { // previous catalog
//					ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
//					double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
//					if(aveLoss==0)
//						continue;
//					System.out.println(rup.getMag()+"\t"+rup.getGeneration()+"\tPrevious catalog"); //  +"\t"+rup.getOriginTimeCal());
//				}
//				for (ObsEqkRupture obsRup : catalog) { // this catalog
//					ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
//					double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
//					if(aveLoss==0)
//						continue;
//					System.out.println(rup.getMag()+"\t"+rup.getGeneration()); //  +"\t"+rup.getOriginTimeCal());
//				}
//				System.exit(0);
//			}

			
		}		
		
		int maxNum = 0, catForMaxNum=-1;
		for(c=0;c<numExceedForCatArray.length;c++ ) {
			int num = numExceedForCatArray[c];
			if(num>maxNum) {
				maxNum=num;
				catForMaxNum=c;
			}
		}
//		System.out.println("maxNum: "+maxNum);
		HistogramFunction hist = new HistogramFunction(0d,(double)maxNum,maxNum+1);
		for(int num:numExceedForCatArray)
			hist.add(num, 1d);
	//	System.out.println("\nhist.calcSumOfY_Vals():\n"+hist.calcSumOfY_Vals());

		hist.normalizeBySumOfY_Vals();
		hist.setInfo("maxNum="+maxNum+", which occurred in subCat="+catForMaxNum);
//		System.out.println("\nnumExceedPDF: for loss="+lossExceedThresh+"\n"+hist);

		return hist;
	}

	
	
	public UncertainArbDiscFunc testGetLossExceedProbCurveFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		for (ObsEqkRupList catalog : subCatList) {

			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(1d);
			double maxLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxLoss<randLoss) maxLoss=randLoss;

				
			}
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < maxLoss)
					lossProbCurveForSubCat.set(i,0.0); // prob of not exceeding
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
			

			
		}		
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		// convert from non-exceet to exceed prob
		for(int i=0;i<lossProbCurve.size();i++) {
			lossProbCurve.set(i,1-lossProbCurve.getY(i));
		}

		
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		
		return uncertFunc;
	}


	
	
	public UncertainArbDiscFunc getLossIncrProbDistFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		// this has offset bin x values by half delta
		HistogramFunction lossHistLogX = getBlankIncrLossCurveLogX();
		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankIncrLossCurve(0);
		double fracZeroLossSubCats = 0;
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		for (ObsEqkRupList catalog : subCatList) {
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankIncrLossCurve(0);
			double totLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				totLoss+=randLoss;
				double log10_Loss = Math.log10(randLoss);
				if(log10_Loss<=lossCurveLogMin) // this is the edge of the first bin
					continue;
//if(randLoss>6.3095734E10)
//	System.out.println("LOSS: "+randLoss);
				int index = lossHistLogX.getClosestXIndex(log10_Loss);
				lossProbCurveForSubCat.set(index, 1); // this ends up 1.0 if one or more occurrences

			}
			if(totLoss==0.0) // this should be good enough to avoid numerical problems
				fracZeroLossSubCats+=1;
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}		
		fracZeroLossSubCats /= (double)subCatList.size();
		lossProbCurve.setInfo("fracZeroLossSubCats="+fracZeroLossSubCats);
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankIncrLossCurve(0);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankIncrLossCurve(0);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		
		return uncertFunc;
	}
	
	
	
	public UncertainArbDiscFunc getAggrLossIncrProbDistFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		// this has offset bin x values by half delta
		HistogramFunction lossHistLogX = getBlankIncrLossCurveLogX();
		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankIncrLossCurve(0);
		double fracZeroLossSubCats = 0;
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		for (ObsEqkRupList catalog : subCatList) {
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankIncrLossCurve(0);
			double totLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				totLoss+=randLoss;
			}
			double log10_Loss = Math.log10(totLoss);
			if(log10_Loss<=lossCurveLogMin) // this is the edge of the first bin
				continue;
			int index = lossHistLogX.getClosestXIndex(log10_Loss);
			lossProbCurve.set(index,lossProbCurve.getY(index)+1); // this ends up 1.0 if one or more occurrences
			if(totLoss==0.0) // this should be good enough to avoid numerical problems
				fracZeroLossSubCats+=1;
		}		
		fracZeroLossSubCats /= (double)subCatList.size();
		lossProbCurve.setInfo("fracZeroLossSubCats="+fracZeroLossSubCats);
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankIncrLossCurve(0);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankIncrLossCurve(0);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		
		return uncertFunc;
	}

	
	
	/**
	 * Make the curve assuming Poisson; this gives the same as above with a randomized catalog
	 * @param catList
	 * @param duration
	 * @return
	 */
	public ArbitrarilyDiscretizedFunc testLossIncrProbDistFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {

		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankIncrLossCurve(0d);
		HistogramFunction lossHistLogX = getBlankIncrLossCurveLogX();
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		for (ObsEqkRupList catalog : subCatList) {
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				double log10_Loss = Math.log10(randLoss);
				if(log10_Loss<=lossCurveLogMin)
					continue;
				lossHistLogX.add(log10_Loss, 1.0);
			}
		}		
		lossHistLogX.scale(1.0/((double)subCatList.size()*duration));
		
		for(int i=0;i<lossHistLogX.size();i++)
			lossProbCurve.set(i,1-Math.exp(-lossHistLogX.getY(i)));
		
		return lossProbCurve;
	}


	
	
	

	public UncertainArbDiscFunc getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(ArrayList<ObsEqkRupList> catList, double duration) {
		
		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		int subCatIndexForMaxLoss = -1;
		int subCatIndexForMaxAggLoss = -1;
		int subCatIndex = -1;
		double maxLoss =  -1;
		double maxAggLoss =  -1;
		double aveAggrLoss = 0;
		for (ObsEqkRupList catalog : subCatList) {
			subCatIndex += 1;
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(0d);
			double totLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get random sample
				totLoss+=randLoss;
				if(randLoss>maxLoss) {
					maxLoss = randLoss;
					subCatIndexForMaxLoss = subCatIndex;
				}
//				totLoss+=aveLoss;
			}
			if(totLoss>maxAggLoss) {
				maxAggLoss = totLoss;
				subCatIndexForMaxAggLoss = subCatIndex;
			}			
			
			aveAggrLoss+=totLoss;
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < totLoss)
					lossProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}	
		aveAggrLoss /= (double)subCatList.size();
		lossProbCurve.scale(1.0/((double)subCatList.size()));
		String info = "aveAggrLoss = "+aveAggrLoss+"\nmaxLoss, subCatIndexForMaxAggLoss: "+maxLoss+", "+subCatIndexForMaxLoss+
				"\nmaxAggLoss, subCatIndexForMaxAggLoss: "+maxAggLoss+", "+subCatIndexForMaxAggLoss;
		lossProbCurve.setInfo(info);
		if(D) System.out.println(info);
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), subCatList.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		return uncertFunc;
	}
	
	
	public UncertainArbDiscFunc getMagExceedProbCurveFromCatalogs(ArrayList<ObsEqkRupList> catList, double duration) {
		
		IncrementalMagFreqDist magExceedProbCurve = makeMFD();
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catList, duration);
		
		for (ObsEqkRupList catalog : subCatList) {
			IncrementalMagFreqDist magExceedProbCurveForSubCat = makeMFD();
			double maxMag = 0;
			for (ObsEqkRupture obsRup : catalog) {
				if(obsRup.getMag()>maxMag)
					maxMag = obsRup.getMag();
			}
			for(int i=0;i<magExceedProbCurveForSubCat.size();i++) {
				if(magExceedProbCurveForSubCat.getX(i) < maxMag)
					magExceedProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<magExceedProbCurve.size();i++) {
				magExceedProbCurve.set(i,magExceedProbCurve.getY(i)+magExceedProbCurveForSubCat.getY(i));
			}
		}	
		magExceedProbCurve.scale(1.0/((double)subCatList.size()));
		IncrementalMagFreqDist low95confCurve = makeMFD();
		IncrementalMagFreqDist high95confCurve = makeMFD();
		for(int i=0;i<magExceedProbCurve.size();i++) {
			double xVal = magExceedProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(magExceedProbCurve.getY(i), subCatList.size());
			if(confArray[0]>magExceedProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+magExceedProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+subCatList.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}

		UncertainArbDiscFunc uncertFunc = new UncertainArbDiscFunc(magExceedProbCurve,low95confCurve,high95confCurve);
		return uncertFunc;
	}



	
	public static ArrayList<ObsEqkRupList> getGK_DeclusteredCatalog(ArrayList<ObsEqkRupList> catalogList) {
		
		if(D) System.out.println("Declustering Catalogs");

		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList declusteredCatalog = GardnerKnopoffDeclustering.getDeclusteredCatalog(rupList);
			declusteredList.add(declusteredCatalog);
		}
		if(D) System.out.println("Done Declustering");
		return declusteredList;
	}
	
	
	/**
	 * This declusters according to whether it was spontaneous in U3ETAS
	 * @param catalogList
	 * @return
	 */
	public static ArrayList<ObsEqkRupList> getSpontaneousEventsCatalog(ArrayList<ObsEqkRupList> catalogList) {
		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList filteredCatalog = new ObsEqkRupList();
			for(ObsEqkRupture rup:rupList) {
				if(((ETAS_EqkRupture)rup).getGeneration() == 0)
							filteredCatalog.add(rup);
			}
			declusteredList.add(filteredCatalog);
		}
		
		return declusteredList;
	}

	
	
	/**
	 * This declusters using the UCERF GK filter
	 * @param catalogList
	 * @return
	 */
	public ArrayList<ObsEqkRupList> getU3_GK_FilteredCatalog(ArrayList<ObsEqkRupList> catalogList) {
		int numGriddedRups=0;
		int numFSS_Rups=0;
		GardnerKnopoffAftershockFilter gkFilter = new GardnerKnopoffAftershockFilter(0.05, 9.95, 100);
		ArrayList<ObsEqkRupList> declusteredList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList rupList: catalogList) {
			ObsEqkRupList filteredCatalog = new ObsEqkRupList();
			for(ObsEqkRupture rup:rupList) {
				if(((ETAS_EqkRupture)rup).getFSSIndex() == -1) {
					numGriddedRups+=1;
					double probKeep = gkFilter.getInterpolatedY(rup.getMag());
						if(randomGen.nextDouble()<probKeep)
							filteredCatalog.add(rup);
				}
				else {
					numFSS_Rups+=1;
					if(randomGen.nextDouble()<0.97) // 3% chance of removal
						filteredCatalog.add(rup);
				}
			}
			declusteredList.add(filteredCatalog);
		}
		if(D)
			System.out.println("numGriddedRups="+numGriddedRups+"\nnumFSS_Rups"+numFSS_Rups);
		
		return declusteredList;
	}

	
	/**
	 * This randomizes the original/full catalogs (creating catalogRandomizedList).  This not only randomizes 
	 * event times within each catalog, but also among catalogs so there isn't a high variance of total rates 
	 * among catalogs.  This clones the ruptures and modifies the origin times.
	 * @param random - supply if reproducibility is desired (null otherwise)
	 * @return
	 */
	public ArrayList<ObsEqkRupList> getRandomizedCatalogs() {
		if(catalogRandomizedList==null) {
			if(D) System.out.println("Making randomized catalogs");
			catalogRandomizedList = new ArrayList<ObsEqkRupList>();
			for(int i=0;i<catalogList.size();i++)
				catalogRandomizedList.add(new ObsEqkRupList());
			for(ObsEqkRupList catalog:catalogList) {
				for(ObsEqkRupture rup:catalog) {
					int ithCatalog = (int)Math.floor(randomGen.nextDouble()*catalogList.size());
					long randTimeMillis = startTimeMillis + (long)(randomGen.nextDouble()*fullDurationMillis);
					ObsEqkRupture newRup = (ObsEqkRupture) rup.clone();
					newRup.setOriginTime(randTimeMillis);
					catalogRandomizedList.get(ithCatalog).add(newRup);
				}
			}
			for(int i=0;i<catalogRandomizedList.size();i++) {
				catalogRandomizedList.get(i).sortByOriginTime();
//				System.out.println(catalogRandomizedList.get(i).size()+" events in randomized catalog "+i);
			}
			if(D) System.out.println("Done making randomized catalogs");

		}
		return catalogRandomizedList;
	}
	
	
	/**
	 * This randomizes a single catalog in case you want it to have the exact same set of events
	 * @param catalog
	 * @return
	 */
	public ObsEqkRupList randomizeCatalog(ObsEqkRupList catalog) {
		ObsEqkRupList randCat =new ObsEqkRupList();
		for(ObsEqkRupture rup:catalog) {
			long randTimeMillis = startTimeMillis + (long)(randomGen.nextDouble()*fullDurationMillis);
			ObsEqkRupture newRup = (ObsEqkRupture) rup.clone();
			newRup.setOriginTime(randTimeMillis);
			randCat.add(newRup);
		}
		randCat.sortByOriginTime();
		return randCat;
	}


	

	/**
	 * First three
	 * @return - first three functions (index 0, 1, 2) are orig, GK declustered, and spontaneous only MDFs,
	 * then come the cumulative MFDs of these in the same order (index 3, 4, 5), index 6 is the ratio of 
	 * GK declustered to original incremental MFD, index 7 is spontaneous over orig incr MFD, and index 
	 * 8 is for the UCERF3 GK filter.
	 */
	public ArrayList<XY_DataSet> old_makeCatalogMFDs(boolean mkPopUpPlots) {
		
		IncrementalMagFreqDist mfd = makeCatalogMFD(catalogList);
		double bVal = Math.log10(mfd.getY(5.05)/mfd.getY(6.95))/(6.95-5.05);
		mfd.setName("Full TD Catalog");
		mfd.setInfo("bVal = "+(float)bVal);
		
		IncrementalMagFreqDist mfdDeclustered = makeCatalogMFD(catalogDeclusteredList);
		bVal = Math.log10(mfdDeclustered.getY(5.05)/mfdDeclustered.getY(6.95))/(6.95-5.05);
		mfdDeclustered.setInfo("GK Declustered");
		mfdDeclustered.setInfo("bVal = "+(float)bVal);
		
//		IncrementalMagFreqDist mfdDeclustered_U3_GKfilter = makeMFD(getU3_GK_FilteredCatalog(catalogList));
//		bVal = Math.log10(mfdDeclustered_U3_GKfilter.getY(5.05)/mfdDeclustered_U3_GKfilter.getY(6.95))/(6.95-5.05);
//		mfdDeclustered_U3_GKfilter.setName("Declustered with U3 GK Filter");
//		mfdDeclustered_U3_GKfilter.setInfo("bVal = "+(float)bVal);

		IncrementalMagFreqDist mfdSpontaneousEventsOnly = makeCatalogMFD(getSpontaneousEventsCatalog(catalogList));
		bVal = Math.log10(mfdSpontaneousEventsOnly.getY(5.05)/mfdSpontaneousEventsOnly.getY(6.95))/(6.95-5.05);
		mfdSpontaneousEventsOnly.setName("Spontaneous Events Only");
		mfdSpontaneousEventsOnly.setInfo("bVal = "+(float)bVal);


				// cumulative MFDs
		EvenlyDiscretizedFunc cumMFD = mfd.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD.getClosestYtoX(5.0)/cumMFD.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD.setInfo("bVal = "+(float)bVal);
		
		EvenlyDiscretizedFunc cumMFD_Declustered = mfdDeclustered.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD_Declustered.getClosestYtoX(5.0)/cumMFD_Declustered.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD_Declustered.setInfo("bVal = "+(float)bVal);

//		EvenlyDiscretizedFunc cumMFD_Declustered_U3_GKfilter = mfdDeclustered_U3_GKfilter.getCumRateDistWithOffset();
//		bVal = Math.log10(cumMFD_Declustered_U3_GKfilter.getClosestYtoX(5.0)/cumMFD_Declustered_U3_GKfilter.getClosestYtoX(6.9))/(6.9-5.0);
//		cumMFD_Declustered_U3_GKfilter.setInfo("bVal = "+(float)bVal);


		EvenlyDiscretizedFunc cumMFD_SpontaneousEventsOnly = mfdSpontaneousEventsOnly.getCumRateDistWithOffset();
		bVal = Math.log10(cumMFD_SpontaneousEventsOnly.getClosestYtoX(5.0)/cumMFD_SpontaneousEventsOnly.getClosestYtoX(6.9))/(6.9-5.0);
		cumMFD_SpontaneousEventsOnly.setInfo("bVal = "+(float)bVal);
		
		// ratio plot
		IncrementalMagFreqDist ratio = makeMFD();
		IncrementalMagFreqDist ratio2 = makeMFD();
		for(int i=0;i<ratio.size();i++) {
			ratio.add(i, mfdDeclustered.getY(i)/mfd.getY(i));
			ratio2.add(i, mfdSpontaneousEventsOnly.getY(i)/mfd.getY(i));
		}
		ratio.setName("MFD ratio of GK declustered over Orig");
		ratio2.setName("MFD ratio of Spontaneous only over Orig");

		GardnerKnopoffAftershockFilter gkFilter = new GardnerKnopoffAftershockFilter(0.05, 9.95, 100);

		
		if(mkPopUpPlots) {
			ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
			funcs.add(mfd);
			funcs.add(cumMFD);
			funcs.add(mfdDeclustered);
			funcs.add(cumMFD_Declustered);
//			funcs.add(mfdDeclustered_U3_GKfilter);
//			funcs.add(cumMFD_Declustered_U3_GKfilter);
			funcs.add(mfdSpontaneousEventsOnly);
			funcs.add(cumMFD_SpontaneousEventsOnly);
			
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GREEN));
			Range xRange = new Range(5, 9);
			Range yRange = new Range(1e-5, 20.);
			PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Num", xRange, yRange, 
					false, true, 3.5, 3.0, null, true);
			
			

			ArrayList<XY_DataSet> funcs2 = new  ArrayList<XY_DataSet>();
			funcs2.add(ratio);
			funcs2.add(gkFilter);
			ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			Range yRange2 = new Range(0, 1.1);
			PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars2, null, "Magnitude", "Num", xRange, yRange2, 
					false, false, 3.5, 3.0, null, true);
		}
		
		ArrayList<XY_DataSet> returnList = new ArrayList<XY_DataSet>();
		returnList.add(mfd);
		returnList.add(mfdDeclustered);
		returnList.add(mfdSpontaneousEventsOnly);
		returnList.add(cumMFD);
		returnList.add(cumMFD_Declustered);
		returnList.add(cumMFD_SpontaneousEventsOnly);
		returnList.add(ratio);
		returnList.add(ratio2);
		returnList.add(gkFilter);
		
		return returnList;
	}
	
	
	public static FaultSystemSolutionERF_ETAS getERF() {
		
		FaultSystemSolution sol=null;
		try {
			sol = FaultSystemSolution.load(new File(fssFileName));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		FaultSystemSolutionERF_ETAS erf = ETAS_Launcher.buildERF(sol,false, 1d, 2012);
		erf.updateForecast();
		return erf;
	}
	
	public ArrayList<ObsEqkRupList> load1yrConditionalCatalogs() throws IOException, DocumentException {
		File[] catalogsFileArray = {new File("/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/2024_10_18_1yrConditionalComCatdata_kCOV1p5_MaxPtSrcM/results_m5_preserve_chain.bin")};
		double minMag=5;
		boolean checkAgainstERF =true;
		ArrayList<ObsEqkRupList> oneYrCats = loadCatalogs(catalogsFileArray, minMag, checkAgainstERF);
		
		// check mags & indices against catalogs
		double maxCatLengthYears=0;
		for (ObsEqkRupList catalog : oneYrCats) {
			if(catalog.size()>1) {
				double catLength = ((double)(catalog.get(catalog.size()-1).getOriginTime()-catalog.get(0).getOriginTime()))/(double)millisPerYr;
				if(maxCatLengthYears<catLength)
					maxCatLengthYears=catLength;
			}
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				int nth = rup.getNthERF_Index();
				double magDiff=Math.abs(rup.getMag()/magForNthRup.get(nth) - 1.0);
				if(magDiff > 1e-4)
					throw new RuntimeException("ERROR: mags differ for nth rup: "+nth+"; cat mag = "+rup.getMag()+"; loss-file mag = "+ magForNthRup.get(nth));		
				if(rup.getFSSIndex() != fssIndexForNthRup.get(nth))
					throw new RuntimeException("ERROR: fss indices differ for nth rup: "+nth+"; cat fssID = "+rup.getFSSIndex()+"; loss-file fssID = "+ fssIndexForNthRup.get(nth));		
				if(rup.getGridNodeIndex() != gridCellIndexForNthRup.get(nth))
					throw new RuntimeException("ERROR: grid cell indices differ for nth rup: "+nth+"; cat cellID = "+rup.getGridNodeIndex()+"; loss-file cellID = "+gridCellIndexForNthRup.get(nth));
			}
		}	
		if(D) System.out.println("maxCatLengthYears="+(float)maxCatLengthYears);

		return oneYrCats;
	}
	
	
	/**
	 * The check against the ERF does not work because indexing between ERF and catalogs somehow got screwed up,
	 * but this shouldn't matter here.  The FSS solution ruptures do not have the finite surfaces populated unless
	 * checkAgainstERF is set as true.
	 * @param fssFile
	 * @param catalogsFile
	 * @param minMag
	 * @param checkAgainstERF - checks things against ERF and adds the finite rupture surfaces
	 * @return
	 * @throws IOException
	 * @throws DocumentException
	 */
	public ArrayList<ObsEqkRupList> loadCatalogs(File[] catalogsFileArray, double minMag, boolean checkAgainstERF) throws IOException, DocumentException {
		
		ArrayList<ETAS_Catalog> catalogs = new ArrayList<ETAS_Catalog>();
		for(File catFile:catalogsFileArray)
			catalogs.addAll(ETAS_CatalogIO.loadCatalogsBinary(catFile, minMag)); 
		
		System.out.println("catalogs.size()="+catalogs.size());
		
		// This checks things against ERF and adds the finite rupture surfaces
		if(checkAgainstERF & D) {
			
			if(erf==null)
				erf = getERF();
			
			int numFSS_Rups = erf.getTotNumRupsFromFaultSystem();
			
			double minERF_Mag = 10000;
			for(int r=0;r<erf.getTotNumRups();r++)
				if(minERF_Mag>erf.getNthRupture(numFSS_Rups).getMag())
					minERF_Mag=erf.getNthRupture(numFSS_Rups).getMag();
			
			System.out.println("ERF totNumRups: "+erf.getTotNumRups());
			System.out.println("ERF getTotNumRupsFromFaultSystem: "+numFSS_Rups);
			System.out.println("ERF getNumFaultSystemSources: "+erf.getNumFaultSystemSources());
			System.out.println("ERF getNumSources: "+erf.getNumSources());
			System.out.println("ERF minMag: "+minERF_Mag);
			System.out.println("ERF numGridCells: "+erf.getGridSourceProvider().getGriddedRegion().getNumLocations());
			
			double[] numProbForGridIndex = new double[erf.getGridSourceProvider().getGriddedRegion().getNumLocations()];

			// set the finite rupture surfaces
			int numProbMags =0;
			int numProbNodeIndices =0;
			double minCatalogMag = 10000;
			int maxCatalogNthIndex = 0;
			int totNumCatRups=0;
			int totFSS_catRups=0;
			int totGridSeisCatRups=0;
			HashMap<Integer,Integer> numProbForGeneration = new HashMap<Integer,Integer>();
			HistogramFunction histFunc = new HistogramFunction(0.5, 10, 1.0);
			for (ETAS_Catalog catalog : catalogs) {
				for (ETAS_EqkRupture rup : catalog) {
					int nth = rup.getNthERF_Index();
					int rupFSS_Index = rup.getFSSIndex();
					
					totNumCatRups+=1;
					if(rupFSS_Index != -1)
						totFSS_catRups += 1;
					else
						totGridSeisCatRups += 1;

					
					if(minCatalogMag>rup.getMag())
						minCatalogMag=rup.getMag();
					if(maxCatalogNthIndex<nth)
						maxCatalogNthIndex=nth;

					EqkRupture nthRupture = erf.getNthRupture(nth);

					//Verify FSS indices are the same
					if(nth < numFSS_Rups) { // gridded rups come after FSS rups; next method should really return -1 for gridded rups
						int fssIndexFromERF = erf.getFltSysRupIndexForNthRup(nth);
						if(fssIndexFromERF != rupFSS_Index) {
//							System.out.println("INFO RupMag="+rup.getMag()+"; ERFmag="+nthRupture.getMag());
							throw new RuntimeException("WARNING: fssIndexFromERF != rupFSS_Index for nth="+nth+"; fssIndexFromERF="+
									fssIndexFromERF+" and rupFSS_Index="+rupFSS_Index);
						}
					}
					else {  // check that gridded rups are in the same grid cell; this fails for known reasons
						int nodeIndex1 = erf.getGridSourceProvider().getGriddedRegion().indexForLocation(rup.getHypocenterLocation());
						int nodeIndex2 = rup.getGridNodeIndex();
						if(nodeIndex1 != nodeIndex2) {
							int gen = rup.getGeneration();
							if(numProbForGeneration.keySet().contains(gen)) {
								int num = numProbForGeneration.get(gen);
								numProbForGeneration.replace(gen, num+1);
							}
							else {
								numProbForGeneration.put(gen, 1);
							}
							numProbNodeIndices+=1;
							if(nodeIndex1>=0 && nodeIndex1<numProbForGridIndex.length) {
//								System.out.println("nodeIndex1="+nodeIndex1+"\tnodeIndex2="+nodeIndex2+"\nHypocenter:\n"+rup.getHypocenterLocation());
								numProbForGridIndex[nodeIndex1]+=1;
								double dist = LocationUtils.horzDistance(rup.getHypocenterLocation(), erf.getGridSourceProvider().getGriddedRegion().locationForIndex(nodeIndex2));
								histFunc.add(dist, 1.0);
							}
							if(numProbNodeIndices==1)
								System.out.println("WARNING: nodeIndex1 != nodeIndex2 for nth="+nth+"; nodeIndex1="+
										nodeIndex1+" and nodeIndex2="+nodeIndex2);
//							throw new RuntimeException("nodeIndex1 != nodeIndex2 for nth="+nth+"; nodeIndex1="+
//							nodeIndex1+" and nodeIndex2="+nodeIndex2);
						}
					}
					
					// verify that magnitudes are the same
					double magDiff = Math.abs(rup.getMag()-nthRupture.getMag());
					if(magDiff>0.001) {
						numProbMags+=1;
						if(numProbMags==1)
//							System.out.println("WARNING: rup.getMag() != nthRupture.getMag() for nth="+nth+"; rup.getMag()="+
//									rup.getMag()+" and nthRupture.getMag()="+nthRupture.getMag());
						throw new RuntimeException("rup.getMag() != nthRupture.getMag() for nth="+nth+"; rup.getMag()="+
								rup.getMag()+" and nthRupture.getMag()="+nthRupture.getMag());
					}

					// replace surface if it's a FSS rupture
					if(rupFSS_Index != -1) {
						rup.setRuptureSurface(nthRupture.getRuptureSurface());
						rup.setAveRake(nthRupture.getAveRake());
					}
					else { // TODO this should be what's in the ERF if we want actual surface info
						rup.setAveRake(0.);
						rup.setPointSurface(rup.getHypocenterLocation(), 0., 90.);
					}
				}
			}
			System.out.println("WARNING: "+numProbNodeIndices+" ruptures have hypocenters in the wrong grid node (1st occrrence exmplified above).");
			System.out.println("\nWARNING: "+(float)((double)numProbNodeIndices/(double)totGridSeisCatRups)+
					"% of gridded seismicity leaked into neighboring grid bin (a known issue with U3ETAS). \nThe distance distribution is (betweeyn hypoc and grid center in km):\n");
			int n=0;
			for(int i=0;i<numProbForGridIndex.length;i++)
				if(numProbForGridIndex[i]>0) {
//					System.out.println("\t"+i+"\t"+numProbForGridIndex[i]); 
					n+=1;
				}
			System.out.println(histFunc);
			System.out.println("\t% of grid nodes with leakage: "+(double)n/(double)numProbForGridIndex.length);
			
			System.out.println("\tThis shows the problem is for triggered events and not spontaneous ones");
			for(int gen:numProbForGeneration.keySet())
				System.out.println("\tgen "+gen+":\t"+numProbForGeneration.get(gen));

			System.out.println("\ntotNumCatRups="+totNumCatRups);
			System.out.println("totFSS_catRups="+totFSS_catRups);
			System.out.println("totGridSeisCatRups="+totGridSeisCatRups);
			System.out.println("num problematic mags="+numProbMags);
			System.out.println("minCatalogMag="+minCatalogMag);
			System.out.println("minERF_Mag="+minERF_Mag);
			System.out.println("maxCatalogNthIndex="+maxCatalogNthIndex);

		}
				
		// convert to a list of ObsEqkRupList objects (and filter M≤5 events out)
		ArrayList<ObsEqkRupList> obsEqkRupListList = new ArrayList<ObsEqkRupList>();
		
		int id = -1; // reset the IDs so I can store random loss samples by this.
		for (ETAS_Catalog catalog : catalogs) {
			ObsEqkRupList obsEqkRupList = new ObsEqkRupList();
			for (ETAS_EqkRupture rup : catalog) {
//if(id<200 && id>100)
//	System.out.println(rup.getEventId()+"\t"+rup.getID()+"\t"+rup.getGeneration()+"\t"+rup.getParentID());
//if(id==199)
//	System.exit(id);
				
				if(rup.getMag()>=minMag) {
					// put id in EventID (string) so save the info (this field is presently null)
					rup.setEventId(Integer.toString(rup.getID()));
					id+=1;
					rup.setID(id);
					obsEqkRupList.add(rup);
				}
			}
			obsEqkRupListList.add(obsEqkRupList);
		}
		totalNumEvents = id+1;
		
		return obsEqkRupListList;
	}
	
	
	private static HistogramFunction getBlankLossExceedCurveLogX() {
		return new HistogramFunction(lossCurveLogMin, lossCurveLogMax, lossCurveNum);
	}
	
	private static ArbitrarilyDiscretizedFunc getBlankLossExceedCurve(double initYvalues) {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		HistogramFunction histFunc = getBlankLossExceedCurveLogX();
		for (int i=0; i<histFunc.size(); i++)
			func.set(Math.pow(10d,histFunc.getX(i)), initYvalues);
		return func;
	}
	
	/**
	 * This as the x-axis bin center offset by half bin width
	 * @return
	 */
	private static HistogramFunction getBlankIncrLossCurveLogX() {
		return new HistogramFunction(lossCurveLogMin+lossCurveDelta/2.0, lossCurveLogMax+lossCurveDelta/2.0, lossCurveNum);
	}

	private static ArbitrarilyDiscretizedFunc getBlankIncrLossCurve(double initYvalues) {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		HistogramFunction histFunc = getBlankIncrLossCurveLogX();
		for (int i=0; i<histFunc.size(); i++)
			func.set(Math.pow(10d,histFunc.getX(i)), initYvalues);
		return func;
	}

	
	
	private static ArbitrarilyDiscretizedFunc junk_convertExceedCurveToPDF(UncertainArbDiscFunc curve) {
		HistogramFunction linearXCurve = getBlankLossExceedCurveLogX(); // for x-axis values
		ArbitrarilyDiscretizedFunc resultCurve = new ArbitrarilyDiscretizedFunc();
		for(int i=0;i<curve.size();i++) {
			if(i==0) {
				double linearXval = linearXCurve.getX(0)-linearXCurve.getDelta()/2d;
				double yVal= 1.0 - curve.getY(i);
				resultCurve.set(Math.pow(10d,linearXval),yVal);
			}
			else {
				double linearXval = linearXCurve.getX(i-1)+linearXCurve.getDelta()/2d;
				double yVal= curve.getY(i-1) - curve.getY(i);
				resultCurve.set(Math.pow(10d,linearXval),yVal);
			}
		}
		resultCurve.scale(1.0/resultCurve.calcSumOfY_Vals());
		return resultCurve;
	}
	
	
	/**
	 * This assumes the given distribution has half x-bin offset relative to what is
	 * returned by getBlankLossCurve()
	 * @return
	 */
	private static ArbitrarilyDiscretizedFunc convertIncrToExceedProbDist(UncertainArbDiscFunc incrDist) {
		ArbitrarilyDiscretizedFunc exceedFunc = getBlankLossExceedCurve(0d);
		for(int i=0;i<incrDist.size();i++) {
			double nonExceed=1;
			for(int j=i; j<incrDist.size();j++) {
				nonExceed *= (1.0-incrDist.getY(j));
			}
			exceedFunc.set(i,1.0-nonExceed);
		}
		return exceedFunc;
	}
	
	/**
	 * This assumes the given distribution has half x-bin offset relative to what is
	 * returned by getBlankLossCurve()
	 * @return
	 */
	private static ArbitrarilyDiscretizedFunc convertAggIncrToExceedProbDist(UncertainArbDiscFunc incrDist) {
		ArbitrarilyDiscretizedFunc exceedFunc = getBlankLossExceedCurve(0d);
		for(int i=0;i<incrDist.size();i++) {
			double sum=0;
			for(int j=i; j<incrDist.size();j++) {
				sum += incrDist.getY(j);
			}
			exceedFunc.set(i,sum);
		}
		return exceedFunc;
	}


	
	
	/**
	 * This assumes the given distribution has half x-bin offset relative to what is
	 * returned by getBlankLossCurve()
	 * @return
	 */
	private static ArbitrarilyDiscretizedFunc convertExceedToIncrProbDist(UncertainArbDiscFunc exceedDist) {
		ArbitrarilyDiscretizedFunc incrFunc = getBlankIncrLossCurve(0d);
		for(int i=0;i<exceedDist.size()-1;i++) {
			double incrProb = 1.0 - (1.0-exceedDist.getY(i))/(1.0-exceedDist.getY(i+1));
			incrFunc.set(i,incrProb);
		}
		return incrFunc;
	}

	
	private static ArbitrarilyDiscretizedFunc convertAggrExceedToIncrProbDist(UncertainArbDiscFunc exceedDist) {
		ArbitrarilyDiscretizedFunc incrFunc = getBlankIncrLossCurve(0d);
		for(int i=0;i<exceedDist.size()-1;i++) {
			double incrProb = exceedDist.getY(i)-exceedDist.getY(i+1);
			incrFunc.set(i,incrProb);
		}
		return incrFunc;
	}

	
	
	/**
	 * 	THIS IS WRONG - CAN'T DIFF THE EXCEED PROB AT BIN EDGES; MUST DISAGG NON-EXCEED PROBS

	 * 
	 * First and last bins contain all the probability below and above those x-axis values, respectively.
	 * 
	 * Not sure the best way to scale the y-axis given varying x-axis bin width.
	 * @param curve
	 * @return
	 */
	private static HistogramFunction convertExceedCurveToLog10_PDF_Hist_WRONG(XY_DataSet curve) {
		HistogramFunction hist = new HistogramFunction(lossCurveLogMin-lossCurveDelta/2d, lossCurveLogMax+lossCurveDelta/2d, lossCurveNum+1);

		// set first value
		hist.set(0,1.0 - curve.getY(0)); // put all probability below in first bin

		// set mid values
		for(int i=1;i<curve.size();i++) {
			hist.set(i,curve.getY(i-1) - curve.getY(i));
		}
		// set last value
		hist.set(curve.size(),curve.getY(curve.size()-1)); // put remaining probability in last bin
		
		hist.setName("Histogram from "+curve.getName());
		
		double mean =0;
		for(int i=0;i<hist.size();i++) {
//			double lowX = Math.pow(10,hist.getX(i)-hist.getDelta()/2.0);
//			double highX = Math.pow(10,hist.getX(i)+hist.getDelta()/2.0);
			mean += Math.pow(10,hist.getX(i))*hist.getY(i);
		}
		hist.setInfo("calc mean = "+mean);

		// test result
		HistogramFunction testHist = hist.getCumulativeDistFunctionWithHalfBinOffset();
		for(int i=0;i<curve.size();i++) {
			double x1=curve.getX(i);
			double y1=curve.getY(i);
			double x2=Math.pow(10, testHist.getX(i+1));
			double y2=(1.0-testHist.getY(i+1));
//			System.out.println(x1+"\t"+x2+"\t"+y1+"\t"+y2);
			double diff1 = Math.abs(1.0-(x1/x2));
			double diff2;
			if(y1>0 && y2>0)
				diff2 = Math.abs(1.0-y1/y2);
			else 
				diff2 = Math.abs(y1-y2);
			if(diff1>0.00001 || diff2>0.00001)
				throw new RuntimeException("Problem with convertExceedCurveToLog10_PDF_Hist()\t"+x1+"\t"+x2+"\t"+y1+"\t"+y2+"\t"+diff1+"\t"+diff2);
		}
//		System.exit(0);
		
//		// TEMP
//		for(int i=0;i<hist.size();i++) {
//			double lowX = Math.pow(10,hist.getX(i)-hist.getDelta()/2.0);
//			double highX = Math.pow(10,hist.getX(i)+hist.getDelta()/2.0);
//			hist.set(i,hist.getY(i)/(highX-lowX));
//		}

		
		return hist;
	}
	
	
	


	/**
	 * 	THIS IS WRONG - CAN'T DIFF THE EXCEED PROB AT BIN EDGES; MUST DISAGG NON-EXCEED PROBS
	 * 
	 * First and last bins contain all the probability below and above those x-axis values, respectively
	 * 
	 * Not sure the best way to scale the y-axis given varying x-axis bin width.

	 * @param curve
	 * @return
	 */
	private static HistogramFunction convertExceedCurveToHist_WRONG(UncertainArbDiscFunc curve, int numXaxisBins) {
		
		
		// find the index of the last no-zero Y value (this will remove a ton of zeros)
		double maxXnonZeroY = 0;
		double firstXzeroY = 0;
		for(int i=0;i<curve.size();i++) {
			if(curve.getY(i)>0) {
				maxXnonZeroY = curve.getX(i);
				firstXzeroY = curve.getX(i+1);
			}
			else
				break;
		}

		
		double delta = Math.pow(10, lossCurveLogMax)/numXaxisBins;
		
		int numBins = (int)Math.ceil(maxXnonZeroY/delta)+10;
		
		HistogramFunction hist = new HistogramFunction(delta/2, numBins, delta);
		
		// set first value
		hist.set(0,1.0-curve.getInterpolatedY_inLogXLogYDomain(delta));
		
		for(int i=1;i<hist.size();i++) {
			double xCenter = hist.getX(i);
			if(xCenter+delta/2.0 < maxXnonZeroY) {// both are nonzero
				double yVal = curve.getInterpolatedY_inLogXLogYDomain(xCenter-delta/2.0)-curve.getInterpolatedY_inLogXLogYDomain(xCenter+delta/2.0);
				if(yVal==0)
					yVal = 1e-16;
				hist.set(i,yVal);
			}
			else {  
				if(xCenter-delta/2.0>firstXzeroY) { // both are zero
					hist.set(i,1e-16);	// zero value			
				}
				else { // use non log interpolation because at least one curve interpolation end point is zero
					double yVal = curve.getInterpolatedY(xCenter-delta/2.0)-curve.getInterpolatedY(xCenter+delta/2.0);
					if(yVal==0)
						yVal = 1e-16;
					hist.set(i,yVal);
				}
			}
		}
		
		hist.setName("Histogram from "+curve.getName());
		hist.setInfo("computed mean = "+hist.computeMean());

		return hist;
	}

	
	
	private static IncrementalMagFreqDist makeMFD() {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(minMag, numMag, deltaMag);
		mfd.setTolerance(deltaMag);
		return mfd;
	}
	
	
	private static IncrementalMagFreqDist makeCatalogMFD(ArrayList<ObsEqkRupList> catalogList) {
		IncrementalMagFreqDist mfd = makeMFD();
		for(ObsEqkRupList rupList: catalogList)
			addToMFD(rupList, mfd);
		mfd.scale(1.0/(catalogList.size()*catalogDuration));
		return mfd;
	}

	
	public ArrayList<XY_DataSet> plotLossRateVsMag() {
		IncrementalMagFreqDist lossRate_cat = makeMFD();
		for(ObsEqkRupList rupList: catalogList) {
			for (ObsEqkRupture rup : rupList) {
				int nth = ((ETAS_EqkRupture)rup).getNthERF_Index();
				if(meanLossForNthRup.containsKey(nth))
					lossRate_cat.add(rup.getMag(), meanLossForNthRup.get(nth));
			}
		}
		lossRate_cat.scale(1.0/(catalogList.size()*catalogDuration));
		lossRate_cat.setName("lossRate_cat");
		
		IncrementalMagFreqDist lossRate_catDecl = makeMFD();
		for(ObsEqkRupList rupList: catalogDeclusteredList) {
			for (ObsEqkRupture rup : rupList) {
				int nth = ((ETAS_EqkRupture)rup).getNthERF_Index();
				if(meanLossForNthRup.containsKey(nth))
					lossRate_catDecl.add(rup.getMag(), meanLossForNthRup.get(nth));
			}
		}
		lossRate_catDecl.scale(1.0/(catalogList.size()*catalogDuration));
		lossRate_catDecl.setName("lossRate_catDecl");
		
		IncrementalMagFreqDist lossRate_spontaneous = makeMFD();
		for(ObsEqkRupList rupList: getSpontaneousEventsCatalog(catalogList)) {
			for (ObsEqkRupture rup : rupList) {
				int nth = ((ETAS_EqkRupture)rup).getNthERF_Index();
				if(meanLossForNthRup.containsKey(nth))
					lossRate_spontaneous.add(rup.getMag(), meanLossForNthRup.get(nth));
			}
		}
		lossRate_spontaneous.scale(1.0/(catalogList.size()*catalogDuration));
		lossRate_spontaneous.setName("lossRate_spontaneous");

		
//		if(erf==null)erf = getERF();
//		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
//		erf.updateForecast();
//		IncrementalMagFreqDist lossRate_erf = makeMFD();
//		lossRate_erf.setName("lossRate_erf");
//		double duration = erf.getTimeSpan().getDuration();
//		for(int nth=0; nth<erf.getTotNumRups();nth++) {
//			if(meanLossForNthRup.containsKey(nth))
//				lossRate_erf.add(erf.getNthRupture(nth).getMag(), 
//						meanLossForNthRup.get(nth)*erf.getNthRupture(nth).getMeanAnnualRate(duration));
//		}
		
		// Scale by delta so y values so the AAL
		lossRate_cat.scale(lossRate_cat.getDelta());
//		lossRate_erf.scale(lossRate_erf.getDelta());
		lossRate_catDecl.scale(lossRate_catDecl.getDelta());
		lossRate_spontaneous.scale(lossRate_spontaneous.getDelta());
		
		
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		mfdList.add(lossRate_cat);
//		mfdList.add(lossRate_erf);
		mfdList.add(lossRate_catDecl);
		mfdList.add(lossRate_spontaneous);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		// new Range(1e3,2e4)
		String fileNamePrefix = dirName+"/LossRateVsMag";

		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, "", "Magnitude", "Loss Rate ($ per year)", 
				new Range(5.0,8.5), new Range(1e5,2e6), false, true, fileNamePrefix, true);
		return mfdList;
	}


	
	private static IncrementalMagFreqDist old_makeMagNumDist(ObsEqkRupList rupList) {
		IncrementalMagFreqDist mfd = makeMFD();
		for(ObsEqkRupture rup:rupList)
			mfd.add(rup.getMag(), 1.0);
		return mfd;
	}

	
	


	private static void addToMFD(ObsEqkRupList rupList, IncrementalMagFreqDist mfd) {
		for (ObsEqkRupture rup : rupList) {
			if(rup.getMag()>1)
				mfd.add(rup.getMag(), 1.0);
		}
	}

	
	public void old_computeNumEventDistribution(ArrayList<ObsEqkRupList> catalogList, double duration, double[] magArray) {
		
		// get subcatalogs
		ArrayList<ObsEqkRupList> eqkRupListList = getSubcatalogList(catalogList, duration);
		
		EvenlyDiscretizedFunc tempFunc = old_makeMagNumDist(eqkRupListList.get(0)).getCumRateDistWithOffset();
		ArbDiscrEmpiricalDistFunc_3D cumMFD_FromAllCatalogsFunc_3D = new ArbDiscrEmpiricalDistFunc_3D(tempFunc.getMinX(), tempFunc.size(), tempFunc.getDelta());
		
		int showProgressAt = eqkRupListList.size()/10;
		for(int i=0;i<eqkRupListList.size();i++) {
			if(i>showProgressAt && D) {
				if(D) System.out.println(i);
				showProgressAt+=eqkRupListList.size()/10;
			}
			IncrementalMagFreqDist func = old_makeMagNumDist(eqkRupListList.get(i));
			cumMFD_FromAllCatalogsFunc_3D.set(func.getCumRateDistWithOffset(), 1.0);
		}		
		
		EvenlyDiscretizedFunc meanExpNumDist = cumMFD_FromAllCatalogsFunc_3D.getMeanCurve();
		
		System.out.println(meanExpNumDist);

		
		for(double mag:magArray) {
//			System.out.println("working on M="+mag);
			int magIndex = meanExpNumDist.getClosestXIndex(mag);
			double expNum = meanExpNumDist.getY(magIndex);
//			System.out.println("Exp Num at M="+mag+" is: "+(float)meanExpNumDist.getY(magIndex));
			ArbDiscrEmpiricalDistFunc numDist = cumMFD_FromAllCatalogsFunc_3D.getArbDiscrEmpDistFuncArray()[magIndex].getNormalizedDist();
//			System.out.println(numDist);

			int maxNum = (int)Math.round(numDist.getMaxX());
			EvenlyDiscretizedFunc poissNumDist = new EvenlyDiscretizedFunc(0.,maxNum+1,1.0);
			PoissonDistribution poisDist = new PoissonDistribution(expNum);
			for(int i=0;i<poissNumDist.size();i++) {
				poissNumDist.set(i,poisDist.probability(i));
			}
			numDist.setName("Catalog Num Distribution for M="+mag);
			numDist.setInfo("ExpNum = : "+expNum);
			poissNumDist.setName("Poisson Num Distribution for M="+mag);
			poissNumDist.setInfo("ExpNum = : "+expNum);
			ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
			plottingFuncsArray.add(numDist);
			plottingFuncsArray.add(poissNumDist);
			this.quickPlot(plottingFuncsArray,"Num","Probability","M="+mag+"; Duration="+duration, false, false);

		}

	}
	
	
	/**
	 * This divides the catalogs up into subcatalogs of the given duration
	 * @param catalogList
	 * @param duration - yrs
	 * @return
	 */
	public static ArrayList<ObsEqkRupList> old_getSubcatalogList(ArrayList<ObsEqkRupList> catalogList, double duration) {
		
		int numSubCatPerCat = (int)Math.floor((catalogDuration+0.003)/duration);	// add a day (0.03) to make sure we get the last window
		if(D) System.out.println("numSubCatPerCat = "+numSubCatPerCat);
		
		int numRupsNotUsed = 0;
		for(ObsEqkRupList catalog : catalogList) 
			numRupsNotUsed += catalog.size();

		ArrayList<ObsEqkRupList> eqkRupListList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList catalog : catalogList) {
			long endEpoch = (long)((catalogStartYear-1970)*millisPerYr) + (long)(duration*millisPerYr);
			ObsEqkRupList currEqkList = new ObsEqkRupList();
			eqkRupListList.add(currEqkList);
			for(ObsEqkRupture rup: catalog) {
				if(rup.getOriginTime()<endEpoch) {
					currEqkList.add(rup);
					numRupsNotUsed -= 1;
				}
				else {
					currEqkList = new ObsEqkRupList();
					eqkRupListList.add(currEqkList);
					currEqkList.add(rup);
					numRupsNotUsed -= 1;
					endEpoch += (long)(duration*millisPerYr);
				}
			}
		}
		if(D) System.out.println("eqkRupListList.size() = "+eqkRupListList.size()+"; it should be "+numSubCatPerCat*catalogList.size()+
				"; numRupsNotUsed = "+numRupsNotUsed);
		
		return eqkRupListList;
	}
	
	
	/**
	 * This divides the catalogs up into subcatalogs of the given duration
	 * @param catalogList
	 * @param subDurationYrs - yrs
	 * @return
	 */
	public ArrayList<ObsEqkRupList> getSubcatalogList(ArrayList<ObsEqkRupList> catalogList, double subDurationYrs) {
		
		long runtime = System.currentTimeMillis();
		long subDurationMillis = (long)(millisPerYr*subDurationYrs);
		
		int numSubCatPerCat = (int)(fullDurationMillis/subDurationMillis); // this should truncate
		if(D) System.out.println("numSubCatPerCat="+numSubCatPerCat);
		long remainderMillis = fullDurationMillis - numSubCatPerCat*subDurationMillis;
		double remainderFraction = (double)remainderMillis/(double)subDurationMillis;
		if(D) System.out.println("remainderFraction="+remainderFraction);

//		numSubCatPerCat = (int)((endTimeMillis-1000-startTimeMillis)/durationMillis);
//		System.out.println("numSubCatPerCat="+numSubCatPerCat);
//		System.exit(0);
				
		int numRupsNotUsed = 0;
		ArrayList<ObsEqkRupList> eqkRupListList = new ArrayList<ObsEqkRupList>();
		for(ObsEqkRupList catalog : catalogList) {
			ArrayList<ObsEqkRupList> subCatList = new ArrayList<ObsEqkRupList>();
			for(int i=0;i<numSubCatPerCat;i++)
				subCatList.add(new ObsEqkRupList());
			for(ObsEqkRupture rup: catalog) {
				int catIndex = (int)((rup.getOriginTime()-startTimeMillis)/subDurationMillis);
				if(catIndex<subCatList.size())
					subCatList.get(catIndex).add(rup);
				else {
					numRupsNotUsed+=1;
//					if(rup.getOriginTime()<endTimeMillis) {
//						throw new RuntimeException("Something wrong in getSubcatalogList()");
//					}
				}
			}
			eqkRupListList.addAll(subCatList);
		}
		runtime = System.currentTimeMillis()-runtime;

		if(D) System.out.println("Subcatalog eqkRupListList.size() = "+eqkRupListList.size()+"; it should be "+numSubCatPerCat*catalogList.size()+
				"; numRupsNotUsed = "+numRupsNotUsed+"; runtime (sec) = "+(runtime/1000));
		
		return eqkRupListList;
	}

	
	
	

	
	
	
	
	
	/**
	 * This plots a 500-year cumulative number of events versus times for the least, most, and median catalogs in terms
	 * total number of events (at 500 years).
	 * @param catalogList
	 * @param dirName
	 * @param popupWindow
	 * @param plotTitle
	 */
	public void old_writeAndOrPlotRateVersusTime(ArrayList<ObsEqkRupList> catalogList, String dirName, boolean popupWindow, String plotTitle) {
		
		double[] totNumArray = new double[catalogList.size()];
		
		ArbitrarilyDiscretizedFunc[] allFuncsArray = new ArbitrarilyDiscretizedFunc[catalogList.size()];

		for(int i=1;i<catalogList.size();i++ ) {
			ArbitrarilyDiscretizedFunc cumNumMgt5 = new ArbitrarilyDiscretizedFunc();
			long firstTime = catalogList.get(i).get(0).getOriginTime();
			int num = 1;
			for(ObsEqkRupture rup : catalogList.get(i)) {
				double yrs = ((double)(rup.getOriginTime()-firstTime))/(double)millisPerYr;
				cumNumMgt5.set(yrs,num);
				num+=1;
			}
			allFuncsArray[i] = cumNumMgt5;
			totNumArray[i] = cumNumMgt5.getY(cumNumMgt5.size()-1);
		}
		
		double median = Quantiles.median().compute(totNumArray);
		double min = Double.MAX_VALUE;
		double max = Double.NEGATIVE_INFINITY;
		int medianIndex=-1, minIndex=-1, maxIndex=-1;
		double minMedianDiff = Double.MAX_VALUE;
		
		for(int i=1;i<catalogList.size();i++ ) {
			double medianDiff = Math.abs(totNumArray[i]-median);
			if(medianDiff < minMedianDiff) {
				minMedianDiff = medianDiff;
				medianIndex=i;
			}
			if(totNumArray[i]<min) {
				min = totNumArray[i];
				minIndex=i;
			}
			if(totNumArray[i]>max) {
				max = totNumArray[i];
				maxIndex=i;
			}
		}		
		
		if(D) {
			double maxMag = 0;
			long ot = 0;
			for(ObsEqkRupture rup : catalogList.get(maxIndex)) {
				if(rup.getMag()>maxMag) {
					maxMag = rup.getMag();
					ot = rup.getOriginTime();
				}
			}
			double yrs = ((double)ot - (catalogStartYear-1970.)*millisPerYr)/millisPerYr;
			System.out.println("maxMag = "+maxMag+" at "+yrs+" yrs");
		}


		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plottingFuncsArray.add(allFuncsArray[medianIndex]);
		plottingFuncsArray.add(allFuncsArray[minIndex]);
		plottingFuncsArray.add(allFuncsArray[maxIndex]);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/cumNumEventsVsTime"+plotTitle;
		String xAxisLabel = "Time (years)";
		String yAxisLabel = "Cumulative Number of Events";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
	}

	
	public static void quickPlot(ArrayList<XY_DataSet> plottingFuncsArray, String xAxisLabel, String yAxisLabel, String title, boolean logX, boolean logY) {
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		
		Range xAxisRange = null;
		Range yAxisRange = null;
//		boolean logX = false;
//		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, title, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, null, true);
	}
	
	
	
	public static void writeAndOrPlotLossCurves(ArrayList<UncertainArbDiscFunc[]> dataSetsArray, ArrayList<XY_DataSet> funcsArray, double duration, 
			double saPeriod, String dirName, boolean popupWindow, String plotTitle) {
		
		Color[] colorArray = {Color.BLUE, Color.RED, Color.BLACK, Color.MAGENTA, Color.CYAN};
		int colorIndex = 0;
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	

		for(UncertainArbDiscFunc[] dataSets:dataSetsArray) {
			plottingFuncsArray.add(dataSets[1]); // for solid line
			plottingFuncsArray.add(dataSets[1]); // for shaded region
			plottingFuncsArray.add(dataSets[0].getLower());
			plottingFuncsArray.add(dataSets[0].getUpper());
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, colorArray[colorIndex]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, colorArray[colorIndex]));
			colorIndex += 1;
		}
		
		
		// add other curves
		for(XY_DataSet func:funcsArray) {
			plottingFuncsArray.add(func);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colorArray[colorIndex]));
			colorIndex += 1;
		}
		
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/"+plotTitle+"_lossCurve";
		String xAxisLabel = "FixLabel";
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = new Range(1e-2,10);
		Range yAxisRange = new Range(1e-4,1.0);
		boolean logX = true;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}
	
	
	

	

	
	
	

	
	
	
	public ArrayList<ObsEqkRupList> getCatalogs() {
		return catalogList;
	}
	
	public ArrayList<ObsEqkRupList> getCatalogsDeclustered() {
		return catalogDeclusteredList;
	}
	
	


	private ArrayList<ObsEqkRupList> old_makeRandomCatalogListFromERF(FaultSystemSolutionERF erf) {
		ArrayList<ObsEqkRupList> catalogList = new ArrayList<ObsEqkRupList>();
		erf.getTimeSpan().setDuration(catalogDuration);
		erf.updateForecast();
		int numCatalogs = 1000;
		for(int i=0;i<numCatalogs; i++) {
			List<EqkRupture> rupList = erf.drawRandomEventSet();
			ObsEqkRupList obsRupList = new ObsEqkRupList();
			for(EqkRupture rup:rupList) {
				if(rup.getMag()<5)
					continue;
				ObsEqkRupture obsRup = new ObsEqkRupture();
				obsRup.setAveRake(rup.getAveRake());
				obsRup.setMag(rup.getMag());
				obsRup.setRuptureSurface(rup.getRuptureSurface());
				double randYrs = (catalogStartYear-1970.) + randomGen.nextDouble()*catalogDuration;
				long ot = (long)(randYrs*this.millisPerYr);
				obsRup.setOriginTime(ot);
				obsRupList.add(obsRup);
			}
			obsRupList.sortByOriginTime();
			catalogList.add(obsRupList);
		}
				
				return catalogList;
	}
	
	
	/**
	 * This demonstrates that when you click x-axis log plot for a PDF (HistogramFunc) the distribution
	 * is distorted because non-even x-axis binning is not being accounted for.  This is exemplified
	 * by the mode of the lognormal distribution not being equal to exp(mu).  It also means a uniform
	 * distribution will not look uniform if corrected for this effect when log plot chosen.
	 */
	public static void testTakingLogOfPDF_Xaxis() {
		double meanLn = 1.0;
		double stdDevLn = 0.5;
		System.out.println("meanLn="+meanLn+"\nstdDevLn="+stdDevLn);
		double mode = Math.exp(meanLn-stdDevLn*stdDevLn);
		double median = Math.exp(meanLn);
		double mean = Math.exp(meanLn+stdDevLn*stdDevLn/2);
		double var = (Math.exp(stdDevLn*stdDevLn)-1)*Math.exp(2*meanLn+stdDevLn*stdDevLn);
		double cov = Math.sqrt(var)/mean;
		System.out.println("Mode = "+mode+"\nMedian = "+median+"\nMean = "+mean+"\nCOV = "+cov);
		
		int numPoints = 1000;
		double maxLn = meanLn*3;
		double minLn = -meanLn*3;
		double deltaLn = 2*maxLn/numPoints;
		
//		org.opensha.commons.calc.GaussianDistCalc gaussDistCalc = new org.opensha.commons.calc.GaussianDistCalc();
		HistogramFunction normEvenPDF = new HistogramFunction(minLn,numPoints,deltaLn);
		for(int i=1;i<normEvenPDF.size();i++) {
			double normVal = (normEvenPDF.getX(i)-meanLn)/stdDevLn;
			double yVal = Math.exp(-0.5*normVal*normVal);
			normEvenPDF.set(i,yVal);
		}
		
		double sum1 = normEvenPDF.calcSumOfY_Vals()*normEvenPDF.getDelta();
		normEvenPDF.scale(1.0/sum1);
		normEvenPDF.setName("normEvenPDF");
		
		ArbitrarilyDiscretizedFunc logNormArbPDF = new ArbitrarilyDiscretizedFunc();
		double sum = 0;
		for(int i=0;i<normEvenPDF.size();i++) {
			double binWidth = Math.exp(normEvenPDF.getX(i)+normEvenPDF.getDelta()/2d) - Math.exp(normEvenPDF.getX(i)-normEvenPDF.getDelta()/2d);
			logNormArbPDF.set(Math.exp(normEvenPDF.getX(i)),normEvenPDF.getY(i)/binWidth);
//			logNormArbPDF.set(Math.exp(normEvenPDF.getX(i)),normEvenPDF.getY(i));
		}
		logNormArbPDF.setName("logNormArbPDF");
		
		double max = logNormArbPDF.getMaxX();
		double min = logNormArbPDF.getMinX();
		double delta = (max-min)/numPoints;
		
		LognormalDistCalc logNormCalc = new LognormalDistCalc();
//		logNormCalc.setAll(mean, cov, delta, numPoints*5);
		logNormCalc.setAllParameters(mean, cov, delta, numPoints*5, 1d, 0d); // numPoints multiplied by 5 to get the right computed mean
		EvenlyDiscretizedFunc logNormEvenPDF = logNormCalc.getPDF();
		logNormEvenPDF.setName("logNormEvenPDF");
//		System.out.println("totalDensityCheck="+logNormEvenPDF.calcSumOfY_Vals()*logNormEvenPDF.getDelta());
		
		ArbitrarilyDiscretizedFunc normArbPDF = new ArbitrarilyDiscretizedFunc();
		for(int i=1;i<logNormEvenPDF.size();i++) {
			double binWidth = Math.log(logNormEvenPDF.getX(i)+logNormEvenPDF.getDelta()/2d) - Math.log(logNormEvenPDF.getX(i)-logNormEvenPDF.getDelta()/2d);
			normArbPDF.set(Math.log(logNormEvenPDF.getX(i)),logNormEvenPDF.getY(i)/binWidth);
		}
		normArbPDF.setName("normArbPDF");
		
		// Make max values the same
		logNormArbPDF.scale(logNormEvenPDF.getMaxY()/logNormArbPDF.getMaxY());
		normArbPDF.scale(normEvenPDF.getMaxY()/normArbPDF.getMaxY());
		
		ArrayList<XY_DataSet> funcList1 = new ArrayList<XY_DataSet>(); 
		funcList1.add(normEvenPDF);
		funcList1.add(normArbPDF);
		quickPlot(funcList1, "X", "Prob", "Normal Dist", false, false);
		
		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); 
		funcList2.add(logNormEvenPDF);
		funcList2.add(logNormArbPDF);
		quickPlot(funcList2, "X", "Prob", "LogNormal Dist", false, false);

	}
	
	
	
	
	/**
	 * 
	 * @param funcList
	 * @param includeUncert - if true, columns with normalized standard deviation (std/loss) are added
	 * @return
	 */
	public static String getTableStringOfCEA_Values(ArrayList<UncertainArbDiscFunc> funcList, boolean includeUncert) {
		String tableString = "";
		for(int row = 0; row<CEA_PROB_ARRAY.length; row++) {
			double prob = CEA_PROB_ARRAY[row];
			for(int col = 0; col<funcList.size()+1; col++) {
				if(col==0) {
					int yrs = (int)(1d/prob);
					if(yrs==350)
						prob = 0.00286;
					tableString +=  prob+" (1/"+yrs+")\t";
				}
				else {
					double loss = funcList.get(col-1).getFirstInterpolatedX_inLogXLogYDomain(prob);
					tableString +=  Float.toString((float)loss);
					if(includeUncert) {
						double std = (funcList.get(col-1).getUpper().getFirstInterpolatedX_inLogXLogYDomain(prob)-
								funcList.get(col-1).getLower().getFirstInterpolatedX_inLogXLogYDomain(prob))/4.0;
						tableString +=  "\t"+Float.toString((float)(std/loss));

					}
					if(col<funcList.size())
						tableString +=  "\t";
				}
			}
			tableString +=  "\n";
		}
		return tableString;
	}

	
	// first column is prob
	public static String getTableStringOfCEA_Ratios(ArrayList<UncertainArbDiscFunc> funcList, boolean includeUncert) {
		String tableString = "";
		for(int row = 0; row<CEA_PROB_ARRAY.length; row++) {
			double prob = CEA_PROB_ARRAY[row];
			double tdLoss = funcList.get(0).getFirstInterpolatedX_inLogXLogYDomain(prob);
			double tdStd = (funcList.get(0).getUpper().getFirstInterpolatedX_inLogXLogYDomain(prob)-
							funcList.get(0).getLower().getFirstInterpolatedX_inLogXLogYDomain(prob))/4.0;
			for(int col = 0; col<funcList.size()+1; col++) {
				if(col==0) {
					int yrs = (int)(1d/prob);
					if(yrs==350)
						prob = 0.00286;
					tableString +=  prob+" (1/"+yrs+")\t";
				}
				else if(col==1) {
					continue; // skip
				}
				else {
					double loss = funcList.get(col-1).getFirstInterpolatedX_inLogXLogYDomain(prob);
					tableString +=  Float.toString((float)(loss/tdLoss));
					if(includeUncert) {
						double std = (funcList.get(col-1).getUpper().getFirstInterpolatedX_inLogXLogYDomain(prob)-
								funcList.get(col-1).getLower().getFirstInterpolatedX_inLogXLogYDomain(prob))/4.0;
						double uncert = Math.sqrt(Math.pow(std/loss,2) + Math.pow(tdStd/tdLoss,2)); // add in quadrature
						tableString +=  "\t"+Float.toString((float)uncert);
					}

					if(col<funcList.size())
						tableString +=  "\t";
				}
			}
			tableString +=  "\n";
		}
		return tableString;
	}
	
	public void make1yrConditionalLossExceedComparisonFig() {
		
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		
		UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1d);
		lossExceedProbFromTD_CatalogsRandFromCOV.setName("lossExceedProbFromTD_CatalogsRandFromCOV");
		
		
		System.out.println("DOING ONE YEAR CATS");
		ArrayList<ObsEqkRupList> oneYrCondCats=null;
		try {
			oneYrCondCats = load1yrConditionalCatalogs();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
		// Need to make exceedance curve by hand because method above subdivides catalogs and not sure it will work
		double aveYearlyLoss=0;
		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		for (ObsEqkRupList catalog : oneYrCondCats) {
			double lossSum=0;
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(0d);
			double maxLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxLoss<randLoss) maxLoss=randLoss;
				lossSum+=randLoss;
			}
			aveYearlyLoss+=lossSum;
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < maxLoss)
					lossProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}		
		aveYearlyLoss /= (double)oneYrCondCats.size();
		lossProbCurve.scale(1.0/((double)oneYrCondCats.size()));
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), oneYrCondCats.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+oneYrCondCats.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}
		UncertainArbDiscFunc oneYrCondExceedUncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		oneYrCondExceedUncertFunc.setName("oneYrCondExceedUncertFunc");
		
		funcList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
		funcList.add(oneYrCondExceedUncertFunc);
		
		// Write tables
		String tableString = "";
		tableString += "\n1yr Loss for CEA probability levels (cond loss PDF for each rup randomly sampled):\n"+"\n";
		tableString += "Probability\tTD_loss\tTD_uncert\tCond2024_loss\tCond2024_uncert"+"\n";
		tableString += getTableStringOfCEA_Values(funcList,true)+"\n";
		tableString += "\n1yr Loss ratios, relative to TD, for CEA probability levels (cond loss PDF for each rup randomly sampled):\n"+"\n";
		tableString += "Probability\tCond2024_loss\tCond2024_uncert"+"\n";
		tableString += getTableStringOfCEA_Ratios(funcList, true)+"\n";
		System.out.println(tableString);
		oneYrCondExceedUncertFunc.setInfo("2024 Expected Loss = "+aveYearlyLoss+"\n\n"+tableString);

		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 180, 0)));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));

		File outputDir = new File(dirName+ "/OneYrCondLossComFigure");
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/OneYrCondLossComFigure/"+"OneYrCondLossComFigure";
		fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Loss ($)", "Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
		fileNamePrefix += "_zoomed";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Loss ($)", "Exceedance Probability", new Range(1e6,1e11), new Range(1e-4,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		UncertainArbDiscFunc denomCurve = funcList.get(0);
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		double xThresh = denomCurve.getClosestXtoY(1e-4);
		for(int i=0;i<funcList.size();i++) {
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			ratioCurve.setName("Ratio for "+funcList.get(i).getName());
			for(int j=0;j<denomCurve.size();j++) {
				if(denomCurve.getX(j)>xThresh)
					break;
				double ratio = funcList.get(i).getY(j)/denomCurve.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(denomCurve.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
		}
		fileNamePrefix += "_Ratio";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Loss ($)", "Exceedance Ratio", new Range(1e4,1e11), new Range(0.0,2.0), true, false, fileNamePrefix, true);

//		// plot incremental distributions
//		ArrayList<XY_DataSet> pdfList = new ArrayList<XY_DataSet>(); 
//		for(int i=0;i<2;i++) { // Only the first two funcs
//			pdfList.add(convertExceedToIncrProbDist(funcList.get(i)));
//		}
//		if(randCOV) {
//			UncertainArbDiscFunc lossIncreProbDistFromTD_CatalogsRandFromCOV = getLossIncrProbDistFromCatalogsRandomFromCOV(getCatalogs(), 1d);
//			lossIncreProbDistFromTD_CatalogsRandFromCOV.setName("lossIncreProbDistFromTD_CatalogsRandFromCOV");
//			UncertainArbDiscFunc lossIcreProbDistFromPoisCatalogsRandFromCOV = getLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1d);
//			lossIcreProbDistFromPoisCatalogsRandFromCOV.setName("lossIcreProbDistFromPoisCatalogsRandFromCOV");
//			pdfList.add(lossIncreProbDistFromTD_CatalogsRandFromCOV);
//			pdfList.add(lossIcreProbDistFromPoisCatalogsRandFromCOV);
//		}
//		ArrayList<PlotCurveCharacterstics> plotChars3 = new ArrayList<PlotCurveCharacterstics>();	
//		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.RED));
//		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
//		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
//		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
//		fileNamePrefix = dirName+"/FigureA/"+durationsYrs+"yr_IncrLoss";
//		if(randCOV) fileNamePrefix += "_randCOV";
//		PlottingUtils.writeAndOrPlotFuncs(pdfList, plotChars3, "", "Loss ($1)", "Incremental Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
	}

	public void make1yrConditionalAggrLossExceedComparisonFig() {
		
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		
		UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1d);
		lossExceedProbFromTD_CatalogsRandFromCOV.setName("aggrLossExceedProbFromTD_CatalogsRandFromCOV");
		
		
		System.out.println("DOING ONE YEAR CATS");
		ArrayList<ObsEqkRupList> oneYrCondCats=null;
		try {
			oneYrCondCats = load1yrConditionalCatalogs();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
		// Need to make exceedance curve by hand because method above subdivides catalogs and not sure it will work
		double aveYearlyLoss=0;
		ArbitrarilyDiscretizedFunc lossProbCurve = getBlankLossExceedCurve(0d);
		for (ObsEqkRupList catalog : oneYrCondCats) {
			double lossSum=0;
			ArbitrarilyDiscretizedFunc lossProbCurveForSubCat = getBlankLossExceedCurve(0d);
//			double maxLoss = 0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double aveLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(aveLoss==0)
					continue;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
//				if(maxLoss<randLoss) maxLoss=randLoss;
				lossSum+=randLoss;
			}
			aveYearlyLoss+=lossSum;
			for(int i=0;i<lossProbCurveForSubCat.size();i++) {
				if(lossProbCurveForSubCat.getX(i) < lossSum)
					lossProbCurveForSubCat.set(i,1.0);
				else
					break;
			}
			for(int i=0;i<lossProbCurve.size();i++) {
				lossProbCurve.set(i,lossProbCurve.getY(i)+lossProbCurveForSubCat.getY(i));
			}
		}		
		aveYearlyLoss /= (double)oneYrCondCats.size();
		lossProbCurve.scale(1.0/((double)oneYrCondCats.size()));
		ArbitrarilyDiscretizedFunc low95confCurve = getBlankLossExceedCurve(0d);
		ArbitrarilyDiscretizedFunc high95confCurve = getBlankLossExceedCurve(0d);
		for(int i=0;i<lossProbCurve.size();i++) {
			double xVal = lossProbCurve.getX(i);
			double[] confArray = ETAS_Utils.getBinomialProportion95confidenceInterval(lossProbCurve.getY(i), oneYrCondCats.size());
			if(confArray[0]>lossProbCurve.getY(i)) {  // numerical problem
				if(D) System.out.println("ISSUE: Low bound above mean:\nmean="+lossProbCurve.getY(i)+"\nlowBound="+confArray[0]+"\nN="+oneYrCondCats.size());
				confArray[0] = 0;
			}
			low95confCurve.set(xVal,confArray[0]);
			high95confCurve.set(xVal,confArray[1]);
		}
		UncertainArbDiscFunc oneYrCondExceedUncertFunc = new UncertainArbDiscFunc(lossProbCurve,low95confCurve,high95confCurve);
		oneYrCondExceedUncertFunc.setName("oneYrCondAggrExceedUncertFunc");
		
		funcList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
		funcList.add(oneYrCondExceedUncertFunc);
		
		// Write tables
		String tableString = "";
		tableString += "\n1yr Aggregate Loss for CEA probability levels (cond loss PDF for each rup randomly sampled):\n"+"\n";
		tableString += "Probability\tTD_loss\tTD_uncert\tCond2024_loss\tCond2024_uncert"+"\n";
		tableString += getTableStringOfCEA_Values(funcList,true)+"\n";
		tableString += "\n1yr Aggregate Loss ratios, relative to TD, for CEA probability levels (cond loss PDF for each rup randomly sampled):\n"+"\n";
		tableString += "Probability\tCond2024_loss\tCond2024_uncert"+"\n";
		tableString += getTableStringOfCEA_Ratios(funcList, true)+"\n";
		System.out.println(tableString);
		oneYrCondExceedUncertFunc.setInfo("2024 Expected Loss = "+aveYearlyLoss+"\n\n"+tableString);

		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(0, 180, 0)));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));

		File outputDir = new File(dirName+ "/OneYrCondAggrLossComFigure");
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/OneYrCondAggrLossComFigure/"+"OneYrCondAggrLossComFigure";
		fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Aggregate Loss ($)", "Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
		fileNamePrefix += "_zoomed";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Aggregate Loss ($)", "Exceedance Probability", new Range(1e6,1e11), new Range(1e-4,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		UncertainArbDiscFunc denomCurve = funcList.get(0);
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		double xThresh = denomCurve.getClosestXtoY(1e-4);
		for(int i=0;i<funcList.size();i++) {
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			ratioCurve.setName("Ratio for "+funcList.get(i).getName());
			for(int j=0;j<denomCurve.size();j++) {
				if(denomCurve.getX(j)>xThresh)
					break;
				double ratio = funcList.get(i).getY(j)/denomCurve.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(denomCurve.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
		}
		fileNamePrefix += "_Ratio";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Aggregate Loss ($)", "Exceedance Ratio", new Range(1e4,1e11), new Range(0.0,2.0), true, false, fileNamePrefix, true);

	}

	
	/**
	 * This is for one or more exceedances
	 * @param randCOV
	 */
	public void makeFigAndTableSetA_Exceedances(double durationsYrs, boolean randCOV) {
				
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		
		String randCOVincluded = "";
		if(randCOV)
			randCOVincluded = " (cond loss PDF for each rup randomly sampled)";

		if(randCOV) {
			UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), durationsYrs);
			lossExceedProbFromTD_CatalogsRandFromCOV.setName("lossExceedProbFromTD_CatalogsRandFromCOV");
			UncertainArbDiscFunc lossExceedProbFromPoisCatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), durationsYrs);
			lossExceedProbFromPoisCatalogsRandFromCOV.setName("lossExceedProbFromPoisCatalogsRandFromCOV");
			UncertainArbDiscFunc lossExceedProbFromDeclusteredCatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogsDeclustered(), durationsYrs);
			lossExceedProbFromDeclusteredCatalogsRandFromCOV.setName("lossExceedProbFromDeclusteredCatalogsRandFromCOV");
			UncertainArbDiscFunc lossExceedProbFromSpontaneousCatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getSpontaneousEventsCatalog(getCatalogs()), durationsYrs);
			lossExceedProbFromSpontaneousCatalogsRandFromCOV.setName("lossExceedProbFromSpontaneousCatalogsRandFromCOV");
			funcList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
			funcList.add(lossExceedProbFromPoisCatalogsRandFromCOV);
			funcList.add(lossExceedProbFromDeclusteredCatalogsRandFromCOV);
			funcList.add(lossExceedProbFromSpontaneousCatalogsRandFromCOV);			
		}
		else {
			UncertainArbDiscFunc lossExceedProbFromTD_Catalogs = getLossExceedProbCurveFromCatalogs(getCatalogs(), durationsYrs);
			lossExceedProbFromTD_Catalogs.setName("lossExceedProbFromTD_Catalogs");
			UncertainArbDiscFunc lossExceedProbFromPoisCatalogs = getLossExceedProbCurveFromCatalogs(getRandomizedCatalogs(), durationsYrs);
			lossExceedProbFromPoisCatalogs.setName("lossExceedProbFromPoisCatalogs");
			UncertainArbDiscFunc lossExceedProbFromDeclusteredCatalogs = getLossExceedProbCurveFromCatalogs(getCatalogsDeclustered(), durationsYrs);
			lossExceedProbFromDeclusteredCatalogs.setName("lossExceedProbFromDeclusteredCatalogs");
			UncertainArbDiscFunc lossExceedProbFromSpontaneousCatalogs = getLossExceedProbCurveFromCatalogs(getSpontaneousEventsCatalog(getCatalogs()), durationsYrs);
			lossExceedProbFromSpontaneousCatalogs.setName("lossExceedProbFromSpontaneousCatalogs");
			funcList.add(lossExceedProbFromTD_Catalogs);
			funcList.add(lossExceedProbFromPoisCatalogs);
			funcList.add(lossExceedProbFromDeclusteredCatalogs);
			funcList.add(lossExceedProbFromSpontaneousCatalogs);
		}
		
		// Write tables
		if(durationsYrs==1.0) { // only for 1-year forecast durations
			System.out.println("\n"+durationsYrs+"yr Loss for CEA probability levels"+randCOVincluded+":\n");
			System.out.println("Probability\tTD_loss\tTD_uncert\tTI_loss\tTI_uncert\tGK_loss\tGK_uncert\tSP_loss\tSP_uncert");
			System.out.println(getTableStringOfCEA_Values(funcList,true));
			System.out.println("\n"+durationsYrs+"yr Loss ratios, relative to TD, for CEA probability levels"+randCOVincluded+":\n");
			System.out.println("Probability\tTIvsTD_loss\tTIvsTD_uncert\tGKvsTD_loss\tGKvsTD_uncert\tSPvsTD_loss\tSPvsTD_uncert");
			System.out.println(getTableStringOfCEA_Ratios(funcList, true));
			System.out.println("TI = Time Dependent\nTI = Time Independent (Poisson)\nGK = Gardner Knopoff declustered\nSP = Spontaneous ETAS events\nuncert = 1-sigma uncertainty normalized by the mean\n");	
		}
		
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));

		File outputDir = new File(dirName+ "/FigureA");
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/FigureA/"+durationsYrs+"yr_LossExceed";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Loss ($)", "Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
		fileNamePrefix += "_zoomed";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Loss ($)", "Exceedance Probability", new Range(1e6,1e10), new Range(1e-2,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		UncertainArbDiscFunc denomCurve = funcList.get(0);
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		double xThresh = denomCurve.getClosestXtoY(1e-4);
		for(int i=0;i<funcList.size();i++) {
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			ratioCurve.setName("Ratio for "+funcList.get(i).getName());
			for(int j=0;j<denomCurve.size();j++) {
				if(denomCurve.getX(j)>xThresh)
					break;
				double ratio = funcList.get(i).getY(j)/denomCurve.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(denomCurve.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
		}
		fileNamePrefix = dirName+"/FigureA/"+durationsYrs+"yr_LossExceedRatio";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Loss ($)", "Exceedance Ratio", new Range(1e4,1e11), new Range(0.0,2.0), true, false, fileNamePrefix, true);

		// plot incremental distributions
		ArrayList<XY_DataSet> pdfList = new ArrayList<XY_DataSet>(); 
		for(int i=0;i<2;i++) { // Only the first two funcs
			pdfList.add(convertExceedToIncrProbDist(funcList.get(i)));
		}
		if(randCOV) {
			UncertainArbDiscFunc lossIncreProbDistFromTD_CatalogsRandFromCOV = getLossIncrProbDistFromCatalogsRandomFromCOV(getCatalogs(), 1d);
			lossIncreProbDistFromTD_CatalogsRandFromCOV.setName("lossIncreProbDistFromTD_CatalogsRandFromCOV");
			UncertainArbDiscFunc lossIcreProbDistFromPoisCatalogsRandFromCOV = getLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1d);
			lossIcreProbDistFromPoisCatalogsRandFromCOV.setName("lossIcreProbDistFromPoisCatalogsRandFromCOV");
			pdfList.add(lossIncreProbDistFromTD_CatalogsRandFromCOV);
			pdfList.add(lossIcreProbDistFromPoisCatalogsRandFromCOV);
		}
		ArrayList<PlotCurveCharacterstics> plotChars3 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.RED));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		fileNamePrefix = dirName+"/FigureA/"+durationsYrs+"yr_IncrLoss";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(pdfList, plotChars3, "", "Loss ($1)", "Incremental Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
	}
	
	
	
	/**
	 * This is for aggregate loss exceedances (only random COV sample supported)
	 * @param randCOV
	 */

	public void makeFigAndTableSetB_1yr_Exceedances() {
		
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		
		boolean randCOV = true;
		String randCOVincluded = "";
		if(randCOV)
			randCOVincluded = " (cond loss PDF for each rup randomly sampled)";

		UncertainArbDiscFunc aggrLossExceedProbFromTD_CatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1d);
		aggrLossExceedProbFromTD_CatalogsRandFromCOV.setName("aggrLossExceedProbFromTD_CatalogsRandFromCOV");
		UncertainArbDiscFunc aggrLossExceedProbFromPoisCatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1d);
		aggrLossExceedProbFromPoisCatalogsRandFromCOV.setName("aggrLossExceedProbFromPoisCatalogsRandFromCOV");
		UncertainArbDiscFunc aggrLossExceedProbFromDeclusteredCatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogsDeclustered(), 1d);
		aggrLossExceedProbFromDeclusteredCatalogsRandFromCOV.setName("aggrLossExceedProbFromDeclusteredCatalogsRandFromCOV");
		UncertainArbDiscFunc aggrLossExceedProbFromSpontaneousCatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getSpontaneousEventsCatalog(getCatalogs()), 1d);
		aggrLossExceedProbFromSpontaneousCatalogsRandFromCOV.setName("aggrLossExceedProbFromSpontaneousCatalogsRandFromCOV");
		funcList.add(aggrLossExceedProbFromTD_CatalogsRandFromCOV);
		funcList.add(aggrLossExceedProbFromPoisCatalogsRandFromCOV);
		funcList.add(aggrLossExceedProbFromDeclusteredCatalogsRandFromCOV);
		funcList.add(aggrLossExceedProbFromSpontaneousCatalogsRandFromCOV);			

		System.out.println("\n1yr Loss for CEA probability levels"+randCOVincluded+":\n");
		System.out.println("Probability\tTD_loss\tTD_uncert\tTI_loss\tTI_uncert\tGK_loss\tGK_uncert\tSP_loss\tSP_uncert");
		System.out.println(getTableStringOfCEA_Values(funcList,true));
		System.out.println("\n1yr Loss ratios, relative to TD, for CEA probability levels"+randCOVincluded+":\n");
		System.out.println("Probability\tTIvsTD_loss\tTIvsTD_uncert\tGKvsTD_loss\tGKvsTD_uncert\tSPvsTD_loss\tSPvsTD_uncert");
		System.out.println(getTableStringOfCEA_Ratios(funcList, true));
		System.out.println("TI = Time Dependent\nTI = Time Independent (Poisson)\nGK = Gardner Knopoff declustered\nSP = Spontaneous ETAS events\nuncert = 1-sigma uncertainty normalized by the mean\n");
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		File outputDir = new File(dirName+ "/FigureB");
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/FigureB/1.0yr_AggrLossExceed";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "", "Aggregate Loss ($)", "Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		UncertainArbDiscFunc denomCurve = funcList.get(0);
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		double xThresh = denomCurve.getClosestXtoY(1e-4);
		for(int i=0;i<funcList.size();i++) {
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			ratioCurve.setName("Ratio for "+funcList.get(i).getName());
			for(int j=0;j<denomCurve.size();j++) {
				if(denomCurve.getX(j)>xThresh)
					break;
				double ratio = funcList.get(i).getY(j)/denomCurve.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(denomCurve.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
		}
		fileNamePrefix = dirName+"/FigureB/1.0yr_AggrLossExceedRatio";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Aggregate Loss ($)", "Exceedance Ratio", new Range(1e4,1e11), new Range(0.0,2.0), true, false, fileNamePrefix, true);

		// plot incremental distributions
		ArrayList<XY_DataSet> pdfList = new ArrayList<XY_DataSet>(); 
		for(int i=0;i<2;i++) {
			pdfList.add(convertExceedToIncrProbDist(funcList.get(i)));
		}
		if(randCOV) {
			UncertainArbDiscFunc lossAggrIncreProbDistFromTD_CatalogsRandFromCOV = getAggrLossIncrProbDistFromCatalogsRandomFromCOV(getCatalogs(), 1d);
			lossAggrIncreProbDistFromTD_CatalogsRandFromCOV.setName("lossAggrIncreProbDistFromTD_CatalogsRandFromCOV");
			UncertainArbDiscFunc lossAggrIcreProbDistFromPoisCatalogsRandFromCOV = getAggrLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1d);
			lossAggrIcreProbDistFromPoisCatalogsRandFromCOV.setName("lossAggrIcreProbDistFromPoisCatalogsRandFromCOV");
			pdfList.add(lossAggrIncreProbDistFromTD_CatalogsRandFromCOV);
			pdfList.add(lossAggrIcreProbDistFromPoisCatalogsRandFromCOV);
		}
		ArrayList<PlotCurveCharacterstics> plotChars3 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.RED));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars3.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		fileNamePrefix = dirName+"/FigureB/1.0yr_AggrIncrLoss";
		if(randCOV) fileNamePrefix += "_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(pdfList, plotChars3, "", "Aggregate Loss ($)", "Incremental Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);
	}
	
	
	public void make1yr_StatsCSV_File() {
//		maxLoss, subCatIndexForMaxAggLoss: 8.253113091481961E10, 244544
//		maxAggLoss, subCatIndexForMaxAggLoss: 1.153828487953038E11, 210663
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(getCatalogs(), 1d);
//		ObsEqkRupList catalog = subCatList.get(244544);
//		MakeFigures.plotEventRateMapForCatalog(catalog, 2.5, "FigRightHere092624", 1d);
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("catindex", "numEvents", "numGT_6pt5", "maxRandLoss", "aggrRandLoss", "maxMag", "origCatIndex", "maxMeanLoss", "aggrMeanLoss");
		int catIndex = -1;
		for (ObsEqkRupList catalog:subCatList) {
			catIndex+=1;
			int origCatIndex = (int)((double)catIndex/500d);
			double maxRandLoss=0, aggrRandLoss=0, maxMeanLoss=0, aggrMeanLoss=0, nGT6pt5=0, maxMag=0;;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				if(rup.getMag()>=6.5)
					nGT6pt5+=1;
				if(rup.getMag()>=maxMag)
					maxMag=rup.getMag();
				double meanLoss = meanLossForNthRup.get(rup.getNthERF_Index());
				if(maxMeanLoss<meanLoss) maxMeanLoss=meanLoss;
				aggrMeanLoss+=meanLoss;

				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxRandLoss<randLoss) maxRandLoss=randLoss;
				aggrRandLoss+=randLoss;
			}
			csv.addLine(catIndex+"", catalog.size()+"", nGT6pt5+"", maxRandLoss+"", aggrRandLoss+"",maxMag+"",origCatIndex+"", maxMeanLoss+"", aggrMeanLoss+"");
		}
		File file = new File(dirName+"/oneYearStats.csv");
		try {
			csv.writeToFile(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void plotMaxAndAggrLossesPerYearInCatalog(int catalogIndex) {
		
		HistogramFunction ratioHistTD = getBlankIncrLossCurveLogX();
		HistogramFunction numHistTD = getBlankIncrLossCurveLogX();
		HistogramFunction ratioHistPoiss = getBlankIncrLossCurveLogX();
		HistogramFunction numHistPoiss = getBlankIncrLossCurveLogX();
		
		ArrayList<ObsEqkRupList> tempCatList = new ArrayList<ObsEqkRupList>();
		tempCatList.add(getCatalogs().get(catalogIndex));
		
		File outputDir = new File(dirName+ "/MaxAndAggLossVsTime");
		if(!outputDir.exists()) outputDir.mkdir();

		
		double aveAggrOverMax_TD=0;
		double aveAggrOverMax_Pois=0;
		double aveAggrOverMax_TD_gt1Billion=0;
		double aveAggrOverMax_Pois_gt1Billion=0;
		double aveAggrOverMax_TD_gt1Billion_num=0;
		double aveAggrOverMax_Pois_gt1Billion_num=0;
		
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(tempCatList, 1d);
		HistogramFunction histTD_MaxLoss = new HistogramFunction(0.05,500,1.0);
		HistogramFunction histTD_AggrLoss = new HistogramFunction(0.05,500,1.0);
		int subCatIndex=-1;
		for (ObsEqkRupList catalog:subCatList) {
			subCatIndex+=1;
			double maxLoss=0, aggrLoss=0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxLoss<randLoss) maxLoss=randLoss;
				aggrLoss+=randLoss;
			}
			if(maxLoss != 0.0) {
				aveAggrOverMax_TD += aggrLoss/maxLoss;
				if(maxLoss>1e9) {
					aveAggrOverMax_TD_gt1Billion+=aggrLoss/maxLoss;
					aveAggrOverMax_TD_gt1Billion_num+=1;
				}
			}
			else {
				aveAggrOverMax_TD += 1.0;
			}
			histTD_MaxLoss.set(subCatIndex,Math.log10(maxLoss));
			histTD_AggrLoss.set(subCatIndex,Math.log10(aggrLoss));
			if(Math.log10(maxLoss)>ratioHistTD.getMinX()-ratioHistTD.getDelta()/2) {
				ratioHistTD.add(Math.log10(maxLoss), aggrLoss/maxLoss);
				numHistTD.add(Math.log10(maxLoss), 1.0);
			}
		}
		aveAggrOverMax_TD_gt1Billion /= aveAggrOverMax_TD_gt1Billion_num;
		aveAggrOverMax_TD /= subCatList.size();
		histTD_AggrLoss.setInfo("aveAggrOverMax_TD="+aveAggrOverMax_TD+"\naveAggrOverMax_TD_gt1Billion="+aveAggrOverMax_TD_gt1Billion);
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
		funcList.add(histTD_AggrLoss);
		funcList.add(histTD_MaxLoss);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLUE));
		String fileNamePrefix = dirName+"/MaxAndAggLossVsTime/maxAndAggLossPerYearForTD_Cat_"+catalogIndex;
		Range xRange = new Range(0,500);
//		Range xRange = new Range(163-50,163+50);
		Range yRange = new Range(4,11.3);
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Full TD", "Year", "Log10(Loss, $)", xRange, yRange, false, false, 6.5, 2.0, fileNamePrefix, true);
		
		// Now do poisson result
		tempCatList = new ArrayList<ObsEqkRupList>();
		tempCatList.add(getRandomizedCatalogs().get(catalogIndex));
		subCatList = getSubcatalogList(tempCatList, 1d);
		HistogramFunction histPois_MaxLoss = new HistogramFunction(0.05,500,1.0);
		HistogramFunction histPois_AggrLoss = new HistogramFunction(0.05,500,1.0);
		subCatIndex=-1;
		for (ObsEqkRupList catalog:subCatList) {
			subCatIndex+=1;
			double maxLoss=0, aggrLoss=0;
			for (ObsEqkRupture obsRup : catalog) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(maxLoss<randLoss) maxLoss=randLoss;
				aggrLoss+=randLoss;
			}
			if(maxLoss != 0.0) {
				aveAggrOverMax_Pois += aggrLoss/maxLoss;
				if(maxLoss>1e9) {
					aveAggrOverMax_Pois_gt1Billion+=aggrLoss/maxLoss;
					aveAggrOverMax_Pois_gt1Billion_num+=1;
				}

			}
			else {
				aveAggrOverMax_Pois += 1.0;
			}

			histPois_MaxLoss.set(subCatIndex,Math.log10(maxLoss));
			histPois_AggrLoss.set(subCatIndex,Math.log10(aggrLoss));
			if(Math.log10(maxLoss)>ratioHistPoiss.getMinX()-ratioHistPoiss.getDelta()/2) {
				ratioHistPoiss.add(Math.log10(maxLoss), aggrLoss/maxLoss);
				numHistPoiss.add(Math.log10(maxLoss), 1.0);
			}

		}
		aveAggrOverMax_Pois_gt1Billion /= aveAggrOverMax_Pois_gt1Billion_num;

		aveAggrOverMax_Pois /= subCatList.size();
		histPois_AggrLoss.setInfo("aveAggrOverMax_Pois="+aveAggrOverMax_Pois+"\naveAggrOverMax_Pois_gt1Billion="+aveAggrOverMax_Pois_gt1Billion);
		ArrayList<XY_DataSet> funcListPoiss = new ArrayList<XY_DataSet>(); 
		funcListPoiss.add(histPois_AggrLoss);
		funcListPoiss.add(histPois_MaxLoss);
		fileNamePrefix = dirName+"/MaxAndAggLossVsTime/maxAndAggLossPerYearForTI_Cat_"+catalogIndex;
		PlottingUtils.writeAndOrPlotFuncs(funcListPoiss, plotChars, "Full Randomized (Poisson)", "Year", "Log10(Loss, $)", xRange, yRange, false, false, 6.5, 2.0, fileNamePrefix, true);

		// plot aggr/max loss histrograms as func of log10(maxLoss):
		for(int i=0;i<ratioHistPoiss.size(); i++) {
			double newVal = ratioHistTD.getY(i)/numHistTD.getY(i);
			ratioHistTD.set(i,newVal);
			newVal = ratioHistPoiss.getY(i)/numHistPoiss.getY(i);
			ratioHistPoiss.set(i,newVal);
		}
		ratioHistTD.setName("ratioHistTD");
		ratioHistPoiss.setName("ratioHistPoiss");
		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); 
		funcList2.add(ratioHistTD);
		funcList2.add(ratioHistPoiss);
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		fileNamePrefix = dirName+"/MaxAndAggLossVsTime/aveRatioVsMaxLoss"+catalogIndex;
		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars2, "", "Log10(MaxLoss ($))", "Ratio", null, null, false, false, fileNamePrefix, true);


			
	}
	
	
	public void plotAccumulatedLossVsTimeForCatalog(int catalogIndex) {
		
		ArbitrarilyDiscretizedFunc accumLossFunc_TD = new ArbitrarilyDiscretizedFunc();
		accumLossFunc_TD.setName("accumLossFunc_TD");
		double accumLoss_TD=0;
		accumLossFunc_TD.set(0.0,0.0);
		long lastMillis = 0;
		for (ObsEqkRupture obsRup : catalogList.get(catalogIndex)) {
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				long timeMillis = rup.getOriginTime()-startTimeMillis;
				if(timeMillis<lastMillis)
					throw new RuntimeException("rupture order proplem");
				lastMillis = timeMillis;
				double timeYears = (double)timeMillis/(double)millisPerYr;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(randLoss == 0.0)
					continue;
				accumLoss_TD += randLoss;
				accumLossFunc_TD.set(timeYears,accumLoss_TD);
		}
		
		
		ArbitrarilyDiscretizedFunc accumLossFunc_Pois = new ArbitrarilyDiscretizedFunc();
		accumLossFunc_Pois.setName("accumLossFunc_Pois");
		double accumLoss_Pois=0;
		accumLossFunc_Pois.set(0.0,0.0);
		lastMillis = 0;
		for (ObsEqkRupture obsRup : randomizeCatalog(catalogList.get(catalogIndex))) { // use special method so catalog has exact same events
				ETAS_EqkRupture rup = (ETAS_EqkRupture)obsRup;
				long timeMillis = rup.getOriginTime()-startTimeMillis;
				if(timeMillis<lastMillis)
					throw new RuntimeException("rupture order proplem");
				lastMillis = timeMillis;
				double timeYears = (double)timeMillis/(double)millisPerYr;
				double randLoss = randomLossForEventID[rup.getID()]; // covModel.getDistribution(aveLoss).sample(); // get randome sample
				if(randLoss == 0.0)
					continue;
				accumLoss_Pois += randLoss;
				accumLossFunc_Pois.set(timeYears,accumLoss_Pois);
		}
		
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
		funcList.add(accumLossFunc_TD);
		funcList.add(accumLossFunc_Pois);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		String fileNamePrefix = dirName+"/JUNK";
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "", "Year", "Loss", new Range(0,500), null, false, false, null, true);
		
	}


	
	public void makeFigSetC_500yrRates(int catIndex) {
		double startTimeYears = 0; // this is the time within the catalog to start, will need to set lower than this as our simulations are 500 years long
		double durationYears = 500;
		int plotWidthPixels = 2000;
		double binWidthDays = 365.25;
		String prefix = null;
		double annotateMinMag = 6;
		double annotateBinWidth = 0.5;
		// use this because it has <≥2.5 events
//		String fullCatFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/testDebugRun/2020_12_03-Start2012_500yr_kCOV1p5_ScaleFactor1p4_Spontaneous_HistoricalCatalog/results_complete.bin"; 
//		ObsEqkRupList rupList=null;
//		try {
//			rupList = U3ETAS_SimulationAnalysis.loadCatalogs(new File(U3ETAS_SimulationAnalysis.fssFileName), new File(fullCatFileName), 2.5).get(0);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		ArrayList<ETAS_EqkRupture> catalog = new ArrayList<ETAS_EqkRupture>();
		for(ObsEqkRupture rup:this.getCatalogs().get(catIndex))
			catalog.add((ETAS_EqkRupture)rup);
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure1";
		File outputDir = new File(dirName+ "/FigureC");
		if(!outputDir.exists()) outputDir.mkdir();
		try {
			ETAS_SimAnalysisTools.plotRateOverTime(catalog, startTimeYears, durationYears, binWidthDays, outputDir, prefix, plotWidthPixels,
					annotateMinMag, annotateBinWidth);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	
	public void makeFig1_100yrRates(int catIndex, double startTimeYears) {
		if(catIndex>999)
			throw new RuntimeException("Only the first set of catalogs is supported");
//		double startTimeYears = 163-50; // time within the catalog to start, this time is 50 yrs before the largest aggregate loss
		double durationYears = 100;
		int plotWidthPixels = 2000;
		double binWidthDays = 30;
		String prefix = null;
		double annotateMinMag = 6;
		double annotateBinWidth = 0.5;
		// use this because it has <≥2.5 events
		String fullCatFileName = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/CEA_WGCEP/UCERF3/DeclusteringOnLossAnalysisPaper2024/Data/results_complete.bin";
//		ObsEqkRupList rupList=null;
		ETAS_Catalog rupList=null;
		try {
			rupList = ETAS_CatalogIO.loadIndividualCatalogBinary(new File(fullCatFileName), catIndex, true);
//			rupList = U3ETAS_SimulationAnalysis.loadCatalogs(new File(U3ETAS_SimulationAnalysis.fssFileName), new File(fullCatFileName), 2.5).get(0);
		} catch (Exception e) {
			e.printStackTrace();
		}
		ArrayList<ETAS_EqkRupture> catalog = new ArrayList<ETAS_EqkRupture>();
		for(ObsEqkRupture rup:rupList)
			catalog.add((ETAS_EqkRupture)rup);
		File outputDir = new File(dirName+ "/Figure1_cat"+catIndex);
		if(!outputDir.exists()) outputDir.mkdir();
		try {
			ETAS_SimAnalysisTools.plotRateOverTime(catalog, startTimeYears, durationYears, binWidthDays, outputDir, prefix, plotWidthPixels,
					annotateMinMag, annotateBinWidth);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	/**
	 * This is for one or more exceedances, randCOV, and for multiple return forecast durations
	 * @param randCOV
	 */
	public void makeFigSetD_Exceedances() throws IOException {
				
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		float[] durationsYrsArray = {1f/365.25f, 7f/365.25f, 1f/12, 1f, 5f, 50f}; // day, week, month, yr, 5yr, 50yr
//		float[] durationsYrsArray = {1f}; // day, week, month, yr, 5yr, 50yr
		
		File outputDir = new File(dirName+ "/FigureD");
		if(!outputDir.exists()) outputDir.mkdir();
		String dataDirName = dirName+ "/FigureD/data";
		File dataDir = new File(dataDirName);
		if(!dataDir.exists()) dataDir.mkdir();
		Gson gson = new GsonBuilder().setPrettyPrinting().create();

		for(float durationsYrs:durationsYrsArray) {
			System.out.println("working on duration (yrs) = "+durationsYrs);
			
			UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV;
			UncertainArbDiscFunc lossExceedProbFromPoisCatalogsRandFromCOV;
			
			String fileName = "TD_Func_"+Float.toString(durationsYrs)+".json";
			File outFile = new File(dataDirName+"/"+fileName);
			String fileName2 = "TI_Func_"+Float.toString(durationsYrs)+".json";
			File outFile2 = new File(dataDirName+"/"+fileName2);
			if(outFile.exists()) {
				FileReader fr = new FileReader(outFile);
				FileReader fr2 = new FileReader(outFile2);
				lossExceedProbFromTD_CatalogsRandFromCOV = gson.fromJson(fr, UncertainArbDiscFunc.class);
				lossExceedProbFromPoisCatalogsRandFromCOV = gson.fromJson(fr2, UncertainArbDiscFunc.class);
			}
			else {
				lossExceedProbFromTD_CatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), durationsYrs);
				lossExceedProbFromTD_CatalogsRandFromCOV.setName("lossExceedProbFromTD_CatalogsRandFromCOV; duration = "+durationsYrs);
				lossExceedProbFromPoisCatalogsRandFromCOV = getLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), durationsYrs);
				lossExceedProbFromPoisCatalogsRandFromCOV.setName("lossExceedProbFromPoisCatalogsRandFromCOV; duration = "+durationsYrs);				

				FileWriter fw = new FileWriter(outFile);
				gson.toJson(lossExceedProbFromTD_CatalogsRandFromCOV, fw);
				FileWriter fw2 = new FileWriter(outFile2);
				gson.toJson(lossExceedProbFromPoisCatalogsRandFromCOV, fw2);
				fw.close();
				fw2.close();
			}
			
			funcList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
			funcList.add(lossExceedProbFromPoisCatalogsRandFromCOV);
			
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			double xThresh = lossExceedProbFromTD_CatalogsRandFromCOV.getClosestXtoY(1e-4);

			ratioCurve.setName("Ratio for "+lossExceedProbFromPoisCatalogsRandFromCOV.getName());
			for(int j=0;j<lossExceedProbFromTD_CatalogsRandFromCOV.size();j++) {
				if(lossExceedProbFromTD_CatalogsRandFromCOV.getX(j)>xThresh)
					break;
				double ratio = lossExceedProbFromPoisCatalogsRandFromCOV.getY(j)/lossExceedProbFromTD_CatalogsRandFromCOV.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(lossExceedProbFromTD_CatalogsRandFromCOV.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
						
		}
		
		
		
//		// Write tables
//		if(durationsYrs==1.0) { // only for 1-year forecast durations
//			System.out.println("\n"+durationsYrs+"yr Loss for CEA probability levels"+randCOVincluded+":\n");
//			System.out.println("Probability\tTD_loss\tTD_uncert\tTI_loss\tTI_uncert\tGK_loss\tGK_uncert\tSP_loss\tSP_uncert");
//			System.out.println(getTableStringOfCEA_Values(funcList,true));
//			System.out.println("\n"+durationsYrs+"yr Loss ratios, relative to TD, for CEA probability levels"+randCOVincluded+":\n");
//			System.out.println("Probability\tTIvsTD_loss\tTIvsTD_uncert\tGKvsTD_loss\tGKvsTD_uncert\tSPvsTD_loss\tSPvsTD_uncert");
//			System.out.println(getTableStringOfCEA_Ratios(funcList, true));
//			System.out.println("TI = Time Dependent\nTI = Time Independent (Poisson)\nGK = Gardner Knopoff declustered\nSP = Spontaneous ETAS events\nuncert = 1-sigma uncertainty normalized by the mean\n");	
//		}
		
		ArrayList<PlotCurveCharacterstics> plotChars1 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));

		String fileNamePrefix = dirName+"/FigureD/LossExceed_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars1, "", "Loss ($)", "Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.MAGENTA));
		fileNamePrefix = dirName+"/FigureD/LossExceedRatio_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Loss ($)", "Exceedance Ratio", new Range(1e4,1e11), new Range(0.5,1.5), true, false, fileNamePrefix, true);

	}
	
	
	/**
	 * This is for aggregate exceedances, randCOV, and for multiple return forecast durations
	 * @param randCOV
	 */
	public void makeFigSetE_Exceedances() throws IOException {
				
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		float[] durationsYrsArray = {1f/365.25f, 7f/365.25f, 1f/12, 1f, 5f, 50f}; // day, week, month, yr, 5yr, 50yr
//		float[] durationsYrsArray = {1f}; // day, week, month, yr, 5yr, 50yr
		
		File outputDir = new File(dirName+ "/FigureE");
		if(!outputDir.exists()) outputDir.mkdir();
		String dataDirName = dirName+ "/FigureE/data";
		File dataDir = new File(dataDirName);
		if(!dataDir.exists()) dataDir.mkdir();
		Gson gson = new GsonBuilder().setPrettyPrinting().create();

		for(float durationsYrs:durationsYrsArray) {
			System.out.println("working on duration (yrs) = "+durationsYrs);
			
			UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV;
			UncertainArbDiscFunc lossExceedProbFromPoisCatalogsRandFromCOV;
			
			String fileName = "TD_AggrFunc_"+Float.toString(durationsYrs)+".json";
			File outFile = new File(dataDirName+"/"+fileName);
			String fileName2 = "TI_AggrFunc_"+Float.toString(durationsYrs)+".json";
			File outFile2 = new File(dataDirName+"/"+fileName2);
			if(outFile.exists()) {
				FileReader fr = new FileReader(outFile);
				FileReader fr2 = new FileReader(outFile2);
				lossExceedProbFromTD_CatalogsRandFromCOV = gson.fromJson(fr, UncertainArbDiscFunc.class);
				lossExceedProbFromPoisCatalogsRandFromCOV = gson.fromJson(fr2, UncertainArbDiscFunc.class);
			}
			else {
				lossExceedProbFromTD_CatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), durationsYrs);
				lossExceedProbFromTD_CatalogsRandFromCOV.setName("aggrLossExceedProbFromTD_CatalogsRandFromCOV; duration = "+durationsYrs);
				lossExceedProbFromPoisCatalogsRandFromCOV = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), durationsYrs);
				lossExceedProbFromPoisCatalogsRandFromCOV.setName("aggrLossExceedProbFromPoisCatalogsRandFromCOV; duration = "+durationsYrs);				

				FileWriter fw = new FileWriter(outFile);
				gson.toJson(lossExceedProbFromTD_CatalogsRandFromCOV, fw);
				FileWriter fw2 = new FileWriter(outFile2);
				gson.toJson(lossExceedProbFromPoisCatalogsRandFromCOV, fw2);
				fw.close();
				fw2.close();
			}
			
			funcList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
			funcList.add(lossExceedProbFromPoisCatalogsRandFromCOV);

			
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			double xThresh = lossExceedProbFromTD_CatalogsRandFromCOV.getClosestXtoY(1e-4);

			ratioCurve.setName("Ratio for "+lossExceedProbFromPoisCatalogsRandFromCOV.getName());
			for(int j=0;j<lossExceedProbFromTD_CatalogsRandFromCOV.size();j++) {
				if(lossExceedProbFromTD_CatalogsRandFromCOV.getX(j)>xThresh)
					break;
				double ratio = lossExceedProbFromPoisCatalogsRandFromCOV.getY(j)/lossExceedProbFromTD_CatalogsRandFromCOV.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(lossExceedProbFromTD_CatalogsRandFromCOV.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
						
		}
		
		ArrayList<PlotCurveCharacterstics> plotChars1 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		String fileNamePrefix = dirName+"/FigureE/AggrLossExceed_randCOV";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars1, "", "Loss ($)", "Aggregate Exceedance Probability", new Range(1e4,1e11), new Range(1e-5,1.0), true, true, fileNamePrefix, true);

		// Make ratio plots
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.MAGENTA));
		fileNamePrefix = dirName+"/FigureE/AggrLossExceedRatio_randCOV";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Loss ($)", "Aggregate Exceedance Ratio", new Range(1e4,1e11), new Range(0.5,1.5), true, false, fileNamePrefix, true);

	}
	
	
	
	
	/**
	 * This is for aggregate exceedances, randCOV, and for multiple return forecast durations
	 * @param randCOV
	 */
	public void makeFigSetG_MagExceedances() throws IOException {
				
		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
		ArrayList<XY_DataSet> funcListRatios = new ArrayList<XY_DataSet>(); 

		float[] durationsYrsArray = {1f/365.25f, 7f/365.25f, 1f/12, 1f, 5f, 50f}; // day, week, month, yr, 5yr, 50yr
//		float[] durationsYrsArray = {1f}; // day, week, month, yr, 5yr, 50yr
		
		File outputDir = new File(dirName+ "/FigureG");
		if(!outputDir.exists()) outputDir.mkdir();
		String dataDirName = dirName+ "/FigureG/data";
		File dataDir = new File(dataDirName);
		if(!dataDir.exists()) dataDir.mkdir();
		Gson gson = new GsonBuilder().setPrettyPrinting().create();

		for(float durationsYrs:durationsYrsArray) {
			System.out.println("working on duration (yrs) = "+durationsYrs);
			
			UncertainArbDiscFunc magExceedProbFromTD_Catalogs;
			UncertainArbDiscFunc magExceedProbFromPoisCatalogs;
			
			String fileName = "TD_MagExceedFunc_"+Float.toString(durationsYrs)+".json";
			File outFile = new File(dataDirName+"/"+fileName);
			String fileName2 = "TI_MagExceedFunc_"+Float.toString(durationsYrs)+".json";
			File outFile2 = new File(dataDirName+"/"+fileName2);
			if(outFile.exists()) {
				FileReader fr = new FileReader(outFile);
				FileReader fr2 = new FileReader(outFile2);
				magExceedProbFromTD_Catalogs = gson.fromJson(fr, UncertainArbDiscFunc.class);
				magExceedProbFromPoisCatalogs = gson.fromJson(fr2, UncertainArbDiscFunc.class);
			}
			else {
				magExceedProbFromTD_Catalogs = getMagExceedProbCurveFromCatalogs(getCatalogs(), durationsYrs);
				magExceedProbFromTD_Catalogs.setName("magExceedProbFromTD_Catalogs; duration = "+durationsYrs);
				magExceedProbFromPoisCatalogs = getMagExceedProbCurveFromCatalogs(getRandomizedCatalogs(), durationsYrs);
				magExceedProbFromPoisCatalogs.setName("magExceedProbFromPoisCatalogs; duration = "+durationsYrs);				

				FileWriter fw = new FileWriter(outFile);
				gson.toJson(magExceedProbFromTD_Catalogs, fw);
				FileWriter fw2 = new FileWriter(outFile2);
				gson.toJson(magExceedProbFromPoisCatalogs, fw2);
				fw.close();
				fw2.close();
			}
			
			funcList.add(magExceedProbFromTD_Catalogs);
			funcList.add(magExceedProbFromPoisCatalogs);

			
			ArbitrarilyDiscretizedFunc ratioCurve = new ArbitrarilyDiscretizedFunc();
			double xThresh = magExceedProbFromTD_Catalogs.getClosestXtoY(1e-4);

			ratioCurve.setName("Ratio for "+magExceedProbFromPoisCatalogs.getName());
			for(int j=0;j<magExceedProbFromTD_Catalogs.size();j++) {
				if(magExceedProbFromTD_Catalogs.getX(j)>xThresh)
					break;
				double ratio = magExceedProbFromPoisCatalogs.getY(j)/magExceedProbFromTD_Catalogs.getY(j);
				if(ratio<1e-9 || ratio>1e9) ratio = 1e-9;
				ratioCurve.set(magExceedProbFromTD_Catalogs.getX(j), ratio);
			}
			funcListRatios.add(ratioCurve);
						
		}
		
		ArrayList<PlotCurveCharacterstics> plotChars1 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars1.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		String fileNamePrefix = dirName+"/FigureG/MagExceed";
		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars1, "", "Magnitude", "Exceedance Probability", new Range(5.0,8.5), new Range(1e-7,1), false, true, fileNamePrefix, true);

		// Make ratio plots
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 1f, Color.MAGENTA));
		fileNamePrefix = dirName+"/FigureG/MagExceedRatio";
		PlottingUtils.writeAndOrPlotFuncs(funcListRatios, plotChars2, "", "Magnitude", "Exceedance Ratio", new Range(5.0,8.5), new Range(0.5,1.5), false, false, fileNamePrefix, true);

	}



	public void makeFigSetF_N_ExceedHists() {
		// For N event histograms:
		File outputDir = new File(dirName+ "/FigureF");
		if(!outputDir.exists()) outputDir.mkdir();
		double[] lossExceedArray = {1e04, 4e7, 1.23e10};
		for(double loss:lossExceedArray) {
			HistogramFunction histN_Dist_TD = getNumLossExceedDistFromCatalogsRandomFromCOV(getCatalogs(), 1.0, loss);
			histN_Dist_TD.setName("TD N exceed dist for loss = "+loss);
			HistogramFunction histN_Dist_TI = getNumLossExceedDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0, loss);
			histN_Dist_TI.setName("TI N exceed dist for loss = "+loss);
			
//			// I need a getRateExceedLoss(loss) method
//			PoissonDistribution poisDist = new PoissonDistribution(expNum);
//			for(int i=0;i<poissDist.size();i++) {
//				poissDist.set(i,poisDist.probability(i));
//			}

			ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
			funcList.add(MakeFigures.getStairStepFunction(histN_Dist_TD));
			funcList.add(MakeFigures.getStairStepFunction(histN_Dist_TI));
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			String fileNamePrefix = dirName+"/FigureF/N_ExceedForLossDist_"+loss;
			PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "", "N Exceedances for loss≥"+loss, "Fraction", null, null, false, true, fileNamePrefix, true);		
			fileNamePrefix = dirName+"/FigureF/N_ExceedForLossDist_"+loss+"_zoomed";
			PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "", "N Exceedances for loss≥"+loss, "Fraction", new Range(-0.5,40), new Range(3e-6,1), false, true, fileNamePrefix, true);		
		}

	}

	
	public void write1yrLargeEvents(int subCatalogIndex) {
//		try {
//			catalogList = loadCatalogs(new File(catalogsFileName), 2.5, false); // this catalog does not have anything below M 5
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (DocumentException e) {
//			e.printStackTrace();
//		}
//		MakeFigures.plotEventRateMapForCatalog(catalogList.get(421), 2.5, "testMapRightHere", 500);

		
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catalogList, 1d);
//		MakeFigures.plotEventRateMapForCatalog(subCatList.get(210663), 2.5, "testMapRightHere", 1d);
		for(ObsEqkRupture rup : subCatList.get(subCatalogIndex)) {
			int fssIndex = ((ETAS_EqkRupture)rup).getFSSIndex();
			int gen = ((ETAS_EqkRupture)rup).getGeneration();
			int parID = ((ETAS_EqkRupture)rup).getParentID();
			int rupID = Integer.parseInt(rup.getEventId()); // orig ID here because in ID replaced in this analysis
// the following event is not found so it's probably M≤5 (filtered out)
if(subCatalogIndex==210663 && rupID == 77552) { // this is the parent of the first M≥6 LA area event in subCatalogIndex=210663 (it started everything)
	GregorianCalendar cal = rup.getOriginTimeCal();
	String calString = "ORIGIN TIME: Yr="+cal.get(cal.YEAR)+", DayOfYr="+cal.get(cal.DAY_OF_YEAR)+", Hr="+cal.get(cal.HOUR_OF_DAY)+", Min="+cal.get(cal.MINUTE)+", Sec="+cal.get(cal.SECOND);
	System.out.println("HERE "+rupID+"\t"+gen+"\t"+parID+"\t"+(float)rup.getMag()+"\t"+calString+"\t"+fssIndex);
}
			if(fssIndex>=0) {
				int srcIndex = erf.getSrcIndexForFltSysRup(fssIndex);
				ProbEqkSource fssSrc = erf.getSource(srcIndex);
				GregorianCalendar cal = rup.getOriginTimeCal();
				String calString = "ORIGIN TIME: Yr="+cal.get(cal.YEAR)+", DayOfYr="+cal.get(cal.DAY_OF_YEAR)+", Hr="+cal.get(cal.HOUR_OF_DAY)+", Min="+cal.get(cal.MINUTE)+", Sec="+cal.get(cal.SECOND);
				System.out.println(rupID+"\t"+gen+"\t"+parID+"\t"+(float)rup.getMag()+"\t"+calString+"\t"+fssSrc.getName());
			}
			else if(rup.getMag()>=6.0) {
				GregorianCalendar cal = rup.getOriginTimeCal();
				String calString = "ORIGIN TIME: Yr="+cal.get(cal.YEAR)+", DayOfYr="+cal.get(cal.DAY_OF_YEAR)+", Hr="+cal.get(cal.HOUR_OF_DAY)+", Min="+cal.get(cal.MINUTE)+", Sec="+cal.get(cal.SECOND);
				System.out.println(rupID+"\t"+gen+"\t"+parID+"\t"+(float)rup.getMag()+"\t"+calString+"\tHypocenter: "+rup.getHypocenterLocation());
			}
		}
	}

	
	
	public void plotCatalogMap(int subCatalogIndex, Region region, boolean includeGriddedSeis) {
		
		File outputDir = new File(dirName+ "/CatalogMaps");
		if(!outputDir.exists()) outputDir.mkdir();

		List<? extends FaultSection> sectionList = erf.getSolution().getRupSet().getFaultSectionDataList();
		
		ArrayList<Double> magForSectionList = new ArrayList<Double>();
		ArrayList<Double> genForSectionList = new ArrayList<Double>();
		for(int s=0;s<sectionList.size();s++) {
			magForSectionList.add(Double.NaN);
			genForSectionList.add(Double.NaN);
		}
		
		CPT cptGen = new CPT();
		cptGen.setBelowMinColor(Color.WHITE);
		cptGen.setNanColor(Color.WHITE);
		cptGen.add(new CPTVal(0f, Color.BLACK, 1f, Color.BLUE));
		cptGen.add(new CPTVal(1f, Color.BLUE, 2f, Color.GREEN));
		cptGen.add(new CPTVal(2f, Color.GREEN, 3f, Color.ORANGE));
		cptGen.add(new CPTVal(3f, Color.ORANGE, 4f, Color.RED));
		cptGen.add(new CPTVal(4f, Color.RED, 5f, Color.MAGENTA));
		cptGen.setAboveMaxColor(Color.MAGENTA);
		System.out.println("cpt.getMaxValue()="+cptGen.getMaxValue());
		System.out.println("cpt.getMinValue()="+cptGen.getMinValue());
		
		CPT cptOrig=null;
		try {
			cptOrig = GMT_CPT_Files.RAINBOW_UNIFORM.instance(); //   RAINBOW_UNIFORM  SEQUENTIAL_BATLOW_UNIFORM
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		CPT cptMag = cptOrig.rescale(5d, 8.5d);
		cptMag.setBelowMinColor(Color.WHITE);
		cptMag.setNanColor(Color.WHITE);

		
		List<Location> hypoList = new ArrayList<Location>();
		List<PlotCurveCharacterstics> charsMag = new  ArrayList<PlotCurveCharacterstics>();
		List<PlotCurveCharacterstics> charsGen = new  ArrayList<PlotCurveCharacterstics>();
		PlotSymbol outlineSymbol = PlotSymbol.CIRCLE;
		String labelMag = "Magnitude";
		String labelGen = "Generation";

		double mMax=0;
		double mMax_gridded=0;
		ArrayList<ObsEqkRupList> subCatList = getSubcatalogList(catalogList, 1d);
		for(ObsEqkRupture rup : subCatList.get(subCatalogIndex)) {
			if(rup.getMag()>mMax) 
				mMax = rup.getMag();
			int fssIndex = ((ETAS_EqkRupture)rup).getFSSIndex();
			if(fssIndex>=0) {
//				if(rup.getMag()>7d) continue;
				int gen = ((ETAS_EqkRupture)rup).getGeneration();
				for(int s:erf.getSolution().getRupSet().getSectionsIndicesForRup(fssIndex)) {
					if(!Double.isNaN(magForSectionList.get(s))) {
						System.out.println("WARNING: Subsection # "+s+" ruptured more than once; previous mag="+magForSectionList.get(s));
						// continue; // this keeps the original
					}
					magForSectionList.set(s, rup.getMag());
					genForSectionList.set(s, (double)gen);
				}				
				int srcIndex = erf.getSrcIndexForFltSysRup(fssIndex);
				ProbEqkSource fssSrc = erf.getSource(srcIndex);
				GregorianCalendar cal = rup.getOriginTimeCal();
				String calString = "ORIGIN TIME: Yr="+cal.get(cal.YEAR)+", DayOfYr="+cal.get(cal.DAY_OF_YEAR)+", Hr="+cal.get(cal.HOUR_OF_DAY)+", Min="+cal.get(cal.MINUTE)+", Sec="+cal.get(cal.SECOND);
				System.out.println(rup.getMag()+"\t"+gen+"\t"+calString+"\t"+fssSrc.getName());

			}
			else { 
				hypoList.add(rup.getHypocenterLocation());  // everything less that 0.5f is the same size (0.7 is larger)
				float size = 0.5f+10f*((float)rup.getMag()-5.0f);
				charsMag.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, size, cptMag.getColor((float)rup.getMag())));  //  (float)rup.getMag()-5f
				int gen = ((ETAS_EqkRupture)rup).getGeneration();
				charsGen.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, size, cptGen.getColor((float)gen)));  //  (float)rup.getMag()-5f
				if(rup.getMag()>mMax_gridded) 
					mMax_gridded = rup.getMag();

			}
		}
		
		System.out.println("mMax="+mMax+"\nmMax_gridded="+mMax_gridded); // 7.76
		
		// make map with color based on mag
		GeographicMapMaker mapMaker = new GeographicMapMaker(sectionList);
		if(region != null)
			mapMaker.setRegion(region);
		mapMaker.setSkipNaNs(true);
		mapMaker.setScalarThickness(8f);
		if(includeGriddedSeis)
			mapMaker.plotScatters(hypoList, charsMag, outlineSymbol, cptMag, labelMag);
		mapMaker.plotSectScalars(magForSectionList, cptMag, null);
		// write mag for each section
		for(int s=0;s<sectionList.size();s++) {
			if(!Double.isNaN(magForSectionList.get(s)))
				System.out.println("here "+s+"\t"+magForSectionList.get(s)+"\t"+sectionList.get(s).getSectionName());
		}
		mapMaker.setFillSurfaces(false);
		mapMaker.setSectOutlineChar(null);
		try {
			String name = "Catalog_"+subCatalogIndex+"_Map_Mag";
			if(!includeGriddedSeis)
				name += "_noGridded";
			mapMaker.plot(outputDir, name, "");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// make map with color based on generation
		GeographicMapMaker mapMakerGen = new GeographicMapMaker(sectionList);
		if(region != null)
			mapMakerGen.setRegion(region);
		mapMakerGen.setSkipNaNs(true);
		mapMakerGen.setScalarThickness(8f);
		if(includeGriddedSeis)
			mapMakerGen.plotScatters(hypoList, charsGen, outlineSymbol, cptGen, labelGen);
		mapMakerGen.plotSectScalars(genForSectionList, cptGen, null);
		// write gen for each section
//		for(int s=0;s<sectionList.size();s++) {
//			if(!Double.isNaN(genForSectionList.get(s)))
//				System.out.println("here "+s+"\t"+genForSectionList.get(s)+"\t"+sectionList.get(s).getSectionName());
//		}
		mapMakerGen.setFillSurfaces(false);
		mapMakerGen.setSectOutlineChar(null);
		try {
			String name = "Catalog_"+subCatalogIndex+"_Map_Gen";
			if(!includeGriddedSeis)
				name += "_noGridded";
			mapMakerGen.plot(outputDir, name, "");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	
	
	
	
	/**
	 * This stores a random loss sample for each event (so we can have the same sampled value across different methods)
	 */
	public void setRandomLossForEvents() {
		randomLossForEventID = new double[totalNumEvents];
		for (ObsEqkRupList catalog : catalogList) {
			for (ObsEqkRupture rup : catalog) {
				ETAS_EqkRupture etasRup = (ETAS_EqkRupture)rup;
				double aveLoss = meanLossForNthRup.get(etasRup.getNthERF_Index());
				if(aveLoss==0)
					randomLossForEventID[etasRup.getID()] = 0;
				else {
					double randLoss = covModel.getDistribution(aveLoss, randomGen).sample(); // get random sample
					randomLossForEventID[etasRup.getID()] = randLoss;
				}
			}
		}
	}
	
	
	/**
	 * 
	 */
	public void testIncrVsExceedCurveTanslations() {

		UncertainArbDiscFunc incrFuncTI = getLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testExceedTI = convertIncrToExceedProbDist(incrFuncTI);
		UncertainArbDiscFunc exceedFuncTI = getLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testIncrFuncTI = convertExceedToIncrProbDist(exceedFuncTI);
		ArbitrarilyDiscretizedFunc testIncrDistFromRates = testLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		
		UncertainArbDiscFunc incrFuncTD = getLossIncrProbDistFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testExceedTD = convertIncrToExceedProbDist(incrFuncTD);
		UncertainArbDiscFunc exceedFuncTD = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testIncrFuncTD = convertExceedToIncrProbDist(exceedFuncTD);

		exceedFuncTI.setName("TI Exceed");
		testExceedTI.setName("TI Exceed converted from Incremental");
		incrFuncTI.setName("TI Incremental");
		testIncrFuncTI.setName("TI Incremental converted from Exceed");
		testIncrDistFromRates.setName("testIncrDistFromRates");
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
		funcList.add(exceedFuncTI);
		funcList.add(testExceedTI);
		funcList.add(incrFuncTI);
		funcList.add(testIncrFuncTI);
		funcList.add(testIncrDistFromRates);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.X, 2f, Color.MAGENTA));
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "TI Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);

		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); 

		exceedFuncTD.setName("TD Exceed");
		testExceedTD.setName("TD Exceed converted from Incremental");
		incrFuncTD.setName("TD Incremental");
		testIncrFuncTD.setName("TD Incremental converted from Exceed");

		funcList2.add(exceedFuncTD);
		funcList2.add(testExceedTD);
		funcList2.add(incrFuncTD);
		funcList2.add(testIncrFuncTD);
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars2, "TD Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);
	}
	
	/**
	 * 
	 */
	public void testAggrIncrVsExceedCurveTanslations() {

		UncertainArbDiscFunc incrFuncTI = getAggrLossIncrProbDistFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testExceedTI = convertAggIncrToExceedProbDist(incrFuncTI);
		UncertainArbDiscFunc exceedFuncTI = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testIncrFuncTI = convertAggrExceedToIncrProbDist(exceedFuncTI);
		
		UncertainArbDiscFunc incrFuncTD = getAggrLossIncrProbDistFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testExceedTD = convertAggIncrToExceedProbDist(incrFuncTD);
		UncertainArbDiscFunc exceedFuncTD = getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		ArbitrarilyDiscretizedFunc testIncrFuncTD = convertAggrExceedToIncrProbDist(exceedFuncTD);

		exceedFuncTI.setName("TI Aggr Exceed");
		testExceedTI.setName("TI Aggr Exceed converted from Incremental");
		incrFuncTI.setName("TI Aggr Incremental");
		testIncrFuncTI.setName("TI Aggr Incremental converted from Exceed");
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
		funcList.add(exceedFuncTI);
		funcList.add(testExceedTI);
		funcList.add(incrFuncTI);
		funcList.add(testIncrFuncTI);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "TI Aggr Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);

		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); 

		exceedFuncTD.setName("TD Aggr Exceed");
		testExceedTD.setName("TD Aggr Exceed converted from Incremental");
		incrFuncTD.setName("TD Aggr Incremental");
		testIncrFuncTD.setName("TD Aggr Incremental converted from Exceed");

		funcList2.add(exceedFuncTD);
		funcList2.add(testExceedTD);
		funcList2.add(incrFuncTD);
		funcList2.add(testIncrFuncTD);
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars2.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars2, "TD Aggr Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);
	}

	
	/**
	 * This tests whether averaging non-exceedcance probabilities from catalogs reaches a different answer
	 * than averaging exceedcance probabilities.  Both TI and TD match here
	 */
	public void testAveragingNonExceedProbs() {
		UncertainArbDiscFunc exceedFuncTI = getLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		UncertainArbDiscFunc exceedFuncTI_Alt = testGetLossExceedProbCurveFromCatalogsRandomFromCOV(getRandomizedCatalogs(), 1.0);
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>(); 
		funcList.add(exceedFuncTI);
		funcList.add(exceedFuncTI_Alt);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		PlottingUtils.writeAndOrPlotFuncs(funcList, plotChars, "Alt Exceed Test TI", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);

		UncertainArbDiscFunc exceedFuncTD = getLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		UncertainArbDiscFunc exceedFuncTD_Alt = testGetLossExceedProbCurveFromCatalogsRandomFromCOV(getCatalogs(), 1.0);
		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); 
		funcList2.add(exceedFuncTD);
		funcList2.add(exceedFuncTD_Alt);
		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars, "Alt Exceed Test TD", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, null, true);
	}


	
	
	
	
	public static void main(String[] args) throws IOException, DocumentException {
		
		
		/**
		 * NOTES:
		 * 
		 * Modal and Median loss will always be zero if we go down to very small earthquakes (or at least it's
		 * always strongly magnitude dependent).
		 * 
		 * AAL is an aggregate metric (summing all losses)
		 * 
		 * 
		 * 
		 */
				
		// make results reproducible
		int seed = 123471;
		U3ETAS_LossSimulationAnalysis analysis = new U3ETAS_LossSimulationAnalysis(seed);
		
		// this shows the MagProbDist discrepancy is less for 1 day than 1 wk or 1 month, as found for losses
		analysis.makeFigSetG_MagExceedances();
		
		// this examines  conditional loss exceedance curves for the year 2024
//		analysis.make1yrConditionalLossExceedComparisonFig();
//		analysis.make1yrConditionalAggrLossExceedComparisonFig();
		
		
		// This is to find interesting catalogs
//		analysis.make1yr_StatsCSV_File();
		// largest aggregate loss (for both rand and ave loss) is for subcatalog 216231 in main catalog=432
		
		// Map & info for above highest aggr loss
//		Location loc1 = new Location(33d,-120d,0d);
//		Location loc2 = new Location(35d,-117d,0d);
//		Region region = new Region(loc1,loc2);
//		analysis.plotCatalogMap(216231, region, true);
//		analysis.plotCatalogMap(216231, region, false);
//		analysis.write1yrLargeEvents(216231);
		// last method implies the sequences starts on year 2243, which is 2243-2012 = 231 years into the main catalog.
		
		// no popup window here
//		analysis.makeFig1_100yrRates(432, 231-50); // put sequence at 50-yr mark
		
		
		// these are long term MFDs so COV not needed
//		analysis.plotMFDs();

//		analysis.writeAAL_ValuesBillions();
		
		// this assumes COV = 0
//		analysis.plotLossRateVsMag();
		
		// Traditional loss exceedance curves:
//		analysis.makeFigAndTableSetA_Exceedances(1d, true);
		
//		analysis.makeFigSetF_N_ExceedHists();

		// Aggregate loss exceedance curves:
//		analysis.makeFigAndTableSetB_1yr_Exceedances();
		
		// This makes/plots TD and TI exceedance curves and ratios for different forecast durations.
//		analysis.makeFigSetD_Exceedances();
		
		// This makes/plots Aggregate TD and TI exceedance curves and ratios for different forecast durations.
//		analysis.makeFigSetE_Exceedances();
		
//		analysis.plotMaxAndAggrLossesPerYearInCatalog(432);


		// this plots the distribution of rupture losses as a function of magnitude (not including rup rates)
//		analysis.plotRupLossVsMagStats(analysis.getCatalogs(), false, true);
//		analysis.plotRupLossVsMagStats(analysis.getCatalogs(), true, true);

		
		// this isn't that informative (comment out part allows to toggle to 1.0-CDF); PDF is more informative
//		analysis.plotCumAggrYearlyLossDistributions(true, true);
		
//		System.exit(0);
		
		

		
//		// Map with greatest number of exceedances of loss=40e6 in subcat = 409649; for both 1000 and 2000 500-yr simulations
//		analysis.plotCatalogMap(409649, null, true);
//		analysis.plotCatalogMap(409649, null, false);
//		analysis.write1yrLargeEvents(409649);


		
		
		
		// OLD:
//		analysis.makeFig1_100yrRates(421, 113); // this start time is 50 yrs before first event
		
		// this doesn't look great (it save to file; no pop up window)
//		analysis.makeFigSetC_500yrRates(421);
		
//		// Map for highest aggr loss year for the 1000 500-yr simulations (first set)
//		Location loc1 = new Location(32.5d,-119.5d,0d);
//		Location loc2 = new Location(34.5d,-116.5d,0d);
//		Region region = new Region(loc1,loc2);
//		analysis.plotCatalogMap(210663, region, true);	
//		// this writes the large events in 1 yr TD subCatalog
//		analysis.write1yrLargeEvents(210663);

		
		
		// this doesn't look great
//		analysis.plotAccumulatedLossVsTimeForCatalog(421);
		
		/**
		 * This tests whether exceed and incremental distributions agree with respect to translations.
		 * The TI results do, but not the TD results, which I don't yet understand.  For TD, values are higher
		 * & closer to Poisson for incremental calculations (less discrepancy), which seems consistent with the
		 * overall rates being lower for incremental calculations (again, less discrepancy).  Why the conversion matches
		 * for Poisson and not for TD is still mysterious.
		 */
//		analysis.testIncrVsExceedCurveTanslations();
		//	This tests whether aggregate incremental distributions from exceedance (and vise versa) are the same.  
		//  They both are identical, but the conversions are different. 
//		analysis.testAggrIncrVsExceedCurveTanslations();
		
		/**
		 * This tests whether averaging non-exceed probs among catalogs (rather than exceed probs) gives a different result.
		 * The answer is no for both TI and TD.
		 */
//		analysis.testAveragingNonExceedProbs();
		
		

		

		

		



		
//		// TEST different time-ind curves - full PDF
//		// Any one of the following three work
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogs(analysis.getCatalogs());
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogs(analysis.getRandomizedCatalogs());
//		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogs(analysis.getSubcatalogList(analysis.getRandomizedCatalogs(), 1.0));
//		lossRateCurveFromCatalogs.scale(500d); // need this for the last option above
//		ArbitrarilyDiscretizedFunc lossExceedProbFromRateCatalogs = analysis.getBlankLossCurve(0.0);
//		lossExceedProbFromRateCatalogs.setName("lossExceedProbFromRateCatalogs");
//		for(int i=0;i<lossExceedProbFromRateCatalogs.size();i++) {
//			lossExceedProbFromRateCatalogs.set(i,1.0-Math.exp(-lossRateCurveFromCatalogs.getY(i)));
//		}
//		UncertainArbDiscFunc lossExceedProbFromCatalogsRandomizedCat = analysis.getLossExceedProbCurveFromCatalogs(analysis.getRandomizedCatalogs(), 1d);
//		lossExceedProbFromCatalogsRandomizedCat.setName("lossExceedProbFromCatalogsRandomizedCat");
//		ArrayList<XY_DataSet> funcProb3 = new ArrayList<XY_DataSet>(); 
//		funcProb3.add(lossExceedProbFromRateCatalogs);
//		funcProb3.add(lossExceedProbFromCatalogsRandomizedCat);
//		analysis.quickPlot(funcProb3, "Loss (thousand $)", "Prob (/yr)", "Test", true, true);

		
		// TEST different time-ind curves - random loss
//		// Any one of the following three work
//		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogsRandomFromCOV(analysis.getCatalogs());
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogsRandomFromCOV(analysis.getRandomizedCatalogs());
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogsRandomFromCOV(analysis.getSubcatalogList(analysis.getRandomizedCatalogs(), 1.0));
////		lossRateCurveFromCatalogs.scale(500d); // need this for the last option above
//		ArbitrarilyDiscretizedFunc lossExceedProbFromRateCatalogs = analysis.getBlankLossCurve(0.0);
//		lossExceedProbFromRateCatalogs.setName("lossExceedProbFromRateCatalogsRandomFromCOV");
//		for(int i=0;i<lossExceedProbFromRateCatalogs.size();i++) {
//			lossExceedProbFromRateCatalogs.set(i,1.0-Math.exp(-lossRateCurveFromCatalogs.getY(i)));
//		}
//		UncertainArbDiscFunc lossExceedProbFromCatalogsRandomizedCat = analysis.getLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getRandomizedCatalogs(), 1d);
//		lossExceedProbFromCatalogsRandomizedCat.setName("lossExceedProbFromCatalogsRandomizedCatRandomFromCOV");
//		ArrayList<XY_DataSet> funcProb3 = new ArrayList<XY_DataSet>(); 
//		funcProb3.add(lossExceedProbFromRateCatalogs);
//		funcProb3.add(lossExceedProbFromCatalogsRandomizedCat);
//		analysis.quickPlot(funcProb3, "Loss (thousand $)", "Prob (/yr)", "Test", true, true);

		
		
//		// trying to figure out what's wrong with the incremental distributions
//		UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV = analysis.getLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getCatalogs(), 1d);
//		lossExceedProbFromTD_CatalogsRandFromCOV.setName("lossExceedProbFromTD_CatalogsRandFromCOV");
//		UncertainArbDiscFunc lossExceedProbFromPoisCatalogsRandFromCOV = analysis.getLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getRandomizedCatalogs(), 1d);
//		lossExceedProbFromPoisCatalogsRandFromCOV.setName("lossExceedProbFromPoisCatalogsRandFromCOV");
//		
//		UncertainArbDiscFunc incrTestTD = analysis.getLossIncrProbDistFromCatalogsRandomFromCOV(analysis.getCatalogs(), 1d);
//		incrTestTD.setName("incrTestTD");
//		UncertainArbDiscFunc incrTestPoiss = analysis.getLossIncrProbDistFromCatalogsRandomFromCOV(analysis.getRandomizedCatalogs(), 1d);
//		incrTestPoiss.setName("incrTestPoiss");
//		ArrayList<UncertainArbDiscFunc> funcList = new ArrayList<UncertainArbDiscFunc>(); 
//		funcList.add(incrTestTD);
//		funcList.add(incrTestPoiss);
//		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
//		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
//		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
//		String fileNamePrefix = null;
//		PlottingUtils.writeAndOrPlotUncertFuncs(funcList, plotChars2, "Incr Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, fileNamePrefix, true);
//		
//		ArbitrarilyDiscretizedFunc exceedTestTD = getBlankLossCurve(0d);
//		ArbitrarilyDiscretizedFunc exceedTestPois = getBlankLossCurve(0d);
//		for(int i=0;i<exceedTestTD.size();i++) {
//			double nonExceed1=1, nonExceed2=1;
//			for(int j=i; j<exceedTestTD.size();j++) {
//				nonExceed1 *= (1.0-incrTestTD.getY(j));
//				nonExceed2 *= (1.0-incrTestPoiss.getY(j));
//			}
//			exceedTestTD.set(i,1.0-nonExceed1);
//			exceedTestPois.set(i,1.0-nonExceed2);
//		}
//		exceedTestTD.setName("exceedTestTD");
//		exceedTestPois.setName("exceedTestPois");
//		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>(); // AbstractDiscretizedFunc
//		funcList2.add(exceedTestTD);
//		funcList2.add(exceedTestPois);
//		funcList2.add(lossExceedProbFromTD_CatalogsRandFromCOV);
//		funcList2.add(lossExceedProbFromPoisCatalogsRandFromCOV);
//		
//		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars2, "Exceed Test", "Loss ($)", "Probability", new Range(1e6,1e12), new Range(1e-6,1.0), true, true, fileNamePrefix, true);

		


	
		
//		// Plotting
//		UncertainArbDiscFunc lossExceedProbFromTD_Catalogs = analysis.getLossExceedProbCurveFromCatalogs(analysis.getCatalogs(), 1d);
//		lossExceedProbFromTD_Catalogs.setName("lossExceedProbFromTD_Catalogs");
//		UncertainArbDiscFunc lossExceedProbFromTD_CatalogsRandFromCOV = analysis.getLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getCatalogs(), 1d);
//		lossExceedProbFromTD_CatalogsRandFromCOV.setName("lossExceedProbFromTD_CatalogsRandFromCOV");
//		UncertainArbDiscFunc aggrLossExceedProbFromTD_CatalogsRandFromCOV = analysis.getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getCatalogs(), 1d);
//		aggrLossExceedProbFromTD_CatalogsRandFromCOV.setName("aggrLossExceedProbFromTD_CatalogsRandFromCOV");
//		UncertainArbDiscFunc aggrLossExceedProbFromTI_CatalogsRandFromCOV = analysis.getAggrLossExceedProbCurveFromCatalogsRandomFromCOV(analysis.getRandomizedCatalogs(random), 1d);
//		aggrLossExceedProbFromTI_CatalogsRandFromCOV.setName("aggrLossExceedProbFromTI_CatalogsRandFromCOV");
//		ArrayList<XY_DataSet> testFuncList = new ArrayList<XY_DataSet>(); 
//		testFuncList.add(lossExceedProbFromTD_Catalogs);
//		testFuncList.add(lossExceedProbFromTD_Catalogs); // for conf bounds
//		testFuncList.add(lossExceedProbFromTD_CatalogsRandFromCOV);
//		testFuncList.add(lossExceedProbFromTD_CatalogsRandFromCOV); // for conf bounds
//		testFuncList.add(aggrLossExceedProbFromTD_CatalogsRandFromCOV);
//		testFuncList.add(aggrLossExceedProbFromTD_CatalogsRandFromCOV); // for conf bounds
//		testFuncList.add(aggrLossExceedProbFromTI_CatalogsRandFromCOV);
//		testFuncList.add(aggrLossExceedProbFromTI_CatalogsRandFromCOV); // for conf bounds
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.BLACK));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.ORANGE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.ORANGE));
//		PlottingUtils.writeAndOrPlotFuncs(testFuncList, plotChars, "test", "loss", "prob", null, null, false, true, null, true);



		
		
		// this demonstrates to be careful when log-plotting PDF x-axes (the distribution is distorted)
//		testTakingLogOfPDF_Xaxis();
//		System.exit(0);
				
		// This was a test showing that LossCOV_Model was not giving the correct distribution.  The two cases produced
		// below are now equivalent because the bug was fixed.
//		double mean = 1e5;
//		LossCOV_Model covModel = LossCOV_Model.PORTER_POWER_LAW_2020_09_01;
//		double cov = covModel.getCOV(mean);
//		double samples[] = 	covModel.getDistribution(mean).sample(1000000);
//		DescriptiveStatistics stats = new DescriptiveStatistics(samples);
//		double meanFromDist = stats.getMean();
//		double covFromDist = stats.getStandardDeviation()/meanFromDist;
//		System.out.println("mean="+mean+";\tmeanFromDist="+meanFromDist+";\tcov="+cov+";\tcovFromDist="+covFromDist+
//				"; meanRatio="+(float)(meanFromDist/mean)+"; covRatio="+(float)(covFromDist/cov));
//		double sigma = Math.sqrt(Math.log(cov*cov+1));
//		double mu = Math.log(mean)-(sigma*sigma/2);
//		LogNormalDistribution dist = new LogNormalDistribution(mu, sigma);
//		double samples2[] = dist.sample(1000000);
//		DescriptiveStatistics stats2 = new DescriptiveStatistics(samples2);
//		double meanFromDist2 = stats2.getMean();
//		double covFromDist2 = stats2.getStandardDeviation()/meanFromDist2;
//		System.out.println("mean="+mean+";\tmeanFromDist2="+meanFromDist2+";\tcov="+cov+";\tcovFromDist2="+covFromDist2+
//				"; meanRatio2="+(float)(meanFromDist2/mean)+"; covRatio2="+(float)(covFromDist2/cov));
//		System.exit(0);

		
		

		
//		// This tests the covModel fix in Aug 2024; blue and black should be identical
//		ArrayList<XY_DataSet> testFuncArray = analysis.testCOV_ModelChangesInAug2024();
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//		PlottingUtils.writeAndOrPlotFuncs(testFuncArray, plotChars, "test covModel changes", "loss", "prob", null, null, false, true, null, true);


		
//		System.exit(0);
//		
//
//		

//		// plot loss PDF in log10 space & compare to pure gaussian (this ignores spike at zero loss)
		// note that this filters out the gazillion little earthquakes that have zero loss (so median and mode are really zero)
//		ArbitrarilyDiscretizedFunc testCurve = analysis.getLossExceedProbCurveFromERF();
//		HistogramFunction testHist = convertExceedCurveToLog10_PDF_Hist(testCurve);
//		testHist.setName("testHist");
////		testHist.set(0,0.0);
//		double mean = testHist.computeMean();
//		double stdDev = testHist.computeStdDev();
//		testHist.setInfo("Mean = "+mean+"\nStdDev = "+stdDev);
//		org.opensha.commons.calc.GaussianDistCalc gaussDistCalc = new org.opensha.commons.calc.GaussianDistCalc();
//		EvenlyDiscretizedFunc testHistFit = testHist.deepClone();
//		for(int i=1;i<testHistFit.size();i++) {
//			double normVal = (testHistFit.getX(i)-mean)/stdDev;
//			double yVal = Math.exp(-0.5*normVal*normVal);
//			testHistFit.set(i,yVal);
//		}
//		testHistFit.scale(testHist.getMaxY()/testHistFit.getMaxY());
//		ArrayList<XY_DataSet> testFuncList2 = new ArrayList<XY_DataSet>(); 
//		testFuncList2.add(testHist);
//		testFuncList2.add(testHistFit);
//		analysis.quickPlot(testFuncList2, "Loss (thousand $)", "Prob", "PDF", false, false);
		

		
		
		
		
//		// plot long-term loss rate curves
//		ArbitrarilyDiscretizedFunc lossRateCurveFromERF = analysis.getLossRateCurveFromERF();
//		lossRateCurveFromERF.setName("lossRateCurveFromERF");
//		ArbitrarilyDiscretizedFunc lossRateCurveFromERF_ZeroCOV = analysis.getLossRateCurveFromERF_zeroCOV();
//		lossRateCurveFromERF_ZeroCOV.setName("lossRateCurveFromERF_ZeroCOV");
//		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogs(analysis.getCatalogs());
//		lossRateCurveFromCatalogs.setName("lossRateCurveFromCatalogs");
//		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>(); 
//		funcs.add(lossRateCurveFromERF);
//		funcs.add(lossRateCurveFromERF_ZeroCOV);
//		funcs.add(lossRateCurveFromCatalogs);
//		analysis.quickPlot(funcs, "Loss (thousand $)", "Rate (/yr)", "Loss Rate from ERF", true, true);
		
		
//		// plot 1-yr loss exceed from ERF
//		ArbitrarilyDiscretizedFunc lossExceedProbFromERF = analysis.getLossExceedProbCurveFromERF();
//		lossExceedProbFromERF.setName("lossExceedProbFromERF");
//		ArbitrarilyDiscretizedFunc testFunc = analysis.getBlankLossCurve(0.0);
//		ArbitrarilyDiscretizedFunc lossExceedProbFromERF_ZeroCOV = analysis.getBlankLossCurve(0.0);
//		ArbitrarilyDiscretizedFunc lossExceedProbFromCatalogs = analysis.getBlankLossCurve(0.0);
//
//		for(int i=0;i<testFunc.size();i++) {
//			testFunc.set(i,1.0-Math.exp(-lossRateCurveFromERF.getY(i)));
//			lossExceedProbFromERF_ZeroCOV.set(i,1.0-Math.exp(-lossRateCurveFromERF_ZeroCOV.getY(i)) );
//			lossExceedProbFromCatalogs.set(i,1.0-Math.exp(-lossRateCurveFromCatalogs.getY(i)) );
//		}
//		testFunc.setName("testFunc (should be same as lossExceedProbFromERF)");
//		lossExceedProbFromERF_ZeroCOV.setName("lossExceedProbFromERF_ZeroCOV");
//		lossExceedProbFromCatalogs.setName("lossExceedProbFromCatalogs");
//		ArrayList<XY_DataSet> funcProb = new ArrayList<XY_DataSet>(); 
//		funcProb.add(lossExceedProbFromERF);
//		funcProb.add(lossExceedProbFromERF_ZeroCOV);
//		funcProb.add(lossExceedProbFromCatalogs);
//		funcProb.add(testFunc);
//		analysis.quickPlot(funcProb, "Loss (thousand $)", "Prob (/yr)", "Loss Exceed Curves", true, true);
//
//
//		UncertainArbDiscFunc lossExceedCurveFromRandCats = analysis.getLossExceedProbCurveFromCatalogs(analysis.getRandomizedCatalogs(random), 1d);
//		lossExceedCurveFromRandCats.setName("lossExceedCurveFromRandCats");
//		ArrayList<XY_DataSet> funcProb2 = new ArrayList<XY_DataSet>(); 
//		funcProb2.add(lossExceedProbFromCatalogs);
//		funcProb2.add(lossExceedCurveFromRandCats);
//		analysis.quickPlot(funcProb2, "Loss (thousand $)", "Prob (/yr)", "Test", true, true);
		
	
		
		
		
		
//		analysis.old_writeAndOrPlotRateVersusTime(analysis.getCatalogs(), null, true, "Orig Catalogs");
//		analysis.old_writeAndOrPlotRateVersusTime(analysis.getRandomizedCatalogs(null), null, true, "Randomized Catalogs");
//		
//		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogs = analysis.getLossRateCurveFromCatalogs(analysis.getCatalogs());
//		ArbitrarilyDiscretizedFunc lossExceedProbFromRateCatalogs = analysis.getBlankLossCurve(0.0);
//		lossExceedProbFromRateCatalogs.setName("lossExceedProbFromRateCatalogs");
//		for(int i=0;i<lossExceedProbFromRateCatalogs.size();i++) {
//			lossExceedProbFromRateCatalogs.set(i,1.0-Math.exp(-lossRateCurveFromCatalogs.getY(i)));
//		}
//		UncertainArbDiscFunc lossExceedProbFromRandCatalogs = analysis.getLossExceedProbCurveFromCatalogs(analysis.getRandomizedCatalogs(null), 1d);
//		lossExceedProbFromRandCatalogs.setName("lossExceedProbFromRandCatalogs");
//		UncertainArbDiscFunc lossExceedProbFromTD_Catalogs = analysis.getLossExceedProbCurveFromCatalogs(analysis.getCatalogs(), 1d);
//		lossExceedProbFromTD_Catalogs.setName("lossExceedProbFromTD_Catalogs");
//		UncertainArbDiscFunc lossExceedProbFromGKdeclusteredCatalogs = analysis.getLossExceedProbCurveFromCatalogs(analysis.getCatalogsDeclustered(), 1d);
//		lossExceedProbFromGKdeclusteredCatalogs.setName("lossExceedProbFromGKdeclusteredCatalogs");
//		ArrayList<XY_DataSet> funcProb3 = new ArrayList<XY_DataSet>(); 
//		funcProb3.add(lossExceedProbFromRateCatalogs);
//		funcProb3.add(lossExceedProbFromRandCatalogs);
//		funcProb3.add(lossExceedProbFromCatalogs);
//		funcProb3.add(lossExceedProbFromTD_Catalogs);
//		funcProb3.add(lossExceedProbFromGKdeclusteredCatalogs);
//		analysis.quickPlot(funcProb3, "Loss (thousand $)", "Prob (/yr)", "Test", true, true);

		
		
		
		
//		// VERIFYING EQUIVALENCE OF VARIOUS TIME-IND APPROACHES FOR ZERO COV CASE (DURING DEBUGGIN)
//
//		int numEvents = 0;
//		for(ObsEqkRupList list:analysis.getCatalogs())
//			numEvents+=list.size();
//		System.out.println("Catalogs numEvents: "+numEvents);
//
//		numEvents = 0;
//		ArrayList<ObsEqkRupList> randCatalogs = analysis.getRandomizedCatalogs(null);
//		for(ObsEqkRupList list:randCatalogs)
//			numEvents+=list.size();
//		System.out.println("randCatalogs numEvents: "+numEvents);
//		ArrayList<ObsEqkRupList> randSubCatalogs = analysis.getSubcatalogList(randCatalogs, 1.0);
//		numEvents = 0;
//		for(ObsEqkRupList list:randSubCatalogs)
//			numEvents+=list.size();
//		System.out.println("randSubCatalogs numEvents: "+numEvents);
//
//		
//		// Any one of the following three work
//		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogsZeroCOV = analysis.getLossRateCurveFromCatalogsZeroCOV(analysis.getCatalogs());
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogsZeroCOV = analysis.getLossRateCurveFromCatalogsZeroCOV(randCatalogs);
////		ArbitrarilyDiscretizedFunc lossRateCurveFromCatalogsZeroCOV = analysis.getLossRateCurveFromCatalogsZeroCOV(randSubCatalogs);
////		lossRateCurveFromCatalogsZeroCOV.scale(500d);
//		ArbitrarilyDiscretizedFunc lossExceedProbFromRateCatalogsZeroCOV = analysis.getBlankLossCurve(0.0);
//		lossExceedProbFromRateCatalogsZeroCOV.setName("lossExceedProbFromRateCatalogsZeroCOV");
//		for(int i=0;i<lossExceedProbFromRateCatalogsZeroCOV.size();i++) {
//			lossExceedProbFromRateCatalogsZeroCOV.set(i,1.0-Math.exp(-lossRateCurveFromCatalogsZeroCOV.getY(i)));
//		}
//		ArbitrarilyDiscretizedFunc lossExceedProbFromCatalogsZeroCOV_RandCat = analysis.getLossExceedProbCurveFromCatalogsZeroCOV(randCatalogs, 1d);
//		lossExceedProbFromCatalogsZeroCOV_RandCat.setName("lossExceedProbFromCatalogsZeroCOV_RandCat");
//		ArbitrarilyDiscretizedFunc lossExceedProbFromCatalogsZeroCOV = analysis.getLossExceedProbCurveFromCatalogsZeroCOV(analysis.getCatalogs(), 1d);
//		lossExceedProbFromCatalogsZeroCOV.setName("lossExceedProbFromCatalogsZeroCOV");
//		ArrayList<XY_DataSet> funcProb3 = new ArrayList<XY_DataSet>(); 
//		funcProb3.add(lossExceedProbFromRateCatalogsZeroCOV);
//		funcProb3.add(lossExceedProbFromCatalogsZeroCOV_RandCat);
////		funcProb3.add(lossExceedProbFromCatalogsZeroCOV);
//		analysis.quickPlot(funcProb3, "Loss (thousand $)", "Prob (/yr)", "Test", true, true);

		
		
		// REALLY OLD STUFF BELOW
		
//		Random random = new Random(102553864); // for reproducibility; change argument to get different results
//		analysis.getRandomizedCatalogs(random);  // make and save the randomized catalogs

		
		// STUFF FROM U3ETAS_SimulationAnalysis
		
		// FIGURE 2 (MFDs)
//		ArrayList<XY_DataSet> funcList = analysis.makeMFDs(false);
//		MakeFigures.makeFigure2_Parts(funcList);

		
		
//		// FIGURE 3 (Hazard Curves)
//		filePrefix = "Figure3";
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		// add random IML for Full TD
//		funcsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, true, imr)[1]);
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure3";
//		MakeFigures.makeFigure3_Parts(dataSetsArray, funcsArray, duration, saPeriod, dirName, true, filePrefix);

		
		
		// FIGURE 9 (1-yr Hazard Curves)
//		duration = 1;
//		filePrefix = "Figure9";
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure9";
//		MakeFigures.makeFigure9_Parts(dataSetsArray, funcsArray, duration, saPeriod, dirName, true, filePrefix);

		
		// FIGURES 4 and 5
//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS USING RANDOM IML SAMPLES
//		double[] imlArray = new double[]{0.423, 1.67}; // the first value is where the probability difference is maximum
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2);
//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, null);

		// THIS IS TO LOOK AT HAZARD CURVES WHERE THE FULL TD AND RANDOMIZED CASES ARE MOST DIFFERENT (SE CORNER OF REGION)
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		loc = new Location(33.9,-113.5);
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/???????";
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, "Lat=33.9, Lon=-113.5");
//		double[] imlArray = new double[]{0.10,0.17};
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2, "maxDiffNumHazExceedanceRandIML_Test");
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, "maxDiffNumHazExceedanceTest");

		
		// THIS IS TO LOOK AT HAZARD CURVES WHERE THE FULL TD AND RANDOMIZED CASES ARE 2nd MOST DIFFERENT (NEAR BIG LAGOON - BALD MTN FAULT)
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		loc = new Location(42.0,-124.6);
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, false, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));
//		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/???????";
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, "Lat=33.9, Lon=-113.5");
//		double[] imlArray = new double[]{0.65,0.80};
//		ArrayList<ArrayList<XY_DataSet>> listList1 = analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats & Rand IML", random);
//		ArrayList<ArrayList<XY_DataSet>> listList2 = analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats & Rand IML", random);
//		MakeFigures.makeNumHazExceedanceRandIML_Figures(imlArray, listList1, listList2, "maxDiffNumHazExceedanceRandIML_Test_BigLagoon");
//		ArrayList<ArrayList<XY_DataSet>> listList3 = analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, imlArray, imr, "Full TD Cats");
//		ArrayList<ArrayList<XY_DataSet>> listList4 = analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(random), loc, duration, saPeriod, imlArray, imr, "Full Randomized Cats");
//		MakeFigures.makeNumHazExceedanceFigures(imlArray, listList3, listList4, "maxDiffNumHazExceedanceTest_BigLagoon");

		
		// REMOVE THIS METHOD ????????????????????? (RUNS TAKE ~90 MINUTES)
//		makeERF_GKfilter_HazMapData(0.2, false);
//		System.exit(0);
		
		// FOR KEVIN **********************
		
		// loop over sites in RELM gridded region
//		Location loc = new Location(34.05,-118.25);
//		Location loc = new Location(34.1,-118.2); // test loc to confirm gridded region curve
		// loop over saPeriods of 0, 0.2, 1, 5
//		double saPeriod = 0; // 0 = PGA
//		
//		ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
//		
//		double duration = 50;
		
//		Random random = new Random();
//		
//		ArrayList<ObsEqkRupList> catalogsList = U3ETAS_SimulationAnalysis.loadCatalogs(new File(fssFileName), new File(catalogsFileName));
//		UncertainArbDiscDataset[] datasetsArray1 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(catalogsList, loc, duration, saPeriod, false, imr);
//		ArbitrarilyDiscretizedFunc poissonCurve = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogsPoisson(catalogsList, loc, duration, saPeriod, imr);
//		ArrayList<ObsEqkRupList> declusteredCatalogsList = U3ETAS_SimulationAnalysis.getGK_DeclusteredCatalog(catalogsList);
//		UncertainArbDiscDataset[] datasetsArray2 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(declusteredCatalogsList, loc, duration, saPeriod, false, imr);
//		ArrayList<ObsEqkRupList> catalogsRandmizedList = U3ETAS_SimulationAnalysis.getRandomizedCatalogs(catalogsList, random);
//		UncertainArbDiscDataset[] datasetsArray3 = U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(catalogsRandmizedList, loc, duration, saPeriod, false, imr);
		
		// write out contents of poissonCurve, datasetsArray1, datasetsArray2, and datasetsArray3,
		
		// END FOR KEVIN **********************
		
		
		
		
//		// TEST PLOT THE ABOVE:
//		ArrayList<UncertainArbDiscDataset[]> dataSetsArray = new ArrayList<UncertainArbDiscDataset[]>();
//		ArrayList<XY_DataSet> funcsArray = new ArrayList<XY_DataSet>();
//		
//		String locName = "Los Angeles";
		
//		dataSetsArray.add(datasetsArray1);
//		dataSetsArray.add(datasetsArray2);
//		dataSetsArray.add(datasetsArray3);
//		funcsArray.add(poissonCurve);
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

	
		
//		U3ETAS_SimulationAnalysis analysis = new U3ETAS_SimulationAnalysis();
//		
//		// TEST OF OUR PRIMARY COMPARISON
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr)); // Magenta
//		funcsArray.add(computeHazardCurvesFromCatalogsPoisson(getU3_GK_FilteredCatalog(analysis.getCatalogs()), loc, duration, saPeriod, imr)); // Magenta
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, false, imr)); // Blue
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getCatalogsDeclustered(), loc, duration, saPeriod, false, imr));	// Red

//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, false, imr)); // Black
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

		
//		analysis.doMFDs();
		
//		// COMPARE ERF-TI CURVES WITH AND WITHOUT GK FILTER
//		FaultSystemSolutionERF erf = getTimeIndERF_Instance(duration, true);
//		ArbitrarilyDiscretizedFunc curveWith = computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr);
//		curveWith.setName("ERF-TI with With GK Filter");
//		funcsArray.add(curveWith); 
//		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
//		erf.updateForecast();
//		ArbitrarilyDiscretizedFunc curveWithOut = computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr);
//		curveWithOut.setName("ERF-TI with Without GK Filter");
//		funcsArray.add(curveWithOut); 
//		writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS COMPARES HAZARD CURVES USING THE TI ERF WITH RANDOM CATALOGS FROM THE ERF
//		// THE MATCH IS PERFECT
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); 
//		System.out.println("starting Rand Cats");
//		ArrayList<ObsEqkRupList> randCats = analysis.makeRandomCatalogListFromERF(erf);
//		System.out.println("done with Rand Cats");
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(randCats, loc, duration, saPeriod, false, imr));
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(randCats, loc, duration, saPeriod, imr));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS COMPARES THE USGS LA HAZARD CURVE TO THAT COMPUTED WITH THE ERF FOR 4 DIFF IMRs
//		// MATCH IS GOOD AT 2 IN 50 ; 10% OFF AT LOWER AND HIGHES END
//		// NO ADDED GMM EPISTEMIC UNCERTAINTY AND ONE SMOOTHED SIES BRANCH HERE SO SOME DIFFS EXPECTED
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, true);
//		erf.updateForecast();
//		imr = AttenRelRef.CB_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.ASK_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.BSSA_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		imr = AttenRelRef.CY_2014.instance(null);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr)); // HAD TO SKIP M≤5 EVENTS IN THE HAZARD CURVE CALCULATOR (HARD CODED) FOR THIS TO WORK
//		funcsArray.add(analysis.getUSGS_LosAngelesPGA_2in50_HazardCurve());
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);

		
//		// THIS COMPARES ERF HAZARD CURVE TO U3ETAS RANDOMIZED CURVE; DIFFERENCE IS 16% AT 2 IN 50, WHICH SEEMS HIGH
//		// BUT IS WITHIN UNCERTAINTIES AT LA FOR U3-TI SUPPLEMENTAL DATA
//		// MFD COMPARISON SHOWS THE PROBLEM THAT U3ETAS IS UNDERPREDICTING (by 26%) EVENTS AROUND M≥7, WHICH IS ABOUT THE SAME
//		// AMOUNT THE BSSA PAPER SHOWS.  NEED TO BUMP UP THE totRateScaleFactor BY 31% (1.14*1.31=1.5) TO 
//		// BALANCE MFD DISCREPENCIES. NoERT CASE SHOWS LESS, BUT STILL A DISCREPANCY (14% ON MFD; 12% ON 2IN50).
//		FaultSystemSolutionERF erf = analysis.getTimeIndERF_Instance(duration);
//		funcsArray.add(analysis.computeHazardCurveFromERF(erf, loc, saPeriod, duration, imr));
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		funcsArray.add(analysis.getUSGS_LosAngelesPGA_2in50_HazardCurve());
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
//		// COMPARE MFDS
//		ArrayList<XY_DataSet> funcsArray2 = new ArrayList<XY_DataSet>();
//		SummedMagFreqDist erfMFD = ERF_Calculator.getTotalMFD_ForERF(erf, minMag, 8.95, numMag, true);
//		funcsArray2.add(erfMFD);
//		funcsArray2.add(erfMFD.getCumRateDistWithOffset());
//		IncrementalMagFreqDist catMFD = analysis.makeMFD(analysis.getCatalogs());
//		funcsArray2.add(catMFD);
//		funcsArray2.add(catMFD.getCumRateDistWithOffset());
//		double aveLnRatio = 0;
//		int num=0;
//		for(int i=0;i<erfMFD.size();i++)
//			if(catMFD.getY(i) > 0.0) {
//				aveLnRatio += Math.log((erfMFD.getY(i)/catMFD.getY(i)));
//				num+=1;
//			}
//		aveLnRatio /=num;
//		System.out.println("aveRatio="+Math.exp(aveLnRatio));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray2, duration, saPeriod, null, true, locName);


//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS USING RANDOM IML SAMPLES
//		analysis.computeNumHazardExceedDistRandomIML(analysis.getCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);
//		analysis.computeNumHazardExceedDistRandomIML(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);

//		// THIS LOOKS AT NUM EXCEEDANCE DISTRIBUTIONS AT THE SPECIFIED IMLS
//		analysis.computeNumHazardExceedDist(analysis.getCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);
//		analysis.computeNumHazardExceedDist(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, new double[]{0.02,0.2,0.7,2.0}, imr);

		
		
		
//		// THIS LOOKS AT THE NUM DISTRIBUTION OF EVENTS ABOVE VARIOUS MAG THRESHOLDS
//		// THINGS LOOK POISSONIAN FOR M≥7 ISH
//		duration = 50;
//		analysis.computeNumEventDistribution(analysis.getCatalogs(), duration, new double[]{5.0,6.0, 7.0, 8.0});
//		analysis.computeNumEventDistribution(analysis.getCatalogsDeclustered(), duration, new double[]{5.0,6.0, 7.0, 8.0});
//		analysis.computeNumEventDistribution(analysis.getRandomizedCatalogs(), duration, new double[]{5.0,6.0, 7.0, 8.0});
		
//		// THIS PLOTS CUMULATIVE RATES VERSUS TIME
//		analysis.writeAndOrPlotRateVersusTime(analysis.getCatalogs(), null, true, "Test U3ETAS");
//		analysis.writeAndOrPlotRateVersusTime(analysis.getRandomizedCatalogs(), null, true, "Test Poisson");
		
//		// THIS TESTS THAT CUTTING SUMULATIONS IN HALF PRODUCES CONSISENT RESULTS
//		// cut catalogs in half:
//		ArrayList<ObsEqkRupList> allCatalogs = analysis.getCatalogs();
//		ArrayList<ObsEqkRupList> firstHalfCats = new ArrayList<ObsEqkRupList> ();
//		ArrayList<ObsEqkRupList> secondHalfCats = new ArrayList<ObsEqkRupList> ();
//		for(int i=0;i<500;i++) {
//			firstHalfCats.add(allCatalogs.get(i));
//			secondHalfCats.add(allCatalogs.get(i+500));
//		}
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(firstHalfCats, loc, duration, saPeriod,false, imr));
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(secondHalfCats, loc, duration, saPeriod,false, imr));
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS TESTS THAT A RANDOMIZED CATALOG PRODUCES THE SAME CURVE AS THAT FROM A POISSON MODEL USING THE AVERAGE NUMBER OF EXCEEDANCES
//		funcsArray.add(analysis.computeHazardCurvesFromCatalogsPoisson(analysis.getCatalogs(), loc, duration, saPeriod, imr));
//		dataSetsArray.add(U3ETAS_SimulationAnalysis.computeHazardCurvesFromCatalogs(analysis.getRandomizedCatalogs(), loc, duration, saPeriod, false, imr));	
//		U3ETAS_SimulationAnalysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
//		// THIS TESTS THAT RANDOM IML SAMPLES PRODUCE THE SAME RESULT
//		boolean randIML = false;
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, randIML, imr));
//		randIML = true;
//		dataSetsArray.add(analysis.computeHazardCurvesFromCatalogs(analysis.getCatalogs(), loc, duration, saPeriod, randIML, imr));
//		analysis.writeAndOrPlotHazardCurves(dataSetsArray, funcsArray, duration, saPeriod, null, true, locName);
		
		// THIS TESTS THAT RANDOM IMLs FROM THE IMR REPRODUCE THE DISTRIBUTION, INCLUDING TRUNCATION
//		analysis.testRandomSamplesFromIMR(analysis.getCatalogs().get(0).get(0), loc, saPeriod, imr);


	}

}
