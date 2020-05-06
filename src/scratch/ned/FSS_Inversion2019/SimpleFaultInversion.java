package scratch.ned.FSS_Inversion2019;


import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;

import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.cybershake.maps.GMT_InterpolationSettings;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;
import org.opensha.sha.earthquake.param.FaultGridSpacingParam;


import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import org.apache.commons.math3.stat.StatUtils;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.simulatedAnnealing.completion.CompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.CompoundCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.EnergyCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.IterationCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.params.CoolingScheduleType;
import scratch.UCERF3.simulatedAnnealing.params.GenerationFunctionType;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.ScalingRelationshipEnum;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.SlipAlongRuptureModelEnum;

/**
 * 
 * This class runs the inversion.
 * 
 * Slip rate and event rate standard deviations are computed as the difference of the 95% confidence bounds divided by 4.
 * 
 * @author 
 *
 */
public class SimpleFaultInversion {

	final static boolean D = true;	// debugging flag
	
	public final static String ROOT_PATH = "src/scratch/ned/FSS_Inversion2019/Results/";
	final static String ROOT_DATA_DIR = "src/scratch/ned/FSS_Inversion2019/data/"; // where to find the data; this really needed?

	
	// These values are the same for all fault sections
	final static double UPPER_SEIS_DEPTH = 0;
	final static double LOWER_SEIS_DEPTH = 14.;
	final static double FAULT_DIP = 90; // FIX MAX_SUBSECT_LENGTH_KM calc below if this is changed
	final static double FAULT_RAKE = 0.0;
	final static double FAULT_SLIP_RATE = 30;	// mm/yr

	double FAULT_LENGTH_KM = 14*29;
	int NUM_SUBSECT_PER_RUP = 4;
	
	final static double hazGridSpacing = 0.05;
		
	ArrayList<FaultSectionPrefData> faultSectionDataList;
	int[][] rupSectionMatrix;
	
	int numSections, numRuptures;
		
	FaultSectionPrefData parentFaultData;
	
	ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
	
	public ArrayList<XY_DataSet> savedHazardCurves = new ArrayList<XY_DataSet>();



	public enum InversionSolutionType {
		FRESH,
		SUNFISH,
		RATES_FROM_MFD,
		NON_NEGATIVE_LEAST_SQUARES,
		SIMULATED_ANNEALING,
		FROM_FILE,
		APPLY_SOLUTION;
	}
	double[] solutionRatesToApplyArray;	

	public enum MFD_TargetType {
		GR_b_0pt32,	// equiv to NSHM combined for EllsworthB
		GR_b_1pt6,
		GR_b_1pt4,
		GR_b_1pt2,
		GR_b_1pt0,
		GR_b_0pt8,
		GR_b_0pt6,
		GR_b_0pt4,
		GR_b_0pt2,
		GR_b_0pt0,
		GR_b_minus1,
		MAX_RATE,  // all in minimum magnitude bin
		MIN_RATE,  // all in maximum magnitude bin
		M7pt25only,
		NSHM_Combined,
		GR_b_0pt8_ForSegmented,
		NONE;
	}

	public enum SlipRateProfileType {
		UNIFORM,
		TAPERED,
		UNI_TRIM;
	}
	
	public enum SA_InitialStateType {
		ALL_ZEROS,
		A_PRIORI_RATES,
		FROM_MFD_CONSTRAINT;
	}

	// the following two should be changed together; WARNING - code below assumes these only have two elements!
	double[] hazardProbArray = {0.02, 0.10};
	String[] hazardProbNameArray = {"2in50", "10in50"};
	double hazardDurationYrs = 50;

	double hazCurveLnMin = Math.log(0.001);
	double hazCurveLnMax = Math.log(10);
	int hazCurveNum = 40;
	double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
	
	// Adjustable Parameters
	String dirName;
	String solutionName; 						// Inversion name
	public boolean shortFault=true;							// whther fault is ~800 or 1600 km
	boolean wtedInversion; 						// Whether to apply data weights
	SlipRateProfileType slipRateProfile;
	SlipAlongRuptureModelEnum slipModelType;	// Average slip along rupture model
	ScalingRelationshipEnum scalingRel; 		// Scaling Relationship
//	boolean applySlipRateSegmentation;			// whether or not to reduce slip rate/std at middle section (hard coded)
	ArrayList<SlipRateSegmentationConstraint> slipRateSegmentationConstraintList;
	ArrayList<SectionRateConstraint> sectionRateConstraintList;
	double relativeSectRateWt; 					// Section Rate Constraints wt
	ArrayList<SegmentationConstraint> segmentationConstrList; 			// 
	double relative_segmentationConstrWt; 		// Segmentation Constraints wt
	double relative_aPrioriRupWt;  				// A priori rupture rates wts
	String aPrioriRupRateFilename;  			// A priori rupture rates file name; e.g., ROOT_DATA_DIR+"aPrioriRupRates.txt"
	boolean setAprioriRatesFromMFD_Constraint;	// this overides the above file option if both are set
	double minRupRate; 							// Minimum Rupture Rate; not applied for NSHMP solution
	boolean applyProbVisible; 					// Apply prob visible
	double moRateReduction;						// reduction due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)		
	MFD_TargetType mfdTargetType; 				// Target MFD Constraint
	double relativeMFD_constraintWt; 			// Target MFD Constraint Wt
	double totalRateConstraint; 				// Total Rate Constraint
	double totalRateSigma; 						// Total Rate Constraint
	double relativeTotalRateConstraintWt;
	ArrayList<int[]> smoothnessConstraintList;	// Section rate smoothness constraints
	double relativeSmoothnessConstraintWt;		// Section rate smoothness weight

	
	InversionSolutionType solutionType;
	long randomSeed;							// for SA reproducibility 
	int numSolutions; 							// this is ignored for NON_NEGATIVE_LEAST_SQUARES which only has one possible solution
	CoolingScheduleType saCooling;
	GenerationFunctionType perturbationFunc;
	CompletionCriteria completionCriteria;
	SA_InitialStateType initialStateType;
	boolean applyRuptureSampler;
//	String rupRatesFromFileDirName;				//if solution is from a file (previous solution); name example: ROOT_PATH+"OutDir/";
	double magAareaAleatoryVariability;			// the gaussian sigma for the mag-area relationship
	// Data to plot and/or save:
	boolean popUpPlots;							// this tells whether to show plots in a window (set null to turn off; e.g., for HPC)
	LocationList hazCurveLocList; 						// to make hazard curve (set loc=null to ignore)
    ArrayList<String> hazCurveLocNameList;
    boolean makeHazardMaps, makeHazardMapsForAllSolutions=false; 					// to make hazard maps
//	double saPeriodForHaz;						// set as 0.0 for PGA


	
	public SimpleFaultInversion() {
		setDefaultParameterValuess();
	}
		
	/**
	 * This read rupture rates from a file
	 * @param fileName - the full path
	 * @return
	 */
	private double[] readRuptureRatesFromFile(String fileName) {
		double rates[] = null;
		File file = new File(fileName);
		List<String> fileLines;
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			rates = new double[fileLines.size()];
			for(int i=0; i<fileLines.size();i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split("\t");
				rates[i] = Double.parseDouble(split[1]);
//				System.out.println(i+"\t"+split[0]+"\t"+rates[i]);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rates;
	}
	
	/**
	 * 
	 */
	private GriddedGeoDataSet readHazardMapDataFromFile(String fileName) {
		GriddedGeoDataSet griddedDataSet = new GriddedGeoDataSet(getGriddedRegion(), true);
		try {
			File file = new File(fileName);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			Preconditions.checkState(fileLines.size()-1 == griddedDataSet.size(), "File and gridded region have different numbers of points (%s vs %s, respectively)", fileLines.size(), griddedDataSet.size());
			for(int i=0;i<griddedDataSet.size();i++) {
				String str = fileLines.get(i+1);	 // skip header
				String[] split = str.split("\t");
				int index = Integer.parseInt(split[0]);
				Preconditions.checkState(index == i, "Incompatible indices (%s vs %s, respectively)", index, i);
				griddedDataSet.set(i, Double.parseDouble(split[1]));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return griddedDataSet;	
	}
	
	
	public void makeFaultSectionDataList() {
		
		double smallCorr = 1e-4; // this is to ensure that there is not an extra subsection (and smaller subsection length)
		if(shortFault) {
			FAULT_LENGTH_KM = 14*29-smallCorr;
			NUM_SUBSECT_PER_RUP = 4;			
		}
		else {
			FAULT_LENGTH_KM = 14*29*2-smallCorr;
			NUM_SUBSECT_PER_RUP = 2;			
		}

		FaultTrace parentFaultTrace = new FaultTrace("parentFaultTrace");
		parentFaultTrace.add(new Location(33d,-117d));
		parentFaultTrace.add(new Location(33d+FAULT_LENGTH_KM/111.195052,-117));
		System.out.println("Fault trace length: "+parentFaultTrace.getTraceLength());
		parentFaultData = new FaultSectionPrefData();
		parentFaultData.setFaultTrace(parentFaultTrace);
		parentFaultData.setAveDip(FAULT_DIP);
		parentFaultData.setAveRake(FAULT_RAKE);
		parentFaultData.setAveUpperDepth(UPPER_SEIS_DEPTH);
		parentFaultData.setAveLowerDepth(LOWER_SEIS_DEPTH);
		parentFaultData.setAveSlipRate(FAULT_SLIP_RATE);
		parentFaultData.setSlipRateStdDev(FAULT_SLIP_RATE/10d);
		if(FAULT_DIP != 90)
			throw new RuntimeException("Fix following code to handle dip != 90");
		double MAX_SUBSECT_LENGTH_KM = (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)/NUM_SUBSECT_PER_RUP;

		faultSectionDataList = parentFaultData.getSubSectionsList(MAX_SUBSECT_LENGTH_KM);
		
		// set all the section ids so we can apply constraints
		for(int i=0; i<faultSectionDataList.size() ;i++) {
			faultSectionDataList.get(i).setSectionId(i);
//			System.out.println(faultSectionDataList.get(i).getFaultTrace().getTraceLength());
		}
//		System.out.println(parentFaultTrace.getTraceLength());
//		System.exit(0);
		
		double moRateSum = 0, lengthSum = 0;
		for(FaultSectionPrefData data : faultSectionDataList) {
			lengthSum += data.getTraceLength();
			moRateSum += data.calcMomentRate(true);
		}
		
		if(D) {
			double totMoRateFromSections = faultSectionDataList.get(0).calcMomentRate(true)*faultSectionDataList.size();
			double totLengthFromSections = faultSectionDataList.get(0).getTraceLength()*faultSectionDataList.size();
			System.out.println("MoRate Tests (parent, subSec sums): "+parentFaultData.calcMomentRate(true)+"\t"+totMoRateFromSections+"\t"+moRateSum);
			System.out.println("Length Tests (parent, subSec sums): "+parentFaultData.getTraceLength()+"\t"+totLengthFromSections+"\t"+lengthSum);
			System.out.println("DDW Tests (parent, subSec): "+parentFaultData.getReducedDownDipWidth()+"\t"+faultSectionDataList.get(0).getReducedDownDipWidth());			
		}
	}
	
		
	private void initData() {

		
		// need to re-make in case tapered was applied before
		makeFaultSectionDataList();
		
		if(D) System.out.println("parentFaultData MoRate: "+parentFaultData.calcMomentRate(true));
		
		numSections = faultSectionDataList.size();
		
		// Apply any taper to slip rates
		switch (slipRateProfile) {
		case TAPERED:
			// sqrtSine Taper
			HistogramFunction taperedSlipCDF = FaultSystemRuptureRateInversion.getTaperedSlipFunc().getCumulativeDistFunction();
			double subSectLength = faultSectionDataList.get(0).getTraceLength();
			double totLength = subSectLength*faultSectionDataList.size();
			double currStartPoint = 0, totWt=0, moRate=0;
			for(int i=0; i<faultSectionDataList.size(); i++) {
				double x1_fract = currStartPoint/totLength;
				double x2_fract = (currStartPoint+subSectLength)/totLength;
				double wt = (taperedSlipCDF.getInterpolatedY(x2_fract) - taperedSlipCDF.getInterpolatedY(x1_fract))*faultSectionDataList.size();
				FaultSectionPrefData data = faultSectionDataList.get(i);
//				System.out.print(i+"\t"+wt+"\t"+data.getOrigAveSlipRate());

				data.setAveSlipRate(wt*data.getOrigAveSlipRate());
				data.setSlipRateStdDev(wt*data.getOrigSlipRateStdDev());
//				System.out.println("\t"+data.getOrigAveSlipRate());
				moRate += data.calcMomentRate(true);
				totWt += wt;
				currStartPoint+=subSectLength;

			}
			if(D)System.out.println("moRate after tapering slip rates = "+moRate);
			break;
		case UNI_TRIM:  // linear taper over first and last 3 subsections
			double moRateOrig = 0;
			for(FaultSectionPrefData data : faultSectionDataList) {
				moRateOrig += data.calcMomentRate(true);
			}
			FaultSectionPrefData data;
			if(NUM_SUBSECT_PER_RUP == 2) {
				data = faultSectionDataList.get(0); 
				data.setAveSlipRate(0.5*data.getOrigAveSlipRate()); 
				data.setSlipRateStdDev(0.5*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(faultSectionDataList.size()-1); 
				data.setAveSlipRate(0.5*data.getOrigAveSlipRate()); 
				data.setSlipRateStdDev(0.5*data.getOrigSlipRateStdDev());				
			}

			else if(NUM_SUBSECT_PER_RUP == 4) {
				data = faultSectionDataList.get(0); data.setAveSlipRate(0.25*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.25*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(1); data.setAveSlipRate(0.50*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.50*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(2); data.setAveSlipRate(0.75*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.75*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(faultSectionDataList.size()-1); data.setAveSlipRate(0.25*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.25*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(faultSectionDataList.size()-2); data.setAveSlipRate(0.50*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.50*data.getOrigSlipRateStdDev());
				data = faultSectionDataList.get(faultSectionDataList.size()-3); data.setAveSlipRate(0.75*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.75*data.getOrigSlipRateStdDev());				
			}
			else
				throw new RuntimeException("NUM_SUBSECT_PER_RUP="+NUM_SUBSECT_PER_RUP+" not supported");
			
			double moRateAfter = 0;
			for(FaultSectionPrefData data2 : faultSectionDataList) {
				moRateAfter += data2.calcMomentRate(true);
			}
			double wt = moRateOrig/moRateAfter;
			double totMoRate = 0;
			for(FaultSectionPrefData data2 : faultSectionDataList) {
				data2.setAveSlipRate(wt*data2.getOrigAveSlipRate());
				data2.setSlipRateStdDev(wt*data2.getOrigSlipRateStdDev());
				totMoRate += data2.calcMomentRate(true);
			}
			if(D)System.out.println("moRate after trimming uniform slip rates = "+totMoRate);
		case UNIFORM:
			break;
		}
		
		numRuptures =0;
		int curNumRup = numSections - (NUM_SUBSECT_PER_RUP-1);
		while(curNumRup>0) {
			numRuptures += curNumRup;
			curNumRup -= 1;
		}

		int rupIndex=0;

//		rupSectionMatrix = new int[numSect][numSect - (NUM_SUBSECT_PER_RUP-1)];
//		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<NUM_SUBSECT_PER_RUP+1; curSectPerRup++)
		rupSectionMatrix = new int[numSections][numRuptures];
		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<numSections+1; curSectPerRup++)
			for(int s=0;s<numSections-curSectPerRup+1;s++) {
				int firstSect = s;
				int lastSect = s+curSectPerRup-1;
//				System.out.println(rupIndex+":\t"+firstSect+"\t"+lastSect);
				for(int col=firstSect; col <= lastSect; col++)
					rupSectionMatrix[col][rupIndex] = 1;
				rupIndex += 1;
			}
		
//		for(int r=0;r<numRup;r++) {
//			System.out.print("\n");
//			for(int s=0;s<numSect;s++)
//				System.out.print(rupSectionMatrix[s][r]+"\t");
//		}
		
//		System.out.print("numSect="+numSect+"\n"+"numRup="+numRup);
//		System.exit(-1);
	}
	
	public ArrayList<FaultSectionPrefData> getFaultSectionDataList() {
		return faultSectionDataList;
	}
	
	public int[][] getRupSectionMatrix() {
		return rupSectionMatrix;
	}
	
	
	
	private void writeApriorRupRatesForMaxRateModel(String fileName) {	
		
		if(rupSectionMatrix == null)
			throw new RuntimeException("must create the rupSectionMatrix before running this method");
		
		try{
			FileWriter fw = new FileWriter(ROOT_DATA_DIR+fileName);
			int numRuptures = rupSectionMatrix[0].length;
			int numSections = rupSectionMatrix.length;

			for(int i=numSections-1;i<numRuptures; i++) {				
				fw.write(i+"\t0.0\t"+1.0+"\n");
			}
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}			
	}


	
	
	

	
	
	/**
	 * This computes and optionally saves and/or plots a hazard curve
	 * @param dirName - set as null if you don't want to save results
	 * @param popupWindow - set as true if you want plot windows to pop up
	 */
	public void writeAndOrPlotHazardCurve(FaultSystemRuptureRateInversion fltSysRupInversion, Location location, 
			double saPeriod, String dirName, boolean popupWindow, String plotTitle) {
		

		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		int numSol = fltSysRupInversion.getNumSolutions();
		if(numSol>1) {
			ArbDiscrEmpiricalDistFunc_3D curvesFromMultRunsFunc_3D = new ArbDiscrEmpiricalDistFunc_3D(hazCurveLnMin,hazCurveNum,hazCurveDelta);
			for(int i=0;i<numSol;i++) {
				EvenlyDiscretizedFunc func = computeHazardCurveLnX(fltSysRupInversion.getFaultSystemSolution(i), location, saPeriod, hazardDurationYrs);
				curvesFromMultRunsFunc_3D.set(func, 1.0);
			}
			
			EvenlyDiscretizedFunc hazCurveMeanLnX = curvesFromMultRunsFunc_3D.getMeanCurve();
			EvenlyDiscretizedFunc hazCurveMinLnX = curvesFromMultRunsFunc_3D.getMinCurve();
			EvenlyDiscretizedFunc hazCurveMaxLnX = curvesFromMultRunsFunc_3D.getMaxCurve();
			UncertainArbDiscDataset hazCurveMean95confLnX = fltSysRupInversion.get95perConfForMultRuns(curvesFromMultRunsFunc_3D);

			ArbitrarilyDiscretizedFunc hazCurveMean = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc hazCurveMin = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc hazCurveMax = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc hazCurveMeanLower95 = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc hazCurveMeanUpper95 = new ArbitrarilyDiscretizedFunc();
			for(int i=0;i<hazCurveNum;i++) {
				hazCurveMean.set(Math.exp(hazCurveMeanLnX.getX(i)), hazCurveMeanLnX.getY(i));
				hazCurveMax.set(Math.exp(hazCurveMaxLnX.getX(i)), hazCurveMaxLnX.getY(i));
				hazCurveMin.set(Math.exp(hazCurveMinLnX.getX(i)), hazCurveMinLnX.getY(i));
				hazCurveMeanLower95.set(Math.exp(hazCurveMean95confLnX.getX(i)), hazCurveMean95confLnX.getLowerY(i));
				hazCurveMeanUpper95.set(Math.exp(hazCurveMean95confLnX.getX(i)), hazCurveMean95confLnX.getUpperY(i));
			}

			hazCurveMean.setName("hazCurveMean");
			UncertainArbDiscDataset hazCurveMinMaxRange = new UncertainArbDiscDataset(hazCurveMean, hazCurveMin, hazCurveMax);
			hazCurveMinMaxRange.setName("hazCurveMinMaxRange");
			UncertainArbDiscDataset hazCurveMean95conf = new UncertainArbDiscDataset(hazCurveMean, hazCurveMeanLower95, hazCurveMeanUpper95);
			hazCurveMean95conf.setName("hazCurveMean95conf");

			plottingFuncsArray.add(hazCurveMinMaxRange);
			plottingFuncsArray.add(hazCurveMean95conf);
			plottingFuncsArray.add(hazCurveMean);

			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,255)));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,255)));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		}

		// now get the result for the mean solution
		EvenlyDiscretizedFunc curveLogXvalues = computeHazardCurveLnX(fltSysRupInversion.getFaultSystemSolution(), location, saPeriod, hazardDurationYrs);
		// convert to linear x valules
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
			for (int i = 0; i < curveLogXvalues.size(); ++i)
				curveLinearXvalues.set(Math.exp(curveLogXvalues.getX(i)), curveLogXvalues.getY(i));

		
		double twoIn50value = Math.exp(curveLogXvalues.getFirstInterpolatedX_inLogYDomain(0.02));
		double tenIn50value = Math.exp(curveLogXvalues.getFirstInterpolatedX_inLogYDomain(0.1));
		curveLinearXvalues.setInfo(imtString+"\n2in50 value: "+twoIn50value+"\n10in50 value: "+tenIn50value+
				"\nLocation: "+location.getLatitude()+", "+location.getLongitude());
		
		savedHazardCurves.add(curveLinearXvalues);
		
		// make the plot
		plottingFuncsArray.add(curveLinearXvalues);
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve"+"_"+imtString+"_"+plotTitle;
		String xAxisLabel = imtString;
		String yAxisLabel = "Probability (in "+hazardDurationYrs+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = new Range(1e-9,1.0);
		boolean logX = true;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}

	
		
	
	/**
	 */
	private EvenlyDiscretizedFunc computeHazardCurveLnX(FaultSystemSolution faultSystemSolution, Location location, 
			double saPeriod, double forecastDuration) {
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(faultSystemSolution);
		erf.setName("Simple Fault ERF");
		erf.getTimeSpan().setDuration(forecastDuration);
		if(magAareaAleatoryVariability != 0) {
			AleatoryMagAreaStdDevParam param = (AleatoryMagAreaStdDevParam) erf.getParameter(AleatoryMagAreaStdDevParam.NAME);
			param.setValue(magAareaAleatoryVariability);
		}
		erf.updateForecast();		// update forecast

		// write out parameter values
		if(D) {
			ParameterList paramList = erf.getAdjustableParameterList();
			for(int i=0;i<paramList.size(); i++) {
				Parameter<?> param = paramList.getByIndex(i);
				System.out.println(param.getName()+"\t"+param.getValue());
			}			
		}
		if(D) System.out.println("NumFaultSystemSources = "+erf.getNumFaultSystemSources());

		
		// chose attenuation relationship (GMPE)
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		imr.setParamDefaults();
		
		EvenlyDiscretizedFunc curveLogXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		
		// make the site object and set values
		Site site = new Site(location);
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
//			System.out.println(param.getName()+"\t"+param.getValue());
		}
		
		// set the IMT
		ArbitrarilyDiscretizedFunc curveLinearXvalues;
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
		}

	
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
		
		return curveLogXvalues;
	}


	/**
	 * This makes the following plots as a function of section index: slip rate, event rate, ave slip/COV,
	 * incremental participation MFD and cumulative participation MFD.  This also combines the first three
	 * in one plot, and the last two in one plot.
	 * @param fltSysRupInversion
	 */
	public void makeResultsVsSectionPlots(FaultSystemRuptureRateInversion fltSysRupInversion) {
		// functions versus subsection index
		ArrayList<PlotSpec> specList = new ArrayList<PlotSpec>();
		Range xAxisRange = new Range(-1,numSections);
		ArrayList<Range> xRanges = new ArrayList<Range>();
		xRanges.add(xAxisRange);
		ArrayList<Range> yRanges = new ArrayList<Range>();

		Range yAxisRange1 = new Range(0,50);
		PlotSpec spec1 = fltSysRupInversion.writeAndOrPlotSlipRates(dirName, popUpPlots, xAxisRange, yAxisRange1, 6.5, 1.7);
		Range yAxisRange2 = new Range(0,0.015);
		PlotSpec spec2 = fltSysRupInversion.writeAndOrPlotEventRates(dirName, popUpPlots, xAxisRange, yAxisRange2, 6.5, 1.7);
		Range yAxisRange3 = new Range(0,5);
		PlotSpec spec3 = fltSysRupInversion.writeAndOrPlotAveSlipAndCOV(dirName, popUpPlots, xAxisRange, yAxisRange3, 6.5, 1.7);
		
		yRanges.add(yAxisRange1);
		yRanges.add(yAxisRange2);
		yRanges.add(yAxisRange3);
		specList.add(spec1);
		specList.add(spec2);
		specList.add(spec3);
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(9);
		gp.setAxisLabelFontSize(11);
		gp.setBackgroundColor(Color.WHITE);
		gp.drawGraphPanel(specList, false, false, xRanges, yRanges);
		int width = (int)(6.85*72.);
		int height = (int)(specList.size()*1.7*72.);
		gp.getChartPanel().setSize(width, height); 
		String fileNamePrefix = dirName+"/allFuncsVsSectionIndex";
		try {
			gp.saveAsPNG(fileNamePrefix+".png");
			gp.saveAsPDF(fileNamePrefix+".pdf");
//			gp.saveAsTXT(fileNamePrefix+".txt");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		// Now do participation MFDs
		double yAddOn = FaultSystemRuptureRateInversion.MAG_DELTA/2;
		Range yAxisRange4 = new Range(fltSysRupInversion.minMagMFD_WithAleatory-yAddOn,fltSysRupInversion.maxMagMFD_WithAleatory+yAddOn);
		ArrayList<XYZPlotSpec> xyzSpecList = fltSysRupInversion.writeAndOrPlotSectPartMFDs(dirName, popUpPlots,6.5, 2.5, xAxisRange, yAxisRange4);
		ArrayList<Range> yRanges2 = new ArrayList<Range>();
		yRanges2.add(yAxisRange4);
		yRanges2.add(yAxisRange4);
		PlotPreferences plotPrefs = XYZGraphPanel.getDefaultPrefs();
		plotPrefs.setTickLabelFontSize(9);
		plotPrefs.setAxisLabelFontSize(11);
		plotPrefs.setBackgroundColor(Color.WHITE);
		XYZGraphPanel gp_xyz = new XYZGraphPanel(plotPrefs);
		gp_xyz.drawPlot(xyzSpecList, false, false, xRanges, yRanges2);
		width = (int)(6.68*72.);
		height = (int)(5.25*72.);
		gp_xyz.getChartPanel().setSize(width, height); 
		fileNamePrefix = dirName+"/allPartMFDvsSectionIndex";
		try {
			gp_xyz.saveAsPNG(fileNamePrefix+".png");
			gp_xyz.saveAsPDF(fileNamePrefix+".pdf");
//			gp_xyz.saveAsTXT(fileNamePrefix+".txt");
		} catch (IOException e) {
			e.printStackTrace();
		}


	}
	
	
	/**
	 * This is makes to hazard maps, one for 10% in 50-year and one for 2% in 50-year ground motions.  
	 * For saPeriod of 1.0 or 0.2, this also make RTGM maps.
	 * @param fltSysRupInversion
	 * @param saPeriod - set as 0 for PGA; this will crash if SA period is not supported
	 * @param dirName
	 * @param popupWindow
	 */
	public void makeHazardMaps(FaultSystemRuptureRateInversion fltSysRupInversion, double saPeriod, String hazDirName, boolean popupWindow) {
		
		if(D) System.out.println("Making hazard map");
		
		File hazDirFile = null;
		if(hazDirName != null) {
			hazDirFile = new File(hazDirName);
			hazDirFile.mkdirs();
		}


		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		// the following makes RTGM if sa period is 1.0 or 0.2
		ArrayList<GriddedGeoDataSet> griddedDataList = makeHazardMapGriddedData(fltSysRupInversion.getFaultSystemSolution(), saPeriod, true);
		
		GriddedGeoDataSet rtgmGriddedDataSetArray=null;
		if(griddedDataList.size()>2)
			rtgmGriddedDataSetArray = griddedDataList.get(2);

		Region region = getGriddedRegion();
		
		try {
			for(int p=0;p<hazardProbNameArray.length;p++) {
				CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(griddedDataList.get(p).getMinZ(), griddedDataList.get(p).getMaxZ());
//				CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(gridData.getMinZ(), gridData.getMaxZ());
				ArrayList<LocationList> faults = FaultBasedMapGen.getTraces(faultSectionDataList);
				double[] values = new double[faults.size()];
				for(int i=0;i<values.length;i++)
					values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;
				
				String label = imtString+"_"+hazardProbNameArray[p];
								
				GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, griddedDataList.get(p), hazGridSpacing, region, true, label);
//				map.setGenerateKML(true); // this tells it to generate a KML file that can be loaded into Google Earthe
				
				// override default trace width
				for(PSXYPolygon fltTracePoly: map.getPolys())
					fltTracePoly.setPenWidth(1);
				
				// remove coast and poitical boundaries
				map.setCoast(null);
				
				try {
					FaultBasedMapGen.SAVE_ZIPS=true;
					FaultBasedMapGen.plotMap(hazDirFile, label, popupWindow, map);	
				} catch (GMT_MapException e) {
					e.printStackTrace();
				}
				
				// write out values to a text file
				if(hazDirName != null) {
					String fileName = hazDirName +"/"+label+".txt";
					writeHazardMapGriddedData(griddedDataList.get(p), fileName);
				}
			}
			
			if(rtgmGriddedDataSetArray != null) {
				CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(rtgmGriddedDataSetArray.getMinZ(), rtgmGriddedDataSetArray.getMaxZ());
				ArrayList<LocationList> faults = FaultBasedMapGen.getTraces(faultSectionDataList);
				double[] values = new double[faults.size()];
				for(int i=0;i<values.length;i++)
					values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;
				
				String label = imtString+"_RTGM";
				
				GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, rtgmGriddedDataSetArray, hazGridSpacing, region, true, label);
//				map.setGenerateKML(true); // this tells it to generate a KML file that can be loaded into Google Earthe
				
				// override default trace width
				for(PSXYPolygon fltTracePoly: map.getPolys())
					fltTracePoly.setPenWidth(1);
				
				// remove coast and poitical boundaries
				map.setCoast(null);

				try {
					FaultBasedMapGen.SAVE_ZIPS=true;
					FaultBasedMapGen.plotMap(hazDirFile, label, popupWindow, map);	
				} catch (GMT_MapException e) {
					e.printStackTrace();
				}
				
				// write out values to a text file
				if(hazDirName != null) {
					String fileName = hazDirName +"/"+label+".txt";
					writeHazardMapGriddedData(rtgmGriddedDataSetArray, fileName);
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void writeHazardMapGriddedData(GriddedGeoDataSet hazData, String fileName) {
		FileWriter fw;
		try {
			fw = new FileWriter(fileName);
			fw.write("index\tvalue\tlatitude\tlongitude\n");
			for(int i=0;i<hazData.size(); i++)	 {
				Location loc = hazData.getLocation(i);
				fw.write(i+"\t"+hazData.get(i)+"\t"+loc.getLatitude()+"\t"+loc.getLongitude()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * 
	 * @param fltSysRupInversion
	 * @param saPeriod
	 * @param hazDirName
	 * @param includeRTGRM
	 * @return index 0 has 2in50, index 1 has 10in50, and index 2 has RTGM if requested (and period =1 or 0.2)
	 */
	public ArrayList<GriddedGeoDataSet[]> makeHazardMapGriddedDataForMultRuns(FaultSystemRuptureRateInversion fltSysRupInversion, 
			double saPeriod, String hazDirName, boolean includeRTGRM) {

		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		String fileNamePrefix2in50=null, fileNamePrefix10in50=null,fileNamePrefixRTGM=null;
		if(hazDirName != null) {
			fileNamePrefix2in50 = hazDirName +"/"+imtString+"_"+hazardProbNameArray[0];
			fileNamePrefix10in50 = hazDirName +"/"+imtString+"_"+hazardProbNameArray[1];
			fileNamePrefixRTGM = hazDirName +"/"+imtString+"_RTGM";
		}
				
		int numSim = fltSysRupInversion.getNumSolutions();
		ArrayList<GriddedGeoDataSet[]> resultList = new ArrayList<GriddedGeoDataSet[]>();
		resultList.add(new GriddedGeoDataSet[numSim]);	// 2in50 data
		resultList.add(new GriddedGeoDataSet[numSim]);	// 10in50 data
		if((saPeriod==1 || saPeriod==0.2) && includeRTGRM)
			resultList.add(new GriddedGeoDataSet[numSim]);  // RTGM

		for(int s=0;s<numSim;s++) {
			if(D) System.out.println("Working on Hazard Map Data #"+s);
			ArrayList<GriddedGeoDataSet> dataList = makeHazardMapGriddedData(fltSysRupInversion.getFaultSystemSolution(s), saPeriod, includeRTGRM);
			resultList.get(0)[s] = dataList.get(0);
			resultList.get(1)[s] = dataList.get(1);
			if(hazDirName != null) {	// save results to file
				writeHazardMapGriddedData(dataList.get(0), fileNamePrefix2in50+"_"+s+".txt");
				writeHazardMapGriddedData(dataList.get(1), fileNamePrefix10in50+"_"+s+".txt");
			}
			if(dataList.size()>2) {
				resultList.get(2)[s] = dataList.get(2);
				if(hazDirName != null)
					writeHazardMapGriddedData(dataList.get(2), fileNamePrefixRTGM+"_"+s+".txt");
			}
		}

		return resultList;
	}
	 
	 
	 public void plotNormPDF_FromMultHazMaps(GriddedGeoDataSet[] griddedHazDataArray, String label, String hazDirName, boolean popupWindow) {
		 
		int numMaps = griddedHazDataArray.length;

		ArrayList<HistogramFunction> ratioHistogramList = new ArrayList<HistogramFunction>();
		for(int j=0;j<numMaps;j++)
			ratioHistogramList.add(new HistogramFunction(0.51,80,0.02));
		HistogramFunction ratioHistogram = new HistogramFunction(0.51,80,0.02);
		ArbDiscrEmpiricalDistFunc ratioData = new ArbDiscrEmpiricalDistFunc();
		int numLocs = griddedHazDataArray[0].size();
		double minRatio = Double.MAX_VALUE;
		double maxRatio = 0.0;
		double fractMoreThan10percASway = 0;
		double[] meanArray = new double[numLocs];
		for(int i=0;i<numLocs;i++) {
			double mean = 0;
			for(int j=0;j<numMaps;j++) 
				mean += griddedHazDataArray[j].get(i);
			mean /= (double)numMaps;
			meanArray[i] = mean;
			for(int j=0;j<numMaps;j++) {
				double ratio = griddedHazDataArray[j].get(i)/mean;
				ratioHistogram.add(ratio, 1.0);
				ratioHistogramList.get(j).add(ratio, 1.0);
				ratioData.set(ratio,1.0);
				if(minRatio>ratio) minRatio=ratio;
				if(maxRatio<ratio) maxRatio=ratio;
				if(ratio>1.1 || ratio<0.9) {
					fractMoreThan10percASway += 1;
				}
			}
		}
		fractMoreThan10percASway /= ratioData.size();
		
		// find which dataset is closest to the mean
		int closestIndex = -1;
		double closestVal = Double.MAX_VALUE;
		for(int j=0;j<numMaps;j++) {
			double sum = 0;
			for(int i=0;i<numLocs;i++) 
				sum += Math.abs((griddedHazDataArray[j].get(i) - meanArray[i])/meanArray[i]);
			sum /= (double)numLocs;
			if(sum<closestVal) {
				closestVal=sum;
				closestIndex = j;
			}
		}

		
		double numData = ratioHistogram.calcSumOfY_Vals();
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		ratioHistogram.setName(label+"Ratio Histogram");
		ratioHistogram.setInfo("\nNum Data = "+(float)numData+"\nCOV from Hist = "+(float)ratioHistogram.computeCOV()+
		"\nCOV from Data = "+(float)ratioData.getCOV()+"\nminRatio = "+(float)minRatio+"\nmaxRatio = "+(float)maxRatio+
		"\nFraction More Than 10% Away From 1.0: "+(float)fractMoreThan10percASway+"\n Data set "+closestIndex+
		" is closest to mean (ave abs norm diff = "+(float)closestVal+")");
		ratioHistogram.normalizeBySumOfY_Vals();
		ratioHistogram.scale(1.0/ratioHistogram.getDelta());
		funcs.add(ratioHistogram);
		
		if(D) System.out.println(label+ ": data set "+closestIndex+" is closest to mean");
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		
		DefaultXY_DataSet tenPercentDiffLines = new DefaultXY_DataSet();
		tenPercentDiffLines.set(0.9,ratioHistogram.getMaxY());
		tenPercentDiffLines.set(0.9,0.0);
		tenPercentDiffLines.set(1.1,0.0);
		tenPercentDiffLines.set(1.1,ratioHistogram.getMaxY());
		PlotCurveCharacterstics tenPercentDiffLinesPlotChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK);
		funcs.add(tenPercentDiffLines);
		plotChars.add(tenPercentDiffLinesPlotChar);

		double minX=0.5;
		double maxX=1.5;

		String fileNamePrefix = null;
		if(hazDirName != null)
			fileNamePrefix = hazDirName+"/"+label+"_NormPDF";
		String plotName =label+"_NormPDF";
		String xAxisLabel = "Ratio";
		String yAxisLabel = "PDF";
		Range xAxisRange = new Range(minX,maxX);
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		ArrayList<XY_DataSet> funcs2 = new ArrayList<XY_DataSet>();
		ArbitrarilyDiscretizedFunc cumDist = ratioData.getCumDist();
		cumDist.scale(1.0/cumDist.getMaxY());
		funcs2.add(cumDist);
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		fileNamePrefix = null;
		if(hazDirName != null)
			fileNamePrefix = hazDirName+"/"+label+"_NormCDF";
		plotName =label+"_NormCDF";
		xAxisLabel = "Ratio";
		yAxisLabel = "CDF";
		PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars2, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		
		
		// stacked histograms
		ArrayList<Color> colorList = new ArrayList<Color>();
		colorList.add(new Color(128,0,128)); colorList.add(Color.BLUE); colorList.add(Color.CYAN); colorList.add(Color.GREEN); colorList.add(Color.LIGHT_GRAY); 
		colorList.add(Color.YELLOW); colorList.add(Color.ORANGE); colorList.add(Color.RED); colorList.add(Color.MAGENTA); 
		colorList.add(Color.DARK_GRAY); colorList.add(Color.GRAY); 
		
		HistogramFunction testRatioHistogram = new HistogramFunction(0.51,80,0.02);
		ArrayList<XY_DataSet> funcs3 = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars3 = new ArrayList<PlotCurveCharacterstics>();
		for(int j=0;j<numMaps;j++) {
			ratioHistogramList.get(j).scale(1.0/(numData));;
			plotChars3.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorList.get(j)));
			for(int k=0;k<testRatioHistogram.size();k++)
				testRatioHistogram.add(k, ratioHistogramList.get(j).getY(k));
		}
		
		for(HistogramFunction func : HistogramFunction.getStackedHists(ratioHistogramList, true)) {
			func.scale(1.0/func.getDelta());
			funcs3.add(func);
		}
		 
			funcs3.add(tenPercentDiffLines);
			plotChars3.add(tenPercentDiffLinesPlotChar);

		
//		for(int k=0;k<testRatioHistogram.size();k++)
//			System.out.println((float)testRatioHistogram.getY(k)+"\t"+(float)ratioHistogram.getY(k));

		fileNamePrefix = null;
		if(hazDirName != null)
			fileNamePrefix = hazDirName+"/"+label+"_NormPDF_Stacked";
		plotName =label+"_NormPDF";
		xAxisLabel = "Ratio";
		yAxisLabel = "PDF";
		xAxisRange = new Range(minX,maxX);
		yAxisRange = null;
		logX = false;
		logY = false;

		PlottingUtils.writeAndOrPlotFuncs(funcs3, plotChars3, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
	 }
	
	
	/**
	 * This is makes to hazard maps, one for 10% in 50-year and one for 2% in 50-year ground motions.  For saPeriod of 1.0 or 0.2,
	 * this also make RTGM maps.
	 * @param fltSysRupInversion
	 * @param saPeriod - set as 0 for PGA; this will crash if SA period is not supported
	 * @param includeRTGRM
	 * @return - array with 2in50 and then 10in50 data, plus RTGM if saPeriod is 1.0 or 0.2
	 */
	public ArrayList<GriddedGeoDataSet> makeHazardMapGriddedData(FaultSystemSolution fltSysSolution, double saPeriod, boolean includeRTGRM) {
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		if(D) System.out.println("Making hazard map data for "+imtString);

		
		//map region
		GriddedRegion region = getGriddedRegion();
		
		// make gridded data sets
		GriddedGeoDataSet[] griddedDataSetArray = new GriddedGeoDataSet[hazardProbArray.length];
		for(int i=0;i<hazardProbArray.length;i++) {
			griddedDataSetArray[i] = new GriddedGeoDataSet(region, true);
		}
		
		GriddedGeoDataSet rtgmGriddedDataSetArray=null;
		if((saPeriod == 1.0 || saPeriod == 0.2) && includeRTGRM) // 1 or 5 Hz SA
			rtgmGriddedDataSetArray = new GriddedGeoDataSet(region, true);
		
		// Create the ERF
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fltSysSolution);
		erf.setName("Simple Fault ERF");
		
		// set the forecast duration
		erf.getTimeSpan().setDuration(hazardDurationYrs);
		
//		FaultGridSpacingParam gridSpaceParam = (FaultGridSpacingParam) erf.getParameter(FaultGridSpacingParam.NAME);
//		gridSpaceParam.setValue(0.2);
		
		// set magAareaAleatoryVariability
		if(magAareaAleatoryVariability != 0) {
			AleatoryMagAreaStdDevParam param = (AleatoryMagAreaStdDevParam) erf.getParameter(AleatoryMagAreaStdDevParam.NAME);
			param.setValue(magAareaAleatoryVariability);
		}

		
		// update forecast
		erf.updateForecast();
		
		// set the attenuation relationship (GMPE)
		imr.setParamDefaults();
		
		// set the IMT & curve x-axis values
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
		
		// adding this to be consistent with computeHazardCurveLnX (and I commented out setting curveLinearXvalues below)
		EvenlyDiscretizedFunc tempCurveLogXvalues = new EvenlyDiscretizedFunc(hazCurveLnMin,hazCurveNum,hazCurveDelta);
		for(int i=0;i<tempCurveLogXvalues.size();i++)
			curveLinearXvalues.set(Math.exp(tempCurveLogXvalues.getX(i)),1.0);

		
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
//			curveLinearXvalues = IMT_Info.getUSGS_PGA_Function();
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
//			curveLinearXvalues = IMT_Info.getUSGS_SA_01_AND_02_Function();
		}
		
		// make the site object and set values
		Site site = new Site(new Location(40.75, -111.90));	// Location will get over written
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
			System.out.println(param.getName()+"\t"+param.getValue());
		}
		
		ArbitrarilyDiscretizedFunc curveLogXvalues = HazardCurveSetCalculator.getLogFunction(curveLinearXvalues); // this is what the calculator expects
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		// loop over sites
		int counter = -1;
		for(int i=0;i<griddedDataSetArray[0].size();i++) {
			counter+=1;
			if(counter == 500) {
				System.out.println(i+" of "+griddedDataSetArray[0].size());
				counter = 0;
			}
			site.setLocation(griddedDataSetArray[0].getLocation(i));
			calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
			// get the points on the hazard curve
			for(int p=0;p<hazardProbArray.length;p++) {
				if(hazardProbArray[p] >= curveLogXvalues.getMinY() && hazardProbArray[p] <= curveLogXvalues.getMaxY()) {
//					double logVal = curveLogXvalues.getFirstInterpolatedX(hazardProbArray[p]);
					// this is more accurate:
					double logVal = curveLogXvalues.getFirstInterpolatedX_inLogYDomain(hazardProbArray[p]);
					griddedDataSetArray[p].set(i, Math.exp(logVal));					
				}
				else {
					griddedDataSetArray[p].set(i,0);					
					System.out.println("prob outside range at i="+i+":\t"+griddedDataSetArray[0].getLocation(i));
//					System.out.println(curveLogXvalues);
//					System.exit(-1);
				}
			}
			// do RTGM is SA is 1 or 5 Hz
			if(rtgmGriddedDataSetArray != null) {
				// get annual rate curve with linear x-axis values
				boolean hasInfValue=false;
				for(int j=0;j<curveLogXvalues.size();j++) {
					// convert probability to annual rate
					double rate = -Math.log(1.0-curveLogXvalues.getY(j))/hazardDurationYrs;
					curveLinearXvalues.set(j, rate);
					if(Double.isInfinite(rate)) {
						hasInfValue = true;
//						System.out.println("Infinite Rate from: "+curveLogXvalues.getY(j));
					}
				}
				// remove Inf values that can occur if above curveLogXvalues.getY(j) = 1.0;
				ArbitrarilyDiscretizedFunc curveLinearXvaluesCleaned;
				if(hasInfValue) {
					curveLinearXvaluesCleaned = new ArbitrarilyDiscretizedFunc();
					for(int j=0;j<curveLinearXvalues.size();j++)
						if(!Double.isInfinite(curveLinearXvalues.getY(j)))
								curveLinearXvaluesCleaned.set(curveLinearXvalues.getX(j),curveLinearXvalues.getY(j));
				}
				else
					curveLinearXvaluesCleaned = curveLinearXvalues;
				// now get RTGM
				Frequency freq;
				if(saPeriod == 1.0)
					freq = Frequency.SA_1P00;
				else
					freq = Frequency.SA_0P20;
				RTGM rtgm;
				
				// this is to write out test value to compare with the USGS on-line calculator
				try {
//					rtgm = RTGM.create(curveLinearXvaluesCleaned, freq, 0.8).call();
					rtgm = RTGM.create(curveLinearXvaluesCleaned, null, null).call();
					rtgmGriddedDataSetArray.set(i, rtgm.get());
					// this is to write one out to test against USGS web calculator; it checked out
//					if(i==0) {
//						System.out.println("TEST VALUES\n\n");
//						for(int j=0;j<curveLinearXvalues.size();j++)
//							System.out.print(curveLinearXvalues.getX(j)+",");
//						System.out.println("\n\n");
//						for(int j=0;j<curveLinearXvalues.size();j++)
//							System.out.print((float)curveLinearXvalues.getY(j)+",");
//						System.out.println("\n\nTest RTGM = "+rtgm.get());
//					}
				} catch (Exception e) {
					System.out.println("RTGM Error; Hazard curve is:\n"+curveLinearXvaluesCleaned);
					System.out.println("Location is:\n"+griddedDataSetArray[0].getLocation(i));
					e.printStackTrace();
					System.exit(0);
				}
			}
		}
		
		ArrayList<GriddedGeoDataSet> returnList = new ArrayList<GriddedGeoDataSet>();
		for(GriddedGeoDataSet dataSet : griddedDataSetArray)
			returnList.add(dataSet);
		if(rtgmGriddedDataSetArray != null) {
			returnList.add(rtgmGriddedDataSetArray);
		}
		return returnList;
	}

	
	private GriddedGeoDataSet[] readMultipleHazardMapGriddedData(ArrayList<String> fileNameList) {
		GriddedGeoDataSet[] dataArray = new GriddedGeoDataSet[fileNameList.size()];
		for(int i=0; i<fileNameList.size();i++)
			dataArray[i] = readHazardMapDataFromFile(fileNameList.get(i));
		return dataArray;
	}
	
	
	/**
	 * 
	 * @param fileName1 - numerator data
	 * @param fileName2 - denominator data
	 * @param label - label for the plot and file name prefix
	 * @param dirName - directory (full path); this cannot be null
	 * @param popupWindow - whether to show results in a pop up window
	 */
	public void makeHazardMapRatio(String fileName1, String fileName2, String label, String dirName, boolean popupWindow) {
		
		GriddedGeoDataSet data1 = readHazardMapDataFromFile(fileName1);
		GriddedGeoDataSet data2 = readHazardMapDataFromFile(fileName2);
		HistogramFunction ratioHistogram = new HistogramFunction(0.0,200,0.05);

		
		double[] ratioArray = new double[data1.size()];
		double fractMoreThan10percASway = 0;
		double aveAbsValDiff = 0;
		double min = Double.MAX_VALUE;
		double max = 0;
		double distFrom1 = Double.MAX_VALUE;
		Location maxLoc=null, minLoc=null, equalLoc=null;
		for(int i=0;i<data1.size();i++) {
			double ratio = data1.get(i)/data2.get(i);
//			if(ratio < ratioHistogram.getMaxX())
			ratioHistogram.add(ratio, 1.0);
			// put log10 ratio in first data set:
			data1.set(i, Math.log10(ratio));
			ratioArray[i] = ratio;
			aveAbsValDiff += Math.abs(ratio-1.0)/ratioArray.length;
			if(ratio>1.1 || ratio<0.9) {
				fractMoreThan10percASway += 1;
			}
			if(min>ratio) {
				min=ratio;
				minLoc = data1.getLocation(i);
			}
			if(max<ratio) {
				max=ratio;
				maxLoc = data1.getLocation(i);
			}
			if(Math.abs(ratio-1.0)<distFrom1) {
				distFrom1 = Math.abs(ratio-1);
				equalLoc = data1.getLocation(i);
			}
		}
		fractMoreThan10percASway /= (double)data1.size();
		double mean = StatUtils.mean(ratioArray);
		double cov = Math.sqrt(StatUtils.populationVariance(ratioArray)) / mean;
		aveAbsValDiff /= mean; // normalize so it's fractional diff
				
		try {
//			CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(-0.48, 0.48); // factor of 3
			CPT cpt = new CPT();
			Color purple = new Color(128,0,128);
			Color orange = new Color(255,165,0); // this looks better than the default
			cpt.add(new CPTVal((float)Math.log10(1.0/5.0), purple, (float)Math.log10(1.0/3.0), Color.BLUE));
			cpt.add(new CPTVal((float)Math.log10(1.0/3.0), Color.BLUE, (float)Math.log10(1.0/2.0), Color.CYAN));
			cpt.add(new CPTVal((float)Math.log10(1.0/2.0), Color.CYAN, (float)Math.log10(1.0/1.1), Color.GREEN));
			cpt.add(new CPTVal((float)Math.log10(1.0/1.1), Color.GREEN, (float)Math.log10(1/1.05), Color.LIGHT_GRAY));
			cpt.add(new CPTVal((float)Math.log10(1/1.05), Color.LIGHT_GRAY, (float)Math.log10(1.05), Color.LIGHT_GRAY));
//			cpt.add(new CPTVal((float)Math.log10(1.0/1.1), Color.GREEN, (float)Math.log10(1.0), Color.LIGHT_GRAY));
//			cpt.add(new CPTVal((float)Math.log10(1.0), Color.LIGHT_GRAY, (float)Math.log10(1.1), Color.YELLOW));
			cpt.add(new CPTVal((float)Math.log10(1.05), Color.LIGHT_GRAY, (float)Math.log10(1.1), Color.YELLOW));
			cpt.add(new CPTVal((float)Math.log10(1.1), Color.YELLOW, (float)Math.log10(2.0), orange));
			cpt.add(new CPTVal((float)Math.log10(2.0), orange, (float)Math.log10(3.0), Color.RED));
			cpt.add(new CPTVal((float)Math.log10(3.0), Color.RED, (float)Math.log10(5.0), Color.MAGENTA));
			cpt.setBelowMinColor(purple);
			cpt.setAboveMaxColor(Color.MAGENTA);
			
			ArrayList<LocationList> faults = FaultBasedMapGen.getTraces(faultSectionDataList);
			double[] values = new double[faults.size()];
			for(int i=0;i<values.length;i++)
				values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;

			GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, data1, hazGridSpacing, getGriddedRegion(), true, "Log10 "+label);

			// override default trace width
			for(PSXYPolygon fltTracePoly: map.getPolys())
				fltTracePoly.setPenWidth(1.1);
			
			// remove coast and political boundaries
			map.setCoast(null);
			
			// I din't see any effect for the following
//			map.setInterpSettings(GMT_InterpolationSettings.getDefaultSettings());

			try {
				FaultBasedMapGen.SAVE_ZIPS=true;
				FaultBasedMapGen.plotMap(new File(dirName), label, popupWindow, map);	
			} catch (GMT_MapException e) {
				e.printStackTrace();
			}

			// write out values to a text file
			FileWriter fw = new FileWriter(dirName+"/"+label+".txt");
			fw.write("index\tvalue\tlatitude\tlongitude\n");
			for(int i=0;i<data1.size(); i++)	 {
				Location loc = data1.getLocation(i);
				fw.write(i+"\t"+data1.get(i)+"\t"+loc.getLatitude()+"\t"+loc.getLongitude()+"\n");
			}
			fw.close();
			
			ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
			ratioHistogram.setName("ratioHistogram");
			ratioHistogram.setInfo("Num Data = "+(float)ratioHistogram.calcSumOfY_Vals()+"\nMean = "+mean+
					"\nFraction More Than 10% Away From 1.0: "+
					fractMoreThan10percASway+"\nCOV = "+(float)cov+
					"\nAveAbsValDiff (mean normalized) = "+(float)aveAbsValDiff+
					"\nminRatio = "+(float)min+" at lat/lon: "+minLoc.getLatitude()+", "+minLoc.getLongitude()+
					"\nmaxRatio = "+(float)max+" at lat/lon: "+maxLoc.getLatitude()+", "+maxLoc.getLongitude()+
					"\nLoc with ratio closest to 1.0: "+equalLoc.getLatitude()+", "+equalLoc.getLongitude());
			ratioHistogram.normalizeBySumOfY_Vals();
			funcs.add(ratioHistogram);
			
			ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
			
			double minX=0.0;
			double maxX=2.0;

			String fileNamePrefix = null;
			if(dirName != null)
				fileNamePrefix = dirName+"/"+label+"_Histogram";
			String plotName =label+"_Histogram";
			String xAxisLabel = "Ratio";
			String yAxisLabel = "PDF";
			Range xAxisRange = new Range(minX,maxX);
			Range yAxisRange = null;
			boolean logX = false;
			boolean logY = false;

			PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, 
					xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
			

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

	// this is the region for hazard maps
	private GriddedRegion getGriddedRegion() {
//		Location loc1 = parentFaultData.getFaultTrace().get(0);
//		Location loc2 = parentFaultData.getFaultTrace().get(parentFaultData.getFaultTrace().size()-1);
//		System.out.println(loc1+"\n"+loc2);
//		System.exit(-1);
		
		// spatial buffer is 1 degree
		Location firstLoc = parentFaultData.getFaultTrace().get(0);
		Location lastLoc = parentFaultData.getFaultTrace().get(parentFaultData.getFaultTrace().size()-1);
		double minLat = firstLoc.getLatitude()-1-hazGridSpacing/2.0;
		double maxLat = lastLoc.getLatitude()+1-hazGridSpacing/2.0;
		double minLon = firstLoc.getLongitude()-1;
		double maxLon = firstLoc.getLongitude()+1;
		Location regionCorner1 = new Location(maxLat, minLon);
		Location regionCorner2 = new Location(minLat, maxLon);
		GriddedRegion region = new GriddedRegion(regionCorner1, regionCorner2, hazGridSpacing, hazGridSpacing, null);
		if(D) System.out.println("num points in gridded region: "+region.getNumLocations());
		return region;
	}
	
	
	public double getMinPossibleRate(ScalingRelationshipEnum scalingRel) {
		double origWidthKm = faultSectionDataList.get(0).getOrigDownDipWidth();
		double maxRupLengthKm = parentFaultData.getTraceLength();
		double maxRupAreaKmSq = parentFaultData.getReducedDownDipWidth()*maxRupLengthKm;
		double mag = scalingRel.getMag(maxRupAreaKmSq*1e6, maxRupLengthKm*1e3, origWidthKm*1e3);
		return parentFaultData.calcMomentRate(true)/MagUtils.magToMoment(mag);
	}
	
	
	public double getMaxPossibleRate(ScalingRelationshipEnum scalingRel) {
		double origWidthKm = faultSectionDataList.get(0).getOrigDownDipWidth();
		double minRupLengthKm = faultSectionDataList.get(0).getTraceLength();
		double minRupAreaKmSq = faultSectionDataList.get(0).getReducedDownDipWidth()*NUM_SUBSECT_PER_RUP*minRupLengthKm;
		double mag = scalingRel.getMag(minRupAreaKmSq*1e6, minRupLengthKm*1e3, origWidthKm*1e3);
		return parentFaultData.calcMomentRate(true)/MagUtils.magToMoment(mag);
	}
	
	
	public void plotRupLengthForEllsworthB() {
		
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();

		double y = 0.0;
		for(double mag=6.4;mag<7.9; mag+=0.1) {
			double length = Math.pow(10, mag-4.2)/14.0;
			y-=1;
			DefaultXY_DataSet func = new DefaultXY_DataSet();
			func.set(0.0,y);
			func.set(length,y);
			func.set(0.0,y);
			func.setName((float)mag+" has length "+length);
			funcs.add(func);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		}
		Range xRange = new Range(0,406);
		Range yRange = new Range(y-1, 1);
		
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Length", "Rup", xRange, yRange, 
					false, false, 3.5, 3.0, ROOT_PATH+"RupLengthForEllsworthB", true);	
	}
	
	
	/**
	 * This plots the classic MFDs for the paper
	 */
	public void plotClassicMFDs() {
		initData();
		IncrementalMagFreqDist grMFD = this.getTargetMFD(ScalingRelationshipEnum.ELLSWORTH_B,  MFD_TargetType.GR_b_1pt0);
		
		GaussianMagFreqDist charMFD = new GaussianMagFreqDist(grMFD.getMinX(), grMFD.size()+3, grMFD.getDelta(), 
				grMFD.getMaxX(), 0.12, grMFD.getTotalMomentRate(), 2.0, 2);
		
		// hard coded from GR results:
		IncrementalMagFreqDist gr_AveSectPartMFD = new IncrementalMagFreqDist(charMFD.getMinX(),charMFD.size(),charMFD.getDelta());
		gr_AveSectPartMFD.set(6.45,7.6359726E-4);	
		gr_AveSectPartMFD.set(6.55,7.581836E-4);
		gr_AveSectPartMFD.set(6.65,7.2269596E-4);
		gr_AveSectPartMFD.set(6.75,7.1735383E-4);
		gr_AveSectPartMFD.set(6.85,7.2180794E-4);
		gr_AveSectPartMFD.set(6.95,6.940851E-4);
		gr_AveSectPartMFD.set(7.05,6.947155E-4);
		gr_AveSectPartMFD.set(7.15,7.041709E-4);
		gr_AveSectPartMFD.set(7.25,6.952328E-4);
		gr_AveSectPartMFD.set(7.35,6.8414275E-4);
		gr_AveSectPartMFD.set(7.45,6.8566634E-4);
		gr_AveSectPartMFD.set(7.55,6.881974E-4);
		gr_AveSectPartMFD.set(7.65,6.8375067E-4);
		gr_AveSectPartMFD.set(7.75,6.837666E-4);
		gr_AveSectPartMFD.set(7.85,6.808676E-4);
		gr_AveSectPartMFD.set(7.95,6.439217E-4);
		
		IncrementalMagFreqDist char_AveSectPartMFD = charMFD.deepClone();
		
		IncrementalMagFreqDist combinedMFD = new IncrementalMagFreqDist(charMFD.getMinX(),charMFD.size(),charMFD.getDelta());
		IncrementalMagFreqDist combined_AveSectPartMFD = new IncrementalMagFreqDist(charMFD.getMinX(),charMFD.size(),charMFD.getDelta());
		for(int i=0;i<combinedMFD.size();i++) {
			if(i<grMFD.size()) {
				combinedMFD.set(i, 0.333*grMFD.getY(i)+0.667*charMFD.getY(i));
				combined_AveSectPartMFD.set(i, 0.333*gr_AveSectPartMFD.getY(i)+0.667*char_AveSectPartMFD.getY(i));
			}
			else {
				combinedMFD.set(i, 0.667*charMFD.getY(i));
				combined_AveSectPartMFD.set(i, 0.667*char_AveSectPartMFD.getY(i));
			}
		}
		
		grMFD.setName("NSHM GR b=1 MFD");
		charMFD.setName("NSHM Char MFD");
		combinedMFD.setName("NSHM Combined MFD");
		
		gr_AveSectPartMFD.setName("gr_AveSectPartMFD");
		char_AveSectPartMFD.setName("char_AveSectPartMFD");
		combined_AveSectPartMFD.setName("combined_AveSectPartMFD");
		
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(gr_AveSectPartMFD);
		funcs.add(grMFD);
		funcs.add(grMFD.getCumRateDistWithOffset());
		funcs.add(char_AveSectPartMFD);
		funcs.add(charMFD);
		funcs.add(charMFD.getCumRateDistWithOffset());
		funcs.add(combined_AveSectPartMFD);
		funcs.add(combinedMFD);
		funcs.add(combinedMFD.getCumRateDistWithOffset());
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		
		Range xRange = new Range(combinedMFD.getMinX()-combinedMFD.getDelta()/2.0, combinedMFD.getMaxX()+combinedMFD.getDelta()/2.0);
		Range yRange = new Range(1e-4, 0.2);
		
		String[] names = {"gr","char","combined"};
		int index = 0;
		String dir = ROOT_PATH+"ClassicMFDs";
		File file = new File(dir);
		file.mkdir();
		for(String name:names) {
			ArrayList<XY_DataSet> funcs2 = new  ArrayList<XY_DataSet>();
			funcs2.add(funcs.get(index));
			funcs2.add(funcs.get(index+1));
			funcs2.add(funcs.get(index+2));
			PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars, null, "Magnitude", "Rate (per yr)", xRange, yRange, 
					false, true, 3.5, 3.0, dir+"/"+name, true);	
			index+=3;
		}
	}
	
	
	public void plotRateVsGRbValue(ScalingRelationshipEnum scalingRel) {
		initData();
		double minB = -2.0;
		double maxB = 4.0;
		int numB = (int)((maxB-minB)*10.0) + 1;
		EvenlyDiscretizedFunc rateVsB_func = new EvenlyDiscretizedFunc(minB, maxB, numB);
		rateVsB_func.setName("rateVsB_func");
		double origWidthKm = faultSectionDataList.get(0).getOrigDownDipWidth();
		double minRupLengthKm = faultSectionDataList.get(0).getTraceLength()*NUM_SUBSECT_PER_RUP;
		double minRupAreaKmSq = faultSectionDataList.get(0).getReducedDownDipWidth()*minRupLengthKm;
		double minMag = scalingRel.getMag(minRupAreaKmSq*1e6, minRupLengthKm*1e3, origWidthKm*1e3);
		double minMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(minMag);
		double maxRupLengthKm = parentFaultData.getTraceLength();
		double maxRupAreaKmSq = parentFaultData.getReducedDownDipWidth()*maxRupLengthKm;
		double maxMag = scalingRel.getMag(maxRupAreaKmSq*1e6, maxRupLengthKm*1e3, origWidthKm*1e3);
		double maxMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(maxMag);
		double toMoRate = parentFaultData.calcMomentRate(true);
		
		int num = (int)Math.round((maxMFD_Mag-minMFD_Mag)/FaultSystemRuptureRateInversion.MAG_DELTA + 1);

		for(int i=0;i<rateVsB_func.size();i++) {
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, rateVsB_func.getX(i));
			rateVsB_func.set(i,gr.getTotalIncrRate());
		}
		
		System.out.println(rateVsB_func);
		
		double minRate = getTargetMFD( scalingRel, MFD_TargetType.MIN_RATE).getTotalIncrRate();
		double maxRate = getTargetMFD( scalingRel, MFD_TargetType.MAX_RATE).getTotalIncrRate();
		System.out.println("minRate="+minRate);
		System.out.println("maxRate="+maxRate);
		
		double NSHM_Rate = 0.333*rateVsB_func.getY(1.0) + 0.6666*minRate*0.93; // 0.93 is for aleatory variability ignored here
		System.out.println("NSHM_Rate="+NSHM_Rate);
		double NSHM_ImpliedB = rateVsB_func.getFirstInterpolatedX(NSHM_Rate);
		System.out.println("NSHM_ImpliedB="+NSHM_ImpliedB);
		rateVsB_func.setInfo("NSHM_Rate="+NSHM_Rate+"\nNSHM_ImpliedB="+NSHM_ImpliedB);
		
		DefaultXY_DataSet minRateFunc = new DefaultXY_DataSet();
		minRateFunc.setName("minRateFunc");
		minRateFunc.set(minB,minRate);
		minRateFunc.set(maxB,minRate);
		
		DefaultXY_DataSet maxRateFunc = new DefaultXY_DataSet();
		maxRateFunc.setName("maxRateFunc");
		maxRateFunc.set(minB,maxRate);
		maxRateFunc.set(maxB,maxRate);
		
		DefaultXY_DataSet nshm_RateFunc = new DefaultXY_DataSet();
		nshm_RateFunc.setName("nshm_RateFunc");
		nshm_RateFunc.set(minB,NSHM_Rate);
		nshm_RateFunc.set(maxB,NSHM_Rate);
		
		
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		funcs.add(rateVsB_func);
		funcs.add(minRateFunc);
		funcs.add(maxRateFunc);
		funcs.add(nshm_RateFunc);
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		
		String dir = ROOT_PATH+"/rateVsGR_bValuePlots";
	    File file = new File(dir);
	    file.mkdirs();			

		
		String fileNamePrefix = dir+"/rateVsGR_bValuePlot_"+scalingRel.getShortName();
		String plotName = scalingRel.getShortName();
		String xAxisLabel = "GR b-value";
		String yAxisLabel = "Total Rate (per yr)";
		Range xAxisRange = new Range(minB,maxB);
		Range yAxisRange = new Range(1e-3,5);
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, 3.5, 3.0, fileNamePrefix, true);

	}
	
	private static void corrMoRateForScalingRelationship(ScalingRelationshipEnum scalingRel, IncrementalMagFreqDist mfd,
			double downDipWidthMeters) {
			
		double moRate = mfd.getTotalMomentRate();
		double logAreaMinMeters = Math.log10(10*1e6);	// 10 km-sq
		double logAreaMaxMeters = Math.log10(50000*1e6);//50,000 km-sq
		EvenlyDiscretizedFunc magAreaFunc = new EvenlyDiscretizedFunc(logAreaMinMeters, logAreaMaxMeters, 200);
		for(int i=0;i<magAreaFunc.size();i++) {
			double areaSqM = Math.pow(10,magAreaFunc.getX(i));
			double length = areaSqM/downDipWidthMeters;
			magAreaFunc.set(i, scalingRel.getMag(areaSqM, length, downDipWidthMeters));
		}
		double moRate2 = 0;
		for(int i=0;i<mfd.size();i++) {
			double mag = mfd.getX(i);
			double area = Math.pow(10, magAreaFunc.getFirstInterpolatedX(mag));
			double slip = scalingRel.getAveSlip(area, area/downDipWidthMeters, downDipWidthMeters);
			double moment = FaultMomentCalc.getMoment(area, slip);
			moRate2+=moment*mfd.getY(i);
		}
		mfd.scale(moRate/moRate2);
//		System.out.println(scalingRel+" correction: "+(float)(moRate2/moRate));
//		System.exit(-1);	
	}
	
	
	public IncrementalMagFreqDist getTargetMFD(ScalingRelationshipEnum scalingRel, MFD_TargetType mfdType) {

		double origWidthKm = faultSectionDataList.get(0).getOrigDownDipWidth();
		
		double minRupLengthKm = faultSectionDataList.get(0).getTraceLength()*NUM_SUBSECT_PER_RUP;
		double minRupAreaKmSq = faultSectionDataList.get(0).getReducedDownDipWidth()*minRupLengthKm;
		double minMag = scalingRel.getMag(minRupAreaKmSq*1e6, minRupLengthKm*1e3, origWidthKm*1e3);
		double minMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(minMag);
		
		double maxRupLengthKm = parentFaultData.getTraceLength();
		double maxRupAreaKmSq = parentFaultData.getReducedDownDipWidth()*maxRupLengthKm;
		double maxMag = scalingRel.getMag(maxRupAreaKmSq*1e6, maxRupLengthKm*1e3, origWidthKm*1e3);
		double maxMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(maxMag);
		
		double toMoRate = parentFaultData.calcMomentRate(true);
		
		int num = (int)Math.round((maxMFD_Mag-minMFD_Mag)/FaultSystemRuptureRateInversion.MAG_DELTA + 1);
		IncrementalMagFreqDist mfd=null;
		switch (mfdType) {
		case GR_b_0pt32:
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.32);
			gr.setName("GR, b=0.32");
			mfd = gr;
			break;
		case GR_b_1pt6:
			GutenbergRichterMagFreqDist gr_b1p6 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b1p6.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.6);
			gr_b1p6.setName("GR, b=1.6");
			mfd = gr_b1p6;
			break;
		case GR_b_1pt4:
			GutenbergRichterMagFreqDist gr_b1p4 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b1p4.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.4);
			gr_b1p4.setName("GR, b=1.4");
			mfd = gr_b1p4;
			break;
		case GR_b_1pt2:
			GutenbergRichterMagFreqDist gr_b1p2 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b1p2.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.2);
			gr_b1p2.setName("GR, b=1.2");
			mfd = gr_b1p2;
			break;
		case GR_b_1pt0:
			GutenbergRichterMagFreqDist gr_b1 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b1.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.0);
			gr_b1.setName("GR, b=1.0");
			mfd = gr_b1;
			break;
		case GR_b_0pt8:
			GutenbergRichterMagFreqDist gr_b0p8 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0p8.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.8);
			gr_b0p8.setName("GR, b=0.8");
			mfd = gr_b0p8;
			break;
		case GR_b_0pt8_ForSegmented:
			GutenbergRichterMagFreqDist tempGR = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			tempGR.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.8);
			mfd = new IncrementalMagFreqDist(minMFD_Mag,num-5,0.1);	// skip M>7.75
			for(int i=0;i<mfd.size();i++)
				if(i != 1) // skip the empty cell
					mfd.set(i, tempGR.getY(i));
			mfd.scaleToTotalMomentRate(toMoRate);
			break;
		case GR_b_0pt6:
			GutenbergRichterMagFreqDist gr_b0p6 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0p6.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.6);
			gr_b0p6.setName("GR, b=0.6");
			mfd = gr_b0p6;
			break;
		case GR_b_0pt4:
			GutenbergRichterMagFreqDist gr_b0p4 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0p4.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.4);
			gr_b0p4.setName("GR, b=0.4");
			mfd = gr_b0p4;
			break;
		case GR_b_0pt2:
			GutenbergRichterMagFreqDist gr_b0p2 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0p2.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.2);
			gr_b0p2.setName("GR, b=0.2");
			mfd = gr_b0p2;
			break;
		case GR_b_0pt0:
			GutenbergRichterMagFreqDist gr_b0 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.0);
			gr_b0.setName("GR, b=0.0");
			mfd = gr_b0;
			break;
		case GR_b_minus1:
			GutenbergRichterMagFreqDist gr_bMinus1 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_bMinus1.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, -1.0);
			gr_bMinus1.setName("GR, b=1");
			mfd = gr_bMinus1;
			break;
		case MAX_RATE:
			IncrementalMagFreqDist maxMFD = new IncrementalMagFreqDist(minMFD_Mag,num,0.1);
			maxMFD.set(0,1.0);
			// this correct for the difference between orig and rounded mag
			double magRoundingCorr = MagUtils.magToMoment(minMFD_Mag)/MagUtils.magToMoment(minMag);
			maxMFD.scaleToTotalMomentRate(toMoRate*magRoundingCorr);
			mfd = maxMFD;
			break;
		case MIN_RATE:
			IncrementalMagFreqDist minMFD = new IncrementalMagFreqDist(minMFD_Mag,num,0.1);
			minMFD.set(maxMFD_Mag,1.0);
			// this correct for the difference between orig and rounded mag
			double magRoundingCorr2 = MagUtils.magToMoment(maxMFD_Mag)/MagUtils.magToMoment(maxMag);
			minMFD.scaleToTotalMomentRate(toMoRate*magRoundingCorr2);
			mfd = minMFD;
			break;
		case M7pt25only:
			mfd = new IncrementalMagFreqDist(minMFD_Mag,num,0.1);
			mfd.set(7.25,1.0);
			// this correct for the difference between orig and rounded mag
			mfd.scaleToTotalMomentRate(toMoRate);
			break;			
		case NSHM_Combined:	// this ignores aleatory variability
			IncrementalMagFreqDist mfd1 = new IncrementalMagFreqDist(minMFD_Mag,num,0.1);
			mfd1.set(maxMFD_Mag,1.0);
			// this correct for the difference between orig and rounded mag
			double magRoundingCorr3 = MagUtils.magToMoment(maxMFD_Mag)/MagUtils.magToMoment(maxMag);
			mfd1.scaleToTotalMomentRate(toMoRate*magRoundingCorr3);
			GutenbergRichterMagFreqDist mfd2 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			mfd2.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.0);
			mfd2.setName("GR, b=1");
			mfd = new IncrementalMagFreqDist(minMFD_Mag,num,0.1);
			for(int i=0;i<mfd.size();i++) {
				mfd.set(i, 0.333*mfd2.getY(i)+0.667*mfd1.getY(i));
			}
			mfd.setName("NSHM_Combined");
			break;
		case NONE:
			mfd = null;
			break;
		}
		
		if(scalingRel == ScalingRelationshipEnum.ELLB_SQRT_LENGTH || scalingRel == ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP)
			corrMoRateForScalingRelationship(scalingRel, mfd, (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*1e3);

		return mfd;
	}
	
	
	/**
	 * This sets/resets default parameter values
	 */
	public void setDefaultParameterValuess() {
		 dirName = ROOT_PATH+"NoNameOutput";
		 solutionName = "No Name Solution"; // Inversion name
//		 shortFault = true;	// change this elsewhere
		 wtedInversion = true; // Whether to apply data weights
		 slipRateProfile = SlipRateProfileType.TAPERED;
		 slipModelType = SlipAlongRuptureModelEnum.TAPERED; // Average slip along rupture model
		 scalingRel = ScalingRelationshipEnum.ELLSWORTH_B; // Scaling Relationship
		 slipRateSegmentationConstraintList = new ArrayList<SlipRateSegmentationConstraint> ();
		 sectionRateConstraintList = new ArrayList<SectionRateConstraint>();
		 relativeSectRateWt=0; // Section Rate Constraints wt
		 segmentationConstrList = new ArrayList<SegmentationConstraint>(); // ;
		 relative_segmentationConstrWt = 0; // Segmentation Constraints wt
		 relative_aPrioriRupWt = 0;  // A priori rupture rates wts
		 aPrioriRupRateFilename = null; // file that has the a prior rupture rates
		 setAprioriRatesFromMFD_Constraint = false;
		 minRupRate = 0.0; // Minimum Rupture Rate; not applied for NSHMP solution
		 applyProbVisible = false; // Apply prob visible
		 moRateReduction = 0.0;	// reduction due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)		
		 mfdTargetType = MFD_TargetType.GR_b_1pt0; // Target MFD Constraint
		 relativeMFD_constraintWt = 0; // Target MFD Constraint Wt
		 totalRateConstraint = 0.005328; // Total Rate Constraint
		 totalRateSigma = 0.1*totalRateConstraint;
		 relativeTotalRateConstraintWt = 0.0;
		 smoothnessConstraintList = null;	// Section rate smoothness constraints
		 relativeSmoothnessConstraintWt = 0.0;		// Section rate smoothness weight
		 applyRuptureSampler = false;
		 
		 solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		 randomSeed = 0;	// zero means use current time in millis as the seed
		 numSolutions = 1; // this is ignored for NON_NEGATIVE_LEAST_SQUARES which only has one possible solution
		 saCooling = CoolingScheduleType.FAST_SA;
		 perturbationFunc = GenerationFunctionType.UNIFORM_NO_TEMP_DEPENDENCE;
		 completionCriteria = new IterationCompletionCriteria((long) 1e5);
//completionCriteria = new EnergyCompletionCriteria(600d);	// number of rows/subsections

		 initialStateType = SA_InitialStateType.ALL_ZEROS;
//		 rupRatesFromFileDirName = null;  //if solution is from a file (previous solution); name example: ROOT_PATH+"OutDir/";
		 magAareaAleatoryVariability = 0.0;
		// Data to plot and/or save:
		 popUpPlots = true;	// this tells whether to show plots in a window (set null to turn off; e.g., for HPC)
//		 hazCurveLoc = new Location(40,-117); // to make hazard curve (set loc=null to ignore)
//	     hazCurveLocName = "Hazard Curve at "+hazCurveLoc.toString();
		 hazCurveLocList = null;
	     hazCurveLocNameList = null;

	     makeHazardMaps = false; //// to make hazard maps
//		 saPeriodForHaz = 0.0;	// set as 0.0 for PGA
	}
	
	
	/**
	 * This generates a solution and diagnostic info for the current parameter settings
	 */
	public FaultSystemRuptureRateInversion getSolution(boolean rePlotOny) {	
		if(randomSeed == 0)
			randomSeed = System.currentTimeMillis();
		initData();
		
		IncrementalMagFreqDist targetMFD = getTargetMFD(scalingRel, mfdTargetType);
		IncrementalMagFreqDist mfdSigma = null;
		if(targetMFD != null && targetMFD.getMinY()>0) {
			mfdSigma = getTargetMFD(scalingRel, mfdTargetType);
			mfdSigma.scale(0.1); // uncertainty is 10%
		}
		if(D) {
			if(targetMFD != null)
				System.out.println(mfdTargetType+" total rate "+targetMFD.getTotalIncrRate());
		}
		
		ArrayList<FaultSectionPrefData> fltSectDataList = getFaultSectionDataList();
		int[][] rupSectionMatrix = getRupSectionMatrix();
		
		// this will be used to keep track of runtimes
		long startTimeMillis = System.currentTimeMillis();
		if(D)
			System.out.println("Starting Inversion");
		
		// create an instance of the inversion class with the above settings
		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
				solutionName,
				slipRateProfile.toString(),
				fltSectDataList, 
				rupSectionMatrix, 
				slipModelType, 
				scalingRel, 
				slipRateSegmentationConstraintList,
				sectionRateConstraintList, 
				relativeSectRateWt, 
				relative_aPrioriRupWt, 
				aPrioriRupRateFilename,
				wtedInversion, 
				minRupRate, 
				applyProbVisible, 
				moRateReduction,
				targetMFD,
				mfdSigma,
				relativeMFD_constraintWt,
				segmentationConstrList,
				relative_segmentationConstrWt,
				totalRateConstraint,
				totalRateSigma,
				relativeTotalRateConstraintWt,
				smoothnessConstraintList,
				relativeSmoothnessConstraintWt,
				magAareaAleatoryVariability);

		
		// make the directory for storing results
		if(!rePlotOny) {
		    File file = new File(dirName);
		    file.mkdirs();			
		}
	    
	    // set a-prior rates from MFD so these can be applied as initial model
	    if(setAprioriRatesFromMFD_Constraint)
	    	fltSysRupInversion.setAprioriRupRatesFromMFD_Constrint();
	    
	    // write the setup info to a file
		if(!rePlotOny)
			fltSysRupInversion.writeInversionSetUpInfoToFile(dirName);

		
		// do this before overriding solutionType for rePlotOnly case (otherwise rupSamplerMFD won't be plotted)
	    IntegerPDF_FunctionSampler rupSampler = null;
	    if(solutionType == InversionSolutionType.SIMULATED_ANNEALING) {
		    // set the rupture sampler; default is GR if targetMFD == null
		    if(applyRuptureSampler) {
		    	double[] rupSampleProbArray;
		    	if(targetMFD == null)	// GR is default
			    	rupSampleProbArray = fltSysRupInversion.getRupRatesForTargetMFD(getTargetMFD(scalingRel, MFD_TargetType.GR_b_1pt0), false);
		    	else
		    		rupSampleProbArray = fltSysRupInversion.getRupRatesForTargetMFD(targetMFD, false);
			    rupSampler = new IntegerPDF_FunctionSampler(rupSampleProbArray.length);
			    for(int r=0; r<rupSampleProbArray.length; r++)
			    	rupSampler.set(r,rupSampleProbArray[r]);	
			    // do this here so re-plot includes this MFD
// System.out.println("targetMFD:\n"+targetMFD);
			    fltSysRupInversion.computeRupSamplerMFD(rupSampler);
		    }
	    }
	    
	    
//	    // ********** TEST **************
	    // test choosing only ruptures not included in a previous run
//	    applyRuptureSampler = true;
//		double[] tempRupRatesArray = readRuptureRatesFromFile(ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B_Alt/ruptureRates.txt");
//	    rupSampler = new IntegerPDF_FunctionSampler(tempRupRatesArray.length);
//	    int numZero=0, numOne=0;
//	    for(int r=0; r<tempRupRatesArray.length; r++)
//	    	if(tempRupRatesArray[r] == 0) {
//	    		rupSampler.set(r,1.0);	
//	    		numOne += 1;
//	    	}
//	    	else {
//	    		rupSampler.set(r,0.0);	
//	    		numZero += 1;
//	    	}
//	    System.out.println("numOne = "+numOne+"; numZero = "+numZero+"; total = "+numZero+numOne);
//	    fltSysRupInversion.computeRupSamplerMFD(rupSampler);


	    
	    
	    
		if(rePlotOny) { // overide solution type
			solutionType = InversionSolutionType.FROM_FILE;
		}


		// set the initial state & rupSampler if SA
	    double[] initialState = null;
	    if(solutionType == InversionSolutionType.SIMULATED_ANNEALING) {
	    	switch(initialStateType) {
	    	case ALL_ZEROS:
	    		initialState = new double[fltSysRupInversion.getNumRuptures()];
	    		break;
	    	case A_PRIORI_RATES:
	    		double tempArray[] = fltSysRupInversion.getAprioriRuptureRates();
	    		if(tempArray.length != fltSysRupInversion.getNumRuptures())
	    			throw new RuntimeException("aPrior rates are not set for each rupture");
	    		initialState = tempArray;
	    		break;
	    	case FROM_MFD_CONSTRAINT:
	    		initialState = fltSysRupInversion.getRupRatesForTargetMFD(targetMFD, false);
	    		break;
	    	}
	    }
	    
	    
	    switch (solutionType) {
	    case NON_NEGATIVE_LEAST_SQUARES:
	    	if(D) System.out.println("NNLS Solution");
	    	fltSysRupInversion.doInversionNNLS();
	    	break;
	    case FRESH:
	    	if(D) System.out.println("FRESH Solution");
	    	fltSysRupInversion.doFRESH_Solution(targetMFD);
	    	break;
	    case SUNFISH:
	    	if(D) System.out.println("SUNFISH Solution");
	    	fltSysRupInversion.doSUNFiSH_Solution();;
	    	break;
	    case RATES_FROM_MFD:
	    	if(D) System.out.println("RATES_FROM_MFD Solution");
	    	fltSysRupInversion.doRatesFromMFD_Solution(targetMFD);
	    	break;
	    case SIMULATED_ANNEALING:
	    	if(numSolutions==1) {
	    		if(D) System.out.println("SIMULATED_ANNEALING; numSolutions=1");
	    		fltSysRupInversion.doInversionSA(completionCriteria, initialState, randomSeed, saCooling, perturbationFunc, rupSampler);
	    	}
	    	else if(numSolutions>1) {
	    		if(D) System.out.println("SIMULATED_ANNEALING; numSolutions="+numSolutions);
	    		fltSysRupInversion.doInversionSA_MultTimes(completionCriteria, initialState, randomSeed, numSolutions, dirName, saCooling,perturbationFunc, rupSampler);
	    	}
	    	else {
	    		throw new RuntimeException("bad numIterations");
	    	}
	    	break;
	    case FROM_FILE:
	    	if(numSolutions==1) {
	    		if(D) System.out.println("FROM_FILE; numSolutions=1");
	    		String rupRatesFileName = dirName+ "/ruptureRatesAlt.txt";
	    		double[] rupRatesArray = readRuptureRatesFromFile(rupRatesFileName);
	    		fltSysRupInversion.setSolution(rupRatesArray, "Solution from file: "+rupRatesFileName);
	    	}
	    	else if(numSolutions>1) {
	    		if(D) System.out.println("FROM_FILE; numSolutions="+numSolutions);
	    		ArrayList<double[]> rupRatesArrayList = new ArrayList<double[]>();
	    		for(int i=0;i<numSolutions;i++) {
	    			String rupRatesFileName = dirName+ "/ruptureRates_"+i+".txt";
	    			rupRatesArrayList.add(readRuptureRatesFromFile(rupRatesFileName));
	    			fltSysRupInversion.setMultipleSolutions(rupRatesArrayList, "Multiple solutions read from files with prefix "+rupRatesFileName, dirName);
	    		}
	    	}
	    	else {
	    		throw new RuntimeException("bad numSolutions");
	    	}
	    	break;
	    case APPLY_SOLUTION:
	    	if(numSolutions==1) {
	    		if(D) System.out.println("APPLY_SOLUTION used");
	    		fltSysRupInversion.setSolution(solutionRatesToApplyArray, "Solution rates applied");
	    	}
	    	else {
	    		throw new RuntimeException("bad numSolutions (can only apply one here)");
	    	}
	    	break;
	    }

		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
		
		String runTimeString = "Done with Inversion after "+(float)runTimeSec+" seconds.";
		if(D) System.out.println(runTimeString);
		
		// Write out info if not replotting
		if(!rePlotOny) { 
			fltSysRupInversion.addToModelRunInfoString("\n"+runTimeString+"\n");
			// write results to file
			fltSysRupInversion.writeInversionRunInfoToFile(dirName);
			fltSysRupInversion.writeRuptureRatesToFile(dirName);
		}
		
		
		// plot various model related histograms
		fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots, 3.5, 3.0);

		// plot solution MFDs
		fltSysRupInversion.writeAndOrPlotMFDs(dirName, popUpPlots, null, new Range(1e-5,1), 3.5, 3.0);
		
		// plot various things as a function of section index
		makeResultsVsSectionPlots(fltSysRupInversion);
		
		// plot normalized residual versus row index
		fltSysRupInversion.writeAndOrPlotNormalizedResiduals(dirName, popUpPlots, 3.5, 3.0);
		
		// plot section boundary rates
		fltSysRupInversion.writeAndOrPlotSectionBoundaryRates(dirName, popUpPlots, null, 6.5, 1.7);
		
		// plot rupture rate versus index (big and not very informative plots)
//			fltSysRupInversion.writeAndOrPlotRupRateVsIndex(dirName, popUpPlots, 6.5, 4.0);
		
		// plot rupture event and slip rate as a function of section for all non-zero ruptures (big and not very informative, except to confirm that tapered slip working)
//			fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots, 9.0, 6.5);

		// plot the rate of each rupture (triangle plot with rate of first section)
		fltSysRupInversion.writeAndOrPlotRupRateOfFirstSection(dirName, popUpPlots, 3.5, 4.0);

		// Plot MFDs for any sections that have constrints
		for(SectionRateConstraint sectRateConstr : sectionRateConstraintList)
			fltSysRupInversion.writeAndOrPlotPartMFD_ForSection(dirName, popUpPlots, sectRateConstr.getSectIndex());
		for(SlipRateSegmentationConstraint srSegConstr : slipRateSegmentationConstraintList)
			fltSysRupInversion.writeAndOrPlotPartMFD_ForSection(dirName, popUpPlots, srSegConstr.getSectIndex());
		for(SegmentationConstraint segConst :segmentationConstrList)
			fltSysRupInversion.writeAndOrPlotJointPartMFD_ForSections(dirName, popUpPlots, segConst.getSect1_Index(), segConst.getSect2_Index());
		
		// Hazard calculations
		double[] saPeriodForHazArray = {0.0, 1.0};
		
		for(double saPeriodForHaz : saPeriodForHazArray) {

			// hazard curve:
			if(hazCurveLocList != null) {
				for(int k=0;k<hazCurveLocList.size();k++)
					writeAndOrPlotHazardCurve(fltSysRupInversion, hazCurveLocList.get(k), saPeriodForHaz, dirName, popUpPlots, hazCurveLocNameList.get(k));
			}

			// second parameter here is SA period; set as 0.0 for PGA:
			if(makeHazardMaps) {
				String hazDirName = null;
				if(dirName != null) {
					hazDirName = dirName+"/hazardMaps";
				}

				// make mean map
				makeHazardMaps(fltSysRupInversion, saPeriodForHaz, hazDirName, popUpPlots);

				if(makeHazardMapsForAllSolutions) {
					if(numSolutions>1) {
						String imtString = "PGA";
						if(saPeriodForHaz != 0)
							imtString = saPeriodForHaz+"secSA";

						ArrayList<GriddedGeoDataSet[]> gridHazDataArray = makeHazardMapGriddedDataForMultRuns(fltSysRupInversion, saPeriodForHaz, hazDirName, true);
						// 2in50
						int probIndex=0;
						String label = imtString+"_"+hazardProbNameArray[probIndex];
						plotNormPDF_FromMultHazMaps(gridHazDataArray.get(probIndex), label, hazDirName, popUpPlots); 
						// 10in50
						probIndex=1;
						label = imtString+"_"+hazardProbNameArray[probIndex];
						plotNormPDF_FromMultHazMaps(gridHazDataArray.get(probIndex), label, hazDirName, popUpPlots); 
						// RTGM
						if(saPeriodForHaz == 1.0) {
							probIndex=2;
							label = imtString+"_RTGM";
							plotNormPDF_FromMultHazMaps(gridHazDataArray.get(probIndex), label, hazDirName, popUpPlots); 
						}
					}
				}

			}
		}
		if(D) System.out.println("Done getting solution");

		return fltSysRupInversion;
		
	}

	
	/**
	 * This applies b=1.0 since aftershocks are included (and to be more general than USGS NSHM)
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 * @param makeHazardMaps
	 */
	public void doNSHMP_GR_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel, 
			boolean makeHazardMaps, LocationList hazCurveLocList, ArrayList<String> hazCurveLocNameList) {
		this.setDefaultParameterValuess();
		this.hazCurveLocList = hazCurveLocList;
		this.hazCurveLocNameList = hazCurveLocNameList;
		this.makeHazardMaps=makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = MFD_TargetType.GR_b_1pt0; // Target MFD Constraint
		solutionName = "NSHMP_GR_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); 
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}
	
	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doNSHMP_Char_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, boolean makeHazardMaps, LocationList hazCurveLocList, ArrayList<String> hazCurveLocNameList) {
		this.setDefaultParameterValuess();
		// do this to get the number of ruptures:
		initData();
		this.hazCurveLocList = hazCurveLocList;
		this.hazCurveLocNameList = hazCurveLocNameList;
		this.makeHazardMaps=makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		magAareaAleatoryVariability = 0.12;	// gaussian distribution sigma
		// make the solution
		solutionType = InversionSolutionType.APPLY_SOLUTION;
		solutionRatesToApplyArray = new double[numRuptures];
		double maxRupAreaKmSq = (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*FAULT_LENGTH_KM;
		double mag = scalingRel.getMag(maxRupAreaKmSq*1e6, FAULT_LENGTH_KM*1e3, (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*1e3);
		double maxMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(mag);
		double toMoRate = FaultMomentCalc.getMoment(maxRupAreaKmSq*1e6, FAULT_SLIP_RATE*1e-3);
		double magRoundingCorr = MagUtils.magToMoment(maxMFD_Mag)/MagUtils.magToMoment(mag);
		solutionRatesToApplyArray[numRuptures-1] = (toMoRate*magRoundingCorr/1.071748)/MagUtils.magToMoment(maxMFD_Mag); // 1.071748 is for aleatory variability in magnitude
		solutionName = "NSHMP_Char_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}
	
	
	// this ignores the aleatory variability on the char part, which has less than a 3% influence on 2in50 values for the pure char haz maps
	// I need to build the MFD etc plots by combining the others
	public void doNSHMP_CombinedSolution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		// do this to get the number of ruptures:
		initData();
		this.makeHazardMaps = makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.APPLY_SOLUTION;
		String gr_solutionName = "NSHMP_GR_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); 
		String gr_fileName = ROOT_PATH+gr_solutionName+"/ruptureRatesAlt.txt";
		double[] grRatesArray = readRuptureRatesFromFile(gr_fileName);
		
		String char_solutionName = "NSHMP_Char_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); 
		String char_fileName = ROOT_PATH+char_solutionName+"/ruptureRatesAlt.txt";
		double[] charRatesArray = readRuptureRatesFromFile(char_fileName);
		
		if(grRatesArray.length != numRuptures || charRatesArray.length != numRuptures)
			throw new RuntimeException("Array length problem (should be 6555): "+grRatesArray.length+", "+charRatesArray.length);
		solutionRatesToApplyArray = new double[6555];
		for(int i=0;i<solutionRatesToApplyArray.length;i++)
			solutionRatesToApplyArray[i] = 0.6666*1.071748*charRatesArray[i] + 0.3333*grRatesArray[i]; // 1.071748 is to correct for no aleatory variability here
		solutionName = "NSHMP_Combined_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}
	
	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doMinRateSolution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, boolean makeHazardMaps, LocationList hazCurveLocList, ArrayList<String> hazCurveLocNameList) {
		this.setDefaultParameterValuess();
		// do this to get the number of ruptures:
		initData();
		this.makeHazardMaps = makeHazardMaps;
		this.hazCurveLocList = hazCurveLocList;
		this.hazCurveLocNameList = hazCurveLocNameList;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		// make the solution
		solutionType = InversionSolutionType.APPLY_SOLUTION;
		solutionRatesToApplyArray = new double[numRuptures];
		double maxRupAreaKmSq = (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*FAULT_LENGTH_KM;
		double mag = scalingRel.getMag(maxRupAreaKmSq*1e6, FAULT_LENGTH_KM*1e3, (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*1e3);
		double toMoRate = FaultMomentCalc.getMoment(maxRupAreaKmSq*1e6, FAULT_SLIP_RATE*1e-3);
//		makeFaultSectionDataList();
//		FaultSectionPrefData sectData = faultSectionDataList.get(0);
//		double maxRulLengthKm = faultSectionDataList.size()*sectData.getTraceLength();	
//		double maxRupAreaKmSq = maxRulLengthKm*sectData.getReducedDownDipWidth();
//		double mag = scalingRel.getMag(maxRupAreaKmSq*1e6, maxRulLengthKm*1e3, sectData.getReducedDownDipWidth()*1e3);
//		double toMoRate = faultSectionDataList.size()*sectData.calcMomentRate(true);	
		
		// round the mag
		double maxMFD_Mag = ((double)Math.round(100*mag))/100.0; // FaultSystemRuptureRateInversion.roundMagTo10thUnit(mag);
		double magRoundingCorr = MagUtils.magToMoment(maxMFD_Mag)/MagUtils.magToMoment(mag);
		solutionRatesToApplyArray[numRuptures-1] = toMoRate*magRoundingCorr/MagUtils.magToMoment(maxMFD_Mag); 
		solutionName = "MinRateSol_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}

	
	/**
	 * The remaining slip-rate discrepancy is from floater end effects, as only one rupture hits the first and last section, whereas 
	 * most sections participate in 4 different ruptures; the average slip-rate discrepancy is zero.
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doMaxRateSolution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel,
			boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = MFD_TargetType.MAX_RATE; // Target MFD Constraint
		solutionName = "MaxRateSol_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}
	


	public void doUnconstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, boolean initState, boolean applyRupSampler) {
				
		this.setDefaultParameterValuess();
		 applyRuptureSampler = applyRupSampler;

//		makeHazardMaps = true;

		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		if(initState) {
			this.mfdTargetType = MFD_TargetType.GR_b_1pt0;
			relativeMFD_constraintWt = 0.0; // Target MFD Constraint Wt			
		}
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		String initStateString = "";
		if(initState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		numSolutions = numSim; 
		solutionName = "UnconstrSA_"+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+"_"+
		slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}
	
	
	public void doSegConstrainedSA_old(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		segmentationConstrList.add(new SegmentationConstraint("Test Seg Const", 58, 59, 5e-4, 5e-6));
		this.relative_segmentationConstrWt = 1.0;
		String namePrefix = "DoubleFltSegConstrSA_";
		String initStateString = "";
		if(mfdTargetType != null) {
			namePrefix += "wMFD_";
			this.mfdTargetType = mfdTargetType;
			relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
			if(asInitState) {
				initStateString = "_initSol";
				initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
			}		
		}
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = namePrefix+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+
						"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	}
	
	public void doSegConstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, 
			double mfdWt, boolean asInitState, ArrayList<SegmentationConstraint> segConstrList, String sectInfo, boolean applyRupSampler, 
			boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		this.segmentationConstrList = segConstrList;
		this.relative_segmentationConstrWt = 1.0;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		String initStateString = "";
		if(asInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		String prefix = "";
		if(!shortFault)
			prefix = "LongFlt_";

		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = prefix+"SegConstrSA_"+numSolutions+"_finalE="+finalEnergy+"_"+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+
				slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString+"_"+sectInfo; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}


	
	
	public void OLDdoSectRateConstrSegSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState,
			ArrayList<SectionRateConstraint> segConstrList, String segInfo, boolean applyRupSampler, boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.relativeSectRateWt = 1;
	//	SectionRateConstraint sectRateConstr = new SectionRateConstraint("Test Sect Constr", 58, 5e-4, 5e-6);
		sectionRateConstraintList = segConstrList;
		this.relative_segmentationConstrWt = 1.0;
		
		
		String namePrefix = "DoubleSectRateConstrSegSA_";
		String initStateString = "";
		if(mfdTargetType != null) {
			namePrefix += "wMFD_";
			this.mfdTargetType = mfdTargetType;
			relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
			if(asInitState) {
				initStateString = "_initSol";
				initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
			}		
		}
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		
		ArrayList<CompletionCriteria> criteriaList = new ArrayList<CompletionCriteria>();
		criteriaList.add(new EnergyCompletionCriteria(finalEnergy));
		long time = 1000*60*10;	// 10 minutes
		criteriaList.add(new TimeCompletionCriteria(time));
		completionCriteria = new CompoundCompletionCriteria(criteriaList);
		
//		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = namePrefix+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	}
	
	public void old_doSectRateConstrSegTwoPtsSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.relativeSectRateWt = 1;
		sectionRateConstraintList.add(new SectionRateConstraint("Test Sect Constr", 58, 5e-4, 5e-6));
		sectionRateConstraintList.add(new SectionRateConstraint("Test Sect Constr", 87, 0.03, 3e-4));
		this.relative_segmentationConstrWt = 1.0;
		String namePrefix = "DoubleSectRateConstrSegTwoPtsSA_";
		String initStateString = "";
		if(mfdTargetType != null) {
			namePrefix += "wMFD_";
			this.mfdTargetType = mfdTargetType;
			relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
			if(asInitState) {
				initStateString = "_initSol";
				initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
			}		
		}
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = namePrefix+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+"_"+
					slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	}
	
	
	public void old_doSectRateConstrSegTwoPtsSmoothSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.relativeSectRateWt = 1;
		sectionRateConstraintList.add(new SectionRateConstraint("Test Sect Constr", 58, 5e-4, 5e-6));
		sectionRateConstraintList.add(new SectionRateConstraint("Test Sect Constr", 87, 0.03, 3e-4));
		this.relative_segmentationConstrWt = 1.0;
		smoothnessConstraintList = new ArrayList<int[]>();
		int[] sectArray = new int[117-2];
		for(int i=0; i<sectArray.length;i++)
			sectArray[i] = i+1;
		smoothnessConstraintList.add(sectArray);
		relativeSmoothnessConstraintWt = 1e3;
		String namePrefix = "DoubleSectRateConstrSegTwoPtsSmoothSA_";
		String initStateString = "";
		if(mfdTargetType != null) {
			namePrefix += "wMFD_";
			this.mfdTargetType = mfdTargetType;
			relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
			if(asInitState) {
				initStateString = "_initSol";
				initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
			}		
		}
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = namePrefix+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+"_"+
					slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	}


	
	public void doSlipRateSegmentedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, 
			double mfdWt, boolean asInitState, ArrayList<SlipRateSegmentationConstraint> segConstrList, String segInfo, boolean applyRupSampler, boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		slipRateSegmentationConstraintList = segConstrList;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		String initStateString = "";
		if(asInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		String prefix = "";
		if(!shortFault)
			prefix = "LongFlt_";

		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = prefix+"SlipRateSegSA_"+numSolutions+"_finalE="+finalEnergy+"_"+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+
				slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString+"_"+segInfo; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}
	
	
	public void doSectRateConstrSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, 
			double mfdWt, boolean asInitState, ArrayList<SectionRateConstraint> sectConstrList, String sectInfo, boolean applyRupSampler, 
			boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		this.sectionRateConstraintList = sectConstrList;
		this.relativeSectRateWt = 1.0;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		String initStateString = "";
		if(asInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		String prefix = "";
		if(!shortFault)
			prefix = "LongFlt_";

		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = prefix+"SectRateConstrSA_"+numSolutions+"_finalE="+finalEnergy+"_"+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+
				slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString+"_"+sectInfo; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}
	
	
	public void doSectRateConstrSmoothedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy, MFD_TargetType mfdTargetType, 
			double mfdWt, boolean asInitState, ArrayList<SectionRateConstraint> sectConstrList, String sectInfo, boolean applyRupSampler, 
			boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		this.sectionRateConstraintList = sectConstrList;
		this.relativeSectRateWt = 1.0;
		
		smoothnessConstraintList = new ArrayList<int[]>();
		int[] sectArray = new int[117-2];
		for(int i=0; i<sectArray.length;i++)
			sectArray[i] = i+1;
		smoothnessConstraintList.add(sectArray);
		relativeSmoothnessConstraintWt = 1e3;

		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		String initStateString = "";
		if(asInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		String prefix = "";
		if(!shortFault)
			prefix = "LongFlt_";

		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = prefix+"SectRateConstrSmoothSA_"+numSolutions+"_finalE="+finalEnergy+"_"+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+
				slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString+"_"+sectInfo; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}




	
	
	public FaultSystemRuptureRateInversion doMFDconstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState, 
			double finalEnergy, boolean applyRupSampler, boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
//		completionCriteria = new IterationCompletionCriteria((long) 1e4);
//		completionCriteria = new EnergyCompletionCriteria(5);
		ArrayList<CompletionCriteria> completionCriteriaList = new ArrayList<CompletionCriteria>();
		completionCriteriaList.add(new IterationCompletionCriteria((long) 1e7));	// this is just a backup
		completionCriteriaList.add(new EnergyCompletionCriteria(finalEnergy));
		completionCriteria = new CompoundCompletionCriteria(completionCriteriaList);

		String initStateString = "";
		if(asInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		numSolutions = numSim; 
		solutionName = "MFDconstrSA_"+numSolutions+"_finalE="+finalEnergy+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+
					slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		return getSolution(rePlotOny);
	
	}
	
	/**
	 * The total rate is obtained here from the supplied mfdTargetType, but the latter is not used as a constraint in the inversion
	 * @param rePlotOny
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 * @param numSim
	 * @param mfdTargetType
	 * @param mfdWt
	 * @param asInitState
	 * @param finalEnergy
	 */
	public void doTotalRateconstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, MFD_TargetType mfdTargetType, double totRateUncertFract, double totRateRelWt, 
			boolean mfdAsInitState, double finalEnergy, boolean applyRupSampler, boolean makeHazardMaps) {
		this.setDefaultParameterValuess();
		this.makeHazardMaps = makeHazardMaps;
		this.applyRuptureSampler = applyRupSampler;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = 0; // Target MFD Constraint Wt
		// I need to do this to make the MFD
		initData();
		totalRateConstraint = getTargetMFD(scalingRel, mfdTargetType).getTotalIncrRate(); // Total Rate Constraint
		if(D) System.out.println("\ntotalRateConstraint="+totalRateConstraint+"\n");
		totalRateSigma = totRateUncertFract*totalRateConstraint;
		relativeTotalRateConstraintWt = totRateRelWt;

		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
//		completionCriteria = new IterationCompletionCriteria((long) 1e4);
//		completionCriteria = new EnergyCompletionCriteria(5);
		ArrayList<CompletionCriteria> completionCriteriaList = new ArrayList<CompletionCriteria>();
		completionCriteriaList.add(new IterationCompletionCriteria((long) 1e7));	// this is just a backup
		completionCriteriaList.add(new EnergyCompletionCriteria(finalEnergy));
		completionCriteria = new CompoundCompletionCriteria(completionCriteriaList);

		String initStateString = "";
		if(mfdAsInitState) {
			initStateString = "_initSol";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		numSolutions = numSim; 
		solutionName = "totRateConstrSA_"+numSolutions+"_finalE="+finalEnergy+"_"+mfdTargetType+"_wt"+Math.round(totRateRelWt)+"_"+slipRateProfile.toString()+"_"+
				slipModelType.toString()+"_"+scalingRel.toString()+initStateString; // Inversion name
		if(applyRuptureSampler)
			solutionName+= "_rupSamp";
		dirName = ROOT_PATH+solutionName;
		getSolution(rePlotOny);
	
	}
	

	/**
	 * This applies the FRESH solution
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 * @param mfdTargetType
	 */
	public void doFRESH_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, MFD_TargetType mfdTargetType, boolean makeHazardMaps, 
			LocationList hazCurveLocList, ArrayList<String> hazCurveLocNameList) {
		this.setDefaultParameterValuess();
		this.hazCurveLocList = hazCurveLocList;
		this.hazCurveLocNameList = hazCurveLocNameList;
		this.makeHazardMaps = makeHazardMaps;
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType; // Target MFD Constraint
		solutionType = InversionSolutionType.FRESH;
		solutionName = "FRESH_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString()+"_"+mfdTargetType.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}
	
	

	/**
	 * This applies the SUNFiSH solution
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 * @param mfdTargetType
	 */
	public void doSUNFiSH_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.SUNFISH;
		solutionName = "SUNFiSH_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}

	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doAppliedMFD_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel,
			MFD_TargetType mfdType) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = mfdType; // Target MFD Constraint
		solutionName = "AppliedMFD_"+mfdType+"_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	}


	/**
	 * give zero wt to MFD constraint to ignore it
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 * @param mfdTargetType
	 * @param mfdWt
	 */
	public void doMFDconstrainedNNLS(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, MFD_TargetType mfdTargetType, double mfdWt) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType;
		relativeMFD_constraintWt = mfdWt; // Target MFD Constraint Wt
		// make the solution
		solutionType = InversionSolutionType.NON_NEGATIVE_LEAST_SQUARES;
		solutionName = "MFDconstrNNLS_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+slipRateProfile.toString()+"_"+
				slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+solutionName;
		getSolution(false);
	
	}


	public static void doLaplacianSmoothingTest() {
		int cols=9;	// num sections
		int rows = cols+1;	// num constraints

		double[][] X = new double[rows][cols];	// inversion matrices
		double[] d = new double[rows];  // the data vector
		d[rows-1] = 1;
		
		for(int r=0;r<cols-2;r++) {
			X[r][r] = -1;
			X[r][r+1] = 2;
			X[r][r+2] = -1;
		}
		X[rows-3][0] = 1;		// set rate of first section to zero
		X[rows-2][cols-1] = 1;	// set rate of last section to zero
		X[rows-1][4] = 1;
		
		for(int r=0;r<rows;r++) {
			System.out.print("\n");
			for(int c=0;c<cols;c++) {
				System.out.print(X[r][c]+"\t");
			}
		}
		
		System.out.print("\n\n");
		for(int r=0;r<rows;r++) {
			System.out.println(d[r]+"\t");
		}
		
//		System.out.println("nRow = "+X.length);
//		System.out.println( "nCol = "+X[0].length);


		double[] result = FaultSystemRuptureRateInversion.getNNLS_solution(X, d);
		System.out.print("\n\n");
		for(int c=0;c<cols;c++) {
			System.out.println(result[c]+"\t");
		}
		
		// compute predicted data
		double[] d_pred = new double[rows];  // predicted data vector
		for(int r=0;r<rows; r++)
			for(int c=0; c <cols; c++)
				d_pred[r] += result[c]*X[r][c];

		System.out.print("\n\n");
		for(int r=0;r<rows;r++) {
			System.out.println(d_pred[r]+"\t");
		}

		
	}
	
	public void computeHazardForGMM() {
		
		// get the reference model (don't recompute and don't compute hazard)
		FaultSystemRuptureRateInversion solution = doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, false);

		ArrayList<ScalarIMR> imrList = new ArrayList<ScalarIMR>();
		imrList.add(AttenRelRef.BSSA_2014.instance(null));
		imrList.add(AttenRelRef.CB_2014.instance(null));
		imrList.add(AttenRelRef.CY_2014.instance(null));
		imrList.add(AttenRelRef.ASK_2014.instance(null));
		
		double[] saPeriodForHazArray = {0.0, 1.0};

		dirName = ROOT_PATH+"/GMM_OptionsHazardVariability";
		File file = new File(dirName);
		file.mkdirs();
		
		for(ScalarIMR imrToApply:imrList) {
			this.imr = imrToApply;
			System.out.println("Working on GMM "+imr.getShortName());
			for(double saPeriodForHaz : saPeriodForHazArray) {
				String subDirName = dirName+"/"+imr.getShortName()+"_HazardMaps";
				makeHazardMaps(solution, saPeriodForHaz, subDirName, popUpPlots);
			}			
		}
	}
	
	public void makeMinMaxRateModelsHazardCurveExamplesPlot() {
		
//		// these are locations for ratios of FRESH max-rate model to min-rate model (tapered slip and slip rate); extremes are for SA 1.0 2in50 (not PGA or 10in50)
		LocationList locList2 = new LocationList();
		ArrayList<String> locNameList2 = new ArrayList<String>();
		locList2.add(new Location(31.975, -118.0));
		locNameList2.add("minRatioLoc");
		locList2.add(new Location(34.825, -117.0));
		locNameList2.add("maxRatioLoc");
		locList2.add(new Location(36.025, -117.2));
		locNameList2.add("unityRatioLoc");

		// this is to put all 1-sec SA hazard curves for the Min/Max ratios of interest on one plot
		shortFault = true;
		doMinRateSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, locList2, locNameList2);
		doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE, false, locList2, locNameList2);
		
		ArrayList<XY_DataSet> curvesToPlot = new ArrayList<XY_DataSet>();
		curvesToPlot.add(savedHazardCurves.get(9)); // Max rate solutions first (solid)
		curvesToPlot.add(savedHazardCurves.get(10));
		curvesToPlot.add(savedHazardCurves.get(11));
		curvesToPlot.add(savedHazardCurves.get(3)); // Min rate solutions (dashed)
		curvesToPlot.add(savedHazardCurves.get(4));
		curvesToPlot.add(savedHazardCurves.get(5));

		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
		
		PlottingUtils.writeAndOrPlotFuncs(curvesToPlot, plotChars, 
				"Hazard Curves", "1 sec SA", "Probability (in 50 yr)", 
				new Range(1e-2,10), new Range(1e-5,1.0), true, true, 3.5, 3.0, ROOT_PATH+"MinMaxRateModelHazardCurves", true);

	}

	
	/**
	 * 
	 * @param args
	 */
	public static void main(String []args) {
		
		// this is to test the laplacian smoothing
//		SimpleFaultInversion.doLaplacianSmoothingTest();
//		System.exit(-1);
		
		SimpleFaultInversion faultInversion = new SimpleFaultInversion();
		
		// following in reverse chronological order (so latest at the top)
		
		// these have the exact same hazard (just different implied slip rates)
//		faultInversion.doNSHMP_Char_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, true, null, null);
//		faultInversion.doNSHMP_Char_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, null, null);
//		faultInversion.doNSHMP_GR_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, false, null, null);
//		faultInversion.doNSHMP_GR_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, null, null);


//		faultInversion.plotClassicMFDs();
//		faultInversion.plotRupLengthForEllsworthB();
		
		
		// Min rate solution
//		faultInversion.doMinRateSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, locList2, locNameList2);

		// Max rate model that honors slip along fault; these two are identical (fit is not perfect, but looks good)
//		faultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE, false, locList2, locNameList2);	

//		// Hazard difference between max-rate and min-rate solutions with tapered slip rate;
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MinRateSol_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps";
//		    String label = name+"_RatioToMinRate";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		// Make hazard curve examples plots
//		faultInversion.makeMinMaxRateModelsHazardCurveExamplesPlot();

		
		// This is to plot tal rate versus b-value
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.ELLSWORTH_B);
		
		// Unconstrained model for Figure ??
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, 5, false, false);
		// show that MFD changes if we use a rupSampler
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, true);
		// show that MFD changes if we use init model
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, true, false);

		// This is the MFD constrained case; GR b=1; this took 17 minutes
//		faultInversion.makeHazardMapsForAllSolutions = true;
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);
		// test above with E=2; this took 28 minutes
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 2, false, true);

		
		// Compare to individual runs that have no overlapping ruptures
		// single run with E=5; this took 1.5 min, and 520 are non-zero
		// I changed the name of the output dir to: MFDconstrSA_1_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B_Alt
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);
		// now I re-ran the run with test code that only allows zero-rate ruptures to be selected in the following (see "********** TEST **************" here); 6.2 min; 584 above zero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);
		// I confirmed that these have none of the same ruptures
		// now look at the hazard difference; as expected, they are very similar; COV=0.015; min =0.910 and max=1.087
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B_Alt/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToNoSimilarRupturesCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		// This is the MFD constrained case; GR b=1, init model; this took 27 minutes
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, true, 5, false, false);
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, true, 2, false, false);
		
		
		// Even-fit solution: Don't include MFD constraint to exclude that as part of the error
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 100, MFD_TargetType.GR_b_1pt0, 0, false, 140, true, false);	

		
		// Do a range of GR bValues
		// 50 min; 2229 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, false, 5, false, true);
		// 35 min; 2270 non zero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt2, 1, false, 5, false, true);
		// 18 min; 2068 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt4, 1, false, 5, false, true);
		// 16 min; 2204 non zero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt6, 1, false, 5, false, true);
		// 14 min; 2265 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, true);
		// 17 min; 2366 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt2, 1, false, 5, false, true);
		// 18 min; 2372 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt4, 1, false, 5, false, true);
		// 23 min; 2432 nonzero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt6, 1, false, 5, false, true);

		// Plost histrogram of hazard diffs for the above range of b-values
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt6_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt4_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt2_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt"); // this one is closest to the mean
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt6_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt4_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt2_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//		    String label = name+"_GR_bValueRange0p0To1p6";
//		    String dirName = ROOT_PATH+"GR_bValueRange0p0To1p6_HazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}

		
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, false, 5, false, true);
		// 17 min; 1557 non-zero
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, false, 5, false, true);
		// 19 min; 2113 nonzero
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt6, 0.02, 1.0, false, 5, false, true);
		//
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt4, 0.02, 1.0, false, 5, false, true);
		// 14 minutes, 2220 nonzero
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt2, 0.02, 1.0, false, 5, false, true);

		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 0.02, 1.0, false, 5, false, true);
		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt6, 0.02, 1.0, false, 5, false, true);
		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt4, 0.02, 1.0, false, 5, false, true);
		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt2, 0.02, 1.0, false, 5, false, true);
	
		
		
		
		
		
		// ****************************************************************************************************
		// OLD SET OF RUNS AND TESTS BEFORE NUM SECT = 116 (RATHER THAN 117); moved to the OldResults050320 dir
		
//		// these are locations for ratios of NSHM Char to GR model (tapered slip and slip rate); extremes are for SA 1.0 2in50 (not PGA or 10in50)
//		LocationList locList1 = new LocationList();
//		ArrayList<String> locNameList1 = new ArrayList<String>();
//		locList1.add(new Location(34.825, -117.0));
//		locNameList1.add("minRatioLoc");
//		locList1.add(new Location(36.725, -117.0));
//		locNameList1.add("maxRatioLoc");
//		locList1.add(new Location(33.525, -117.85));
//		locNameList1.add("unityRatioLoc");
//
//		
//		// these are locations for ratios of FRESH max-rate model to min-rate model (tapered slip and slip rate); extremes are for SA 1.0 2in50 (not PGA or 10in50)
//		LocationList locList2 = new LocationList();
//		ArrayList<String> locNameList2 = new ArrayList<String>();
//		locList2.add(new Location(31.975, -117.0));
//		locNameList2.add("minRatioLoc");
//		locList2.add(new Location(34.875, -117.0));
//		locNameList2.add("maxRatioLoc");
//		locList2.add(new Location(35.325, -117.35));
//		locNameList2.add("unityRatioLoc");

//		
//		SimpleFaultInversion faultInversion = new SimpleFaultInversion();
		
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.ELLSWORTH_B);
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.HANKS_BAKUN_08);
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.SHAW_2009_MOD);
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP);
//		faultInversion.plotRateVsGRbValue(ScalingRelationshipEnum.ELLB_SQRT_LENGTH);
//		System.exit(0);
//		faultInversion.initData();
//		faultInversion.getGriddedRegion();
		
		// these have the exact same hazard (just different implied slip rates)
//		faultInversion.doNSHMP_Char_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, true, null, null);
//		faultInversion.doNSHMP_Char_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, locList1, locNameList1);
////		
//		// test that hazard is the same:
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//	    String fileName1 = ROOT_PATH+"NSHMP_Char_UNIFORM_UNIFORM_ELLSWORTH_B/hazardMaps/PGA_2in50.txt";
//	    String fileName2 = ROOT_PATH+"NSHMP_Char_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/PGA_2in50.txt";
//	    String dirName = ROOT_PATH+"NSHMP_Char_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//	    String label = "PGA_2in50_test_ratio";
//		faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);

		
		// these also have the exact same hazard (just different implied slip rates)
//		faultInversion.doNSHMP_GR_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, true, null, null);
//		faultInversion.doNSHMP_GR_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, false, locList1, locNameList1);
		
//		// compute ratios of Char to GR hazard
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"NSHMP_Char_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"NSHMP_GR_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"NSHMP_Char_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToGR";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		// these two are identical (only implied slip rate etc differ):
//		faultInversion.doNSHMP_CombinedSolution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);
//		faultInversion.doNSHMP_CombinedSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, true);
		
		// these two are identical (only implied slip rate etc differ); rates are 6% higher than char model above:
//		faultInversion.doMinRateSolution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, true, null, null);
//		faultInversion.doMinRateSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, true, locList2, locNameList2);

//		// compute ratios of Char to Min-rate model hazard; 10in50 for both SA01 and PGA have the biggest diff (~11%) due to curve having lower slope
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"NSHMP_Char_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MinRateSol_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MinRateSol_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToChar";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

//		//This misfits the very ends and a little elsewhere (2.7% higher away from ends)
//		faultInversion.doMaxRateSolution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, true);
//		// This fixes that (slip rates trimmed down at ends)
//		faultInversion.doMaxRateSolution(SlipRateProfileType.UNI_TRIM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, false);
//
//		// this still applies uniform dist along strike (because only two-section ruptures)
//		faultInversion.doMaxRateSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, true);
		
		// Max rate model that honors slip along fault; these two are identical (fit is not perfect, but looks good)
//		faultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE, false, null, null);	
//		faultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE, true, locList2, locNameList2);	
		
		// this shows the ratio max rate with slip-rate taper to uniform case; AveAbsValDiff (mean normalized) is less than 0.1 for all cases
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MaxRateSol_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";  // this has uniform slip distribution
//		    String dirName = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps";
//		    String label = name+"_RatioToUniformSlip";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		// Hazard difference between max-rate and min-rate solutions with tapered slip rate;
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MinRateSol_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"FRESH_TAPERED_TAPERED_ELLSWORTH_B_MAX_RATE/hazardMaps";
//		    String label = name+"_RatioToMinRate";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

//		// Hazard difference between min and max-rate solutions with uniform slip rate;
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MaxRateSol_UNIFORM_UNIFORM_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MinRateSol_UNIFORM_UNIFORM_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MaxRateSol_UNIFORM_UNIFORM_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToMinRate";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		
		// for the paper
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, 5, false, false);

		// this took 5.6 minutes (final energy E = 2)
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, 2, false, false);
		// this E=5 goes faster (3.7 min),and corresponds to RMS = sqrt(E/N) = sqrt(5/117) = 0.2 (ave Norm abs diff is less than this)
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, 5, false, false);
		// now do 10 runs (took 33 minutes); this looks more GR than TAPERED slip rate case
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, false);
		// this still produces edge effects (only use this slip rate distribution for max rate model)
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNI_TRIM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, false);

		// this took 10 min:
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, false);
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, false);

		// Apply default b=1 initial solution; this took 22 min; 
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, true, false);

		// this took <5 min; using the rup sampler leads to a different final solution (uses default b=1)
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, false, true);


		// this is to not over fit the data
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, 140, false);
		// this took 10 min
//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, 140, false, false);
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 25, 220, false, false);
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 25, 500, false, false);
//		faultInversion.doUnconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 25, 1000, false, false);

//		faultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, 140, true);

		// this should be the same as equiv unconstrainted version above (since wt=0 and zero initial state); looks good
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, false, 5, false);	

		// Use MFD as starting model only:
		// This hugely influences the final model, but slip rates are biased in the middle  It took 20 minutes
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, true, 5, false);
		// this fits it more closely; it took 33 minutes
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, true, 1, false);
		// this took 14 minutes; slip rates are a little worse
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, true, 1, true);
		
		// this took about an hour; slip rates slightly discrepancy at center and MFD half way between b =1 and b=0
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0, true, 1, false);
		// this show discrepancy at center
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0, true, 5, false);
		// OLD (GR rup sampler); this took twice as long for some reason; MFD more b=1ish; maybe sample from the target?
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0, true, 1, true);
			
		// this took 40 minutes; MFDs don't get close to starting model
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0, true, 1, false);	
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0, true, 5, false);	
		// this took 46 minutes; MFD fit no better
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0, true, 1, true);	

	
		// MFD WEIGHTED; target E higher because more data
		
		// b=1:
		
		// this took 22 minutes (2503 out of 6555 non zero)
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 2, false, true);
		// this took 50 minutes (48 minutes on second run) (2437 out of 6555 non zero); it's very strange that this is not faster, and that it's slower than the init model below
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 2, true);	
		// this took 35 min and all ruptures are nonzero (I commented out maps for all the different runs, only the mean)
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, true, 2, false, true);	
		// this took 34 minutes and all ruptures are nonzero; not any faster than not using rupSampler
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, true, 2, true);	
		
		// test higher energy case; took 14.5 minutes; 2320 are non zero
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);

//		// Hazard difference between E=5 over E=2; these are the same
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToE2case";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

//		// Hazard difference between zero init model and GR init model; only slight diff for 10in50 near ends
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_initSol_10_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToInitSolcase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		// Compare to individual runs that have no overlapping ruptures
		// single run with E=5; took 85 sec; 491 are non zero
		// I changed the name of the output dir to: MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B_Alt
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);
		// now I re-ran the run with test code that only allows zero-rate ruptures to be selected in the following (see "********** TEST **************" here);
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 1, false, 5, false, true);
		// I confirmed that these have none of the same ruptures
		// now look at the hazard difference; as expected, they are very similar
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B_Alt/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_1_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToNoSimilarRupturesCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		// b=0:
		
		// this took 4.5 hrs (3978 out of 6555 are non zero)
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, false, 2,false, true);

		// this took 4.8 hrs (3423 are non zero); rup sampler does not speed things up
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, false, 2,true);	

		// this took 9.7 hours; no zero ruptures; slight slip rate discrepancy
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, true, 2, false, false);	

		// try one run, higher energy and new rupSampler; this took 44 minutes (7 hrs for 10 runs); slip rates and event rates are a bit biased
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_0pt0, 1, true, 5, true, false);	
		// now try without rup sampler; 40 minutes; rup sampler does not speed things up and results the same
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_0pt0, 1, true, 5, false, false);	

		// b=-1 case:
		
		// gave up after ~10 hrs:
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1, false, 2, false);	
		// this took 2.4 hrs; 1739 rups are non zero
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1, false, 10, false, true);
		// this too 1.6 hrs, 33% faster
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1, false, 10, true, false);

		// this took 1.8 hours; 7 ruptures are zero
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1, true, 10, false);	
		// This took 2.0 hrs, no zero ruptures
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1, true, 10, true);	

		
		// test even fitting with MFD constraint; this also over-estimates slip rates (biased)
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 160, false, false);	
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, MFD_TargetType.GR_b_1pt0, 1, false, 160, false, false);	
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 25, MFD_TargetType.GR_b_1pt0, 1, false, 1000, false, false);	

		// USE THIS ONE IN PAPER: Don't include MFD constraint to exclude that as part of the error
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 100, MFD_TargetType.GR_b_1pt0, 0, false, 140, true, false);	

		
		// NNLS Solutions
		// this is instantaneous; 133 out of 6555 are non zero (2%)
//		faultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_minus1, 1);	
		// 133 non-zero here too
//		faultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt0, 1);	
		// 133 non-zero here too
//		faultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0, 1);	
		// removing MFD produces the minimum rate model
//		faultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0, 0);	

		
		// TOT RATE CONSTRAINED RUNS
		
		// b==1 rate:
			
		// this took ~13 minutes, 2225 are non zero
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, false, 5, false, true);
		// add initSol; this took 21 minutes; very close to GR; no zeros; this over estimates slip rates near fault center
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 5, false);
		// add sampler; this took less than 5 minutes with 1764 non-zero; closer to GR and less biased by sampling due to sampler
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, false, 5, true);
		// add sampler and initSol; this took about 10 minutes; no zeros; very close to GR, but with one dip at low M; this over estimates slip rates near fault center too
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 5, true);
		// Lower E; this took 16 minutes; slip rates are a bit better
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 2, true);
		// Even lower E; this took 15 minutes; and now slip rates are even better
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 1, true);
		
		// initSol make it very close to GR but with a weird dip at low mag and it biases slip rates (need lower E)
		// rupSampler speeds it up; makes it closer to GR and less biased by default sampler
		// initSol makes all non-zero but result is basically same as applying GR?
		
		
		// b=0 rate:
			
		// this took 15 minutes; slight negative b-value (increase with mag); biased by default sampler?
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, false, 5, false, true);
		// add initSol; 3.8 hrs; no zeros; about these same final MFD shape, but less biased by default sampler:
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, true, 5, false);
		// this took 12 minutes; closer to b=0; this seems to remove biased of the default sampler
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, false, 5, true);
		// 2.7 hrs; no zero rups; MFD is strongly rainbow with peak ~M7.4; this seems very strange
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, true, 5, true);
		
		// intoSol makes runs much longer
		// rupSampler speeds it up a bit and makes it closer to b=0 if no initSol
		
		// b=-1 rate:
		
		// this took 1.14 hours; 997 non-zero ruptures; looks pretty close to b=-1
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.02, 1.0, false, 5, false, true);
		// add intiSol; this took 1.9 hrs (19 hrs for 10 runs) so I gave up; 5 rups are zero
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0.02, 1.0, true, 5, false);
		// add rupSampler; this is to see whether sampling peaks in MFD are removed (they look about the same, maybe worse here); 8 hrs; 1080 non-zero, and took way longer
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.02, 1.0, false, 5, true);
		// add both; this took 11 hours; a bit faster?; no zero ruptures; MFD has inexplicable drop and jump at last mags, the single initSol model above is consistent with this.
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.02, 1.0, true, 5, true);

		// initSol slows things way down, and adds weird drop and jump at last mags
		// rupSampler neither improved MFD nor sped things up; not sure why
		
//		// ratio of tot-rate to MFD constrained; b=1; min=0.96 & max=1.04; ave norm abs diff is 0.01
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToMFD_ConstrCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
//		// ratio of tot-rate to MFD constrained; b=0; min=0.931 & max=1.066; max ave norm abs diff is ~0.026
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_0pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_0pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_0pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToMFD_ConstrCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		// ratio of tot-rate to MFD constrained; b=-1; min=0.955 & max=1.051; max ave norm abs diff is 0.018; COV = 0.02
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_minus1_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_10_finalE=10.0_GR_b_minus1_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"totRateConstrSA_10_finalE=5.0_GR_b_minus1_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToMFD_ConstrCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		
		
		// GR with b=0.32 case to compare with NSHM Combined
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt32, 1, false, 5, false, true);
		// ratio of b=0.32 case to NSHM Combined; min=0. & max=1.; max ave norm abs diff is 1.0; COV = 
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt32_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"NSHMP_Combined_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt32_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToNSHMcombinedCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		// the above NSHM Combined model does not fit slip rates perfectly, and has 10% and 30% hazard points > 10% away for PGA and 1secSA 2in50:
		
		// the following is NSHM Combined MFD and fitting slip rates; hazard ratios are significantly smaller (less than 2% of points are more than 10% away, on ave)
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.NSHM_Combined, 1, false, 5, true, true);
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt32_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_NSHM_Combined_wt1_TAPERED_TAPERED_ELLSWORTH_B_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt32_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps";
//		    String label = name+"_RatioToNSHMcombinedMFDCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		
		// Do a range of GR bValues
		// 26 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt2, 1, false, 5, false, true);
		// 20 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt4, 1, false, 5, false, true);
		// 16 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt6, 1, false, 5, false, true);
		// 15 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, true);
		// 16 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt2, 1, false, 5, false, true);
		// 19 min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt4, 1, false, 5, false, true);
		//  min
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt6, 1, false, 5, false, true);

//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt6_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt4_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt2_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_1pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt"); // this one is closest to the mean
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt6_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt4_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt2_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=2.0_GR_b_0pt0_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//		    String label = name+"_GR_bValueRange0p0To1p6";
//		    String dirName = ROOT_PATH+"GR_bValueRange0p0To1p6_HazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}

		
		// Do different scaling relationships for b=0.8
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.HANKS_BAKUN_08, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, true);
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLB_SQRT_LENGTH, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, true);
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP, 10, MFD_TargetType.GR_b_0pt8, 1, false, 105, false, true);
		// this skips a mag bin so E cannot get below 100
//		faultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.SHAW_2009_MOD, 10, MFD_TargetType.GR_b_0pt8, 1, false, 105, false, true);
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=105.0_GR_b_0pt8_wt1_TAPERED_TAPERED_SHAW_2009_MOD/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=105.0_GR_b_0pt8_wt1_TAPERED_TAPERED_SHAW_CONST_STRESS_DROP/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLB_SQRT_LENGTH/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_HANKS_BAKUN_08/hazardMaps/"+name+".txt"); // this one is closest to the mean
//		    String label = name+"_ScalingRelationshipOptions";
//		    String dirName = ROOT_PATH+"ScalingRelationshipOptionsHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}

		// Do alternate SlipAlongRuptureModel
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 1, false, 5, false, true);
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_UNIFORM_ELLSWORTH_B/hazardMaps/"+name+".txt");
//		    String label = name+"_SlipAlongRuptureOptions";
//		    String dirName = ROOT_PATH+"SlipAlongRuptureOptionsHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}


		
		// Do alternate SlipRateProfile (this needs higher energy because of necessary slip-rate misfits at ends)
//		faultInversion.doMFDconstrainedSA(false, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt8, 1, false, 50, true, true);
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=5.0_GR_b_0pt8_wt1_TAPERED_TAPERED_ELLSWORTH_B/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"MFDconstrSA_10_finalE=50.0_GR_b_0pt8_wt1_UNIFORM_TAPERED_ELLSWORTH_B_rupSamp/hazardMaps/"+name+".txt");
//		    String label = name+"_SlipRateProfileOptions";
//		    String dirName = ROOT_PATH+"SlipRateProfileOptionsHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}

		// This compares different GMMs for the GR b=0.8 case
//		faultInversion.computeHazardForGMM();
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"GMM_OptionsHazardVariability/ASK2014_HazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"GMM_OptionsHazardVariability/CY2014_HazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"GMM_OptionsHazardVariability/CB2014_HazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"GMM_OptionsHazardVariability/BSSA2014_HazardMaps/"+name+".txt");
//		    String label = name+"_GMM_OptionsHazardVariability";
//		    String dirName = ROOT_PATH+"GMM_OptionsHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}
		// here's a test ratio to explore why ASK is an outlier; it's at far distances, which I confirmed for median values with the GMM gui
//		faultInversion.shortFault=true; faultInversion.makeFaultSectionDataList();
//		File file = new File(ROOT_PATH+"TestRatioASKvsCY"); file.mkdir();
//		faultInversion.makeHazardMapRatio(ROOT_PATH+"GMM_OptionsHazardVariability/ASK2014_HazardMaps/1.0secSA_2in50.txt"
//				, ROOT_PATH+"GMM_OptionsHazardVariability/CY2014_HazardMaps/1.0secSA_2in50.txt"
//				, "TestRatioASKvsCY", ROOT_PATH+"TestRatioASKvsCY", true);

		
		
		// max Rate case: (THIS WOULD NEVER BE USED?)
		
		// this took 8.5 hrs; 2242 are nonzero
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, false, 5, false);
		// add initSol; 90 minutes
//		faultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 5, false);
		// this took 1 min; it just replicates the true max-rate solution
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MAX_RATE, 0.02, 1.0, false, 5, true);
		// this took ~4 min
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 5, true);

		// this is a check, as expected, has a 10e-22 initial energy (no iterations necessary)
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.UNI_TRIM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 2, false);

		// min Rate case: (THIS WOULD NEVER BE USED?)
		
		// this takes forever (as expected); gave up after an hour (only got to E=24) 
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, false, 10, false);
		// this takes a long time too, and not sure why initial state doesn't have lower E; following test implies there is biased moment rate due to amny rups in last mag bin
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, true, 2, false);
//		faultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MIN_RATE);
		// 7 minutes here
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, false, 10, true);
		// 5 min
//		faultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, true, 10, true);

		
		// NEED TO UNDERSTAND WHY THESE FRESH SOLUTION ARE NOT BETTER FITS
		// E=1293 for this:
//		faultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt8, false, null, null);	
		// try uniform slip; E=483
//		faultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt8, false, null, null);	
		// E=747 (better than FRESH)
//		faultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt8);
		// E=184; this is because the standard floater model has a taper that happens to match the slip-rate taper here.
//		faultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt8);

		// SUNFISH MODELS
		// this has E=1908; this has a negative b-value (rates increase with magnitude)
//		faultInversion.doSUNFiSH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B);	
		// this has E=2231; looks about the same
//		faultInversion.doSUNFiSH_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);	

		
		
		// SEGMENTATION CONSTRAINTS
		faultInversion.shortFault = false;	// use longer fault for the rest of these
		
		ArrayList<SlipRateSegmentationConstraint> segConstr = new ArrayList<SlipRateSegmentationConstraint>();
		segConstr.add(new SlipRateSegmentationConstraint("sect58", 58, 0.1, 0.001));	
//		 this has high energy because of a zero MFD bin and highest mag rates conflict with segmentation
//		faultInversion.doSlipRateSegmentedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 1, 600, MFD_TargetType.GR_b_0pt8, 1.0, false, segConstr, "s58_0.1", false, false);
		// inti model:
//		faultInversion.doSlipRateSegmentedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//				ScalingRelationshipEnum.ELLSWORTH_B, 1, 600, MFD_TargetType.GR_b_0pt8, 1.0, true, segConstr, "s58_0.1", false, false);
//		faultInversion.doSlipRateSegmentedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//				ScalingRelationshipEnum.ELLSWORTH_B, 10, 100, MFD_TargetType.GR_b_0pt8, 0.0, true, segConstr, "s58_0.1", false, false);
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "s58_0.1", true);
//		faultInversion.doSlipRateSegmentedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8_ForSegmented, 1.0, false, segConstr, "s58_0.1");
		// PREFERRED MODEL:
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "s58_0.1", true, true);

//		ArrayList<SlipRateSegmentationConstraint> segConstr = new ArrayList<SlipRateSegmentationConstraint>();
//		segConstr.add(new SlipRateSegmentationConstraint("sect58", 58, 0.0, 0.000001));	
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "s58_0.0", true, true);

		// 50% segmentation
//		ArrayList<SlipRateSegmentationConstraint> segConstr = new ArrayList<SlipRateSegmentationConstraint>();
//		segConstr.add(new SlipRateSegmentationConstraint("sect58", 58, 0.5, 0.001));	
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "s58_0.5", true, true);


		// no constraint:
//		ArrayList<SlipRateSegmentationConstraint> segConstr = new ArrayList<SlipRateSegmentationConstraint>();
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//				ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "noSeg", true);
//		faultInversion.doSlipRateSegmentedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//				ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, MFD_TargetType.GR_b_0pt8, 0.0, true, segConstr, "noSeg");
		// PREFERRED MODEL:
//		faultInversion.doSlipRateSegmentedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstr, "noSeg", true, true);

		
		// Ratio of no RupSampler cases for total segmentation
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.0_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.0_rupSamp/hazardMaps";
//		    String label = name+"_RatioToNoSegCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
//		// The treats the three segmentation levels as branches
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.0_rupSamp/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg_rupSamp/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt");
//		    String label = name+"_SlipRateSegmentationHazardVariability";
//		    String dirName = ROOT_PATH+"SlipRateSegmentationHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}
		
		// Ratio of no RupSampler cases
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1/hazardMaps";
//		    String label = name+"_RatioToNoSegCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		// Ratio of RupSampler cases; PREFERRED MODEL
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1_rupSamp/hazardMaps";
//		    String label = name+"_RatioToNoSegCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		// Test Ratio (influence of RupSampler); hazard results are close enough
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.1/hazardMaps";
//		    String label = name+"_RatioToRupSamplerCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_noSeg/hazardMaps";
//		    String label = name+"_RatioToRupSamplerCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		
		// do the 50% segmentation case using section rate constraint; need higher energy because slip rate not fit at sect 58
//		ArrayList<SectionRateConstraint> sectConstrList = new ArrayList<SectionRateConstraint>();
//		sectConstrList.add(new SectionRateConstraint("50percSeg", 58, 0.0037, 0.0037/500.)); // the rate from the slip rate case
//		faultInversion.doSectRateConstrSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 40, MFD_TargetType.GR_b_0pt8, 0.0, false, sectConstrList, "s58_0.5", true, true);
		
		// confirm that the two different ways of doing the 50% segmentation produce the same hazard
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SectRateConstrSA_10_finalE=40.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps";
//		    String label = name+"_RatioToSectRateSegCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}

		
		// do the 50% segmentation case using segmentation constraint
//		ArrayList<SegmentationConstraint> segConstrList = new ArrayList<SegmentationConstraint>();
//		segConstrList.add(new SegmentationConstraint("50percSeg", 58, 59, 0.0037, 0.0037/500.)); // the rate from the slip rate case
//		faultInversion.doSegConstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 10, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstrList, "s58_0.5", true, false);

//		// confirm that this third way of doing the 50% segmentation produce the same hazard; they do
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SectRateConstrSA_10_finalE=40.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps";
//		    String label = name+"_RatioToSectRateConstrCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
	

		// do the 50% segmentation case using segmentation constraint
//		ArrayList<SegmentationConstraint> segConstrList = new ArrayList<SegmentationConstraint>();
//		segConstrList.add(new SegmentationConstraint("50percSeg", 58, 59, 0.02, 0.02/500., true)); // the rate from the slip rate case
//		faultInversion.makeHazardMapsForAllSolutions=true;
//		faultInversion.doSegConstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 5, MFD_TargetType.GR_b_0pt8, 0.0, false, segConstrList, "s58_slipRate0.5", true, true);

		
//		// confirm that this fourth way of doing the 50% segmentation produces the same hazard; it does
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps";
//		    String label = name+"_RatioToSegConstrCase";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
//		// compare hazard maps for all four ways of doing the 50% segmentation
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//			ArrayList<String> fileNameList = new ArrayList<String>();
//			fileNameList.add(ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"LongFlt_SlipRateSegSA_10_finalE=10.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt");
//			fileNameList.add(ROOT_PATH+"LongFlt_SectRateConstrSA_10_finalE=40.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_0.5_rupSamp/hazardMaps/"+name+".txt");
//		    String label = name+"_4_Diff_0.5_SegmentationHazardVariability";
//		    String dirName = ROOT_PATH+"4_Diff_0.5_SegmentationHazardVariability";
//			File tempFile = new File(dirName); tempFile.mkdirs();
//		    GriddedGeoDataSet[] dataArray = faultInversion.readMultipleHazardMapGriddedData(fileNameList);
//		    faultInversion.plotNormPDF_FromMultHazMaps(dataArray, label, dirName, true);
//		}

		
		// For comparison, I then computed hazard map PDFs for all 10 runs of LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp
		
		// I also did another run of LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp to compare the ave solutions:
		// this has about the same variability as the above 4 50-percent segmentation approaches, implying the latter produce identical results
//		faultInversion.shortFault=false; faultInversion.makeFaultSectionDataList();
//		String[] nameArray = {"PGA_2in50", "PGA_10in50", "1.0secSA_2in50", "1.0secSA_10in50", "1.0secSA_RTGM"};
//		for(String name: nameArray) {
//		    String fileName1 = ROOT_PATH+"altLongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String fileName2 = ROOT_PATH+"LongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps/"+name+".txt";
//		    String dirName = ROOT_PATH+"altLongFlt_SegConstrSA_10_finalE=5.0__GR_b_0pt8_wt0_TAPERED_TAPERED_ELLSWORTH_B_s58_slipRate0.5_rupSamp/hazardMaps";
//		    String label = name+"_RatioToAnotherIdenticalRun";
//			faultInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName, true);			
//		}
		
		
		
		// SECTION RATES CONSTRAINTS AT TWO POINTS
//		ArrayList<SectionRateConstraint> sectConstrList = new ArrayList<SectionRateConstraint>();
//		double rate1=0.008/2.0, rate2=0.008*2.0; // 0.03 is the rate at these sections in the unconstrained result
//		sectConstrList.add(new SectionRateConstraint("Test Sect Constr 1", 29, rate1, rate1/100));
//		sectConstrList.add(new SectionRateConstraint("Test Sect Constr 2", 87, rate2, rate2/100));

//		faultInversion.doSectRateConstrSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 60, MFD_TargetType.GR_b_0pt8, 0.0, false, sectConstrList, "TwoPoints", true, false);

//		faultInversion.doSectRateConstrSmoothedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, 
//		ScalingRelationshipEnum.ELLSWORTH_B, 10, 90, MFD_TargetType.GR_b_0pt8, 0.0, false, sectConstrList, "TwoPoints", true, false);

		
		
		
		


		
		
		

		
	}
}
