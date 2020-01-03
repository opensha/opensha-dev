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
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
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
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.simulatedAnnealing.completion.CompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.CompoundCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.EnergyCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.IterationCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.params.CoolingScheduleType;
import scratch.UCERF3.simulatedAnnealing.params.GenerationFunctionType;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.ScalingRelationshipEnum;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.SlipAlongRuptureModelEnum;

/**
 * This class runs the inversion.
 * 
 * Slip rate and event rate standard deviations are computed as the difference of the 95% confidence bounds divided by 4.
 * 
 * @author 
 *
 */
public class SingleFaultInversion {

	final static boolean D = true;	// debugging flag
	
	public final static String ROOT_PATH = "src/scratch/ned/FSS_Inversion2019/";
	final static String ROOT_DATA_DIR = "src/scratch/ned/FSS_Inversion2019/data/"; // where to find the data

	
	// These values are the same for all fault sections
	final static double UPPER_SEIS_DEPTH = 0;
	final static double LOWER_SEIS_DEPTH = 14.;
	final static double FAULT_DIP = 90; // FIX MAX_SUBSECT_LENGTH_KM calc below if this is changed
	final static double FAULT_RAKE = 0.0;
	final static double FAULT_SLIP_RATE = 30;	// mm/yr
	final static double FAULT_LENGTH_KM = 14*29;
	
	final static int NUM_SUBSECT_PER_RUP = 4;
	final static double MAX_SUBSECT_LENGTH_KM = (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)/NUM_SUBSECT_PER_RUP;
	
	final static double hazGridSpacing = 0.05;
		
	ArrayList<FaultSectionPrefData> faultSectionDataList;
	ArrayList<SegRateConstraint> sectionRateConstraints;
	int[][] rupSectionMatrix;
		
	FaultSectionPrefData parentFaultData;

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
		GR_b_1pt0,
		GR_b_0pt8,
		GR_b_0pt0,
		GR_b_minus1,
		MAX_RATE,  // all in minimum magnitude bin
		MIN_RATE,  // all in maximum magnitude bin
		M7pt25only,
		NONE;
	}

	public enum SlipRateProfileType {
		UNIFORM,
		TAPERED,
		UNIFORM_TRIMMED;
	}
	
	public enum SA_InitialStateType {
		ALL_ZEROS,
		A_PRIORI_RATES,
		FROM_MFD_CONSTRAINT;
	}

	
	double hazCurveLnMin = Math.log(0.001);
	double hazCurveLnMax = Math.log(10);
	int hazCurveNum = 20;
	double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);
	
	// Adjustable Parameters
	String dirName;
	String solutionName; 						// Inversion name
	boolean wtedInversion; 						// Whether to apply data weights
	SlipRateProfileType slipRateProfile;
	SlipAlongRuptureModelEnum slipModelType;	// Average slip along rupture model
	ScalingRelationshipEnum scalingRel; 		// Scaling Relationship
	double relativeSectRateWt; 					// Section Rate Constraints wt
	String segmentationConstrFilename; 			// = ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE;
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
	InversionSolutionType solutionType;
	long randomSeed;							// for SA reproducibility 
	int numSolutions; 							// this is ignored for NON_NEGATIVE_LEAST_SQUARES which only has one possible solution
	CoolingScheduleType saCooling;
	GenerationFunctionType perturbationFunc;
	CompletionCriteria completionCriteria;
	SA_InitialStateType initialStateType;
//	String rupRatesFromFileDirName;				//if solution is from a file (previous solution); name example: ROOT_PATH+"OutDir/";
	double magAareaAleatoryVariability;			// the gaussian sigma for the mag-area relationship
	// Data to plot and/or save:
	boolean popUpPlots;							// this tells whether to show plots in a window (set null to turn off; e.g., for HPC)
	boolean doDataFits;
	boolean doMagHistograms;
	boolean doNonZeroRateRups;
    boolean doSectPartMFDs;
    boolean doSectSlipDists;
	Location hazCurveLoc; 						// to make hazard curve (set loc=null to ignore)
    String hazCurveLocName;
    boolean makeHazardMaps; 					// to make hazard maps
	double saPeriodForHaz;						// set as 0.0 for PGA


	
	public SingleFaultInversion() {
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
	private GriddedGeoDataSet readHazardDataFromFile(String fileName) {
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
	
	
	private void tempSetNSHMP_GR_Model(ScalingRelationshipEnum scalingRel) {
		IncrementalMagFreqDist targetMFD = getTargetMFD(scalingRel, MFD_TargetType.GR_b_1pt0);	// this should be b=0.8
		
		// make magLogAreaFunc to interpolate from because scalingRel doesn't (yet) give area from mag
		double minLogArea = Math.log10(10e6); // 10 km squared)
		double maxLogArea = Math.log10(10e3 *10e6); // 10 km by 10k km
		EvenlyDiscretizedFunc magLogAreaFunc = new EvenlyDiscretizedFunc(minLogArea, maxLogArea, 101);
		for(int i=0;i<magLogAreaFunc.size();i++) {
			double mag = scalingRel.getMag(Math.pow(10, magLogAreaFunc.getX(i)), Double.NaN, Double.NaN);
			magLogAreaFunc.set(i,mag);
		}
		System.out.println(magLogAreaFunc.toString());
			
		
//		WC1994_MagLengthRelationship wcScRel = new WC1994_MagLengthRelationship();
//		wcScRel.setRake(0.0);
		for(int i=0;i<targetMFD.size();i++) {
//			System.out.println((float)targetMFD.getX(i)+"\t"+(float)wcScRel.getMedianLength(targetMFD.getX(i)));
			double areaKm = Math.pow(10.0, magLogAreaFunc.getFirstInterpolatedX(targetMFD.getX(i)))/1e6;
			double lengthKm = areaKm/(LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH);
			System.out.println((float)targetMFD.getX(i)+"\t"+(float)lengthKm+"\t"+(float)(lengthKm/3.5));

		}
		System.exit(0);
	}
	
	
	private void makeFaultSectionDataList() {
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
		faultSectionDataList = parentFaultData.getSubSectionsList(MAX_SUBSECT_LENGTH_KM);
		
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
		
		System.out.println("parentFaultData MoRate: "+parentFaultData.calcMomentRate(true));
		
		int numSect = faultSectionDataList.size();
		
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
		case UNIFORM_TRIMMED:  // linear taper over first and last 3 subsections
			double moRateOrig = 0;
			for(FaultSectionPrefData data : faultSectionDataList) {
				moRateOrig += data.calcMomentRate(true);
			}
			FaultSectionPrefData data;
			data = faultSectionDataList.get(0); data.setAveSlipRate(0.25*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.25*data.getOrigSlipRateStdDev());
			data = faultSectionDataList.get(1); data.setAveSlipRate(0.50*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.50*data.getOrigSlipRateStdDev());
			data = faultSectionDataList.get(2); data.setAveSlipRate(0.75*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.75*data.getOrigSlipRateStdDev());
			data = faultSectionDataList.get(faultSectionDataList.size()-1); data.setAveSlipRate(0.25*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.25*data.getOrigSlipRateStdDev());
			data = faultSectionDataList.get(faultSectionDataList.size()-2); data.setAveSlipRate(0.50*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.50*data.getOrigSlipRateStdDev());
			data = faultSectionDataList.get(faultSectionDataList.size()-3); data.setAveSlipRate(0.75*data.getOrigAveSlipRate()); data.setSlipRateStdDev(0.75*data.getOrigSlipRateStdDev());
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
		
		int numRup =0;
		int curNumRup = numSect - (NUM_SUBSECT_PER_RUP-1);
		while(curNumRup>0) {
			numRup += curNumRup;
			curNumRup -= 1;
		}

		int rupIndex=0;

//		rupSectionMatrix = new int[numSect][numSect - (NUM_SUBSECT_PER_RUP-1)];
//		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<NUM_SUBSECT_PER_RUP+1; curSectPerRup++)
		rupSectionMatrix = new int[numSect][numRup];
		for(int curSectPerRup = NUM_SUBSECT_PER_RUP; curSectPerRup<numSect+1; curSectPerRup++)
			for(int s=0;s<numSect-curSectPerRup+1;s++) {
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
		
		sectionRateConstraints = new ArrayList<SegRateConstraint>();
	}
	
	public ArrayList<FaultSectionPrefData> getFaultSectionDataList() {
		return faultSectionDataList;
	}
	
	public ArrayList<SegRateConstraint> getSectionRateConstraints() {
		return sectionRateConstraints;
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
		
		// set the forecast duration
		double duration = 50;

		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		int numSol = fltSysRupInversion.getNumSolutions();
		if(numSol>1) {
			ArbDiscrEmpiricalDistFunc_3D curvesFromMultRunsFunc_3D = new ArbDiscrEmpiricalDistFunc_3D(hazCurveLnMin,hazCurveNum,hazCurveDelta);
			for(int i=0;i<numSol;i++) {
				EvenlyDiscretizedFunc func = computeHazardCurveLnX(fltSysRupInversion.getFaultSystemSolution(i), location, saPeriod, duration);
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
		EvenlyDiscretizedFunc curveLogXvalues = computeHazardCurveLnX(fltSysRupInversion.getFaultSystemSolution(), location, saPeriod, duration);
		// convert to linear x valules
		ArbitrarilyDiscretizedFunc curveLinearXvalues = new ArbitrarilyDiscretizedFunc();
			for (int i = 0; i < curveLogXvalues.size(); ++i)
				curveLinearXvalues.set(Math.exp(curveLogXvalues.getX(i)), curveLogXvalues.getY(i));

		
		double twoIn50value = Math.exp(curveLogXvalues.getFirstInterpolatedX(0.02));
		double tenIn50value = Math.exp(curveLogXvalues.getFirstInterpolatedX(0.1));
		curveLinearXvalues.setInfo("2in50 value: "+twoIn50value+"\n"+"10in50 value: "+tenIn50value+
				"\nLocation: "+location.getLatitude()+", "+location.getLongitude());
		
		
		// make the plot
		plottingFuncsArray.add(curveLinearXvalues);
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve";
		String xAxisLabel = imtString;
		String yAxisLabel = "Probability (in "+duration+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = true;
		boolean logY = true;

		fltSysRupInversion.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, plotTitle, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}

	
	
	/**
	 */
	private ArbitrarilyDiscretizedFunc old_computeHazardCurve(FaultSystemSolution faultSystemSolution, Location location, 
			double saPeriod, double forecastDuration) {
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(faultSystemSolution);
		erf.setName("WasatchERF");
		erf.getTimeSpan().setDuration(forecastDuration);
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
		
		// set the IMT & curve x-axis values
		ArbitrarilyDiscretizedFunc curveLinearXvalues;
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
			curveLinearXvalues = IMT_Info.getUSGS_PGA_Function();
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
			curveLinearXvalues = IMT_Info.getUSGS_SA_01_AND_02_Function();
		}
		
		// make the site object and set values
		Site site = new Site(location);
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
			System.out.println(param.getName()+"\t"+param.getValue());
		}
					
		ArbitrarilyDiscretizedFunc curveLogXvalues = HazardCurveSetCalculator.getLogFunction(curveLinearXvalues);
	
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
		
		return curveLogXvalues;
	}
	
	
	/**
	 */
	private EvenlyDiscretizedFunc computeHazardCurveLnX(FaultSystemSolution faultSystemSolution, Location location, 
			double saPeriod, double forecastDuration) {
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(faultSystemSolution);
		erf.setName("WasatchERF");
		erf.getTimeSpan().setDuration(forecastDuration);
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
	 * This is makes to hazard maps, one for 10% in 50-year and one for 2% in 50-year ground motions
	 * @param fltSysRupInversion
	 * @param saPeriod - set as 0 for PGA; this will crash if SA period is not supported
	 * @param dirName
	 * @param popupWindow
	 */
	public void makeHazardMaps(FaultSystemRuptureRateInversion fltSysRupInversion, double saPeriod, String dirName, boolean popupWindow) {
		
		File subDirFile = null;
		if(dirName != null) {
		    subDirFile = new File(dirName+"/hazardMaps");
		    subDirFile.mkdirs();
		}

		double duration = 50;
		double[] probabilityArray = {0.02, 0.10};
		String[] probabilityNameArray = {"2in50", "10in50"};
		
		String imtString = "PGA";
		if(saPeriod != 0)
			imtString = saPeriod+"secSA";
		
		//map region
		GriddedRegion region = getGriddedRegion();
		
		// make gridded data sets
		GriddedGeoDataSet[] griddedDataSetArray = new GriddedGeoDataSet[probabilityArray.length];
		for(int i=0;i<probabilityArray.length;i++) {
			griddedDataSetArray[i] = new GriddedGeoDataSet(region, true);
		}
		
		// Create the ERF
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fltSysRupInversion.getFaultSystemSolution());
		erf.setName("WasatchERF");
		
		// set the forecast duration
		erf.getTimeSpan().setDuration(duration);
		
		// update forecast
		erf.updateForecast();
		
		// set the attenuation relationship (GMPE)
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		imr.setParamDefaults();
		
		// set the IMT & curve x-axis values
		ArbitrarilyDiscretizedFunc curveLinearXvalues;
		if(saPeriod == 0) {
			imr.setIntensityMeasure(PGA_Param.NAME);
			curveLinearXvalues = IMT_Info.getUSGS_PGA_Function();
		}
		else {
			SA_Param saParam = (SA_Param)imr.getParameter(SA_Param.NAME);
			saParam.getPeriodParam().setValue(saPeriod);
			imr.setIntensityMeasure(saParam);
			curveLinearXvalues = IMT_Info.getUSGS_SA_01_AND_02_Function();
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
		for(int i=0;i<griddedDataSetArray[0].size();i++) {
			System.out.println(i+" of "+griddedDataSetArray[0].size());
			site.setLocation(griddedDataSetArray[0].getLocation(i));
			calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
			// get the points on the hazard curve
			for(int p=0;p<probabilityArray.length;p++) {
				double logVal = curveLogXvalues.getFirstInterpolatedX(probabilityArray[p]);
				griddedDataSetArray[p].set(i, Math.exp(logVal));
			}
		}
		
		try {
			for(int p=0;p<probabilityArray.length;p++) {
				CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(griddedDataSetArray[p].getMinZ(), griddedDataSetArray[p].getMaxZ());
//				CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(gridData.getMinZ(), gridData.getMaxZ());
				ArrayList<LocationList> faults = FaultBasedMapGen.getTraces(faultSectionDataList);
				double[] values = new double[faults.size()];
				for(int i=0;i<values.length;i++)
					values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;
				
				String label = imtString+"_"+probabilityNameArray[p];
				
				GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, griddedDataSetArray[p], hazGridSpacing, region, true, label);
//				map.setGenerateKML(true); // this tells it to generate a KML file that can be loaded into Google Earthe
				
				try {
					FaultBasedMapGen.SAVE_ZIPS=true;
					FaultBasedMapGen.plotMap(subDirFile, label, popupWindow, map);	
				} catch (GMT_MapException e) {
					e.printStackTrace();
				}
				
				// write out values to a text file
				FileWriter fw = new FileWriter(dirName+"/hazardMaps/"+label+".txt");
				fw.write("index\tvalue\tlatitude\tlongitude\n");
				for(int i=0;i<griddedDataSetArray[p].size(); i++)	 {
					Location loc = griddedDataSetArray[p].getLocation(i);
					fw.write(i+"\t"+griddedDataSetArray[p].get(i)+"\t"+loc.getLatitude()+"\t"+loc.getLongitude()+"\n");
				}
				fw.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
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
		
		GriddedGeoDataSet data1 = readHazardDataFromFile(fileName1);
		GriddedGeoDataSet data2 = readHazardDataFromFile(fileName2);
		
		// put log10 ratio in first data set:
		for(int i=0;i<data1.size();i++)
			data1.set(i, Math.log10(data1.get(i)/data2.get(i)));
				
		try {
			CPT cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance().rescale(-1.0, 1.0);
			ArrayList<LocationList> faults = FaultBasedMapGen.getTraces(faultSectionDataList);
			double[] values = new double[faults.size()];
			for(int i=0;i<values.length;i++)
				values[i] = FaultBasedMapGen.FAULT_HIGHLIGHT_VALUE;

			GMT_Map map = FaultBasedMapGen.buildMap(cpt, faults, values, data1, hazGridSpacing, getGriddedRegion(), true, "Log10 "+label);

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

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}


	private GriddedRegion getGriddedRegion() {
		double minLat = 39.5;
		double maxLat = 41.8;
		double minLon = -113.1;
		double maxLon = -110.6;
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
		case GR_b_1pt0:
			GutenbergRichterMagFreqDist gr_b1 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b1.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 1.0);
			gr_b1.setName("GR, b=1");
			mfd = gr_b1;
			break;
		case GR_b_0pt8:
			GutenbergRichterMagFreqDist gr_b0p8 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0p8.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.8);
			gr_b0p8.setName("GR, b=0.8");
			mfd = gr_b0p8;
			break;
		case GR_b_0pt0:
			GutenbergRichterMagFreqDist gr_b0 = new GutenbergRichterMagFreqDist(minMFD_Mag,num,0.1);
			gr_b0.setAllButTotCumRate(minMFD_Mag, maxMFD_Mag, toMoRate, 0.0);
			gr_b0.setName("GR, b=1");
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
		case NONE:
			mfd = null;
			break;
		}

		return mfd;
	}
	
	
	/**
	 * This sets/resets default parameter values
	 */
	public void setDefaultParameterValuess() {
		 dirName = ROOT_PATH+"NoNameOutput";
		 solutionName = "No Name Solution"; // Inversion name
		 wtedInversion = true; // Whether to apply data weights
		 slipRateProfile = SlipRateProfileType.TAPERED;
		 slipModelType = SlipAlongRuptureModelEnum.TAPERED; // Average slip along rupture model
		 scalingRel = ScalingRelationshipEnum.ELLSWORTH_B; // Scaling Relationship
		 relativeSectRateWt=0; // Section Rate Constraints wt
		 segmentationConstrFilename = null; // = ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE;
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
		 doDataFits=true;
		 doMagHistograms=true;
		 doNonZeroRateRups=true;
	     doSectPartMFDs=true;
	     doSectSlipDists=true;
		 hazCurveLoc = null; // to make hazard curve (set loc=null to ignore)
	     hazCurveLocName = "Hazard Curve";
	     makeHazardMaps = false; //// to make hazard maps
		 saPeriodForHaz = 1.0;	// set as 0.0 for PGA
	}
	
	
	/**
	 * This generates a solution and diagnostic info for the current parameter settings
	 */
	public void getSolution(boolean rePlotOny) {	
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
		ArrayList<SegRateConstraint> sectionRateConstraints = getSectionRateConstraints();
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
				sectionRateConstraints, 
				rupSectionMatrix, 
				slipModelType, 
				scalingRel, 
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
				segmentationConstrFilename,
				relative_segmentationConstrWt,
				totalRateConstraint,
				totalRateSigma,
				relativeTotalRateConstraintWt,
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
	    
		if(rePlotOny) { // overide solution type
			solutionType = InversionSolutionType.FROM_FILE;
		}


		// set the intial state if SA
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
	    		initialState = fltSysRupInversion.getRupRatesForTargetMFD(targetMFD);
//	    		for(int i=0;i<initialState.length;i++)
//	    			initialState[i] *= 0.25;
	    		break;
	    	}
	    }


	    switch (solutionType) {
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
	    		fltSysRupInversion.doInversionSA(completionCriteria, initialState, randomSeed, saCooling, perturbationFunc);
	    	}
	    	else if(numSolutions>1) {
	    		if(D) System.out.println("SIMULATED_ANNEALING; numSolutions="+numSolutions);
	    		fltSysRupInversion.doInversionSA_MultTimes(completionCriteria, initialState, randomSeed, numSolutions, dirName, saCooling,perturbationFunc);
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
		
		if(!rePlotOny) { // don't overwrite info about the run
			fltSysRupInversion.addToModelRunInfoString("\n"+runTimeString+"\n");
			// write results to file
			fltSysRupInversion.writeInversionRunInfoToFile(dirName);			
		}
					
		if(doDataFits) fltSysRupInversion.writeAndOrPlotDataFits(dirName, popUpPlots);
		if(doMagHistograms) fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots);
		if(doNonZeroRateRups) fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots);
		if(doSectPartMFDs) fltSysRupInversion.writeAndOrPlotSectPartMFDs(dirName, popUpPlots);
		if(doSectSlipDists) fltSysRupInversion.writeAndOrPlotSectSlipDistributions(dirName, popUpPlots);
	    
	    // hazard curve:
		if(hazCurveLoc != null) {
			writeAndOrPlotHazardCurve(fltSysRupInversion, hazCurveLoc, saPeriodForHaz, dirName, popUpPlots, hazCurveLocName);
		}
	    
	    // second parameter here is SA period; set as 0.0 for PGA:
		if(makeHazardMaps)
			makeHazardMaps(fltSysRupInversion, saPeriodForHaz, dirName, popUpPlots);
	}

	
	
	public void doNSHMP_GR_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = MFD_TargetType.GR_b_0pt8; // Target MFD Constraint
		solutionName = "NSHMP_GR_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}
	
	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doNSHMP_Char_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		magAareaAleatoryVariability = 0.12;	// gaussian distribution sigma
		// make the solution
		solutionType = InversionSolutionType.APPLY_SOLUTION;
		solutionRatesToApplyArray = new double[6555];
		double maxRupAreaKmSq = (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*FAULT_LENGTH_KM;
		double mag = scalingRel.getMag(maxRupAreaKmSq*1e6, FAULT_LENGTH_KM*1e3, (LOWER_SEIS_DEPTH-UPPER_SEIS_DEPTH)*1e3);
		double maxMFD_Mag = FaultSystemRuptureRateInversion.roundMagTo10thUnit(mag);
		double toMoRate = FaultMomentCalc.getMoment(maxRupAreaKmSq*1e6, FAULT_SLIP_RATE*1e-3);
		double magRoundingCorr = MagUtils.magToMoment(maxMFD_Mag)/MagUtils.magToMoment(mag);
		solutionRatesToApplyArray[6554] = (toMoRate*magRoundingCorr/1.071748)/MagUtils.magToMoment(maxMFD_Mag); // 1.071748 is for aleatory variability in magnitude
		solutionName = "NSHMP_Char_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}
	
	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doMinRateSolution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		// make the solution
		solutionType = InversionSolutionType.APPLY_SOLUTION;
		solutionRatesToApplyArray = new double[6555];
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
		
		solutionRatesToApplyArray[6554] = toMoRate*magRoundingCorr/MagUtils.magToMoment(maxMFD_Mag); 
		solutionName = "MinRateSol_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString()+"_"; // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}

	
	/**
	 * The remaining slip-rate discrepancy is from floater end effects, as only one rupture hits the first and last section, whereas 
	 * most sections participate in 4 different ruptures; the average slip-rate discrepancy is zero.
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doMaxRateSolution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = MFD_TargetType.MAX_RATE; // Target MFD Constraint
		solutionName = "MaxRateSol_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}
	


	public void doUnconstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, double finalEnergy) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		// make the solution
		solutionType = InversionSolutionType.SIMULATED_ANNEALING;
		completionCriteria = new EnergyCompletionCriteria(finalEnergy);	// number of rows/subsections
		numSolutions = numSim; 
		solutionName = "UnconstrainedSA_"+numSolutions+"_finalE="+finalEnergy+"_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(rePlotOny);
	
	}

	
	
	public void doMFDconstrainedSA(boolean rePlotOny, SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, int numSim, MFD_TargetType mfdTargetType, double mfdWt, boolean asInitState, double finalEnergy) {
		this.setDefaultParameterValuess();
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
			initStateString = "initSt_";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		numSolutions = numSim; 
		solutionName = "MFDconstrSA_"+initStateString+numSolutions+"_finalE="+finalEnergy+"_"+mfdTargetType+"_wt"+Math.round(mfdWt)+"_"+slipRateProfile.toString()+"_"+
				slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(rePlotOny);
	
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
			ScalingRelationshipEnum scalingRel, int numSim, MFD_TargetType mfdTargetType, double totRateUncertFract, double totRateRelWt, boolean mfdAsInitState, double finalEnergy) {
		this.setDefaultParameterValuess();
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
			initStateString = "initSt_";
			initialStateType = SA_InitialStateType.FROM_MFD_CONSTRAINT;
		}
		numSolutions = numSim; 
		solutionName = "totRateConstrSA_"+initStateString+numSolutions+"_finalE="+finalEnergy+"_"+mfdTargetType+"_wt"+Math.round(totRateRelWt)+"_"+slipRateProfile.toString()+"_"+
				slipModelType.toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
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
			ScalingRelationshipEnum scalingRel, MFD_TargetType mfdTargetType) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		this.mfdTargetType = mfdTargetType; // Target MFD Constraint
		solutionType = InversionSolutionType.FRESH;
		solutionName = "FRESH_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
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
		solutionName = "SUNFiSH_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}

	
	/**
	 * 
	 * @param slipRateProfile
	 * @param slipModelType
	 * @param scalingRel
	 */
	public void doAppliedMFD_Solution(SlipRateProfileType slipRateProfile, SlipAlongRuptureModelEnum slipModelType, ScalingRelationshipEnum scalingRel,
			MFD_TargetType mfdTylpe) {
		this.setDefaultParameterValuess();
		this.slipRateProfile = slipRateProfile;
		this.slipModelType = slipModelType;
		this.scalingRel = scalingRel;
		solutionType = InversionSolutionType.RATES_FROM_MFD;
		mfdTargetType = mfdTylpe; // Target MFD Constraint
		solutionName = "AppliedMFD_"+mfdTylpe+"_"+slipRateProfile.toString()+"_"+slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	}



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
				slipModelType.toString()+"_".toString()+"_"+scalingRel.toString(); // Inversion name
		dirName = ROOT_PATH+"Output_"+solutionName;
		getSolution(false);
	
	}



	
	/**
	 * As written, the FaultSystemRuptureRateInversion constructor used here always sets the a-priori rupture rates and MFD constraint
	 * from a best-fitting GR distribution (whether these are applied depends on the weights given and the Simulated Annealing initial 
	 * solution)
	 * @param args
	 */
	public static void main(String []args) {
		
		SingleFaultInversion singleFaultInversion = new SingleFaultInversion();
		
//		singleFaultInversion.doNSHMP_GR_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);
		// this is qualitatively similar, but higher final energy
//		singleFaultInversion.doNSHMP_GR_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.HANKS_BAKUN_08);
		
		// this is qualitatively similar, but lower energy
//		singleFaultInversion.doNSHMP_GR_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B);
		
//		singleFaultInversion.doNSHMP_Char_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);
//		singleFaultInversion.doNSHMP_Char_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B);
		
//		singleFaultInversion.doMinRateSolution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);
//		singleFaultInversion.doMinRateSolution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B);
		
		//This misfits the very ends and a little elsewhere
//		singleFaultInversion.doMaxRateSolution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);
		// This fixes that (slip rates trimmed down at ends)
//		singleFaultInversion.doMaxRateSolution(SlipRateProfileType.UNIFORM_TRIMMED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);

		// this took six minutes (final energy E = 2)
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, 2);
		// this E=5 goes fast,and corresponds to RMS = sqrt(E/N) = sqrt(5/117) = 0.2 (ave Norm abs diff is less than this)
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, 5);
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5);
		// this still produces edge effects (only use this slip rate distribution for max rate model)
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.UNIFORM_TRIMMED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5);
		
		
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5);
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, 5);

		// this is to not over fit the data
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, 137);
		// on average this over-fits the data (also the error on the average is 94 versus 137, whereas it should be ~16 if it's not systematically over fitting)
//		singleFaultInversion.doUnconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 50, 137);

		// this should be the same as 3rd above (since wt=0 and zero initial state); looks good
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, false, 5);	
		// Use MFD as starting model only
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, true, 1);
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0, true, 1);	
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0, true, 1);	
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0, true, 1);	
	
		
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, false, 2);	
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1, true, 2);	
		
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, false, 2);	
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1, true, 2);	

		// MFD not perfectly fit in following; what happens if it's is (given high weight)?  This also took 3 hrs
//		singleFaultInversion.doMFDconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1, false, 2);	

		// NNLS Solutions
//		singleFaultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_minus1, 1);	
//		singleFaultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_0pt0, 1);	
//		singleFaultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0, 1);	
		// removing MFD produces the minimum rate model
//		singleFaultInversion.doMFDconstrainedNNLS(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0, 0);	

		// this took 118 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0.02, 1.0, true, 2);
		singleFaultInversion.doTotalRateconstrainedSA(true, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.02, 1.0, true, 2);
		
		// this took 16 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0.02, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.02, 1.0, false, 2);

		// this took 25 min
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, true, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, true, 2);
		// this took 2.5 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.02, 1.0, false, 2);
		// this took 2.5 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, true, 2);
		// this took <1.5 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.02, 1.0, false, 2);
		// this took 15 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 2);
		// this, as expected, has a 10e-22 initial energy (no iterations necessary)
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.UNIFORM_TRIMMED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.UNIFORM_TRIMMED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, true, 2);
		// this took 86 minutes
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MAX_RATE, 0.02, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MAX_RATE, 0.02, 1.0, false, 2);

		// this takes forever and never got below energy of ~10
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, true, 2);
		// this takes forever and never got below energy of ~10
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.MIN_RATE, 0.02, 1.0, false, 2);

		
		// NEED TO UNDERSTAND WHY THESE FRESH SOLUTION ARE NOT BETTER FITS
		// E=1170 for this:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);	
		// E = 543 for this:
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);
		// E=488 for this:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);	
		// E=131 for this:
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);
		// E=1348 for this:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);	
		// E=1336 for this (should be same as above):
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_1pt0);

		// this has E=86:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE);	
		// this has E=2424:
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MAX_RATE);

		// this has E=437:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MIN_RATE);	
		// this has E=435:
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.MIN_RATE);

		// this has E=880:
//		singleFaultInversion.doFRESH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.M7pt25only);	
		// this has E=456:
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.M7pt25only);

		// this has E=895
//		singleFaultInversion.doSUNFiSH_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);	
		// this has E=2231
//		singleFaultInversion.doSUNFiSH_Solution(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B);	

		
		
		
		
		
		
		
		
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 0.1, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.1, 1.0, false, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0.1, 1.0, true, 2);
		
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.1, 10.0, false, 5);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.1, 1.0, true, 2);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0.1, 10.0, true, 5);

//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0.1, 10.0, false, 5);
//		singleFaultInversion.doTotalRateconstrainedSA(false, SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0.1, 10.0, true, 5);

				
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1, true);	
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1, false);	

//		singleFaultInversion.doMFDconstrainedNNLS(SlipRateProfileType.UNIFORM, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_minus1, 1);	


//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1, true);	
//		singleFaultInversion.doAppliedMFD_Solution(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, MFD_TargetType.GR_b_minus1);

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 2, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.UNIFORM_TRIMMED, SlipAlongRuptureModelEnum.UNIFORM, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_0pt0, 1, false);
		
		

		
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1e4, false);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_1pt0, 1e6, false); // I also increased the number of iterations here
		
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1e4, false);

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1e4, false);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1e4, false);

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1e4, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 1e4, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 1e2, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 1, MFD_TargetType.GR_b_minus1, 0, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_minus1, 0, true);
		

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 1e4, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_0pt0, 0, true);

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 1e4, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.GR_b_1pt0, 0, true);

//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MIN_RATE, 0, true);
//		singleFaultInversion.doMFDconstrainedSA(SlipRateProfileType.TAPERED, SlipAlongRuptureModelEnum.TAPERED, ScalingRelationshipEnum.ELLSWORTH_B, 10, MFD_TargetType.MIN_RATE, 1e4, true);



		//DELETE THE FOLLOWING EVENTULLY)

//		// THE FOLLOWING SETS ALL THE INVERSION ATTRIBUTES:
//		//------------------------------------------------
//		
//		String dirName = ROOT_PATH+"Output";
//	
//		// Inversion name
//		String solutionName = "Single Fault Inversion";
//		
//		// Weighted inversion (whether to apply data weights)
//		boolean wtedInversion = true;
//		
//		// fault slip rates:
//		SlipRateProfileType slipRateProfile = SlipRateProfileType.TAPERED;
//		
//		// Average slip along rupture model
//		SlipAlongRuptureModelEnum slipModelType = SlipAlongRuptureModelEnum.UNIFORM;
//		
//		// Scaling Relationship
//		ScalingRelationshipEnum scalingRel = ScalingRelationshipEnum.ELLSWORTH_B;
//		
//		// Section Rate Constraints
//		double relativeSectRateWt=0;
//		
//		// Segmentation Constraints
////		String segmentationConstrFilename = ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE;
//		String segmentationConstrFilename = null;
//		double relative_segmentationConstrWt = 0;
//		
//		// A priori rupture rates
//		double relative_aPrioriRupWt = 0;
//		
//		// Minimum Rupture Rate
////		double minRupRate = 1e-8;
//		double minRupRate = 0.0;
//		
//		// Apply prob visible
//		boolean applyProbVisible = false;
//		
//		// Moment Rate Reduction
//		double moRateReduction = 0.0;	// this is the amount to reduce due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)
//		
//		// Target MFD Constraint
//		MFD_TargetType mfdTargetType = MFD_TargetType.GR_b_1pt0;
//		double relativeMFD_constraintWt = 0; // 
//		
//		// Total Rate Constraint
//		double totalRateConstraint = 0.005328;
//		double relativeTotalRateConstraintWt = 100/totalRateConstraint;
//
//		
//		// Inversion Solution Type:
//		
////		InversionSolutionType solutionType = InversionSolutionType.NSHMP_SOLUTION;
//
////		InversionSolutionType solutionType = InversionSolutionType.FROM_FILE;
//		String rupRatesFileDirName = ROOT_PATH+"OutputDataAndFigs_SA10_Segmented_Uniform_4/";	// don't comment out; this is only used if solutionType = InversionSolutionType.FROM_FILE
//
////		InversionSolutionType solutionType = InversionSolutionType.NON_NEGATIVE_LEAST_SQUARES;
//		
//		InversionSolutionType solutionType = InversionSolutionType.SIMULATED_ANNEALING;
//		int numSolutions = 1; // this is ignored for NON_NEGATIVE_LEAST_SQUARES which only has one possible solution
//		CoolingScheduleType saCooling = CoolingScheduleType.LINEAR;
//		GenerationFunctionType perturbationFunc = GenerationFunctionType.UNIFORM_NO_TEMP_DEPENDENCE;
//		// set the SA completion criteria (choose one of these):
//		long numIterations = (long) 1e5;
//		IterationCompletionCriteria completionCriteria = new IterationCompletionCriteria(numIterations);
////		double finalEnergy = 117d;	// this should be the number of rows to match data uncertinaties
////		EnergyCompletionCriteria completionCriteria = new EnergyCompletionCriteria(finalEnergy);
//		boolean initStateFromAprioriRupRates = true;
//		long randomSeed = System.currentTimeMillis();
////		long randomSeed = 1525892588112l; // for reproducibility; note that the last character here is the letter "l" to indicated a long value
//		
//		
//		// Data to plot and/or save:
//		boolean popUpPlots = true;	// this tells whether to show plots in a window (set null to turn off; e.g., for HPC)
//		boolean doDataFits=true;
//		boolean doMagHistograms=true;
//		boolean doNonZeroRateRups=true;
//	    boolean doSectPartMFDs=true;
//	    
//	    // to make hazard curve (set loc=null to ignore)
//		Location hazCurveLoc = null; 
//	    String hazCurveLocName = "Hazard Curve";
//	    
//	    // to make hazard maps
//	    boolean makeHazardMaps = false;
//		double saPeriodForHaz = 1.0;	// set as 0.0 for PGA
//		

		// THIS IS THE END OF THE INVERSION SETTINGS

		// Create an instance of this inversion class
//		SingleFaultInversion singleFaultInversion = new SingleFaultInversion();
		
//		singleFaultInversion.tempSetNSHMP_GR_Model(scalingRel);
		
//		// get the target MFD
//		IncrementalMagFreqDist targetMFD = singleFaultInversion.getTargetMFD(scalingRel, mfdTargetType);
//		if(D) System.out.println(mfdTargetType+" total rate "+targetMFD.getTotalIncrRate());
		
//		// This writes out a range of total rates
//		System.out.println("Rates for "+scalingRel);
//		double rate = singleFaultInversion.getTargetMFD(scalingRel, MFD_TargetType.GR_b_1pt0).getTotalIncrRate();
//		System.out.printf("\tGR b=1 Rate: %f\n",rate);
//		rate = singleFaultInversion.getTargetMFD(scalingRel, MFD_TargetType.GR_b_0pt0).getTotalIncrRate();
//		System.out.printf("\tGR b=0 Rate: %f\n",rate);
//		rate = singleFaultInversion.getTargetMFD(scalingRel, MFD_TargetType.BR_b_minus1).getTotalIncrRate();
//		System.out.printf("\tGR b=-1 Rate: %f\n",rate);
//		System.out.printf("\tMax Rate: %f\n",singleFaultInversion.getMaxPossibleRate(scalingRel));
//		System.out.printf("\tMin Rate: %f\n",singleFaultInversion.getMinPossibleRate(scalingRel));
//		System.exit(0);
//		// Result:
//		Rates for ELLSWORTH_B
//			GR b=1 Rate: 	0.088471
//			GR b=0 Rate: 	0.021428
//			GR b=-1 Rate: 	0.009476
//			Max Rate: 		0.842805
//			Min Rate: 		0.005328

		
//		singleFaultInversion.writeApriorRupRatesForMaxRateModel();


//		//Write section data (e.g., so it can be checked)
//		if(D) {
//			String str = "sectID\tslipRate\tslipRateStd\tlength\tdip\trake\tupDepth\tlowDepth\tDDW+\n";
//			for(int s=0;s<fltSectDataList.size();s++) {
//				FaultSectionPrefData fltSectData = fltSectDataList.get(s);
//				str += fltSectData.getSectionId()+"\t";
//				str += fltSectData.getOrigAveSlipRate()+"\t";
//				str += fltSectData.getOrigSlipRateStdDev()+"\t";
//				str += (float)fltSectData.getTraceLength()+"\t";
//				str += fltSectData.getAveDip()+"\t";
//				str += fltSectData.getAveRake()+"\t";
//				str += fltSectData.getOrigAveUpperDepth()+"\t";
//				str += fltSectData.getAveLowerDepth()+"\t";
//				str += fltSectData.getOrigDownDipWidth()+"\t";
//				str += (float)fltSectData.getStirlingGriddedSurface(1.0).getArea()+"\n";
//			}
//			System.out.println(str);			
//		}
		
//	    // THE FOLLOWING IS TO MAKE HAZARD MAP RATIOS FOR DATA IN FILES
//	    String fileName1 = ROOT_PATH+"OutputDataAndFigs_SA10_Unsegmented_Uniform"+"/hazardMaps/1.0secSA_2in50.txt";
//	    String fileName2 = ROOT_PATH+"OutputDataAndFigs_SA10_Segmented_Uniform"+"/hazardMaps/1.0secSA_2in50.txt";
//	    String label = "1.0secSA_2in50_2in50_ratio";
//	    wasatchInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName+"/hazardMaps", true);
//		System.exit(-1);
//		


//		System.exit(0);
		
//		// this will be used to keep track of runtimes
//		long startTimeMillis = System.currentTimeMillis();
//		if(D)
//			System.out.println("Starting Inversion");
//		
//		// create an instance of the inversion class with the above settings
//		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
//				solutionName,
//				"default",
//				fltSectDataList, 
//				sectionRateConstraints, 
//				rupSectionMatrix, 
//				slipModelType, 
//				scalingRel, 
//				relativeSectRateWt, 
//				relative_aPrioriRupWt, 
//				ROOT_DATA_DIR+APRIORI_RUP_RATE_FILENAME,
//				wtedInversion, 
//				minRupRate, 
//				applyProbVisible, 
//				moRateReduction,
//				targetMFD,
//				relativeMFD_constraintWt,
//				segmentationConstrFilename,
//				relative_segmentationConstrWt,
//				totalRateConstraint,
//				relativeTotalRateConstraintWt);
//
//		
//		// make the directory for storing results
//	    File file = new File(dirName);
//	    file.mkdirs();
//	    
//	    // set a-prior rates from MFD so these can be applied as initial model
//	    fltSysRupInversion.setAprioriRupRatesFromMFD_Constrint();
//	    
//	    // write the setup info to a file
//	    fltSysRupInversion.writeInversionSetUpInfoToFile(dirName);
//	    
//	    switch (solutionType) {
//	    		case NON_NEGATIVE_LEAST_SQUARES:
//	    			if(D) System.out.println("NON_NEGATIVE_LEAST_SQUARES");
//	    			fltSysRupInversion.doInversionNNLS();
//	    			break;
//	    		case NSHMP_SOLUTION:
//	    			if(D) System.out.println("NSHMP_SOLUTION");
//	    			fltSysRupInversion.doNSHMP_Solution(targetMFD);
//	    			break;
//	    		case SIMULATED_ANNEALING:
//	    			if(numSolutions==1) {
//		    			if(D) System.out.println("SIMULATED_ANNEALING; numSolutions=1");
//	    				fltSysRupInversion.doInversionSA(completionCriteria, initStateFromAprioriRupRates, randomSeed, saCooling, perturbationFunc);
//	    			}
//	    			else if(numSolutions>1) {
//		    			if(D) System.out.println("SIMULATED_ANNEALING; numSolutions="+numSolutions);
//	    				fltSysRupInversion.doInversionSA_MultTimes(completionCriteria, initStateFromAprioriRupRates, randomSeed, numSolutions, dirName, saCooling,perturbationFunc);
//	    			}
//	    			else {
//	    				throw new RuntimeException("bad numIterations");
//	    			}
//	    			break;
//	    		case FROM_FILE:
//			    	if(numSolutions==1) {
//		    			if(D) System.out.println("FROM_FILE; numSolutions=1");
//		    			String rupRatesFileName = rupRatesFileDirName+ "ruptureRatesAlt.txt";
//			    		double[] rupRatesArray = singleFaultInversion.readRuptureRatesFromFile(rupRatesFileName);
//			    		fltSysRupInversion.setSolution(rupRatesArray, "Solution from file: "+rupRatesFileName);
//			    	}
//			    	else if(numSolutions>1) {
//		    			if(D) System.out.println("FROM_FILE; numSolutions="+numSolutions);
//		    			ArrayList<double[]> rupRatesArrayList = new ArrayList<double[]>();
//		    			for(int i=0;i<numSolutions;i++) {
//			    			String rupRatesFileName = rupRatesFileDirName+ "ruptureRates_"+i+".txt";
//			    			rupRatesArrayList.add(singleFaultInversion.readRuptureRatesFromFile(rupRatesFileName));
//			    			fltSysRupInversion.setMultipleSolutions(rupRatesArrayList, "Multiple solutions read from files with prefix "+rupRatesFileName, dirName);
//		    			}
//			    	}
//			    	else {
//			    		throw new RuntimeException("bad numIterations");
//			    	}
//			    	break;
//	    }
//
//		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
//		if(D) System.out.println("Done with Inversion after "+(float)runTimeSec+" seconds.");
//				
//		// write results to file
//		fltSysRupInversion.writeInversionRunInfoToFile(dirName);
//			
//		if(doDataFits) fltSysRupInversion.writeAndOrPlotDataFits(dirName, popUpPlots);
//		if(doMagHistograms) fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots);
//		if(doNonZeroRateRups) fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots);
//		if(doSectPartMFDs) fltSysRupInversion.writeAndOrPlotSectPartMFDs(dirName, popUpPlots);
//	    
//	    // hazard curve:
//		if(hazCurveLoc != null) {
//			singleFaultInversion.writeAndOrPlotHazardCurve(fltSysRupInversion, hazCurveLoc, saPeriodForHaz, dirName, popUpPlots, hazCurveLocName);
//		}
//	    
//	    // second parameter here is SA period; set as 0.0 for PGA:
//		if(makeHazardMaps)
//			singleFaultInversion.makeHazardMaps(fltSysRupInversion, saPeriodForHaz, dirName, popUpPlots);
	}
}
