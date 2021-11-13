/**
 * 
 */
package scratch.ned.FSS_Inversion2019;


import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.nnls.NNLSWrapper;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotColorAndLineTypeSelectorControlPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.RunScript;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ThreadedSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.EnergyCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.CoolingScheduleType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.finalReferenceFaultParamDb.DeformationModelPrefDataFinal;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import scratch.UCERF3.U3FaultSystemRupSet;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.ScalingRelationshipEnum;
import scratch.ned.FSS_Inversion2019.logicTreeEnums.SlipAlongRuptureModelEnum;

/**
 * 
 * This class does an inversion for the rate of events in an unsegmented fault model:
 * 
 * TO DO:
 * 
 * 1) Improve the segmentation constraint implementation?
 * 2) Input prob visible model rather than computing here  (already in U3?)
 * 3) sample MRIs via monte carlo simulations (same for slip rates?) for more epistemic 
 *    uncertainty (or do this outside with zero errors); or just use SA with fewer iterations?
 *
 */
public class FaultSystemRuptureRateInversion {
	
	final static boolean D = true;	// debugging flag
	
	public CoolingScheduleType sa_coolingSchedule = CoolingScheduleType.VERYFAST_SA;	// this is the default

	double GAUSS_MFD_SIGMA = 0.0; // typically otherwise 0.12;
	public final static double GAUSS_MFD_TRUNCATION = 2;
	public final static double MAG_DELTA = 0.1;
	
	public final static double MIN_RUP_RATE = 1e-9;
	
	String modelName;
	
	// Strings to keep track of results
	String modelSetUpInfoString, modelRunInfoString;
	
	// input data
	private ArrayList<FaultSectionPrefData> fltSectionDataList;
	private ArrayList<SectionRateConstraint> sectionRateConstraints; // using old class for "segments"
	private int[][] rupSectionMatrix;
	
	// section attributes
	private int numSections;
	private double[] sectArea, sectLength, sectSlipRate, sectSlipRateStdDev, sectMoRate, sectRake;

	// rupture attributes
	int numRuptures;
	private String[] rupNameShort;
	private double[] rupLength, rupArea, rupMeanMag, rupMeanMo, rupMoRateMax, rupAveSlip, rupAveRake; // rupMoRateMax is the moment rate available if the rupture is the only one on the sections
	private double minRupMag, maxRupMag,minRupMagWithAleatory, maxRupMagWithAleatory;
	double[] rupRateSolution; // these are the rates from the inversion (not total rate of MFD)
	
	// these are for MFD plots (rounded to nearest 0.1 mag unit)
	public double minMagMFD, maxMagMFD,minMagMFD_WithAleatory, maxMagMFD_WithAleatory;


	// section-rupture attributes
	private int[] numSectInRup, firstSectOfRup;
	private int minNumSectInRup;  // this sets the number of sections for the smallest ruptures
	private double[][] sectSlipInRup;
	private double[] rateOfThroughGoingRupsAtSectBoudary;

	double totMoRate;
	
	IncrementalMagFreqDist mfdConstraint, mfdSigma; 

	private SummedMagFreqDist aveOfSectPartMFDs, rupSamplerMFD;
	
	
	// this accounts for fact that ave slip from Gauss MFD is greater than the slip of the average mag
	private double gaussMFD_slipCorr;
	
	// a-priori rate constraints
	int num_aPriori_constraints;
	int[] aPriori_rupIndex;
	double[] aPriori_rate, aPriori_wt;
	String aPrioriRupRatesFilename;
	
	// segmentation (rup rate = 0) constraints
	int num_segConstraints;
	ArrayList<SegmentationConstraint> segmentationConstraintList;
	ArrayList<SlipRateSegmentationConstraint> slipRateSegmentationConstraintList;
	
	// the following specifies fault section indices where a Laplacian smoothness constraint should be applied
	ArrayList<int[]> smoothnessConstraintList;


	
	private static boolean MATLAB_TEST = false;
	double[][] C_wted, C;	// inversion matrices
	double[] d, d_wted, data_wt, full_wt, d_pred;  // the data vector
	
	private double minRupRateArray[]; // the minimum rate constraint for each rupture
	private double minRupRate;
	private boolean wtedInversion, applyProbVisible;	// weight the inversion according to slip rate and segment rate uncertainties
	private double relativeSectRateWt, relative_aPrioriRupWt, relative_segConstraintWt; 
	private double relativeMFD_constraintWt;
	private double totalRateConstraint, totalRateSigma, relativeTotalRateConstraintWt;
	private double relativeSmoothnessConstraintWt;

	private int  firstRowSectSlipRateData=-1,firstRowSectEventRateData=-1, firstRowAprioriData=-1, firstRowSegConstraint=-1;
	private int  lastRowSectSlipRateData=-1,lastRowSectEventRateData=-1, lastRowAprioriData=-1, lastRowSegConstraint=-1;
	private int firstRowMFD_constraintData=-1, lastRowMFD_constraintData=-1;
	private int firstRowTotalRateConstraint=-1, lastRowTotalRateConstraint=-1;
	private int firstRowSmoothnessConstraint=-1, lastRowSmoothnessConstraint=-1;
	private int totNumRows;
	
	// average slip along rupture model
	private SlipAlongRuptureModelEnum slipModelType;

	private static HistogramFunction taperedSlipPDF, taperedSlipCDF;
	
	// MFDs
	SummedMagFreqDist magFreqDist, meanMagHistorgram, magHistorgram;
	ArrayList<SummedMagFreqDist> sectNucleationMFDs;
	ArrayList<SummedMagFreqDist> sectParticipationMFDs;
	
	
	
	private double[] finalSectEventRate, finalPaleoVisibleSectEventRate, finalSectSlipRate, finalSectMeanSlip, finalSectSlipCOV;
	
	// the following is the total moment-rate reduction, including that which goes to the  
	// background, afterslip, events smaller than the min mag here, and aftershocks and foreshocks.
	private double moRateReduction;  
	
	private ScalingRelationshipEnum magAreaRel;
	
	// NNLS inversion solver - static to save time and memory
	private static NNLSWrapper nnls = new NNLSWrapper();

	
	// These contain data from multiple SA runs:
	ArrayList<double[]> rupRatesFromMultRunsArrayList;
	ArbDiscrEmpiricalDistFunc_3D rupRatesFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D mfdsFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D cumMfdsFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D finalSectSlipRateFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D finalPaleoVisibleSectEventRateFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D rateOfThroughGoingRupsAtSectBoudaryFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D finalSectMeanSlipFromMultRuns;
	ArbDiscrEmpiricalDistFunc_3D finalSectSlipCOV_FromMultRuns;
	

	/**
	 * This returns the number of solutions
	 * @return
	 */
	public int getNumSolutions() {
		if(rupRatesFromMultRunsArrayList == null)
			return 1;
		else
			return rupRatesFromMultRunsArrayList.size();
	}
	
	public int getNumRuptures() {
		return numRuptures;
	}
	
	
	
	/**
	 * This writes the 2D segPartMFDs to a file and/or makes a plot of the
	 * results. 
	 * @param dirName - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 */
	public ArrayList<XYZPlotSpec> writeAndOrPlotSectPartMFDs(String dirName, boolean popUpWindow, double widthInches, double heightInches,
			Range xRange, Range yRange, Range zRange) {
		
		ArrayList<XYZPlotSpec> specList = new ArrayList<XYZPlotSpec>();
		
		// this writes out 
		try{
			
			String fileName = dirName+"/sectPartMFDsData";
			String fileNameCum = dirName+"/sectPartCumMFDsData";
			FileWriter fw=null, cfw=null;
			if(dirName != null) {
				fw = new FileWriter(fileName+".txt");
				cfw = new FileWriter(fileNameCum+".txt");
				fw.write("seg_index\tmag\tpart_rate\n");
				cfw.write("seg_index\tmag\tpart_cum_rate\n");				
			}
			SummedMagFreqDist mfd = sectParticipationMFDs.get(0);
			EvenlyDiscretizedFunc cmfd = mfd.getCumRateDist();
			EvenlyDiscrXYZ_DataSet xyzDataSectPart = new EvenlyDiscrXYZ_DataSet(sectParticipationMFDs.size(), mfd.size(), 0, mfd.getMinX(), 1.0 , mfd.getDelta());
			EvenlyDiscrXYZ_DataSet xyzDataSectPartCum = new EvenlyDiscrXYZ_DataSet(sectParticipationMFDs.size(), cmfd.size(), 0, cmfd.getMinX(), 1, cmfd.getDelta());
			for(int s=0; s<sectParticipationMFDs.size(); s++) {
				mfd = sectParticipationMFDs.get(s);
				cmfd = mfd.getCumRateDist();
				for(int m=0; m<mfd.size();m++) {
					if( mfd.getY(m) != 0.0) {
						if(dirName != null)
							fw.write(s+"\t"+(float)mfd.getX(m)+"\t"+(float)Math.log10(mfd.getY(m))+"\n");
						xyzDataSectPart.set(s, mfd.getX(m), Math.log10(mfd.getY(m)));
					}
					else {
						xyzDataSectPart.set(s, mfd.getX(m), Double.NaN);
					}

					if( cmfd.getY(m) != 0.0) {
						if(dirName != null)
							cfw.write(s+"\t"+(float)cmfd.getX(m)+"\t"+(float)Math.log10(cmfd.getY(m))+"\n");
						xyzDataSectPartCum.set(s, cmfd.getX(m), Math.log10(cmfd.getY(m)));
					}
					else {
						xyzDataSectPartCum.set(s, cmfd.getX(m), Double.NaN);
					}
				}
			}
			
			fw.close();
			cfw.close();
			
			specList.add(PlottingUtils.make2D_plot(xyzDataSectPart, null, "Subsection", "Magnitude", "Incremental Rate", 
					fileName, popUpWindow, xRange, yRange, zRange, widthInches, heightInches));
			specList.add(PlottingUtils.make2D_plot(xyzDataSectPartCum, null, "Subsection", "Magnitude", "Cumulative Rate", 
					fileNameCum, popUpWindow, xRange, yRange, zRange, widthInches, heightInches));		

		}catch(Exception e) {
			e.printStackTrace();
		}
		
		return specList;
	}

	
	
	
	
	
	/**
	 * This writes and/or plots the rates and slip-rate contribution of each rupture that has a 
	 * final rate above the minimum.  
	 * @param dirName - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 */
	public void writeAndOrPlotNonZeroRateRups(String dirName, boolean popUpWindow, double widthInches, double heightInches) {
				
		// get number of ruptures above minRupRate
		int numAboveMinRate = 0;
		for(int rup=0; rup<numRuptures;rup++)
			if(rupRateSolution[rup]>minRupRate)
				numAboveMinRate += 1;
		
		if(numAboveMinRate==1) {
			System.out.println("numAboveMinRate==1, so aborting writeAndOrPlotNonZeroRateRups");
			return;
		}
		
		EvenlyDiscrXYZ_DataSet xyzDataRupSlipRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);
		EvenlyDiscrXYZ_DataSet xyzDataRupRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);
		
		int index = 0;
		try{
			FileWriter fw=null,fw2=null;
			String fileName_sr=null, fileName_rr=null;
			if(dirName != null) {
				fileName_sr = dirName+"/rupSlipRateOnSectData";
				fileName_rr = dirName+"/rupRateOnSectData";
				fw = new FileWriter(fileName_sr+".txt");
				fw2 = new FileWriter(fileName_rr+".txt");
				fw.write("rup_index\tseg_index\tslip_rate\n");
				fw2.write("rup_index\tseg_index\trate\n");				
			}
			for(int rup=0; rup<numRuptures;rup++) {
				if(rupRateSolution[rup]>minRupRate) {
					double logRate = Math.log10(rupRateSolution[rup]);
					for(int s=0; s<numSections; s++) {
						if(rupSectionMatrix[s][rup] == 1) {
							double logSlipRate = Math.log10(sectSlipInRup[s][rup]*rupRateSolution[rup]);
							if(dirName != null) {
								fw.write(index+"\t"+s+"\t"+(float)logSlipRate+"\n");
								fw2.write(index+"\t"+s+"\t"+(float)logRate+"\n");
							}
							xyzDataRupSlipRate.set(s, index,logSlipRate);
							xyzDataRupRate.set(s, index,logRate);
						}
						else {
							xyzDataRupSlipRate.set(s, index, Double.NaN);
							xyzDataRupRate.set(s, index, Double.NaN);		
						}
					}
					index+= 1;
				}
			}	

			fw.close();
			fw2.close();
			
			Range xAxisRange = new Range(0-0.5, numSections+0.5);
			Range yAxisRange = new Range(xyzDataRupSlipRate.getMinY()-0.5, xyzDataRupSlipRate.getMaxY()+0.5);

			PlottingUtils.make2D_plot(xyzDataRupSlipRate, "Slip Rate from Ruptures", "Section", "Rupture", "Log10 Slip Rate", 
					fileName_sr, popUpWindow, xAxisRange, yAxisRange, null, widthInches, heightInches);
			PlottingUtils.make2D_plot(xyzDataRupRate, "Rate of Ruptures", "Section", "Rupture", "Log10 Rate", 
					fileName_rr, popUpWindow, xAxisRange, yAxisRange, null, widthInches, heightInches);


		}catch(Exception e) {
			e.printStackTrace();
		}			
	}
	
	
	
	/**
	 *   
	 * @param dirName - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 */
	public void writeAndOrPlotRupRateOfFirstSection(String dirName, boolean popUpWindow, double widthInches, double heightInches, Range zRange) {
		
		EvenlyDiscrXYZ_DataSet xyzDataRupRate = new EvenlyDiscrXYZ_DataSet(numSections,numSections+1, 0, 0, 1, 1);
		for(int x=0; x<numSections;x++)
			for(int y=0; y<numSections+1;y++)
				xyzDataRupRate.set(x, y, Double.NaN);
		
		try{
			FileWriter fw=null;
			String fileName_rr=null;
			if(dirName != null) {
				fileName_rr = dirName+"/rupRateOnFirstSectData";
				fw = new FileWriter(fileName_rr+".txt");
				fw.write("firstSect\tnumSectInRup\tlog10_rate\n");				
			}
			for(int rup=0; rup<numRuptures;rup++) {
				int numSect = numSectInRup[rup];
					double logRate = Double.NaN;
					if(rupRateSolution[rup]>0)
						logRate = Math.log10(rupRateSolution[rup]);
					for(int s=0; s<numSections; s++) {
						if(rupSectionMatrix[s][rup] == 1) {
							if(dirName != null) {
								fw.write(s+"\t"+numSect+"\t"+(float)logRate+"\n");
							}
							xyzDataRupRate.set(s, numSect,logRate);
							break;
						}
					}
			}	

			Range rangeX = new Range(-0.5,numSections+0.5);
			Range rangeY = new Range(-0.5,numSections+1.5);
			fw.close();
			PlottingUtils.make2D_plot(xyzDataRupRate, null, "First Subsection", "Num Subsections In Rupture", "Log10 Rupture Rate", 
					fileName_rr, popUpWindow, rangeX, rangeY, zRange, widthInches, heightInches);

		}catch(Exception e) {
			e.printStackTrace();
		}			
	}


	/**
	 * Write the inversion-run info to a file
	 * @param dirName
	 */
	public void writeInversionRunInfoToFile(String dirName) {
		FileWriter fw;
		try {
			fw = new FileWriter(dirName+"/InversionRunInfo.txt");
			fw.write(modelRunInfoString);
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Write the inversion-run info to a file
	 * @param dirName
	 */
	public void writeInversionSetUpInfoToFile(String dirName) {
		FileWriter fw;
		try {
			fw = new FileWriter(dirName+"/InversionSetUpInfo.txt");
			fw.write(modelSetUpInfoString);
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}




	
	/**
	 * @param modelName - any name the user wants to give the inversion model
	 * @param slipRateModelName - the name of the slip rate model (only used for metadata)
	 * @param fltSectionDataList
	 * @param sectionRateConstraints
	 * @param rupSectionMatrix
	 * @param slipModelType
	 * @param scalingRel
	 * @param relativeSectRateWt - weight on segment rate equations (relative to slip rate)
	 * @param relative_aPrioriRupWt - weight on a-priori rates (relative to slip rate)
	 * @param aPrioriRupRatesFilename - filename for aPriori rupture rates
	 * @param wtedInversion - apply data uncertainties?
	 * @param minRupRate - constrain all rupture rates to be greater than this value
	 * @param applyProbVisible - account for likelihood that Paleoseismology will see the rupture
	 * @param moRateReduction - fraction reduction from smaller events (and not aseismicity or coupling coefficient, which are set in fltSectionDataList)
	 * @param mfdConstraint - an IncrementalMagFreqDist constraint
	 * @param mfdSigma - the 1 standard deviation uncertainty on mfdConstraint
	 * @param relativeMFD_constraintWt - weight for MFD constraint
	 * @param segConstraintFilename - each line has sect1_ID, sect2_ID, and wt (rups that utilize these sections are given a prior rate of zero)
	 * @param relative_segConstraintWt - 
	 * @param totalRateConstraint
	 * @param totalRateSigma
	 * @param relativeTotalRateConstraintWt,
	 * @param smoothnessConstraintList - ArrayList<int[]>, Laplacian smoothing applied to each list of fault sections
	 * @param relativeSmoothnessConstraintWt,
	 * @param magAareaAleatoryVariability - 
	 */
	public FaultSystemRuptureRateInversion( String modelName,
			String slipRateModelName,
			ArrayList<FaultSectionPrefData> fltSectionDataList, 
			int[][] rupSectionMatrix, 
			SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, 
			ArrayList<SlipRateSegmentationConstraint> slipRateSegmentationConstraintList,
			ArrayList<SectionRateConstraint> sectionRateConstraints,
			double relativeSectRateWt, 
			double relative_aPrioriRupWt, 
			String aPrioriRupRatesFilename,
			boolean wtedInversion, 
			double minRupRate,
			boolean applyProbVisible, 
			double moRateReduction, 
			IncrementalMagFreqDist mfdConstraint,
			IncrementalMagFreqDist mfdSigma, // uncertainty of the MFD (1 sigma)
			double relativeMFD_constraintWt,
			ArrayList<SegmentationConstraint> segmentationConstraintList,
			double relative_segConstraintWt,
			double totalRateConstraint,
			double totalRateSigma,
			double relativeTotalRateConstraintWt,
			ArrayList<int[]> smoothnessConstraintList,
			double relativeSmoothnessConstraintWt,
			double magAareaAleatoryVariability) {
		
		this.modelName = modelName;
		this.fltSectionDataList = fltSectionDataList;
		this.sectionRateConstraints = sectionRateConstraints;
		this.rupSectionMatrix = rupSectionMatrix;
		this.slipModelType = slipModelType;
		this.magAreaRel = scalingRel;
		this.slipRateSegmentationConstraintList = slipRateSegmentationConstraintList;
		this.relativeSectRateWt = relativeSectRateWt;
		this.relative_aPrioriRupWt = relative_aPrioriRupWt;
		this.aPrioriRupRatesFilename = aPrioriRupRatesFilename;
		this.wtedInversion = wtedInversion;
		this.minRupRate = minRupRate;
		this.applyProbVisible = applyProbVisible;
		this.moRateReduction = moRateReduction;
		this.mfdConstraint = mfdConstraint;
		this.mfdSigma = mfdSigma;
		this.relativeMFD_constraintWt = relativeMFD_constraintWt;
		this.segmentationConstraintList = segmentationConstraintList;
		this.relative_segConstraintWt = relative_segConstraintWt;
		this.totalRateConstraint = totalRateConstraint;
		this.totalRateSigma = totalRateSigma;
		this.relativeTotalRateConstraintWt = relativeTotalRateConstraintWt;
		this.smoothnessConstraintList = smoothnessConstraintList;
		this.relativeSmoothnessConstraintWt = relativeSmoothnessConstraintWt;
		this.GAUSS_MFD_SIGMA = magAareaAleatoryVariability;

		
		if(modelName == null) {
			modelName = this.getClass().toString();
		}

		// set info string
		modelSetUpInfoString = modelName+" with:\n\n";
		modelSetUpInfoString += "\tslipRateModelName = "+slipRateModelName+"\n";
		modelSetUpInfoString += "\tslipModelType = "+slipModelType+"\n";
		modelSetUpInfoString += "\tmagAreaRel = "+magAreaRel+"\n";
		modelSetUpInfoString += "\trelativeSectRateWt = "+relativeSectRateWt+"\n";
		modelSetUpInfoString += "\trelative_aPrioriRupWt = "+relative_aPrioriRupWt+"\n";
		modelSetUpInfoString += "\taPrioriRupRatesFilename = "+aPrioriRupRatesFilename+"\n";
		modelSetUpInfoString += "\twtedInversion = "+wtedInversion+"\n";
		modelSetUpInfoString += "\tminRupRate = "+minRupRate+"\n";
		modelSetUpInfoString += "\tapplyProbVisible = "+applyProbVisible+"\n";
		modelSetUpInfoString += "\tmoRateReduction = "+moRateReduction+"\n";
		modelSetUpInfoString += "\trelativeMFD_constraintWt = "+relativeMFD_constraintWt+"\n";
		if(segmentationConstraintList != null)
			modelSetUpInfoString += "\tsegmentationConstraintList.size() = "+segmentationConstraintList.size()+"\n";
		if(slipRateSegmentationConstraintList != null)
			modelSetUpInfoString += "\tslipRateSegmentationConstraintList.size() = "+slipRateSegmentationConstraintList.size()+"\n";
		modelSetUpInfoString += "\trelative_segConstraintWt = "+relative_segConstraintWt+"\n";
		modelSetUpInfoString += "\ttotalRateConstraint = "+totalRateConstraint+"\n";
		modelSetUpInfoString += "\ttotalRateSigma = "+totalRateSigma+"\n";
		modelSetUpInfoString += "\trelativeTotalRateConstraintWt = "+relativeTotalRateConstraintWt+"\n";
		modelSetUpInfoString += "\tGAUSS_MFD_SIGMA = "+GAUSS_MFD_SIGMA+"\n";
		if(smoothnessConstraintList != null)
			modelSetUpInfoString += "\tsmoothnessConstraintList.size() = "+smoothnessConstraintList.size()+"\n";
		modelSetUpInfoString += "\tGAUSS_MFD_TRUNCATION = "+GAUSS_MFD_TRUNCATION+"\n";
		
		
		// apply slip rate segmentation constraints
		if(slipRateSegmentationConstraintList != null) {
			for(SlipRateSegmentationConstraint srConstr: slipRateSegmentationConstraintList) {
				for(FaultSectionPrefData fltSect:fltSectionDataList)
					if(srConstr.getSectIndex() == fltSect.getSectionId()) {
						fltSect.setAveSlipRate(fltSect.getOrigAveSlipRate()*srConstr.getSlipRateReductionFactor());
						fltSect.setSlipRateStdDev(fltSect.getOrigSlipRateStdDev()*srConstr.getSlipRateStdevReductionFactor());
					}
			}
		}

		
		// initialize section and rupture attributes
		initSectAndRupAttributes();
		
		
		// compute matrix of Dsr (slip on each segment in each rupture)
		computeSectSlipInRupMatrix();
		
		// read the a-priori rup rates & wts
		num_aPriori_constraints = 0;
		if(aPrioriRupRatesFilename != null) {
			readApriorRupRateConstraintsFromFile(aPrioriRupRatesFilename);
			num_aPriori_constraints = aPriori_rupIndex.length;
		}
		
		num_segConstraints = 0;
		if(relative_segConstraintWt>0 && segmentationConstraintList !=null ) {
			num_segConstraints = segmentationConstraintList.size();
		}
		
		initMatricesEtc();
	}
	
	
	private void initMatricesEtc() {
		
		// write out rupture attributes
		if(D) {
			System.out.println("Rupture attributes:");
			System.out.println("index\tmeanMag\tlength\tarea\taveSlip\tname");
			for(int r=0;r<numRuptures;r++) {
				double aveSlip = Math.round(rupAveSlip[r]*100.0)/100.0;	// cm precision
				System.out.println(r+"\t"+rupMeanMag[r]+"\t"+Math.round(rupLength[r]*1e-3)+"\t"+
						Math.round(rupArea[r]*1e-6)+"\t"+(float)aveSlip+"\t"+rupNameShort[r]);
			}
		}
					
		// get the number of section rate constraints
		int numSectRateConstraints = 0;
		if(sectionRateConstraints != null)
			numSectRateConstraints = sectionRateConstraints.size();
		
		// set the minimum rupture rate constraints
		if(minRupRate >0.0) {
			minRupRateArray = new double[numRuptures]; // this sets them all to zero
			for(int rup=0; rup<numRuptures; rup++) 
				minRupRateArray[rup] = minRupRate;			
		}
		
		// SET NUMBER OF ROWS AND IMPORTANT INDICES	
		// all section slip-rates used
		firstRowSectSlipRateData = 0;
		totNumRows = numSections;	// assumes all sections have slip-rates
		lastRowSectSlipRateData = totNumRows-1;
		
		// add section rate constrains if needed
		if(relativeSectRateWt > 0.0) {
			firstRowSectEventRateData = totNumRows;
			totNumRows += numSectRateConstraints;
			lastRowSectEventRateData = totNumRows-1;
		}
		
		// add a-priori rate constrains if needed
		if(relative_aPrioriRupWt > 0.0) {
			firstRowAprioriData  = totNumRows;
			totNumRows += num_aPriori_constraints;
			lastRowAprioriData = totNumRows-1;
		}
		
		// add segmentation constrains if needed
		if(relative_segConstraintWt > 0.0) {
			firstRowSegConstraint  = totNumRows;
			totNumRows += num_segConstraints;
			lastRowSegConstraint = totNumRows-1;
		}

		// add number of MFD constaints
		int numMFD_constraints=0;
		if(relativeMFD_constraintWt>0) {
			numMFD_constraints = mfdConstraint.size();
			firstRowMFD_constraintData = totNumRows;
			totNumRows += numMFD_constraints;
			lastRowMFD_constraintData = totNumRows-1;
		}
		
		// add total rate constraint
		int numTotalRateConstraints=0;
		if(relativeTotalRateConstraintWt>0) {
			numTotalRateConstraints = 1;
			firstRowTotalRateConstraint = totNumRows;
			totNumRows += numTotalRateConstraints;
			lastRowTotalRateConstraint = totNumRows-1;
		}
		// smoothness constraints
		int numSmoothnessConstraints=0;
		if(relativeSmoothnessConstraintWt>0 && smoothnessConstraintList != null) {
			for(int[] sectArray : smoothnessConstraintList)
				numSmoothnessConstraints += sectArray.length-2;
			firstRowSmoothnessConstraint = totNumRows;
			totNumRows += numSmoothnessConstraints;
			lastRowSmoothnessConstraint = totNumRows-1;
		}

				
		String tempString = "firstRowSegEventRateData = "+firstRowSectEventRateData+
				"; firstRowAprioriData = "+firstRowAprioriData+
				"; firstRowSegConstraint = "+firstRowSegConstraint+
				"; firstRowMFD_constraintData = "+firstRowMFD_constraintData+
				"; firstRowTotalRateConstraint = "+firstRowTotalRateConstraint+
				"; firstRowSmoothnessConstraint = "+firstRowSmoothnessConstraint+
				"; totNumRows = "+totNumRows;
		if(D) System.out.println("\n"+tempString+"\n");
		modelSetUpInfoString += "\n"+tempString+"\n";
			
		C = new double[totNumRows][numRuptures];
		d = new double[totNumRows];  // data vector
		C_wted = new double[totNumRows][numRuptures]; // wted version
		d_wted = new double[totNumRows];  // wted data vector

		data_wt = new double[totNumRows];  // data weights
		full_wt = new double[totNumRows];  // data weights
		d_pred = new double[totNumRows];  // predicted data vector
		
		// initialize wts to 1.0
		for(int i=0;i<data_wt.length;i++)  data_wt[i]=1.0;
			
		// CREATE THE MODEL AND DATA MATRICES
		
		// first fill in the slip-rate constraints & wts
		for(int row = firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row ++) {
			d[row] = sectSlipRate[row]*(1-moRateReduction);
			if(wtedInversion)
				data_wt[row] = 1/((1-moRateReduction)*sectSlipRateStdDev[row]);
			for(int col=0; col<numRuptures; col++)
				C[row][col] = sectSlipInRup[row][col];
		}
		
		// now fill in the section event rate constraints if requested
		if(relativeSectRateWt > 0.0) {
			SectionRateConstraint constraint;
			for(int i = 0; i < numSectRateConstraints; i ++) {
				constraint = sectionRateConstraints.get(i);
				int sect = constraint.getSectIndex();
				int row = i+firstRowSectEventRateData;
				d[row] = constraint.getMeanRate(); // this is the average sub-section rate
				if(wtedInversion)
					data_wt[row] = 1/constraint.getStdDevOfMean();
				for(int col=0; col<numRuptures; col++)
					if(applyProbVisible)
						C[row][col] = rupSectionMatrix[sect][col]*getProbVisible(rupMeanMag[col]);
					else
						C[row][col] = rupSectionMatrix[sect][col];
			}
		}
		
		// now fill in the a-priori rates if needed
		if(relative_aPrioriRupWt > 0.0) {
			for(int i=0; i < num_aPriori_constraints; i++) {
				int row = i+firstRowAprioriData;
				int col = aPriori_rupIndex[i];
				d[row] = aPriori_rate[i];
				if(wtedInversion)
					data_wt[row] = aPriori_wt[i];
				C[row][col]=1.0;
			}
		}
		
		
		// now fill in the segmentation constraint if needed
		if(relative_segConstraintWt > 0.0) {
			for(int i=0; i < num_segConstraints; i++) {
				SegmentationConstraint segConst = this.segmentationConstraintList.get(i);
				int sect1 = segConst.getSect1_Index();
				int sect2 = segConst.getSect2_Index();
				int row = i+firstRowSegConstraint;
				d[row] = segConst.getMeanJointRate(); 
				if(wtedInversion) {
					if(segConst.isSlipRateConstraint())
						data_wt[row] = 1.0/(segConst.getStdDevOfMean()*(1-moRateReduction));
					else
						data_wt[row] = 1.0/segConst.getStdDevOfMean();
				}
				for(int col=0;col<numRuptures; col++) {
					if(rupSectionMatrix[sect1][col]==1 && rupSectionMatrix[sect2][col]==1) {
						if(segConst.isSlipRateConstraint()) {
							double aveSlip = (sectSlipInRup[sect1][col]+sectSlipInRup[sect2][col])/2.0;
							C[row][col] = aveSlip*(1-moRateReduction);
						}
						else {
							C[row][col]=1.0;
				// System.out.println("HERE: "+segConstraint_rupIndex[i]+"\t\t"+ segConstraint_rupRate[i] +"\t\t"+segConstraint_RupWt[i]);
						}
					}
				}
			}
		}
		
		// now fill in the MFD constraint if needed
		if(relativeMFD_constraintWt > 0.0) {
			for(int i=0; i < numMFD_constraints; i++) {
				int row = i+firstRowMFD_constraintData;
				double mag = mfdConstraint.getX(i);
				d[row] = mfdConstraint.getY(mag);
				if(wtedInversion && mfdSigma != null) {
					data_wt[row] = 1.0/mfdSigma.getY(mag);
				}
				for(int col=0; col<numRuptures; col++)
					if(mfdConstraint.getClosestXIndex(rupMeanMag[col]) == i)
						C[row][col]=1.0;
			}
		}
		
		// now fill in total rate constraint if needed
		if(relativeTotalRateConstraintWt>0) {
			int row = firstRowTotalRateConstraint;
			d[row] = totalRateConstraint;
			if(wtedInversion && totalRateSigma != 0)
				data_wt[row] = 1.0/totalRateSigma;
			for(int col=0; col<numRuptures; col++)
				C[row][col]=1.0;
		}
		
		// now fill in the smoothness constraints if requested
		if(relativeSmoothnessConstraintWt > 0.0 && smoothnessConstraintList != null) {
			int rowIncrement = 0;
			for(int[] sectArray : smoothnessConstraintList) {
				for(int s = 1; s<sectArray.length-1; s++) {
					int row = rowIncrement+firstRowSmoothnessConstraint;
					d[row] = 0.0; 
					for(int col=0; col<numRuptures; col++)
						C[row][col] = -rupSectionMatrix[s-1][col] + 2*rupSectionMatrix[s][col] - rupSectionMatrix[s+1][col];
					rowIncrement += 1;
				}
			}
		}

		
		
		// copy un-wted data to wted versions (wts added below)
		for(int row=0;row<totNumRows; row++){
			d_wted[row] = d[row];
			for(int col=0;col<numRuptures; col++)
				C_wted[row][col] = C[row][col];
		}
			

		// CORRECT IF MINIMUM RATE CONSTRAINT DESIRED
		if(minRupRate >0.0) {
			double[] Cmin = new double[totNumRows];  // the data vector
			// correct the data vector
			for(int row=0; row <totNumRows; row++) {
				for(int col=0; col < numRuptures; col++)
					Cmin[row]+=minRupRateArray[col]*C_wted[row][col];
				d_wted[row] -= Cmin[row];
			}
		}
		
		// APPLY WEIGHTS

		// segment slip rates first (no equation-set weight because others are relative)
		for(int row = firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row ++) {
			if(wtedInversion) 
				full_wt[row] = data_wt[row];
			else
				full_wt[row] = 1.0;
			d_wted[row] *= full_wt[row];
			for(int col=0; col<numRuptures; col++)
				C_wted[row][col] *= full_wt[row];
		}
		// sect event rate wts
		if(relativeSectRateWt > 0.0) {
			for(int i = 0; i < numSectRateConstraints; i ++) {
				int row = i+firstRowSectEventRateData;
				full_wt[row] = relativeSectRateWt;
				if(wtedInversion) full_wt[row] *= data_wt[row];
				d_wted[row] *= full_wt[row];
				for(int col=0; col<numRuptures; col++)
					C_wted[row][col] *= full_wt[row];
			}
		}
		// a-priori rate wts
		if(relative_aPrioriRupWt > 0.0) {
			for(int i=0; i < num_aPriori_constraints; i++) {
				int row = i+firstRowAprioriData;
				int col = aPriori_rupIndex[i];
				full_wt[row] = relative_aPrioriRupWt;
				if(wtedInversion) full_wt[row] *= data_wt[row];
				d_wted[row] *= full_wt[row];
				C_wted[row][col]*=full_wt[row];
			}
		}
		// seg const wts
		if(relative_segConstraintWt > 0.0) {
			for(int i=0; i < num_segConstraints; i++) {
				SegmentationConstraint segConst = segmentationConstraintList.get(i);
				int row = i+firstRowSegConstraint;
				full_wt[row] = relative_segConstraintWt;
				if(wtedInversion) full_wt[row] *= data_wt[row];
				d_wted[row] *= full_wt[row];
				for(int col=0;col<numRuptures; col++)
					if(rupSectionMatrix[segConst.getSect1_Index()][col]==1 && rupSectionMatrix[segConst.getSect2_Index()][col]==1)
						C_wted[row][col]*=full_wt[row];
			}
		}
		// MFD constraint wts
		if(relativeMFD_constraintWt > 0.0) {
			for(int i=0; i < numMFD_constraints; i++) {
				int row = i+firstRowMFD_constraintData;
				full_wt[row] = relativeMFD_constraintWt;
				if(wtedInversion) full_wt[row] *= data_wt[row];
				d_wted[row] *= full_wt[row];
				for(int rup=0; rup < numRuptures; rup++) 
					C_wted[row][rup] *= full_wt[row];
			}
		}
		
		// total rate constraint
		if(relativeTotalRateConstraintWt>0) {
			int row = firstRowTotalRateConstraint;
			full_wt[row] = relativeTotalRateConstraintWt;
			if(wtedInversion) full_wt[row] *= data_wt[row];
			d_wted[row] *= full_wt[row];
			for(int rup=0; rup<numRuptures; rup++)
				C_wted[row][rup] *= full_wt[row];

		}
		
		// smoothness constraint
		if(relativeSmoothnessConstraintWt>0) {
			for(int i=0; i<numSmoothnessConstraints; i++) {
				int row = i+firstRowSmoothnessConstraint;
				full_wt[row] = relativeSmoothnessConstraintWt;
//				if(wtedInversion) full_wt[row] *= data_wt[row]; // not data wt for smoothness constraint
				d_wted[row] *= full_wt[row];
				for(int rup=0; rup<numRuptures; rup++)
					C_wted[row][rup] *= full_wt[row];
			}
		}
	}
	
	
	/**
	 * This does the inversion using non-negative least squares
	 */
	public void doInversionNNLS() {
		
		modelRunInfoString = "\nInversion Type = Non-Negative Least Squares\n";
		
/*
		// manual check of matrices
			int nRow = C.length;
			int nCol = C[0].length;
			System.out.println("C = [");
			for(int i=0; i<nRow;i++) {
				for(int j=0;j<nCol;j++) 
					System.out.print(C[i][j]+"   ");
				System.out.print("\n");
			}
			System.out.println("];");
			System.out.println("d = [");
			for(int i=0; i<nRow;i++)
				System.out.println(d[i]);
			System.out.println("];");
*/

		// SOLVE THE INVERSE PROBLEM
		rupRateSolution = getNNLS_solution(C_wted, d_wted);

		// CORRECT FINAL RATES IF MINIMUM RATE CONSTRAINT APPLIED
		if(minRupRate >0.0)
			for(int rup=0; rup<numRuptures;rup++) rupRateSolution[rup] += minRupRateArray[rup];
		
		// set rates below the minimum to zero
		setRupRatesBelowMinToZero();

		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
		
		
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();
		
	}
	
	/**
	 * This allows one to set the solution by hand, such as the average of a bunch of SA solutions (minRupRate is ignored and potentially violated)
	 * @param rupRatesArray
	 */
	public void setSolution(double[] rupRatesArray, String info) {
		Preconditions.checkState(rupRatesArray.length == numRuptures, "rupRatesArray does not have the correct number of rupture (%s vs %s", rupRatesArray.length, numRuptures);
		
		if(info !=null)
			modelRunInfoString = "\n"+info+"\n";

		rupRateSolution = rupRatesArray;
		
		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();

	}
		
	
	/**
	 * This does the inversion using simulated annealing
	 */
	public void setMultipleSolutions(ArrayList<double[]> rupRatesArrayList, String info, String dirName) {
				
		// set these to null in case this method was already called
		rupRatesFromMultRunsArrayList = null;
		rupRatesFromMultRuns = null;
		mfdsFromMultRuns = null; 
		cumMfdsFromMultRuns = null; 
		finalSectSlipRateFromMultRuns = null;
		finalPaleoVisibleSectEventRateFromMultRuns = null;
		rateOfThroughGoingRupsAtSectBoudaryFromMultRuns = null;
		finalSectMeanSlipFromMultRuns = null;
		finalSectSlipCOV_FromMultRuns = null;

		if(info !=null)
			modelRunInfoString = "\n"+info+"\n";
		
		for(int solNum=0; solNum<rupRatesArrayList.size();solNum++) {
			
			rupRateSolution = rupRatesArrayList.get(solNum);
			Preconditions.checkState(rupRateSolution.length == numRuptures, "input rupture rates does not have the correct number of ruptures (%s vs %s", rupRateSolution.length, numRuptures);

			// compute predicted data
			d_pred = new double[totNumRows];  // predicted data vector
			for(int row=0;row<totNumRows; row++)
				for(int col=0; col <numRuptures; col++)
					d_pred[row] += rupRateSolution[col]*C[row][col];
					
			String solNumString = "\nFOR SOLUTION NUMBER "+solNum+":\n----------------------------\n";
			modelRunInfoString += solNumString;
			if(D) System.out.println(solNumString);
			
			// Compute stuff
			computeFinalStuff();
			computeSectMFDs();
			setMiscRunInfo();
									
			// create mult run data objects if currently null
			if(rupRatesFromMultRuns==null) {
				rupRatesFromMultRunsArrayList = rupRatesArrayList;
				rupRatesFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numRuptures, 1.0);
				mfdsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(magFreqDist.getMinX(), magFreqDist.size(), magFreqDist.getDelta()); 
				EvenlyDiscretizedFunc cumTemp = magFreqDist.getCumRateDistWithOffset();
				cumMfdsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(cumTemp.getMinX(), cumTemp.size(), cumTemp.getDelta()); 
				finalSectSlipRateFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				finalPaleoVisibleSectEventRateFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				rateOfThroughGoingRupsAtSectBoudaryFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(-0.5, rateOfThroughGoingRupsAtSectBoudary.length, 1.0);
				finalSectMeanSlipFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				finalSectSlipCOV_FromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
			}
			for(int i=0;i<numRuptures;i++)
				rupRatesFromMultRuns.set(i, rupRateSolution[i], 1.0);

			mfdsFromMultRuns.set(magFreqDist, 1.0);
			cumMfdsFromMultRuns.set(magFreqDist.getCumRateDistWithOffset(), 1.0);
			
			for(int s=0;s<numSections;s++) {
				finalSectSlipRateFromMultRuns.set(s, finalSectSlipRate[s], 1.0);
				finalPaleoVisibleSectEventRateFromMultRuns.set(s, finalPaleoVisibleSectEventRate[s], 1.0);
				finalSectMeanSlipFromMultRuns.set(s, finalSectMeanSlip[s], 1.0);
				finalSectSlipCOV_FromMultRuns.set(s, finalSectSlipCOV[s], 1.0);
			}
			
			for(int s=0;s<rateOfThroughGoingRupsAtSectBoudary.length;s++) {
				rateOfThroughGoingRupsAtSectBoudaryFromMultRuns.set(s, rateOfThroughGoingRupsAtSectBoudary[s], 1.0);
			}
		}
		
		modelRunInfoString += "\nFOR MEAN INVERSION:\n---------------------------------\n";
		EvenlyDiscretizedFunc meanRupRateFunc = rupRatesFromMultRuns.getMeanCurve();
		double[] rupRatesArray = new double[numRuptures];
		for(int i=0;i<numRuptures;i++)
			rupRatesArray[i] = meanRupRateFunc.getY(i);
		
		setSolution(rupRatesArray, null);

	}
	
	/**
	 * This does the inversion using simulated annealing
	 */
	public void doInversionSA(CompletionCriteria completionCriteria, double[] initialState, long randomSeed, CoolingScheduleType sa_coolingSchedule, 
			GenerationFunctionType perturbationFunc, IntegerPDF_FunctionSampler rupSampler) {
		
		this.sa_coolingSchedule = sa_coolingSchedule;
		
		
		modelRunInfoString = "\nInversion Type = Simulated Annearling with\n\n\tcompletionCriteria = "+completionCriteria.toString()+
				"\n\trandomSeed = "+randomSeed+
				"\n\tsa_CoolingSchedule = "+sa_coolingSchedule+
				"\n\tperturbationFunc = "+perturbationFunc.toString();
		if(initialState != null)
			modelRunInfoString += "\n\tinitialState is NOT null";
		else
			modelRunInfoString += "\n\tinitialState is null";
		if(rupSampler != null)
			modelRunInfoString += "\n\trupSampler is NOT null";
		else
			modelRunInfoString += "\n\trupSampler is null\n";

		// SOLVE THE INVERSE PROBLEM
//		rupRateSolution = getSimulatedAnnealingSolution(C_wted, d_wted, initialState, completionCriteria, randomSeed, perturbationFunc);

		// Not exactly sure how the following works, and doesn't appear to take a seed
		rupRateSolution = getSimulatedAnnealingThreadedSolution(C_wted, d_wted, initialState, completionCriteria, randomSeed, rupSampler);
		
		// CORRECT FINAL RATES IF MINIMUM RATE CONSTRAINT APPLIED
		if(minRupRate >0.0)
			for(int rup=0; rup<numRuptures;rup++) rupRateSolution[rup] += minRupRateArray[rup];
		
		// set rates below the minimum to zero
		setRupRatesBelowMinToZero();

		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();
		
	}

	
	private void setRupRatesBelowMinToZero() {
		for(int r=0;r<numRuptures;r++)
			if(rupRateSolution[r]<MIN_RUP_RATE)
				rupRateSolution[r]=0;
	}

	
	/**
	 * This does the inversion using simulated annealing
	 */
	public void doInversionSA_MultTimes(CompletionCriteria completionCriteria, double[] initialState, long randomSeed, int numInversions, String dirName,
			CoolingScheduleType sa_coolingSchedule, GenerationFunctionType perturbationFunc, IntegerPDF_FunctionSampler rupSampler) {
		
		this.sa_coolingSchedule = sa_coolingSchedule;
		
		// set these to null in case this method was already called
		rupRatesFromMultRunsArrayList=null;
		rupRatesFromMultRuns = null;
		mfdsFromMultRuns = null; 
		cumMfdsFromMultRuns = null; 
		finalSectSlipRateFromMultRuns = null;
		finalPaleoVisibleSectEventRateFromMultRuns = null;
		rateOfThroughGoingRupsAtSectBoudaryFromMultRuns = null;
		finalSectMeanSlipFromMultRuns = null;
		finalSectSlipCOV_FromMultRuns = null;

				
		modelRunInfoString = "\nInversion Type = Simulated Annearling with\n\n\tcompletionCriteria = "+
				completionCriteria.toString()+
				"\n\trandomSeed = "+randomSeed+
				"\n\tnumInversions = "+numInversions+
				"\n\tsa_CoolingSchedule = "+sa_coolingSchedule+
				"\n\tperturbationFunc = "+perturbationFunc.toString();
		if(initialState != null)
			modelRunInfoString += "\n\tinitialState is NOT null";
		else
			modelRunInfoString += "\n\tinitialState is null";
		if(rupSampler != null)
			modelRunInfoString += "\n\trupSampler is NOT null";
		else
			modelRunInfoString += "\n\trupSampler is null\n";

		
		for(int invNum=0; invNum<numInversions;invNum++) {
			
			randomSeed += invNum;
			
			// SOLVE THE INVERSE PROBLEM
//			rupRateSolution = getSimulatedAnnealingSolution(C_wted, d_wted, initialState, completionCriteria, randomSeed, perturbationFunc);
			
			// Not exactly sure how the following works, and doesn't appear to take a seed
			rupRateSolution = getSimulatedAnnealingThreadedSolution(C_wted, d_wted, initialState, completionCriteria, randomSeed, rupSampler);


			// CORRECT FINAL RATES IF MINIMUM RATE CONSTRAINT APPLIED
			if(minRupRate >0.0)
				for(int rup=0; rup<numRuptures;rup++) rupRateSolution[rup] += minRupRateArray[rup];
			
			// set rates below the minimum to zero
			setRupRatesBelowMinToZero();

			// compute predicted data
			d_pred = new double[totNumRows];  // predicted data vector
			for(int row=0;row<totNumRows; row++)
				for(int col=0; col <numRuptures; col++)
					d_pred[row] += rupRateSolution[col]*C[row][col];
					
			String invNumString = "\nFOR INVERSION NUMBER "+invNum+":\n----------------------------\n";
			modelRunInfoString += invNumString;
			if(D) System.out.println(invNumString);
			
			// Compute stuff
			computeFinalStuff();
			computeSectMFDs();
			setMiscRunInfo();
						
			// write these out now in case of crash
			String fileNamePrefix = null;
			if(dirName != null) {
				fileNamePrefix = dirName+"/ruptureRates_"+invNum;
				try{
					FileWriter fw = new FileWriter(fileNamePrefix+".txt");
					for(int i=0;i<numRuptures; i++) {				
						fw.write(i+"\t"+rupRateSolution[i]+"\n");
					}
					fw.close();
				}catch(Exception e) {
					e.printStackTrace();
				}		
			}
			
			// create mult run data objects if currently null
			if(rupRatesFromMultRunsArrayList==null) {
				rupRatesFromMultRunsArrayList = new ArrayList<double[]>();
				rupRatesFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numRuptures, 1.0);
				mfdsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(magFreqDist.getMinX(), magFreqDist.size(), magFreqDist.getDelta()); 
				EvenlyDiscretizedFunc cumTemp = magFreqDist.getCumRateDistWithOffset();
				cumMfdsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(cumTemp.getMinX(), cumTemp.size(), cumTemp.getDelta()); 
				finalSectSlipRateFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				finalPaleoVisibleSectEventRateFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				rateOfThroughGoingRupsAtSectBoudaryFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(-0.5, rateOfThroughGoingRupsAtSectBoudary.length, 1.0);
				finalSectMeanSlipFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
				finalSectSlipCOV_FromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(0, numSections, 1.0);
			}
			rupRatesFromMultRunsArrayList.add(rupRateSolution);
			for(int i=0;i<numRuptures;i++)
				rupRatesFromMultRuns.set(i, rupRateSolution[i], 1.0);

			mfdsFromMultRuns.set(magFreqDist, 1.0);
			cumMfdsFromMultRuns.set(magFreqDist.getCumRateDistWithOffset(), 1.0);
			
			for(int s=0;s<numSections;s++) {
				finalSectSlipRateFromMultRuns.set(s, finalSectSlipRate[s], 1.0);
				finalPaleoVisibleSectEventRateFromMultRuns.set(s, finalPaleoVisibleSectEventRate[s], 1.0);
				finalSectMeanSlipFromMultRuns.set(s, finalSectMeanSlip[s], 1.0);
				finalSectSlipCOV_FromMultRuns.set(s, finalSectSlipCOV[s], 1.0);
			}
			
			for(int s=0;s<rateOfThroughGoingRupsAtSectBoudary.length;s++) {
				rateOfThroughGoingRupsAtSectBoudaryFromMultRuns.set(s, rateOfThroughGoingRupsAtSectBoudary[s], 1.0);
			}

		}
		
		modelRunInfoString += "\nFOR MEAN INVERSION:\n---------------------------------\n";
		EvenlyDiscretizedFunc meanRupRateFunc = rupRatesFromMultRuns.getMeanCurve();
		double[] rupRatesArray = new double[numRuptures];
		for(int i=0;i<numRuptures;i++)
			rupRatesArray[i] = meanRupRateFunc.getY(i);
		
		setSolution(rupRatesArray, null);

	}

	
	
	/**
	 * This creates the segSlipInRup (Dsr) matrix based on the value of slipModelType.
	 * This slips are in meters.
	 *
	 */
	private void computeSectSlipInRupMatrix() {
		sectSlipInRup = new double[numSections][numRuptures];
		FaultSectionPrefData segData;
				
		// for case where ave slip computed from mag & area, and is same on all segments 
		if (slipModelType == SlipAlongRuptureModelEnum.UNIFORM) {
			for(int rup=0; rup<numRuptures; ++rup) {
				for(int seg=0; seg<numSections; seg++) {
					sectSlipInRup[seg][rup] = rupSectionMatrix[seg][rup]*rupAveSlip[rup];
				}
			}
		}
		else if (slipModelType == SlipAlongRuptureModelEnum.TAPERED) {
			// note that the ave slip is partitioned by area, not length; this is so the final model is moment balanced.
			mkTaperedSlipFuncs();
			for(int rup=0; rup<numRuptures; ++rup) {
				double normBegin=0, normEnd, scaleFactor;
				for(int seg=0; seg<numSections; seg++) {
					if(rupSectionMatrix[seg][rup]==1) {
						normEnd = normBegin + sectArea[seg]/rupArea[rup];
						if(normEnd > 1 && normEnd < 1.00001) normEnd = 1.0;
						scaleFactor = taperedSlipCDF.getInterpolatedY(normEnd)-taperedSlipCDF.getInterpolatedY(normBegin);
						scaleFactor /= (normEnd-normBegin);
						sectSlipInRup[seg][rup] = rupAveSlip[rup]*scaleFactor;
						normBegin = normEnd;
					}
				}
				/*
				if(rup == num_rup-1) { // check results
					double d_aveTest=0;
					for(int seg=0; seg<num_seg; seg++)
						d_aveTest += segSlipInRup[seg][rup]*segArea[seg]/rupArea[rup];
					System.out.println("AveSlipCheck: " + (float) (d_aveTest/aveSlip));
				}
				*/
			}
		}
		else throw new RuntimeException("slip model not supported");
		
//		for(int seg=0; seg<num_seg; seg++) System.out.println(seg+"\t"+segSlipInRup[seg][num_rup-1]);
		
	}
	
	
	/**
	 * This makes a tapered slip function based on the [Sin(x)]^0.5 fit of 
	 * Biasi & Weldon (in prep; pesonal communication), which is based on  
	 * the data comilation of Biasi & Weldon (2006, "Estimating Surface  
	 * Rupture Length and Magnitude of Paleoearthquakes from Point 
	 * Measurements of Rupture Displacement", Bull. Seism. Soc. Am. 96, 
	 * 1612-1623, doi: 10.1785/0120040172 E)
	 *
	 */
	private static void mkTaperedSlipFuncs() {
		
		// only do if another instance has not already done this
		if(taperedSlipCDF != null) return;
		
		taperedSlipPDF = getTaperedSlipFunc();
		taperedSlipCDF = taperedSlipPDF.getCumulativeDistFunction();	

	}
	
	
	/**
	 * This makes a tapered slip function based on the [Sin(x)]^0.5 fit of 
	 * Biasi & Weldon (in prep; pesonal communication), which is based on  
	 * the data comilation of Biasi & Weldon (2006, "Estimating Surface  
	 * Rupture Length and Magnitude of Paleoearthquakes from Point 
	 * Measurements of Rupture Displacement", Bull. Seism. Soc. Am. 96, 
	 * 1612-1623, doi: 10.1785/0120040172 E)
	 *
	 */
	public static HistogramFunction getTaperedSlipFunc() {
		HistogramFunction taperedSlipPDF = new HistogramFunction(0, 5001, 0.0002);
		double x,y;
		for(int i=0; i<taperedSlipPDF.size();i++) {
			x = taperedSlipPDF.getX(i);
			y = Math.pow(Math.sin(x*Math.PI), 0.5);
			taperedSlipPDF.set(i,y);
		}
		taperedSlipPDF.normalizeBySumOfY_Vals();
		return taperedSlipPDF;
	}
	
	
	/**
	 * This computes the total event rate prediction error (ignoring wt given to this equation set),
	 * meaning it will give a result even if the equation set wt was zero.
	 * @return
	 */
	private double getTotEventRatePredErr() {
		SectionRateConstraint constraint;
		double totPredErr=0;
		for(int i = 0; i < sectionRateConstraints.size(); i ++) {
			constraint = sectionRateConstraints.get(i);
			int sect = constraint.getSectIndex();
			double normResid = (finalPaleoVisibleSectEventRate[sect]-constraint.getMeanRate())/constraint.getStdDevOfMean();
			totPredErr += normResid*normResid;
		}
		return totPredErr;
	}
	
		
	private double[] getSimulatedAnnealingSolution(double[][] C, double[] d, double[] initialState, CompletionCriteria completionCriteria,  long randomSeed,
			GenerationFunctionType perturbationFunc, IntegerPDF_FunctionSampler rupSampler) {
		SparseDoubleMatrix2D matrixC = new SparseDoubleMatrix2D(C); //
		SerialSimulatedAnnealing simulatedAnnealing =new SerialSimulatedAnnealing(matrixC, d, initialState);
		simulatedAnnealing.setCoolingFunc(sa_coolingSchedule);
		simulatedAnnealing.setRandom(new Random(randomSeed));
		simulatedAnnealing.setPerturbationFunc(perturbationFunc);
		simulatedAnnealing.setRuptureSampler(rupSampler);
		simulatedAnnealing.iterate(completionCriteria);
		return simulatedAnnealing.getBestSolution();
	}



	private double[] getSimulatedAnnealingThreadedSolution(double[][] C, double[] d, double[] initialState, CompletionCriteria completionCriteria,  
			long randomSeed, IntegerPDF_FunctionSampler rupSampler) {
		SparseDoubleMatrix2D matrixC = new SparseDoubleMatrix2D(C); //
		//this is the "sub completion criteria" - the amount of time (or iterations) between synchronization
		CompletionCriteria subCompetionCriteria = TimeCompletionCriteria.getInSeconds(1); // 1 second;
		// this will use all available processors
		int numThreads = Runtime.getRuntime().availableProcessors();

		ThreadedSimulatedAnnealing simulatedAnnealing = new ThreadedSimulatedAnnealing(
				matrixC, d, initialState, numThreads, subCompetionCriteria);
		simulatedAnnealing.setRuptureSampler(rupSampler);
//		simulatedAnnealing.setRandom(new Random(randomSeed));
		simulatedAnnealing.iterate(completionCriteria);
		return simulatedAnnealing.getBestSolution();
	}


	
	
	/**
	 * This gets the non-negative least squares solution for the matrix C
	 * and data vector d.
	 * @param C
	 * @param d
	 * @return
	 */
	public static double[] getNNLS_solution(double[][] C, double[] d) {

		int nRow = C.length;
		int nCol = C[0].length;
		
		double[] A = new double[nRow*nCol];
		double[] x = new double[nCol];
		
		int i,j,k=0;
	
		if(MATLAB_TEST) {
			System.out.println("display "+"SSAF Inversion test");
			System.out.println("C = [");
			for(i=0; i<nRow;i++) {
				for(j=0;j<nCol;j++) 
					System.out.print(C[i][j]+"   ");
				System.out.print("\n");
			}
			System.out.println("];");
			System.out.println("d = [");
			for(i=0; i<nRow;i++)
				System.out.println(d[i]);
			System.out.println("];");
		}
/////////////////////////////////////
		
		for(j=0;j<nCol;j++) 
			for(i=0; i<nRow;i++)	{
				A[k]=C[i][j];
				k+=1;
			}
		nnls.update(A,nRow,nCol);
		
		boolean converged = nnls.solve(d,x);
		if(!converged)
			throw new RuntimeException("ERROR:  NNLS Inversion Failed");
		
		if(MATLAB_TEST) {
			System.out.println("x = [");
			for(i=0; i<x.length;i++)
				System.out.println(x[i]);
			System.out.println("];");
			System.out.println("max(abs(x-lsqnonneg(C,d)))");
		}
		
		return x;
	}
	
	/**
	 * Computer Final Slip Rate for each segment (& aPrioriSegSlipRate)
	 *
	 */
	private void computeFinalStuff() {
		
		// compute segment slip and event rates
		finalSectSlipRate = new double[numSections];
		finalSectEventRate = new double[numSections];
		finalPaleoVisibleSectEventRate = new double[numSections];
		for(int s=0; s < numSections; s++) {
			finalSectSlipRate[s] = 0;
			finalSectEventRate[s] = 0;
			for(int rup=0; rup < numRuptures; rup++) 
				if(rupSectionMatrix[s][rup]==1) {
					finalSectSlipRate[s] += rupRateSolution[rup]*sectSlipInRup[s][rup];
					finalSectEventRate[s]+=rupRateSolution[rup];
					if(applyProbVisible)
						finalPaleoVisibleSectEventRate[s]+=rupRateSolution[rup]*getProbVisible(rupMeanMag[rup]);
					else
						finalPaleoVisibleSectEventRate[s]+=rupRateSolution[rup];
				}
		}
		
		
		// compute section mean slip and cov
		finalSectMeanSlip = new double[numSections];
		finalSectSlipCOV = new double[numSections];
		for(int s=0; s < numSections; s++) {
			ArbDiscrEmpiricalDistFunc dist = new ArbDiscrEmpiricalDistFunc();
			for(int r=0; r < numRuptures; r++) {
				if(rupRateSolution[r] > 0 && rupSectionMatrix[s][r] == 1) {
					dist.set(sectSlipInRup[s][r]*gaussMFD_slipCorr, rupRateSolution[r]);
				}
			}
			finalSectMeanSlip[s] = dist.getMean();
			finalSectSlipCOV[s] = dist.getCOV();
		}

		
		// Compute the total Mag Freq Dist
		int num = (int)Math.round((maxMagMFD_WithAleatory-minMagMFD_WithAleatory)/MAG_DELTA + 1);
		if(num==1) num=2;
		magFreqDist = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
		double totRupMoRate = 0;
		for(int rup=0; rup<numRuptures;rup++) {
			if(GAUSS_MFD_SIGMA==0d) {
				int index = magFreqDist.getClosestXIndex(rupMeanMag[rup]);
				magFreqDist.add(index, rupRateSolution[rup]);
				// the following produces different results
//				magFreqDist.addResampledMagRate(rupMeanMag[rup], rupRateSolution[rup], true);
				totRupMoRate += rupRateSolution[rup]*rupMeanMo[rup];
			}
			else {
				GaussianMagFreqDist gDist = new GaussianMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA,rupMeanMag[rup],GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
				gDist.scaleToCumRate(0, rupRateSolution[rup]);
				magFreqDist.addIncrementalMagFreqDist(gDist);
//				totRupMoRate += magFreqDist.getTotalMomentRate();
				totRupMoRate += rupRateSolution[rup]*rupMeanMo[rup]*gaussMFD_slipCorr;
			}
		}
		magFreqDist.setInfo("Incremental Mag Freq Dist");
		
		
		// check moment rate
		double origMoRate = totMoRate*(1-moRateReduction);
		double mfdRatio = origMoRate/magFreqDist.getTotalMomentRate();
		String tempString = "Moment Rates: \n\n\tFrom Fault Sections (possibly reduced) = "+ (float)origMoRate+
				"\n\tSum From Rups = "+(float)totRupMoRate+",  ratio = "+(float)(totRupMoRate/origMoRate)+
				"\n\tFrom total MFD = "+(float)magFreqDist.getTotalMomentRate()+
				",  ratio = "+(float)mfdRatio;
		if(D)
			System.out.println(tempString+"\n");
		modelRunInfoString += "\n"+ tempString+"\n";
		
		
		// COMPUTE RATE AT WHICH SECTION BOUNDARIES CONSTITUTE RUPTURE ENDPOINTS
		rateOfThroughGoingRupsAtSectBoudary = new double[numSections+1];  // there is one more boundary than sections
		for(int rup=0; rup<numRuptures;rup++) {
			int beginBoundary = firstSectOfRup[rup];
			int endBoundary = firstSectOfRup[rup]+numSectInRup[rup];
			for(int s=beginBoundary+1;s<endBoundary;s++)
				rateOfThroughGoingRupsAtSectBoudary[s] += rupRateSolution[rup];
		}
	}
	
	
	
	/**
	 * This computes both participation and nucleation MFDs for each sub-section
	 */
	private void computeSectMFDs() {
		int num = (int)Math.round((maxMagMFD_WithAleatory-minMagMFD_WithAleatory)/MAG_DELTA + 1);
		if(num==1) num = 2;
		sectNucleationMFDs = new ArrayList<SummedMagFreqDist>();
		sectParticipationMFDs = new ArrayList<SummedMagFreqDist>();
		SummedMagFreqDist sumOfSegPartMFDs = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
		aveOfSectPartMFDs = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
		
		SummedMagFreqDist segPartMFD, segNuclMFD;
//		double mag, rate;
		for(int seg=0; seg < numSections; seg++) {
			segPartMFD = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
			segNuclMFD = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
			for(int rup=0; rup < numRuptures; rup++) {
				if(this.rupSectionMatrix[seg][rup] == 1) {
					if(rupRateSolution[rup] > 0) {
						double mag;
						if(GAUSS_MFD_SIGMA==0d) {
							int index = segPartMFD.getClosestXIndex(rupMeanMag[rup]);
							mag = segPartMFD.getX(index);
							// not sure the following will always produce the same result as above
//							mag = this.roundMagTo10thUnit(rupMeanMag[rup]);
						}
						else
							mag = rupMeanMag[rup];
						GaussianMagFreqDist mfd = new GaussianMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
						mfd.scaleToCumRate(0, rupRateSolution[rup]);
						segPartMFD.addIncrementalMagFreqDist(mfd);
						mfd.scaleToCumRate(0, rupRateSolution[rup]/numSectInRup[rup]); // assume uniform distribution of nucleations
						segNuclMFD.addIncrementalMagFreqDist(mfd);						
					}
				}
			}
			sectNucleationMFDs.add(segNuclMFD);
			sectParticipationMFDs.add(segPartMFD);
			sumOfSegPartMFDs.addIncrementalMagFreqDist(segPartMFD);
		}
		// compute aveOfSegPartMFDs from sumOfSegPartMFDs
		for(int m=0; m<sumOfSegPartMFDs.size();m++) aveOfSectPartMFDs.add(m, sumOfSegPartMFDs.getY(m)/numSections);
		aveOfSectPartMFDs.setInfo("Average Seg Participation MFD");
		
	}
	
	
	
	
	private void setMiscRunInfo() {

		// write out rupture rates and mags
//		System.out.println("Final Rupture Rates & Mags:");
//		for(int rup=0; rup < num_rup; rup++)
//		System.out.println(rup+"\t"+(float)rupRateSolution[rup]+"\t"+(float)rupMeanMag[rup]);
		
		modelRunInfoString += "\nMisc. Run Info:\n\n";
		
		// WRITE OUT MAXIMUM EVENT RATE & INDEX
		double maxRate=0;
		int index=-1;
		for(int seg = 0; seg < numSections; seg ++)
			if(finalSectEventRate[seg] > maxRate) {
				maxRate = finalSectEventRate[seg];
				index = seg;
			}
		int mri = (int)Math.round(1.0/maxRate);
		if(D) System.out.println("\nMax sect rate (MRI): "+(float)maxRate+" ("+mri+" yrs) at index "+index);
		modelRunInfoString += "\tMax sect rate (MRI): "+(float)maxRate+" ("+mri+" yrs) at index "+index+"\n";

		
		//write out number of ruptures that have rates above minRupRate
		int numAbove = 0;
		for(int rup=0; rup<this.rupRateSolution.length; rup++)
			if(rupRateSolution[rup] > minRupRate) numAbove += 1;
		if(D) System.out.println("\nNum Ruptures above minRupRate = "+numAbove+"\t(out of "+rupRateSolution.length+")\n");
		modelRunInfoString += "\tNum Ruptures above minRupRate = "+numAbove+"\t(out of "+rupRateSolution.length+")\n";

		// write out final segment slip rates
//		System.out.println("\nSegment Slip Rates: index, final, orig, and final/orig (orig is corrected for moRateReduction)");
		double aveRatio=0;
		for(int seg = 0; seg < numSections; seg ++) {
			double slipRate = sectSlipRate[seg]*(1-this.moRateReduction);
			aveRatio += finalSectSlipRate[seg]/slipRate;
//			System.out.println(seg+"\t"+(float)finalSegSlipRate[seg]+"\t"+(float)slipRate+"\t"+(float)(finalSegSlipRate[seg]/slipRate));
		}
		aveRatio /= numSections;
		if(D)System.out.println("Ave final/orig slip rate = "+(float)aveRatio);
		modelRunInfoString += "\tAve final/orig slip rate = "+(float)aveRatio+"\n";


		// write out final segment event rates
		if(relativeSectRateWt > 0.0) {
//			System.out.println("\nSegment Rates: index, final, orig, and final/orig");
			aveRatio = 0;
			SectionRateConstraint constraint;
			for(int i = 0; i < sectionRateConstraints.size(); i ++) {
				int row = firstRowSectEventRateData+i;
				constraint = sectionRateConstraints.get(i);
				int seg = constraint.getSectIndex();
//				this checks that finalSegEventRate[seg] and d_pred[row} are the same (RECONSIDER PROB VISIBLE IF I USE THIS AGAIN)
//				System.out.println(seg+"\t"+(float)finalSegEventRate[seg]+"\t"+(float)d_pred[row]+"\t"+(float)(finalSegEventRate[seg]/d_pred[row]));
				if(applyProbVisible) {
					aveRatio += finalPaleoVisibleSectEventRate[seg]/constraint.getMeanRate();
//					System.out.println(seg+"\t"+(float)finalPaleoReducedSegEventRate[seg]+"\t"+(float)constraint.getMean()+"\t"+(float)aveRatio);				
				}
				else {
					aveRatio += finalSectEventRate[seg]/constraint.getMeanRate();
//					System.out.println(seg+"\t"+(float)finalSegEventRate[seg]+"\t"+(float)constraint.getMean()+"\t"+(float)aveRatio);					
				}
			}
			aveRatio /= sectionRateConstraints.size();
			if(D) System.out.println("Ave final/orig segment rate = "+(float)aveRatio);
			modelRunInfoString += "\tAve final/orig segment rate = "+(float)aveRatio+"\n";
		}

		// write out final rates for ruptures with an a-priori constraint
		if(this.relative_aPrioriRupWt >0.0) {
			aveRatio = 0;
			//		System.out.println("\nA Priori Rates: index, final, orig, and final/orig");
			for(int i=0; i<num_aPriori_constraints;i++) {
				double ratio;
				if(rupRateSolution[aPriori_rupIndex[i]] > 1e-14 && aPriori_rate[i] > 1e-14)  // if both are not essentially zero
					ratio = (rupRateSolution[aPriori_rupIndex[i]]/aPriori_rate[i]);
				else
					ratio = 1;
				aveRatio += ratio;
				//			System.out.println(aPriori_rupIndex[i]+"\t"+(float)rupRateSolution[aPriori_rupIndex[i]]+"\t"+aPriori_rate[i]+"\t"+(float)ratio);				
			}
			aveRatio /= num_aPriori_constraints;
			if(D) System.out.println("Ave final/orig a-priori rate = "+(float)aveRatio);
			modelRunInfoString += "\tAve final/orig a-priori rate = "+(float)aveRatio+"\n";
		}

		// Now do the Model error metrics
		
		// First without equation set weights
		String resultsString = getErrorTableString(false);
		if(D) System.out.println(resultsString);
		modelRunInfoString += "\n"+resultsString+"\n";
				
		// Now with equation weights
		resultsString = getErrorTableString(true);
		if(D) System.out.println(resultsString);
		modelRunInfoString += "\n"+resultsString+"\n";
			
	}
	
	
	private String getErrorTableString(boolean includeEquationSetWts) {
		// First without equation set weights
		double totPredErr=0, slipRateErr=0, eventRateErr=0, aPrioriErr=0, segConstrErr=0, mfdErr=0, totRateErr=0, totSmoothnessErr=0;
		double totRMS=0, slipRateRMS=0, eventRateRMS=0, aPrioriRMS=0, segConstrRMS=0, mfdRMS=0, totRateRMS=0, totSmoothnessRMS=0;
		double totAveAbsDiff=0, slipRateDiff=0, eventRateDiff=0, aPrioriDiff=0, segConstrDiff=0, mfdDiff=0, totRateDiff=0, totSmoothnessDiff=0; 	// these are for ave absolute value of diffs
		double testNumRows=0, slipRateNumRows=0, eventRateNumRows=0, aPrioriNumRows=0, segConstrNumRows=0, mfdNumRows=0, totRateNumRows=0, totSmoothnessNumRows=0;

		double[] wt;
		if(includeEquationSetWts)
			wt = full_wt;
		else
			wt = data_wt;

		
		for(int row=firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row++) {
			slipRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
			slipRateDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
			slipRateNumRows += 1;
		}
		if(relativeSectRateWt >0)
			for(int row=firstRowSectEventRateData; row <= lastRowSectEventRateData; row++) {
				eventRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
				eventRateDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
				eventRateNumRows += 1;
			}
		if(relative_aPrioriRupWt > 0)
			for(int row=firstRowAprioriData; row <= lastRowAprioriData; row++) {
				aPrioriErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
				aPrioriDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
				aPrioriNumRows += 1;
			}
		if(relative_segConstraintWt > 0)
			for(int row=firstRowSegConstraint; row <= lastRowSegConstraint; row++) {
//System.out.println("segConstraint HERE:"+d[row]+"\t"+d_pred[row]+"\t"+wt[row]);
				segConstrErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
				segConstrDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
				segConstrNumRows += 1;
			}
		if(relativeMFD_constraintWt>0)
			for(int row=firstRowMFD_constraintData; row <= lastRowMFD_constraintData; row++){
				mfdErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
				mfdDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
				mfdNumRows += 1;
			}
		if(relativeTotalRateConstraintWt>0) {
			int row=firstRowTotalRateConstraint;
			totRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
			totRateDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
			totRateNumRows += 1;
		}
		if(relativeSmoothnessConstraintWt>0) {
			for(int row=firstRowSmoothnessConstraint; row <= lastRowSmoothnessConstraint; row++) {
				totSmoothnessErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*wt[row]*wt[row];
				totSmoothnessDiff += Math.abs((d[row]-d_pred[row])*wt[row]);
				totSmoothnessNumRows += 1;		
			}
		}

		testNumRows = slipRateNumRows+eventRateNumRows+aPrioriNumRows+segConstrNumRows+mfdNumRows+totRateNumRows+totSmoothnessNumRows;
		if(testNumRows != totNumRows) {
			throw new RuntimeException("Problem with number of rows");
		}

		totPredErr = slipRateErr + eventRateErr + aPrioriErr + segConstrErr + mfdErr + totRateErr + totSmoothnessErr;
		totRMS = Math.sqrt(totPredErr/totNumRows);
		totAveAbsDiff = (slipRateDiff+eventRateDiff+aPrioriDiff+segConstrDiff+mfdDiff+totRateDiff+totSmoothnessDiff)/totNumRows;
		
		slipRateRMS = Math.sqrt(slipRateErr/slipRateNumRows);
		eventRateRMS = Math.sqrt(eventRateErr/eventRateNumRows);
		aPrioriRMS = Math.sqrt(aPrioriErr/aPrioriNumRows);
		segConstrRMS = Math.sqrt(segConstrErr/segConstrNumRows);
		mfdRMS = Math.sqrt(mfdErr/mfdNumRows); 
		totRateRMS = Math.sqrt(totRateErr/totRateNumRows);
		totSmoothnessRMS = Math.sqrt(totSmoothnessErr/totSmoothnessNumRows);
		
		slipRateDiff /= slipRateNumRows;
		eventRateDiff /= eventRateNumRows;
		aPrioriDiff /= aPrioriNumRows;
		segConstrDiff /= segConstrNumRows;
		mfdDiff /= mfdNumRows;
		totRateDiff /= totRateNumRows;
		totSmoothnessDiff /= totSmoothnessNumRows;
		
		// get alt versions in case equation wts zero
		double eventRateErrAlt = getTotEventRatePredErr();
		
		String resultsString = "";
		if(includeEquationSetWts)
			resultsString = 	"\nFinal Energy/Error Including Equation Set Wts:"+
								"\n---------------------------------------------";
		else
			resultsString = 	"\nFinal Energy/Error Ignoring Equation Set Wts:"+
								"\n--------------------------------------------";

		resultsString += "\n         \tError    \tRMS      \tAveAbsDiff\tRelWt    \n"+
								"Total    \t"+(float)totPredErr+"\t"+(float)totRMS+"\t"+(float)totAveAbsDiff+"\tNA\n"+
								"SlipRate \t"+(float)slipRateErr+"\t"+(float)slipRateRMS+"\t"+(float)slipRateDiff+"\t1.0\n";
		if(eventRateNumRows>0)
			resultsString += "EventRate \t"+(float)eventRateErr+"\t"+(float)eventRateRMS+"\t"+(float)eventRateDiff+"\t"+relativeSectRateWt+"\n";
		if(aPrioriNumRows>0)
			resultsString += "A Priori  \t"+(float)aPrioriErr+"\t"+(float)aPrioriRMS+"\t"+(float)aPrioriDiff+"\t"+relative_aPrioriRupWt+"\n";
		if(segConstrNumRows>0)
			resultsString += "Seg Constr  \t"+(float)segConstrErr+"\t"+(float)segConstrRMS+"\t"+(float)segConstrDiff+"\t"+relative_segConstraintWt+"\n";
		if(mfdNumRows>0) 
			resultsString += "MFD Rates  \t"+(float)mfdErr+"\t"+(float)mfdRMS+"\t"+(float)mfdDiff+"\t"+relativeMFD_constraintWt+"\n";
		if(totRateNumRows>0)
			resultsString += "Total Rate \t"+(float)totRateErr+"\t"+(float)totRateRMS+"\t"+(float)totRateDiff+"\t"+relativeTotalRateConstraintWt+"\n";
		if(totSmoothnessNumRows>0)
			resultsString += "Total Smoothness \t"+(float)totSmoothnessErr+"\t"+(float)totSmoothnessRMS+"\t"+(float)totSmoothnessDiff+"\t"+relativeSmoothnessConstraintWt+"\n";
		
		return resultsString;
	}
	
	
	/**
	 * This plots/writes out the participation MFD for the given section
	 * @param dirName
	 * @param popupWindow
	 * @param targetSect
	 */
	public void writeAndOrPlotPartMFD_ForSection(String dirName, boolean popupWindow, int targetSect) {
		
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		SummedMagFreqDist avePartMFD = sectParticipationMFDs.get(targetSect);


		if(rupRatesFromMultRuns != null) {
			
			ArbDiscrEmpiricalDistFunc_3D partMFDsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(avePartMFD.getMinX(), avePartMFD.size(), avePartMFD.getDelta());
			
			for(double[] rupRatesArray : rupRatesFromMultRunsArrayList) {
				SummedMagFreqDist sectPartMFD = new SummedMagFreqDist(avePartMFD.getMinX(),avePartMFD.size(),MAG_DELTA);
				
				for(int rup=0; rup < rupRatesArray.length; rup++) {
					if(rupSectionMatrix[targetSect][rup] == 1) {
						double rupRate = rupRatesArray[rup];
						if(rupRate > 0) {
							double mag;
							if(GAUSS_MFD_SIGMA==0d) {
								int index = sectPartMFD.getClosestXIndex(rupMeanMag[rup]);
								mag = sectPartMFD.getX(index);
							}
							else
								mag = rupMeanMag[rup];
							GaussianMagFreqDist mfd = new GaussianMagFreqDist(avePartMFD.getMinX(),avePartMFD.size(),MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
							mfd.scaleToCumRate(0, rupRate);
							sectPartMFD.addIncrementalMagFreqDist(mfd);
						}
					}
				}
				partMFDsFromMultRuns.set(sectPartMFD, 1.0);
			}
			
			UncertainArbDiscFunc mfdMeanMinMaxRange = new UncertainArbDiscFunc(partMFDsFromMultRuns.getMeanCurve(), 
					partMFDsFromMultRuns.getMinCurve(), partMFDsFromMultRuns.getMaxCurve());
			mfdMeanMinMaxRange.setName("partMFD_MeanMinMaxRange");
			mfdList.add(mfdMeanMinMaxRange);
			
			UncertainArbDiscFunc mfdMean95conf = get95perConfForMultRuns(partMFDsFromMultRuns);
			mfdMean95conf.setName("partMFD_Mean95conf");
			mfdList.add(mfdMean95conf);
			
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,200)));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,120)));

		}
		
		mfdList.add(avePartMFD);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		//
		SummedMagFreqDist sectMagProbMFD = new SummedMagFreqDist(avePartMFD.getMinX(),avePartMFD.size(),MAG_DELTA);
		for(int rup=0; rup < numRuptures; rup++) {
			if(rupSectionMatrix[targetSect][rup] == 1) {
				double mag;
				if(GAUSS_MFD_SIGMA==0d) {
					int index = sectMagProbMFD.getClosestXIndex(rupMeanMag[rup]);
					mag = sectMagProbMFD.getX(index);
				}
				else
					mag = rupMeanMag[rup];
				GaussianMagFreqDist mfd = new GaussianMagFreqDist(avePartMFD.getMinX(),avePartMFD.size(),MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
				mfd.scaleToCumRate(0, 1.0);
				sectMagProbMFD.addIncrementalMagFreqDist(mfd);
			}
		}
		sectMagProbMFD.scaleToCumRate(0, avePartMFD.calcSumOfY_Vals());
		sectMagProbMFD.setName("sectMagProbMFD");
		sectMagProbMFD.setInfo("Relative probability of sampling magnitude");
		mfdList.add(sectMagProbMFD);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sect"+targetSect+"partMFD";
		String plotName = "Section "+targetSect+" Part. MFD";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Rate (per yr)";
		Range xAxisRange = new Range(6.0, 9.0);
		Range yAxisRange = new Range(1e-8, 1e-1);
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		// write alternative rupture rates file to get rid of big header
		try{
			FileWriter fw = new FileWriter(fileNamePrefix+"Alt.txt");
			for(int i=0;i<numRuptures; i++) {				
				fw.write(i+"\t"+rupRateSolution[i]+"\n");
			}
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}		
	}
	
	/**
	 * This plots/writes out the participation MFD for the given section
	 * @param dirName
	 * @param popupWindow
	 * @param targetSect
	 */
	public void writeAndOrPlotJointPartMFD_ForSections(String dirName, boolean popupWindow, int sect1, int sect2) {
		
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		SummedMagFreqDist tempMFD = sectParticipationMFDs.get(sect1);
		double tempSum=0;

		if(rupRatesFromMultRuns != null) {
			
			ArbDiscrEmpiricalDistFunc_3D partMFDsFromMultRuns = new ArbDiscrEmpiricalDistFunc_3D(tempMFD.getMinX(), tempMFD.size(), tempMFD.getDelta());
			
			for(double[] rupRatesArray : rupRatesFromMultRunsArrayList) {
				SummedMagFreqDist sectPartMFD = new SummedMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA);
				
				for(int rup=0; rup < rupRatesArray.length; rup++) {
					if(rupSectionMatrix[sect1][rup] == 1 && rupSectionMatrix[sect2][rup] == 1) {
						double rupRate = rupRatesArray[rup];
						if(rupRate > 0) {
							double mag;
							if(GAUSS_MFD_SIGMA==0d) {
								int index = sectPartMFD.getClosestXIndex(rupMeanMag[rup]);
								mag = sectPartMFD.getX(index);
							}
							else
								mag = rupMeanMag[rup];
							GaussianMagFreqDist mfd = new GaussianMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
							mfd.scaleToCumRate(0, rupRate);
							sectPartMFD.addIncrementalMagFreqDist(mfd);
						}
					}
				}
				partMFDsFromMultRuns.set(sectPartMFD, 1.0);
			}
			tempSum = partMFDsFromMultRuns.getMeanCurve().calcSumOfY_Vals();
			UncertainArbDiscFunc mfdMeanMinMaxRange = new UncertainArbDiscFunc(partMFDsFromMultRuns.getMeanCurve(), 
					partMFDsFromMultRuns.getMinCurve(), partMFDsFromMultRuns.getMaxCurve());
			mfdMeanMinMaxRange.setName("jointPartMFD_MeanMinMaxRange");
			mfdList.add(mfdMeanMinMaxRange);
			
			UncertainArbDiscFunc mfdMean95conf = get95perConfForMultRuns(partMFDsFromMultRuns);
			mfdMean95conf.setName("jointPartMFD_Mean95conf");
			mfdList.add(mfdMean95conf);
			
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,200)));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,120)));
			
			EvenlyDiscretizedFunc meanMFD = partMFDsFromMultRuns.getMeanCurve();
			meanMFD.setName("Joint Part. MFD");
			mfdList.add(meanMFD);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		}
		else {
			SummedMagFreqDist sectPartMFD = new SummedMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA);
			
			for(int rup=0; rup < numRuptures; rup++) {
				if(rupSectionMatrix[sect1][rup] == 1 && rupSectionMatrix[sect2][rup] == 1) {
					double rupRate = this.rupRateSolution[rup];
					if(rupRate > 0) {
						double mag;
						if(GAUSS_MFD_SIGMA==0d) {
							int index = sectPartMFD.getClosestXIndex(rupMeanMag[rup]);
							mag = sectPartMFD.getX(index);
						}
						else
							mag = rupMeanMag[rup];
						GaussianMagFreqDist mfd = new GaussianMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
						mfd.scaleToCumRate(0, rupRate);
						sectPartMFD.addIncrementalMagFreqDist(mfd);
					}
				}
			}
			sectPartMFD.setName("Joint Part. MFD");
			mfdList.add(sectPartMFD);
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		}
		

		
		//
		SummedMagFreqDist sectMagProbMFD = new SummedMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA);
		for(int rup=0; rup < numRuptures; rup++) {
			if(rupSectionMatrix[sect1][rup] == 1 && rupSectionMatrix[sect2][rup] == 1) {
				double mag;
				if(GAUSS_MFD_SIGMA==0d) {
					int index = sectMagProbMFD.getClosestXIndex(rupMeanMag[rup]);
					mag = sectMagProbMFD.getX(index);
				}
				else
					mag = rupMeanMag[rup];
				GaussianMagFreqDist mfd = new GaussianMagFreqDist(tempMFD.getMinX(),tempMFD.size(),MAG_DELTA,mag,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
				mfd.scaleToCumRate(0, 1.0);
				sectMagProbMFD.addIncrementalMagFreqDist(mfd);
			}
		}
		
		sectMagProbMFD.scaleToCumRate(0, tempSum);
		sectMagProbMFD.setName("sectMagProbMFD");
		sectMagProbMFD.setInfo("Relative probability of sampling magnitude");
		mfdList.add(sectMagProbMFD);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sect"+sect1+"and"+sect2+"_JointPartMFD";
		String plotName = "Subsection "+sect1+" & "+sect2+" Joint Part. MFD";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Rate (per yr)";
		Range xAxisRange = new Range(6.0, 9.0);
		Range yAxisRange = new Range(1e-8, 1e-1);
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		// write alternative rupture rates file to get rid of big header
		try{
			FileWriter fw = new FileWriter(fileNamePrefix+"Alt.txt");
			for(int i=0;i<numRuptures; i++) {				
				fw.write(i+"\t"+rupRateSolution[i]+"\n");
			}
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}		
	}
	
	
	public void computeRupSamplerMFD(IntegerPDF_FunctionSampler rupSampler) {
		
		rupSamplerMFD = new SummedMagFreqDist(meanMagHistorgram.getMinX(),meanMagHistorgram.size(),meanMagHistorgram.getDelta());
		for(int rup=0; rup<numRuptures;rup++) {
			int index = rupSamplerMFD.getClosestXIndex(rupMeanMag[rup]);
			rupSamplerMFD.add(index, rupSampler.getY(rup));
		}
		rupSamplerMFD.setName("RupSamplerMFD");
		rupSamplerMFD.normalizeByTotalRate();
		
//System.out.println("\nrupSamplerMFD:\n"+rupSamplerMFD);
//		System.out.println(rupSamplerMFD);
//		System.exit(0);
	}

	
	public void writeAndOrPlotMFDs(String dirName, boolean popupWindow, Range xAxisRange, Range yAxisRange, 
			double widthInches, double heightInches) {

		ArrayList<XY_DataSet> mfd_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> mfd_plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		// If multiple SA runs have been generated
		if(mfdsFromMultRuns != null) {
			// Incremental Distributions:
			EvenlyDiscretizedFunc meanMfdFromMultipleRuns = mfdsFromMultRuns.getMeanCurve();
			meanMfdFromMultipleRuns.setName("meanMfdFromMultipleRuns");
			UncertainArbDiscFunc mfdMinMaxRange = new UncertainArbDiscFunc(meanMfdFromMultipleRuns, 
					mfdsFromMultRuns.getMinCurve(), mfdsFromMultRuns.getMaxCurve());
			mfdMinMaxRange.setName("mfdMinMaxRange");
			mfd_funcs.add(mfdMinMaxRange);
			UncertainArbDiscFunc mfdMean95conf = get95perConfForMultRuns(mfdsFromMultRuns);
			mfdMean95conf.setName("mfdMean95conf");
			mfd_funcs.add(mfdMean95conf);
			mfd_funcs.add(meanMfdFromMultipleRuns);
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,255)));
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,255)));
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));

			// Cumulative Distributions:
			EvenlyDiscretizedFunc meanCumMfdFromMultipleRuns = cumMfdsFromMultRuns.getMeanCurve();
			meanCumMfdFromMultipleRuns.setName("meanCumMfdFromMultipleRuns");
			UncertainArbDiscFunc cumMfdMinMaxRange = new UncertainArbDiscFunc(meanCumMfdFromMultipleRuns, 
					cumMfdsFromMultRuns.getMinCurve(), cumMfdsFromMultRuns.getMaxCurve());
			cumMfdMinMaxRange.setName("cumMfdMinMaxRange");
			mfd_funcs.add(cumMfdMinMaxRange);
			UncertainArbDiscFunc cumMfdMean95conf = get95perConfForMultRuns(cumMfdsFromMultRuns);
			cumMfdMean95conf.setName("cumMfdMean95conf");
			mfd_funcs.add(cumMfdMean95conf);
			mfd_funcs.add(meanCumMfdFromMultipleRuns);
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(255,200,200)));
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(255,120,120)));
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		}
		else { // just single run value
			mfd_funcs.add(magFreqDist);
			EvenlyDiscretizedFunc cumMagFreqDist = magFreqDist.getCumRateDistWithOffset();
			cumMagFreqDist.setInfo("Cumulative Mag Freq Dist");
			mfd_funcs.add(cumMagFreqDist);
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
			mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		}
		
		// add average/single section participation MFD
		mfd_funcs.add(aveOfSectPartMFDs);
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));

		// Add target or a comparison MFD
		IncrementalMagFreqDist mfdComparison;
		EvenlyDiscretizedFunc mfdCumComparison;
		if(mfdConstraint != null) {
			mfdComparison = mfdConstraint;
			mfdComparison.setName("Target MFD Constraint Incr. Dist.");
			mfdCumComparison = mfdComparison.getCumRateDistWithOffset();
			mfdCumComparison.setName("Target Cumulative MFD Constraint");
		}
		else {
			mfdComparison = getGR_DistFit();
			mfdComparison.setName("GR Incr. Dist. with same Mmin, Mmax and moment rate, and b=1");
			mfdCumComparison = mfdComparison.getCumRateDistWithOffset();
			mfdCumComparison.setName("GR Cumulative Dist. with same Mmin, Mmax and moment rate, and b=1");
		}
		mfd_funcs.add(mfdComparison);
		mfd_funcs.add(mfdCumComparison);
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.RED));

		if(yAxisRange==null) {
			double minY=Double.POSITIVE_INFINITY, maxY=Double.NEGATIVE_INFINITY;
			for(XY_DataSet func : mfd_funcs) {
				// loop over all y values to avoid zero (which will screw up log plotting)
				for(int i=0; i<func.size();i++) {
					double val = func.getY(i);
					if(val>0) 
						if(minY>val) minY=val;
					if(maxY<val) maxY=val;

				}
			}
			if(minY==maxY) {
				minY /= 2;
				maxY *= 2;
			}
			// round minY and maxY to tick marks in log space
			double temp = Math.pow(10, Math.floor(Math.log10(minY)));
			double ratio = Math.floor(minY/temp);
			minY = temp*ratio;
			temp = Math.pow(10, Math.floor(Math.log10(maxY)));
			ratio = Math.ceil(maxY/temp);
			maxY = temp*ratio;
				
			yAxisRange = new Range(minY, maxY);
		}
		
		if(xAxisRange==null) {
			double minX=minMagMFD_WithAleatory, maxX=maxMagMFD_WithAleatory;
			if(minX==maxX) {
				minX += 2;
				maxX += 2;
			}
			xAxisRange = new Range(minX-0.05, maxX+0.05);
		}

		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/MFDs";
		String plotName = "Magnitude Frequency Distributions";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Rate (per yr)";
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(mfd_funcs, mfd_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
	}

	
	public void writeAndOrPlotNormalizedResiduals(String dirName, boolean popupWindow, double widthInches, double heightInches) {
		
		DefaultXY_DataSet normResids = new DefaultXY_DataSet();
		normResids.setName("normResidsArray");
		normResids.setInfo("(p-o)/std for each row (data constraint), where p = predicted value, o = obs value, and std is the standard deviation of the observation (assuming wt = 1/std)");
		double nr_min = Double.POSITIVE_INFINITY;
		double nr_max = Double.NEGATIVE_INFINITY;
		for(int row=0; row < totNumRows; row++) {
			double nr = (d_pred[row]-d[row])*data_wt[row];
			normResids.set((double)row,nr);
			if(nr<nr_min)
				nr_min = nr;
			if(nr>nr_max)
				nr_max = nr;
		}
		//System.out.println("nr_min="+nr_min+"\nnr_max="+nr_max);
		
		ArrayList<XY_DataSet> nr_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> nr_plotChars = new ArrayList<PlotCurveCharacterstics>();
		nr_funcs.add(normResids);
		nr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 2f, Color.BLUE));
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/normalizedResiduals";

		String plotName = "";
		String xAxisLabel = "Data Constraint Index";
		String yAxisLabel = "Normalized Residuals";
		if(nr_max<0) nr_max *= -1;
		if(nr_min>0) nr_min *= -1;
		double maxNR = Math.max(-nr_min, nr_max);
		Range xAxisRange = new Range(-1,totNumRows);
		Range yAxisRange = new Range(-maxNR,maxNR);
		boolean logX = false;
		boolean logY = false;
		
		PlottingUtils.writeAndOrPlotFuncs(nr_funcs, nr_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}
	
	public PlotSpec writeAndOrPlotSlipRates(String dirName, boolean popupWindow, Range xAxisRange, Range yAxisRange, 
			double widthInches, double heightInches) {

		double min = 0, max = numSections-1;
		EvenlyDiscretizedFunc origSlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc origUpper68_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc origLower68_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc finalSlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		ArrayList<DefaultXY_DataSet> slipRate95confFuncsList = new ArrayList<DefaultXY_DataSet>();
		DefaultXY_DataSet slipRate95confFunc;
		for(int s=0; s<numSections;s++) {
			origSlipRateFunc.set(s,sectSlipRate[s]*(1-moRateReduction));
			slipRate95confFunc = new DefaultXY_DataSet();
			slipRate95confFunc.set((double)s,(sectSlipRate[s]-1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
			slipRate95confFunc.set((double)s,(sectSlipRate[s]+1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
			slipRate95confFunc.setName("95% conf for orig slip rate on section "+s);
			slipRate95confFuncsList.add(slipRate95confFunc);
			origUpper68_SlipRateFunc.set(s,(sectSlipRate[s]+sectSlipRateStdDev[s])*(1-moRateReduction));
			origLower68_SlipRateFunc.set(s,(sectSlipRate[s]-sectSlipRateStdDev[s])*(1-moRateReduction));
			finalSlipRateFunc.set(s,finalSectSlipRate[s]);
		}
		origSlipRateFunc.setName("Target Slip Rates (mm/yr)");
		origUpper68_SlipRateFunc.setName("Target Slip Rate Upper (+stdev)");
		origLower68_SlipRateFunc.setName("Target Slip Rate Lower (-stdev)");
		finalSlipRateFunc.setName("Final Slip Rates (mm/yr)");
		
		// convert to mm/yr:
		origSlipRateFunc.scale(1000);;
		origUpper68_SlipRateFunc.scale(1000);
		origLower68_SlipRateFunc.scale(1000);
		finalSlipRateFunc.scale(1000);
	

		ArrayList<XY_DataSet> sr_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> sr_plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		if(finalSectSlipRateFromMultRuns != null) {
			EvenlyDiscretizedFunc meanSlipRateFromMultipleRuns = finalSectSlipRateFromMultRuns.getMeanCurve();
			meanSlipRateFromMultipleRuns.setName("meanSlipRateFromMultipleRuns");
			meanSlipRateFromMultipleRuns.scale(1000);	// convert to mm/yr
			
			DiscretizedFunc minCurve = finalSectSlipRateFromMultRuns.getMinCurve();;
			DiscretizedFunc maxCurve = finalSectSlipRateFromMultRuns.getMaxCurve();
			minCurve.scale(1000);	// convert to mm/yr
			maxCurve.scale(1000);	// convert to mm/yr
			UncertainArbDiscFunc slipRatesMinMaxRange = new UncertainArbDiscFunc(meanSlipRateFromMultipleRuns, minCurve, maxCurve);
			slipRatesMinMaxRange.setName("slipRatesMinMaxRange");
			sr_funcs.add(slipRatesMinMaxRange);
			
//			UncertainArbDiscDataset slipRatesStdDevRange = new UncertainArbDiscDataset(meanSlipRateFromMultipleRuns, 
//					finalSectSlipRateFromMultRuns.getMeanPlusXstdDevCurve(-1.0), finalSectSlipRateFromMultRuns.getMeanPlusXstdDevCurve(1.0));
//			slipRatesStdDevRange.setName("slipRatesPlusMinus1stdDevRange");
//			sr_funcs.add(slipRatesStdDevRange);


			// need the following to scale to mm/yr
			EvenlyDiscretizedFunc meanCurve = finalSectSlipRateFromMultRuns.getMeanCurve();
			EvenlyDiscretizedFunc stdevCurve = finalSectSlipRateFromMultRuns.getStdDevCurve();
			EvenlyDiscretizedFunc upper95 = stdevCurve.deepClone();
			EvenlyDiscretizedFunc lower95 = stdevCurve.deepClone();
			double numRuns = finalSectSlipRateFromMultRuns.getArbDiscrEmpDistFuncArray()[0].calcSumOfY_Vals();
			double sqrtNum = Math.sqrt(numRuns);
			for(int i=0;i<meanCurve.size();i++) {
				double mean = meanCurve.getY(i);
				double stdom = stdevCurve.getY(i)/sqrtNum;
				upper95.set(i,mean+1.96*stdom);
				lower95.set(i,mean-1.96*stdom);
			}
			meanCurve.scale(1000);
			lower95.scale(1000);
			upper95.scale(1000);
			UncertainArbDiscFunc slipRatesMean95conf = new UncertainArbDiscFunc(meanCurve,lower95,upper95);

			// this cannot be scaled to mm/yr:
//			UncertainArbDiscDataset slipRatesMean95conf = get95perConfForMultRuns(finalSectSlipRateFromMultRuns);
			slipRatesMean95conf.setName("slipRatesMean95conf");
			sr_funcs.add(slipRatesMean95conf);
			sr_funcs.add(meanSlipRateFromMultipleRuns);
	
			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,255)));
			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,255)));
//			// if stdDev included:
//			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(220,220,255)));
//			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(150,150,255)));
//			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(80,80,255)));
			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		}
		else {
			sr_funcs.add(finalSlipRateFunc);
			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		}

		sr_funcs.add(origSlipRateFunc);
		sr_funcs.add(origUpper68_SlipRateFunc);
		sr_funcs.add(origLower68_SlipRateFunc);
//		sr_funcs.addAll(slipRate95confFuncsList);
		
//		sr_plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLUE));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.BLACK));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 0.5f, Color.BLACK));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 0.5f, Color.BLACK));
//		for(int i=0;i<slipRate95confFuncsList.size();i++) {
//			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 4f, Color.BLUE));
//		}
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sectionSlipRates";

		String plotName = "";
		String xAxisLabel = "Subsection";
		String yAxisLabel = "Slip Rate (mm/yr)";
		boolean logX = false;
		boolean logY = false;
		
		return PlottingUtils.writeAndOrPlotFuncs(sr_funcs, sr_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
	}
	
	
	public PlotSpec writeAndOrPlotEventRates(String dirName, boolean popupWindow, Range xAxisRange, Range yAxisRange, 
			double widthInches, double heightInches) {

		double min = 0, max = numSections-1;
		ArrayList<XY_DataSet> er_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> er_plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		EvenlyDiscretizedFunc finalEventRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc finalPaleoVisibleEventRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		for(int s=0;s < numSections; s++) {
			finalPaleoVisibleEventRateFunc.set(s, finalPaleoVisibleSectEventRate[s]);
			finalEventRateFunc.set(s,finalSectEventRate[s]);
		}
		finalPaleoVisibleEventRateFunc.setName("Final Paleoseismically Visible Event Rates");
		finalEventRateFunc.setName("Final Event Rates (dashed)");


		if(finalPaleoVisibleSectEventRateFromMultRuns != null) {
			EvenlyDiscretizedFunc meanPaleoVisEventRateFromMultipleRuns = finalPaleoVisibleSectEventRateFromMultRuns.getMeanCurve();
			meanPaleoVisEventRateFromMultipleRuns.setName("meanPaleoVisEventRateFromMultipleRuns");

			UncertainArbDiscFunc paleoVisSlipRatesMinMaxRange = new UncertainArbDiscFunc(meanPaleoVisEventRateFromMultipleRuns, 
					finalPaleoVisibleSectEventRateFromMultRuns.getMinCurve(), finalPaleoVisibleSectEventRateFromMultRuns.getMaxCurve());
			paleoVisSlipRatesMinMaxRange.setName("paleoVisSlipRatesMinMaxRange");
			er_funcs.add(paleoVisSlipRatesMinMaxRange);
			
			UncertainArbDiscFunc paleoVisEventRatesMean95conf = get95perConfForMultRuns(finalPaleoVisibleSectEventRateFromMultRuns);
			paleoVisEventRatesMean95conf.setName("paleoVisEventRatesMean95conf");
			er_funcs.add(paleoVisEventRatesMean95conf);
			
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(255,200,200)));
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(255,120,120)));
			er_funcs.add(finalEventRateFunc);
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
			er_funcs.add(meanPaleoVisEventRateFromMultipleRuns);
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		}
		else {
			er_funcs.add(finalEventRateFunc);
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
			er_funcs.add(finalPaleoVisibleEventRateFunc);
			er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		}
		
		int num = sectionRateConstraints.size();
		ArrayList<AbstractXY_DataSet> obs_er_funcs = new ArrayList<AbstractXY_DataSet>();

		if(num>0) {
			DefaultXY_DataSet func;
			ArbitrarilyDiscretizedFunc meanER_Func = new ArbitrarilyDiscretizedFunc();
			er_funcs.add(meanER_Func);
			obs_er_funcs.add(meanER_Func);
			SectionRateConstraint constraint;
			for(int c=0;c<num;c++) {
				func = new DefaultXY_DataSet();
				constraint = sectionRateConstraints.get(c);
				int s = constraint.getSectIndex();
				meanER_Func.set((double)s, constraint.getMeanRate());
				func.set((double)s, constraint.getLower95Conf());
				func.set((double)s, constraint.getUpper95Conf());
				func.setName(constraint.getFaultName());
				obs_er_funcs.add(func);
				er_funcs.add(func);
			}			
		}
		er_plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
		for(int c=0;c<num;c++)
			er_plotChars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 4f, Color.RED));

		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sectionEventRates";
		String plotName = null;
		String xAxisLabel = "Subsection";
		String yAxisLabel = "Event Rate (per yr)";
		boolean logX = false;
		boolean logY = false;

		return PlottingUtils.writeAndOrPlotFuncs(er_funcs, er_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
	}


	
	public PlotSpec writeAndOrPlotAveSlipAndCOV(String dirName, boolean popupWindow, Range xAxisRange, Range yAxisRange, 
			double widthInches, double heightInches) {

		double min = 0, max = numSections-1;

		ArrayList<XY_DataSet> slip_mean_cov_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> slip_mean_cov_plotChars = new ArrayList<PlotCurveCharacterstics>();
		
		if(finalSectMeanSlipFromMultRuns != null) {
			EvenlyDiscretizedFunc meanSlipFromMultipleRuns = finalSectMeanSlipFromMultRuns.getMeanCurve();
			meanSlipFromMultipleRuns.setName("meanSlipFromMultipleRuns");
			UncertainArbDiscFunc meanSlipMinMaxRange = new UncertainArbDiscFunc(meanSlipFromMultipleRuns, 
					finalSectMeanSlipFromMultRuns.getMinCurve(), finalSectMeanSlipFromMultRuns.getMaxCurve());
			meanSlipMinMaxRange.setName("meanSlipMinMaxRange");
			slip_mean_cov_funcs.add(meanSlipMinMaxRange);
			
			UncertainArbDiscFunc meanSlip95conf = get95perConfForMultRuns(finalSectMeanSlipFromMultRuns);
			meanSlip95conf.setName("meanSlip95conf");
			slip_mean_cov_funcs.add(meanSlip95conf);
			slip_mean_cov_funcs.add(meanSlipFromMultipleRuns);
	
//			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,255)));
//			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,255)));
//			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(225,195,225)));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(180,128,180)));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(128,0,128)));
			
			
			EvenlyDiscretizedFunc slipCOV_FromMultipleRuns = finalSectSlipCOV_FromMultRuns.getMeanCurve();
			slipCOV_FromMultipleRuns.setName("slipCOV_FromMultipleRuns");
			UncertainArbDiscFunc slipCOV_MinMaxRange = new UncertainArbDiscFunc(slipCOV_FromMultipleRuns, 
					finalSectSlipCOV_FromMultRuns.getMinCurve(), finalSectSlipCOV_FromMultRuns.getMaxCurve());
			slipCOV_MinMaxRange.setName("slipCOV_MinMaxRange");
			slip_mean_cov_funcs.add(slipCOV_MinMaxRange);
			
			UncertainArbDiscFunc slipCOV_95conf = get95perConfForMultRuns(finalSectSlipCOV_FromMultRuns);
			slipCOV_95conf.setName("slipCOV_95conf");
			slip_mean_cov_funcs.add(slipCOV_95conf);
			slip_mean_cov_funcs.add(slipCOV_FromMultipleRuns);
	
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,255,200)));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,255,120)));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		}
		else {
			EvenlyDiscretizedFunc finalMeanSlipFunc = new EvenlyDiscretizedFunc(min, max, numSections);
			EvenlyDiscretizedFunc finalSlipCOVFunc = new EvenlyDiscretizedFunc(min, max, numSections);
			for(int s=0; s<numSections;s++) {
				finalMeanSlipFunc.set(s, finalSectMeanSlip[s]);
				finalSlipCOVFunc.set(s,finalSectSlipCOV[s]);
			}
			finalMeanSlipFunc.setName("Mean Slip");
			finalSlipCOVFunc.setName("Slip COV");
			slip_mean_cov_funcs.add(finalMeanSlipFunc);
			slip_mean_cov_funcs.add(finalSlipCOVFunc);
			
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(128,0,128)));
			slip_mean_cov_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		}
		

		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sectionMeanSlipAndCOV";

		String plotName = null;
		String xAxisLabel = "Subsection";
		String yAxisLabel = "Slip (m)";
		boolean logX = false;
		boolean logY = false;
		
		return PlottingUtils.writeAndOrPlotFuncs(slip_mean_cov_funcs, slip_mean_cov_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
	}


	public void writeAndOrPlotRupRateVsIndex(String dirName, boolean popupWindow, double widthInches, double heightInches) {
		
		double min = -1;
		double max = numRuptures-1;

		ArrayList<XY_DataSet> rup_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> rup_plotChars = new ArrayList<PlotCurveCharacterstics>();

		if(rupRatesFromMultRuns != null) {
			UncertainArbDiscFunc rupRatesMinMaxRange = new UncertainArbDiscFunc(rupRatesFromMultRuns.getMeanCurve(), 
					rupRatesFromMultRuns.getMinCurve(), rupRatesFromMultRuns.getMaxCurve());
			rupRatesMinMaxRange.setName("rupRatesMinMaxRange");
			rup_funcs.add(rupRatesMinMaxRange);
			
			UncertainArbDiscFunc rupRatesMean95conf = get95perConfForMultRuns(rupRatesFromMultRuns);
			rupRatesMean95conf.setName("rupRatesMean95conf");
			rup_funcs.add(rupRatesMean95conf);
			
			rup_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,200)));
			rup_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,120)));

		}
		
		EvenlyDiscretizedFunc rupRateFunc = new EvenlyDiscretizedFunc(min, max, numRuptures);
		for(int rup=0; rup<numRuptures;rup++) {
			rupRateFunc.set(rup,rupRateSolution[rup]);
		}
		rupRateFunc.setName("Rupture Rates");
		rup_funcs.add(rupRateFunc);
		rup_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/ruptureRatesVsIndex";
		String plotName = "Rupture Rates";
		String xAxisLabel = "Rupture Index";
		String yAxisLabel = "Rate (per yr)";
		Range xAxisRange = new Range(-1, numRuptures);
		Range yAxisRange = new Range(1e-10, rupRateFunc.getMaxY());
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(rup_funcs, rup_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);

	}
	
	/**
	 * This only writes out the single or average solution (multiple SA solutions are
	 * written out in the simulation loop in case of crash)
	 * @param fileName
	 */
	public void writeRuptureRatesToFile(String dirName) {
		// write alternative rupture rates file to get rid of big header
		try{
			FileWriter fw = new FileWriter(dirName+ "/ruptureRates.txt");	// fileNamePrefix+"Alt.txt"
			for(int i=0;i<numRuptures; i++) {				
				fw.write(i+"\t"+rupRateSolution[i]+"\n");
			}
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}			
	}

	
		
	public void writeAndOrPlotSectionBoundaryRates(String dirName, boolean popupWindow, Range yAxisRange, 
			double widthInches, double heightInches) {

		EvenlyDiscretizedFunc rateOfThroughgoingRupsFunc = new EvenlyDiscretizedFunc(-0.5, numSections+1, 1.0);
		for(int s=0; s<numSections+1;s++) {
			rateOfThroughgoingRupsFunc.set(s,rateOfThroughGoingRupsAtSectBoudary[s]);
		}
		rateOfThroughgoingRupsFunc.setName("Rate of throughgoing ruptures at each section boundary");

		ArrayList<XY_DataSet> sect_funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();

		if(rateOfThroughGoingRupsAtSectBoudaryFromMultRuns != null) {
			EvenlyDiscretizedFunc meanRateOfThroughGoingRupsFromMultipleRuns = rateOfThroughGoingRupsAtSectBoudaryFromMultRuns.getMeanCurve();
			meanRateOfThroughGoingRupsFromMultipleRuns.setName("meanRateOfThroughGoingRupsFromMultipleRuns");

			UncertainArbDiscFunc rateOfThroughGoingRupsMinMaxRange = new UncertainArbDiscFunc(meanRateOfThroughGoingRupsFromMultipleRuns, 
					rateOfThroughGoingRupsAtSectBoudaryFromMultRuns.getMinCurve(), rateOfThroughGoingRupsAtSectBoudaryFromMultRuns.getMaxCurve());
			rateOfThroughGoingRupsMinMaxRange.setName("rateOfThroughGoingRupsMinMaxRange");
			sect_funcs.add(rateOfThroughGoingRupsMinMaxRange);

			UncertainArbDiscFunc rateOfThroughGoingRupsMean95conf = get95perConfForMultRuns(rateOfThroughGoingRupsAtSectBoudaryFromMultRuns);
			rateOfThroughGoingRupsMean95conf.setName("rateOfThroughGoingRupsMean95conf");
			sect_funcs.add(rateOfThroughGoingRupsMean95conf);

			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(200,200,200)));
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120,120,120)));
		}

		sect_funcs.add(rateOfThroughgoingRupsFunc);
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sectionBoundaryRates";
		String plotName =  "Subsection Boundary Rates";
		String xAxisLabel = "Subaection Boundary";
		String yAxisLabel = "Rate (per yr)";
		Range xAxisRange = new Range(-1,numSections);
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(sect_funcs, plotChars2, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);


	}

	
	public 	UncertainArbDiscFunc get95perConfForMultRuns(ArbDiscrEmpiricalDistFunc_3D arbDiscrEmpiricalDistFunc_3D) {
		EvenlyDiscretizedFunc meanCurve = arbDiscrEmpiricalDistFunc_3D.getMeanCurve();
		EvenlyDiscretizedFunc stdevCurve = arbDiscrEmpiricalDistFunc_3D.getStdDevCurve();
		EvenlyDiscretizedFunc upper95 = stdevCurve.deepClone();
		EvenlyDiscretizedFunc lower95 = stdevCurve.deepClone();
		double numRuns = arbDiscrEmpiricalDistFunc_3D.getArbDiscrEmpDistFuncArray()[0].calcSumOfY_Vals();
		
//		ArbDiscrEmpiricalDistFunc[] funcArray = arbDiscrEmpiricalDistFunc_3D.getArbDiscrEmpDistFuncArray();
//		for(int i=0;i<funcArray.length;i++) {
//			ArbDiscrEmpiricalDistFunc func = funcArray[i];
//			System.out.println("HERE: "+func.size()+"\t"+func.calcSumOfY_Vals());
//		}
		
		double sqrtNum = Math.sqrt(numRuns);
		for(int i=0;i<meanCurve.size();i++) {
			double mean = meanCurve.getY(i);
			double stdom = stdevCurve.getY(i)/sqrtNum;
			upper95.set(i,mean+1.96*stdom);
			lower95.set(i,mean-1.96*stdom);
		}
		return new UncertainArbDiscFunc(meanCurve,lower95,upper95);
	}
	
	

	

	/**
	 * This makes the following histograms: 1) the number of ruptures having N subsection, the number
	 * of ruptures in each magnitude bin, the default likelihood that a subsection will be included in
	 * a randomly chosen rupture, and the ruptureSampler implied MFD (if ruptureSampler != null)
	 * @param dirName - set as null if you don't want to save results
	 * @param popupWindow - set as true if you want plot windows to pop up
	 */
	public void writeAndOrPlotMagHistograms(String dirName, boolean popupWindow, double widthInches, double heightInches) {
		
		ArrayList<XY_DataSet> funcs1 = new ArrayList<XY_DataSet>();
		funcs1.add(meanMagHistorgram);
		funcs1.add(magHistorgram);

		
		// make a numSegInRupHistogram
		HistogramFunction numSectInRupHistogram = new HistogramFunction(1.0,numSections,1.0);
		for(int r=0;r<numSectInRup.length;r++) numSectInRupHistogram.add((double)numSectInRup[r], 1.0);
		numSectInRupHistogram.setName("Num Segments In Rupture Histogram");
		ArrayList<XY_DataSet> funcs2 = new ArrayList<XY_DataSet>();
		funcs2.add(numSectInRupHistogram);	
		
		
		// make prob of section being randomly selected
		HistogramFunction probSectSelectedHistogram = new HistogramFunction(0.0,numSections,1.0);
		for(int r=0;r<numRuptures;r++) {
			for(int s=0; s<numSections; s++) {
				if(rupSectionMatrix[s][r] == 1)
					probSectSelectedHistogram.add(s, 1);
			}
		}
//		probSectSelectedHistogram.normalizeBySumOfY_Vals();
		probSectSelectedHistogram.setName("Prob Section Randomly Selected Histogram (or num rups that touch sect)");
		ArrayList<XY_DataSet> funcs3 = new ArrayList<XY_DataSet>();
		funcs3.add(probSectSelectedHistogram);

		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.gray));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		double minX=magHistorgram.getMinX();
		double maxX=magHistorgram.getMaxX();
		if(minX==maxX) {
			minX += 1;
			maxX -= 1;
		}

		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/magHistograms";
		String plotName ="Mag Histograms";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Num Ruptures";
		Range xAxisRange = new Range(minX-MAG_DELTA/2.0,maxX+MAG_DELTA/2.0);
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(funcs1, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/numSectInRuptHistogram";
		plotName = "Num Sect In Rup Histograms";
		xAxisLabel = "Num Subsections";
		yAxisLabel = "Num Ruptures";
		xAxisRange = new Range(-1,numSections);
		yAxisRange = null;
		logX = false;
		logY = false;

		PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/probSegSelectedHistogram";
		plotName = "Prob Section Randomly Selected Histogram";
		xAxisLabel = "Subsection ID";
		yAxisLabel = "Num Ruptures";

		PlottingUtils.writeAndOrPlotFuncs(funcs3, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);
		
		
		if(rupSamplerMFD != null) {
			// plot the SA rupture sampler
			ArrayList<XY_DataSet> funcs4 = new ArrayList<XY_DataSet>();
			funcs4.add(rupSamplerMFD);
			ArrayList<PlotCurveCharacterstics> plotChars4 = new ArrayList<PlotCurveCharacterstics>();
			plotChars4.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));

			fileNamePrefix = null;
			if(dirName != null)
				fileNamePrefix = dirName+"/rupSamplerMFD";
			plotName ="SA RupSamplerMFD";
			xAxisLabel = "Magnitude";
			yAxisLabel = "PDF";
			xAxisRange = new Range(minX-0.05,maxX+0.05);
			yAxisRange = null;
			logX = false;
			logY = true;

			PlottingUtils.writeAndOrPlotFuncs(funcs4, plotChars4, plotName, xAxisLabel, yAxisLabel, 
					xAxisRange, yAxisRange, logX, logY, widthInches, heightInches, fileNamePrefix, popupWindow);	
		
		}
	}
	
	

	
	
	private GutenbergRichterMagFreqDist getGR_DistFit() {
		int num = (int)Math.round((maxMagMFD-minMagMFD)/0.1 + 1);
		if(num==1) num=2;
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(minMagMFD,num,0.1);
		double moRate = totMoRate*(1-moRateReduction);
//		double altMoRate = magFreqDist.getTotalMomentRate();
		gr.setAllButTotCumRate(minMagMFD, maxMagMFD, moRate, 1.0);
		gr.setName("GR fit");
//		if(D) {
//			System.out.println("GR MFD FIT:\n\tMmin="+minRupMag+"\n\tMmax="+maxRupMag+
//					"\n\tnum="+num+"\n\tdelta="+gr.getDelta()+"\n\tmoRate="+moRate);
//		}
		return gr;
	}
	
	
	/**
	 * This writes a-priori rupture rates computed from the MFD constraint to the
	 * specified file.
	 * @param fileName
	 */
	public void writeAprioriRupRatesFromMFD_Constraint(String fileName) {
		try{
			double aPriori_rate, aPriori_wt;
			FileWriter fw = new FileWriter(fileName);
			for(int r=0;r<numRuptures;r++) {
				aPriori_rate = mfdConstraint.getY(rupMeanMag[r])/meanMagHistorgram.getY(rupMeanMag[r]);
				aPriori_wt = 1;
				fw.write(r+"\t"+aPriori_rate+"\t"+aPriori_wt+"\n");
			}
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}			
	}
	
	/**
	 * This sets a-priori rupture rates from the MFD constraint.  Each rate is given a weight of 1.0, even if mfdSigma is non null.
	 */
	public void setAprioriRupRatesFromMFD_Constrint() {
		aPriori_rupIndex = new int[numRuptures];
		aPriori_rate  = new double[numRuptures];
		aPriori_wt  = new double[numRuptures];
		for(int r=0;r<numRuptures;r++) {
			aPriori_rupIndex[r] = r;
			aPriori_rate[r] = mfdConstraint.getClosestYtoX(rupMeanMag[r])/meanMagHistorgram.getClosestYtoX(rupMeanMag[r]);
			aPriori_wt[r] = 1;
		}
	}

	
	/**
	 * This reads the a-priori rupture rates from a file, where the first column is
	 * rupture index, the second column is rate, and third is weight
	 * @param fileName
	 */
	private void readApriorRupRateConstraintsFromFile(String fileName) {
		try {
			File file = new File(fileName);
			List<String> fileLines;
			fileLines = Files.readLines(file, Charset.defaultCharset());
			int numLines = fileLines.size();
			aPriori_rupIndex = new int[numLines];
			aPriori_rate  = new double[numLines];
			aPriori_wt  = new double[numLines];
			int index=0;
			for (String line : fileLines) {
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 3, "Expected 3 items, got %s", split.length);
				int rupIndex = Integer.valueOf(split[0]);
				aPriori_rupIndex[index] = rupIndex;
				aPriori_rate[index] = Double.valueOf(split[1]);	
				aPriori_wt[index] = Double.valueOf(split[2]);
				index+=1;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public double[] getAprioriRuptureRates() {
		return aPriori_rate;
	}
	
	
	/**
	 * This returns the probability that the given magnitude 
	 * event will be observed at the ground surface.  This is based
	 * on equation 4 of Youngs et al. [2003, A Methodology for Probabilistic Fault
	 * Displacement Hazard Analysis (PFDHA), Earthquake Spectra 19, 191-219] using 
	 * the coefficients they list in their appendix for "Data from Wells and 
	 * Coppersmith (1993) 276 worldwide earthquakes".  Their function has the following
	 * probabilities:
	 * 
	 * mag	prob
	 * 5		0.10
	 * 6		0.45
	 * 7		0.87
	 * 8		0.98
	 * 9		1.00
	 * 
	 * @return
	 */
	private double getProbVisible(double mag) {
		return Math.exp(-12.51+mag*2.053)/(1.0 + Math.exp(-12.51+mag*2.053));
		/* Ray & Glenn's equation
		 if(mag <= 5) 
		 return 0.0;
		 else if (mag <= 7.6)
		 return -0.0608*mag*mag + 1.1366*mag + -4.1314;
		 else 
		 return 1.0;
		 */
	}
	
	/**
	 * This initializes the section and rupture attributes
	 */
	private void initSectAndRupAttributes() {
		
		// compute slip correction for Gaussian MFD
		GaussianMagFreqDist gDist1 = new GaussianMagFreqDist(5.0,9.0,41,7,GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2);
		GaussianMagFreqDist gDist2 = new GaussianMagFreqDist(5.0,9.0,41,7,0.01,1.0,0.01,2);
		gDist1.scaleToCumRate(0, 1.0);
		gDist2.scaleToCumRate(0, 1.0);
		gaussMFD_slipCorr = gDist1.getTotalMomentRate()/gDist2.getTotalMomentRate();
		System.out.println("gaussMFD_slipCorr="+(float)gaussMFD_slipCorr+"\n");

		// set section data values
		numSections = fltSectionDataList.size();
		if(D) System.out.println("numSections="+numSections);
		modelSetUpInfoString += "\nnumSections = "+numSections+"\n";
		
		sectSlipRate = new double[numSections];
		sectSlipRateStdDev = new double[numSections];
		sectLength = new double[numSections];
		sectArea = new double[numSections];
		sectMoRate = new double[numSections];
		sectRake = new double[numSections];
		totMoRate=0;
		for(int s=0;s<numSections;s++) {
			FaultSectionPrefData fltSectData = fltSectionDataList.get(s);
			// Note that moRateReduction is not applied in the following, but later when setting up matrix
			sectSlipRate[s] = fltSectData.getReducedAveSlipRate()*1e-3;	// convert to meters/yr
			sectSlipRateStdDev[s] = fltSectData.getReducedSlipRateStdDev()*1e-3;	// convert to meters/yr
			sectLength[s] = fltSectData.getTraceLength()*1e3; // convert to meters
			sectArea[s] = sectLength[s] * fltSectData.getReducedDownDipWidth()*1e3; // convert latter to meters
			sectMoRate[s] = fltSectData.calcMomentRate(true);
			sectRake[s] = fltSectData.getAveRake();
			totMoRate += sectMoRate[s];
		}

		// set rupture attributes
		numRuptures = rupSectionMatrix[0].length;
		if(D) System.out.println("numRuptures="+numRuptures+"\n");
		modelSetUpInfoString += "\nnumRuptures = "+numRuptures+"\n";


		rupNameShort = new String[numRuptures];
		rupLength = new double[numRuptures];
		rupArea = new double[numRuptures];
		rupMoRateMax = new double[numRuptures];
		rupAveRake = new double[numRuptures];
		numSectInRup = new int[numRuptures];
		firstSectOfRup = new int[numRuptures];
		minNumSectInRup = Integer.MAX_VALUE;

		for(int r=0;r<numRuptures;r++) {
			boolean foundFirstSect = false;
			for(int s=0; s<numSections; s++) {
				if(rupSectionMatrix[s][r] == 1) {
					rupLength[r] += sectLength[s];
					rupArea[r] += sectArea[s];
					rupMoRateMax[r] += sectMoRate[s];
					rupAveRake[r] += sectRake[s];
					numSectInRup[r] += 1;
					if(!foundFirstSect) {
						firstSectOfRup[r] = s;
						rupNameShort[r] = ""+s;
						foundFirstSect=true;
					}
					else {
						rupNameShort[r] += "+"+s;
					}
				}
				if(minNumSectInRup>numSectInRup[r])
					minNumSectInRup=numSectInRup[r];
			}
			rupAveRake[r] /= rupArea[r];
		}

		if(D) System.out.println("minNumSectInRup="+minNumSectInRup);
		modelSetUpInfoString += "\nminNumSectInRup = "+minNumSectInRup+"\n\n"+
				"gaussMFD_slipCorr="+(float)gaussMFD_slipCorr+"\n";

		// compute rupture mean mags etc
		rupMeanMag = new double[numRuptures];
		rupAveSlip = new double[numRuptures];
		rupMeanMo = new double[numRuptures];
		minRupMag = Double.MAX_VALUE;
		maxRupMag = 0;
		
		double origWidth = fltSectionDataList.get(0).getOrigDownDipWidth()*1e3;	// assuming this is the same for all sections
		for(int r=0; r <numRuptures; r++) {
			double mag = magAreaRel.getMag(rupArea[r], rupLength[r], origWidth);	
			
			//round this to nearest 100th unit
			rupMeanMag[r] = ((double)Math.round(100*mag))/100.0;
			rupMeanMo[r] = MagUtils.magToMoment(rupMeanMag[r])*gaussMFD_slipCorr;   // increased if magSigma >0
//			rupAveSlip[r] = rupMeanMo[r]/(rupArea[r]*FaultMomentCalc.SHEAR_MODULUS);  // inlcudes aveSlipCor in rupMeanMo
			rupAveSlip[r] = magAreaRel.getAveSlip(rupArea[r], rupLength[r], rupArea[r]/rupLength[r])*gaussMFD_slipCorr;
			
			if(minRupMag>rupMeanMag[r])
				minRupMag=rupMeanMag[r];
			if(maxRupMag<rupMeanMag[r])
				maxRupMag=rupMeanMag[r];
		}
		
		double tempMag = minRupMag-GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		minRupMagWithAleatory = ((double)Math.round(100*tempMag))/100.0;
		tempMag = maxRupMag+GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		maxRupMagWithAleatory = ((double)Math.round(100*tempMag))/100.0;
		String tempString = "minRupMag = "+minRupMag+"; maxRupMag = "+maxRupMag+
				"; minRupMagWithAleatory = "+minRupMagWithAleatory+
				"; maxRupMagWithAleatory = "+maxRupMagWithAleatory;
		
		minMagMFD = roundMagTo10thUnit(minRupMag);
		maxMagMFD = roundMagTo10thUnit(maxRupMag);
		minMagMFD_WithAleatory = roundMagTo10thUnit(minRupMagWithAleatory);
		maxMagMFD_WithAleatory = roundMagTo10thUnit(maxRupMagWithAleatory);
		
		tempString += "\nminMagMFD = "+minMagMFD+"; maxMagMFD = "+maxMagMFD+
				"; minMagMFD_WithAleatory = "+minMagMFD_WithAleatory+
				"; maxMagMFD_WithAleatory = "+maxMagMFD_WithAleatory;

		
		if(D) System.out.println(tempString);
		modelSetUpInfoString += "\n"+tempString+"\n";

		// compute meanMagHistorgram
		int num = (int)Math.round((maxMagMFD-minMagMFD)/MAG_DELTA + 1);
		if(num==1) num=2;
		meanMagHistorgram = new SummedMagFreqDist(minMagMFD,num,MAG_DELTA);
		String[] numSectInBinArray = new String[meanMagHistorgram.size()];
		ArrayList<Integer> numBeenDone = new ArrayList<Integer>();
		for(int i=0;i<numSectInBinArray.length;i++)
			numSectInBinArray[i] = (float)meanMagHistorgram.getX(i)+": ";
		for(int rup=0; rup<numRuptures;rup++) {
			int index = meanMagHistorgram.getClosestXIndex(rupMeanMag[rup]);
			meanMagHistorgram.add(index, 1.0);
			int numSect = numSectInRup[rup];
			if(!numBeenDone.contains(numSect)) {
				numSectInBinArray[index] += ", "+numSect;
				numBeenDone.add(numSect);
			}
//			meanMagHistorgram.addResampledMagRate(rupMeanMag[rup], 1.0, true);
		}
		String info = "Mean Mag Histogram\n\nNum Sect in Bin:\n\n";
		for(String st:numSectInBinArray)
			info += st+"\n";
		meanMagHistorgram.setInfo(info);
		
		// compute mag historgram with aleatory variability
		num = (int)Math.round((maxMagMFD_WithAleatory-minMagMFD_WithAleatory)/MAG_DELTA + 1);
		if(num==1) num=2;
		magHistorgram = new SummedMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA);
		for(int rup=0; rup<numRuptures;rup++) {
			if(GAUSS_MFD_SIGMA == 0d) {
				int index = magHistorgram.getClosestXIndex(rupMeanMag[rup]);
				magHistorgram.add(index, 1.0);
//				magHistorgram.addResampledMagRate(rupMeanMag[rup], 1.0, true);
			}
			else {
				GaussianMagFreqDist gDist = new GaussianMagFreqDist(minMagMFD_WithAleatory,num,MAG_DELTA,rupMeanMag[rup],GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2);
				gDist.scaleToCumRate(0, 1.0); // this makes it a PDF
				magHistorgram.addIncrementalMagFreqDist(gDist);				
			}
		}
		magHistorgram.setInfo("Mag Historgram (including aleatory variability)");
	}
	
	
	public static double roundMagTo10thUnit(double mag) {
		return ((double)Math.round(10*(mag-0.05)))/10.0 +0.05;
	}
	

	
	
	// The returns the fault system solution for current solution (or the mean from multiple solutions)
	public U3FaultSystemSolution getFaultSystemSolution() {
		
		List<List<Integer>> sectionsForRups = Lists.newArrayList();;
		for(int r=0;r<numRuptures;r++) {
			ArrayList<Integer> sectList = new ArrayList<Integer>();
			for(int s=0;s<numSections;s++) {
				if(rupSectionMatrix[s][r] == 1)
					sectList.add(s);
			}
			sectionsForRups.add(sectList);
		}
		
		U3FaultSystemRupSet rupSet = new U3FaultSystemRupSet(
				fltSectionDataList,
				sectSlipRate,
				sectSlipRateStdDev,
				sectArea,
				sectionsForRups,
				rupMeanMag,
				rupAveRake,
				rupArea,
				rupLength,
				modelName);
		
		return new U3FaultSystemSolution(rupSet, rupRateSolution);
		
	}
	
	public U3FaultSystemSolution getFaultSystemSolution(int solutionIndex) {
		
		if(rupRatesFromMultRunsArrayList == null)
			throw new RuntimeException("Multiple solutions do not exist");
		
		List<List<Integer>> sectionsForRups = Lists.newArrayList();;
		for(int r=0;r<numRuptures;r++) {
			ArrayList<Integer> sectList = new ArrayList<Integer>();
			for(int s=0;s<numSections;s++) {
				if(rupSectionMatrix[s][r] == 1)
					sectList.add(s);
			}
			sectionsForRups.add(sectList);
		}
		
		U3FaultSystemRupSet rupSet = new U3FaultSystemRupSet(
				fltSectionDataList,
				sectSlipRate,
				sectSlipRateStdDev,
				sectArea,
				sectionsForRups,
				rupMeanMag,
				rupAveRake,
				rupArea,
				rupLength,
				modelName);
		
		return new U3FaultSystemSolution(rupSet, this.rupRatesFromMultRunsArrayList.get(solutionIndex));
		
	}
	
	
	/**
	 * This returns an array of rupture rates computed from the given MFD, where the rate in each MFD mag bin is
	 * distributed equally among the number of ruptures that fall in that mag bin
	 * @param targetMFD
	 * @return
	 */
	public double[] getRupRatesForTargetMFD(IncrementalMagFreqDist targetMFD, boolean correctForSegConstraints) {
		double[] rupRateArray = new double[this.numRuptures];
		for(int r=0;r<numRuptures;r++) {
			int index = targetMFD.getClosestXIndex(rupMeanMag[r]);
			rupRateArray[r] = targetMFD.getY(index)/meanMagHistorgram.getY(index);
//			rupRateArray[r] = targetMFD.getClosestYtoX(rupMeanMag[r])/meanMagHistorgram.getClosestYtoX(rupMeanMag[r]);
		}
		
		if(correctForSegConstraints) {
		    // adjust initial state for various segmentation type constraints
			// for slip rate segmentation constraints:
			if(slipRateSegmentationConstraintList !=null && slipRateSegmentationConstraintList.size()>0) {
				for(SlipRateSegmentationConstraint srSegConstr : slipRateSegmentationConstraintList) {		
					int sectID = srSegConstr.getSectIndex();
					double sectSlipRateTemp = 0;
					for(int rup=0; rup < numRuptures; rup++) {
						if(rupSectionMatrix[sectID][rup]==1)
							sectSlipRateTemp += rupRateArray[rup]*sectSlipInRup[sectID][rup];
					}
					double scaleFactor = srSegConstr.getSlipRateReductionFactor()*sectSlipRate[sectID]/sectSlipRateTemp;
					for(int r=0;r<rupRateArray.length;r++) {
						if(rupSectionMatrix[srSegConstr.getSectIndex()][r] == 1.0)
							rupRateArray[r] = scaleFactor*rupRateArray[r];
					}
				}
			}
			// for segmentation constraints:
			if(segmentationConstraintList !=null && segmentationConstraintList.size()>0) {
				for(SegmentationConstraint segConstr : segmentationConstraintList) {		
					int sect1 = segConstr.getSect1_Index();
					int sect2 = segConstr.getSect1_Index();
					double rateTemp = 0;
					for(int r=0; r < numRuptures; r++) {
						if(rupSectionMatrix[sect1][r]==1 && rupSectionMatrix[sect2][r]==1)
							if(segConstr.isSlipRateConstraint()) {
								double aveSlip = (sectSlipInRup[sect1][r]+sectSlipInRup[sect2][r])/2.0;
								rateTemp += rupRateArray[r]*aveSlip;
							}
							else {
								rateTemp += rupRateArray[r];
							}
					}
					double scaleFactor = segConstr.getMeanJointRate()/rateTemp;
					for(int r=0;r<rupRateArray.length;r++) {
						if(rupSectionMatrix[sect1][r]==1 && rupSectionMatrix[sect2][r]==1)
							rupRateArray[r] = scaleFactor*rupRateArray[r];
					}
				}
			}
		    // for section rate constraints:
			if(sectionRateConstraints !=null && sectionRateConstraints.size()>0) {
				for(SectionRateConstraint sectRateConstr : sectionRateConstraints) {		
					int sectID = sectRateConstr.getSectIndex();
					double sectRateTemp = 0;
					for(int r=0; r < numRuptures; r++) {
						if(rupSectionMatrix[sectID][r]==1)
							sectRateTemp += rupRateArray[r];
					}
					double scaleFactor = sectRateConstr.getMeanRate()/sectRateTemp;
	//scaleFactor=0;
					for(int r=0;r<rupRateArray.length;r++) {
						if(rupSectionMatrix[sectID][r] == 1.0)
							rupRateArray[r] = scaleFactor*rupRateArray[r];
					}
				}
			}
		}
		

		return rupRateArray;
	}
	
	/**
	 * This evenly distributes ruptures of each magnitude along the fault
	 * @param targetMFD
	 */
	public void doRatesFromMFD_Solution(IncrementalMagFreqDist targetMFD) {
		
		modelRunInfoString = "\nInversion Type = NSHMP Solution\n";

		// set rates from MFD
		rupRateSolution = getRupRatesForTargetMFD(targetMFD, false);

		// MINIMUM RATE CONSTRAINT IS INGORED

		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();
	}
	
	/**
	 * This sets the rate of each rupture as the rupMoRateMax*totMoRate/totMaxMoRate, where rupMoRateMax is the sum of 
	 * the moment rate of sections involved, and totMaxMoRate is the sum of rupMoRateMax over all ruptures.  This model
	 * is described in Visini et al., 2019, Pure Appl. Geophys, https://doi.org/10.1007/s00024-019-02114-6.
	 * @param
	 */
	public void doSUNFiSH_Solution() {
		
		modelRunInfoString = "\nInversion Type = SUNFiSH Solution\n";
		
		double totMaxMoRate = 0;
		for(int r=0;r<numRuptures;r++) {
			totMaxMoRate += rupMoRateMax[r];
		}
		double[] rupMoRateArray = new double[numRuptures];
		for(int r=0;r<numRuptures;r++) {
			rupMoRateArray[r] = rupMoRateMax[r]*totMoRate/totMaxMoRate;
		}
		rupRateSolution = new double[numRuptures];
		for(int r=0;r<numRuptures;r++) {
			rupRateSolution[r] = rupMoRateArray[r]/rupMeanMo[r];
		}

		// MINIMUM RATE CONSTRAINT IS INGORED

		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();
	}

	
	
	/**
	 * This sets the rate of each rupture from the supplied MFD, weighted by rupMoRateMax (the sum of the
	 * the moment rate of sections involved in each rupture.  This model is described in 
	 * Visini et al., 2019, Pure Appl. Geophys, https://doi.org/10.1007/s00024-019-02114-6.
	 * @param targetMFD
	 */
	public void doFRESH_Solution(IncrementalMagFreqDist targetMFD) {
		
		modelRunInfoString = "\nInversion Type = FRESH Solution\n";
		
		SummedMagFreqDist mfdWt = new SummedMagFreqDist(minMagMFD,targetMFD.size(),MAG_DELTA);

		for(int r=0;r<numRuptures;r++) {
			int index = targetMFD.getClosestXIndex(rupMeanMag[r]);
			mfdWt.add(index, rupMoRateMax[r]);
		}
		rupRateSolution = new double[numRuptures];
		for(int r=0;r<numRuptures;r++) {
			int index = targetMFD.getClosestXIndex(rupMeanMag[r]);
			rupRateSolution[r] = targetMFD.getY(index)*rupMoRateMax[r]/mfdWt.getY(index);
//System.out.println(r+"\t"+index+"\t"+rupMeanMag[r]+"\t"+rupRateSolution[r]+"\t"+rupMoRateMax[r]);
		}
//System.exit(1);


		// MINIMUM RATE CONSTRAINT IS INGORED

		// compute predicted data
		d_pred = new double[totNumRows];  // predicted data vector
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSectMFDs();
		
		setMiscRunInfo();
	}

	
	/**
	 * This enables one to add info to the modelRunInfoString (such as runtime)
	 * @param str
	 */
	public void addToModelRunInfoString(String str) {
		modelRunInfoString += str;
	}

	
	/**
	 * @param args
	 */
	public static void main(String []args) {
		
		
	}

}
