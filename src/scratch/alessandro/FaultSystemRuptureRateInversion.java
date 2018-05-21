/**
 * 
 */
package scratch.alessandro;

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
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotColorAndLineTypeSelectorControlPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.RunScript;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
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
import com.google.common.io.Files;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.simulatedAnnealing.SerialSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.ThreadedSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.completion.CompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;
import scratch.alessandro.logicTreeEnums.ScalingRelationshipEnum;
import scratch.alessandro.logicTreeEnums.SlipAlongRuptureModelEnum;

/**
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

	public final static double GAUSS_MFD_SIGMA = 0.12;
	public final static double GAUSS_MFD_TRUNCATION = 2;
	public final static double MAG_DELTA = 0.1;
	
	String modelName;
	
	// Strings to keep track of results
	String modelSetUpInfoString, modelRunInfoString;
	
	// input data
	private ArrayList<FaultSectionPrefData> fltSectionDataList;
	private ArrayList<SegRateConstraint> sectionRateConstraints; // using old class for "segments"
	private int[][] rupSectionMatrix;
	
	// section attributes
	private int numSections;
	private double[] sectArea, sectLength, sectSlipRate, sectSlipRateStdDev, sectMoRate;

	// rupture attributes
	int numRuptures;
	private String[] rupNameShort;
	private double[] rupLength, rupArea, rupMeanMag, rupMeanMo, rupMoRate, totRupRate,rupAveSlip;
	private double minRupMag, maxRupMag,minRupMagWithAleatory, maxRupMagWithAleatory;
	double[] rupRateSolution; // these are the rates from the inversion (not total rate of MFD)

	// section-rupture attributes
	private int[] numSectInRup, firstSectOfRup;
	private int minNumSectInRup;  // this sets the number of sections for the smallest ruptures
	private double[][] sectSlipInRup;
	private double[] rateOfThroughGoingRupsAtSectBoudary;

	double totMoRate;
	
	IncrementalMagFreqDist mfdConstraint; 

	private SummedMagFreqDist aveOfSectPartMFDs;
	
	// this accounts for fact that ave slip from Gauss MFD is greater than the slip of the average mag
	private double gaussMFD_slipCorr;
	
	// a-priori rate constraints
	int num_aPriori_constraints;
	int[] aPriori_rupIndex;
	double[] aPriori_rate, aPriori_wt;
	String aPrioriRupRatesFilename;

	
	private static boolean MATLAB_TEST = false;
	double[][] C_wted, C;	// inversion matrices
	double[] d, d_wted, data_wt, full_wt, d_pred;  // the data vector
	
	private double minRupRateArray[]; // the minimum rate constraint for each rupture
	private double minRupRate;
	private boolean wtedInversion, applyProbVisible;	// weight the inversion according to slip rate and segment rate uncertainties
	private double relativeSectRateWt, relative_aPrioriRupWt; //, relative_aPrioriSegRateWt;
	private double relativeMFD_constraintWt;

	private int  firstRowSectSlipRateData=-1,firstRowSectEventRateData=-1, firstRowAprioriData=-1;
	private int  lastRowSectSlipRateData=-1,lastRowSectEventRateData=-1, lastRowAprioriData=-1;
	private int firstRowGR_constraintData=-1, lastRowGR_constraintData=-1; //, firstRowParkSegConstraint=-1, lastRowParkSegConstraint=-1;
	private int totNumRows;
	
	// slip model:  CHANGE TO ENUM
	private SlipAlongRuptureModelEnum slipModelType;

	private static EvenlyDiscretizedFunc taperedSlipPDF, taperedSlipCDF;
	
	// MFDs
	SummedMagFreqDist magFreqDist, meanMagHistorgram, magHistorgram;
	ArrayList<SummedMagFreqDist> sectNucleationMFDs;
	ArrayList<SummedMagFreqDist> sectParticipationMFDs;
	
	
	
	private double[] finalSectEventRate, finalPaleoVisibleSectEventRate, finalSectSlipRate;
	
	// the following is the total moment-rate reduction, including that which goes to the  
	// background, afterslip, events smaller than the min mag here, and aftershocks and foreshocks.
	private double moRateReduction;  
	
	private ScalingRelationshipEnum magAreaRel;
	
	// NNLS inversion solver - static to save time and memory
	private static NNLSWrapper nnls = new NNLSWrapper();

	
	
	
	
	

	/**
	 * This writes the segPartMFDs to a file and/or makes a plot of the
	 * results. 
	 * @param dirName - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 */
	public void writeAndOrPlotSegPartMFDs(String dirName, boolean popUpWindow) {
		
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
			
			make2D_plot(xyzDataSectPart, "Incr. Participation", "Section", "Magnitude", "Rate", fileName, popUpWindow);
			make2D_plot(xyzDataSectPartCum, "Cum. Participation", "Section", "Magnitude", "CumRate", fileNameCum, popUpWindow);		

		}catch(Exception e) {
			e.printStackTrace();
		}
	}

	
	
	
	
	
	/**
	 * This writes and/or plots the rates and slip-rate contribution of each rupture that has a 
	 * final rate above the minimum.  
	 * @param dirName - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 */
	public void writeAndOrPlotNonZeroRateRups(String dirName, boolean popUpWindow) {
		
		// get number of ruptures above minRupRate
		int numAboveMinRate = 0;
		for(int rup=0; rup<numRuptures;rup++)
			if(rupRateSolution[rup]>minRupRate)
				numAboveMinRate += 1;

		EvenlyDiscrXYZ_DataSet xyzDataRupSlipRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);
		EvenlyDiscrXYZ_DataSet xyzDataRupRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);
		

		int index = 0;
		try{
			FileWriter fw=null,fw2=null;
			String fileName_sr=null, fileName_rr=null;
			if(dirName != null) {
				fileName_sr = dirName+"/rupSlipRateData";
				fileName_rr = dirName+"/rupRateData";
				fw = new FileWriter(fileName_sr+".txt");
				fw2 = new FileWriter(fileName_rr+".txt");
				fw.write("index\tseg_index\tslip_rate\n");
				fw2.write("index\tseg_index\trate\n");				
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
			make2D_plot(xyzDataRupSlipRate, "Slip Rate from Ruptures", "Section", "Rupture", "Log10 Slip Rate", fileName_sr, popUpWindow);
			make2D_plot(xyzDataRupRate, "Rate of Ruptures", "Section", "Rupture", "Log10 Rate", fileName_rr, popUpWindow);


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
	 *
	 * @param fltSectionDataList - any name the user wants to give the inversion model
	 * @param fltSectionDataList
	 * @param sectionRateConstraints
	 * @param rupSectionMatrix
	 * @param slipModelType - TAPERED_SLIP_MODEL, UNIFORM_SLIP_MODEL, or WG02_SLIP_MODEL
	 * @param magAreaRel
	 * @param relativeSectRateWt - weight on segment rate equations (relative to slip rate)
	 * @param relative_aPrioriRupWt - weight on a-priori rates (relative to slip rate)
	 * @param aPrioriRupRatesFilename - filename for aPriori rupture rates
	 * @param wtedInversion - apply data uncertainties?
	 * @param minRupRate - constrain all rupture rates to be greater than this value
	 * @param applyProbVisible - account for likelihood that Paleoseismology will see the rupture
	 * @param moRateReduction - fraction reduction from smaller events (and not aseismicity or coupling coefficient, which are set in fltSectionDataList)
	 * @param relativeMFD_constraintWt - weight for MFD constraint
	 */
	public FaultSystemRuptureRateInversion( String modelName,
			ArrayList<FaultSectionPrefData> fltSectionDataList, 
			ArrayList<SegRateConstraint> sectionRateConstraints,
			int[][] rupSectionMatrix, 
			SlipAlongRuptureModelEnum slipModelType, 
			ScalingRelationshipEnum scalingRel, 
			double relativeSectRateWt, 
			double relative_aPrioriRupWt, 
			String aPrioriRupRatesFilename,
			boolean wtedInversion, 
			double minRupRate,
			boolean applyProbVisible, 
			double moRateReduction, 
			IncrementalMagFreqDist mfdConstraint,
			double relativeMFD_constraintWt) {
		
		this.modelName = modelName;
		this.fltSectionDataList = fltSectionDataList;
		this.sectionRateConstraints = sectionRateConstraints;
		this.rupSectionMatrix = rupSectionMatrix;
		this.slipModelType = slipModelType;
		this.magAreaRel = scalingRel;
		this.relativeSectRateWt = relativeSectRateWt;
		this.relative_aPrioriRupWt = relative_aPrioriRupWt;
		this.aPrioriRupRatesFilename = aPrioriRupRatesFilename;
		this.wtedInversion = wtedInversion;
		this.minRupRate = minRupRate;
		this.applyProbVisible = applyProbVisible;
		this.moRateReduction = moRateReduction;
		this.mfdConstraint = mfdConstraint;
		this.relativeMFD_constraintWt = relativeMFD_constraintWt;
		
		if(modelName == null) {
			modelName = this.getClass().toString();
		}

		// set info string
		modelSetUpInfoString = modelName+" with:\n\n";
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
		modelSetUpInfoString += "\tGAUSS_MFD_SIGMA = "+GAUSS_MFD_SIGMA+"\n";
		modelSetUpInfoString += "\tGAUSS_MFD_TRUNCATION = "+GAUSS_MFD_TRUNCATION+"\n";

		
		// compute slip correction for Gaussian MFD
		GaussianMagFreqDist gDist1 = new GaussianMagFreqDist(5.0,9.0,41,7,this.GAUSS_MFD_SIGMA,1.0,this.GAUSS_MFD_TRUNCATION,2);
		GaussianMagFreqDist gDist2 = new GaussianMagFreqDist(5.0,9.0,41,7,0.01,1.0,0.01,2);
		gDist1.scaleToCumRate(0, 1.0);
		gDist2.scaleToCumRate(0, 1.0);
		gaussMFD_slipCorr = gDist1.getTotalMomentRate()/gDist2.getTotalMomentRate();
//		System.out.println("gaussMFD_slipCorr="+(float)gaussMFD_slipCorr+"\n");
		
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
		int numSectRateConstraints = sectionRateConstraints.size();
		
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
		
		
		// add number of MFD constaints
		int numMFD_constraints=0;
		if(relativeMFD_constraintWt>0) {
			numMFD_constraints = mfdConstraint.size();
			firstRowGR_constraintData = totNumRows;
			totNumRows += numMFD_constraints;
			lastRowGR_constraintData = totNumRows-1;
		}
				
		String tempString = "firstRowSegEventRateData = "+firstRowSectEventRateData+
				"; firstRowAprioriData = "+firstRowAprioriData+
				"; firstRowGR_constraintData = "+firstRowGR_constraintData+
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
		
		// now fill in the segment event rate constraints if requested
		if(relativeSectRateWt > 0.0) {
			SegRateConstraint constraint;
			for(int i = 0; i < numSectRateConstraints; i ++) {
				constraint = sectionRateConstraints.get(i);
				int seg = constraint.getSegIndex();
				int row = i+firstRowSectEventRateData;
				d[row] = constraint.getMean(); // this is the average sub-section rate
				if(wtedInversion)
					data_wt[row] = 1/constraint.getStdDevOfMean();
				for(int col=0; col<numRuptures; col++)
					if(applyProbVisible)
						C[row][col] = rupSectionMatrix[seg][col]*getProbVisible(rupMeanMag[col]);
					else
						C[row][col] = rupSectionMatrix[seg][col];
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
		
		
		// now fill in the MFD constraint if needed
		if(relativeMFD_constraintWt > 0.0) {
			for(int i=0; i < numMFD_constraints; i++) {
				int row = i+firstRowGR_constraintData;
				double mag = mfdConstraint.getX(i);
				d[row] = mfdConstraint.getY(mag);
				for(int col=0; col<numRuptures; col++)
					if(mfdConstraint.getClosestXIndex(rupMeanMag[col]) == i)
						C[row][col]=1.0;
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
		// segment event rate wts
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
				C_wted[row][col]=full_wt[row];
			}
		}

		// MFD constraint wts
		if(relativeMFD_constraintWt > 0.0) {
			for(int i=0; i < numMFD_constraints; i++) {
				int row = i+firstRowGR_constraintData;
				full_wt[row] = relativeMFD_constraintWt;
				d_wted[row] *= full_wt[row];
				for(int rup=0; rup < numRuptures; rup++) C_wted[row][rup] *= full_wt[row];
			}
		}
		
	}
	
	/**
	 * This deos the inversion using non-negative least squares
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

//		// set initial state from MFD constraint
////		double[] initialState = aPriori_rate;
//		// or set initial state to zero
//		double[] initialState = new double[numRuptures];
//		rupRateSolution = getSimulatedAnnealingSolution(C_wted, d_wted, initialState);

		// CORRECT FINAL RATES IF MINIMUM RATE CONSTRAINT APPLIED
		if(minRupRate >0.0)
			for(int rup=0; rup<numRuptures;rup++) rupRateSolution[rup] += minRupRateArray[rup];

		// compute predicted data
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSegMFDs();
		
		setMiscRunInfo();
		
	}
	
		
	/**
	 * This does the inversion using simulated annealing
	 */
	public void doInversionSA(long numIterations, boolean initStateFromAprioriRupRates, long randomSeed) {
		
		
		modelRunInfoString = "\nInversion Type = Simulated Annearling with\n\n\tnumIterations = "+
				numIterations+"\n\tinitStateFromAprioriRupRates = "+initStateFromAprioriRupRates+
				"\n\trandomSeed = "+randomSeed+"\n";

		// set the intial state
		double[] initialState;
		if(initStateFromAprioriRupRates) {
			// check that there is a rate for each rupture
			if(aPriori_rate.length != numRuptures)
				throw new RuntimeException("Can't use a priori rupture rates as initialState because sizes are different");
			for(int r=0;r<numRuptures;r++) {
				if(this.aPriori_rupIndex[r] != r)
					throw new RuntimeException("Can't use a priori rupture rates as initialState because there is an index problem");
			}
			initialState = aPriori_rate;
		}
		else { // inital state is all zeros
			initialState = new double[numRuptures];	
		}
		
		// SOLVE THE INVERSE PROBLEM
		rupRateSolution = getSimulatedAnnealingSolution(C_wted, d_wted, initialState, numIterations, randomSeed);
		
		// Not sure how the following works (not faster but lower prediction error); and doesn't appear to take a seed
//		rupRateSolution = getSimulatedAnnealingThreadedSolution(C_wted, d_wted, initialState, numIterations, randomSeed);

		// CORRECT FINAL RATES IF MINIMUM RATE CONSTRAINT APPLIED
		if(minRupRate >0.0)
			for(int rup=0; rup<numRuptures;rup++) rupRateSolution[rup] += minRupRateArray[rup];

		// compute predicted data
		for(int row=0;row<totNumRows; row++)
			for(int col=0; col <numRuptures; col++)
				d_pred[row] += rupRateSolution[col]*C[row][col];
				
		// Compute final segment slip rates and event rates
		computeFinalStuff();
		
		computeSegMFDs();
		
		setMiscRunInfo();
		
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
		
		taperedSlipCDF = new EvenlyDiscretizedFunc(0, 5001, 0.0002);
		taperedSlipPDF = new EvenlyDiscretizedFunc(0, 5001, 0.0002);
		double x,y, sum=0;
		int num = taperedSlipPDF.size();
		for(int i=0; i<num;i++) {
			x = taperedSlipPDF.getX(i);
			// y = Math.sqrt(1-(x-0.5)*(x-0.5)/0.25);
			y = Math.pow(Math.sin(x*Math.PI), 0.5);
			taperedSlipPDF.set(i,y);
			sum += y;
		}

		// now make final PDF & CDF
		y=0;
		for(int i=0; i<num;i++) {
				y += taperedSlipPDF.getY(i);
				taperedSlipCDF.set(i,y/sum);
				taperedSlipPDF.set(i,taperedSlipPDF.getY(i)/sum);
//				System.out.println(taperedSlipCDF.getX(i)+"\t"+taperedSlipPDF.getY(i)+"\t"+taperedSlipCDF.getY(i));
		}
	}
	
	
	/**
	 * This computes the total event rate prediction error (ignoring wt given to this equation set),
	 * meaning it will give a result even if the equation set wt was zero.
	 * @return
	 */
	private double getTotEventRatePredErr() {
		SegRateConstraint constraint;
		double totPredErr=0;
		for(int i = 0; i < sectionRateConstraints.size(); i ++) {
			constraint = sectionRateConstraints.get(i);
			int seg = constraint.getSegIndex();
			double normResid = (finalPaleoVisibleSectEventRate[seg]-constraint.getMean())/constraint.getStdDevOfMean();
			totPredErr += normResid*normResid;
		}
		return totPredErr;
	}
	
	

	
	private static double[] getSimulatedAnnealingSolution(double[][] C, double[] d, double[] initialState, long numIterations,  long randomSeed) {
		SparseDoubleMatrix2D matrixC = new SparseDoubleMatrix2D(C); //
		SerialSimulatedAnnealing simulatedAnnealing =new SerialSimulatedAnnealing(matrixC, d, initialState);
		simulatedAnnealing.setRandom(new Random(randomSeed));
		simulatedAnnealing.iterate(numIterations);
		return simulatedAnnealing.getBestSolution();
	}



	private static double[] getSimulatedAnnealingThreadedSolution(double[][] C, double[] d, double[] initialState, long numIterations,  long randomSeed) {
		SparseDoubleMatrix2D matrixC = new SparseDoubleMatrix2D(C); //
		//this is the "sub completion criteria" - the amount of time (or iterations) between synchronization
		CompletionCriteria subCompetionCriteria = TimeCompletionCriteria.getInSeconds(1); // 1 second;
		// this will use all available processors
		int numThreads = Runtime.getRuntime().availableProcessors();

		ThreadedSimulatedAnnealing simulatedAnnealing = new ThreadedSimulatedAnnealing(
				matrixC, d, initialState, numThreads, subCompetionCriteria);
//		simulatedAnnealing.setRandom(new Random(randomSeed));
		simulatedAnnealing.iterate(numIterations);
		return simulatedAnnealing.getBestSolution();
	}


	
	
	/**
	 * This gets the non-negative least squares solution for the matrix C
	 * and data vector d.
	 * @param C
	 * @param d
	 * @return
	 */
	private static double[] getNNLS_solution(double[][] C, double[] d) {

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
		for(int seg=0; seg < numSections; seg++) {
			finalSectSlipRate[seg] = 0;
			finalSectEventRate[seg] = 0;
			for(int rup=0; rup < numRuptures; rup++) 
				if(rupSectionMatrix[seg][rup]==1) {
					finalSectSlipRate[seg] += rupRateSolution[rup]*sectSlipInRup[seg][rup];
					finalSectEventRate[seg]+=rupRateSolution[rup];
					if(applyProbVisible)
						finalPaleoVisibleSectEventRate[seg]+=rupRateSolution[rup]*getProbVisible(rupMeanMag[rup]);
					else
						finalPaleoVisibleSectEventRate[seg]+=rupRateSolution[rup];
				}
		}
		
		// Compute the total Mag Freq Dist
		int num = (int)Math.round((maxRupMagWithAleatory-minRupMagWithAleatory)/MAG_DELTA + 1);
		magFreqDist = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
		for(int rup=0; rup<numRuptures;rup++) {
			// add MFD of this rupture
			GaussianMagFreqDist gDist = new GaussianMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA,rupMeanMag[rup],GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
			gDist.scaleToCumRate(0, rupRateSolution[rup]);
			magFreqDist.addIncrementalMagFreqDist(gDist);
		}
		magFreqDist.setInfo("Incremental Mag Freq Dist");
		
		
		// check moment rate
		double origMoRate = totMoRate*(1-moRateReduction);
		double ratio = origMoRate/magFreqDist.getTotalMomentRate();
		String tempString = "Moment Rates: from Fault Sections (possible reduced) = "+ (float)origMoRate+
				";  from total MFD = "+(float)magFreqDist.getTotalMomentRate()+
				",  ratio = "+(float)ratio;
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
	private void computeSegMFDs() {
		int num = (int)Math.round((maxRupMagWithAleatory-minRupMagWithAleatory)/MAG_DELTA + 1);
		sectNucleationMFDs = new ArrayList<SummedMagFreqDist>();
		sectParticipationMFDs = new ArrayList<SummedMagFreqDist>();
		SummedMagFreqDist sumOfSegPartMFDs = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
		aveOfSectPartMFDs = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
		
		SummedMagFreqDist segPartMFD, segNuclMFD;
//		double mag, rate;
		for(int seg=0; seg < numSections; seg++) {
			segPartMFD = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
			segNuclMFD = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
			for(int rup=0; rup < numRuptures; rup++) {
				if(this.rupSectionMatrix[seg][rup] == 1) {
					if(rupRateSolution[rup] > 0) {
						GaussianMagFreqDist mfd = new GaussianMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA,rupMeanMag[rup],GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2); // dist w/ unit moment rate
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
			SegRateConstraint constraint;
			for(int i = 0; i < sectionRateConstraints.size(); i ++) {
				int row = firstRowSectEventRateData+i;
				constraint = sectionRateConstraints.get(i);
				int seg = constraint.getSegIndex();
//				this checks that finalSegEventRate[seg] and d_pred[row} are the same (RECONSIDER PROB VISIBLE IF I USE THIS AGAIN)
//				System.out.println(seg+"\t"+(float)finalSegEventRate[seg]+"\t"+(float)d_pred[row]+"\t"+(float)(finalSegEventRate[seg]/d_pred[row]));
				if(applyProbVisible) {
					aveRatio += finalPaleoVisibleSectEventRate[seg]/constraint.getMean();
//					System.out.println(seg+"\t"+(float)finalPaleoReducedSegEventRate[seg]+"\t"+(float)constraint.getMean()+"\t"+(float)aveRatio);				
				}
				else {
					aveRatio += finalSectEventRate[seg]/constraint.getMean();
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
		
		// First without equation weights
		double totPredErr=0, slipRateErr=0, eventRateErr=0, aPrioriErr=0, grErr=0;
		
		for(int row=firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row++)
			slipRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relativeSectRateWt >0)
			for(int row=firstRowSectEventRateData; row <= lastRowSectEventRateData; row++)
				eventRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relative_aPrioriRupWt > 0)
			for(int row=firstRowAprioriData; row <= lastRowAprioriData; row++)
				aPrioriErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relativeMFD_constraintWt>0)
			for(int row=firstRowGR_constraintData; row <= lastRowGR_constraintData; row++){
				grErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
			}

		totPredErr = slipRateErr+eventRateErr+aPrioriErr+grErr;
		
		// get alt versions in case equation wts zero
		double eventRateErrAlt = getTotEventRatePredErr();
		
		double aveSRnormResid = Math.sqrt(slipRateErr/numSections);
		double aveERnormResid = Math.sqrt(eventRateErrAlt/sectionRateConstraints.size());
		
		String resultsString = "\nTotal Prediction Error =\t"+(float)totPredErr+"\n\t"+
				               "Slip Rate Err =\t\t"+(float)slipRateErr+"\trel. wt = 1.0"+
				               "\tnorm resid rms = "+(float)aveSRnormResid+"\n\t"+
				               "Event Rate Err =\t"+(float)eventRateErrAlt+"\trel. wt = "+relativeSectRateWt+
				               "\tnorm resid rms = "+(float)aveERnormResid+"\n\t"+
				               "A Priori Err =\t\t"+(float)aPrioriErr+"\trel. wt = "+relative_aPrioriRupWt+"\n\t"+
        					   "GR Err =\t"+(float)grErr+"\trel. wt = "+relativeMFD_constraintWt+"\n";

		if(D) System.out.println(resultsString);
		modelRunInfoString += "\n"+resultsString+"\n";
				
		// Now with equation weights
		totPredErr=0; slipRateErr=0; eventRateErr=0; aPrioriErr=0; grErr=0;
		for(int row=firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row++)
			slipRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relativeSectRateWt >0)
			for(int row=firstRowSectEventRateData; row <= lastRowSectEventRateData; row++)
				eventRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relative_aPrioriRupWt > 0)
			for(int row=firstRowAprioriData; row <= lastRowAprioriData; row++)
				aPrioriErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relativeMFD_constraintWt>0)
			for(int row=firstRowGR_constraintData; row <= lastRowGR_constraintData; row++)
				grErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];

		totPredErr = slipRateErr+eventRateErr+aPrioriErr+grErr;

		resultsString = "\nTotal Pred Err w/ Eq Wts =\t"+(float)totPredErr+"\n\t"+
				"Slip Rate Err =\t\t"+(float)slipRateErr+"\trel. wt = 1.0\n\t"+
				"Event Rate Err =\t"+(float)eventRateErr+"\trel. wt = "+relativeSectRateWt+"\n\t"+
				"A Priori Err =\t\t"+(float)aPrioriErr+"\trel. wt = "+relative_aPrioriRupWt+"\n\t"+
				"GR Err =\t"+(float)grErr+"\trel. wt = "+relativeMFD_constraintWt+"\n";
		
		if(D) System.out.println(resultsString);
		modelRunInfoString += "\n"+resultsString+"\n";

			
	}
	
	/**
	 * This writes and/or plot various data fits.
	 * @param dirName - set as null if you don't want to save results
	 * @param popupWindow - set as true if you want plot windows to pop up
	 */
	public void writeAndOrPlotDataFits(String dirName, boolean popupWindow) {
		
		// SLIP RATES:
		double min = 0, max = numSections-1;
		EvenlyDiscretizedFunc origSlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
//		EvenlyDiscretizedFunc origUpper95_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
//		EvenlyDiscretizedFunc origLower95_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
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
//			origUpper95_SlipRateFunc.set(s,(sectSlipRate[s]+1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
//			origLower95_SlipRateFunc.set(s,(sectSlipRate[s]-1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
			finalSlipRateFunc.set(s,finalSectSlipRate[s]);
		}
		ArrayList<XY_DataSet> sr_funcs = new ArrayList<XY_DataSet>();
		origSlipRateFunc.setName("Orig Slip Rates");
//		origUpper95_SlipRateFunc.setName("Orig Upper 95% Confidence for Slip Rates");
//		origLower95_SlipRateFunc.setName("Orig Lower 95% Confidence for Slip Rates");
		finalSlipRateFunc.setName("Final Slip Rates");
		sr_funcs.add(finalSlipRateFunc);
		sr_funcs.add(origSlipRateFunc);
		sr_funcs.addAll(slipRate95confFuncsList);
		
		ArrayList<PlotCurveCharacterstics> sr_plotChars = new ArrayList<PlotCurveCharacterstics>();
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLUE));
		for(int i=0;i<slipRate95confFuncsList.size();i++) {
			sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 4f, Color.BLUE));
		}
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/slipRates";

		String plotName = "";
		String xAxisLabel = "Section";
		String yAxisLabel = "Slip Rate (m/sec)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		writeAndOrPlotFuncs(sr_funcs, sr_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		
		
		
		// EVENT RATES:
		ArrayList<XY_DataSet> er_funcs = new ArrayList<XY_DataSet>();
		// now fill in final event rates
		EvenlyDiscretizedFunc finalEventRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc finalPaleoVisibleEventRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		for(int s=0;s < numSections; s++) {
			finalEventRateFunc.set(s,finalSectEventRate[s]);
			finalPaleoVisibleEventRateFunc.set(s, finalPaleoVisibleSectEventRate[s]);
		}
		finalPaleoVisibleEventRateFunc.setName("Final Paleoseismically Visible Event Rates");
		finalEventRateFunc.setName("Final Event Rates (dashed)");
		er_funcs.add(finalPaleoVisibleEventRateFunc);
		er_funcs.add(finalEventRateFunc);
		int num = sectionRateConstraints.size();
		DefaultXY_DataSet func;
		ArbitrarilyDiscretizedFunc meanER_Func = new ArbitrarilyDiscretizedFunc();
		er_funcs.add(meanER_Func);
		ArrayList<AbstractXY_DataSet> obs_er_funcs = new ArrayList<AbstractXY_DataSet>();
		obs_er_funcs.add(meanER_Func);
		SegRateConstraint constraint;
		for(int c=0;c<num;c++) {
			func = new DefaultXY_DataSet();
			constraint = sectionRateConstraints.get(c);
			int s = constraint.getSegIndex();
			meanER_Func.set((double)s, constraint.getMean());
			func.set((double)s, constraint.getLower95Conf());
			func.set((double)s, constraint.getUpper95Conf());
			func.setName(constraint.getFaultName());
			obs_er_funcs.add(func);
			er_funcs.add(func);
		}
		ArrayList<PlotCurveCharacterstics> er_plotChars = new ArrayList<PlotCurveCharacterstics>();
		er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		er_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		er_plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
		for(int c=0;c<num;c++)
			er_plotChars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 4f, Color.RED));

		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/eventRates";
		plotName = "";
		xAxisLabel = "Section";
		yAxisLabel = "Event Rate (per yr)";
		xAxisRange = null;
		yAxisRange = null;
		logX = false;
		logY = true;

		writeAndOrPlotFuncs(er_funcs, er_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		
				
		// RUPTURE RATES
		// plot the final rupture rates
		max = numRuptures-1;
		EvenlyDiscretizedFunc rupRateFunc = new EvenlyDiscretizedFunc(min, max, numRuptures);
		for(int rup=0; rup<numRuptures;rup++) {
			rupRateFunc.set(rup,rupRateSolution[rup]);
		}
		rupRateFunc.setName("Rupture Rates");
		ArrayList<XY_DataSet> rup_funcs = new ArrayList<XY_DataSet>();
		rup_funcs.add(rupRateFunc);
		ArrayList<PlotCurveCharacterstics> rup_plotChars = new ArrayList<PlotCurveCharacterstics>();
		rup_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/ruptureRates";
		plotName = "Rupture Rates";
		xAxisLabel = "Rupture Index";
		yAxisLabel = "Rate (per yr)";
		xAxisRange = new Range(0, numRuptures);
		yAxisRange = new Range(minRupRate/2.0, rupRateFunc.getMaxY());
		logX = false;
		logY = true;

		writeAndOrPlotFuncs(rup_funcs, rup_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		

		// PLOT MFDs
		ArrayList<XY_DataSet> mfd_funcs = new ArrayList<XY_DataSet>();
		mfd_funcs.add(magFreqDist);
		EvenlyDiscretizedFunc cumMagFreqDist = magFreqDist.getCumRateDistWithOffset();
		cumMagFreqDist.setInfo("Cumulative Mag Freq Dist");
		mfd_funcs.add(cumMagFreqDist);
		// add average seg participation MFD
		mfd_funcs.add(aveOfSectPartMFDs);
		EvenlyDiscretizedFunc cumAveOfSegPartMFDs = aveOfSectPartMFDs.getCumRateDistWithOffset();
		cumAveOfSegPartMFDs.setInfo("cumulative "+aveOfSectPartMFDs.getInfo());
		
		GutenbergRichterMagFreqDist grFit = getGR_DistFit();
		EvenlyDiscretizedFunc cumGR_fit = grFit.getCumRateDistWithOffset();
		cumGR_fit.setName("Cumulative GR Fit");
		mfd_funcs.add(grFit);
		mfd_funcs.add(cumGR_fit);
		
		ArrayList<PlotCurveCharacterstics> mfd_plotChars = new ArrayList<PlotCurveCharacterstics>();
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.MAGENTA));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		mfd_plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK));
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/MFDs";
		plotName = "Magnitude Frequency Distributions";
		xAxisLabel = "Magnitude";
		yAxisLabel = "Rate (per yr)";
		xAxisRange = new Range(minRupMagWithAleatory, maxRupMagWithAleatory);
		yAxisRange = new Range(1e-7, 0.1);
		logX = false;
		logY = true;

		writeAndOrPlotFuncs(mfd_funcs, mfd_plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);




		// PLOT RATE AT WHICH RUPTURES PASS THROUGH EACH SECTION BOUNDARY
		min = 0; max = numSections;
		EvenlyDiscretizedFunc rateOfThroughgoingRupsFunc = new EvenlyDiscretizedFunc(min, max, numSections+1);
		for(int s=0; s<numSections+1;s++) {
			rateOfThroughgoingRupsFunc.set(s,rateOfThroughGoingRupsAtSectBoudary[s]);
		}
		ArrayList<XY_DataSet> sect_funcs = new ArrayList<XY_DataSet>();
		rateOfThroughgoingRupsFunc.setName("Rate of throughgoing ruptures at each section boundary");
		sect_funcs.add(finalSlipRateFunc);
		sect_funcs.add(finalPaleoVisibleEventRateFunc);
		// add paleoseismic obs
		for(int i=0; i<obs_er_funcs.size();i++) sect_funcs.add(obs_er_funcs.get(i));
		sect_funcs.add(rateOfThroughgoingRupsFunc);
		
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
		for(int c=0;c<num;c++)
			plotChars2.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 1f, PlotSymbol.CROSS, 4f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/sectionBoundaryRates";
		plotName =  "Section Boundary Rates";
		xAxisLabel = "Section Boundary";
		yAxisLabel = "Rate (per yr)";
		xAxisRange = null;
		yAxisRange = null;
		logX = false;
		logY = false;

		writeAndOrPlotFuncs(sect_funcs, plotChars2, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);


	}

	
	

	/**
	 * This makes histograms of the number of ruptures in each magnitude interval and the number
	 * of sections in each
	 * @param dirName - set as null if you don't want to save results
	 * @param popupWindow - set as true if you want plot windows to pop up
	 */
	public void writeAndOrPlotMagHistograms(String dirName, boolean popupWindow) {
		
		ArrayList<XY_DataSet> funcs1 = new ArrayList<XY_DataSet>();
		funcs1.add(meanMagHistorgram);
		funcs1.add(magHistorgram);

		
		// make a numSegInRupHistogram
		SummedMagFreqDist numSegInRupHistogram = new SummedMagFreqDist(1.0,numSections,1.0);
		for(int r=0;r<numSectInRup.length;r++) numSegInRupHistogram.add((double)numSectInRup[r], 1.0);
		numSegInRupHistogram.setName("Num Segments In Rupture Histogram");
		ArrayList<XY_DataSet> funcs2 = new ArrayList<XY_DataSet>();
		funcs2.add(numSegInRupHistogram);		
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/magHistograms";
		String plotName ="Mag Histograms";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Num Ruptures";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		writeAndOrPlotFuncs(funcs1, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/numSectInRuptHistogram";
		plotName = "Num Sect In Rup Histograms";
		xAxisLabel = "Num Section";
		yAxisLabel = "Num Ruptures";
		xAxisRange = null;
		yAxisRange = null;
		logX = false;
		logY = false;

		writeAndOrPlotFuncs(funcs2, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

	}
	
	
	/**
	 * The general x-y plotting method
	 * @param funcs
	 * @param plotChars
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	protected void writeAndOrPlotFuncs(
			ArrayList<XY_DataSet> funcs, 
			ArrayList<PlotCurveCharacterstics> plotChars, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		
		if(popupWindow) {
			GraphWindow graph = new GraphWindow(funcs, plotName);
			if(xAxisRange != null)
				graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
			if(yAxisRange != null)
				graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
			graph.setXLog(logX);
			graph.setYLog(logY);
			graph.setPlotChars(plotChars);
			graph.setX_AxisLabel(xAxisLabel);
			graph.setY_AxisLabel(yAxisLabel);
			graph.setTickLabelFontSize(18);
			graph.setAxisLabelFontSize(20);
			
			if(fileNamePrefix != null) {
				try {
					graph.saveAsPDF(fileNamePrefix+".pdf");
					graph.saveAsPNG(fileNamePrefix+".png");
					graph.saveAsTXT(fileNamePrefix+".txt");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}		
		}
		else if (fileNamePrefix != null){	// no pop-up window; just save plot
			PlotSpec spec = new PlotSpec(funcs, plotChars, plotName, xAxisLabel, yAxisLabel);
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(24);
			gp.setAxisLabelFontSize(26);
			gp.drawGraphPanel(spec, logX, logY);
			gp.getChartPanel().setSize(1000, 800);
			
			try {
				gp.saveAsPNG(fileNamePrefix+".png");
				gp.saveAsPDF(fileNamePrefix+".pdf");
				gp.saveAsTXT(fileNamePrefix+".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	
	
	private GutenbergRichterMagFreqDist getGR_DistFit() {
		int num = (int)Math.round((maxRupMag-minRupMag)/0.1 + 1);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(minRupMag,num,0.1);
		double moRate = totMoRate*(1-moRateReduction);
//		double altMoRate = magFreqDist.getTotalMomentRate();
		gr.setAllButTotCumRate(minRupMag, maxRupMag, moRate, 1.0);
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
	public void writeApriorRupRatesFromMFD_Constrint(String fileName) {
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

		// set section data values
		numSections = fltSectionDataList.size();
		if(D) System.out.println("numSections="+numSections);
		modelSetUpInfoString += "\nnumSections = "+numSections+"\n";
		
		sectSlipRate = new double[numSections];
		sectSlipRateStdDev = new double[numSections];
		sectLength = new double[numSections];
		sectArea = new double[numSections];
		sectMoRate = new double[numSections];
		totMoRate=0;
		for(int s=0;s<numSections;s++) {
			FaultSectionPrefData fltSectData = fltSectionDataList.get(s);
			// Note that moRateReduction is not applied in the following, but later when setting up matrix
			sectSlipRate[s] = fltSectData.getReducedAveSlipRate()*1e-3;	// convert to meters/yr
			sectSlipRateStdDev[s] = fltSectData.getReducedSlipRateStdDev()*1e-3;	// convert to meters/yr
			sectLength[s] = fltSectData.getTraceLength()*1e3; // convert to meters
			sectArea[s] = sectLength[s] * fltSectData.getReducedDownDipWidth()*1e3; // convert latter to meters
			sectMoRate[s] = fltSectData.calcMomentRate(true);
			totMoRate += sectMoRate[s];
		}

		// set rupture attributes
		numRuptures = rupSectionMatrix[0].length;
		if(D) System.out.println("numRuptures="+numRuptures+"\n");
		modelSetUpInfoString += "\nnumRuptures = "+numRuptures+"\n";


		rupNameShort = new String[numRuptures];
		rupLength = new double[numRuptures];
		rupArea = new double[numRuptures];
		numSectInRup = new int[numRuptures];
		firstSectOfRup = new int[numRuptures];
		minNumSectInRup = Integer.MAX_VALUE;

		for(int r=0;r<numRuptures;r++) {
			boolean foundFirstSect = false;
			for(int s=0; s<numSections; s++) {
				if(rupSectionMatrix[s][r] == 1) {
					rupLength[r] += sectLength[s];
					rupArea[r] += sectArea[s];
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
		}

		if(D) System.out.println("minNumSectInRup="+minNumSectInRup);
		modelSetUpInfoString += "\nminNumSectInRup = "+minNumSectInRup+"\n";

		
		// compute rupture mean mags etc
		rupMeanMag = new double[numRuptures];
		rupAveSlip = new double[numRuptures];
		rupMeanMo = new double[numRuptures];
		minRupMag = Double.MAX_VALUE;
		maxRupMag = 0;
		
		for(int r=0; r <numRuptures; r++) {
			double mag = magAreaRel.getMag(rupArea[r], rupLength[r], Double.NaN);	// orig down dip with set to NaN because scaling relationships used here don't depende on this
			//round this to nearest 10th unit
			rupMeanMag[r] = ((double)Math.round(10*mag))/10.0;
			rupMeanMo[r] = MagUtils.magToMoment(rupMeanMag[r])*gaussMFD_slipCorr;   // increased if magSigma >0
//			rupAveSlip[r] = rupMeanMo[r]/(rupArea[r]*FaultMomentCalc.SHEAR_MODULUS);  // inlcudes aveSlipCor in rupMeanMo
			rupAveSlip[r] = magAreaRel.getAveSlip(rupArea[r], rupLength[r], Double.NaN)*gaussMFD_slipCorr;

			if(minRupMag>rupMeanMag[r])
				minRupMag=rupMeanMag[r];
			if(maxRupMag<rupMeanMag[r])
				maxRupMag=rupMeanMag[r];
		}

		double tempMag = minRupMag-GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		minRupMagWithAleatory = ((double)Math.round(10*tempMag))/10.0;
		tempMag = maxRupMag+GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		maxRupMagWithAleatory = ((double)Math.round(10*tempMag))/10.0;
		String tempString = "maxRupMag = "+maxRupMag+"; tminRupMag = "+minRupMag+
				"; minRupMagWithAleatory = "+minRupMagWithAleatory+
				"; maxRupMagWithAleatory = "+maxRupMagWithAleatory;
		
		if(D) System.out.println(tempString);
		modelSetUpInfoString += "\n"+tempString+"\n";

		
		// compute meanMagHistorgram
		int num = (int)Math.round((maxRupMag-minRupMag)/MAG_DELTA + 1);
		meanMagHistorgram = new SummedMagFreqDist(minRupMag,num,MAG_DELTA);
		for(int rup=0; rup<numRuptures;rup++) {
			meanMagHistorgram.add(rupMeanMag[rup], 1.0);
		}
		meanMagHistorgram.setInfo("Mean Mag Histogram");
		
		// compute mag historgram with aleatory variability
		num = (int)Math.round((maxRupMagWithAleatory-minRupMagWithAleatory)/MAG_DELTA + 1);
		magHistorgram = new SummedMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA);
		for(int rup=0; rup<numRuptures;rup++) {
			// add to mag historgrams
			GaussianMagFreqDist gDist = new GaussianMagFreqDist(minRupMagWithAleatory,num,MAG_DELTA,rupMeanMag[rup],GAUSS_MFD_SIGMA,1.0,GAUSS_MFD_TRUNCATION,2);
			gDist.scaleToCumRate(0, 1.0); // this makes it a PDF
			magHistorgram.addIncrementalMagFreqDist(gDist);
		}
		magHistorgram.setInfo("Mag Historgram (including aleatory variability)");

		
		
	}
	
	/**
	 * This is the general 2D plotting method.
	 * @param xyzData
	 * @param title
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param zAxisLabel
	 * @param fileNamePrefix - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it

	 */
	protected static void make2D_plot(EvenlyDiscrXYZ_DataSet xyzData, String title,
			String xAxisLabel, String yAxisLabel, String zAxisLabel, String fileNamePrefix, boolean popupWindow) {
		CPT cpt=null;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(xyzData.getMinZ(), xyzData.getMaxZ());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		XYZPlotSpec spec = new XYZPlotSpec(xyzData, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
		
		try {
			if(popupWindow) {
				XYZPlotWindow window = new XYZPlotWindow(spec, new Range(xyzData.getMinX(),xyzData.getMaxX()), new Range(xyzData.getMinY(),xyzData.getMaxY()));
				XYZGraphPanel panel = window.getXYZPanel();
				if(fileNamePrefix != null) {
					panel.saveAsPDF(fileNamePrefix+".pdf");
					panel.saveAsPNG(fileNamePrefix+".png");
				}
			}
			else if (fileNamePrefix != null) {
				XYZGraphPanel xyzGP = new XYZGraphPanel();
				xyzGP.drawPlot(spec, false, false, new Range(xyzData.getMinX(),xyzData.getMaxX()), new Range(xyzData.getMinY(),xyzData.getMaxY()));
				xyzGP.getChartPanel().setSize(700, 800);
				xyzGP.saveAsPNG(fileNamePrefix+".png");
				xyzGP.saveAsPDF(fileNamePrefix+".pdf");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}



	
	/**
	 * @param args
	 */
	public static void main(String []args) {
		
	}

}
