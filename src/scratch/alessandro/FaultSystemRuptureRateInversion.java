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
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.PlotColorAndLineTypeSelectorControlPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
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
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import scratch.UCERF3.simulatedAnnealing.SerialSimulatedAnnealing;

/**
 * This class does an inversion for the rate of events in an unsegmented fault model:
 * 
 * TO DO:
 * 
 * 1) Add a "redoInverion() method that uses existing constraints?
 * 3) Input a-priori rup-rate constraints
 * 4) Make slip model and enum (already in U3?)
 * 5) Input prob visible model rather than computing here  (already in U3?)
 * 6) sample MRIs via monte carlo simulations (same for slip rates?) for more epistemic uncertainty (or do this outside with zero errors)
 *
 */
public class FaultSystemRuptureRateInversion {
	
	final static boolean D = true;	// debugging flag

	public final static double GAUSS_MFD_SIGMA = 0.12;
	public final static double GAUSS_MFD_TRUNCATION = 2;
	public final static double MAG_DELTA = 0.1;
	
	
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
	private double[] rateOfRupEndsOnSect;

	double totMoRate;
	
	IncrementalMagFreqDist mfdConstraint; 

	private SummedMagFreqDist aveOfSectPartMFDs;
	
	// this accounts for fact that ave slip from Gauss MFD is greater than the slip of the average mag
	private double gaussMFD_slipCorr;
	
	// a-priori rate constraints
	int[] aPriori_rupIndex;
	double[] aPriori_rate, aPriori_wt;
	
	private static boolean MATLAB_TEST = false;
	double[][] C_wted, C;	// inversion matrices
	double[] d, d_wted, data_wt, full_wt, d_pred;  // the data vector
	
	private double minRupRateArray[]; // the minimum rate constraint for each rupture
	private double minRupRate;
	private boolean wtedInversion, applyProbVisible;	// weight the inversion according to slip rate and segment rate uncertainties
	private double relativeSectRateWt, relative_aPrioriRupWt, relative_smoothnessWt; //, relative_aPrioriSegRateWt;
	private double relativeMFD_constraintWt;

	private int  firstRowSectSlipRateData=-1,firstRowSectEventRateData=-1, firstRowAprioriData=-1, firstRowSmoothnessData=-1;
	private int  lastRowSectSlipRateData=-1,lastRowSectEventRateData=-1, lastRowAprioriData=-1, lastRowSmoothnessData=-1;
	private int firstRowGR_constraintData=-1, lastRowGR_constraintData=-1; //, firstRowParkSegConstraint=-1, lastRowParkSegConstraint=-1;
	private int totNumRows;
	
	// slip model:  CHANGE TO ENUM
	private String slipModelType;
	public final static String CHAR_SLIP_MODEL = "Characteristic (Dsr=Ds)";
	public final static String UNIFORM_SLIP_MODEL = "Uniform/Boxcar (Dsr=Dr)";
	public final static String WG02_SLIP_MODEL = "WGCEP-2002 model (Dsr prop to Vs)";
	public final static String TAPERED_SLIP_MODEL = "Tapered Ends ([Sin(x)]^0.5)";
	
	private static EvenlyDiscretizedFunc taperedSlipPDF, taperedSlipCDF;
	
	// MFDs
	SummedMagFreqDist magFreqDist, meanMagHistorgram, magHistorgram;
	ArrayList<SummedMagFreqDist> sectNucleationMFDs;
	ArrayList<SummedMagFreqDist> sectParticipationMFDs;
	
	
	
	private double[] finalSectEventRate, finalPaleoVisibleSectEventRate, finalSectSlipRate;
	private double[] segSlipRateResids, segEventRateResids;
	
	// the following is the total moment-rate reduction, including that which goes to the  
	// background, afterslip, events smaller than the min mag here, and aftershocks and foreshocks.
	private double moRateReduction;  
	
	private MagAreaRelationship magAreaRel;
	
	// NNLS inversion solver - static to save time and memory
	private static NNLSWrapper nnls = new NNLSWrapper();

	public FaultSystemRuptureRateInversion() {
		
		// compute slip correction for Gaussian MFD
		GaussianMagFreqDist gDist1 = new GaussianMagFreqDist(5.0,9.0,41,7,this.GAUSS_MFD_SIGMA,1.0,this.GAUSS_MFD_TRUNCATION,2);
		GaussianMagFreqDist gDist2 = new GaussianMagFreqDist(5.0,9.0,41,7,0.01,1.0,0.01,2);
		gDist1.scaleToCumRate(0, 1.0);
		gDist2.scaleToCumRate(0, 1.0);
		gaussMFD_slipCorr = gDist1.getTotalMomentRate()/gDist2.getTotalMomentRate();
//		System.out.println("gaussMFD_slipCorr="+(float)gaussMFD_slipCorr+"\n");

	}

	
	
	
	
	

	/**
	 * This writes the segPartMFDs to a file and, optionally, makes a plot of the
	 * result (plot not yet saved)
	 * @param dirName
	 * @param gmt_plot
	 */
	public void writeAndPlotSegPartMFDs(String dirName, boolean makePlot) {
		
		// this writes out 
		try{
			FileWriter fw = new FileWriter(dirName+"/segPartMFDsData.txt");
			FileWriter cfw = new FileWriter(dirName+"/segPartCumMFDsData.txt");
			fw.write("seg_index\tmag\tpart_rate\n");
			cfw.write("seg_index\tmag\tpart_cum_rate\n");
			SummedMagFreqDist mfd = sectParticipationMFDs.get(0);
			EvenlyDiscretizedFunc cmfd = mfd.getCumRateDist();
			EvenlyDiscrXYZ_DataSet xyzDataSectPart = new EvenlyDiscrXYZ_DataSet(sectParticipationMFDs.size(), mfd.size(), 0, mfd.getMinX(), 1.0 , mfd.getDelta());
			EvenlyDiscrXYZ_DataSet xyzDataSectPartCum = new EvenlyDiscrXYZ_DataSet(sectParticipationMFDs.size(), cmfd.size(), 0, cmfd.getMinX(), 1, cmfd.getDelta());
			for(int s=0; s<sectParticipationMFDs.size(); s++) {
				mfd = sectParticipationMFDs.get(s);
				cmfd = mfd.getCumRateDist();
				for(int m=0; m<mfd.size();m++) {
					if( mfd.getY(m) != 0.0) {
						fw.write(s+"\t"+(float)mfd.getX(m)+"\t"+(float)Math.log10(mfd.getY(m))+"\n");
						xyzDataSectPart.set(s, mfd.getX(m), Math.log10(mfd.getY(m)));
					}
					else {
						xyzDataSectPart.set(s, mfd.getX(m), Double.NaN);
					}

					if( cmfd.getY(m) != 0.0) {
						cfw.write(s+"\t"+(float)cmfd.getX(m)+"\t"+(float)Math.log10(cmfd.getY(m))+"\n");
						xyzDataSectPartCum.set(s, cmfd.getX(m), Math.log10(cmfd.getY(m)));
					}
					else {
						xyzDataSectPartCum.set(s, cmfd.getX(m), Double.NaN);
					}
				}
			}	
			
			if(makePlot) {
				make2D_plot(xyzDataSectPart, "Incr. Participation", "Section", "Magnitude", "Rate");
				make2D_plot(xyzDataSectPartCum, "Cum. Participation", "Section", "Magnitude", "CumRate");				
			}

			fw.close();
			cfw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		
		
//		}
			
	}

	
	
	
	
	
	/**
	 * This writes and optionally plots the rates and slip-rate contribution of each rupture that has a final rate above the minimum.
	 * @param dirName
	 * @param gmt_plot
	 */
	public void writeAndPlotNonZeroRateRups(String dirName, boolean makePlot) {
		
		// get number of ruptures above minRupRate
		int numAboveMinRate = 0;
		for(int rup=0; rup<numRuptures;rup++)
			if(rupRateSolution[rup]>minRupRate)
				numAboveMinRate += 1;

		EvenlyDiscrXYZ_DataSet xyzDataRupSlipRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);
		EvenlyDiscrXYZ_DataSet xyzDataRupRate = new EvenlyDiscrXYZ_DataSet(numSections,numAboveMinRate, 0, 0, 1, 1);

		int index = 0;
		try{
			FileWriter fw = new FileWriter(dirName+"/rupSlipRateData.txt");
			FileWriter fw2 = new FileWriter(dirName+"/rupRateData.txt");
			fw.write("index\tseg_index\tslip_rate\n");
			fw2.write("index\tseg_index\trate\n");
			for(int rup=0; rup<numRuptures;rup++) {
				if(rupRateSolution[rup]>minRupRate) {
					double logRate = Math.log10(rupRateSolution[rup]);
					for(int s=0; s<numSections; s++) {
						if(rupSectionMatrix[s][rup] == 1) {
							double logSlipRate = Math.log10(sectSlipInRup[s][rup]*rupRateSolution[rup]);
							fw.write(index+"\t"+s+"\t"+(float)logSlipRate+"\n");
							fw2.write(index+"\t"+s+"\t"+(float)logRate+"\n");
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

			if(makePlot) {
				make2D_plot(xyzDataRupSlipRate, "Slip Rate from Ruptures", "Section", "Rupture", "Log10 Slip Rate");
				make2D_plot(xyzDataRupRate, "Rate of Ruptures", "Section", "Rupture", "Log10 Rate");	
			}

			fw.close();
			fw2.close();

		}catch(Exception e) {
			e.printStackTrace();
		}			
	}

	
	

	/**
	 *
	 * @param fltSectionDataList
	 * @param sectionRateConstraints
	 * @param rupSectionMatrix
	 * @param slipModelType - TAPERED_SLIP_MODEL, UNIFORM_SLIP_MODEL, or WG02_SLIP_MODEL
	 * @param magAreaRel
	 * @param relativeSectRateWt - weight on segment rate equations (relative to slip rate)
	 * @param relative_aPrioriRupWt - weight on a-priori rates (relative to slip rate)
	 * @param relative_smoothnessWt - weight on smoothness equations (relative to slip rate)
	 * @param wtedInversion - apply data uncertainties?
	 * @param minRupRate - constrain all rupture rates to be greater than this value
	 * @param applyProbVisible - account for likelihood that Paleoseismology will see the rupture
	 * @param moRateReduction - fraction reduction from smaller events (and not aseismicity or coupling coefficient, which are set in fltSectionDataList)
	 * @param relativeMFD_constraintWt - weight for MFD constraint
	 */
	public void doInversion(
			ArrayList<FaultSectionPrefData> fltSectionDataList, 
			ArrayList<SegRateConstraint> sectionRateConstraints,
			int[][] rupSectionMatrix, 
			String slipModelType, 
			MagAreaRelationship magAreaRel, 
			double relativeSectRateWt, 
			double relative_aPrioriRupWt, 
			double relative_smoothnessWt, 
			boolean wtedInversion, 
			double minRupRate,
			boolean applyProbVisible, 
			double moRateReduction, 
			IncrementalMagFreqDist mfdConstraint,
			double relativeMFD_constraintWt) {
		
		this.fltSectionDataList = fltSectionDataList;
		this.sectionRateConstraints = sectionRateConstraints;
		this.rupSectionMatrix = rupSectionMatrix;
		this.slipModelType = slipModelType;
		this.magAreaRel = magAreaRel;
		this.relativeSectRateWt = relativeSectRateWt;
		this.relative_aPrioriRupWt = relative_aPrioriRupWt;
		this.relative_smoothnessWt = relative_smoothnessWt;
		this.wtedInversion = wtedInversion;
		this.minRupRate = minRupRate;
		this.applyProbVisible = applyProbVisible;
		this.moRateReduction = moRateReduction;
		this.mfdConstraint = mfdConstraint;
		this.relativeMFD_constraintWt = relativeMFD_constraintWt;

		
		// initialize section and rupture attributes
		initSectAndRupAttributes();
		
		// set the a-priori rup rates & wts now that mags are set THESE NEED TO BE INPUT
		aPriori_rupIndex = new int[0];
		aPriori_rate = new double[0];
		

		// compute matrix of Dsr (slip on each segment in each rupture)
		computeSectSlipInRupMatrix();
		
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
		
		// get the number of a-priori rate constraints
		int num_aPriori_constraints = aPriori_rupIndex.length;
		
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
			setApriorRupRatesFromMFD_Constrint();
			num_aPriori_constraints = aPriori_rupIndex.length;
			firstRowAprioriData  = totNumRows;
			totNumRows += num_aPriori_constraints;
			lastRowAprioriData = totNumRows-1;
		}
		
		// add number of smoothness constraints
		if(relative_smoothnessWt > 0) {
			firstRowSmoothnessData=totNumRows;
			if(minNumSectInRup==1)
				totNumRows += numRuptures-numSections;
			else if(minNumSectInRup==2) // the case where numSegForSmallestRups=2
				totNumRows += numRuptures-(numSections-1);
			else
				throw new RuntimeException("numSegForSmallestRups="+minNumSectInRup+" not supported");
			lastRowSmoothnessData = totNumRows-1;
		}
		
		// add number of MFD constaints
		int numMFD_constraints=0;
		if(relativeMFD_constraintWt>0) {
			numMFD_constraints = mfdConstraint.size();
			firstRowGR_constraintData = totNumRows;
			totNumRows += numMFD_constraints;
			lastRowGR_constraintData = totNumRows-1;
		}
				
		System.out.println("\nfirstRowSegEventRateData="+firstRowSectEventRateData+
				";\tfirstRowAprioriData="+firstRowAprioriData+
				";\tfirstRowSmoothnessData="+firstRowSmoothnessData+
				";\tfirstRowGR_constraintData="+firstRowGR_constraintData+
				";\ttotNumRows="+totNumRows+"\n");
			
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
		
		// add the smoothness constraint
		if(relative_smoothnessWt > 0.0) {
			int row = firstRowSmoothnessData;
			int counter = 0;
			for(int rup=0; rup < numRuptures; rup++) {
				// check to see if the last segment is used by the rupture (skip this last rupture if so)
				if(rupSectionMatrix[numSections-1][rup] != 1) {
					d[row] = 0;
					C[row][rup]=1.0;
					C[row][rup+1]=-1.0;
					row += 1;
					counter +=1;
				}
//				else
//					System.out.println("REJECT: row="+row+"\trup="+rup+"\tcounter="+counter);
			}
		}
//		System.out.println("num_smooth_constrints="+num_smooth_constrints);
		
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
		// smoothness constraint wts
		if(relative_smoothnessWt > 0.0) {
			int row = firstRowSmoothnessData;
			for(int rup=0; rup < numRuptures; rup++) {
				// check to see if the last segment is used by the rupture (skip this last rupture if so)
				if(rupSectionMatrix[numSections-1][rup] != 1) {
					full_wt[row] = relative_smoothnessWt;
					d_wted[row] *= full_wt[row];
					C_wted[row][rup] *= full_wt[row];
					C_wted[row][rup+1] *= full_wt[row];
					row += 1;
				}
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
		

		
//		for(int row=0;row<totNumRows; row++)
//			System.out.println(row+"\t"+(float)d[row]);

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
//		rupRateSolution = getNNLS_solution(C_wted, d_wted);

		// set initial state from MFD constraint
		setApriorRupRatesFromMFD_Constrint();
		double[] initialState = aPriori_rate;
		// or set initial state to zero
//		double[] initialState = new double[numRuptures];
		rupRateSolution = getSimulatedAnnealingSolution(C_wted, d_wted, initialState);

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
		
	}
	
		
	
	

	
	/**
	 * This creates the segSlipInRup (Dsr) matrix based on the value of slipModelType.
	 * This slips are in meters.
	 *
	 */
	private void computeSectSlipInRupMatrix() {
		sectSlipInRup = new double[numSections][numRuptures];
		FaultSectionPrefData segData;
		
		// for case segment slip is independent of rupture (constant), and equal to slip-rate * MRI
		if(slipModelType.equals(CHAR_SLIP_MODEL)) {
			throw new RuntimeException(CHAR_SLIP_MODEL+ " not yet supported");
		}
		
		// the rest get average slip from mag-area relationship
		rupAveSlip = new double[numRuptures];
		for(int rup=0; rup<numRuptures; ++rup)
			rupAveSlip[rup] = rupMeanMo[rup]/(rupArea[rup]*FaultMomentCalc.SHEAR_MODULUS);  // inlcudes aveSlipCor
		
		// for case where ave slip computed from mag & area, and is same on all segments 
		if (slipModelType.equals(UNIFORM_SLIP_MODEL)) {
			for(int rup=0; rup<numRuptures; ++rup) {
				for(int seg=0; seg<numSections; seg++) {
					sectSlipInRup[seg][rup] = rupSectionMatrix[seg][rup]*rupAveSlip[rup];
				}
			}
		}
		// this is the model where seg slip is proportional to segment slip rate 
		// (bumped up or down based on ratio of seg slip rate over wt-ave slip rate (where wts are seg areas)
		else if (slipModelType.equals(WG02_SLIP_MODEL)) {
			for(int rup=0; rup<numRuptures; ++rup) {
				double totMoRate = 0;	// a proxi for slip-rate times area
				double totArea = 0;
				for(int seg=0; seg<numSections; seg++) {
					if(rupSectionMatrix[seg][rup]==1) {
						totMoRate += sectMoRate[seg]; // a proxi for Vs*As
						totArea += sectArea[seg];
					}
				}
				for(int seg=0; seg<numSections; seg++) {
					sectSlipInRup[seg][rup] = rupAveSlip[rup]*rupSectionMatrix[seg][rup]*sectMoRate[seg]*totArea/(totMoRate*sectArea[seg]);
				}
			}
		}
		else if (slipModelType.equals(TAPERED_SLIP_MODEL)) {
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
	
	
	/**
	 * This gets the smoothness prediction error regardless of whether it was used in the inversion
	 * @return
	 */
	private double getTotSmoothnessPredErr() {
		double predErr=0;	
		for(int rup=0; rup < numRuptures; rup++) {
			// check to see if the last segment is used by the rupture (skip this last rupture if so)
			if(rupSectionMatrix[numSections-1][rup] != 1) {
				predErr += (rupRateSolution[rup]-rupRateSolution[rup+1])*(rupRateSolution[rup]-rupRateSolution[rup+1]);
			}
		}
		return predErr;
	}

	
private static double[] getSimulatedAnnealingSolution(double[][] C, double[] d, double[] initialState) {
	SparseDoubleMatrix2D matrixC = new SparseDoubleMatrix2D(C); //
	SerialSimulatedAnnealing simulatedAnnealing =new SerialSimulatedAnnealing(matrixC, d, initialState);
	long numIterations = (long)1000000;
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
					finalPaleoVisibleSectEventRate[seg]+=rupRateSolution[rup]*getProbVisible(rupMeanMag[rup]);
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
		System.out.println("Moment Rates: from Fault Sections (possible reduced) = "+ (float)origMoRate+
				";  from total MFD = "+(float)magFreqDist.getTotalMomentRate()+
				",  ratio = "+(float)ratio+"\n");
		
		
		// COMPUTE RATE AT WHICH SECTION BOUNDARIES CONSTITUTE RUPTURE ENDPOINTS
		rateOfRupEndsOnSect = new double[numSections+1];  // there is one more boundary than segments
		for(int rup=0; rup<numRuptures;rup++) {
			int beginBoundary = firstSectOfRup[rup];
			int endBoundary = firstSectOfRup[rup]+numSectInRup[rup];
			rateOfRupEndsOnSect[beginBoundary] += rupRateSolution[rup];
			rateOfRupEndsOnSect[endBoundary] += rupRateSolution[rup];
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
	
	public void writeFinalStuff() {

		// write out rupture rates and mags
//		System.out.println("Final Rupture Rates & Mags:");
//		for(int rup=0; rup < num_rup; rup++)
//		System.out.println(rup+"\t"+(float)rupRateSolution[rup]+"\t"+(float)rupMeanMag[rup]);
		
		// WRITE OUT MAXIMUM EVENT RATE & INDEX
		double maxRate=0;
		int index=-1;
		for(int seg = 0; seg < numSections; seg ++)
			if(finalSectEventRate[seg] > maxRate) {
				maxRate = finalSectEventRate[seg];
				index = seg;
			}
		int mri = (int)Math.round(1.0/maxRate);
		System.out.println("\nMax seg rate (MRI): "+(float)maxRate+" ("+mri+" yrs) at index "+index);



		//write out number of ruptures that have rates above minRupRate
		int numAbove = 0;
		for(int rup=0; rup<this.rupRateSolution.length; rup++)
			if(rupRateSolution[rup] > minRupRate) numAbove += 1;
		System.out.println("\nNum Ruptures above minRupRate = "+numAbove+"\t(out of "+rupRateSolution.length+")\n");

		// write out final segment slip rates
//		System.out.println("\nSegment Slip Rates: index, final, orig, and final/orig (orig is corrected for moRateReduction)");
		double aveRatio=0;
		for(int seg = 0; seg < numSections; seg ++) {
			double slipRate = sectSlipRate[seg]*(1-this.moRateReduction);
			aveRatio += finalSectSlipRate[seg]/slipRate;
//			System.out.println(seg+"\t"+(float)finalSegSlipRate[seg]+"\t"+(float)slipRate+"\t"+(float)(finalSegSlipRate[seg]/slipRate));
		}
		aveRatio /= numSections;
		System.out.println("Ave final/orig slip rate = "+(float)aveRatio);


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
			System.out.println("Ave final/orig segment rate = "+(float)aveRatio);
		}

		// write out final rates for ruptures with an a-priori constraint
		if(this.relative_aPrioriRupWt >0.0)
			aveRatio = 0;
//		System.out.println("\nA Priori Rates: index, final, orig, and final/orig");
		for(int i=0; i<this.aPriori_rate.length;i++) {
			double ratio;
			if(rupRateSolution[aPriori_rupIndex[i]] > 1e-14 && aPriori_rate[i] > 1e-14)  // if both are not essentially zero
				ratio = (rupRateSolution[aPriori_rupIndex[i]]/aPriori_rate[i]);
			else
				ratio = 1;
			aveRatio += ratio;
//			System.out.println(aPriori_rupIndex[i]+"\t"+(float)rupRateSolution[aPriori_rupIndex[i]]+"\t"+aPriori_rate[i]+"\t"+(float)ratio);				
		}
		aveRatio /= aPriori_rate.length;
		System.out.println("Ave final/orig a-priori rate = "+(float)aveRatio);

	}
	
	public void writePredErrorInfo() {
		
		// First without equation weights
		double totPredErr=0, slipRateErr=0, eventRateErr=0, aPrioriErr=0, smoothnessErr=0, grErr=0;
		
		for(int row=firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row++)
			slipRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relativeSectRateWt >0)
			for(int row=firstRowSectEventRateData; row <= lastRowSectEventRateData; row++)
				eventRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relative_aPrioriRupWt > 0)
			for(int row=firstRowAprioriData; row <= lastRowAprioriData; row++)
				aPrioriErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relative_smoothnessWt>0)
			for(int row=firstRowSmoothnessData; row <= lastRowSmoothnessData; row++)
				smoothnessErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
		if(relativeMFD_constraintWt>0)
			for(int row=firstRowGR_constraintData; row <= lastRowGR_constraintData; row++){
				grErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
//double err = (d[row]-d_pred[row])*(d[row]-d_pred[row])*data_wt[row]*data_wt[row];
//System.out.println(row+"\t"+(float)d_pred[row]);
			}

		totPredErr = slipRateErr+eventRateErr+aPrioriErr+smoothnessErr+grErr;
		
		// get alt versions in case equation wts zero
		double eventRateErrAlt = getTotEventRatePredErr();
		double smoothnessErrAlt = getTotSmoothnessPredErr();
		
		double aveSRnormResid = Math.sqrt(slipRateErr/numSections);
		double aveERnormResid = Math.sqrt(eventRateErrAlt/sectionRateConstraints.size());
		
		String resultsString = "\nTotal Prediction Error =\t"+(float)totPredErr+"\n\t"+
				               "Slip Rate Err =\t\t"+(float)slipRateErr+"\trel. wt = 1.0"+
				               "\tnorm resid rms = "+(float)aveSRnormResid+"\n\t"+
				               "Event Rate Err =\t"+(float)eventRateErrAlt+"\trel. wt = "+relativeSectRateWt+
				               "\tnorm resid rms = "+(float)aveERnormResid+"\n\t"+
				               "A Priori Err =\t\t"+(float)aPrioriErr+"\trel. wt = "+relative_aPrioriRupWt+"\n\t"+
				               "Smoothness Err =\t"+(float)smoothnessErrAlt+"\trel. wt = "+relative_smoothnessWt+"\n\t"+
        					   "GR Err =\t"+(float)grErr+"\trel. wt = "+relativeMFD_constraintWt+"\n";

		System.out.println(resultsString);
				
		// Now with equation weights
		totPredErr=0; slipRateErr=0; eventRateErr=0; aPrioriErr=0; smoothnessErr=0; grErr=0;
		for(int row=firstRowSectSlipRateData; row <= lastRowSectSlipRateData; row++)
			slipRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relativeSectRateWt >0)
			for(int row=firstRowSectEventRateData; row <= lastRowSectEventRateData; row++)
				eventRateErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relative_aPrioriRupWt > 0)
			for(int row=firstRowAprioriData; row <= lastRowAprioriData; row++)
				aPrioriErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relative_smoothnessWt>0)
			for(int row=firstRowSmoothnessData; row <= lastRowSmoothnessData; row++)
				smoothnessErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];
		if(relativeMFD_constraintWt>0)
			for(int row=firstRowGR_constraintData; row <= lastRowGR_constraintData; row++)
				grErr += (d[row]-d_pred[row])*(d[row]-d_pred[row])*full_wt[row]*full_wt[row];

		totPredErr = slipRateErr+eventRateErr+aPrioriErr+smoothnessErr+grErr;

		System.out.println("\nTotal Pred Err w/ Eq Wts =\t"+(float)totPredErr+"\n\t"+
				"Slip Rate Err =\t\t"+(float)slipRateErr+"\trel. wt = 1.0\n\t"+
				"Event Rate Err =\t"+(float)eventRateErr+"\trel. wt = "+relativeSectRateWt+"\n\t"+
				"A Priori Err =\t\t"+(float)aPrioriErr+"\trel. wt = "+relative_aPrioriRupWt+"\n\t"+
				"Smoothness Err =\t"+(float)smoothnessErr+"\trel. wt = "+relative_smoothnessWt+"\n\t"+
				"GR Err =\t"+(float)grErr+"\trel. wt = "+relativeMFD_constraintWt+"\n");
			
	}
	

	public void plotStuff(String dirName) {
		
		// plot orig and final slip rates	
		double min = 0, max = numSections-1;
		EvenlyDiscretizedFunc origSlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc origUpper95_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc origLower95_SlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		EvenlyDiscretizedFunc finalSlipRateFunc = new EvenlyDiscretizedFunc(min, max, numSections);
		for(int s=0; s<numSections;s++) {
			origSlipRateFunc.set(s,sectSlipRate[s]*(1-moRateReduction));
			origUpper95_SlipRateFunc.set(s,(sectSlipRate[s]+1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
			origLower95_SlipRateFunc.set(s,(sectSlipRate[s]-1.96*sectSlipRateStdDev[s])*(1-moRateReduction));
			finalSlipRateFunc.set(s,finalSectSlipRate[s]);
		}
		ArrayList sr_funcs = new ArrayList();
		origSlipRateFunc.setName("Orig Slip Rates");
		origUpper95_SlipRateFunc.setName("Orig Upper 95% Confidence for Slip Rates");
		origLower95_SlipRateFunc.setName("Orig Lower 95% Confidence for Slip Rates");
		finalSlipRateFunc.setName("Final Slip Rates");
		sr_funcs.add(finalSlipRateFunc);
		sr_funcs.add(origSlipRateFunc);
		sr_funcs.add(origUpper95_SlipRateFunc);
		sr_funcs.add(origLower95_SlipRateFunc);
		GraphWindow sr_graph = new GraphWindow(sr_funcs, "");  
		ArrayList<PlotCurveCharacterstics> sr_plotChars = new ArrayList<PlotCurveCharacterstics>();
		sr_plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		sr_plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		sr_graph.setPlotChars(sr_plotChars);
		sr_graph.setX_AxisLabel("Section");
		sr_graph.setY_AxisLabel("Slip Rate (m/sec)");
//		sr_graph.setY_AxisRange(0.0, 0.04);
		sr_graph.setTickLabelFontSize(12);
		sr_graph.setAxisLabelFontSize(14);
		if(dirName != null) {
			String filename = dirName+"/slipRates";
			try {
				sr_graph.saveAsPDF(filename+".pdf");
				sr_graph.saveAsPNG(filename+".png");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		

		// plot orig and final section event rates	
		ArrayList er_funcs = new ArrayList();
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
		ArbitrarilyDiscretizedFunc func;
		ArrayList obs_er_funcs = new ArrayList();
		SegRateConstraint constraint;
		for(int c=0;c<num;c++) {
			func = new ArbitrarilyDiscretizedFunc();
			constraint = sectionRateConstraints.get(c);
			int seg = constraint.getSegIndex();
			func.set((double)seg-0.0001, constraint.getLower95Conf());
			func.set((double)seg, constraint.getMean());
			func.set((double)seg+0.0001, constraint.getUpper95Conf());
			func.setName(constraint.getFaultName());
			obs_er_funcs.add(func);
			er_funcs.add(func);
		}
//		er_funcs.add(obs_er_funcs);
		GraphWindow er_graph = new GraphWindow(er_funcs, ""); 
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		for(int c=0;c<num;c++)
			plotChars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 1f, PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
		er_graph.setPlotChars(plotChars);
		er_graph.setX_AxisLabel("Section");
		er_graph.setY_AxisLabel("Event Rate (per yr)");
		er_graph.setYLog(true);
		er_graph.setTickLabelFontSize(12);
		er_graph.setAxisLabelFontSize(14);
		if(dirName != null) {
			String filename = dirName+"/eventRates";
			try {
				er_graph.saveAsPDF(filename+".pdf");
				er_graph.saveAsPNG(filename+".png");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		// plot the final rupture rates
		max = numRuptures-1;
		EvenlyDiscretizedFunc rupRateFunc = new EvenlyDiscretizedFunc(min, max, numRuptures);
		for(int rup=0; rup<numRuptures;rup++) {
			rupRateFunc.set(rup,rupRateSolution[rup]);
		}
		ArrayList rup_funcs = new ArrayList();
		rupRateFunc.setName("Rupture Rates");
		rup_funcs.add(rupRateFunc);
		GraphWindow rup_graph = new GraphWindow(rup_funcs, "Rupture Rates");   


		// PLOT MFDs
		ArrayList mfd_funcs = new ArrayList();
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
		

		GraphWindow mfd_graph = new GraphWindow(mfd_funcs, "Magnitude Frequency Distributions");   
		mfd_graph.setYLog(true);
		mfd_graph.setY_AxisRange(1e-7, 0.2);
		mfd_graph.setX_AxisRange(minRupMagWithAleatory, maxRupMagWithAleatory);
		mfd_graph.setX_AxisLabel("Magnitude");
		mfd_graph.setY_AxisLabel("Rate (per yr)");

		ArrayList<PlotCurveCharacterstics> plotMFD_Chars = new ArrayList<PlotCurveCharacterstics>();
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.MAGENTA));
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		plotMFD_Chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK));
		mfd_graph.setPlotChars(plotMFD_Chars);
		mfd_graph.setTickLabelFontSize(12);
		mfd_graph.setAxisLabelFontSize(14);
		if(dirName != null) {
			String filename = dirName+"/MFDs";
			try {
				mfd_graph.saveAsPDF(filename+".pdf");
				mfd_graph.saveAsPNG(filename+".png");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}



		
		
		// PLOT RATE AT WHICH SECTIONS ENDS REPRESENT RUPTURE ENDS
		min = 0; max = numSections;
		EvenlyDiscretizedFunc rateOfRupEndsOnSegFunc = new EvenlyDiscretizedFunc(min, max, numSections+1);
		for(int seg=0; seg<numSections+1;seg++) {
			rateOfRupEndsOnSegFunc.set(seg,rateOfRupEndsOnSect[seg]);
		}
		ArrayList seg_funcs = new ArrayList();
		rateOfRupEndsOnSegFunc.setName("Rate that section ends represent rupture ends");
		seg_funcs.add(finalSlipRateFunc);
		seg_funcs.add(finalPaleoVisibleEventRateFunc);
		// add paleoseismic obs
		for(int i=0; i<obs_er_funcs.size();i++) seg_funcs.add(obs_er_funcs.get(i));
		seg_funcs.add(rateOfRupEndsOnSegFunc);
		
		GraphWindow seg_graph = new GraphWindow(seg_funcs, ""); 
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		for(int c=0;c<num;c++)
			plotChars2.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 1f, PlotSymbol.FILLED_CIRCLE, 4f, Color.RED));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		seg_graph.setPlotChars(plotChars2);
		seg_graph.setX_AxisLabel("Subsection");
		seg_graph.setY_AxisLabel("Rates");
		seg_graph.setYLog(true);
		seg_graph.setTickLabelFontSize(12);
		seg_graph.setAxisLabelFontSize(14);
		if(dirName != null) {
			String filename = dirName+"/endpointRates";
			try {
				seg_graph.saveAsPDF(filename+".pdf");
				seg_graph.saveAsPNG(filename+".png");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}


	}
	
	public void plotMagHistograms() {
		
		ArrayList funcs = new ArrayList();
		funcs.add(meanMagHistorgram);
		funcs.add(magHistorgram);
		GraphWindow mHist_graph = new GraphWindow(funcs, "Mag Histograms");   
//		mfd_graph.setYLog(true);
//		mfd_graph.setY_AxisRange(1e-5, 1);
		mHist_graph.setX_AxisRange(minRupMagWithAleatory,maxRupMagWithAleatory);
		
		// make a numSegInRupHistogram
//		SummedMagFreqDist numSegInRupHistogram = new SummedMagFreqDist(1.0,numSections,1.0);
//		for(int r=0;r<numSectInRup.length;r++) numSegInRupHistogram.add((double)numSectInRup[r], 1.0);
//		numSegInRupHistogram.setName("Num Segments In Rupture Histogram");
//		ArrayList funcs2 = new ArrayList();
//		funcs2.add(numSegInRupHistogram);
//		GraphWindow graph = new GraphWindow(funcs2, "Num Segments In Rupture Histogram");   

	}
	
	
	private GutenbergRichterMagFreqDist getGR_DistFit() {
		int num = (int)Math.round((maxRupMag-minRupMag)/0.1 + 1);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(minRupMag,num,0.1);
		double moRate = totMoRate*(1-moRateReduction);
//		double altMoRate = magFreqDist.getTotalMomentRate();
		gr.setAllButTotCumRate(minRupMag, maxRupMag, moRate, 1.0);
		gr.setName("GR fit");
		if(D) {
			System.out.println("GR MFD FIT:\n\tMmin="+minRupMag+"\n\tMmax="+maxRupMag+
					"\n\tnum="+num+"\n\tdelta="+gr.getDelta()+"\n\tmoRate="+moRate);
		}
		return gr;
	}
	
	private void setApriorRupRatesFromMFD_Constrint() {
		aPriori_rupIndex = new int[numRuptures];
		aPriori_rate  = new double[numRuptures];
		aPriori_wt  = new double[numRuptures];
		for(int r=0;r<numRuptures;r++) {
			aPriori_rupIndex[r] = r;
			aPriori_rate[r] = mfdConstraint.getY(rupMeanMag[r])/meanMagHistorgram.getY(rupMeanMag[r]);
			aPriori_wt[r]=1;
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
	 * 5	0.10
	 * 6	0.45
	 * 7	0.87
	 * 8	0.98
	 * 9	1.00
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


		// compute rupture mean mags etc
		rupMeanMag = new double[numRuptures];
		rupMeanMo = new double[numRuptures];
		minRupMag = Double.MAX_VALUE;
		maxRupMag = 0;
		if(slipModelType.equals(CHAR_SLIP_MODEL))
			// getRupMeanMagsAssumingCharSlip();
			throw new RuntimeException(CHAR_SLIP_MODEL+" is not yet supported");
		else {
			// compute from mag-area relationship
			for(int r=0; r <numRuptures; r++) {
				double mag = magAreaRel.getMedianMag(rupArea[r]*1e-6);
				//round this to nearst 10th unit
				rupMeanMag[r] = ((double)Math.round(10*mag))/10.0;
				rupMeanMo[r] = MagUtils.magToMoment(rupMeanMag[r])*gaussMFD_slipCorr;   // increased if magSigma >0
				if(minRupMag>rupMeanMag[r])
					minRupMag=rupMeanMag[r];
				if(maxRupMag<rupMeanMag[r])
					maxRupMag=rupMeanMag[r];

			}
		}
		double tempMag = minRupMag-GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		minRupMagWithAleatory = ((double)Math.round(10*tempMag))/10.0;
		tempMag = maxRupMag+GAUSS_MFD_SIGMA*GAUSS_MFD_TRUNCATION;
		maxRupMagWithAleatory = ((double)Math.round(10*tempMag))/10.0;
		if(D) System.out.println("maxRupMag="+maxRupMag+"\tminRupMag="+minRupMag+
				"\tminRupMagWithAleatory="+minRupMagWithAleatory+
				"\tmaxRupMagWithAleatory="+maxRupMagWithAleatory);
		
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
	
	public static void make2D_plot(EvenlyDiscrXYZ_DataSet xyzData, String title,
			String xAxisLabel, String yAxisLabel, String zAxisLabel) {
		CPT cpt=null;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(xyzData.getMinZ(), xyzData.getMaxZ());
		} catch (IOException e) {
			e.printStackTrace();
		}
		XYZPlotSpec logLikeSpec = new XYZPlotSpec(xyzData, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
		XYZPlotWindow window_logLikeSpec = new XYZPlotWindow(logLikeSpec, new Range(xyzData.getMinX(),xyzData.getMaxX()), new Range(xyzData.getMinY(),xyzData.getMaxY()));
	}



	
	/**
	 * @param args
	 */
	public static void main(String []args) {
		
	}

}
