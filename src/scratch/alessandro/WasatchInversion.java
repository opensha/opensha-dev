package scratch.alessandro;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.alessandro.logicTreeEnums.ScalingRelationshipEnum;
import scratch.alessandro.logicTreeEnums.SlipAlongRuptureModelEnum;

/**
 * This class reads Wasatch inversion data from files, provides methods for getting the various constraints, and runs the inversion.
 * 
 * Slip rate and event rate standard deviations are computed as the difference of the 95% confidence bounds divided by 4.
 * 
 * @author 
 *
 */
public class WasatchInversion {

	final static boolean D = false;	// debugging flag
	
	public final static String ROOT_PATH = "src/scratch/alessandro/";
	final static String ROOT_DATA_DIR = "src/scratch/alessandro/data/"; // where to find the data

	
	// These values are the same for all fault sections
	final static double UPPER_SEIS_DEPTH = 0;
	final static double LOWER_SEIS_DEPTH = 15;
	final static double FAULT_DIP = 50;
	final static double FAULT_RAKE = -90;
	
	
	ArrayList<FaultSectionPrefData> faultSectionDataList;
	ArrayList<SegRateConstraint> sectionRateConstraints;
	int[][] rupSectionMatrix;
	
	final static String SLIP_RATE_FILENAME = "sliprate_wasatch_final.txt";
	final static String PALEO_RATE_FILENAME = "paleorate_wasatch_and_95p.txt";
	final static String SEGMENT_BOUNDARY_DATA_FILE = "segmentBoundaryData.txt";
	final static String APRIORI_RUP_RATE_FROM_SECT_CONSTR_FILENAME = "aPrioriRupRatesFromSegmentationConstraints.txt";

	final static String FAULT_TRACE_DIR_NAME = "subsections_traces/";

	public WasatchInversion() {
		
		readData();
		
	}
		
	/**
	 * This read rupture rates from a file
	 * @param fileName - the full path
	 * @return
	 */
	private double[] readRuptureRatesFromFile(String fileName) {
		double rates[] = null;
		int numPts=-1;
		int startLineNum=-1;
		int counter = -1;
		File file = new File(fileName);
		List<String> fileLines;
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			for(String str:fileLines) {
				counter += 1;
				if(str.contains("Points")) {
					String[] split = str.split(": ");
					numPts = Integer.valueOf(split[1]);
				}
				if(str.contains("Data[x,y]"))
						startLineNum = counter +1;
			}
			rates = new double[numPts];
			counter =0;
			for(int i=startLineNum; i<startLineNum+numPts;i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split("\t");
				rates[counter] = Double.parseDouble(split[1]);
//				System.out.println(rates[counter]+"\t"+split[0]);
				counter+=1;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		System.out.println("numPts="+numPts);
//		System.out.println("startLineNum="+startLineNum);
//		System.exit(0);;
		
		return rates;
	}
	
	
		
	private void readData() {
	
		int numSections, numRuptures;
		double[] sectSlipRate, sectSlipRateStdDev;

		// Read Data:

		try {
			File file = new File(ROOT_DATA_DIR+SLIP_RATE_FILENAME);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			numSections = fileLines.size();
			if(D) System.out.println("numSections="+numSections);
			int sectIndex = 0;
			sectSlipRate = new double[numSections];
			sectSlipRateStdDev = new double[numSections];
			for (String line : fileLines) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 4, "Expected 4 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				int testIndex = Integer.valueOf(split[0]);
				if(sectIndex != testIndex)
					throw new RuntimeException("Bad section index; "+sectIndex+" != "+testIndex);
				sectSlipRate[sectIndex] = Double.valueOf(split[1]);
				double low95 = Double.valueOf(split[2]);
				double upp95 = Double.valueOf(split[3]);
				sectSlipRateStdDev[sectIndex] = (upp95-low95)/(2*1.96);	
				if(D) System.out.println(sectIndex+"\t"+sectSlipRate[sectIndex]+"\t"+sectSlipRateStdDev[sectIndex]);
				sectIndex+=1;
			}


			// Now read section rate constraints
			sectionRateConstraints   = new ArrayList<SegRateConstraint>();
			file = new File(ROOT_DATA_DIR+PALEO_RATE_FILENAME);
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 5, "Expected 3 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				sectIndex = Integer.valueOf(split[0]);
				double meanRate = Double.valueOf(split[1]);
				double stdDev = Double.valueOf(split[2]);
				double upp95 = Double.valueOf(split[3]);
				double low95 = Double.valueOf(split[4]);
//				stdDev = (upp95-low95)/4.0; 
				SegRateConstraint sectionRateConstraint = new SegRateConstraint("Section "+split[0]); // Names are not unique!
				sectionRateConstraint.setSegRate(sectIndex, meanRate, stdDev, low95, upp95);
				sectionRateConstraints.add(sectionRateConstraint);
				if (D) System.out.println(sectionRateConstraint.getFaultName()+"\t"+
						sectionRateConstraint.getSegIndex()+"\t"+
						sectionRateConstraint.getMean()+"\t"+
						sectionRateConstraint.getStdDevOfMean()+"\t"+
						sectionRateConstraint.getLower95Conf()+"\t"+
						sectionRateConstraint.getUpper95Conf());
			}


			// Now read section rupture matrix
			file = new File("src/scratch/alessandro/data/Gsr_matrix_new.txt");
			fileLines = Files.readLines(file, Charset.defaultCharset());
			numRuptures = fileLines.size()-1;
			rupSectionMatrix = new int[numSections][numRuptures];
			int rupIndex=-1;
			for (String line : fileLines) {
				if(rupIndex==-1) {
					rupIndex+=1;
					continue;
				}
				//					System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length-1 == numSections, "Number of columns (%s) not consistent with numSections (%s)", split.length, numSections);

				int testRupIndex = Integer.valueOf(split[0]);
				Preconditions.checkState(testRupIndex == rupIndex, "Rup index problem on input file (%s vs %s)", testRupIndex, rupIndex);

				for(int s=0; s<numSections; s++) {
					int sectInRup = Integer.valueOf(split[s+1]);
					if(sectInRup == 1) {
						rupSectionMatrix[s][rupIndex] = 1;
					}
				}

				if(D) {
					System.out.println("Rupture "+rupIndex);
					for(int s=0; s<numSections; s++)
						System.out.print(rupSectionMatrix[s][rupIndex]+" ");
					System.out.print("\n");
				}

				rupIndex+=1;
			}


			// Now read traces and make FaultSectionPrefData for each section
			
			faultSectionDataList = new ArrayList<FaultSectionPrefData>();
			for(int s=0; s<numSections; s++) {
				String traceFileName;
				FaultTrace fltTrace = new FaultTrace("Trace "+s);
				traceFileName = s+".txt";
				
				// read fault trace from file
				file = new File(ROOT_DATA_DIR+FAULT_TRACE_DIR_NAME+traceFileName);
				for (String line : Files.readLines(file, Charset.defaultCharset())) {
					//			System.out.println(line);
					line = line.trim();
					String[] split = line.split("\t");	// tab delimited
					Preconditions.checkState(split.length == 2, "Expected 2 items, got %s", split.length);
					//			System.out.println(split[0]+"\t"+split[1]);
					double lon = Double.valueOf(split[0]);
					double lat = Double.valueOf(split[1]);
					fltTrace.add(new Location(lat,lon,UPPER_SEIS_DEPTH));
				}

				FaultSectionPrefData fltSectData = new FaultSectionPrefData();
				fltSectData.setSectionId(s);
				fltSectData.setSectionName("Section "+s);
				fltSectData.setShortName("Sect"+s);
				fltSectData.setFaultTrace(fltTrace);
				fltSectData.setAveDip(FAULT_DIP);
				fltSectData.setAveUpperDepth(UPPER_SEIS_DEPTH);
				fltSectData.setAveLowerDepth(LOWER_SEIS_DEPTH);
				fltSectData.setAveSlipRate(sectSlipRate[s]);
				fltSectData.setSlipRateStdDev(sectSlipRateStdDev[s]);
				fltSectData.setAveRake(FAULT_RAKE);;
				
				if(D) {
					String str = new String();
					str += "sectionId = "+fltSectData.getSectionId()+"\n";
					str += "sectionName = "+fltSectData.getSectionName()+"\n";
					str += "shortName = "+fltSectData.getShortName()+"\n";
					str += "aveLongTermSlipRate = "+fltSectData.getOrigAveSlipRate()+"\n";
					str += "slipRateStdDev = "+fltSectData.getOrigSlipRateStdDev()+"\n";
					str += "aveDip = "+fltSectData.getAveDip()+"\n";
					str += "aveRake = "+fltSectData.getAveRake()+"\n";
					str += "aveUpperDepth = "+fltSectData.getOrigAveUpperDepth()+"\n";
					str += "aveLowerDepth = "+fltSectData.getAveLowerDepth()+"\n";
					str += "aseismicSlipFactor = "+fltSectData.getAseismicSlipFactor()+"\n";
					str += "couplingCoeff = "+fltSectData.getCouplingCoeff()+"\n";
					str += "dipDirection = "+fltSectData.getDipDirection()+"\n";
					str += "dateOfLastEventMillis = "+fltSectData.getDateOfLastEvent()+"\n";
					str += "slipInLastEvent = "+fltSectData.getSlipInLastEvent()+"\n";
					str += "traceLength = "+fltSectData.getTraceLength()+"\n";
					str += "downDipWidth = "+fltSectData.getOrigDownDipWidth()+"\n";
					str += "area (stirling surface) = "+fltSectData.getStirlingGriddedSurface(1.0).getArea()+"\n";
					str += "faultTrace:\n";
					for(int i=0; i <fltSectData.getFaultTrace().size();i++) {
						Location loc = fltSectData.getFaultTrace().get(i);
						str += "\t"+loc.getLatitude()+", "+loc.getLongitude()+", "+loc.getDepth()+"\n";
					}
					System.out.println(str);
				}
				
				faultSectionDataList.add(fltSectData);
				
			}
			

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

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
	
	
	/**
	 * This writes a-priori rupture rate constraints from "segmentation" boundaries to the
	 * specified file.
	 * @param fileName
	 */
	private void writeApriorRupRatesForSegmentationConstrints() {	
		
		if(rupSectionMatrix == null)
			throw new RuntimeException("must create the rupSectionMatrix before running this method");
		
		// Read Segment Boundary Data:

		int numBoundaries=0;
		int sect1_array[]=null;
		int sect2_array[]=null;
		double wt_array[]=null;

		try {
			File file = new File(ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			numBoundaries = fileLines.size();
			sect1_array = new int[numBoundaries];
			sect2_array = new int[numBoundaries];
			wt_array = new double[numBoundaries];
			if(D) System.out.println("numBoundaries="+numBoundaries);
			for (int i=0; i<numBoundaries; i++) {
				//			System.out.println(line);
				String line = fileLines.get(i).trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 3, "Expected 3 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				sect1_array[i] = Integer.valueOf(split[0]);
				sect2_array[i] = Integer.valueOf(split[1]);
				wt_array[i] = Double.valueOf(split[2]);
			}
		}catch(Exception e) {
			e.printStackTrace();
		}			


		try{
			FileWriter fw = new FileWriter(ROOT_DATA_DIR+APRIORI_RUP_RATE_FROM_SECT_CONSTR_FILENAME);
			int numRuptures = rupSectionMatrix[0].length;
			int numSections = rupSectionMatrix.length;
			System.out.println(numRuptures+"\t"+numSections);

			for(int i=0;i<numBoundaries; i++) {				
				for(int rup=0; rup<numRuptures; ++rup) {
					if(rupSectionMatrix[sect1_array[i]][rup]==1 && rupSectionMatrix[sect2_array[i]][rup]==1) {
						fw.write(rup+"\t0.0\t"+wt_array[i]+"\n");
					}
				}
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
	public void writeAndOrPlotHazardCurve(FaultSystemRuptureRateInversion fltSysRupInversion, String dirName, boolean popupWindow) {
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fltSysRupInversion.getFaultSystemSolution());
		erf.setName("WasatchERF");
		
		// set the forecast duration
		double duration = 1;
		erf.getTimeSpan().setDuration(duration);
		
		// write out parameter values
		if(D) {
			ParameterList paramList = erf.getAdjustableParameterList();
			for(int i=0;i<paramList.size(); i++) {
				Parameter<?> param = paramList.getByIndex(i);
				System.out.println(param.getName()+"\t"+param.getValue());
			}			
		}
		
		// update forecast
		erf.updateForecast();
		
		if(D) System.out.println("NumFaultSystemSources = "+erf.getNumFaultSystemSources());
		
//		for(int r=0;r<erf.getNumFaultSystemSources();r++) {
//			System.out.println(r+"\t"+erf.getNthRupture(r).getProbability());
//		}
		
		// compute hazard curve
		// chose attenuation relationship (GMPE)
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGA_Param.NAME);
		// make the site object and set values
		Site site = new Site(new Location(40.75, -111.90));
		for (Parameter<?> param : imr.getSiteParams()) {
			site.addParameter(param);
			System.out.println(param.getName()+"\t"+param.getValue());
		}
					
		ArbitrarilyDiscretizedFunc curveLinearXvalues = IMT_Info.getUSGS_PGA_Function();
//		System.out.println("curveLinearXvalues:\n "+curveLinearXvalues);
		ArbitrarilyDiscretizedFunc curveLogXvalues = HazardCurveSetCalculator.getLogFunction(curveLinearXvalues);
//		System.out.println("curveLogXvalues:\n "+curveLogXvalues);
	
		HazardCurveCalculator calc = new HazardCurveCalculator();
		calc.getHazardCurve(curveLogXvalues, site, imr, erf); // result is put into curveLogXvalues
		
		curveLinearXvalues = HazardCurveSetCalculator.unLogFunction(curveLinearXvalues, curveLogXvalues);
		
		// make the plot
		ArrayList<XY_DataSet> funcs1 = new ArrayList<XY_DataSet>();
		funcs1.add(curveLinearXvalues);
		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		
		String fileNamePrefix = null;
		if(dirName != null)
			fileNamePrefix = dirName+"/hazardCurve";
		String plotName ="Salt Lake City Hazard Curve";
		String xAxisLabel = "PGA";
		String yAxisLabel = "Probability ("+duration+" yr)";
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = true;
		boolean logY = true;

		fltSysRupInversion.writeAndOrPlotFuncs(funcs1, plotChars, plotName, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		

	}

	
	
	/**
	 * @param args
	 */
	public static void main(String []args) {
		
		
		WasatchInversion wasatchInversion = new WasatchInversion();
		// ONLY NEED TO DO ThE FOLLOWING ONCE:
//		wasatchInversion.writeApriorRupRatesForSegmentationConstrints();
//		System.exit(-1);
		ArrayList<FaultSectionPrefData> fltSectDataList = wasatchInversion.getFaultSectionDataList();
		ArrayList<SegRateConstraint> sectionRateConstraints = wasatchInversion.getSectionRateConstraints();
		int[][] rupSectionMatrix = wasatchInversion.getRupSectionMatrix();

		//Write section data (e.g., so it can be checked)
		if(D) {
			String str = "sectID\tslipRate\tslipRateStd\tlength\tdip\trake\tupDepth\tlowDepth\tDDW+\n";
			for(int s=0;s<fltSectDataList.size();s++) {
				FaultSectionPrefData fltSectData = fltSectDataList.get(s);
				str += fltSectData.getSectionId()+"\t";
				str += fltSectData.getOrigAveSlipRate()+"\t";
				str += fltSectData.getOrigSlipRateStdDev()+"\t";
				str += (float)fltSectData.getTraceLength()+"\t";
				str += fltSectData.getAveDip()+"\t";
				str += fltSectData.getAveRake()+"\t";
				str += fltSectData.getOrigAveUpperDepth()+"\t";
				str += fltSectData.getAveLowerDepth()+"\t";
				str += fltSectData.getOrigDownDipWidth()+"\t";
				str += (float)fltSectData.getStirlingGriddedSurface(1.0).getArea()+"\n";
			}
			System.out.println(str);			
		}

//		System.exit(0);
		
		// this will be used to keep track of runtimes
		long startTimeMillis = System.currentTimeMillis();
		if(D)
			System.out.println("Starting Inversion");

		// set inversion attributes
		String name = "Wasatch Inversion";
		SlipAlongRuptureModelEnum slipModelType = SlipAlongRuptureModelEnum.UNIFORM;
		ScalingRelationshipEnum scalingRel = ScalingRelationshipEnum.WC94_SRL_ALL;
		double relativeSectRateWt=1;
		
		double relative_aPrioriRupWt = 0;	//

		boolean wtedInversion = true;
		double minRupRate = 1e-8;
//		double minRupRate = 0.0;
		boolean applyProbVisible = true;
		double moRateReduction = 0.1;	// this is the amount to reduce due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)
		double relativeMFD_constraintWt = 0; // setting this to 1e6
		
		// Segmentation constraint:
		String segmentationConstrFilename = ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE;
		double relative_segmentationConstrWt = 0;
		
		// create an instance of the inversion class with the above settings
		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
				name,
				fltSectDataList, 
				sectionRateConstraints, 
				rupSectionMatrix, 
				slipModelType, 
				scalingRel, 
				relativeSectRateWt, 
				relative_aPrioriRupWt, 
				wtedInversion, 
				minRupRate, 
				applyProbVisible, 
				moRateReduction,
				relativeMFD_constraintWt,
				segmentationConstrFilename,
				relative_segmentationConstrWt);
		
		// make the directory for storing results (set as null if you don't want to save anything)
		String dirName = ROOT_PATH+"OutputFigsAndData";
	    File file = new File(dirName);
	    file.mkdirs();
	    
	    // write the setup info to a file
	    fltSysRupInversion.writeInversionSetUpInfoToFile(dirName);
		
	    
		// NON_NEGATIVE LEAST SQUARES:
//		fltSysRupInversion.doInversionNNLS();
		
//		// SIMULATED ANNEALING
//		long numIterations = (long) 1e4;
//		boolean initStateFromAprioriRupRates = false;
//		long randomSeed = System.currentTimeMillis();
////		long randomSeed = 1525892588112l; // not that the last character here is the letter "l" to indicated a long value
//		fltSysRupInversion.doInversionSA(numIterations, initStateFromAprioriRupRates, randomSeed);
		
		
		// SOLUTION FROM FILE:
		String ratesFileName = ROOT_PATH+"OutputFigsAndData/ruptureRates.txt"; // assumed consistent with values above
		double[] rupRatesArray = wasatchInversion.readRuptureRatesFromFile(ratesFileName);
		fltSysRupInversion.setSolution(rupRatesArray);

		

		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
		if(D) System.out.println("Done with Inversion after "+(float)runTimeSec+" seconds.");
				
		// write results to file
		fltSysRupInversion.writeInversionRunInfoToFile(dirName);
		
		// Now make plots if desired
		boolean popUpPlots = true;	// this tells whether to show plots in a window (turn off for HPC)
//		dirName = null;	// set as null if you don't want to save to file
		fltSysRupInversion.writeAndOrPlotDataFits(dirName, popUpPlots);
		fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots);
		fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots);
	    fltSysRupInversion.writeAndOrPlotSegPartMFDs(dirName, popUpPlots);
	    
	    // hazard curve:
	    wasatchInversion.writeAndOrPlotHazardCurve(fltSysRupInversion, dirName, popUpPlots);

		
	}



}
