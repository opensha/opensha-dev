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
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc_3D;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
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
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.CoolingScheduleType;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.alessandro.logicTreeEnums.ScalingRelationshipEnum;
import scratch.alessandro.logicTreeEnums.SlipAlongRuptureModelEnum;
import scratch.alessandro.logicTreeEnums.WasatchSlipRatesEnum;

/**
 * This class reads Wasatch inversion data from files, provides methods for getting the various constraints, and runs the inversion.
 * 
 * Slip rate and event rate standard deviations are computed as the difference of the 95% confidence bounds divided by 4.
 * 
 * @author 
 *
 */
public class WasatchInversion {

	final static boolean D = true;	// debugging flag
	
	public final static String ROOT_PATH = "src/scratch/alessandro/";
	final static String ROOT_DATA_DIR = "src/scratch/alessandro/data/"; // where to find the data

	
	// These values are the same for all fault sections
	final static double UPPER_SEIS_DEPTH = 0;
	final static double LOWER_SEIS_DEPTH = 15;
	final static double FAULT_DIP = 50;
	final static double FAULT_RAKE = -90;
	
	final static double hazGridSpacing = 0.05;
	
	WasatchSlipRatesEnum wasatchSlipRatesEnum;

	
	
	ArrayList<FaultSectionPrefData> faultSectionDataList;
	ArrayList<SegRateConstraint> sectionRateConstraints;
	int[][] rupSectionMatrix;
	
	final static String SLIP_RATE_FILENAME = "sliprate_wasatch_final.txt";
	final static String PALEO_RATE_FILENAME = "paleorate_wasatch_and_95p_updated.txt";
	final static String SEGMENT_BOUNDARY_DATA_FILE = "segmentBoundaryData.txt";
	final static String APRIORI_RUP_RATE_FROM_SECT_CONSTR_FILENAME = "aPrioriRupRatesFromSegmentationConstraints.txt";

	final static String FAULT_TRACE_DIR_NAME = "subsections_traces/";
	
	public enum InversionSolutionType {
		NON_NEGATIVE_LEAST_SQUARES,
		SIMULATED_ANNEALING,
		FROM_FILE;
	}
	
	double hazCurveLnMin = Math.log(0.001);
	double hazCurveLnMax = Math.log(10);
	int hazCurveNum = 20;
	double hazCurveDelta = (hazCurveLnMax-hazCurveLnMin)/(double)(hazCurveNum-1);

	
	public WasatchInversion(WasatchSlipRatesEnum wasatchSlipRatesEnum) {
		this.wasatchSlipRatesEnum = wasatchSlipRatesEnum;
		readData();
		
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
	
	
		
	private void readData() {
	
		int numSections, numRuptures;
		
		
		
		double[] sectSlipRate = wasatchSlipRatesEnum.getSlipRateArray();
		double[] sectSlipRateStdDev = wasatchSlipRatesEnum.getSlipRateStdomArray();
		numSections = sectSlipRate.length;
		

		// Read Data:

// NO LONGER NEEDED:
//		double[] sectSlipRate, sectSlipRateStdDev;
//		try {
//			File file = new File(ROOT_DATA_DIR+SLIP_RATE_FILENAME);
//			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
//			numSections = fileLines.size();
//			if(D) System.out.println("numSections="+numSections);
//			int sectIndex = 0;
//			sectSlipRate = new double[numSections];
//			sectSlipRateStdDev = new double[numSections];
//			for (String line : fileLines) {
//				//			System.out.println(line);
//				line = line.trim();
//				String[] split = line.split("\t");	// tab delimited
//				Preconditions.checkState(split.length == 4, "Expected 4 items, got %s", split.length);
//				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
//				int testIndex = Integer.valueOf(split[0]);
//				if(sectIndex != testIndex)
//					throw new RuntimeException("Bad section index; "+sectIndex+" != "+testIndex);
//				sectSlipRate[sectIndex] = Double.valueOf(split[1]);
//				double low95 = Double.valueOf(split[2]);
//				double upp95 = Double.valueOf(split[3]);
//				sectSlipRateStdDev[sectIndex] = (upp95-low95)/(2*1.96);	
//				if(D) System.out.println(sectIndex+"\t"+sectSlipRate[sectIndex]+"\t"+sectSlipRateStdDev[sectIndex]);
//				sectIndex+=1;
//			}

		try {
			
			// Read section rate constraints
			sectionRateConstraints   = new ArrayList<SegRateConstraint>();
			File file = new File(ROOT_DATA_DIR+PALEO_RATE_FILENAME);
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 5, "Expected 3 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				int sectIndex = Integer.valueOf(split[0]);
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
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
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
			UncertainArbDiscFunc hazCurveMean95confLnX = fltSysRupInversion.get95perConfForMultRuns(curvesFromMultRunsFunc_3D);

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
			UncertainArbDiscFunc hazCurveMinMaxRange = new UncertainArbDiscFunc(hazCurveMean, hazCurveMin, hazCurveMax);
			hazCurveMinMaxRange.setName("hazCurveMinMaxRange");
			UncertainArbDiscFunc hazCurveMean95conf = new UncertainArbDiscFunc(hazCurveMean, hazCurveMeanLower95, hazCurveMeanUpper95);
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
	
	/**
	 * As written, the FaultSystemRuptureRateInversion constructor used here always sets the a-priori rupture rates and MFD constraint
	 * from a best-fitting GR distribution (whether these are applied depends on the weights given and the Simulated Annealing initial 
	 * solution)
	 * @param args
	 */
	public static void main(String []args) {
		
		// THE FOLLOWING SETS ALL THE INVERSION ATTRIBUTES:
		//------------------------------------------------
		
		String dirName = ROOT_PATH+"OutputDataAndFigs_SA10_Uniform_noER_boxcar_fromZero_segmented";
	
		// Inversion name
		String name = "Wasatch Inversion";
		
		// Slip Rate
		WasatchSlipRatesEnum slipRates = WasatchSlipRatesEnum.UNIFORM_MEAN;
		
		// Average slip along rupture model
//		SlipAlongRuptureModelEnum slipModelType = SlipAlongRuptureModelEnum.TAPERED;
		SlipAlongRuptureModelEnum slipModelType = SlipAlongRuptureModelEnum.UNIFORM;
		
		// Scaling Relationship
		ScalingRelationshipEnum scalingRel = ScalingRelationshipEnum.WC94_SRL_ALL;
		double relativeSectRateWt=1;
		
		// Segmentation constraint filename and weight (file created by running the following method once):
//		wasatchInversion.writeApriorRupRatesForSegmentationConstrints();
//		System.exit(-1);
		String segmentationConstrFilename = ROOT_DATA_DIR+SEGMENT_BOUNDARY_DATA_FILE;
		double relative_segmentationConstrWt = 0;
		
		// Misc settings
		double relative_aPrioriRupWt = 0;	//
		boolean wtedInversion = true;
//		double minRupRate = 1e-8;
		double minRupRate = 0.0;
		boolean applyProbVisible = true;
		double moRateReduction = 0.1;	// this is the amount to reduce due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)
		double relativeMFD_constraintWt = 0; // 
		
		// Inversion Solution Type:
//		InversionSolutionType solutionType = InversionSolutionType.SIMULATED_ANNEALING;

//		InversionSolutionType solutionType = InversionSolutionType.FROM_FILE;
		// the following is the directory where to find this file - CANNOT COMMENT THIS OUT
		String rupRatesFileDirName = ROOT_PATH+"OutputDataAndFigs_SA10_Segmented_Uniform_4/";	// this is only used if solutionType = InversionSolutionType.FROM_FILE

		InversionSolutionType solutionType = InversionSolutionType.NON_NEGATIVE_LEAST_SQUARES;

		int numSolutions = 10; // this is ignored for NON_NEGATIVE_LEAST_SQUARES which only has one possible solution
		
		// Simulated Annealing Parameters (ignored for NON_NEGATIVE_LEAST_SQUARES)
		CoolingScheduleType saCooling = CoolingScheduleType.VERYFAST_SA;
		long numIterations = (long) 1e6;
		boolean initStateFromAprioriRupRates = false;
		long randomSeed = System.currentTimeMillis();
//		long randomSeed = 1525892588112l; // for reproducibility; note that the last character here is the letter "l" to indicated a long value
		
		// Data to plot and/or save:
		boolean popUpPlots = true;	// this tells whether to show plots in a window (set null to turn off; e.g., for HPC)
		boolean doDataFits=true;
		boolean doMagHistograms=true;
		boolean doNonZeroRateRups=true;
	    boolean doSectPartMFDs=true;
	    
	    // to make hazard curve (set loc=null to ignore)
		Location hazCurveLoc = new Location(40.75, -111.90);	// Salt Lake City
	    String hazCurveLocName = "Salt Lake City Hazard Curve";
	    
	    // to make hazard maps
	    boolean makeHazardMaps = true;
		double saPeriodForHaz = 1.0;	// set as 0.0 for PGA

		// THIS IS THE END OF THE INVERSION SETTINGS

		// Create an instance of this inversion class
		WasatchInversion wasatchInversion = new WasatchInversion(slipRates);
				

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
		
//	    // THE FOLLOWING IS TO MAKE HAZARD MAP RATIOS FOR DATA IN FILES
//	    String fileName1 = ROOT_PATH+"OutputDataAndFigs_SA10_Unsegmented_Uniform"+"/hazardMaps/1.0secSA_2in50.txt";
//	    String fileName2 = ROOT_PATH+"OutputDataAndFigs_SA10_Segmented_Uniform"+"/hazardMaps/1.0secSA_2in50.txt";
//	    String label = "1.0secSA_2in50_2in50_ratio";
//	    wasatchInversion.makeHazardMapRatio(fileName1, fileName2, label, dirName+"/hazardMaps", true);
//		System.exit(-1);
//		


//		System.exit(0);
		
		// this will be used to keep track of runtimes
		long startTimeMillis = System.currentTimeMillis();
		if(D)
			System.out.println("Starting Inversion");
		
		// create an instance of the inversion class with the above settings
		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
				name,
				slipRates.getName(),
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
		
		// make the directory for storing results
	    File file = new File(dirName);
	    file.mkdirs();
	    
	    // write the setup info to a file
	    fltSysRupInversion.writeInversionSetUpInfoToFile(dirName);
	    
	    switch (solutionType) {
	    		case NON_NEGATIVE_LEAST_SQUARES:
	    			if(D) System.out.println("NON_NEGATIVE_LEAST_SQUARES");
	    			fltSysRupInversion.doInversionNNLS();
	    			break;
	    		case SIMULATED_ANNEALING:
	    			if(numSolutions==1) {
		    			if(D) System.out.println("SIMULATED_ANNEALING; numSolutions=1");
	    				fltSysRupInversion.doInversionSA(numIterations, initStateFromAprioriRupRates, randomSeed, saCooling);
	    			}
	    			else if(numSolutions>1) {
		    			if(D) System.out.println("SIMULATED_ANNEALING; numSolutions="+numSolutions);
	    				fltSysRupInversion.doInversionSA_MultTimes(numIterations, initStateFromAprioriRupRates, randomSeed, numSolutions, dirName, saCooling);
	    			}
	    			else {
	    				throw new RuntimeException("bad numIterations");
	    			}
	    			break;
	    		case FROM_FILE:
			    	if(numSolutions==1) {
		    			if(D) System.out.println("FROM_FILE; numSolutions=1");
		    			String rupRatesFileName = rupRatesFileDirName+ "ruptureRatesAlt.txt";
			    		double[] rupRatesArray = wasatchInversion.readRuptureRatesFromFile(rupRatesFileName);
			    		fltSysRupInversion.setSolution(rupRatesArray, "Solution from file: "+rupRatesFileName);
			    	}
			    	else if(numSolutions>1) {
		    			if(D) System.out.println("FROM_FILE; numSolutions="+numSolutions);
		    			ArrayList<double[]> rupRatesArrayList = new ArrayList<double[]>();
		    			for(int i=0;i<numSolutions;i++) {
			    			String rupRatesFileName = rupRatesFileDirName+ "ruptureRates_"+i+".txt";
			    			rupRatesArrayList.add(wasatchInversion.readRuptureRatesFromFile(rupRatesFileName));
			    			fltSysRupInversion.setMultipleSolutions(rupRatesArrayList, "Multiple solutions read from files with prefix "+rupRatesFileName, dirName);
		    			}
			    	}
			    	else {
			    		throw new RuntimeException("bad numIterations");
			    	}
			    	break;
	    }

		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
		if(D) System.out.println("Done with Inversion after "+(float)runTimeSec+" seconds.");
				
		// write results to file
		fltSysRupInversion.writeInversionRunInfoToFile(dirName);
			
		if(doDataFits) fltSysRupInversion.writeAndOrPlotDataFits(dirName, popUpPlots);
		if(doMagHistograms) fltSysRupInversion.writeAndOrPlotMagHistograms(dirName, popUpPlots);
		if(doNonZeroRateRups) fltSysRupInversion.writeAndOrPlotNonZeroRateRups(dirName, popUpPlots);
		if(doSectPartMFDs) fltSysRupInversion.writeAndOrPlotSectPartMFDs(dirName, popUpPlots);
	    
	    // hazard curve:
		if(hazCurveLoc != null) {
				wasatchInversion.writeAndOrPlotHazardCurve(fltSysRupInversion, hazCurveLoc, saPeriodForHaz, dirName, popUpPlots, hazCurveLocName);
		}
	    
	    // second parameter here is SA period; set as 0.0 for PGA:
		if(makeHazardMaps)
			wasatchInversion.makeHazardMaps(fltSysRupInversion, saPeriodForHaz, dirName, popUpPlots);
	}
}
