package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.Declustering.GardnerKnopoffDeclustering;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.UCERF3.erf.ETAS.ETAS_CubeDiscretizationParams;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class MakeFigures {
	
	public static final String outputDirString = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse";

	public static void makeNumHazExceedanceRandIML_Figures(double[] imlArray, ArrayList<ArrayList<XY_DataSet>> listList1, 
			ArrayList<ArrayList<XY_DataSet>> listList2, String subDirName) {
		
		String iml1 = Double.toString(imlArray[0]);
		String iml2 = Double.toString(imlArray[1]);
		
		if(subDirName==null)
			subDirName = "NumHazExceedanceRandIMLFigures";
		String dirName = outputDirString+"/"+subDirName;
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdir();
		
		double tsWidth = 6.5;
		double tsHeight = 3.0;
		
		// plot time series
		// there are 10,000 points for 50-year windows
		
		// 30,000 years is what can be resolved with line width of 1.0 (??)
		double startYear = 10000;
		double duration = 20000;
		int index1 = (int)startYear/50;
		int index2 = (int)(startYear+duration)/50;
		Range xRange = new Range(0,duration);
		Range yRange = new Range(-0.5,25);
//		double lineWidth = 72*5.62/(duration/50.); // this is ~1.0 for 20,000 yr duration
		float lineWidth = 1f;
		PlotCurveCharacterstics pcc = new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, Color.BLACK);

		PlottingUtils.writeAndOrPlotFuncs(getVerticalBarFunction(listList1.get(0).get(0), index1, index2),
				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml1, true);
		PlottingUtils.writeAndOrPlotFuncs(getVerticalBarFunction(listList2.get(0).get(0), index1, index2),
				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml1, true);
		pcc = new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.BLACK);

//		yRange = new Range(-0.05,3);
//		
//		PlottingUtils.writeAndOrPlotFuncs(getVerticalBarFunction(listList1.get(0).get(1), index1, index2),
//				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
//				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml2, true);
//		PlottingUtils.writeAndOrPlotFuncs(getVerticalBarFunction(listList2.get(0).get(1), index1, index2),
//				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
//				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml2, true);
		
		
		
//		PlottingUtils.writeAndOrPlotFuncs(listList1.get(0).get(0),
//				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
//				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml1, true);
//		PlottingUtils.writeAndOrPlotFuncs(listList2.get(0).get(0),
//				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
//				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml1, true);
//		pcc = new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.BLACK);
		xRange = new Range(0,5e5);
		yRange = new Range(-0.05,3);
		
		PlottingUtils.writeAndOrPlotFuncs(listList1.get(0).get(1),
				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml2, true);
		PlottingUtils.writeAndOrPlotFuncs(listList2.get(0).get(1),
				 pcc, "", "Time (yrs)", "Obs Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml2, true);

		

		// plot pdfs
		ArrayList<XY_DataSet> pdfList1 = new ArrayList<XY_DataSet>();
		pdfList1.add(getStairStepFunction(listList2.get(1).get(0))); // Randomized PDF for first IML
		pdfList1.add(getStairStepFunction(listList2.get(2).get(0))); // Poisson expectation PDF for first IML
		pdfList1.add(getStairStepFunction(listList1.get(1).get(0))); // TD PDF for first IML
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		PlottingUtils.writeAndOrPlotFuncs(pdfList1, plotChars, "", "Num Exceedances", "Probability (per 50 yrs)", 
				new Range(-0.6,18), new Range(0.9e-4,1.1), false, true, dirName+"/numPDF_IML"+iml1, true);
		ArrayList<XY_DataSet> pdfList2 = new ArrayList<XY_DataSet>();
		pdfList2.add(getStairStepFunction(listList2.get(1).get(1))); // Randomized PDF for second IML
		pdfList2.add(getStairStepFunction(listList2.get(2).get(1))); // Poisson expectation PDF for second IML
		pdfList2.add(getStairStepFunction(listList1.get(1).get(1))); // TD PDF for second IML
		PlottingUtils.writeAndOrPlotFuncs(pdfList2, plotChars, "", "Num Exceedances", "Probability (per 50 yrs)", 
				new Range(-0.52,4), new Range(0.9e-4,1.1), false, true, dirName+"/numPDF_IML"+iml2, true);

	}
	
	
	
	public static void makeNumHazExceedanceFigures(double[] imlArray, ArrayList<ArrayList<XY_DataSet>> listList1, 
			ArrayList<ArrayList<XY_DataSet>> listList2, String subDirName) {
		
		String iml1 = Double.toString(imlArray[0]);
		String iml2 = Double.toString(imlArray[1]);
		if(subDirName==null)
			subDirName = "NumHazExceedanceFigures";
		String dirName = outputDirString+"/"+subDirName;
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdir();
		
		double tsWidth = 6.5;
		double tsHeight = 1.5;
		
		// plot time series
		// there are 10,000 points for 50-year windows
		
		// 30,000 years is what can be resolved with line width of 1.0 (??)
		double startYear = 200000;
		double duration = 250000;
//		int index1 = (int)startYear/50;
//		int index2 = (int)(startYear+duration)/50;
		Range xRange = new Range(startYear,duration);
		Range yRange = new Range(-4,1.5);
//		double lineWidth = 72*5.62/(duration/50.); // this is ~1.0 for 20,000 yr duration
		float lineWidth = 0.5f;
		PlotCurveCharacterstics pcc = new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, Color.BLACK);

		PlottingUtils.writeAndOrPlotFuncs(convertYvalsToLog10(listList1.get(0).get(0)),
				 pcc, "", "Time (yrs)", "Log10 Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml1, true);
		PlottingUtils.writeAndOrPlotFuncs(convertYvalsToLog10(listList2.get(0).get(0)),
				 pcc, "", "Time (yrs)", "Log10 Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml1, true);
//		pcc = new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.BLACK);

		
//		xRange = new Range(0,5e5);
		yRange = new Range(-8,0.5);
		
		PlottingUtils.writeAndOrPlotFuncs(convertYvalsToLog10(listList1.get(0).get(1)),
				 pcc, "", "Time (yrs)", "Log10 Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesFullTD_IML"+iml2, true);
		PlottingUtils.writeAndOrPlotFuncs(convertYvalsToLog10(listList2.get(0).get(1)),
				 pcc, "", "Time (yrs)", "Log10 Num Exceed per 50 yrs", 
				 xRange, yRange, false, false, tsWidth, tsHeight, dirName+"/timeSeriesPoiss_IML"+iml2, true);

		

		// plot pdfs
		ArrayList<XY_DataSet> pdfList1 = new ArrayList<XY_DataSet>();
		pdfList1.add(getStairStepFunction(listList2.get(1).get(0))); // Randomized PDF for first IML
		pdfList1.add(getStairStepFunction(listList1.get(1).get(0))); // TD PDF for first IML
		pdfList1.add(listList2.get(3).get(0)); // Randomized CDF for first IML
		pdfList1.add(listList1.get(3).get(0)); // TD CDF for first IML
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.RED));
		PlottingUtils.writeAndOrPlotFuncs(pdfList1, plotChars, "", "Num Exceedances", "Probability (per 50 yrs)", 
				new Range(-0.2,15), new Range(1e-4,10), false, true, dirName+"/numPDF_IML"+iml1, true);
		ArrayList<XY_DataSet> pdfList2 = new ArrayList<XY_DataSet>();
		pdfList2.add(getStairStepFunction(listList2.get(1).get(1))); // Randomized PDF for second IML
		pdfList2.add(getStairStepFunction(listList1.get(1).get(1))); // TD PDF for second IML
		pdfList2.add(getStairStepFunction(listList2.get(3).get(1))); // Randomized CDF for second IML
		pdfList2.add(getStairStepFunction(listList1.get(3).get(1))); // TD CDF for second IML
		PlottingUtils.writeAndOrPlotFuncs(pdfList2, plotChars, "", "Num Exceedances", "Probability (per 50 yrs)", 
				new Range(-0.2,3), new Range(1e-4,10), false, true, dirName+"/numPDF_IML"+iml2, true);

	}

	
	/**
	 * This turns a function into a stair-step histogram
	 * @param func
	 * @return
	 */
	static public XY_DataSet getStairStepFunction(XY_DataSet func) {
		// first points
		DefaultXY_DataSet newFunc = new DefaultXY_DataSet(); 
		double firstX = func.getX(0) - 0.5*(func.getX(1)-func.getX(0));
		double secondX = func.getX(0) + 0.5*(func.getX(1)-func.getX(0));
		newFunc.set(firstX,1e-9);
		newFunc.set(firstX,func.getY(0));
		newFunc.set(secondX,func.getY(0));
		for(int i=1;i<func.size()-1; i++) {
			firstX = 0.5*(func.getX(i) + func.getX(i-1));
			secondX = 0.5*(func.getX(i) + func.getX(i+1));
			newFunc.set(firstX,func.getY(i));
			newFunc.set(secondX,func.getY(i));
		}
		// last points
		int lastIndex = func.size()-1;
		firstX = func.getX(lastIndex) - 0.5*(func.getX(lastIndex)-func.getX(lastIndex-1));
		secondX = func.getX(lastIndex) + 0.5*(func.getX(lastIndex)-func.getX(lastIndex-1));
		newFunc.set(firstX,func.getY(lastIndex));
		newFunc.set(secondX,func.getY(lastIndex));
		newFunc.set(secondX,1e-9);

		newFunc.setName(func.getName());
		newFunc.setInfo(func.getInfo());
		return newFunc;
	}
	
	
	/**
	 * This turns a function into a stair-step histogram
	 * @param func
	 * @return
	 */
	static public XY_DataSet getVerticalBarFunction(XY_DataSet func, int firstIndex, int lastIndex) {
		
		DefaultXY_DataSet newFunc = new DefaultXY_DataSet(); 
		double xAtFirstIndex = func.getX(firstIndex);

		for(int i=firstIndex;i<=lastIndex; i++) {
			if(func.getY(i)>0.0) {
				newFunc.set(func.getX(i)-xAtFirstIndex, 0.0);
				newFunc.set(func.getX(i)-xAtFirstIndex, func.getY(i));
				newFunc.set(func.getX(i)-xAtFirstIndex, 0.0);				
			}
			else {
				newFunc.set(func.getX(i)-xAtFirstIndex, func.getY(i));			
			}
		}
		newFunc.setName(func.getName());
		newFunc.setInfo(func.getInfo());
		return newFunc;
	}
	
	public static XY_DataSet convertYvalsToLog10(XY_DataSet func) {
		for(int i=0;i<func.size();i++)
			func.set(i,Math.log10(func.getY(i)));
		return func;
	}
	
	
	public static String plotLongTermEventRateMapFromFile(String dirName) {
		String dataFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/Data/U3ETAS_SimulationData/gridded_nucleation_m2.5.xyz";

		GriddedRegion gr = new GriddedRegion(new CaliforniaRegions.RELM_TESTING(), 0.02, GriddedRegion.ANCHOR_0_0);
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gr, true);

		File file = new File(dataFileName);
		List<String> fileLines;
		
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			for(int i=0; i<fileLines.size();i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split("\t");
				double lon = Double.parseDouble(split[0]);
				double lat = Double.parseDouble(split[1]);
				double rate = Double.parseDouble(split[2]);
				int indexTest = gr.indexForLocation(new Location(lat,lon));
				if(indexTest != i)
					throw new RuntimeException("Problem with index");
				xyzDataSet.set(i, rate);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return plotEventRateMap(xyzDataSet, gr, 2.5, dirName);

	}
	
	
	
	public static String plotTest() {
		
		String dirName = "testMap";
		CaliforniaRegions.RELM_TESTING_GRIDDED griddedRegion = RELM_RegionUtils.getGriddedRegionInstance();
		ETAS_CubeDiscretizationParams  cubeParams = new ETAS_CubeDiscretizationParams(griddedRegion);
		GriddedRegion gridRegForCubes = cubeParams.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);

//		int index = gridRegForCubes.indexForLocation(new Location(37d,-120d));
//		xyzDataSet.set(index, 1);
		
		for(double lat = gridRegForCubes.getMinGridLat(); lat<gridRegForCubes.getMaxGridLat(); lat+=0.06)
			for(double lon = gridRegForCubes.getMinGridLon(); lon<gridRegForCubes.getMaxGridLon(); lon+=0.06) {
				int index = gridRegForCubes.indexForLocation(new Location(lat,lon));
				if(index != -1)
					xyzDataSet.set(index, 1);
			}		
		
		for(int i=0;i<xyzDataSet.size();i++)
			if(xyzDataSet.get(i)==0)
				xyzDataSet.set(i, 1e-10);
		
		return plotEventRateMap(xyzDataSet, gridRegForCubes, 2.5, dirName);
	}

	
	
	/**
	 * This plots the hypocenter rates above the specified magnitude for the given catalog
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public static String plotEventRateMapForCatalog(ObsEqkRupList catalog, double magThresh, 
			String dirName, double duration) {
		
		CaliforniaRegions.RELM_TESTING_GRIDDED griddedRegion = RELM_RegionUtils.getGriddedRegionInstance();
		ETAS_CubeDiscretizationParams  cubeParams = new ETAS_CubeDiscretizationParams(griddedRegion);
		GriddedRegion gridRegForCubes = cubeParams.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);
//		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(griddedRegion, true);

//		int index = gridRegForCubes.indexForLocation(new Location(40d,-125d));
//		System.out.println(gridRegForCubes.getLocation(index));
//		System.exit(0);
		
		System.out.println(catalog.size()+" events for \n\t"+dirName);

		for(ObsEqkRupture rup : catalog) {
			if(rup.getMag()<magThresh)
				continue;
			int locIndex = xyzDataSet.indexOf(rup.getHypocenterLocation());
			if(locIndex == -1)
				continue;
			double lastValue = xyzDataSet.get(locIndex);
			xyzDataSet.set(locIndex, lastValue+1.0/duration);
//System.out.println(rup.getHypocenterLocation());
		}
//System.exit(0);
		for(int i=0;i<xyzDataSet.size();i++)
			if(xyzDataSet.get(i)==0)
				xyzDataSet.set(i, 1e-10);
		
		return plotEventRateMap(xyzDataSet, gridRegForCubes, magThresh, dirName);

	}
	
	
	
	/**
	 * This plots an event rate maps
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public static String plotEventRateMap(GriddedGeoDataSet xyzDataSet, GriddedRegion gridRegForCubes,
			double magThresh, String dirName) {
		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same

		// this draws each cell as a box to avoid weird aliasing issues
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME, true);

		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);

//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-2d);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1d);			
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3.5d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,0.5d);	
		cptParam.getValue().setNanColor(Color.WHITE);
		
		System.out.println("\tminZ="+(float)xyzDataSet.getMinZ()+"\n\tmaxZ="+(float)xyzDataSet.getMaxZ()+
				"\n\tnumGridPts="+xyzDataSet.size());
		
		// these have zero influence
//		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
//		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, 0.1);
//		mapGen.setParameter(GMT_MapGenerator.DPI_PARAM_NAME, 300);

		String metadata = "Map from calling plotEventRateMapForCatalog(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "M>="+magThresh+" Rates", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(GMT_CA_Maps.GMT_DIR, dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "For rates map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}

	
	
	
	public static void makeFigure1_Parts() {
		//
		double startTimeYears = 100; // this is the time within the catalog to start, will need to set lower than this as our simulations are 500 years long
		double durationYears = 100;
		int plotWidthPixels = 2000;
		double binWidthDays = 30;
		String prefix = null;
		double annotateMinMag = 6;
		double annotateBinWidth = 0.5;
		// use this because it has <≥2.5 events
		String fullCatFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/testDebugRun/2020_12_03-Start2012_500yr_kCOV1p5_ScaleFactor1p4_Spontaneous_HistoricalCatalog/results_complete.bin"; 
		ObsEqkRupList rupList=null;
		try {
			rupList = U3ETAS_SimulationAnalysis.loadCatalogs(new File(U3ETAS_SimulationAnalysis.fssFileName), new File(fullCatFileName), 2.5).get(0);
		} catch (Exception e) {
			e.printStackTrace();
		}
		ArrayList<ETAS_EqkRupture> catalog = new ArrayList<ETAS_EqkRupture>();
		for(ObsEqkRupture rup:rupList)
			catalog.add((ETAS_EqkRupture)rup);
		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure1";
		File outputDir = new File(dirName);
		if(!outputDir.exists()) outputDir.mkdir();
		try {
			ETAS_SimAnalysisTools.plotRateOverTime(catalog, startTimeYears, durationYears, binWidthDays, outputDir, prefix, plotWidthPixels,
					annotateMinMag, annotateBinWidth);
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		long actualOT = catalog.get(0).getOriginTime();
		long actualMaxOT = catalog.get(catalog.size()-1).getOriginTime();
		long startOT;
		if (startTimeYears > 0)
			startOT = actualOT + (long)(startTimeYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		else
			startOT = actualOT;
		Preconditions.checkState(startOT < actualMaxOT, "Start time is after end of catalog");
		long maxOT = (long)(startOT + durationYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		
		ObsEqkRupList rupListCut = new ObsEqkRupList();
		for(ObsEqkRupture rup:rupList) {
			long ot = rup.getOriginTime();
			if(ot>startOT && ot<maxOT)
				rupListCut.add(rup);
		}
		
		ObsEqkRupList declusteredCatalog = GardnerKnopoffDeclustering.getDeclusteredCatalog(rupList);
		ObsEqkRupList rupListDeclusteredCut = new ObsEqkRupList();
		for(ObsEqkRupture rup:declusteredCatalog) {
			long ot = rup.getOriginTime();
			if(ot>startOT && ot<maxOT)
				rupListDeclusteredCut.add(rup);
		}



		MakeFigures.plotEventRateMapForCatalog(rupListCut, 2.5, "Figure1b_Data", durationYears);

		MakeFigures.plotEventRateMapForCatalog(rupListDeclusteredCut, 2.5, "DeclusteredMap_Data", durationYears);

		plotLongTermEventRateMapFromFile("Figure1c_Data");
	}
	
	
	
	
	public static void makeFigure2_Parts(ArrayList<XY_DataSet> funcList) {
		ArrayList<XY_DataSet> funcList2 = new ArrayList<XY_DataSet>();
		funcList2.add(funcList.get(3)); // full cum MFD
		funcList2.add(funcList.get(4)); // GK declustered cum MFD
		funcList2.add(funcList.get(5)); // Spontaneous events only cum MFD
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.GRAY));
		Range xRange = new Range(5, 8.5);
		Range yRange = new Range(4e-5, 10.);
		String dirName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/FiguresFromEclipse/Figure2";
		File outputDir = new File(dirName);
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/Figure2a";
		PlottingUtils.writeAndOrPlotFuncs(funcList2, plotChars, null, "Magnitude", "Cumulative Rate (per year)", xRange, yRange, 
				false, true, 3.5, 3.0, fileNamePrefix, true);
		
		ArrayList<XY_DataSet> funcList3 = new  ArrayList<XY_DataSet>();
		funcList3.add(funcList.get(6)); // 
		funcList3.add(funcList.get(7)); // 
		funcList3.add(funcList.get(8)); // 
		fileNamePrefix = dirName+"/Figure2b";
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.BLUE));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.GRAY));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.BLACK));
		Range yRange2 = new Range(0, 1.1);
		PlottingUtils.writeAndOrPlotFuncs(funcList3, plotChars2, null, "Magnitude", "Fraction Main Shocks", xRange, yRange2, 
				false, false, 3.5, 3.0, fileNamePrefix, true);

	}

	
	public static String getMaxRatioInfo(DiscretizedFunc func1, DiscretizedFunc func2) {
		
		double maxProbRatio = Double.NEGATIVE_INFINITY;
		double maxIMLratio = Double.NEGATIVE_INFINITY;
		int indexProb=-1;
		int indexIML=-1;
		double imlAtMaxIML_Ratio=-1;
		for(int i=0;i<func1.size();i++) {
			if(func1.getY(i)<1e-3)
				break;
			double probRatio = func1.getY(i)/func2.getY(i);
			if(maxProbRatio<probRatio) {
				maxProbRatio=probRatio;
				indexProb = i;
			}
			double iml1 = func1.getX(i);
			double prob = func1.getY(i);
			double iml2 = func2.getFirstInterpolatedX_inLogXLogYDomain(prob);
			double imlRatio = iml1/iml2;
			if(maxIMLratio<imlRatio) {
				maxIMLratio=imlRatio;
				indexIML=i;
				imlAtMaxIML_Ratio=(iml1+iml2)/2;;
			}
			System.out.println((float)((iml1+iml2)/2.)+"\t"+(float)imlRatio);
			
		}
		
		double imlAtMaxProbRatio = func1.getX(indexProb);
		String str1 = "Max Prob Ratio is "+(float)maxProbRatio+" at IML = "+imlAtMaxProbRatio;

		String str2 = "Max IML Ratio is "+(float)maxIMLratio+" at ave IML = "+imlAtMaxIML_Ratio;

		return str1+"\n"+str2;
	}

	
	public static void makeFigure3_Parts(ArrayList<UncertainArbDiscDataset[]> dataSetsArray, ArrayList<XY_DataSet> funcsArray, double duration, 
			double saPeriod, String dirName, boolean popupWindow, String plotTitle) {
		
		dataSetsArray.get(0)[0].getLower().setName("Full TD min");
		dataSetsArray.get(0)[0].getUpper().setName("Full TD max");
		dataSetsArray.get(0)[1].setName("Full TD; mean & 95% conf");
		dataSetsArray.get(1)[0].getLower().setName("GK Declustered min");
		dataSetsArray.get(1)[0].getUpper().setName("GK Declustered max");
		dataSetsArray.get(1)[1].setName("GK Declustered; mean & 95% conf");
		dataSetsArray.get(2)[0].getLower().setName("Full Randomized min");
		dataSetsArray.get(2)[0].getUpper().setName("Full Randomized max");
		dataSetsArray.get(2)[1].setName("Full Randomized; mean & 95% conf");
		
		System.out.println("MAX RATIO INFO FOR RANDOMIZED VS FULL TD:\n"+getMaxRatioInfo(dataSetsArray.get(2)[0], dataSetsArray.get(0)[0]));
		String newInfo = dataSetsArray.get(2)[1].getInfo()+"\n\nMAX RATIO INFO FOR RANDOMIZED VS FULL TD:\n"+
		getMaxRatioInfo(dataSetsArray.get(2)[0], dataSetsArray.get(0)[0]);
		dataSetsArray.get(2)[1].setInfo(newInfo);
		
		
		Color[] colorArray = {Color.RED, Color.BLUE, Color.BLACK, Color.MAGENTA, Color.CYAN};
		float lineWidth = 0.5f;
//		int colorIndex = 0;
		
		String imtString = U3ETAS_SimulationAnalysis.getIMT_String(saPeriod);

		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	

//		for(UncertainArbDiscDataset[] dataSets:dataSetsArray) {
//			plottingFuncsArray.add(dataSets[1]); // for solid line
//			plottingFuncsArray.add(dataSets[1]); // for shaded region
//			plottingFuncsArray.add(dataSets[0].getLower());
//			plottingFuncsArray.add(dataSets[0].getUpper());
//			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, lineWidth, colorArray[colorIndex]));
//			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, colorArray[colorIndex]));
//			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth, colorArray[colorIndex]));
//			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth, colorArray[colorIndex]));
//			colorIndex += 1;
//		}
		
		plottingFuncsArray.add(funcsArray.get(1));	// Full TD Random IML
		plottingFuncsArray.add(funcsArray.get(0));	// Poisson calc
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1.2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1.2f, Color.BLACK));
		
		for(int i = dataSetsArray.size()-1; i>=0;i--) {  // reverse order for better layering
			UncertainArbDiscDataset[] dataSets = dataSetsArray.get(i);
			plottingFuncsArray.add(dataSets[1]); // for solid line
			plottingFuncsArray.add(dataSets[1]); // for shaded region
			plottingFuncsArray.add(dataSets[0].getLower());
			plottingFuncsArray.add(dataSets[0].getUpper());
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, lineWidth, colorArray[i]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, colorArray[i]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth*2, colorArray[i]));
			plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth*2, colorArray[i]));
		}


		
		File outputDir = new File(dirName);
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/Figure3a";
		String xAxisLabel = imtString;
		String yAxisLabel = ((int)Math.round(duration))+" yr Exceedance Probability";
		boolean logX = true;
		boolean logY = true;
		Range xAxisRange = new Range(3e-2,10);
		Range yAxisRange = new Range(1e-4,1.1);
		
		plotChars.get(0).setSymbolWidth(1.2f);
		plotChars.get(1).setSymbolWidth(1.2f);
		
		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);

		
		fileNamePrefix = dirName+"/Figure3b";
		xAxisRange = new Range(1e-1,3);
		yAxisRange = new Range(1e-2,1.1);
		

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
		
		ArrayList<XY_DataSet> plottingFuncsArray2 = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();	
		plottingFuncsArray2.add(funcsArray.get(1));			// Full TD Random IML
		plottingFuncsArray2.add(funcsArray.get(1));			// Full TD Random IML
		plottingFuncsArray2.add(dataSetsArray.get(0)[1]);	// Full TD

		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, lineWidth, new Color(0f,102f/255f,0f)));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, new Color(0f,102f/255f,0f)));
		plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, Color.red));
		fileNamePrefix = dirName+"/Figure3c";
		xAxisRange = new Range(3e-2,10);
		yAxisRange = new Range(1e-4,1.1);
		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray2, plotChars2, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);



	}
	
	public static void makeFigure9_Parts(ArrayList<UncertainArbDiscDataset[]> dataSetsArray, ArrayList<XY_DataSet> funcsArray, double duration, 
			double saPeriod, String dirName, boolean popupWindow, String plotTitle) {
		
		dataSetsArray.get(0)[0].getLower().setName("Full TD min");
		dataSetsArray.get(0)[0].getUpper().setName("Full TD max");
		dataSetsArray.get(0)[1].setName("Full TD; mean & 95% conf");
			
		Color[] colorArray = {Color.RED, Color.BLUE, Color.BLACK, Color.MAGENTA, Color.CYAN};
		float lineWidth = 0.5f;
		
		String imtString = U3ETAS_SimulationAnalysis.getIMT_String(saPeriod);

		ArrayList<XY_DataSet> plottingFuncsArray = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	

		plottingFuncsArray.add(funcsArray.get(0));	// Poisson calc
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, Color.BLACK));
		
		UncertainArbDiscDataset[] dataSets = dataSetsArray.get(0);
		plottingFuncsArray.add(dataSets[1]); // for solid line
		plottingFuncsArray.add(dataSets[1]); // for shaded region
		plottingFuncsArray.add(dataSets[0].getLower());
		plottingFuncsArray.add(dataSets[0].getUpper());
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, lineWidth, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, lineWidth, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth*2, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, lineWidth*2, Color.RED));
		
		File outputDir = new File(dirName);
		if(!outputDir.exists()) outputDir.mkdir();
		String fileNamePrefix = dirName+"/Figure9";
		String xAxisLabel = imtString;
		String yAxisLabel = ((int)Math.round(duration))+" yr Exceedance Probability";
		boolean logX = true;
		boolean logY = true;
		Range xAxisRange = new Range(1e-3,10);
		Range yAxisRange = new Range(1e-4,1.1);
		
		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, null, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, fileNamePrefix, popupWindow);
	}


	
	public static void main(String[] args) {
		plotTest();
//		makeFigure1_Parts();
//		plotLongTermEventRateMapFromFile("Figure1c_Data");
	}

}
