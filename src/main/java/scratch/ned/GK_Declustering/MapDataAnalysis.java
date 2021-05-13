package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.math.Quantiles;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class MapDataAnalysis {
	
	final static boolean D = true;
	
	File file = new File("/Users/field/Field_Other/CEA_WGCEP/UCERF3/DeclusteringAnalysis/Data/KevinsCurvesForMapsRuns/results.zip");
	
	GriddedRegion region;
	int numLocs;
	
	double duration = 50;
	
	ArbitrarilyDiscretizedFunc[] pga_full_mean_curves, pga_full_min_curves, pga_full_max_curves, pga_full_low95_curves, pga_full_upp95_curves;
	ArbitrarilyDiscretizedFunc[] pga_poisson_curves;
	ArbitrarilyDiscretizedFunc[] pga_declustered_mean_curves, pga_declustered_min_curves, pga_declustered_max_curves, pga_declustered_low95_curves, pga_declustered_upp95_curves;
	ArbitrarilyDiscretizedFunc[] pga_randomized_mean_curves, pga_randomized_min_curves, pga_randomized_max_curves, pga_randomized_low95_curves, pga_randomized_upp95_curves;

	ArbitrarilyDiscretizedFunc[] sa0pt2_full_mean_curves, sa0pt2_full_min_curves, sa0pt2_full_max_curves, sa0pt2_full_low95_curves, sa0pt2_full_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa0pt2_poisson_curves;
	ArbitrarilyDiscretizedFunc[] sa0pt2_declustered_mean_curves, sa0pt2_declustered_min_curves, sa0pt2_declustered_max_curves, sa0pt2_declustered_low95_curves, sa0pt2_declustered_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa0pt2_randomized_mean_curves, sa0pt2_randomized_min_curves, sa0pt2_randomized_max_curves, sa0pt2_randomized_low95_curves, sa0pt2_randomized_upp95_curves;
	
	ArbitrarilyDiscretizedFunc[] sa1pt0_full_mean_curves, sa1pt0_full_min_curves, sa1pt0_full_max_curves, sa1pt0_full_low95_curves, sa1pt0_full_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa1pt0_poisson_curves;
	ArbitrarilyDiscretizedFunc[] sa1pt0_declustered_mean_curves, sa1pt0_declustered_min_curves, sa1pt0_declustered_max_curves, sa1pt0_declustered_low95_curves, sa1pt0_declustered_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa1pt0_randomized_mean_curves, sa1pt0_randomized_min_curves, sa1pt0_randomized_max_curves, sa1pt0_randomized_low95_curves, sa1pt0_randomized_upp95_curves;

	ArbitrarilyDiscretizedFunc[] sa5pt0_full_mean_curves, sa5pt0_full_min_curves, sa5pt0_full_max_curves, sa5pt0_full_low95_curves, sa5pt0_full_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa5pt0_poisson_curves;
	ArbitrarilyDiscretizedFunc[] sa5pt0_declustered_mean_curves, sa5pt0_declustered_min_curves, sa5pt0_declustered_max_curves, sa5pt0_declustered_low95_curves, sa5pt0_declustered_upp95_curves;
	ArbitrarilyDiscretizedFunc[] sa5pt0_randomized_mean_curves, sa5pt0_randomized_min_curves, sa5pt0_randomized_max_curves, sa5pt0_randomized_low95_curves, sa5pt0_randomized_upp95_curves;

	
	public MapDataAnalysis() {
		
		region = new CaliforniaRegions.RELM_TESTING_GRIDDED(0.1);
		numLocs = region.getNodeCount();
		
		
		System.out.println("Loading data");
		long startTime = System.currentTimeMillis();
		instantiateCurveArrays();
		try {
			loadZippedCSVs(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
		long runtime = System.currentTimeMillis() - startTime;
		System.out.println("That took "+runtime/1000+" seconds");
		
		// test for LA location (have to use the exact loc of 2in50 values to match)
		Location loc = new Location(34.05,-118.25);
		int indexForLA = region.indexForLocation(loc);
		double twoIn50 = this.getIML_ForProb(pga_full_mean_curves[indexForLA], 0.02);
		System.out.println("LA PGA 2in50:  "+twoIn50);
		System.out.println(pga_full_mean_curves[indexForLA]);
		
		Location locTest = region.getLocation(indexForLA);
		System.out.println("\n"+locTest+"\n"+loc);

	}
	
	
	private void instantiateCurveArrays() {
		pga_full_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_full_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_full_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_full_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_full_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_poisson_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_declustered_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_declustered_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_declustered_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_declustered_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_declustered_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_randomized_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_randomized_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_randomized_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_randomized_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		pga_randomized_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];

		sa0pt2_full_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_full_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_full_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_full_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_full_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_poisson_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_declustered_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_declustered_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_declustered_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_declustered_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_declustered_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_randomized_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_randomized_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_randomized_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_randomized_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa0pt2_randomized_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];

			
		sa1pt0_full_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_full_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_full_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_full_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_full_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_poisson_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_declustered_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_declustered_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_declustered_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_declustered_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_declustered_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_randomized_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_randomized_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_randomized_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_randomized_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa1pt0_randomized_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];

		sa5pt0_full_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_full_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_full_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_full_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_full_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_poisson_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_declustered_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_declustered_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_declustered_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_declustered_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_declustered_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_randomized_mean_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_randomized_min_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_randomized_max_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_randomized_low95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
		sa5pt0_randomized_upp95_curves = new ArbitrarilyDiscretizedFunc[numLocs];
	}

	private GriddedGeoDataSet getGriddedGeoDataSet() {
		return new GriddedGeoDataSet(region,false);
	}
	
	private void loadZippedCSVs(File file) throws IOException {
		ZipFile zip = new ZipFile(file);
		
		Enumeration<? extends ZipEntry> entries = zip.entries();
		
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			String name = entry.getName();
			if (!name.endsWith(".csv"))
				continue;
//			System.out.println("Loading "+entry);
			String dirName = name.substring(0, name.indexOf('/'));
			String fileName = name.substring(name.indexOf('/')+1);
			
			// parse the name for node/location
			String[] split = dirName.split("_");
			int node = Integer.parseInt(split[1]);
			double lat = Double.parseDouble(split[2]);
			double lon = Double.parseDouble(split[3]);
			Location loc = new Location(lat, lon);
			
			String imt;
			double period;
			if (fileName.startsWith("sa_")) {
				imt = SA_Param.NAME;
				String periodStr = fileName.substring(3, fileName.indexOf(".csv"));
				period = Double.parseDouble(periodStr);
			} else {
				imt = PGA_Param.NAME;
				period = 0d;
			}
			
//			System.out.println(fileName);
//			System.out.println("Loading "+name+": node="+node+", loc="+loc+", "+imt+" ("+(float)period+" s)");
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			
			int col = 1;
			
			if(fileName.equals("pga.csv")) {
//				if(node != pga_full_mean_curves.size())
//					throw new RuntimeException("node="+node+" and size = "+pga_full_mean_curves.size());
//				System.out.println("node="+node);
				pga_full_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_full_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_full_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_full_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_full_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_poisson_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_declustered_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_declustered_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_declustered_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_declustered_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_declustered_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_randomized_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_randomized_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_randomized_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_randomized_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				pga_randomized_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
			}
			else if(fileName.equals("sa_0.2.csv")) {
				sa0pt2_full_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_full_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_full_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_full_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_full_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_poisson_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_declustered_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_declustered_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_declustered_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_declustered_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_declustered_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_randomized_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_randomized_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_randomized_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_randomized_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa0pt2_randomized_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
			}
			else if(fileName.equals("sa_1.0.csv")) {
				sa1pt0_full_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_full_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_full_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_full_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_full_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_poisson_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_declustered_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_declustered_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_declustered_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_declustered_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_declustered_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_randomized_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_randomized_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_randomized_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_randomized_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa1pt0_randomized_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);				
			}
			else if(fileName.equals("sa_5.0.csv")) {
				sa5pt0_full_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_full_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_full_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_full_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_full_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_poisson_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_declustered_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_declustered_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_declustered_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_declustered_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_declustered_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_randomized_mean_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_randomized_min_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_randomized_max_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_randomized_low95_curves[node] = loadFuncFromCSV(csv, 0, col++);
				sa5pt0_randomized_upp95_curves[node] = loadFuncFromCSV(csv, 0, col++);				
			}
			else {
				throw new RuntimeException("Problem");

			}
			
//			ArbitrarilyDiscretizedFunc fullMean = loadFuncFromCSV(csv, 0, col++);
//			UncertainArbDiscDataset fullMinMax = new UncertainArbDiscDataset(
//					fullMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			UncertainArbDiscDataset fullConf95 = new UncertainArbDiscDataset(
//					fullMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			ArbitrarilyDiscretizedFunc poisson = loadFuncFromCSV(csv, 0, col++);
//			ArbitrarilyDiscretizedFunc declusteredMean = loadFuncFromCSV(csv, 0, col++);
//			UncertainArbDiscDataset declusteredMinMax = new UncertainArbDiscDataset(
//					declusteredMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			UncertainArbDiscDataset declusteredConf95 = new UncertainArbDiscDataset(
//					declusteredMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			ArbitrarilyDiscretizedFunc randomMean = loadFuncFromCSV(csv, 0, col++);
//			UncertainArbDiscDataset randomMinMax = new UncertainArbDiscDataset(
//					randomMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			UncertainArbDiscDataset randomConf95 = new UncertainArbDiscDataset(
//					randomMean, loadFuncFromCSV(csv, 0, col++), loadFuncFromCSV(csv, 0, col++));
//			if (node == 0)
//				System.out.println(fullMinMax);
			Preconditions.checkState(col == csv.getNumCols());
		}
		
		zip.close();
		
		System.out.println("pga_full_mean_curves.length = "+pga_full_mean_curves.length);
		System.out.println("getGriddedGeoDataSet().size() = "+getGriddedGeoDataSet().size());
	}
	
	private static ArbitrarilyDiscretizedFunc loadFuncFromCSV(CSVFile<String> csv, int xCol, int yCol) {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		for (int row=1; row<csv.getNumRows(); row++)
			func.set(csv.getDouble(row, xCol), csv.getDouble(row, yCol));
		return func;
	}
	
	
	public double getIML_ForProb(ArbitrarilyDiscretizedFunc func, double prob) {
		return func.getFirstInterpolatedX_inLogXLogYDomain(prob);
	}
	
	public GriddedGeoDataSet getIML_ForProbGeoDataset(ArbitrarilyDiscretizedFunc[] funcArray, double prob) {
		GriddedGeoDataSet dataset = getGriddedGeoDataSet();
		for(int i=0;i<funcArray.length;i++)
			dataset.set(i, funcArray[i].getFirstInterpolatedX_inLogXLogYDomain(prob));
		return dataset;
	}
	

	public void quickPlot(ArrayList<XY_DataSet> plottingFuncsArray, String xAxisLabel, String yAxisLabel, String title) {
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		Range xAxisRange = null;
		Range yAxisRange = null;
		boolean logX = false;
		boolean logY = false;

		PlottingUtils.writeAndOrPlotFuncs(plottingFuncsArray, plotChars, title, xAxisLabel, yAxisLabel, 
				xAxisRange, yAxisRange, logX, logY, null, true);

	}

	private String compareCases(ArbitrarilyDiscretizedFunc[] case1_FuncsArray, 
			ArbitrarilyDiscretizedFunc[] case2_FuncsArray, double prob, String title, String imtString) {
		GriddedGeoDataSet dataset1 = getIML_ForProbGeoDataset(case1_FuncsArray, prob);
		GriddedGeoDataSet dataset2 = getIML_ForProbGeoDataset(case2_FuncsArray, prob);
		return compareCases(dataset1, dataset2, title, imtString);
	}

	
	
	private String compareCases(GriddedGeoDataSet dataset1, GriddedGeoDataSet dataset2, String title, String imtString) {
		
		String dirName = MakeFigures.outputDirString+"/MapRatioHistograms";
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdir();
		
		String fileNamePrefix = dirName+"/"+title.replace(" ", "_");


		
		HistogramFunction ratioHist = new HistogramFunction(0.0,2.0,201);
		double mean=0, min=Double.MAX_VALUE, max = 0d; 
		int minLocIndex = -1, maxLocIndex = -1;
		double[] ratioArray = new double[dataset1.size()];
		for(int i=0;i<dataset1.size();i++) {
			double ratio = dataset1.get(i)/dataset2.get(i);
			ratioArray[i] = ratio;
			if(ratio < ratioHist.getMaxX()+ratioHist.getDelta()/2.0)
				ratioHist.add(ratio, 1.0);
			mean += ratio;
			if(min>ratio) {
				min = ratio;
				minLocIndex=i;
			}
			if(max<ratio) {
				max = ratio;
				maxLocIndex=i;
			}

		}
		mean /= dataset1.size();
		Location minLoc = region.getLocation(minLocIndex);
		Location maxLoc = region.getLocation(maxLocIndex);
		
		double median = Quantiles.median().compute(ratioArray);
		
		ratioHist.scale(1.0/(dataset1.size()*ratioHist.getDelta()));
		double maxYval = ratioHist.getMaxY();
		
		String infoString =
				"mean = "+(float)mean+"   ("+((int)Math.round((mean-1.0)*100.))+ " %)\n" +
				"median = "+(float)median+"   ("+((int)Math.round((median-1.0)*100.))+ " %)\n" +
				"min = "+(float)min+"   ("+((int)Math.round((min-1.0)*100.))+ " %)\n" +
				"max = "+(float)max+"   ("+((int)Math.round((max-1.0)*100.))+ " %)\n"+
				"min loc: i="+minLocIndex+" ("+(float)minLoc.getLatitude()+", "+(float)minLoc.getLongitude()+")\n"+
				"max loc: i="+maxLocIndex+" ("+(float)maxLoc.getLatitude()+", "+(float)maxLoc.getLongitude()+")\n";
		
		String stringToPassBack = title+"\t"+(float)mean+"\t"+(float)min+"\t"+(float)max+"\t"+(float)maxYval+"\n";
//		System.out.println(title+"\t"+(float)mean+"\t"+(float)min+"\t"+(float)max);
		
		float roundedMean = (float)(Math.round(mean*100d)/100d);
		float roundedMin = (float)(Math.round(min*100d)/100d);
		float roundedMax = (float)(Math.round(max*100d)/100d);
		String annotation = imtString+" "+roundedMean+" ("+roundedMin+", "+roundedMax+")";
		
		ratioHist.setName(title);
		ratioHist.setInfo(infoString);
		
		PlotCurveCharacterstics plotChar =new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK);
		
		Range xAxisRange = new Range(0.7,1.3);
		Range yAxisRange = new Range(0.0,60.0);

//		PlottingUtils.writeAndOrPlotFuncs(ratioHist, plotChar, title, "Ratio", "Density", 
//				xAxisRange, yAxisRange, false, false, 2.2, 2.2, fileNamePrefix, true);
		writeAndOrPlotFuncsWithAnnotation(ratioHist, plotChar, null, null, null, 
				xAxisRange, yAxisRange, false, false, 2.0, 1.75, fileNamePrefix, true, annotation);
		
//		// Make ratio map
//		try {
//			plotAndWriteRatioMap(dataset1, dataset2, title+" Ratio", "test", "testDir");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		
//		System.out.println(ratioHist);
		
		return stringToPassBack;
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
	 * @param widthInches
	 * @param heightInches
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 * @param annotationString - annotation string
	 */
	public static void writeAndOrPlotFuncsWithAnnotation(
			XY_DataSet func, 
			PlotCurveCharacterstics plotChar, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			double widthInches,
			double heightInches,
			String fileNamePrefix, 
			boolean popupWindow,
			String annotationString) {
		
		
		ArrayList<XY_DataSet> funcList = new ArrayList<XY_DataSet>();
		funcList.add(func);
		ArrayList<PlotCurveCharacterstics> plotCharList = new ArrayList<PlotCurveCharacterstics>();
		plotCharList.add(plotChar);

		if(popupWindow) {
			
			GraphWindow graph = new GraphWindow(funcList, plotName);

			if(xAxisRange != null)
				graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
			if(yAxisRange != null)
				graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
			graph.setXLog(logX);
			graph.setYLog(logY);
			graph.setPlotChars(plotCharList);
			graph.setX_AxisLabel(xAxisLabel);
			graph.setY_AxisLabel(yAxisLabel);
			graph.setTickLabelFontSize(18);
			graph.setAxisLabelFontSize(20);

		}
		
		PlotSpec spec = new PlotSpec(funcList, plotCharList, plotName, xAxisLabel, yAxisLabel);
		
		if (fileNamePrefix != null){
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xAxisRange, yAxisRange);
			gp.setTickLabelFontSize(9);
			gp.setAxisLabelFontSize(11);
			gp.setPlotLabelFontSize(9);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec, logX, logY); // spec can be a list
			int width = (int)(widthInches*72.);
			int height = (int)(heightInches*72.);
			gp.getChartPanel().setSize(width, height); 

//			double annotationX = xAxisRange.getLowerBound()+0.25*(xAxisRange.getUpperBound()-xAxisRange.getLowerBound());
//			double annotationY = yAxisRange.getLowerBound()+0.75*(yAxisRange.getUpperBound()-yAxisRange.getLowerBound());
			XYTextAnnotation annotation = new XYTextAnnotation(annotationString,xAxisRange.getCentralValue(),yAxisRange.getCentralValue());
			Font font = new Font("Helvetica", Font.PLAIN, 8);
			annotation.setFont(font);
			gp.getChartPanel().getChart().getXYPlot().addAnnotation(annotation);	
			try {
				gp.saveAsPNG(fileNamePrefix+".png");
				gp.saveAsPDF(fileNamePrefix+".pdf");
				gp.saveAsTXT(fileNamePrefix+".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	

	public static void plotAndWriteRatioMap(GeoDataSet geoDataSet1, GeoDataSet geoDataSet2, String scaleLabel,
			String metadata, String dirName) throws IOException {
		
		GeoDataSet ratioGeoDataSet = GeoDataSetMath.divide(geoDataSet1, geoDataSet2);
		String cptFileString = GMT_CPT_Files.GMT_POLAR.getFileName();
		
//		if(D) {
//			System.out.println("MIN RATIO = "+ratioGeoDataSet.getMinZ());
//			System.out.println("MAX RATIO = "+ratioGeoDataSet.getMaxZ());			
//		}

		GMT_MapGenerator gmt_MapGenerator = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		//override default scale
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME, 0.5);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME, 1.5);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.LOG_PLOT_NAME, false);

		
		
		// must set this parameter this way because the setValue(CPT) method takes a CPT object, and it must be the
		// exact same object as in the constraint (same instance); the setValue(String) method was added for convenience
		// but it won't succeed for the isAllowed(value) call.
		CPTParameter cptParam = (CPTParameter )gmt_MapGenerator.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(cptFileString);
		cptParam.getValue().setNanColor(Color.WHITE);

		GMT_CA_Maps.makeMap(ratioGeoDataSet, scaleLabel, metadata, dirName, gmt_MapGenerator);	
		
	}
	
	
	public static void plotAndWriteHapMap(GeoDataSet geoDataSet, String scaleLabel,
			String metadata, String dirName) throws IOException {
		
		String cptFileString = GMT_CPT_Files.NSHMP_1hz.getFileName();
		
//		if(D) {
//			System.out.println("MIN RATIO = "+ratioGeoDataSet.getMinZ());
//			System.out.println("MAX RATIO = "+ratioGeoDataSet.getMaxZ());			
//		}

		GMT_MapGenerator gmt_MapGenerator = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		//override default scale
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME, 0.0);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME, 1.0);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.LOG_PLOT_NAME, false);

		
		
		// must set this parameter this way because the setValue(CPT) method takes a CPT object, and it must be the
		// exact same object as in the constraint (same instance); the setValue(String) method was added for convenience
		// but it won't succeed for the isAllowed(value) call.
		CPTParameter cptParam = (CPTParameter )gmt_MapGenerator.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(cptFileString);

		GMT_CA_Maps.makeMap(geoDataSet, scaleLabel, metadata, dirName, gmt_MapGenerator);	
		
	}

	
	public GriddedGeoDataSet getRTGM_Dataset(ArbitrarilyDiscretizedFunc[] curvesArray, double saPeriod) {
		
		GriddedGeoDataSet rtgmGriddedDataSetArray = getGriddedGeoDataSet();
		
		for(int i=0;i<curvesArray.length;i++) {
			
			ArbitrarilyDiscretizedFunc hazCurveLinearXvalues = curvesArray[i];
			ArbitrarilyDiscretizedFunc rateCurveLinearXvalues = hazCurveLinearXvalues.deepClone();
			
			// get annual rate curve with linear x-axis values
			boolean hasInfValue=false;
			for(int j=0;j<hazCurveLinearXvalues.size();j++) {
				// convert probability to annual rate
				double rate = -Math.log(1.0-hazCurveLinearXvalues.getY(j))/duration;
				rateCurveLinearXvalues.set(j, rate);
				if(Double.isInfinite(rate)) {
					hasInfValue = true;
//					System.out.println("Infinite Rate from: "+curveLogXvalues.getY(j));
				}
			}
			// remove Inf values that can occur if above curveLogXvalues.getY(j) = 1.0;
			ArbitrarilyDiscretizedFunc rateCurveLinearXvaluesCleaned;
			if(hasInfValue) {
				rateCurveLinearXvaluesCleaned = new ArbitrarilyDiscretizedFunc();
				for(int j=0;j<rateCurveLinearXvalues.size();j++)
					if(!Double.isInfinite(rateCurveLinearXvalues.getY(j)))
						rateCurveLinearXvaluesCleaned.set(rateCurveLinearXvalues.getX(j),rateCurveLinearXvalues.getY(j));
			}
			else
				rateCurveLinearXvaluesCleaned = rateCurveLinearXvalues;

			// now get RTGM
			Frequency freq;
			if(saPeriod == 1.0)
				freq = Frequency.SA_1P00;
			else if (saPeriod == 0.2)
				freq = Frequency.SA_0P20;
			else
				throw new RuntimeException("saPeriod not supported");
			
			
			RTGM rtgm;
			
			// this is to write out test value to compare with the USGS on-line calculator
			try {
//				rtgm = RTGM.create(curveLinearXvaluesCleaned, freq, 0.8).call();
				rtgm = RTGM.create(rateCurveLinearXvaluesCleaned, null, null).call();
				rtgmGriddedDataSetArray.set(i, rtgm.get());
				// this is to write one out to test against USGS web calculator; it checked out
//				if(i==0) {
//					System.out.println("TEST VALUES\n\n");
//					for(int j=0;j<curveLinearXvalues.size();j++)
//						System.out.print(curveLinearXvalues.getX(j)+",");
//					System.out.println("\n\n");
//					for(int j=0;j<curveLinearXvalues.size();j++)
//						System.out.print((float)curveLinearXvalues.getY(j)+",");
//					System.out.println("\n\nTest RTGM = "+rtgm.get());
//				}
			} catch (Exception e) {
				System.out.println("RTGM Error; Hazard curve is:\n"+rateCurveLinearXvaluesCleaned);
				System.out.println("Location is:\n"+rtgmGriddedDataSetArray.getLocation(i));
				e.printStackTrace();
				System.exit(0);
			}
		}
		
		return rtgmGriddedDataSetArray;
	}
	
	
	public void mkRatioVsHazardScatterPlotTest() {
		GriddedGeoDataSet dataSet1 = getIML_ForProbGeoDataset(sa0pt2_poisson_curves, 0.02); 
		GriddedGeoDataSet dataSet2 = getIML_ForProbGeoDataset(sa0pt2_full_mean_curves, 0.02); 
		
		DefaultXY_DataSet scatterData = new DefaultXY_DataSet();
		for(int i=0;i<dataSet1.size();i++)
			scatterData.set(dataSet2.get(i), dataSet1.get(i)/dataSet2.get(i));
		
		String dirName = MakeFigures.outputDirString+"/RatioVsHazardScatterPlotTest";
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdir();
		
		String fileNamePrefix = dirName+"/RatioVsHazardScatterPlotTest";

		
		PlottingUtils.writeAndOrPlotFuncs(scatterData, new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.RED), 
				"Ratio vs 2in50 IML", "5-Hz SA", "Ratio", null, null, true, false, fileNamePrefix, true);
	}
	
	/**
	 * This compares Poisson to Full TD 2in50 5 Hz SA ratios with the 
	 * CharFactor plotted in Fig 9c or the UCERF3-ETAS paper.
	 */
	public void mkRatioVsCharFactorScatterPlotTest() {
		
		String dataFileName = "/Users/field/Field_Other/CEA_WGCEP/UCERF3/UCERF3-ETAS/BSSA_PaperStuff/Figures/Fig09_CubeRatesEtcFigs/CharFactorAtDepth7km_Poiss_withCorr/map_data.txt";

		GriddedGeoDataSet charFactorDataset = getGriddedGeoDataSet();
		File file = new File(dataFileName);
		List<String> fileLines;
		
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
			for(int i=0; i<fileLines.size();i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split("\t");
				double lat = Double.parseDouble(split[0]);
				double lon = Double.parseDouble(split[1]);
				double charFactor = Double.parseDouble(split[2]);
				int index = charFactorDataset.indexOf(new Location(lat,lon));
//				if(index != i)
//					throw new RuntimeException("Problem with index; index="+index+"; i="+i);
				if(index == -1) {
					continue;
				}
				charFactorDataset.set(index, Math.pow(10,charFactor));
			}
//			System.out.println("numFilled="+numFilled+"; charFactorDataset.size()="+charFactorDataset.size()+"\nfileLines.size()="+fileLines.size());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		GriddedGeoDataSet dataSet1 = getIML_ForProbGeoDataset(sa0pt2_poisson_curves, 0.02); 
		GriddedGeoDataSet dataSet2 = getIML_ForProbGeoDataset(sa0pt2_full_mean_curves, 0.02); 
		
//		System.out.println("dataSet1.size()="+dataSet1.size());

		DefaultXY_DataSet scatterData = new DefaultXY_DataSet();
		for(int i=0;i<dataSet1.size();i++)
			scatterData.set(charFactorDataset.get(i), dataSet1.get(i)/dataSet2.get(i));
		Range xRange = new Range(1e-3,1e3);
		Range yRange = new Range(0.95,1.6);
		
		String dirName = MakeFigures.outputDirString+"/RatioVsCharFactorScatterPlotTest";
		File dir = new File(dirName);
		if(!dir.exists())
			dir.mkdir();
		
		String fileNamePrefix = dirName+"/RatioVsCharFactorScatterPlotTest";

		
		PlottingUtils.writeAndOrPlotFuncs(scatterData, new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.RED), 
				null, "Char Factor", "Ratio", xRange, yRange, true, false, fileNamePrefix, true);
	}

	
	
	public void doAnalysis() {
		
		// verify ratio plots 
//		compareCases(pga_poisson_curves, pga_poisson_curves, 0.02, "Same Data Test, 2in50 PGA");

		// this is to compare a hazard map with that in Figure 3b of Powers and Field (2015, Eqk Spectra, Vol 31, pages S177-S200)
		// it looks as close as expected
//		try {
//			plotAndWriteHapMap(getIML_ForProbGeoDataset(sa1pt0_declustered_mean_curves, 0.02), "2% in 50 1-Hz SA", "", "test");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
//		// This is for Figure 6 Map,
//		// which has to be fetched from:
//		// /Users/field/workspace/git/opensha-dev/src/scratch/UCERF3/data/scratch/GMT/Figure6_Map
//		String title = "Ratio of Full Poisson to Full TD, 2% in 50 yr 5-Hz SA";
//		try {
//			plotAndWriteRatioMap(
//					getIML_ForProbGeoDataset(sa0pt2_poisson_curves, 0.02), 
//					getIML_ForProbGeoDataset(sa0pt2_full_mean_curves, 0.02), 
//					title, title, "Figure6_Map");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//
//		try {
//			plotAndWriteRatioMap(
//					getIML_ForProbGeoDataset(sa0pt2_declustered_mean_curves, 0.02), 
//					getIML_ForProbGeoDataset(sa0pt2_poisson_curves, 0.02), 
//					title, title, "GK_DeclusteredVsFullPois_Map");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		
//		System.exit(0);
		
		
//		// FOR 2% in 50-years: *********************************
//		
//		String infoString="";
//		
//		String imtString = "PGA";
////		infoString += compareCases(pga_randomized_mean_curves, pga_poisson_curves, 0.02, "Randomized vs Full Poisson, 2in50 PGA", imtString);
//		infoString += compareCases(pga_poisson_curves, pga_full_mean_curves, 0.02, "Full Poisson vs Full TD, 2in50 PGA", imtString);
////		System.exit(0);
//		infoString += compareCases(pga_declustered_mean_curves, pga_full_mean_curves, 0.02, "GK Declustered vs Full TD, 2in50 PGA", imtString);
//		infoString += compareCases(pga_declustered_mean_curves, pga_poisson_curves, 0.02, "GK Declustered vs Full Poisson, 2in50 PGA", imtString);
//
//		imtString = "5-Hz SA";
////		infoString += compareCases(sa0pt2_randomized_mean_curves, sa0pt2_poisson_curves, 0.02, "Randomized vs Full Poisson, 2in50 0.2secSA", imtString);
//		infoString += compareCases(sa0pt2_poisson_curves, sa0pt2_full_mean_curves, 0.02, "Full Poisson vs Full TD, 2in50 0.2secSA", imtString);
//		infoString += compareCases(sa0pt2_declustered_mean_curves, sa0pt2_full_mean_curves, 0.02, "GK Declustered vs Full TD, 2in50 0.2secSA", imtString);
//		infoString += compareCases(sa0pt2_declustered_mean_curves, sa0pt2_poisson_curves, 0.02, "GK Declustered vs Full Poisson, 2in50 0.2secSA", imtString);
//		
////		infoString += compareCases(getRTGM_Dataset(sa0pt2_randomized_mean_curves, 0.2),getRTGM_Dataset(sa0pt2_poisson_curves, 0.2, imtString),
////				"Randomized vs Full Poisson, RTGM 0.2secSA");
//		infoString += compareCases(getRTGM_Dataset(sa0pt2_poisson_curves, 0.2),getRTGM_Dataset(sa0pt2_full_mean_curves, 0.2),
//				"Full Poisson vs Full TD, RTGM 0.2secSA", imtString);
//		infoString += compareCases(getRTGM_Dataset(sa0pt2_declustered_mean_curves, 0.2),getRTGM_Dataset(sa0pt2_full_mean_curves, 0.2),
//				"GK Declustered vs Full TD, RTGM 0.2secSA", imtString);
//		infoString += compareCases(getRTGM_Dataset(sa0pt2_declustered_mean_curves, 0.2),
//				getRTGM_Dataset(sa0pt2_poisson_curves, 0.2),"GK Declustered vs Full Poisson, RTGM 0.2secSA", imtString);
//
//		imtString = "1-Hz SA";
////		infoString += compareCases(sa1pt0_randomized_mean_curves, sa1pt0_poisson_curves, 0.02, "Randomized vs Full Poisson, 2in50 1.0secSA", imtString);
//		infoString += compareCases(sa1pt0_poisson_curves, sa1pt0_full_mean_curves, 0.02, "Full Poisson vs Full TD, 2in50 1.0secSA", imtString);
//		infoString += compareCases(sa1pt0_declustered_mean_curves, sa1pt0_full_mean_curves, 0.02, "GK Declustered vs Full TD, 2in50 1.0secSA", imtString);
//		infoString += compareCases(sa1pt0_declustered_mean_curves, sa1pt0_poisson_curves, 0.02, "GK Declustered vs Full Poisson, 2in50 1.0secSA", imtString);
//
////		infoString += compareCases(getRTGM_Dataset(sa1pt0_randomized_mean_curves, 0.2),getRTGM_Dataset(sa1pt0_poisson_curves, 0.2, imtString),
////				"Randomized vs Full Poisson, RTGM 1.0secSA");
//		infoString += compareCases(getRTGM_Dataset(sa1pt0_poisson_curves, 1.0),getRTGM_Dataset(sa1pt0_full_mean_curves, 1.0),
//				"Full Poisson vs Full TD, RTGM 1.0secSA", imtString);
//		infoString += compareCases(getRTGM_Dataset(sa1pt0_declustered_mean_curves, 1.0),getRTGM_Dataset(sa1pt0_full_mean_curves, 1.0),
//				"GK Declustered vs Full TD, RTGM 1.0secSA", imtString);
//		infoString += compareCases(getRTGM_Dataset(sa1pt0_declustered_mean_curves, 1.0), getRTGM_Dataset(sa1pt0_poisson_curves, 1.0),
//				"GK Declustered vs Full Poisson, RTGM 1.0secSA", imtString);
//
//		imtString = "0.2-Hz SA";
////		infoString += compareCases(sa5pt0_randomized_mean_curves, sa5pt0_poisson_curves, 0.02, "Randomized vs Full Poisson, 2in50 5.0secSA", imtString);
//		infoString += compareCases(sa5pt0_poisson_curves, sa5pt0_full_mean_curves, 0.02, "Full Poisson vs Full TD, 2in50 5.0secSA", imtString);
//		infoString += compareCases(sa5pt0_declustered_mean_curves, sa5pt0_full_mean_curves, 0.02, "GK Declustered vs Full TD, 2in50 5.0secSA", imtString);
//		infoString += compareCases(sa5pt0_declustered_mean_curves, sa5pt0_poisson_curves, 0.02, "GK Declustered vs Full Poisson, 2in50 5.0secSA", imtString);
//		
//		System.out.println("2%in 50:\nName\tmean\tmin\tmax\tmaxY");
//		System.out.println(infoString);
		
		
		
		
		
		// FOR 40% in 50-years: *********************************
		
		String infoString="";
		
		String imtString = "PGA";
//		infoString += compareCases(pga_randomized_mean_curves, pga_poisson_curves, 0.4, "Randomized vs Full Poisson, 40in50 PGA", imtString);
		infoString += compareCases(pga_poisson_curves, pga_full_mean_curves, 0.4, "Full Poisson vs Full TD, 40in50 PGA", imtString);
//		System.exit(0);
		infoString += compareCases(pga_declustered_mean_curves, pga_full_mean_curves, 0.4, "GK Declustered vs Full TD, 40in50 PGA", imtString);
		infoString += compareCases(pga_declustered_mean_curves, pga_poisson_curves, 0.4, "GK Declustered vs Full Poisson, 40in50 PGA", imtString);

		imtString = "5-Hz SA";
//		infoString += compareCases(sa0pt2_randomized_mean_curves, sa0pt2_poisson_curves, 0.4, "Randomized vs Full Poisson, 40in50 0.2secSA", imtString);
		infoString += compareCases(sa0pt2_poisson_curves, sa0pt2_full_mean_curves, 0.4, "Full Poisson vs Full TD, 40in50 0.2secSA", imtString);
		infoString += compareCases(sa0pt2_declustered_mean_curves, sa0pt2_full_mean_curves, 0.4, "GK Declustered vs Full TD, 40in50 0.2secSA", imtString);
		infoString += compareCases(sa0pt2_declustered_mean_curves, sa0pt2_poisson_curves, 0.4, "GK Declustered vs Full Poisson, 40in50 0.2secSA", imtString);
		

		imtString = "1-Hz SA";
//		infoString += compareCases(sa1pt0_randomized_mean_curves, sa1pt0_poisson_curves, 0.4, "Randomized vs Full Poisson, 40in50 1.0secSA", imtString);
		infoString += compareCases(sa1pt0_poisson_curves, sa1pt0_full_mean_curves, 0.4, "Full Poisson vs Full TD, 40in50 1.0secSA", imtString);
		infoString += compareCases(sa1pt0_declustered_mean_curves, sa1pt0_full_mean_curves, 0.4, "GK Declustered vs Full TD, 40in50 1.0secSA", imtString);
		infoString += compareCases(sa1pt0_declustered_mean_curves, sa1pt0_poisson_curves, 0.4, "GK Declustered vs Full Poisson, 40in50 1.0secSA", imtString);

		imtString = "0.2-Hz SA";
//		infoString += compareCases(sa5pt0_randomized_mean_curves, sa5pt0_poisson_curves, 0.4, "Randomized vs Full Poisson, 40in50 5.0secSA", imtString);
		infoString += compareCases(sa5pt0_poisson_curves, sa5pt0_full_mean_curves, 0.4, "Full Poisson vs Full TD, 40in50 5.0secSA", imtString);
		infoString += compareCases(sa5pt0_declustered_mean_curves, sa5pt0_full_mean_curves, 0.4, "GK Declustered vs Full TD, 40in50 5.0secSA", imtString);
		infoString += compareCases(sa5pt0_declustered_mean_curves, sa5pt0_poisson_curves, 0.4, "GK Declustered vs Full Poisson, 40in50 5.0secSA", imtString);
		
		System.out.println("40% in 50:\nName\tmean\tmin\tmax\tmaxY");
		System.out.println(infoString);

		
		
		
	}
	



	
	public static void main(String[] args) {
		
		MapDataAnalysis analysis = new MapDataAnalysis();
		// this generates the figures for the paper
		analysis.doAnalysis();
		
		// With respect to high values in the map of Ratio of Full Poisson to Full TD, 2% in 50 yr 5−Hz SA (Figure 6)
		// This demonstrates that low hazard areas have higher discrepancies (low hazard is dominated by smaller events):
//		analysis.mkRatioVsHazardScatterPlotTest();
		// And this demonstrates that high ratios correlate with high char-factor values
//		analysis.mkRatioVsCharFactorScatterPlotTest();

	}
	
	

}
