package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;

import com.google.common.base.Preconditions;

public class EALConvergenceResults {
	
	private static void makePortfolioMap(Portfolio portfolio, File dir) throws IOException, GMT_MapException, RuntimeException {
		double inc = 0.1;
		GriddedRegion region = new CaliforniaRegions.RELM_TESTING_GRIDDED(inc);
		GriddedGeoDataSet numMap = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet valueMap = new GriddedGeoDataSet(region, true);
		
		double max = 0;
		
		for (Asset asset : portfolio.getAssetList()) {
			Location loc = asset.getLocation();
			Preconditions.checkNotNull(loc);
			int index = numMap.indexOf(loc);
			double val = numMap.get(index);
			if (val > max)
				max = val;
			val++;
			numMap.set(index, val);
			
			val = valueMap.get(index);
			valueMap.set(index, val+inMillions(asset.getValue()));
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(1, max);
		cpt.setBelowMinColor(Color.GRAY);
		
		GMT_Map map = new GMT_Map(region, numMap, inc, cpt);
		map.setLogPlot(false);
		
		GMT_MapGenerator gmt = new GMT_MapGenerator();
		gmt.setParameter(GMT_MapGenerator.LOG_PLOT_NAME, false);
		gmt.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME, 1e-5);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME, max);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME, GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		map.setCustomScaleMin(1d);
		map.setCustomScaleMax(max);
		map.setUseGMTSmoothing(false);
		map.setCustomLabel("# Assets Per Cell");
		
		String url = gmt.makeMapUsingServlet(map, "", null);
		
		String baseURL = url.substring(0, url.lastIndexOf('/')+1);
		
		File pngFile = new File(dir, "assets_map.png");
		FileUtils.downloadURL(baseURL+"map.png", pngFile);
		
		map = new GMT_Map(region, valueMap, inc, cpt);
		map.setLogPlot(false);
		map.setCustomScaleMin(1d);
		map.setCustomScaleMax(valueMap.getMaxZ());
		map.setUseGMTSmoothing(false);
		map.setCustomLabel("Asset Value Per Cell ($ million)");
		
		url = gmt.makeMapUsingServlet(map, "", null);
		
		baseURL = url.substring(0, url.lastIndexOf('/')+1);
		
		pngFile = new File(dir, "assets_value_map.png");
		FileUtils.downloadURL(baseURL+"map.png", pngFile);
	}
	
	private static void makeAssetVariabilityMap(Portfolio portfolio, List<PortfolioEAL> eals, File dir) throws IOException, GMT_MapException, RuntimeException {
		double inc = 0.1;
		GriddedRegion region = new CaliforniaRegions.RELM_TESTING_GRIDDED(inc);
//		GriddedGeoDataSet[] data = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet[] datas = new GriddedGeoDataSet[eals.size()];
		for (int i=0; i<eals.size(); i++)
			datas[i] = new GriddedGeoDataSet(region, true);
		
		ArrayList<Asset> assets = portfolio.getAssetList();
		
		for (int i=0; i<assets.size(); i++) {
			Asset asset = assets.get(i);
			Location loc = asset.getLocation();
			Preconditions.checkNotNull(loc);
			int index = datas[0].indexOf(loc);
			for (int j=0; j<eals.size(); j++) {
				double val = datas[j].get(index);
				datas[j].set(index, val+inMillions(eals.get(j).assetEALs.get(i)));
			}
		}
		
		double max = 0;
		
		GriddedGeoDataSet meanData = new GriddedGeoDataSet(region, true);
		GriddedGeoDataSet meanOverStdDev = new GriddedGeoDataSet(region, true);
		for (int i=0; i<meanOverStdDev.size(); i++) {
			double[] vals = new double[datas.length];
			for (int j=0; j<datas.length; j++)
				vals[j] = datas[j].get(i);
			double mean = StatUtils.mean(vals);
			if (mean == 0) {
				meanOverStdDev.set(i, Double.NaN);
				meanData.set(i, Double.NaN);
				continue;
			}
			double stdDev = Math.sqrt(StatUtils.variance(vals, mean));
			double val = mean / stdDev;
			if (val > max)
				max = val;
			meanData.set(i, mean);
			meanOverStdDev.set(i, val);
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(1, max);
		cpt.setBelowMinColor(Color.GRAY);
		
		GMT_Map map = new GMT_Map(region, meanOverStdDev, inc, cpt);
		map.setLogPlot(true);
		
		GMT_MapGenerator gmt = new GMT_MapGenerator();
		gmt.setParameter(GMT_MapGenerator.LOG_PLOT_NAME, true);
		gmt.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME, 1e-5);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME, max);
		gmt.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME, GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		map.setCustomScaleMin(0d);
//		map.setCustomScaleMax(600d);
		map.setCustomScaleMin(1d);
		map.setCustomScaleMax(3d);
		map.setUseGMTSmoothing(false);
		map.setCustomLabel("Mean / Std Dev of EAL Per Cell");
		
		String url = gmt.makeMapUsingServlet(map, "", null);
		
		String baseURL = url.substring(0, url.lastIndexOf('/')+1);
		
		File pngFile = new File(dir, "assets_variability_map.png");
		FileUtils.downloadURL(baseURL+"map.png", pngFile);
		
		max = meanData.getMaxZ();
		
		gmt.setParameter(GMT_MapGenerator.LOG_PLOT_NAME, false);
		map = new GMT_Map(region, meanData, inc, cpt);
		map.setLogPlot(false);
		map.setCustomScaleMin(0d);
		map.setCustomScaleMax(13d);
		map.setUseGMTSmoothing(false);
		map.setCustomLabel("Mean EAL Per Cell ($ Millions)");
		
		url = gmt.makeMapUsingServlet(map, "", null);
		
		baseURL = url.substring(0, url.lastIndexOf('/')+1);
		
		pngFile = new File(dir, "assets_mean_map.png");
		FileUtils.downloadURL(baseURL+"map.png", pngFile);
	}
	
	private static List<PortfolioEAL> loadEALs(File dir) throws IOException {
		File[] dirFiles = dir.listFiles();
		Arrays.sort(dirFiles, new FileNameComparator());
		
		ArrayList<PortfolioEAL> eals = new ArrayList<PortfolioEAL>();
		
		for (File file : dirFiles) {
			if (!file.getName().endsWith(".txt"))
				continue;
			
			ArrayList<String> lines = FileUtils.loadFile(file.getAbsolutePath());
			
			String line0 = lines.get(0).trim();
			line0 = line0.substring(line0.indexOf(':')+1).trim();
			double portEAL = Double.parseDouble(line0);
			
			ArrayList<IndexedEAL> assetEALs = new ArrayList<IndexedEAL>();
			
			for (String line : lines.subList(1, lines.size())) {
				line = line.trim();
				if (line.isEmpty())
					continue;
				String[] splt = line.split(",");
				int index = Integer.parseInt(splt[0]);
				double eal = Double.parseDouble(splt[1]);
				assetEALs.add(new IndexedEAL(index, eal));
			}
			
			Collections.sort(assetEALs);
			
			ArrayList<Double> ealDoubles = new ArrayList<Double>();
			for (IndexedEAL eal : assetEALs)
				ealDoubles.add(eal.eal);
			
			eals.add(new PortfolioEAL(portEAL, ealDoubles));
		}
		
		return eals;
	}
	
	private static class IndexedEAL implements Comparable<IndexedEAL> {
		private int index;
		private double eal;
		public IndexedEAL(int index, double eal) {
			this.index = index;
			this.eal = eal;
		}
		@Override
		public int compareTo(IndexedEAL o) {
			return Double.compare(index, o.index);
		}
	}
	
	private static class PortfolioEAL {
		private double portfolioEAL;
		private List<Double> assetEALs;
		public PortfolioEAL(double portfolioEAL, List<Double> assetEALs) {
			this.portfolioEAL = portfolioEAL;
			this.assetEALs = assetEALs;
		}
	}
	
	private static void plotCumDist(List<PortfolioEAL> eals, File dir) throws IOException {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (PortfolioEAL eal : eals)
			track.addValue(eal.portfolioEAL);
		
		double min = inMillions(track.getMin());
		min = ((double)(int)(min * 100d)) / 100d;
		double delta = 0.5;
		double max = inMillions(track.getMax());
		int num = 0;
		for (double val=min; val<max; val+=delta)
			num++;
		HistogramFunction hist = new HistogramFunction(min, num, delta);
		
		System.out.println(track);
		System.out.println("Min: "+min+", max: "+hist.getMaxX()+", num: "+num);
		
		for (PortfolioEAL eal : eals) {
			double billions = inMillions(eal.portfolioEAL);
//			System.out.println("val: "+billions);
			double prevY = hist.getY(billions);
			hist.set(billions, prevY+1d);
		}
		
		ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
		funcs.add(hist);
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		
		gp.drawGraphPanel("EAL ($ million)", "Num", funcs, chars, "EAL Histogram");
		gp.getChartPanel().setSize(1000, 800);
		gp.setBackground(Color.WHITE);
		gp.saveAsPNG(new File(dir, "eal_hist.png").getAbsolutePath());
		gp.getChartPanel().setSize(500, 400);
		gp.setBackground(Color.WHITE);
		gp.saveAsPNG(new File(dir, "eal_hist_small.png").getAbsolutePath());
		
		hist.normalizeBySumOfY_Vals();
		HistogramFunction cumDist = hist.getCumulativeDistFunction();
		
		funcs.set(0, cumDist);
		chars.set(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		gp.setAutoRange();
		gp.drawGraphPanel("EAL ($ million)", "Probability [EAL <= X]", funcs, chars, "EAL Prob Dist");
		gp.getChartPanel().setSize(1000, 800);
		gp.setBackground(Color.WHITE);
		gp.saveAsPNG(new File(dir, "eal_prob_dist.png").getAbsolutePath());
		gp.getChartPanel().setSize(500, 400);
		gp.setBackground(Color.WHITE);
		gp.saveAsPNG(new File(dir, "eal_prob_dist_small.png").getAbsolutePath());
	}
	
	private static double inMillions(double eal) {
		// it's already in thousands
		return eal / 1000d;
	}
	
	public static void main(String[] args) throws IOException, GMT_MapException, RuntimeException {
//		File dir = new File("/home/kevin/OpenSHA/UCERF3/eal/2012_05_03-eal-tests-apriori-1000/CB2008");
		File dir = new File("/home/kevin/OpenSHA/UCERF3/eal/2012_05_17-rtgm-tests-apriori-1000/CB2008");
		
		List<PortfolioEAL> eals = loadEALs(dir);
		
//		plotCumDist(eals, dir);
		
		File portfolioFile = new File(dir.getParentFile(), "portfolio.csv");
		Portfolio port = Portfolio.createPortfolio(portfolioFile);
		
		int numAssets = port.getAssetList().size();
		for (int i=0; i<eals.size(); i++) {
			PortfolioEAL eal = eals.get(i);
			Preconditions.checkState(eal.assetEALs.size() == numAssets, "oops ("+i+")! "+eal.assetEALs.size()+" != "+numAssets);
		}
		
		makePortfolioMap(port, dir);
		makeAssetVariabilityMap(port, eals, dir);
	}

}
