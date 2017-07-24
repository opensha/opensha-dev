package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class HazardMapComparePlotter {
	
	private static final boolean map_parallel = true;

	public static void main(String[] args) throws Exception {
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_23-bruce2142-vs-ucerf3-m7.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m5.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m6.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-bruce2142-vs-ucerf3-m6.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_1-vs-ucerf3-m6.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_1-vs-ucerf3-m7.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_2-vs-ucerf3-m6.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_2-vs-ucerf3-m7.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_shortTestCatalog-vs-ucerf3-m6.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_shortTestCatalog-vs-ucerf3-m7.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_05-bruce2194-vs-ucerf3-m6.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_05-bruce2194-vs-ucerf3-m7.0");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2142-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2142-m6.5-sectArea0.3");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2142-m6.5-sectArea0.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2194-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2194-m6.5-sectArea0.3");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_13-bruce2194-m6.5-sectArea0.5");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-m6.5-sectArea0.2");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-m6.5-sectArea0.4");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-matchU3supra-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-matchU3supra-sectArea0.2");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-matchU3supra-sectArea0.4");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-m6.5-sectArea0.2");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-m6.5-sectArea0.4");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-matchU3supra-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-matchU3supra-sectArea0.2");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-matchU3supra-sectArea0.4");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_18-jacqui_shortTestCatalog-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_18-jacqui_slipWeakening_calibrated_2-m6.5-sectArea0.1");
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_07_20-jacqui_slipWeakening_calibrated_1-m6.5-sectArea0.1");
		
		List<File> jobDirs = Lists.newArrayList();
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_18-jacqui_slipWeakening_calibrated_2-m6.5-sectArea0.1"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_20-jacqui_slipWeakening_calibrated_1-m6.5-sectArea0.1"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2142-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_14-bruce2194-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2230-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2231-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2232-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2233-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2234-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2240-m6.5-sectArea0.2"));
		jobDirs.add(new File("/home/kevin/Simulators/hazard/2017_07_21-bruce2241-m6.5-sectArea0.2"));
//		
		// fallback
		Map<Double, File> u3Files = Maps.newHashMap();
		u3Files.put(5.5, new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m5.5/ucerf3/curves/imrs1.bin"));
		u3Files.put(6.0, new File("/home/kevin/Simulators/hazard/2017_06_28-bruce2142-vs-ucerf3-m6.0/ucerf3/curves/imrs1.bin"));
		u3Files.put(6.5, new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m6.5/ucerf3/curves/imrs1.bin"));
		u3Files.put(7.0, new File("/home/kevin/Simulators/hazard/2017_06_23-bruce2142-vs-ucerf3-m7.0/ucerf3/curves/imrs1.bin"));
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		File u3FullFile = new File("/home/kevin/OpenSHA/UCERF3/maps/2017_07_14-ucerf3-gridded-tests/full/curves/imrs1.bin");
		
		boolean plotMaps = false;
		boolean plotHist = true;
		boolean plotCurves = false;
		boolean plotDisagg = false;
		
		int[] mapRPs = { 1000, 2500, 10000};
		
		int[] histRPs = { 1000, 2500, 5000, 10000 };
		int histHighlightIndex = 1;
		
		Map<String, Location> curveLocs = new HashMap<>();
		curveLocs.put("Pasadena", new Location(34.148426, -118.17119));
		curveLocs.put("San Bernardino", new Location(34.064986, -117.29201));
		curveLocs.put("USC", new Location(34.0192, -118.286));
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			curveLocs.put(city.toString(), city.location());
		}
		
		String imt = "PGA (g)";
		
		for (int d=0; d<jobDirs.size(); d++) {
			File jobDir = jobDirs.get(d);
			System.out.println("Processing "+(d+1)+"/"+jobDirs.size()+": "+jobDir);
			File rsqsimFile = new File(jobDir, "rsqsim/curves/imrs1.bin");
			File u3File = new File(jobDir, "ucerf3/curves/imrs1.bin");
			if (!u3File.exists()) {
				u3File = u3Files.get(6.5);
				System.out.println("Using fallback UCERF3 file: "+u3File.getAbsolutePath());
			}
			
			BinaryHazardCurveReader rsqsimReader = new BinaryHazardCurveReader(rsqsimFile.getAbsolutePath());
			BinaryHazardCurveReader u3Reader = new BinaryHazardCurveReader(u3File.getAbsolutePath());
			
			Map<Location, ArbitrarilyDiscretizedFunc> rsqsimCurves = rsqsimReader.getCurveMap();
			Map<Location, ArbitrarilyDiscretizedFunc> u3Curves = u3Reader.getCurveMap();
			
			Map<Location, ArbitrarilyDiscretizedFunc> u3FullCurves = null;
			if (plotCurves && u3FullFile != null) {
				BinaryHazardCurveReader u3FullReader = new BinaryHazardCurveReader(u3FullFile.getAbsolutePath());
				u3FullCurves = u3FullReader.getCurveMap();
			}
			
			CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
			hazardCPT = hazardCPT.rescale(0d, 1.2d);
			
			for (int mapRP : mapRPs) {
				boolean isProbAtIML = false;
				double level = 1d/(double)mapRP;
				
				String durationLabel = mapRP+"yr";
				
				GriddedGeoDataSet rsqsimData = loadFromBinary(gridReg, rsqsimCurves, isProbAtIML, level);
				GriddedGeoDataSet u3Data = loadFromBinary(gridReg, u3Curves, isProbAtIML, level);
				
				if (plotMaps) {
					plotMaps(jobDir, "map_"+durationLabel+"_rsqsim", rsqsimData, region,
							(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "RSQSim, "+durationLabel, hazardCPT, false);
					plotMaps(jobDir, "map_"+durationLabel+"_u3"+durationLabel, u3Data, region,
							(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "UCERF3, "+durationLabel, hazardCPT, false);
				}
				
				GriddedGeoDataSet ratioData = new GriddedGeoDataSet(gridReg, false);
				for (int i=0; i<gridReg.getNodeCount(); i++)
					ratioData.set(i, rsqsimData.get(i) / u3Data.get(i));
				ratioData.log();
				
				if (plotMaps) {
					CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance();
					ratioCPT = ratioCPT.rescale(-0.2d, 0.2d);
					plotMaps(jobDir, "map_"+durationLabel+"_ratio_log_tight", ratioData, region, (double)ratioCPT.getMinValue(),
							(double)ratioCPT.getMaxValue(), "Ln(RSQSim / UCERF3), "+durationLabel, ratioCPT, false);
					ratioCPT = ratioCPT.rescale(-0.5d, 0.5d);
					plotMaps(jobDir, "map_"+durationLabel+"_ratio_log", ratioData, region, (double)ratioCPT.getMinValue(),
							(double)ratioCPT.getMaxValue(), "Ln(RSQSim / UCERF3), "+durationLabel, ratioCPT, false);
				}
			}
			
			if (plotHist) {
				plotHists(u3Curves, rsqsimCurves, gridReg, histRPs, histHighlightIndex, jobDir, true);
				plotHists(u3Curves, rsqsimCurves, gridReg, histRPs, histHighlightIndex, jobDir, false);
				plotMeanStdDevTrend(500, 30000, u3Curves, rsqsimCurves, gridReg, jobDir);
			}
			
			if (plotCurves) {
				File curveDir = new File(jobDir, "curves");
				Preconditions.checkState(curveDir.exists() || curveDir.mkdir());
				
				for (String siteName : curveLocs.keySet()) {
					Location siteLoc = curveLocs.get(siteName);
//					Location gridLoc = rsqsimData.getLocation(rsqsimData.indexOf(siteLoc));
					Location gridLoc = gridReg.locationForIndex(gridReg.indexForLocation(siteLoc));
					double dist = LocationUtils.horzDistance(siteLoc, gridLoc);
					System.out.println(siteName+" is "+(float)dist+" km away from nearest grid point");
					
					DiscretizedFunc u3Curve = u3Curves.get(gridLoc);
					DiscretizedFunc rsqsimCurve = rsqsimCurves.get(gridLoc);
					DiscretizedFunc u3FullCurve = null;
					if (u3FullCurves != null)
						u3FullCurve = u3FullCurves.get(gridLoc);
					
					plotCurves(u3Curve, u3FullCurve, rsqsimCurve, curveDir, siteName, imt, mapRPs);
				}
			}
		}
		
		waitOnMaps();
	}
	
	private static GriddedGeoDataSet loadFromBinary(GriddedRegion gridReg, Map<Location, ArbitrarilyDiscretizedFunc> curves,
			boolean isProbAtIML, double level) {
		GriddedGeoDataSet data = new GriddedGeoDataSet(gridReg, false);
		
		for (Location loc : curves.keySet()) {
			DiscretizedFunc curve = curves.get(loc);
			double val = HazardDataSetLoader.getCurveVal(curve, isProbAtIML, level);
			data.set(loc, val);
		}
		
		return data;
	}
	
	private static void plotMaps(File outputDir, String prefix, GriddedGeoDataSet data, Region region,
			Double customMin, Double customMax, String label, CPT cpt, boolean rescaleCPT)
					throws GMT_MapException, IOException {
		
		System.out.println("Creating map instance...");
		GMT_Map map = new GMT_Map(region, data, data.getRegion().getLatSpacing(), cpt);
		
		map.setCustomLabel(label);
//		map.setTopoResolution(TopographicSlopeFile.US_SIX);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setCoast(new CoastAttributes(Color.GRAY, 1));
		map.setLogPlot(false);
		map.setDpi(300);
		map.setCustomScaleMin(customMin);
		map.setCustomScaleMax(customMax);
		map.setRescaleCPT(rescaleCPT);
		map.setBlackBackground(false);
		
		System.out.println("Making map...");
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		
		Runnable run = new Runnable() {
			
			@Override
			public void run() {
				try {
					FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		};
		
		if (map_parallel) {
			if (mapExec == null) {
				mapExec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
				mapFutures = new ArrayList<>();
			}
			mapFutures.add(mapExec.submit(run));
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			run.run();
		}
	}
	
	private static void waitOnMaps() {
		if (mapFutures == null)
			return;
		System.out.println("Waiting on "+mapFutures.size()+" maps...");
		for (Future<?> future : mapFutures) {
			try {
				future.get();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		System.out.print("DONE");
		mapExec.shutdown();
	}
	
	private static ExecutorService mapExec;
	private static List<Future<?>> mapFutures;
	
	private static void plotHists(Map<Location, ArbitrarilyDiscretizedFunc> ucerf3Curves,
			Map<Location, ArbitrarilyDiscretizedFunc> rsqsimCurves, GriddedRegion gridReg,
			int[] returnPeriods, int hightlightIndex, File outputDir, boolean log) throws IOException {
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		List<List<Double>> logRatioValsList = new ArrayList<>();
		
		for (int rp : returnPeriods) {
			GriddedGeoDataSet ucerf3 = loadFromBinary(gridReg, ucerf3Curves, false, 1d/(double)rp);
			GriddedGeoDataSet rsqsim = loadFromBinary(gridReg, rsqsimCurves, false, 1d/(double)rp);
			
			// min non zero
			double minZ = Double.POSITIVE_INFINITY;
			for (int i=0; i<ucerf3.size(); i++) {
				if (ucerf3.get(i) > 0)
					minZ = Math.min(minZ, ucerf3.get(i));
				if (rsqsim.get(i) > 0)
					minZ = Math.min(minZ, rsqsim.get(i));
			}
			double maxZ = Math.max(ucerf3.getMaxZ(), rsqsim.getMaxZ());
			
			double delta = 0.02;
			if (log) {
				rsqsim = rsqsim.copy();
				ucerf3 = ucerf3.copy();
				rsqsim.log();
				ucerf3.log();
				minZ = Math.log(minZ);
				maxZ = Math.log(maxZ);
				delta = 0.05;
			}
			
			HistogramFunction sampleHist = HistogramFunction.getEncompassingHistogram(minZ, maxZ, delta);
			
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(sampleHist.size(), sampleHist.size(), minZ, minZ, delta);
			
			List<Double> logRatioVals = new ArrayList<>();
			logRatioValsList.add(logRatioVals);
			
			for (int i=0; i<ucerf3.size(); i++) {
				double u3 = ucerf3.get(i);
				double rs = rsqsim.get(i);
				if (log && (Double.isNaN(u3) || Double.isNaN(rs)))
					continue;
				if (!log && (u3 == 0 || rs == 0))
					continue;
				int xInd = sampleHist.getClosestXIndex(u3);
				int yInd = sampleHist.getClosestXIndex(rs);
				
				xyz.set(xInd, yInd, xyz.get(xInd, yInd) + 1d);
				
				if (log)
					// already in log space, subract
					logRatioVals.add(rs - u3);
			}
			
			xyz.log10();
			
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
			cpt = cpt.rescale(Math.log10(1d), xyz.getMaxZ());
			cpt.setNanColor(Color.WHITE);
			cpt.setBelowMinColor(Color.WHITE);
			
			String xAxisLabel = "UCERF3 1/"+rp+"yr";
			String yAxisLabel = "RSQSim 1/"+rp+"yr";
			if (log) {
				xAxisLabel = "Ln "+xAxisLabel;
				yAxisLabel = "Ln "+yAxisLabel;
			}
			
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, "Hazard Histogram", xAxisLabel, yAxisLabel, "Log10(Num)");
			
			Range plotRange = new Range(minZ - 0.5*delta, maxZ + 0.5*delta);
			
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
			oneToOne.set(plotRange.getLowerBound(), plotRange.getLowerBound());
			oneToOne.set(plotRange.getUpperBound(), plotRange.getUpperBound());
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			xyzSpec.setXYElems(funcs);
			xyzSpec.setXYChars(chars);
			
			XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
			xyzGP.drawPlot(xyzSpec, false, false, plotRange, plotRange);
			// write plot
			xyzGP.getChartPanel().setSize(1000, 800);
			String prefix = "hist_2d_"+rp+"yr";
			xyzGP.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			xyzGP.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		}
		
		CPT rpCPT = getRPlogCPT(returnPeriods);
		
		if (log) {
			boolean logY = true;
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double highlightMean = 0;
			double highlightStdDev = 0;
			
			double histMinNonZero = Double.POSITIVE_INFINITY;
			double histMax = 0;
			
			for (int i=0; i<returnPeriods.length; i++) {
				int rp = returnPeriods[i];
				Color color = rpCPT.getColor((float)Math.log10(rp));
				
				double[] logRatioArray = Doubles.toArray(logRatioValsList.get(i));
				HistogramFunction oneDHist = HistogramFunction.getEncompassingHistogram(
						StatUtils.min(logRatioArray), StatUtils.max(logRatioArray), 0.05);
				for (double val : logRatioArray)
					oneDHist.add(val, 1d);
				double mean = StatUtils.mean(logRatioArray);
				double stdDev = Math.sqrt(StatUtils.variance(logRatioArray));
				
				oneDHist.normalizeBySumOfY_Vals();
				oneDHist.setName(rp+"yr");
				
				for (Point2D pt : oneDHist)
					if (pt.getY() > 0)
						histMinNonZero = Math.min(histMinNonZero, pt.getY());
				histMax = Math.max(histMax, oneDHist.getMaxY());
				
				if (i == hightlightIndex) {
					funcs.add(0, oneDHist);
					chars.add(0, new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
					highlightMean = mean;
					highlightStdDev = stdDev;
				} else {
					funcs.add(oneDHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
				}
			}
			
			Range yRange = null;
			if (logY) {
				yRange = new Range(histMinNonZero*0.7, histMax*1.5);
			}
			
			if (hightlightIndex >= 0) {
				DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
				meanLine.set(highlightMean, 0d);
				if (logY) {
					meanLine.set(highlightMean, yRange.getLowerBound());
					meanLine.set(highlightMean, yRange.getUpperBound());
				}
				meanLine.set(highlightMean, histMax);
				meanLine.setName(returnPeriods[hightlightIndex]+"yr Mean: "+fourDigits.format(highlightMean)
					+", Std Dev: "+fourDigits.format(highlightStdDev));
				
				funcs.add(meanLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Hazard 1D Histogram", "Ln(RSQSim/UCERF3)", "Fraction");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, false, logY, null, yRange);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(new File(outputDir, "hist_1d.png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, "hist_1d.pdf").getAbsolutePath());
		}
	}
	
	private static final DecimalFormat fourDigits = new DecimalFormat("0.0000");
	
	private static void plotMeanStdDevTrend(int minRP, int maxRP, Map<Location, ArbitrarilyDiscretizedFunc> ucerf3Curves,
			Map<Location, ArbitrarilyDiscretizedFunc> rsqsimCurves, GriddedRegion gridReg, File outputDir) throws IOException {
		int numRPs = 20;
		
		EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc((double)minRP, (double)maxRP, numRPs);
		EvenlyDiscretizedFunc stdDevFunc = new EvenlyDiscretizedFunc((double)minRP, (double)maxRP, numRPs);
		
		double[] pga_thresholds = {0.0, 0.1, 0.2, 0.4, 0.8, 1.2};
		CPT pgaWeightCPT = new CPT(pga_thresholds[0], pga_thresholds[pga_thresholds.length-1], Color.GRAY, Color.BLACK);
		EvenlyDiscretizedFunc[] pgaWeightFuncs = new EvenlyDiscretizedFunc[pga_thresholds.length];
		for (int p=0; p<pga_thresholds.length; p++)
			pgaWeightFuncs[p] = new EvenlyDiscretizedFunc((double)minRP, (double)maxRP, numRPs);
		EvenlyDiscretizedFunc[] pgaWeightStdDevFuncs = new EvenlyDiscretizedFunc[pga_thresholds.length];
		for (int p=0; p<pga_thresholds.length; p++)
			pgaWeightStdDevFuncs[p] = new EvenlyDiscretizedFunc((double)minRP, (double)maxRP, numRPs);
		
		EvenlyDiscretizedFunc nehrpMeanFunc = new EvenlyDiscretizedFunc((double)minRP, (double)maxRP, numRPs);
		
		HashSet<Integer> nehrpGridIndexes = new HashSet<>();
		for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {
			Location loc = city.getSite().getLocation();
			nehrpGridIndexes.add(gridReg.indexForLocation(loc));
		}
		
		for (int i=0; i<numRPs; i++) {
			double rp = meanFunc.getX(i);
			GriddedGeoDataSet ucerf3 = loadFromBinary(gridReg, ucerf3Curves, false, 1d/rp);
			GriddedGeoDataSet rsqsim = loadFromBinary(gridReg, rsqsimCurves, false, 1d/rp);
			
			List<Double> ratioVals = new ArrayList<>();
			
			List<List<Double>> pgaWeights = new ArrayList<>();
			for (int p=0; p<pga_thresholds.length; p++)
				pgaWeights.add(new ArrayList<>());
			
			List<Double> nehrpRatios = new ArrayList<>();
			
			for (int j=0; j<ucerf3.size(); j++) {
				double u3 = ucerf3.get(j);
				double rs = rsqsim.get(j);
				if (u3 == 0 || rs == 0 || Double.isNaN(u3) || Double.isNaN(rs))
					continue;
				
				double ratio = Math.log(rs/u3);
				Preconditions.checkState(Double.isFinite(ratio), "Bad ratio: Ln(%s/%s) = %s", rs, u3, ratio);
				
				if (nehrpGridIndexes.contains(j))
					nehrpRatios.add(ratio);
				
				ratioVals.add(Math.log(rs/u3));
				for (int p=0; p<pga_thresholds.length; p++) {
					double weight = u3 - pga_thresholds[p];
					if (weight < 0)
						weight = 0;
					pgaWeights.get(p).add(weight);
				}
			}
			
			double[] ratioArray = Doubles.toArray(ratioVals);
			
			double mean = StatUtils.mean(ratioArray);
			double stdDev = Math.sqrt(StatUtils.variance(ratioArray));
			
			meanFunc.set(i, mean);
			stdDevFunc.set(i, stdDev);
			
			double nehrpMean = StatUtils.mean(Doubles.toArray(nehrpRatios));
			nehrpMeanFunc.set(i, nehrpMean);
			
			for (int p=0; p<pga_thresholds.length; p++) {
				ArbDiscrEmpiricalDistFunc func = new ArbDiscrEmpiricalDistFunc();
				for (int j=0; j<ratioVals.size(); j++)
					func.set(ratioVals.get(j), pgaWeights.get(p).get(j));
				pgaWeightFuncs[p].set(i, func.getMean());
				pgaWeightStdDevFuncs[p].set(i, func.getStdDev());
			}
			
//			System.out.println(rp+" "+mean+" "+stdDev+" "+ratioVals.size());
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int p=0; p<pga_thresholds.length; p++) {
			funcs.add(pgaWeightFuncs[p]);
			String name;
			if (p == 0 || p == pga_thresholds.length-1)
				name = "PGA_0="+(float)pga_thresholds[p];
			else
				name = (float)pga_thresholds[p]+"";
			pgaWeightFuncs[p].setName(name);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f+2f*(float)pga_thresholds[p],
					pgaWeightCPT.getColor((float)pga_thresholds[p])));
		}
		
		for (int p=0; p<pga_thresholds.length; p++) {
			funcs.add(pgaWeightStdDevFuncs[p]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f+2f*(float)pga_thresholds[p],
					pgaWeightCPT.getColor((float)pga_thresholds[p])));
		}
		
		funcs.add(nehrpMeanFunc);
		nehrpMeanFunc.setName("NEHRP City Mean");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		funcs.add(meanFunc);
		meanFunc.setName("Mean");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		funcs.add(stdDevFunc);
		stdDevFunc.setName("StdDev");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Mean/StdDev Trend", "Return Period (yr)", "Ln(RSQSim/UCERF3)");
		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		gp.drawGraphPanel(spec, false, false, null, null);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(new File(outputDir, "hist_0d.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "hist_0d.pdf").getAbsolutePath());
	}
	
	private static CPT getRPlogCPT(int[] rps) {
		int minRP = Integer.MAX_VALUE;
		int maxRP = 0;
		
		for (int rp : rps) {
			if (rp < minRP)
				minRP = rp;
			if (rp > maxRP)
				maxRP = rp;
		}
		return new CPT(Math.log10(minRP), Math.log10(maxRP), Color.LIGHT_GRAY, Color.BLACK);
	}
	
	private static void plotCurves(DiscretizedFunc ucerf3Compare, DiscretizedFunc ucerf3Full,
			DiscretizedFunc rsqsim, File outputDir, String siteName,
			String imt, int[] rps) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = new Range(1e-3, 1e1);
		Range yRange = new Range(1e-8, 1e0);
		
		if (ucerf3Full != null) {
			ucerf3Full.setName("UCERF3 Full");
			funcs.add(ucerf3Full);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		}
		
		ucerf3Compare.setName("UCERF3");
		funcs.add(ucerf3Compare);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		rsqsim.setName("RSQSim");
		funcs.add(rsqsim);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		if (rps != null && rps.length > 0) {
			CPT rpCPT = getRPlogCPT(rps);
			for (int rp : rps) {
				Color color = rpCPT.getColor((float)Math.log10(rp));
				double probLevel = 1d/(double)rp;
				DiscretizedFunc probLine = new ArbitrarilyDiscretizedFunc();
				probLine.set(xRange.getLowerBound(), probLevel);
				probLine.set(xRange.getUpperBound(), probLevel);
				probLine.setName(rp+"yr");
				funcs.add(probLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color));
			}
		}

		PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Hazard Curves", imt, "Annual Probabilitiy");
		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		String prefix = "curves_"+siteName.replaceAll(" ", "_");
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(new File(outputDir, prefix).getAbsolutePath()+".png");
		gp.saveAsPDF(new File(outputDir, prefix).getAbsolutePath()+".pdf");
		gp.saveAsTXT(new File(outputDir, prefix).getAbsolutePath()+".txt");
	}

}
