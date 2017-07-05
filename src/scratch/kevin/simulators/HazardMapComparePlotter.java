package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class HazardMapComparePlotter {

	public static void main(String[] args) throws Exception {
		double minMag;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_23-bruce2142-vs-ucerf3-m7.0"); minMag = 7.0;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m5.5"); minMag = 5.5;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m6.5"); minMag = 6.5;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-bruce2142-vs-ucerf3-m6.0"); minMag = 6.0;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_1-vs-ucerf3-m6.5"); minMag = 6.5;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_1-vs-ucerf3-m7.0"); minMag = 7.0;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_2-vs-ucerf3-m6.5"); minMag = 6.5;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_slipWeakening_calibrated_2-vs-ucerf3-m7.0"); minMag = 7.0;
//		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_shortTestCatalog-vs-ucerf3-m6.5"); minMag = 6.5;
		File jobDir = new File("/home/kevin/Simulators/hazard/2017_06_28-jacqui_shortTestCatalog-vs-ucerf3-m7.0"); minMag = 7.0;
		
		// fallback
		Map<Double, File> u3Files = Maps.newHashMap();
		u3Files.put(5.5, new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m5.5/ucerf3/curves/imrs1.bin"));
		u3Files.put(6.0, new File("/home/kevin/Simulators/hazard/2017_06_28-bruce2142-vs-ucerf3-m6.0/ucerf3/curves/imrs1.bin"));
		u3Files.put(6.5, new File("/home/kevin/Simulators/hazard/2017_06_27-bruce2142-vs-ucerf3-m6.5/ucerf3/curves/imrs1.bin"));
		u3Files.put(7.0, new File("/home/kevin/Simulators/hazard/2017_06_23-bruce2142-vs-ucerf3-m7.0/ucerf3/curves/imrs1.bin"));
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		
		boolean plotMaps = false;
		boolean plotHist = true;
		
		// the point on the hazard curve we are plotting
		boolean isProbAtIML = false;
		double level = 0.0004;
		String durationLabel = "2% in 50 yrs";
		String fileLabel = "2pin50";
		
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
		
		GriddedGeoDataSet rsqsimData = loadFromBinary(gridReg, rsqsimCurves, isProbAtIML, level);
		GriddedGeoDataSet u3Data = loadFromBinary(gridReg, u3Curves, isProbAtIML, level);
		
		CPT hazardCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		hazardCPT = hazardCPT.rescale(0d, 1.2d);
		if (plotMaps) {
			plot(jobDir, "rsqsim_"+fileLabel+"_"+(float)minMag, rsqsimData, region,
					(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "RSQSim M"+(float)minMag+", "+durationLabel, hazardCPT, false);
			plot(jobDir, "u3_"+fileLabel+"_"+(float)minMag, u3Data, region,
					(double)hazardCPT.getMinValue(), (double)hazardCPT.getMaxValue(), "UCERF3 M"+(float)minMag+", "+durationLabel, hazardCPT, false);
		}
		
		GriddedGeoDataSet ratioData = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<gridReg.getNodeCount(); i++)
			ratioData.set(i, rsqsimData.get(i) / u3Data.get(i));
		ratioData.log10();
		
		if (plotMaps) {
			CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance();
			ratioCPT = ratioCPT.rescale(-0.2d, 0.2d);
			plot(jobDir, "ratio_log_tight_"+fileLabel+"_"+(float)minMag, ratioData, region, (double)ratioCPT.getMinValue(), (double)ratioCPT.getMaxValue(),
				"Log10(RSQSim / UCERF3), M"+(float)minMag, ratioCPT, false);
			ratioCPT = ratioCPT.rescale(-0.5d, 0.5d);
			plot(jobDir, "ratio_log_"+fileLabel+"_"+(float)minMag, ratioData, region, (double)ratioCPT.getMinValue(), (double)ratioCPT.getMaxValue(),
					"Log10(RSQSim / UCERF3), M"+(float)minMag, ratioCPT, false);
		}
		
		if (plotHist) {
			plotHists(u3Data, rsqsimData, jobDir, "hist_2d_log_"+fileLabel, true, durationLabel);
			plotHists(u3Data, rsqsimData, jobDir, "hist_2d_"+fileLabel, false, durationLabel);
		}
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
	
	private static void plot(File outputDir, String prefix, GriddedGeoDataSet data, Region region,
			Double customMin, Double customMax, String label, CPT cpt, boolean rescaleCPT)
					throws GMT_MapException, IOException {
		
		System.out.println("Creating map instance...");
		GMT_Map map = new GMT_Map(region, data, data.getRegion().getLatSpacing(), cpt);
		
		map.setCustomLabel(label);
//		map.setTopoResolution(TopographicSlopeFile.US_SIX);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setLogPlot(false);
		map.setDpi(300);
		map.setCustomScaleMin(customMin);
		map.setCustomScaleMax(customMax);
		map.setRescaleCPT(rescaleCPT);
		map.setBlackBackground(false);
		
		System.out.println("Making map...");
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
	}
	
	private static void plotHists(GriddedGeoDataSet ucerf3, GriddedGeoDataSet rsqsim, File outputDir, String prefix, boolean log,
			String label) throws IOException {
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
			rsqsim.log10();
			ucerf3.log10();
			minZ = Math.log10(minZ);
			maxZ = Math.log10(maxZ);
			delta = 0.05;
		}
		
		HistogramFunction sampleHist = HistogramFunction.getEncompassingHistogram(minZ, maxZ, delta);
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(sampleHist.size(), sampleHist.size(), minZ, minZ, delta);
		
		List<Double> logRatioVals = Lists.newArrayList();
		
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
		
		String xAxisLabel = "UCERF3 "+label;
		String yAxisLabel = "RSQSim "+label;
		if (log) {
			xAxisLabel = "Log10 "+xAxisLabel;
			yAxisLabel = "Log10 "+yAxisLabel;
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
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
		xyzGP.drawPlot(xyzSpec, false, false, plotRange, plotRange);
		// write plot
		xyzGP.getChartPanel().setSize(1000, 800);
		xyzGP.saveAsPNG(new File(outputDir, prefix+"_hist2D.png").getAbsolutePath());
		xyzGP.saveAsPDF(new File(outputDir, prefix+"_hist2D.pdf").getAbsolutePath());
		
		if (log) {
			double[] logRatioArray = Doubles.toArray(logRatioVals);
			HistogramFunction oneDHist = HistogramFunction.getEncompassingHistogram(
					StatUtils.min(logRatioArray), StatUtils.max(logRatioArray), 0.025);
			for (double val : logRatioArray)
				oneDHist.add(val, 1d);
			double mean = StatUtils.mean(logRatioArray);
			double stdDev = Math.sqrt(StatUtils.variance(logRatioArray));
			funcs = Lists.newArrayList();
			chars = Lists.newArrayList();
			
			oneDHist.normalizeBySumOfY_Vals();
			oneDHist.setName("Histogram");
			
			funcs.add(oneDHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.DARK_GRAY));
			
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			meanLine.set(mean, 0d);
			meanLine.set(mean, oneDHist.getMaxY());
			meanLine.setName("Mean: "+(float)mean+", Std Dev: "+(float)stdDev);
			
			funcs.add(meanLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Hazard 1D Histogram", "Log10(RSQSim/UCERF3)", "Fraction");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			gp.drawGraphPanel(spec, false, false);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(new File(outputDir, "hist_1d.png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, "hist_1d.pdf").getAbsolutePath());
		}
	}

}
