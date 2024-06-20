package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.ArbDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.base.Preconditions;

public class MultiCheckerboardPlot {

	public static void main(String[] args) throws IOException {
		List<File> refDirs = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		List<String> outPrefixes = new ArrayList<>();
		List<Boolean> stdDevs = new ArrayList<>();
		List<Boolean> doIndvs = new ArrayList<>();
		List<String> titles = new ArrayList<>();
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_residuals");
		outPrefixes.add("r5652_residuals_ask2014");
		stdDevs.add(false);
		doIndvs.add(true);
		titles.add(" ");
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_residuals");
		outPrefixes.add("r5652_std_devs_ask2014");
		stdDevs.add(true);
		doIndvs.add(false);
		titles.add(" ");
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_rake_ss_residuals");
		outPrefixes.add("r5652_ss_residuals_ask2014");
		titles.add("Strike-Slip");
		stdDevs.add(false);
		doIndvs.add(false);
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_rake_reverse_residuals");
		outPrefixes.add("r5652_reverse_residuals_ask2014");
		titles.add("Reverse");
		stdDevs.add(false);
		doIndvs.add(false);
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_rake_normal_residuals");
		outPrefixes.add("r5652_normal_residuals_ask2014");
		titles.add("Normal");
		stdDevs.add(false);
		doIndvs.add(false);
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_ASK2014_GriddedSites/resources"));
		prefixes.add("detrend_rake_oblique_residuals");
		outPrefixes.add("r5652_oblique_residuals_ask2014");
		titles.add("Oblique");
		stdDevs.add(false);
		doIndvs.add(false);
		
		refDirs.add(new File("/home/kevin/markdown/rsqsim-analysis/catalogs/rundir5652/bbp_LA_BASIN_500/"
				+ "gmpe_bbp_comparisons_NGAWest_2014_NoIdr_GriddedSites/resources"));
		prefixes.add("detrend_residuals");
		outPrefixes.add("r5652_residuals_ngaW2");
		stdDevs.add(false);
		doIndvs.add(false);
		titles.add(" ");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_Bruce_GMM_Multifault/figures/checkerboards");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double[] periods = {-1d, 2, 3, 5, 10};
		
		CPT residualCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1.5, 1.5);
		residualCPT.setNanColor(Color.WHITE);
		CPT stdDevCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0, 1d);
		stdDevCPT.setNanColor(Color.WHITE);
		CPT magCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(6.5d, 8d);
		magCPT.setBelowMinColor(magCPT.getMinColor());
		magCPT.setAboveMaxColor(magCPT.getMaxColor());
		
		Range scatterYRange = new Range(-1.5, 1.5);
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().setPlotLabelFontSize(36);
		gp.getPlotPrefs().setAxisLabelFontSize(28);
		gp.getPlotPrefs().setTickLabelFontSize(26);
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 36);
		
		PaintScaleLegend magSubtitle = GraphPanel.getLegendForCPT(magCPT, "Magnitude", gp.getPlotPrefs().getAxisLabelFontSize(),
				gp.getPlotPrefs().getTickLabelFontSize(), 0.5, RectangleEdge.LEFT);
		
		for (int n=0; n<refDirs.size(); n++) {
			File refDir = refDirs.get(n);
			String prefix = prefixes.get(n);
			String outPrefix = outPrefixes.get(n);
			String title = titles.get(n);
			boolean stdDev = stdDevs.get(n);
			boolean doIndv = doIndvs.get(n);
			
			CPT cpt = stdDev ? stdDevCPT : residualCPT;
			
			List<PlotSpec> specs = new ArrayList<>();
			Range xRange = null;
			Range yRange = null;
			
			List<PlotSpec> scatterSpecs = new ArrayList<>();
			
			for (int p=0; p<periods.length; p++) {
				String perPrefix;
				String perLabel;
				if (periods[p] == -1) {
					perPrefix = "pgv";
					perLabel = "PGV";
				} else {
					Preconditions.checkState(periods[p] > 0d);
					perPrefix = oDF.format(periods[p])+"s";
					perLabel = oDF.format(periods[p])+"s SA";
				}
				File csvFile = new File(refDir, prefix+"_"+perPrefix+".csv");
				CSVFile<String> csv = CSVFile.readFile(csvFile, true);
				
				Map<Double, XY_DataSet> magScatters = new HashMap<>();
				
				ArbDiscrXYZ_DataSet xyz = new ArbDiscrXYZ_DataSet();
				
				double minDist = Double.POSITIVE_INFINITY;
				double maxDist = Double.NEGATIVE_INFINITY;
				double minMag = Double.POSITIVE_INFINITY;
				double maxMag = Double.NEGATIVE_INFINITY;
				
				Map<Double, Double> distAverages = new HashMap<>();
				Map<Double, Integer> distCounts = new HashMap<>();
				
				for (int row=1; row<csv.getNumRows(); row++) {
					double distMin = csv.getDouble(row, 0);
					double distMax = csv.getDouble(row, 1);
					double magMin = csv.getDouble(row, 2);
					double magMax = csv.getDouble(row, 3);
					
					minDist = Math.min(minDist, distMin);
					maxDist = Math.max(maxDist, distMax);
					minMag = Math.min(minMag, magMin);
					maxMag = Math.max(maxMag, magMax);
					
					int count = csv.getInt(row, 4);
					
					double val;
					if (count == 0)
						val = Double.NaN;
					else if (stdDev && count == 1)
						val = Double.NaN;
					else if (stdDev)
						val = csv.getDouble(row, 10);
					else
						val = csv.getDouble(row, 5);
					
					double x = 0.5*(Math.log10(distMin) + Math.log10(distMax));
					double y = 0.5*(magMin + magMax);
					xyz.set(x, y, val);
					if (Double.isFinite(val) && !stdDev) {
						double mag = y;
						double dist = Math.pow(10, x);
						double resid = val;
						XY_DataSet scatter = magScatters.get(mag);
						if (scatter == null) {
							scatter = new DefaultXY_DataSet();
							magScatters.put(mag, scatter);
						}
						scatter.set(dist, resid);
						
						Double distAvg = distAverages.get(dist);
						Integer distCount = distCounts.get(dist);
						if (distAvg == null) {
							distAvg = 0d;
							distCount = 0;
						}
						distAverages.put(dist, distAvg + count*val);
						distCounts.put(dist, distCount + count);
					}
				}
				
				if (xRange == null) {
//					xRange = new Range(Math.log10(minDist), Math.log10(maxDist));
					xRange = new Range(minDist, maxDist);
					yRange = new Range(minMag, maxMag);
				}
				
				String xLabel = "Distance Rup  (km)";
				String zLabel = stdDev ? "Residual Std. Dev." : "Mean Residual  (ln)";
				
				XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, title, xLabel, "Magnitude", zLabel);
				spec.setCPTPosition(RectangleEdge.LEFT);
				spec.setCPTTickUnit(stdDev ? 0.2 : 0.5);
				
				XYTextAnnotation perAnn = new XYTextAnnotation("  "+perLabel, xRange.getLowerBound(),
						yRange.getLowerBound() + 0.98*yRange.getLength());
				perAnn.setFont(annFont);
				perAnn.setTextAnchor(TextAnchor.TOP_LEFT);
				spec.addPlotAnnotation(perAnn);
				
				specs.add(spec);
				
				List<XY_DataSet> scatterFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> scatterChars = new ArrayList<>();
				List<Double> mags = new ArrayList<>(magScatters.keySet());
				Collections.sort(mags);
				DefaultXY_DataSet scatterZeroXY = new DefaultXY_DataSet();
				scatterZeroXY.set(xRange.getLowerBound(), 0d);
				scatterZeroXY.set(xRange.getUpperBound(), 0d);
				scatterFuncs.add(scatterZeroXY);
				scatterChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.GRAY));
				
				// average on bottom
				XY_DataSet distAvgFunc = new DefaultXY_DataSet();
				for (Double dist : distAverages.keySet()) {
					double avg = distAverages.get(dist);
					int count = distCounts.get(dist);
					avg /= count;
					distAvgFunc.set(dist, avg);
				}
				scatterFuncs.add(distAvgFunc);
				scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 8f, Color.DARK_GRAY));
				
				for (double mag : mags) {
					XY_DataSet scatter = magScatters.get(mag);
					Color color = magCPT.getColor((float)mag);
					scatterFuncs.add(scatter);
					scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, color));
				}
				PlotSpec scatterSpec = new PlotSpec(scatterFuncs, scatterChars, title, xLabel, zLabel);
				scatterSpec.addSubtitle(magSubtitle);
				scatterSpecs.add(scatterSpec);
				
				if (doIndv) {
					gp.drawGraphPanel(spec, true, false, xRange, yRange);
					PlotUtils.setYTick(gp, 0.2d);
					
					PlotUtils.writePlots(outputDir, outPrefix+"_"+perPrefix, gp, 800, 600, true, true, false);
					
					if (!stdDev) {
						gp.drawGraphPanel(scatterSpec, true, false, xRange, scatterYRange);
						PlotUtils.setYTick(gp, 0.2d);
						
						PlotUtils.writePlots(outputDir, outPrefix+"_"+perPrefix+"_scatter", gp, 800, 600, true, true, false);
					}
				}
				if (p > 0)
					scatterSpec.setSubtitles(null);
			}
			
			int width = 400 + 600*specs.size();
			List<Range> xRanges = new ArrayList<>();
			for (int i=0; i<specs.size(); i++)
				xRanges.add(xRange);
			List<Range> yRanges = List.of(yRange);
			gp.drawGraphPanel(specs, true, false, xRanges, yRanges);
//			PlotUtils.setXTick(gp, 0.25d);
			PlotUtils.setYTick(gp, 0.2d);
			PlotUtils.setSubplotGap(gp, 80);
			
			PlotUtils.writePlots(outputDir, outPrefix+"_combined", gp, width, 600, true, true, false);
			
			if (!stdDev) {
				yRanges = List.of(scatterYRange);
				gp.drawGraphPanel(scatterSpecs, true, false, xRanges, yRanges);
				PlotUtils.setYTick(gp, 0.5d);
				PlotUtils.setSubplotGap(gp, 80);
				
				PlotUtils.writePlots(outputDir, outPrefix+"_combined_scatter", gp, width, 600, true, true, false);
			}
		}
	}

}
