package scratch.kevin.sampling;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;

import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceMetric;
import scratch.kevin.sampling.HazardConvergenceCalcs.HazardStatistics;
import scratch.kevin.sampling.HazardConvergenceCalcs.ModelHazardMaps;

public class HazardMapPlots {

	private static final ReturnPeriods RP = ReturnPeriods.TWO_IN_50;

	public static void main(String[] args) throws IOException {
		File convergenceDir = new File(PaperPaths.FIGURES_DIR, "hazard_convergence");
		int[] sizes = { 512, 1024, 2048, 4096, 8192 };
		plotPeriod(new File(convergenceDir, "pga_two_in_50"), 0d, "PGA, "+RP.label, sizes);
		plotPeriod(new File(convergenceDir, "1s_sa_two_in_50"), 1d, "1s SA, "+RP.label, sizes);
	}

	private static void plotPeriod(File periodDir, double period, String perLabel, int... indvSizes)
			throws IOException {
		File mapDir = new File(periodDir, "hazard_maps");
		Preconditions.checkState(mapDir.exists() || mapDir.mkdir());

		File firstMCSDir = HazardConvergenceCalcs.runDirs.row(SamplingMethod.MONTE_CARLO)
				.values().iterator().next().get(0);
		GriddedRegion gridReg = GriddedRegion.fromFeature(
				Feature.read(new File(firstMCSDir, "gridded_region.geojson")));

		List<File> mcsPoolDirs = flatten(HazardConvergenceCalcs.runDirs.row(SamplingMethod.MONTE_CARLO));
		List<File> sobolPoolDirs;
		if (HazardConvergenceCalcs.FIXED_SOBOL_CONSENSUS_SIZE == null) {
			sobolPoolDirs = flatten(HazardConvergenceCalcs.runDirs.row(SamplingMethod.OWEN_SCRAMBLED_SOBOL));
		} else {
			sobolPoolDirs = Preconditions.checkNotNull(HazardConvergenceCalcs.runDirs.get(
					SamplingMethod.OWEN_SCRAMBLED_SOBOL,
					HazardConvergenceCalcs.FIXED_SOBOL_CONSENSUS_SIZE));
		}
		
		DecimalFormat groupedDF = new DecimalFormat("0");
		groupedDF.setGroupingSize(3);
		groupedDF.setGroupingUsed(true);

		PooledHazardData mcsPool = loadPool(periodDir, "pooled_mcs", mcsPoolDirs, gridReg, period);
		PooledHazardData sobolPool = loadPool(periodDir, "pooled_sobol", sobolPoolDirs, gridReg, period);
		List<File> poLHSPoolDirs = largestRuns(HazardConvergenceCalcs.runDirs.row(
				SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE));
		PooledHazardData poLHSPool = poLHSPoolDirs.isEmpty() ? null
				: loadPool(periodDir, "pooled_po_lhs", poLHSPoolDirs, gridReg, period);

		plotMeanHazard(gridReg, sobolPool.data(), perLabel, mapDir, "pooled_sobol");
		plotMeanHazard(gridReg, mcsPool.data(), perLabel, mapDir, "pooled_mcs");
		plotComparisons(gridReg, sobolPool.data(), mcsPool.data(), perLabel, mapDir, "pooled_sobol_vs_mcs",
				"Pooled Sobol' ("+nStr(sobolPool.data)+") vs MCS ("+nStr(mcsPool.data)+")");
		if (poLHSPool != null) {
			plotMeanHazard(gridReg, poLHSPool.data(), perLabel, mapDir, "pooled_po_lhs");
			plotComparisons(gridReg, poLHSPool.data(), mcsPool.data(), perLabel, mapDir,
					"pooled_po_lhs_vs_mcs", "Pooled Pairwise-Optimized LHS ("+nStr(poLHSPool.data)
					+") vs MCS ("+nStr(mcsPool.data)+")");
			plotComparisons(gridReg, poLHSPool.data(), sobolPool.data(), perLabel, mapDir,
					"pooled_po_lhs_vs_sobol", "Pooled Pairwise-Optimized LHS ("+nStr(poLHSPool.data)
					+") vs Sobol' ("+nStr(sobolPool.data)+")");
		}

		// Plot the first realization available for each method and requested sample count.
		for (int size : indvSizes) {
			File sizeDir = new File(mapDir, size+"_samples");
			Preconditions.checkState(sizeDir.exists() || sizeDir.mkdir());
			for (SamplingMethod method : SamplingMethod.values()) {
				HazardData data;
				HazardData refMCS;
				HazardData refSobol;
				if (method == SamplingMethod.MONTE_CARLO) {
					// MCS comparisons use the first span of the first run and remove that span from the MCS pool.
					data = loadRunPrefix(firstMCSDir, size, gridReg, period, mcsPool);
					refMCS = mcsPool.without(firstMCSDir, size, data.meanCurves(), gridReg);
					refSobol = sobolPool.data();
				} else {
					List<File> runDirs = HazardConvergenceCalcs.runDirs.get(method, size);
					if (runDirs == null || runDirs.isEmpty())
						continue;
					File runDir = runDirs.get(0);
					data = loadRunPrefix(runDir, size, gridReg, period,
							method == SamplingMethod.OWEN_SCRAMBLED_SOBOL ? sobolPool : null);
					refMCS = mcsPool.data();
					if (method == SamplingMethod.OWEN_SCRAMBLED_SOBOL && sobolPool.contains(runDir))
						refSobol = sobolPool.without(runDir, size, data.meanCurves(), gridReg);
					else
						refSobol = sobolPool.data();
				}
				
				String name = "Individual "+HazardConvergencePlots.getMethodName(method)+" ("+nStr(data)+")";
				
				plotComparisons(gridReg, data, refMCS, perLabel, sizeDir,
						method.name().toLowerCase()+"_vs_pooled_mcs", name+" vs pooled MCS ("+nStr(refMCS)+")");
				plotComparisons(gridReg, data, refSobol, perLabel, sizeDir,
						method.name().toLowerCase()+"_vs_pooled_sobol", name+" vs pooled Sobol' ("+nStr(refSobol)+")");
			}
		}
	}
	
	private static final DecimalFormat groupedDF = new DecimalFormat("0");
	static {
		groupedDF.setGroupingSize(3);
		groupedDF.setGroupingUsed(true);
	}
	
	private static String nStr(HazardData data) {
		return nStr(data.numBranches());
	}
	
	private static String nStr(int size) {
		return "N="+groupedDF.format(size);
	}

	private static List<File> flatten(Map<Integer, List<File>> dirsBySize) {
		List<File> dirs = new ArrayList<>();
		dirsBySize.entrySet().stream().sorted(Map.Entry.comparingByKey())
				.forEach(entry -> dirs.addAll(entry.getValue()));
		return dirs;
	}

	private static List<File> largestRuns(Map<Integer, List<File>> dirsBySize) {
		return dirsBySize.entrySet().stream().max(Map.Entry.comparingByKey())
				.map(Map.Entry::getValue).orElse(List.of());
	}

	private static PooledHazardData loadPool(File periodDir, String poolName, List<File> runDirs,
			GriddedRegion gridReg, double period) throws IOException {
		Preconditions.checkState(!runDirs.isEmpty(), "No runs available for %s", poolName);
		DiscretizedFunc[] meanCurves = loadCurves(
				new File(new File(periodDir, poolName),
						SolHazardMapCalc.getCSV_FileName("mean_curves", period)+".gz"), gridReg);
		List<RunBlock> blocks = new ArrayList<>();
		for (File runDir : runDirs) {
			LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree_analysis.json"));
			ModelHazardMaps maps = HazardConvergenceCalcs.loadMaps(
					new File(runDir, "results_hazard.zip"), tree, gridReg, period, RP);
			blocks.add(new RunBlock(runDir.getAbsoluteFile(), HazardConvergenceCalcs.copyValues(maps.individual())));
		}
		return new PooledHazardData(buildHazardData(meanCurves, concatenate(blocks), gridReg), blocks);
	}

	private static HazardData loadRunPrefix(File runDir, int sampleCount, GriddedRegion gridReg,
			double period, PooledHazardData cachedPool) throws IOException {
		LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree_analysis.json"));
		Preconditions.checkState(sampleCount <= tree.size(), "Requested %s of %s branches from %s",
				sampleCount, tree.size(), runDir.getName());
		double[][] allBranchMaps = cachedPool == null ? null : cachedPool.branchMaps(runDir);
		if (allBranchMaps == null) {
			ModelHazardMaps maps = HazardConvergenceCalcs.loadMaps(
					new File(runDir, "results_hazard.zip"), tree, gridReg, period, RP);
			allBranchMaps = HazardConvergenceCalcs.copyValues(maps.individual());
		}
		double[][] branchMaps = Arrays.copyOf(allBranchMaps, sampleCount);
		DiscretizedFunc[] meanCurves = loadMeanCurves(runDir, tree, sampleCount, gridReg, period);
		return buildHazardData(meanCurves, branchMaps, gridReg);
	}

	private static DiscretizedFunc[] loadMeanCurves(File runDir, LogicTree<?> tree, int sampleCount,
			GriddedRegion gridReg, double period) throws IOException {
		double[] xValues = null;
		double[][] sums = null;
		try (HazardConvergenceCalcs.BranchCurveLoader curveLoader =
				new HazardConvergenceCalcs.BranchCurveLoader(runDir)) {
			for (int b=0; b<sampleCount; b++) {
				DiscretizedFunc[] curves = curveLoader.load(tree.getBranch(b), gridReg, period);
				if (sums == null) {
					xValues = new double[curves[0].size()];
					for (int i=0; i<xValues.length; i++)
						xValues[i] = curves[0].getX(i);
					sums = new double[curves.length][xValues.length];
				}
				Preconditions.checkState(curves.length == sums.length);
				for (int n=0; n<curves.length; n++) {
					Preconditions.checkState(curves[n].size() == xValues.length);
					for (int i=0; i<xValues.length; i++) {
						Preconditions.checkState((float)curves[n].getX(i) == (float)xValues[i]);
						sums[n][i] += curves[n].getY(i);
					}
				}
			}
		}
		return scaledCurves(sums, xValues, 1d/sampleCount);
	}

	private static DiscretizedFunc[] loadCurves(File file, GriddedRegion gridReg) throws IOException {
		Preconditions.checkState(file.isFile(), "Hazard curves file doesn't exist: %s", file.getAbsolutePath());
		return SolHazardMapCalc.loadCurvesCSV(CSVFile.readFile(file, true), gridReg);
	}

	private static HazardData buildHazardData(DiscretizedFunc[] meanCurves, double[][] branchMaps,
			GriddedRegion gridReg) {
		GriddedGeoDataSet meanMap = SolHazardMapCalc.buildMap(meanCurves, gridReg, RP);
		double[] meanValues = new double[meanMap.size()];
		for (int i=0; i<meanValues.length; i++)
			meanValues[i] = meanMap.get(i);
		HazardStatistics statistics = HazardConvergenceCalcs.calcHazardStatistics(
				branchMaps, branchMaps.length, meanValues);
		return new HazardData(meanCurves, branchMaps, statistics);
	}

	private static DiscretizedFunc[] subtractMeanCurves(DiscretizedFunc[] fullMean, int fullCount,
			DiscretizedFunc[] excludedMean, int excludedCount) {
		Preconditions.checkState(excludedCount < fullCount);
		Preconditions.checkState(fullMean.length == excludedMean.length);
		DiscretizedFunc[] curves = new DiscretizedFunc[fullMean.length];
		for (int n=0; n<curves.length; n++) {
			Preconditions.checkState(fullMean[n].size() == excludedMean[n].size());
			double[] xValues = new double[fullMean[n].size()];
			double[] yValues = new double[xValues.length];
			for (int i=0; i<xValues.length; i++) {
				xValues[i] = fullMean[n].getX(i);
				Preconditions.checkState((float)xValues[i] == (float)excludedMean[n].getX(i));
				yValues[i] = (fullCount*fullMean[n].getY(i)-excludedCount*excludedMean[n].getY(i))
						/(fullCount-excludedCount);
			}
			curves[n] = new LightFixedXFunc(xValues, yValues);
		}
		return curves;
	}

	private static DiscretizedFunc[] scaledCurves(double[][] values, double[] xValues, double scale) {
		DiscretizedFunc[] curves = new DiscretizedFunc[values.length];
		for (int n=0; n<curves.length; n++) {
			double[] yValues = values[n].clone();
			for (int i=0; i<yValues.length; i++)
				yValues[i] *= scale;
			curves[n] = new LightFixedXFunc(xValues, yValues);
		}
		return curves;
	}

	private static double[][] concatenate(List<RunBlock> blocks) {
		int size = blocks.stream().mapToInt(block -> block.branchMaps().length).sum();
		double[][] values = new double[size][];
		int index = 0;
		for (RunBlock block : blocks)
			for (double[] row : block.branchMaps())
				values[index++] = row;
		return values;
	}
	
	private static Region getMapReagion(GriddedRegion gridReg) {
		if (!gridReg.isRectangular())
			return gridReg;
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (Location loc : gridReg.getNodeList()) {
			latTrack.addValue(loc.lat);
			lonTrack.addValue(loc.lon);
		}
		double halfLat = gridReg.getLatSpacing()*0.5;
		double halfLon = gridReg.getLonSpacing()*0.5;
		return new Region(new Location(latTrack.getMin()-halfLat, lonTrack.getMin()-halfLon),
				new Location(latTrack.getMax()+halfLat, lonTrack.getMax()+halfLon));
	}

	private static void plotMeanHazard(GriddedRegion gridReg, HazardData data,
			String perLabel, File outputDir, String prefix) throws IOException {
		GeographicMapMaker mapMaker = new GeographicMapMaker(getMapReagion(gridReg));
		CPT hazCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(1e-2, 3d).asLog10();
		GriddedGeoDataSet map = SolHazardMapCalc.buildMap(data.meanCurves(), gridReg, RP);
		mapMaker.plotXYZData(map, hazCPT, "Mean hazard, "+perLabel+" (g)");
		mapMaker.plot(outputDir, prefix, "", PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 300);
	}

	private static void plotComparisons(GriddedRegion gridReg, HazardData data, HazardData reference,
			String perLabel, File outputDir, String prefix, String title) throws IOException {
		Region region = getMapReagion(gridReg);
		GeographicMapMaker mapMaker = new GeographicMapMaker(region);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-5d, 5d);
		pDiffCPT.setPreferredTickInterval(1d);
		
		DecimalFormat df = new DecimalFormat("0.00");
		Font statsFont = new Font(Font.SANS_SERIF, Font.PLAIN, 8);
		Font metricFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
//		Color bgPaint = new Color(255, 255, 255, 60);
//		Color bgPaint = new Color(200, 200, 200, 120);
		Color bgPaint = new Color(220, 220, 220, 200);
		
		List<PlotSpec> plots = new ArrayList<>();
		
		for (ConvergenceMetric metric : HazardConvergencePlots.PLOT_METRICS) {
			GriddedGeoDataSet xyz = asGeoDataSet(gridReg, data.statistics().values(metric));
			GriddedGeoDataSet refXYZ = asGeoDataSet(gridReg, reference.statistics().values(metric));
			GriddedGeoDataSet pDiff = pDiff(xyz, refXYZ);
			mapMaker.plotXYZData(pDiff, pDiffCPT, perLabel+", % change");

			mapMaker.clearAnnotations();
			
			double mean = 0d;
			double meanAbs = 0d;
			double min=Double.MAX_VALUE;
			double max=Double.MIN_VALUE;
			for (int i=0; i<pDiff.size(); i++) {
				double v = pDiff.get(i);
				mean += v;
				meanAbs += Math.abs(v);
				min = Math.min(min, v);
				max = Math.max(max, v);
			}
			mean /= xyz.size();
			meanAbs /= xyz.size();
			String label = "mean="+df.format(mean)+"%; meanAbs="+df.format(meanAbs)
					+"%; range=["+df.format(min)+", "+df.format(max)+"]%";
			double annX = 0.5*(region.getMinLon() + region.getMaxLon());
			double annY = region.getMinLat()+0.25;
			XYTextAnnotation ann = new XYTextAnnotation(" "+label+" ", annX, annY);
			ann.setBackgroundPaint(bgPaint);
			ann.setFont(statsFont);
			ann.setTextAnchor(TextAnchor.CENTER);
			mapMaker.addAnnotation(ann);
			
			label = metric.shortLabel;
			annX = region.getMaxLon()-0.1;
			annY = region.getMaxLat()-0.1;
			ann = new XYTextAnnotation(" "+label+" ", annX, annY);
			ann.setBackgroundPaint(bgPaint);
			ann.setFont(metricFont);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			mapMaker.addAnnotation(ann);
			
			plots.add(mapMaker.buildPlot(title, true));
			
//			mapMaker.plot(outputDir, metric.name()+"_"+prefix, "",
//					PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 300);
		}
		Range xRange = mapMaker.getXRange();
		Range yRange = mapMaker.getYRange();
		
		if (plots.size() == 1) {
			HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
			
			gp.drawGraphPanel(plots.getFirst(), false, false, xRange, yRange);
			
			PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, true, 300, true, true, false);
		} else if (plots.size() == 2) {
			HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
			
			gp.drawGraphPanel(plots, false, false, List.of(xRange, xRange), List.of(yRange));
			
			PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH, true, 300, true, true, false);
		} else if (plots.size() % 2 != 0) {
			// vertical
			List<Range> yRanges = new ArrayList<>(plots.size());
			for (int i=0; i<plots.size(); i++)
				yRanges.add(yRange);
			
			HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
			
			gp.drawGraphPanel(plots, false, false, List.of(xRange), yRanges);
			
			PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, true, 300, true, true, false);
		} else {
			// arrange into rows
			List<List<PlotSpec>> rows = new ArrayList<>();
			List<PlotSpec> curRow = null;
			for (int i=0; i<plots.size(); i++) {
				if (i % 2 == 0) {
					curRow = new ArrayList<>();
					rows.add(curRow);
				}
				curRow.add(plots.get(i));
			}
			List<HeadlessGraphPanel> gps = new ArrayList<>(rows.size());
			for (List<PlotSpec> row : rows) {
				HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
				gp.drawGraphPanel(row, false, false, List.of(xRange, xRange), List.of(yRange));
				gps.add(gp);
			}
			
			PlotUtils.stitchPlotRows(outputDir, prefix, gps, false, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH, -1d, 300, true, true, true);
		}
//		List<Range> yRanges = new ArrayList<>(plots.size());
//		for (int i=0; i<plots.size(); i++)
//			yRanges.add(yRange);
//		
//		HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
//		
//		gp.drawGraphPanel(plots, false, false, List.of(xRange), yRanges);
//		
//		PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, true, 300, true, true, false);
	}

	private static GriddedGeoDataSet asGeoDataSet(GriddedRegion gridReg, double[] values) {
		Preconditions.checkState(values.length == gridReg.getNodeCount());
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<values.length; i++)
			xyz.set(i, values[i]);
		return xyz;
	}

	private static GriddedGeoDataSet pDiff(GriddedGeoDataSet xyz, GriddedGeoDataSet refXYZ) {
		Preconditions.checkState(refXYZ.size() == xyz.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(xyz.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double val = xyz.get(i);
			double refVal = refXYZ.get(i);
			if (!Double.isFinite(val) || !Double.isFinite(refVal) || refVal == 0d)
				ret.set(i, Double.NaN);
			else
				ret.set(i, 100d*(val-refVal)/refVal);
		}
		return ret;
	}

	private record HazardData(DiscretizedFunc[] meanCurves, double[][] branchMaps,
			HazardStatistics statistics) {
		
		public int numBranches() {
			return branchMaps.length;
		}
	}

	private record RunBlock(File directory, double[][] branchMaps) {
		boolean matches(File other) {
			return directory.equals(other.getAbsoluteFile());
		}
	}

	private record PooledHazardData(HazardData data, List<RunBlock> blocks) {
		boolean contains(File runDir) {
			return blocks.stream().anyMatch(block -> block.matches(runDir));
		}

		double[][] branchMaps(File runDir) {
			return blocks.stream().filter(block -> block.matches(runDir)).findFirst()
					.map(RunBlock::branchMaps).orElse(null);
		}

		HazardData without(File runDir, int excludedCount, DiscretizedFunc[] excludedMean,
				GriddedRegion gridReg) {
			Preconditions.checkState(contains(runDir), "Run is not in pooled data: %s", runDir.getName());
			List<RunBlock> retained = new ArrayList<>();
			for (RunBlock block : blocks) {
				if (block.matches(runDir)) {
					Preconditions.checkState(excludedCount <= block.branchMaps().length);
					if (excludedCount < block.branchMaps().length)
						retained.add(new RunBlock(block.directory(),
								Arrays.copyOfRange(block.branchMaps(), excludedCount, block.branchMaps().length)));
				} else {
					retained.add(block);
				}
			}
			int fullCount = data.branchMaps().length;
			DiscretizedFunc[] meanCurves = subtractMeanCurves(
					data.meanCurves(), fullCount, excludedMean, excludedCount);
			return buildHazardData(meanCurves, concatenate(retained), gridReg);
		}
	}
}
