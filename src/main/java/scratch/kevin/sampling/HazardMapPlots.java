package scratch.kevin.sampling;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
		int[] sizes = { 512, 1024, 2048, 4096, 8192, 16384 };
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
		File mcsReserveDir = HazardConvergenceCalcs.MCS_RESERVE_DIR;
		LogicTree<?> mcsReserveTree = LogicTree.read(new File(mcsReserveDir, "logic_tree_analysis.json"));
		Preconditions.checkState(mcsReserveTree.getSamplingMethod() == SamplingMethod.MONTE_CARLO,
				"Expected MCS reserve, have %s", mcsReserveTree.getSamplingMethod());
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
		Map<Integer, PooledHazardData> exactSobolPools = new TreeMap<>();
		for (int sampleCount : HazardConvergenceCalcs.SOBOL_POOL_OUTPUT_SIZES) {
			List<File> exactSizeDirs = HazardConvergenceCalcs.runDirs.get(
					SamplingMethod.OWEN_SCRAMBLED_SOBOL, sampleCount);
			if (exactSizeDirs == null || exactSizeDirs.isEmpty())
				continue;
			PooledHazardData pool = HazardConvergenceCalcs.FIXED_SOBOL_CONSENSUS_SIZE != null
					&& sampleCount == HazardConvergenceCalcs.FIXED_SOBOL_CONSENSUS_SIZE
					&& exactSizeDirs.equals(sobolPoolDirs) ? sobolPool
							: loadPool(periodDir, HazardConvergenceCalcs.sobolPoolDirName(sampleCount),
									exactSizeDirs, gridReg, period);
			exactSobolPools.put(sampleCount, pool);
		}
		List<File> poLHSPoolDirs = largestRuns(HazardConvergenceCalcs.runDirs.row(
				SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE));
		PooledHazardData poLHSPool = poLHSPoolDirs.isEmpty() ? null
				: loadPool(periodDir, "pooled_po_lhs", poLHSPoolDirs, gridReg, period);

		plotMeanHazard(gridReg, sobolPool.data(), perLabel, mapDir, "pooled_sobol");
		plotMeanHazard(gridReg, mcsPool.data(), perLabel, mapDir, "pooled_mcs");
		plotComparisons(gridReg, sobolPool.data(), mcsPool.data(), perLabel, mapDir, "pooled_sobol_vs_mcs",
				"Pooled Sobol' ("+nStr(sobolPool.data)+") vs MCS ("+nStr(mcsPool.data)+")");
		PooledHazardData sobol8192 = exactSobolPools.get(8192);
		PooledHazardData sobol16384 = exactSobolPools.get(16384);
		if (sobol8192 != null && sobol16384 != null)
			plotComparisons(gridReg, sobol16384.data(), sobol8192.data(), perLabel, mapDir,
					"pooled_sobol_16384_vs_8192", "Pooled 16,384-sample Sobol' ("+nStr(sobol16384.data)
					+") vs 8,192-sample Sobol' ("+nStr(sobol8192.data)+")");
		if (poLHSPool != null) {
			plotMeanHazard(gridReg, poLHSPool.data(), perLabel, mapDir, "pooled_po_lhs");
			plotComparisons(gridReg, poLHSPool.data(), mcsPool.data(), perLabel, mapDir,
					"pooled_po_lhs_vs_mcs", "Pooled Pairwise-Optimized LHS ("+nStr(poLHSPool.data)
					+") vs MCS ("+nStr(mcsPool.data)+")");
			plotComparisons(gridReg, poLHSPool.data(), sobolPool.data(), perLabel, mapDir,
					"pooled_po_lhs_vs_sobol", "Pooled Pairwise-Optimized LHS ("+nStr(poLHSPool.data)
					+") vs Sobol' ("+nStr(sobolPool.data)+")");
		}

		CSVFile<String> selectionCSV = initSelectionCSV();
		CSVFile<String> referenceComparisons = CSVFile.readFile(
				new File(periodDir, "reference_comparisons.csv"), true);
		Map<File, HazardData> fullSobolRunCache = new HashMap<>();
		Map<File, double[][]> runMapCache = new HashMap<>();
		for (int size : indvSizes) {
			File sizeDir = new File(mapDir, size+"_samples");
			Preconditions.checkState(sizeDir.exists() || sizeDir.mkdir());
			for (SamplingMethod method : SamplingMethod.values()) {
				HazardData data;
				HazardData refMCS;
				HazardData refSobol;
				RepresentativeSelection selection;
				if (method == SamplingMethod.MONTE_CARLO) {
					if (size > mcsReserveTree.size())
						continue;
					// The reserve is independent of, and deliberately excluded from, the primary MCS pool. Use
					// only its nonoverlapping spans as candidates; production-pool spans are not considered.
					List<RepresentativeCandidate> candidates = new ArrayList<>();
					for (int start=0; start+size<=mcsReserveTree.size(); start+=size) {
						HazardData candidateData = loadRunSpan(mcsReserveDir, start, size, gridReg, period,
								null, runMapCache);
						candidates.add(new RepresentativeCandidate(mcsReserveDir, start, start+size,
								candidateData, meanAbsErrors(candidateData, mcsPool.data())));
					}
					selection = selectRepresentative(candidates);
					data = selection.candidate().data();
					refMCS = mcsPool.data();
					refSobol = sobolPool.data();
				} else {
					List<RepresentativeCandidate> candidates = loadCandidateErrors(
							referenceComparisons, method, size);
					if (candidates.isEmpty())
						continue;
					selection = selectRepresentative(candidates);
					RepresentativeCandidate selected = selection.candidate();
					File runDir = selected.runDir();
					data = loadRunSpan(runDir, selected.startIndex(), size, gridReg, period,
							method == SamplingMethod.OWEN_SCRAMBLED_SOBOL ? sobolPool : null, runMapCache);
					selection = selection.withData(data);
					refMCS = mcsPool.data();
					if (method == SamplingMethod.OWEN_SCRAMBLED_SOBOL && sobolPool.contains(runDir)) {
						int fullRunSize = sobolPool.runSize(runDir);
						HazardData fullRun = data;
						if (size < fullRunSize) {
							fullRun = fullSobolRunCache.get(runDir.getAbsoluteFile());
							if (fullRun == null) {
								fullRun = loadRunSpan(runDir, 0, fullRunSize, gridReg, period,
										sobolPool, runMapCache);
								fullSobolRunCache.put(runDir.getAbsoluteFile(), fullRun);
							}
						}
						// Prefixes from one scramble are not independent of its suffix, so exclude the entire parent run.
						refSobol = sobolPool.without(runDir, fullRunSize, fullRun.meanCurves(), gridReg);
					} else {
						refSobol = sobolPool.data();
					}
				}
				appendSelection(selectionCSV, method, size, selection);
				System.out.println("Representative "+method.getShortName()+", N="+size+": "
						+selection.candidate().description()+", candidates="+selection.candidateCount()
						+", mean rank distance="+selection.meanRankDistance());
				
				String name = "Individual "+HazardConvergencePlots.getMethodName(method)+" ("+nStr(data)+")";
				
				plotComparisons(gridReg, data, refMCS, perLabel, sizeDir,
						method.name().toLowerCase()+"_vs_pooled_mcs", name+" vs pooled MCS ("+nStr(refMCS)+")");
				plotComparisons(gridReg, data, refSobol, perLabel, sizeDir,
						method.name().toLowerCase()+"_vs_pooled_sobol", name+" vs pooled Sobol' ("+nStr(refSobol)+")");
			}
		}
		selectionCSV.writeToFile(new File(mapDir, "representative_selections.csv"));
	}

	private static List<RepresentativeCandidate> loadCandidateErrors(CSVFile<String> csv,
			SamplingMethod method, int sampleCount) {
		List<String> header = csv.getLine(0);
		int runCol = header.indexOf("Run");
		int methodCol = header.indexOf("Sampling method");
		int countCol = header.indexOf("Sample count");
		int referenceCol = header.indexOf("Reference");
		int metricCol = header.indexOf("Metric");
		int meanAbsCol = header.indexOf("Spatial mean absolute % change");
		int startCol = header.indexOf("Start index");
		int endCol = header.indexOf("End index (exclusive)");
		Preconditions.checkState(runCol >= 0 && methodCol >= 0 && countCol >= 0 && referenceCol >= 0
				&& metricCol >= 0 && meanAbsCol >= 0 && startCol >= 0 && endCol >= 0,
				"Unexpected reference comparison CSV header");

		Map<CandidateID, double[]> errorsByCandidate = new LinkedHashMap<>();
		for (int row=1; row<csv.getNumRows(); row++) {
			if (!csv.get(row, methodCol).equals(method.name())
					|| Integer.parseInt(csv.get(row, countCol)) != sampleCount
					|| !csv.get(row, referenceCol).equals(HazardConvergenceCalcs.MCS_REFERENCE_NAME))
				continue;
			int metricIndex = metricIndex(csv.get(row, metricCol));
			if (metricIndex < 0)
				continue;
			CandidateID id = new CandidateID(csv.get(row, runCol),
					Integer.parseInt(csv.get(row, startCol)), Integer.parseInt(csv.get(row, endCol)));
			double[] errors = errorsByCandidate.computeIfAbsent(id, key -> {
				double[] values = new double[HazardConvergencePlots.PLOT_METRICS.length];
				Arrays.fill(values, Double.NaN);
				return values;
			});
			Preconditions.checkState(Double.isNaN(errors[metricIndex]),
					"Duplicate comparison for %s, %s", id, csv.get(row, metricCol));
			errors[metricIndex] = Double.parseDouble(csv.get(row, meanAbsCol));
		}

		List<RepresentativeCandidate> candidates = new ArrayList<>();
		for (Map.Entry<CandidateID, double[]> entry : errorsByCandidate.entrySet()) {
			for (double value : entry.getValue())
				Preconditions.checkState(Double.isFinite(value), "Incomplete comparison for %s", entry.getKey());
			CandidateID id = entry.getKey();
			File runDir = findRunDir(method, id.runName());
			candidates.add(new RepresentativeCandidate(runDir, id.startIndex(), id.endIndex(),
					null, entry.getValue()));
		}
		return candidates;
	}

	private static int metricIndex(String label) {
		for (int i=0; i<HazardConvergencePlots.PLOT_METRICS.length; i++)
			if (HazardConvergencePlots.PLOT_METRICS[i].label.equals(label))
				return i;
		return -1;
	}

	private static File findRunDir(SamplingMethod method, String runName) {
		for (List<File> dirs : HazardConvergenceCalcs.runDirs.row(method).values())
			for (File dir : dirs)
				if (dir.getName().equals(runName))
					return dir;
		throw new IllegalStateException("No configured "+method+" run named "+runName);
	}

	private static double[] meanAbsErrors(HazardData data, HazardData reference) {
		double[] errors = new double[HazardConvergencePlots.PLOT_METRICS.length];
		for (int m=0; m<errors.length; m++) {
			double[] values = data.statistics().values(HazardConvergencePlots.PLOT_METRICS[m]);
			double[] refValues = reference.statistics().values(HazardConvergencePlots.PLOT_METRICS[m]);
			double sum = 0d;
			int count = 0;
			for (int i=0; i<values.length; i++) {
				if (Double.isFinite(values[i]) && Double.isFinite(refValues[i]) && refValues[i] != 0d) {
					sum += Math.abs(100d*(values[i]-refValues[i])/refValues[i]);
					count++;
				}
			}
			Preconditions.checkState(count > 0);
			errors[m] = sum/count;
		}
		return errors;
	}

	static RepresentativeSelection selectRepresentative(List<RepresentativeCandidate> candidates) {
		Preconditions.checkState(!candidates.isEmpty());
		int count = candidates.size();
		double[][] ranks = new double[count][HazardConvergencePlots.PLOT_METRICS.length];
		for (int m=0; m<HazardConvergencePlots.PLOT_METRICS.length; m++) {
			final int metricIndex = m;
			List<Integer> order = new ArrayList<>(count);
			for (int i=0; i<count; i++)
				order.add(i);
			order.sort(Comparator.comparingDouble(i -> candidates.get(i).errors()[metricIndex]));
			for (int start=0; start<count;) {
				int end = start+1;
				double value = candidates.get(order.get(start)).errors()[m];
				while (end < count && Double.compare(candidates.get(order.get(end)).errors()[m], value) == 0)
					end++;
				double rank = count == 1 ? 0.5 : (0.5*(start+end-1))/(count-1d);
				for (int i=start; i<end; i++)
					ranks[order.get(i)][m] = rank;
				start = end;
			}
		}

		int bestIndex = -1;
		double bestMean = Double.POSITIVE_INFINITY;
		double bestMax = Double.POSITIVE_INFINITY;
		for (int i=0; i<count; i++) {
			double sum = 0d;
			double max = 0d;
			for (double rank : ranks[i]) {
				double distance = Math.abs(rank-0.5);
				sum += distance;
				max = Math.max(max, distance);
			}
			double mean = sum/ranks[i].length;
			if (mean < bestMean || (mean == bestMean && max < bestMax)) {
				bestIndex = i;
				bestMean = mean;
				bestMax = max;
			}
		}
		return new RepresentativeSelection(candidates.get(bestIndex), count,
				ranks[bestIndex], bestMean, bestMax);
	}

	private static CSVFile<String> initSelectionCSV() {
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>(List.of("Sampling method", "Sample count", "Run",
				"Start index", "End index (exclusive)", "Candidate count",
				"Mean absolute quantile-rank distance", "Maximum absolute quantile-rank distance"));
		for (ConvergenceMetric metric : HazardConvergencePlots.PLOT_METRICS) {
			header.add(metric.label+" spatial mean absolute % change");
			header.add(metric.label+" quantile rank");
		}
		csv.addLine(header);
		return csv;
	}

	private static void appendSelection(CSVFile<String> csv, SamplingMethod method, int sampleCount,
			RepresentativeSelection selection) {
		RepresentativeCandidate candidate = selection.candidate();
		List<String> line = new ArrayList<>(List.of(method.name(), Integer.toString(sampleCount),
				candidate.runDir().getName(), Integer.toString(candidate.startIndex()),
				Integer.toString(candidate.endIndex()), Integer.toString(selection.candidateCount()),
				Double.toString(selection.meanRankDistance()), Double.toString(selection.maxRankDistance())));
		for (int m=0; m<candidate.errors().length; m++) {
			line.add(Double.toString(candidate.errors()[m]));
			line.add(Double.toString(selection.ranks()[m]));
		}
		csv.addLine(line);
	}

	private record CandidateID(String runName, int startIndex, int endIndex) {}

	record RepresentativeCandidate(File runDir, int startIndex, int endIndex,
			HazardData data, double[] errors) {
		String description() {
			return runDir.getName()+"["+startIndex+","+endIndex+")";
		}
	}

	record RepresentativeSelection(RepresentativeCandidate candidate, int candidateCount,
			double[] ranks, double meanRankDistance, double maxRankDistance) {
		RepresentativeSelection withData(HazardData data) {
			return new RepresentativeSelection(new RepresentativeCandidate(candidate.runDir(),
					candidate.startIndex(), candidate.endIndex(), data, candidate.errors()),
					candidateCount, ranks, meanRankDistance, maxRankDistance);
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

	private static HazardData loadRunSpan(File runDir, int startIndex, int sampleCount,
			GriddedRegion gridReg, double period, PooledHazardData cachedPool,
			Map<File, double[][]> runMapCache) throws IOException {
		LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree_analysis.json"));
		Preconditions.checkState(startIndex >= 0 && startIndex+sampleCount <= tree.size(),
				"Requested [%s,%s) of %s branches from %s",
				startIndex, startIndex+sampleCount, tree.size(), runDir.getName());
		double[][] allBranchMaps = cachedPool == null ? null : cachedPool.branchMaps(runDir);
		File cacheKey = runDir.getAbsoluteFile();
		if (allBranchMaps == null)
			allBranchMaps = runMapCache.get(cacheKey);
		if (allBranchMaps == null) {
			ModelHazardMaps maps = HazardConvergenceCalcs.loadMaps(
					new File(runDir, "results_hazard.zip"), tree, gridReg, period, RP);
			allBranchMaps = HazardConvergenceCalcs.copyValues(maps.individual());
			runMapCache.put(cacheKey, allBranchMaps);
		}
		double[][] branchMaps = Arrays.copyOfRange(allBranchMaps, startIndex, startIndex+sampleCount);
		DiscretizedFunc[] meanCurves = loadMeanCurves(
				runDir, tree, startIndex, sampleCount, gridReg, period);
		return buildHazardData(meanCurves, branchMaps, gridReg);
	}

	private static DiscretizedFunc[] loadMeanCurves(File runDir, LogicTree<?> tree,
			int startIndex, int sampleCount, GriddedRegion gridReg, double period) throws IOException {
		double[] xValues = null;
		double[][] sums = null;
		try (HazardConvergenceCalcs.BranchCurveLoader curveLoader =
				new HazardConvergenceCalcs.BranchCurveLoader(runDir)) {
			for (int b=startIndex; b<startIndex+sampleCount; b++) {
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

		int runSize(File runDir) {
			double[][] maps = branchMaps(runDir);
			Preconditions.checkState(maps != null, "Run is not in pooled data: %s", runDir.getName());
			return maps.length;
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
