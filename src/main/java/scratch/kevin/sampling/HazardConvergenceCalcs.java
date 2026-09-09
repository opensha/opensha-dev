package scratch.kevin.sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SplittableRandom;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.IntStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.xyz.AbstractXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.commons.util.RandomSeedUtils;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

public class HazardConvergenceCalcs {
	static final String MCS_REFERENCE_NAME = "Pooled MCS";
	static final String LOO_MCS_REFERENCE_NAME = "Pooled MCS, leave span out";
	static final String POOLED_SOBOL_REFERENCE_NAME = "Pooled Sobol";
	static final String LOO_SOBOL_REFERENCE_NAME = "Pooled Sobol, leave one out";
	private static final int MAX_RUN_LOAD_THREADS = 4;

	public static void main(String[] args) throws IOException {
		File outputDir = new File(PaperPaths.FIGURES_DIR, "hazard_convergence");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Couldn't create output directory: %s", outputDir.getAbsolutePath());

		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;

		Table<SamplingMethod, Integer, List<File>> runDirs = HashBasedTable.create();

		/*
		 * MCS runs
		 */
		runDirs.put(SamplingMethod.MONTE_CARLO, 20000, List.of(
				new File(PaperPaths.INVS_DIR, "2026_07_17-nshm27-AMSAM-20000samples-mcs"),
				new File(PaperPaths.INVS_DIR, "2026_09_03-nshm27-AMSAM-20000samples-mcs-unique_seed"),
				new File(PaperPaths.INVS_DIR, "2026_09_05-nshm27-AMSAM-20000samples-mcs-unique_seed-2")
				));

		/*
		 * Sobol runs
		 */
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 512, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-512samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-512samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 1024, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-1024samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-1024samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 2048, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-2048samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-2048samples-sobol_scrambled-unique_seed")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 4096, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-4096samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_25-nshm27-AMSAM-4096samples-sobol_scrambled-unique_seed"),
				new File(PaperPaths.INVS_DIR, "2026_08_27-nshm27-AMSAM-4096samples-sobol_scrambled-unique_seed-2"),
				new File(PaperPaths.INVS_DIR, "2026_08_27-nshm27-AMSAM-4096samples-sobol_scrambled-unique_seed-3")
				));
		runDirs.put(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 8192, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_28-nshm27-AMSAM-8192samples-sobol_scrambled"),
				new File(PaperPaths.INVS_DIR, "2026_08_28-nshm27-AMSAM-8192samples-sobol_scrambled-unique_seed"),
				new File(PaperPaths.INVS_DIR, "2026_08_29-nshm27-AMSAM-8192samples-sobol_scrambled-unique_seed-2"),
				new File(PaperPaths.INVS_DIR, "2026_08_29-nshm27-AMSAM-8192samples-sobol_scrambled-unique_seed-3")
				));
		
		/*
		 * LHS runs
		 */
		runDirs.put(SamplingMethod.LATIN_HYPERCUBE, 4096, List.of(
				new File(PaperPaths.INVS_DIR, "2026_09_08-nshm27-AMSAM-4096samples-lhs"),
				new File(PaperPaths.INVS_DIR, "2026_09_08-nshm27-AMSAM-4096samples-lhs-unique_seed")
//				new File(PaperPaths.INVS_DIR, "2026_09_08-nshm27-AMSAM-4096samples-lhs-unique_seed-2"),
//				new File(PaperPaths.INVS_DIR, "2026_09_08-nshm27-AMSAM-4096samples-lhs-unique_seed-3")
				));
		
		/*
		 * Pairwise-LHS runs
		 */
		runDirs.put(SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE, 4096, List.of(
				new File(PaperPaths.INVS_DIR, "2026_08_28-nshm27-AMSAM-4096samples-lhs_pairwise"),
				new File(PaperPaths.INVS_DIR, "2026_08_28-nshm27-AMSAM-4096samples-lhs_pairwise-unique_seed"),
				new File(PaperPaths.INVS_DIR, "2026_08_29-nshm27-AMSAM-4096samples-lhs_pairwise-unique_seed-2"),
				new File(PaperPaths.INVS_DIR, "2026_08_29-nshm27-AMSAM-4096samples-lhs_pairwise-unique_seed-3")
				));

		List<RunSpec> sobolRuns = loadRunSpecs(SamplingMethod.OWEN_SCRAMBLED_SOBOL,
				runDirs.row(SamplingMethod.OWEN_SCRAMBLED_SOBOL));
		List<RunSpec> fullDesignRuns = new ArrayList<>();
		for (SamplingMethod method : runDirs.rowKeySet().stream().sorted().toList())
			if (method != SamplingMethod.MONTE_CARLO && method != SamplingMethod.OWEN_SCRAMBLED_SOBOL)
				fullDesignRuns.addAll(loadRunSpecs(method, runDirs.row(method)));
		List<RunSpec> mcsRuns = loadRunSpecs(SamplingMethod.MONTE_CARLO,
				runDirs.row(SamplingMethod.MONTE_CARLO));
		System.out.println("MCS reference has "+mcsRuns.stream().mapToInt(RunSpec::maxSamples).sum()
				+" branches across "+mcsRuns.size()+" runs");

		GriddedRegion gridReg = GriddedRegion.fromFeature(
				Feature.read(new File(mcsRuns.get(0).directory, "gridded_region.geojson")));

		double[] periods = { 0d, 1d };
		String[] periodNames = { "PGA", "1 s SA" };
		String[] periodPrefixes = { "pga", "1s_sa" };
		for (int p=0; p<periods.length; p++) {
			File periodDir = new File(outputDir, periodPrefixes[p]+"_"+rp.name().toLowerCase());
			Preconditions.checkState(periodDir.exists() || periodDir.mkdir(),
					"Couldn't create output directory: %s", periodDir.getAbsolutePath());
			System.out.println("\n========== "+periodNames[p]+", "+rp+" ==========");
			runSamplingConvergence(sobolRuns, fullDesignRuns, mcsRuns,
					gridReg, periods[p], periodNames[p], rp, periodDir);
		}
	}

	private static List<RunSpec> loadRunSpecs(SamplingMethod method,
			Map<Integer, List<File>> runDirs) throws IOException {
		List<RunSpec> runs = new ArrayList<>();
		for (Map.Entry<Integer, List<File>> entry : new TreeMap<>(runDirs).entrySet()) {
			for (File dir : entry.getValue()) {
				Preconditions.checkState(dir.isDirectory(), "Run directory doesn't exist: %s", dir.getAbsolutePath());
				LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree_analysis.json"));
				Preconditions.checkState(tree.size() == entry.getKey(), "Expected %s branches in %s, have %s",
						entry.getKey(), dir.getName(), tree.size());
				Preconditions.checkState(tree.getSamplingMethod() == method,
						"Expected %s tree, have %s: %s", method, tree.getSamplingMethod(), dir.getName());
				runs.add(new RunSpec(dir.getName(), dir, tree, method,
						tree.getSamplingRandomSeed(), tree.size()));
			}
		}
		return runs;
	}

	private static List<RunSpec> loadRunSpecs(SamplingMethod method, List<File> runDirs) throws IOException {
		List<RunSpec> runs = new ArrayList<>();
		for (File dir : runDirs) {
			Preconditions.checkState(dir.isDirectory(), "Run directory doesn't exist: %s", dir.getAbsolutePath());
			LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree_analysis.json"));
			Preconditions.checkState(tree.getSamplingMethod() == method,
					"Expected %s tree, have %s: %s", method, tree.getSamplingMethod(), dir.getName());
			runs.add(new RunSpec(dir.getName(), dir, tree, method,
					tree.getSamplingRandomSeed(), tree.size()));
		}
		return runs;
	}

	private static void runSamplingConvergence(List<RunSpec> sobolRuns, List<RunSpec> fullDesignRuns,
			List<RunSpec> mcsRuns,
			GriddedRegion gridReg, double period, String periodName, ReturnPeriods rp, File outputDir) throws IOException {
		List<RunPeriodData> sobolData = loadRunPeriodData(sobolRuns, gridReg, period, rp);
		List<RunPeriodData> fullDesignData = loadRunPeriodData(fullDesignRuns, gridReg, period, rp);
		int[] sampleCounts = sobolData.stream().flatMap(data -> data.checkpoints().keySet().stream())
				.mapToInt(Integer::intValue).distinct().sorted().toArray();
		List<RunPeriodData> mcsData = loadRunPeriodData(mcsRuns, gridReg, period, rp, sampleCounts);

		PooledHazardData pooledMCSData = buildPooledHazardData(mcsData, null, gridReg, rp,
				MCS_REFERENCE_NAME, MCS_REFERENCE_NAME);
		ReferenceStatistics mcsReference = pooledMCSData.reference();
		Map<RunSpec, ReferenceStatistics> leaveOneOut = new LinkedHashMap<>();
		for (RunPeriodData data : sobolData)
			leaveOneOut.put(data.run(), buildPooledHazardData(sobolData, data.run(), gridReg, rp,
					POOLED_SOBOL_REFERENCE_NAME, LOO_SOBOL_REFERENCE_NAME).reference());
		PooledHazardData pooledSobolData = buildPooledHazardData(sobolData, null, gridReg, rp,
				POOLED_SOBOL_REFERENCE_NAME, LOO_SOBOL_REFERENCE_NAME);
		ReferenceStatistics pooledSobol = pooledSobolData.reference();

		writePooledHazardFiles(new File(outputDir, "pooled_mcs"), pooledMCSData, gridReg, period, rp);
		writePooledHazardFiles(new File(outputDir, "pooled_sobol"), pooledSobolData, gridReg, period, rp);

		List<ReferenceComparison> comparisons = new ArrayList<>();
		for (RunPeriodData data : sobolData) {
			ReferenceStatistics sobolReference = leaveOneOut.get(data.run());
			for (Map.Entry<Integer, HazardStatistics> entry : data.checkpoints().entrySet()) {
				appendComparisons(comparisons, data.run(), entry.getKey(), entry.getValue(), sobolReference, gridReg);
				appendComparisons(comparisons, data.run(), entry.getKey(), entry.getValue(), mcsReference, gridReg);
			}
		}
		// LHS and other full designs are evaluated at their generated sizes, without taking prefixes.
		for (RunPeriodData data : fullDesignData) {
			HazardStatistics statistics = data.checkpoints().get(data.run().maxSamples());
			appendComparisons(comparisons, data.run(), data.run().maxSamples(),
					statistics, pooledSobol, gridReg);
			appendComparisons(comparisons, data.run(), data.run().maxSamples(),
					statistics, mcsReference, gridReg);
		}
		List<Realization> mcsSpans = appendMCSSpanComparisons(comparisons, mcsData, sampleCounts,
				pooledSobol, gridReg, rp);
		writeReferenceComparisons(new File(outputDir, "reference_comparisons.csv"), comparisons);
		writeReferenceComparisonSummary(new File(outputDir, "reference_comparison_summary.csv"), comparisons);

		List<DoublingComparison> doublings = new ArrayList<>();
		for (RunPeriodData data : sobolData) {
			for (Map.Entry<Integer, HazardStatistics> entry : data.checkpoints().entrySet()) {
				int lowerCount = entry.getKey();
				HazardStatistics upper = data.checkpoints().get(2*lowerCount);
				if (upper == null)
					continue;
				for (ConvergenceMetric metric : ConvergenceMetric.values()) {
					MapComparison comparison = compare(upper.values(metric), entry.getValue().values(metric));
					doublings.add(new DoublingComparison(data.run(), lowerCount, 2*lowerCount,
							metric, comparison, gridReg.getLocation(comparison.maximumAbsoluteIndex())));
				}
			}
		}
		writeDoublingComparisons(new File(outputDir, "paired_doubling_comparisons.csv"), doublings);
		writeDoublingComparisonSummary(new File(outputDir, "paired_doubling_summary.csv"), doublings);

		List<RunPeriodData> allData = new ArrayList<>(sobolData);
		allData.addAll(fullDesignData);
		allData.addAll(mcsData);
		List<Realization> realizations = new ArrayList<>(mcsSpans);
		for (RunPeriodData data : allData)
			data.checkpoints().forEach((count, statistics) -> realizations.add(new Realization(data.run(), count, statistics)));
		List<RealizationPairComparison> realizationPairs = buildRealizationPairComparisons(realizations, gridReg);
		writeRealizationPairComparisons(new File(outputDir, "realization_pair_comparisons.csv"), realizationPairs);
		writeRealizationPairComparisonSummary(
				new File(outputDir, "realization_pair_summary.csv"), realizationPairs);

		HazardConvergencePlots.plotPeriod(outputDir, periodName);
	}

	private static List<RunPeriodData> loadRunPeriodData(List<RunSpec> runs,
			GriddedRegion gridReg, double period, ReturnPeriods rp, int... spanCounts) throws IOException {
		if (runs.isEmpty())
			return List.of();
		int threads = Math.min(MAX_RUN_LOAD_THREADS, runs.size());
		if (threads == 1)
			return List.of(loadRunPeriodData(runs.get(0), gridReg, period, rp, spanCounts));

		ExecutorService executor = Executors.newFixedThreadPool(threads);
		try {
			List<Future<RunPeriodData>> futures = new ArrayList<>(runs.size());
			for (RunSpec run : runs)
				futures.add(executor.submit(() -> loadRunPeriodData(run, gridReg, period, rp, spanCounts)));
			List<RunPeriodData> data = new ArrayList<>(runs.size());
			// Retrieve in input order so downstream CSV and plot ordering remains stable.
			for (Future<RunPeriodData> future : futures) {
				try {
					data.add(future.get());
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
					throw new IOException("Interrupted while loading hazard runs", e);
				} catch (ExecutionException e) {
					Throwable cause = e.getCause();
					if (cause instanceof IOException ioException)
						throw ioException;
					if (cause instanceof RuntimeException runtimeException)
						throw runtimeException;
					if (cause instanceof Error error)
						throw error;
					throw new IOException("Exception while loading hazard runs", cause);
				}
			}
			return data;
		} finally {
			executor.shutdownNow();
		}
	}

	private static RunPeriodData loadRunPeriodData(RunSpec run, GriddedRegion gridReg,
			double period, ReturnPeriods rp, int[] spanCounts) throws IOException {
		System.out.println("\nLoading "+run.method().getShortName()+" run: "+run.id());
		ModelHazardMaps maps = loadMaps(new File(run.directory(), "results_hazard.zip"),
				run.tree(), gridReg, period, rp);
		double[][] branchMaps = copyValues(maps.individual());
		File hazardResultsDir = new File(run.directory(), "results");

		double[] curveX = null;
		double[][] curveSums = null;
		Map<Integer, HazardStatistics> checkpoints = new TreeMap<>();
		Map<Integer, double[][]> curveBoundaries = new TreeMap<>();
		for (int b=0; b<run.maxSamples(); b++) {
			DiscretizedFunc[] curves = loadBranchCurves(hazardResultsDir, run.tree().getBranch(b), gridReg, period);
			if (curveSums == null) {
				Preconditions.checkState(curves.length == gridReg.getNodeCount());
				curveX = new double[curves[0].size()];
				for (int i=0; i<curveX.length; i++)
					curveX[i] = curves[0].getX(i);
				curveSums = new double[curves.length][curveX.length];
			}
			for (int n=0; n<curves.length; n++) {
				DiscretizedFunc curve = curves[n];
				Preconditions.checkState(curve.size() == curveX.length);
				for (int i=0; i<curveX.length; i++) {
					Preconditions.checkState((float)curve.getX(i) == (float)curveX[i]);
					curveSums[n][i] += curve.getY(i);
				}
			}
			int count = b+1;
			boolean fullRun = count == run.maxSamples();
			boolean sobolCheckpoint = run.method() == SamplingMethod.OWEN_SCRAMBLED_SOBOL
					&& count >= 512 && Integer.bitCount(count) == 1;
			if (Arrays.stream(spanCounts).anyMatch(size -> count % size == 0))
				curveBoundaries.put(count, Arrays.stream(curveSums).map(double[]::clone).toArray(double[][]::new));
			if (fullRun || sobolCheckpoint) {
				double[] curveMean = buildCurveMeanMap(curveSums, curveX, count, rp);
				checkpoints.put(count, calcHazardStatistics(branchMaps, count, curveMean));
				System.out.println("\tBuilt "+count+"-sample checkpoint");
			}
		}
		Preconditions.checkState(checkpoints.containsKey(run.maxSamples()));
		MapComparison archivedComparison = compare(
				checkpoints.get(run.maxSamples()).values(ConvergenceMetric.MEAN_HAZARD), copyValues(maps.mean()));
		System.out.println("\tCurve mean versus archived mean: "+archivedComparison);
		return new RunPeriodData(run, branchMaps, curveX, curveSums, checkpoints, curveBoundaries);
	}

	private static PooledHazardData buildPooledHazardData(List<RunPeriodData> allRuns,
			RunSpec excluded, GriddedRegion gridReg, ReturnPeriods rp,
			String pooledName, String leaveOneOutName) {
		int sampleCount = 0;
		int curveSize = -1;
		double[] curveX = null;
		for (RunPeriodData data : allRuns) {
			if (data.run().equals(excluded))
				continue;
			sampleCount += data.run().maxSamples();
			if (curveX == null) {
				curveX = data.curveX();
				curveSize = curveX.length;
			} else {
				Preconditions.checkState(curveSize == data.curveX().length);
				for (int i=0; i<curveSize; i++)
					Preconditions.checkState((float)curveX[i] == (float)data.curveX()[i]);
			}
		}
		Preconditions.checkState(sampleCount > 0);
		double[][] curveSums = new double[gridReg.getNodeCount()][curveSize];
		double[][] branchMaps = new double[sampleCount][];
		int destBranch = 0;
		for (RunPeriodData data : allRuns) {
			if (data.run().equals(excluded))
				continue;
			for (int n=0; n<curveSums.length; n++)
				for (int i=0; i<curveSize; i++)
					curveSums[n][i] += data.curveSums()[n][i];
			for (double[] branchMap : data.branchMaps())
				branchMaps[destBranch++] = branchMap;
		}
		Preconditions.checkState(destBranch == sampleCount);
		double[] curveMean = buildCurveMeanMap(curveSums, curveX, sampleCount, rp);
		String name = excluded == null ? pooledName : leaveOneOutName;
		ReferenceStatistics reference = new ReferenceStatistics(name,
				excluded == null ? null : excluded.id(), sampleCount,
				calcHazardStatistics(branchMaps, sampleCount, curveMean));
		return new PooledHazardData(reference, curveX, curveSums, curveMean);
	}

	private static void writePooledHazardFiles(File outputDir, PooledHazardData pooled,
			GriddedRegion gridReg, double period, ReturnPeriods rp) throws IOException {
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Couldn't create pooled hazard directory: %s", outputDir.getAbsolutePath());
		int sampleCount = pooled.reference().sampleCount();
		double[] curveX = pooled.curveX();
		double[][] curveSums = pooled.curveSums();
		DiscretizedFunc[] meanCurves = new DiscretizedFunc[gridReg.getNodeCount()];
		for (int n=0; n<meanCurves.length; n++) {
			double[] meanY = new double[curveX.length];
			for (int i=0; i<meanY.length; i++)
				meanY[i] = curveSums[n][i]/sampleCount;
			meanCurves[n] = new LightFixedXFunc(curveX, meanY);
		}

		String curvesName = SolHazardMapCalc.getCSV_FileName("mean_curves", period)+".gz";
		SolHazardMapCalc.writeCurvesCSV(new File(outputDir, curvesName),
				meanCurves, gridReg.getNodeList());
		String mapName = "mean_"+MPJ_LogicTreeHazardCalc.mapPrefix(period, rp)+".txt";
		AbstractXYZ_DataSet.writeXYZFile(
				new GriddedGeoDataSet(gridReg, pooled.meanMap()), new File(outputDir, mapName));
		System.out.println("Wrote "+pooled.reference().name()+" hazard curves and map to "
				+outputDir.getAbsolutePath());
	}

	/** Uses disjoint full spans within each run; leftover branches still contribute to every reference. */
	static List<Realization> appendMCSSpanComparisons(List<ReferenceComparison> comparisons,
			List<RunPeriodData> runs, int[] sampleCounts, ReferenceStatistics sobolReference,
			GriddedRegion gridReg, ReturnPeriods rp) {
		List<Realization> spans = new ArrayList<>();
		double[][] allMaps = runs.stream().flatMap(data -> Arrays.stream(data.branchMaps()))
				.toArray(double[][]::new);
		double[] curveX = runs.get(0).curveX();
		double[][] totalCurves = new double[gridReg.getNodeCount()][curveX.length];
		for (RunPeriodData data : runs) {
			Preconditions.checkState(Arrays.equals(curveX, data.curveX()), "MCS curve grids differ");
			for (int n=0; n<totalCurves.length; n++)
				for (int i=0; i<curveX.length; i++)
					totalCurves[n][i] += data.curveSums()[n][i];
		}
		int runOffset = 0;
		for (RunPeriodData data : runs) {
			for (int count : sampleCounts) {
				for (int start=0; start <= data.branchMaps().length-count; start+=count) {
					int end = start+count;
					Preconditions.checkState(allMaps.length-count > 1, "Too few MCS reference samples after exclusion");
					double[][] spanCurves = new double[totalCurves.length][curveX.length];
					double[][] remainingCurves = new double[totalCurves.length][curveX.length];
					double[][] before = start == 0 ? null : data.curveBoundaries().get(start);
					double[][] after = data.curveBoundaries().get(end);
					for (int n=0; n<totalCurves.length; n++) {
						for (int i=0; i<curveX.length; i++) {
							spanCurves[n][i] = after[n][i] - (before == null ? 0d : before[n][i]);
							remainingCurves[n][i] = totalCurves[n][i]-spanCurves[n][i];
						}
					}
					// Only copy row references; the maps themselves remain shared and read-only.
					double[][] spanMaps = Arrays.copyOfRange(data.branchMaps(), start, end);
					double[][] remainingMaps = excludeSpan(allMaps, runOffset+start, runOffset+end);
					HazardStatistics statistics = calcHazardStatistics(spanMaps, count,
							buildCurveMeanMap(spanCurves, curveX, count, rp));
					// Retain the small statistics arrays for pair comparisons; curve and map arrays are unnecessary there.
					RunSpec spanRun = new RunSpec(data.run().id()+"["+start+","+end+")",
							data.run().directory(), data.run().tree(), data.run().method(), data.run().seed(), count);
					spans.add(new Realization(spanRun, count, statistics));
					ReferenceStatistics reference = new ReferenceStatistics(LOO_MCS_REFERENCE_NAME,
							data.run().id()+"["+start+","+end+")", remainingMaps.length,
							calcHazardStatistics(remainingMaps, remainingMaps.length,
									buildCurveMeanMap(remainingCurves, curveX, remainingMaps.length, rp)));
					appendComparisons(comparisons, data.run(), count, statistics, reference, gridReg, start);
					appendComparisons(comparisons, data.run(), count, statistics, sobolReference, gridReg, start);
				}
				System.out.println("Compared "+data.run().id()+": "+count+"-sample MCS spans");
			}
			runOffset += data.branchMaps().length;
		}
		return spans;
	}

	static double[][] excludeSpan(double[][] rows, int start, int end) {
		Preconditions.checkArgument(start >= 0 && start < end && end <= rows.length);
		double[][] remaining = new double[rows.length-(end-start)][];
		System.arraycopy(rows, 0, remaining, 0, start);
		System.arraycopy(rows, end, remaining, start, rows.length-end);
		return remaining;
	}

	static double[] buildCurveMeanMap(double[][] curveSums, double[] curveX,
			int sampleCount, ReturnPeriods rp) {
		double[] map = new double[curveSums.length];
		double[] meanY = new double[curveX.length];
		for (int n=0; n<curveSums.length; n++) {
			for (int i=0; i<meanY.length; i++)
				meanY[i] = curveSums[n][i]/sampleCount;
			LightFixedXFunc curve = new LightFixedXFunc(curveX, meanY);
			if (rp.oneYearProb > curve.getMaxY())
				map[n] = 0d;
			else if (rp.oneYearProb < curve.getMinY())
				map[n] = curve.getMaxX();
			else
				map[n] = curve.getFirstInterpolatedX_inLogXLogYDomain(rp.oneYearProb);
		}
		return map;
	}

	static HazardStatistics calcHazardStatistics(double[][] branchMaps, int sampleCount,
			double[] curveMean) {
		Preconditions.checkArgument(sampleCount > 1 && sampleCount <= branchMaps.length);
		Preconditions.checkArgument(curveMean.length == branchMaps[0].length);
		int numSites = curveMean.length;
		double[] standardDeviation = new double[numSites];
		double[] iqr = new double[numSites];
		double[] central68 = new double[numSites];
		double[] central95 = new double[numSites];
		int workers = Math.min(numSites, Runtime.getRuntime().availableProcessors());
		// Partition sites rather than allocating a sample array for every parallel-stream element.
		IntStream.range(0, workers).parallel().forEach(worker -> {
			double[] values = new double[sampleCount];
			for (int n=worker; n<numSites; n+=workers) {
				double sum = 0d;
				double sumSquares = 0d;
				for (int b=0; b<sampleCount; b++) {
					double value = branchMaps[b][n];
					values[b] = value;
					sum += value;
					sumSquares += value*value;
				}
				double mapMean = sum/sampleCount;
				double varianceNumerator = sumSquares - sum*mapMean;
				if (varianceNumerator < 0d && varianceNumerator > -1e-12*sumSquares)
					varianceNumerator = 0d;
				Preconditions.checkState(varianceNumerator >= 0d,
						"Negative variance numerator at site "+n+": "+varianceNumerator);
				standardDeviation[n] = Math.sqrt(varianceNumerator/sampleCount);
				Arrays.sort(values);
				iqr[n] = empiricalFractile(values, 0.75)-empiricalFractile(values, 0.25);
				central68[n] = empiricalFractile(values, 0.84)-empiricalFractile(values, 0.16);
				central95[n] = empiricalFractile(values, 0.975)-empiricalFractile(values, 0.025);
			}
		});
		Map<ConvergenceMetric, double[]> metricValues = new EnumMap<>(ConvergenceMetric.class);
		metricValues.put(ConvergenceMetric.MEAN_HAZARD, curveMean);
		metricValues.put(ConvergenceMetric.STANDARD_DEVIATION, standardDeviation);
		metricValues.put(ConvergenceMetric.IQR, iqr);
		metricValues.put(ConvergenceMetric.CENTRAL_68_RANGE, central68);
		metricValues.put(ConvergenceMetric.CENTRAL_95_RANGE, central95);
		return new HazardStatistics(metricValues);
	}

	private static double empiricalFractile(double[] sortedValues, double fractile) {
		Preconditions.checkArgument(sortedValues.length > 0 && fractile >= 0d && fractile <= 1d);
		double previousValue = sortedValues[0];
		int index = 1;
		while (index < sortedValues.length && (float)sortedValues[index] == (float)previousValue)
			index++;
		double previousCDF = (double)index/sortedValues.length;
		if (fractile <= previousCDF)
			return previousValue;
		while (index < sortedValues.length) {
			double value = sortedValues[index++];
			while (index < sortedValues.length && (float)sortedValues[index] == (float)value)
				index++;
			double cdf = (double)index/sortedValues.length;
			if (fractile == cdf)
				return value;
			if (fractile < cdf) {
				double relative = (fractile-previousCDF)/(cdf-previousCDF);
				return previousValue + relative*(value-previousValue);
			}
			previousValue = value;
			previousCDF = cdf;
		}
		return sortedValues[sortedValues.length-1];
	}

	private static void appendComparisons(List<ReferenceComparison> comparisons, RunSpec run,
			int sampleCount, HazardStatistics statistics, ReferenceStatistics reference, GriddedRegion gridReg) {
		appendComparisons(comparisons, run, sampleCount, statistics, reference, gridReg, 0);
	}

	private static void appendComparisons(List<ReferenceComparison> comparisons, RunSpec run,
			int sampleCount, HazardStatistics statistics, ReferenceStatistics reference, GriddedRegion gridReg,
			int startIndex) {
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			MapComparison comparison = compare(statistics.values(metric), reference.statistics().values(metric));
			comparisons.add(new ReferenceComparison(run, sampleCount, reference.name(), reference.sampleCount(),
					metric, comparison, gridReg.getLocation(comparison.maximumAbsoluteIndex()), startIndex));
		}
	}

	private static void writeReferenceComparisons(File file, List<ReferenceComparison> comparisons)
			throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Run", "Sampling method", "Seed", "Maximum run size", "Sample count",
				"Reference", "Reference sample count",
				"Metric", "Spatial mean % change", "Spatial mean absolute % change",
				"Spatial P95 absolute % change", "Maximum absolute % change", "Minimum % change",
				"Maximum % change", "Worst longitude", "Worst latitude", "Start index", "End index (exclusive)");
		for (ReferenceComparison row : comparisons) {
			MapComparison comparison = row.comparison();
			csv.addLine(row.run().id(), row.run().method().name(), row.run().seed()+"",
					row.run().maxSamples()+"", row.sampleCount()+"",
					row.referenceName(), row.referenceSampleCount()+"", row.metric().label,
					comparison.meanPercentChange()+"", comparison.meanAbsolutePercentChange()+"",
					comparison.p95AbsolutePercentChange()+"", comparison.maximumAbsolutePercentChange()+"",
					comparison.minimumPercentChange()+"", comparison.maximumPercentChange()+"",
					row.worstLocation().lon+"", row.worstLocation().lat+"",
					row.startIndex()+"", (row.startIndex()+row.sampleCount())+"");
		}
		csv.writeToFile(file);
	}

	private static void writeReferenceComparisonSummary(File file, List<ReferenceComparison> comparisons)
			throws IOException {
		Map<ReferenceComparisonGroup, List<ReferenceComparison>> groups = new LinkedHashMap<>();
		for (ReferenceComparison comparison : comparisons) {
			ReferenceComparisonGroup group = new ReferenceComparisonGroup(comparison.run().method(),
					comparison.sampleCount(),
					comparison.referenceName(), comparison.metric());
			groups.computeIfAbsent(group, key -> new ArrayList<>()).add(comparison);
		}
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sampling method", "Sample count", "Reference", "Metric", "Spatial summary", "Realizations",
				"Mean", "Minimum", "P2.5", "P16", "P50", "P84", "P97.5", "Maximum");
		for (Map.Entry<ReferenceComparisonGroup, List<ReferenceComparison>> entry : groups.entrySet()) {
			ReferenceComparisonGroup group = entry.getKey();
			for (ConvergenceSummary summary : ConvergenceSummary.values()) {
				double[] values = new double[entry.getValue().size()];
				for (int i=0; i<values.length; i++)
					values[i] = summary.value(entry.getValue().get(i).comparison());
				addSummaryLine(csv, List.of(group.method().name(), group.sampleCount()+"",
						group.referenceName(), group.metric().label,
						summary.label, values.length+""), values);
			}
		}
		csv.writeToFile(file);
	}

	private static void writeDoublingComparisons(File file, List<DoublingComparison> comparisons)
			throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Run", "Sampling method", "Seed", "Maximum run size", "Lower sample count",
				"Upper sample count", "Metric",
				"Spatial mean % change", "Spatial mean absolute % change", "Spatial P95 absolute % change",
				"Maximum absolute % change", "Minimum % change", "Maximum % change",
				"Worst longitude", "Worst latitude");
		for (DoublingComparison row : comparisons) {
			MapComparison comparison = row.comparison();
			csv.addLine(row.run().id(), row.run().method().name(), row.run().seed()+"",
					row.run().maxSamples()+"", row.lowerCount()+"",
					row.upperCount()+"", row.metric().label, comparison.meanPercentChange()+"",
					comparison.meanAbsolutePercentChange()+"", comparison.p95AbsolutePercentChange()+"",
					comparison.maximumAbsolutePercentChange()+"", comparison.minimumPercentChange()+"",
					comparison.maximumPercentChange()+"", row.worstLocation().lon+"", row.worstLocation().lat+"");
		}
		csv.writeToFile(file);
	}

	private static void writeDoublingComparisonSummary(File file, List<DoublingComparison> comparisons)
			throws IOException {
		Map<DoublingComparisonGroup, List<DoublingComparison>> groups = new LinkedHashMap<>();
		for (DoublingComparison comparison : comparisons) {
			DoublingComparisonGroup group = new DoublingComparisonGroup(comparison.run().method(), comparison.lowerCount(),
					comparison.upperCount(), comparison.metric());
			groups.computeIfAbsent(group, key -> new ArrayList<>()).add(comparison);
		}
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sampling method", "Lower sample count", "Upper sample count", "Metric",
				"Spatial summary", "Realizations",
				"Mean", "Minimum", "P2.5", "P16", "P50", "P84", "P97.5", "Maximum");
		for (Map.Entry<DoublingComparisonGroup, List<DoublingComparison>> entry : groups.entrySet()) {
			DoublingComparisonGroup group = entry.getKey();
			for (ConvergenceSummary summary : ConvergenceSummary.values()) {
				double[] values = new double[entry.getValue().size()];
				for (int i=0; i<values.length; i++)
					values[i] = summary.value(entry.getValue().get(i).comparison());
				addSummaryLine(csv, List.of(group.method().name(), group.lowerCount()+"", group.upperCount()+"",
						group.metric().label,
						summary.label, values.length+""), values);
			}
		}
		csv.writeToFile(file);
	}

	static List<RealizationPairComparison> buildRealizationPairComparisons(
			List<Realization> allData, GriddedRegion gridReg) {
		Map<MethodCount, List<Realization>> groups = new LinkedHashMap<>();
		for (Realization data : allData) {
			MethodCount group = new MethodCount(data.run().method(), data.sampleCount());
			groups.computeIfAbsent(group, unused -> new ArrayList<>()).add(data);
		}
		List<RealizationPairComparison> comparisons = new ArrayList<>();
		for (Map.Entry<MethodCount, List<Realization>> entry : groups.entrySet()) {
			List<Realization> data = entry.getValue();
			for (int i=0; i<data.size(); i++) {
				HazardStatistics first = data.get(i).statistics();
				for (int j=i+1; j<data.size(); j++) {
					HazardStatistics second = data.get(j).statistics();
					for (ConvergenceMetric metric : ConvergenceMetric.values()) {
						MapComparison comparison = compare(first.values(metric), second.values(metric));
						comparisons.add(new RealizationPairComparison(entry.getKey().method(),
								entry.getKey().sampleCount(), data.get(i).run(), data.get(j).run(), metric,
								comparison, gridReg.getLocation(comparison.maximumAbsoluteIndex())));
					}
				}
			}
		}
		return comparisons;
	}

	private static void writeRealizationPairComparisons(File file,
			List<RealizationPairComparison> comparisons) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sampling method", "Sample count", "First run", "First seed", "Second run", "Second seed",
				"Metric", "Spatial mean % change", "Spatial mean absolute % change",
				"Spatial P95 absolute % change", "Maximum absolute % change", "Minimum % change",
				"Maximum % change", "Worst longitude", "Worst latitude");
		for (RealizationPairComparison row : comparisons) {
			MapComparison comparison = row.comparison();
			csv.addLine(row.method().name(), row.sampleCount()+"", row.first().id(), row.first().seed()+"",
					row.second().id(), row.second().seed()+"", row.metric().label,
					comparison.meanPercentChange()+"", comparison.meanAbsolutePercentChange()+"",
					comparison.p95AbsolutePercentChange()+"", comparison.maximumAbsolutePercentChange()+"",
					comparison.minimumPercentChange()+"", comparison.maximumPercentChange()+"",
					row.worstLocation().lon+"", row.worstLocation().lat+"");
		}
		csv.writeToFile(file);
	}

	private static void writeRealizationPairComparisonSummary(File file,
			List<RealizationPairComparison> comparisons) throws IOException {
		Map<RealizationPairComparisonGroup, List<RealizationPairComparison>> groups = new LinkedHashMap<>();
		for (RealizationPairComparison comparison : comparisons) {
			RealizationPairComparisonGroup group = new RealizationPairComparisonGroup(comparison.method(),
					comparison.sampleCount(), comparison.metric());
			groups.computeIfAbsent(group, unused -> new ArrayList<>()).add(comparison);
		}
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sampling method", "Sample count", "Metric", "Spatial summary", "Realization pairs",
				"Mean", "Minimum", "P2.5", "P16", "P50", "P84", "P97.5", "Maximum");
		for (Map.Entry<RealizationPairComparisonGroup, List<RealizationPairComparison>> entry : groups.entrySet()) {
			RealizationPairComparisonGroup group = entry.getKey();
			for (ConvergenceSummary summary : ConvergenceSummary.values()) {
				double[] values = new double[entry.getValue().size()];
				for (int i=0; i<values.length; i++)
					values[i] = summary.value(entry.getValue().get(i).comparison());
				addSummaryLine(csv, List.of(group.method().name(), group.sampleCount()+"",
						group.metric().label, summary.label, values.length+""), values);
			}
		}
		csv.writeToFile(file);
	}

	private static void addSummaryLine(CSVFile<String> csv, List<String> prefix, double[] values) {
		List<String> line = new ArrayList<>(prefix);
		line.add(StatUtils.mean(values)+"");
		line.add(StatUtils.min(values)+"");
		line.add(percentile(values, 2.5)+"");
		line.add(percentile(values, 16d)+"");
		line.add(percentile(values, 50d)+"");
		line.add(percentile(values, 84d)+"");
		line.add(percentile(values, 97.5)+"");
		line.add(StatUtils.max(values)+"");
		csv.addLine(line);
	}

	/**
	 * Resamples the reference MCS branches with replacement. Smaller resample sizes describe the error expected from
	 * ordinary MCS at those sample counts; the full-size resample describes uncertainty in the reference itself.
	 * Results are deterministic for a given seed even though replicates are evaluated in parallel.
	 */
	static void runBootstrapTests(ModelHazardMaps reference, int[] sampleCounts, int numReplicates, long seed,
			File outputDir, String outputPrefix) throws IOException {
		Preconditions.checkArgument(numReplicates > 1);
		Preconditions.checkArgument(sampleCounts.length > 0);

		double[][] branchValues = copyValues(reference.individual());
		StatisticMaps referenceStats = calcStatistics(branchValues, null, branchValues.length);
		MapComparison archivedMeanComparison = compare(referenceStats.mean(), copyValues(reference.mean()));
		System.out.println("Branch-map arithmetic mean versus archived curve-derived mean (expected to differ): "
				+archivedMeanComparison);

		CSVFile<String> replicateCSV = new CSVFile<>(true);
		replicateCSV.addLine("Sample count", "Replicate", "Metric", "Spatial mean % change",
				"Spatial mean absolute % change", "Minimum % change", "Maximum % change");
		CSVFile<String> summaryCSV = new CSVFile<>(true);
		summaryCSV.addLine("Sample count", "Metric", "Spatial summary", "Mean", "Minimum", "P2.5", "P16", "P50",
				"P84", "P97.5", "Maximum");
		CSVFile<String> siteCSV = new CSVFile<>(true);
		siteCSV.addLine("Sample count", "Metric", "Longitude", "Latitude", "Reference value",
				"Mean % change", "Std. dev. % change", "P2.5", "P16", "P50", "P84", "P97.5");

		for (int sampleCount : sampleCounts) {
			Preconditions.checkArgument(sampleCount > 1);
			System.out.println("\nBootstrapping "+sampleCount+" branches x "+numReplicates+" replicates");
			BootstrapReplicate[] replicates = new BootstrapReplicate[numReplicates];
			IntStream.range(0, numReplicates).parallel().forEach(r -> {
				long replicateSeed = RandomSeedUtils.uniqueSeedCombination(seed, sampleCount, r);
				int[] counts = bootstrapCounts(branchValues.length, sampleCount, replicateSeed);
				StatisticMaps statistics = calcStatistics(branchValues, counts, sampleCount);
				Map<HazardMetric, MapComparison> comparisons = new EnumMap<>(HazardMetric.class);
				for (HazardMetric metric : HazardMetric.values())
					comparisons.put(metric, compare(metric.values(statistics), metric.values(referenceStats)));
				replicates[r] = new BootstrapReplicate(statistics, comparisons);
			});

			appendReplicateCSVs(replicateCSV, summaryCSV, siteCSV, reference.mean(), referenceStats,
					sampleCount, replicates);
			for (HazardMetric metric : HazardMetric.values()) {
				double[] meanAbs = new double[numReplicates];
				for (int r=0; r<numReplicates; r++)
					meanAbs[r] = replicates[r].comparisons().get(metric).meanAbsolutePercentChange();
				System.out.println("\t"+metric.label+" mean absolute spatial % change: "+formatDistribution(meanAbs));
			}
		}

		File replicateFile = new File(outputDir, outputPrefix+"_bootstrap_replicates.csv");
		File summaryFile = new File(outputDir, outputPrefix+"_bootstrap_summary.csv");
		File siteFile = new File(outputDir, outputPrefix+"_bootstrap_sites.csv");
		replicateCSV.writeToFile(replicateFile);
		summaryCSV.writeToFile(summaryFile);
		siteCSV.writeToFile(siteFile);
		System.out.println("\nWrote bootstrap results:");
		System.out.println("\t"+replicateFile.getAbsolutePath());
		System.out.println("\t"+summaryFile.getAbsolutePath());
		System.out.println("\t"+siteFile.getAbsolutePath());
	}

	private static void appendReplicateCSVs(CSVFile<String> replicateCSV, CSVFile<String> summaryCSV,
			CSVFile<String> siteCSV, GriddedGeoDataSet referenceMap, StatisticMaps referenceStats,
			int sampleCount, BootstrapReplicate[] replicates) {
		for (HazardMetric metric : HazardMetric.values()) {
			for (int r=0; r<replicates.length; r++) {
				MapComparison comparison = replicates[r].comparisons().get(metric);
				replicateCSV.addLine(sampleCount+"", r+"", metric.label,
						comparison.meanPercentChange()+"", comparison.meanAbsolutePercentChange()+"",
						comparison.minimumPercentChange()+"", comparison.maximumPercentChange()+"");
			}
			for (ComparisonQuantity quantity : ComparisonQuantity.values()) {
				double[] values = new double[replicates.length];
				for (int r=0; r<replicates.length; r++)
					values[r] = quantity.value(replicates[r].comparisons().get(metric));
				appendDistribution(summaryCSV, sampleCount, metric.label, quantity.label, values);
			}

			double[] referenceValues = metric.values(referenceStats);
			double[] changes = new double[replicates.length];
			for (int n=0; n<referenceValues.length; n++) {
				for (int r=0; r<replicates.length; r++)
					changes[r] = percentChange(metric.values(replicates[r].statistics())[n], referenceValues[n]);
				Location location = referenceMap.getLocation(n);
				siteCSV.addLine(sampleCount+"", metric.label, location.lon+"", location.lat+"",
						referenceValues[n]+"", StatUtils.mean(changes)+"", Math.sqrt(StatUtils.variance(changes))+"",
						percentile(changes, 2.5)+"", percentile(changes, 16d)+"", percentile(changes, 50d)+"",
						percentile(changes, 84d)+"", percentile(changes, 97.5)+"");
			}
		}
	}

	private static void appendDistribution(CSVFile<String> csv, int sampleCount, String metric,
			String quantity, double[] values) {
		csv.addLine(sampleCount+"", metric, quantity, StatUtils.mean(values)+"", StatUtils.min(values)+"",
				percentile(values, 2.5)+"",
				percentile(values, 16d)+"", percentile(values, 50d)+"", percentile(values, 84d)+"",
				percentile(values, 97.5)+"", StatUtils.max(values)+"");
	}

	private static String formatDistribution(double[] values) {
		return "median="+(float)percentile(values, 50d)+", 95%=["+(float)percentile(values, 2.5)
				+", "+(float)percentile(values, 97.5)+"]";
	}

	private static double percentile(double[] values, double percentile) {
		// StatUtils does not modify its input, which lets us reuse the per-replicate arrays for each percentile.
		return StatUtils.percentile(values, percentile);
	}

	private static int[] bootstrapCounts(int numBranches, int sampleCount, long seed) {
		int[] counts = new int[numBranches];
		SplittableRandom random = new SplittableRandom(seed);
		for (int i=0; i<sampleCount; i++)
			counts[random.nextInt(numBranches)]++;
		return counts;
	}

	private static double[][] copyValues(GriddedGeoDataSet[] maps) {
		double[][] values = new double[maps.length][];
		for (int i=0; i<maps.length; i++)
			values[i] = copyValues(maps[i]);
		return values;
	}

	private static double[] copyValues(GriddedGeoDataSet map) {
		double[] values = new double[map.size()];
		for (int i=0; i<values.length; i++)
			values[i] = map.get(i);
		return values;
	}

	private static StatisticMaps calcStatistics(double[][] branchValues, int[] counts, int sampleCount) {
		Preconditions.checkArgument(branchValues.length > 0);
		Preconditions.checkArgument(sampleCount > 1);
		if (counts != null)
			Preconditions.checkArgument(counts.length == branchValues.length);
		int numSites = branchValues[0].length;
		double[] sums = new double[numSites];
		double[] sumSquares = new double[numSites];
		for (int b=0; b<branchValues.length; b++) {
			int count = counts == null ? 1 : counts[b];
			if (count == 0)
				continue;
			double[] values = branchValues[b];
			Preconditions.checkState(values.length == numSites);
			for (int n=0; n<numSites; n++) {
				double value = values[n];
				sums[n] += count*value;
				sumSquares[n] += count*value*value;
			}
		}

		double[] means = new double[numSites];
		double[] standardDeviations = new double[numSites];
		double[] coefficientsOfVariation = new double[numSites];
		for (int n=0; n<numSites; n++) {
			double mean = sums[n]/sampleCount;
			// LogicTreeHazardCompare treats the sampled tree as the full population, so use the population variance.
			double varianceNumerator = sumSquares[n] - sums[n]*mean;
			if (varianceNumerator < 0d && varianceNumerator > -1e-12*sumSquares[n])
				varianceNumerator = 0d;
			Preconditions.checkState(varianceNumerator >= 0d,
					"Negative variance numerator at site "+n+": "+varianceNumerator);
			double standardDeviation = Math.sqrt(varianceNumerator/sampleCount);
			means[n] = mean;
			standardDeviations[n] = standardDeviation;
			coefficientsOfVariation[n] = standardDeviation/mean;
		}
		return new StatisticMaps(means, standardDeviations, coefficientsOfVariation);
	}

	private static MapComparison compare(double[] testValues, double[] referenceValues) {
		Preconditions.checkArgument(testValues.length == referenceValues.length);
		double sum = 0d;
		double sumAbsolute = 0d;
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		double maxAbsolute = Double.NEGATIVE_INFINITY;
		int maxAbsoluteIndex = -1;
		double[] absoluteChanges = new double[testValues.length];
		for (int i=0; i<testValues.length; i++) {
			double change = percentChange(testValues[i], referenceValues[i]);
			double absolute = Math.abs(change);
			sum += change;
			sumAbsolute += absolute;
			absoluteChanges[i] = absolute;
			min = Math.min(min, change);
			max = Math.max(max, change);
			if (absolute > maxAbsolute) {
				maxAbsolute = absolute;
				maxAbsoluteIndex = i;
			}
		}
		Arrays.sort(absoluteChanges);
		return new MapComparison(sum/testValues.length, sumAbsolute/testValues.length,
				empiricalFractile(absoluteChanges, 0.95), maxAbsolute, maxAbsoluteIndex, min, max);
	}

	private static double percentChange(double testValue, double referenceValue) {
		Preconditions.checkState(Double.isFinite(testValue));
		Preconditions.checkState(Double.isFinite(referenceValue) && referenceValue > 0d,
				"Reference values must be finite and positive: %s", referenceValue);
		return 100d*(testValue/referenceValue - 1d);
	}

	private static String mapFilePrefix(double period, ReturnPeriods rp) {
		String perStr = period == 0d ? "pga" : (float)period+"s";
		return perStr+"_"+rp.name();
	}

	private static GriddedGeoDataSet readMap(GriddedRegion gridReg, InputStream is) throws IOException {
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		BufferedReader bRead = new BufferedReader(new InputStreamReader(is));
		String line = bRead.readLine();
		int index = 0;
		while (line != null) {
			line = line.trim();
			if (!line.startsWith("#")) {
				StringTokenizer tok = new StringTokenizer(line);
				double lon = Double.parseDouble(tok.nextToken());
				double lat = Double.parseDouble(tok.nextToken());
				double val = Double.parseDouble(tok.nextToken());
				Location loc = new Location(lat, lon);
				Preconditions.checkState(LocationUtils.areSimilar(loc, gridReg.getLocation(index)));
				xyz.set(index++, val);
			}
			line = bRead.readLine();
		}
		Preconditions.checkState(index == gridReg.getNodeCount());
		bRead.close();
		return xyz;
	}

	private static ModelHazardMaps loadMaps(File hazardZip, LogicTree<?> tree, GriddedRegion gridReg,
			double period, ReturnPeriods rp) throws ZipException, IOException {
		System.out.println("Loading maps from "+hazardZip.getAbsolutePath());
		try (ZipFile zip = new ZipFile(hazardZip)) {
			String suffix = mapFilePrefix(period, rp)+".txt";
			String meanEntryName = "mean_map_"+suffix;
			ZipEntry meanEntry = zip.getEntry(meanEntryName);
			Preconditions.checkNotNull(meanEntry, "Entry doesn't exist in %s: %s", hazardZip.getAbsolutePath(), meanEntryName);
			GriddedGeoDataSet meanMap = readMap(gridReg, zip.getInputStream(meanEntry));

			GriddedGeoDataSet[] individual = new GriddedGeoDataSet[tree.size()];
			for (int i=0; i<individual.length; i++) {
				LogicTreeBranch<?> branch = tree.getBranch(i);
				String mapName = branch.buildFileName()+"/map_"+suffix;
				ZipEntry mapEntry = zip.getEntry(mapName);
				Preconditions.checkNotNull(mapEntry, "Entry doesn't exist in %s: %s", hazardZip.getAbsolutePath(), mapName);
				individual[i] = readMap(gridReg, zip.getInputStream(mapEntry));
			}

			System.out.println("\tLoaded mean & "+individual.length+" individual");

			return new ModelHazardMaps(meanMap, individual);
		}
	}

	private static DiscretizedFunc[] loadBranchCurves(File hazardResultsDir, LogicTreeBranch<?> branch,
			GriddedRegion gridReg, double period) throws IOException {
		File runDir = branch.getBranchDirectory(hazardResultsDir, false);
		File hazardDir = new File(runDir, "hazard_"+(float)gridReg.getSpacing()+"deg_grid_seis_INCLUDE");
		Preconditions.checkState(hazardDir.exists(), "Hazard directory doesn't exist: %s", hazardDir.getAbsolutePath());
		File hazardFile = new File(hazardDir, SolHazardMapCalc.getCSV_FileName("curves", period));
		if (!hazardFile.exists())
			hazardFile = new File(hazardDir, SolHazardMapCalc.getCSV_FileName("curves", period)+".gz");
		Preconditions.checkState(hazardFile.exists(), "Hazard curves file doesn't exist: %s", hazardFile.getAbsolutePath());
		// this will detect that it's gzipped
		CSVFile<String> csv = CSVFile.readFile(hazardFile, true);
		return SolHazardMapCalc.loadCurvesCSV(csv, gridReg);
	}

	private enum HazardMetric {
		MEAN("Mean hazard") {
			@Override double[] values(StatisticMaps maps) { return maps.mean(); }
		},
		STANDARD_DEVIATION("SD of hazard") {
			@Override double[] values(StatisticMaps maps) { return maps.standardDeviation(); }
		},
		COEFFICIENT_OF_VARIATION("CV of hazard") {
			@Override double[] values(StatisticMaps maps) { return maps.coefficientOfVariation(); }
		};

		final String label;
		private HazardMetric(String label) {
			this.label = label;
		}
		abstract double[] values(StatisticMaps maps);
	}

	enum ConvergenceMetric {
		MEAN_HAZARD("Mean hazard", "Mean"),
		STANDARD_DEVIATION("SD of hazard", "SD"),
		IQR("Interquartile range", "IQR"),
		CENTRAL_68_RANGE("Central 68% range", "68%"),
		CENTRAL_95_RANGE("Central 95% range", "95%");

		final String label;
		final String shortLabel;
		private ConvergenceMetric(String label, String shortLabel) {
			this.label = label;
			this.shortLabel = shortLabel;
		}

		static ConvergenceMetric fromLabel(String label) {
			for (ConvergenceMetric metric : values())
				if (metric.label.equals(label))
					return metric;
			throw new IllegalArgumentException("Unknown convergence metric: "+label);
		}
	}

	enum ConvergenceSummary {
		MEAN_SIGNED("Spatial mean % change") {
			@Override double value(MapComparison comparison) { return comparison.meanPercentChange(); }
		},
		MEAN_ABSOLUTE("Spatial mean absolute % change") {
			@Override double value(MapComparison comparison) { return comparison.meanAbsolutePercentChange(); }
		},
		P95_ABSOLUTE("Spatial P95 absolute % change") {
			@Override double value(MapComparison comparison) { return comparison.p95AbsolutePercentChange(); }
		},
		MAXIMUM_ABSOLUTE("Maximum absolute % change") {
			@Override double value(MapComparison comparison) { return comparison.maximumAbsolutePercentChange(); }
		};

		final String label;
		private ConvergenceSummary(String label) {
			this.label = label;
		}

		static ConvergenceSummary fromLabel(String label) {
			for (ConvergenceSummary summary : values())
				if (summary.label.equals(label))
					return summary;
			throw new IllegalArgumentException("Unknown convergence summary: "+label);
		}

		abstract double value(MapComparison comparison);
	}

	private enum ComparisonQuantity {
		MEAN("Spatial mean % change") {
			@Override double value(MapComparison comparison) { return comparison.meanPercentChange(); }
		},
		MEAN_ABSOLUTE("Spatial mean absolute % change") {
			@Override double value(MapComparison comparison) { return comparison.meanAbsolutePercentChange(); }
		},
		MINIMUM("Minimum % change") {
			@Override double value(MapComparison comparison) { return comparison.minimumPercentChange(); }
		},
		MAXIMUM("Maximum % change") {
			@Override double value(MapComparison comparison) { return comparison.maximumPercentChange(); }
		};

		final String label;
		private ComparisonQuantity(String label) {
			this.label = label;
		}
		abstract double value(MapComparison comparison);
	}

	record ModelHazardMaps(GriddedGeoDataSet mean, GriddedGeoDataSet[] individual) {}

	private record StatisticMaps(double[] mean, double[] standardDeviation, double[] coefficientOfVariation) {}

	record MapComparison(double meanPercentChange, double meanAbsolutePercentChange,
			double p95AbsolutePercentChange, double maximumAbsolutePercentChange, int maximumAbsoluteIndex,
			double minimumPercentChange, double maximumPercentChange) {
		@Override
		public String toString() {
			return "mean="+(float)meanPercentChange+"%, abs="+(float)meanAbsolutePercentChange
					+"%, p95abs="+(float)p95AbsolutePercentChange+"%, maxAbs="
					+(float)maximumAbsolutePercentChange+"%, range=["+(float)minimumPercentChange+"%, "
					+(float)maximumPercentChange+"%]";
		}
	}

	private record BootstrapReplicate(StatisticMaps statistics, Map<HazardMetric, MapComparison> comparisons) {}

	record RunSpec(String id, File directory, LogicTree<?> tree, SamplingMethod method,
			long seed, int maxSamples) {}

	record RunPeriodData(RunSpec run, double[][] branchMaps, double[] curveX,
			double[][] curveSums, Map<Integer, HazardStatistics> checkpoints,
			Map<Integer, double[][]> curveBoundaries) {}

	record HazardStatistics(Map<ConvergenceMetric, double[]> metricValues) {
		double[] values(ConvergenceMetric metric) {
			return Preconditions.checkNotNull(metricValues.get(metric));
		}
	}

	record ReferenceStatistics(String name, String excludedRun, int sampleCount,
			HazardStatistics statistics) {}

	private record PooledHazardData(ReferenceStatistics reference, double[] curveX,
			double[][] curveSums, double[] meanMap) {}

	record ReferenceComparison(RunSpec run, int sampleCount, String referenceName,
			int referenceSampleCount, ConvergenceMetric metric, MapComparison comparison,
			Location worstLocation, int startIndex) {}

	private record ReferenceComparisonGroup(SamplingMethod method, int sampleCount, String referenceName,
			ConvergenceMetric metric) {}

	private record DoublingComparison(RunSpec run, int lowerCount, int upperCount,
			ConvergenceMetric metric, MapComparison comparison, Location worstLocation) {}

	private record DoublingComparisonGroup(SamplingMethod method, int lowerCount, int upperCount,
			ConvergenceMetric metric) {}

	private record MethodCount(SamplingMethod method, int sampleCount) {}

	record Realization(RunSpec run, int sampleCount, HazardStatistics statistics) {}

	record RealizationPairComparison(SamplingMethod method, int sampleCount,
			RunSpec first, RunSpec second, ConvergenceMetric metric, MapComparison comparison,
			Location worstLocation) {}

	private record RealizationPairComparisonGroup(SamplingMethod method, int sampleCount,
			ConvergenceMetric metric) {}

}
