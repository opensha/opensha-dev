package scratch.kevin.sampling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.TreeSet;
import java.util.random.RandomGenerator;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.sampling.ContinuousSamplingDimension;
import org.opensha.commons.data.sampling.PermutedPointSet;
import org.opensha.commons.data.sampling.PointSet;
import org.opensha.commons.data.sampling.SamplingDimension;
import org.opensha.commons.data.sampling.optimization.PointSetHillClimber;
import org.opensha.commons.data.sampling.optimization.PointSetHillClimber.Result;
import org.opensha.commons.data.sampling.optimization.PointSetObjective.SwapSession;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.commons.util.RandomSeedUtils;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

/**
 * Small diagnostic for deciding whether the production PO-LHS and CDO-LHS iteration limits are long enough. Each
 * trial follows one optimization trajectory through logarithmic checkpoints, continuing to four times the production
 * limit. The most useful summary value is "Further improvement through terminal": if it is already small at the
 * production limit, additional iterations have little left to gain on that trajectory.
 */
public class OptimizationIterationConvergenceCalc {

	private static final int[] SAMPLE_COUNTS = { 256, 512, 1024, 2048, 4096, 8192, 16384 };
	private static final SamplingMethod[] METHODS = {
			SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE,
			SamplingMethod.CENTERED_DISCREPANCY_OPTIMIZED_LATIN_HYPERCUBE
	};
	private static final int DEFAULT_TRIALS = 1;
	private static final int TERMINAL_ITERATION_MULTIPLIER = 4;

	private enum DimensionConfiguration {
		CONTINUOUS_10D("continuous_10d") {
			@Override List<SamplingDimension> dimensions() {
				return Collections.nCopies(10, ContinuousSamplingDimension.INSTANCE);
			}
		},
		NSHM27_AMSAM("nshm27_amsam") {
			@Override List<SamplingDimension> dimensions() {
				return SamplingScoreFigures.getDimsNSHM27_AmSam();
			}
		};

		private final String prefix;

		DimensionConfiguration(String prefix) {
			this.prefix = prefix;
		}

		abstract List<SamplingDimension> dimensions();
	}

	private record Checkpoint(long iterations, long acceptedInterval, long acceptedCumulative,
			double objective, double intervalSeconds, double cumulativeSeconds) {}

	private record TrialResult(SamplingMethod method, int sampleCount, int trial, long productionIterations,
			double initialObjective, List<Checkpoint> checkpoints) {

		Checkpoint productionCheckpoint() {
			return checkpoint(productionIterations);
		}

		Checkpoint terminalCheckpoint() {
			return checkpoints.getLast();
		}

		private Checkpoint checkpoint(long iterations) {
			return checkpoints.stream().filter(checkpoint -> checkpoint.iterations() == iterations).findFirst()
					.orElseThrow(() -> new IllegalStateException("Missing checkpoint at "+iterations+" iterations"));
		}
	}

	public static void main(String[] args) throws IOException {
		DimensionConfiguration configuration = args.length > 0
				? DimensionConfiguration.valueOf(args[0].toUpperCase(Locale.US))
				: DimensionConfiguration.NSHM27_AMSAM;
		int trials = args.length > 1 ? Integer.parseInt(args[1]) : DEFAULT_TRIALS;
		Preconditions.checkArgument(trials > 0, "Trials must be positive");

		List<SamplingDimension> dimensions = configuration.dimensions();
		File outputDir = new File(PaperPaths.FIGURES_DIR,
				"scores/optimization_convergence/"+configuration.prefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs(),
				"Couldn't create output directory: %s", outputDir.getAbsolutePath());

		System.out.println("Testing "+dimensions.size()+" "+configuration.prefix+" dimensions with "
				+trials+" trials per sample count and method");
		System.out.println("Each trajectory continues through "+TERMINAL_ITERATION_MULTIPLIER
				+"x the production iteration count\n");

		List<TrialResult> results = new ArrayList<>();
		for (int sampleCount : SAMPLE_COUNTS) {
			for (int trial=0; trial<trials; trial++) {
				long lhsSeed = RandomSeedUtils.seedForStrings("optimization-convergence", configuration.prefix,
						"lhs", Integer.toString(sampleCount), Integer.toString(trial));
				PointSet lhs = SamplingMethod.LATIN_HYPERCUBE.prepare(sampleCount, dimensions, lhsSeed);
				for (SamplingMethod method : METHODS)
					results.add(runTrial(configuration, lhs, method, sampleCount, trial));
			}
		}

		writeCheckpointCSV(new File(outputDir, "checkpoints.csv"), configuration, dimensions.size(), results);
		writeProductionSummaryCSV(new File(outputDir, "production_summary.csv"), results);
		printProductionSummary(results);
		System.out.println("\nWrote results to "+outputDir.getAbsolutePath());
	}

	private static TrialResult runTrial(DimensionConfiguration configuration, PointSet lhs, SamplingMethod method,
			int sampleCount, int trial) {
		long productionIterations = method.pairwiseIterations(sampleCount);
		long[] checkpoints = buildCheckpoints(productionIterations);
		PermutedPointSet permuted = PermutedPointSet.independentDimensions(lhs);
		SwapSession session = method.getObjective().prepare(permuted);
		long optimizationSeed = RandomSeedUtils.seedForStrings("optimization-convergence", configuration.prefix,
				method.name(), Integer.toString(sampleCount), Integer.toString(trial));
		RandomGenerator random = new Random(optimizationSeed);
		double initialObjective = session.getCurrentValue();
		List<Checkpoint> values = new ArrayList<>(checkpoints.length+1);
		values.add(new Checkpoint(0L, 0L, 0L, initialObjective, 0d, 0d));
		long previousIterations = 0L;
		long acceptedCumulative = 0L;
		double cumulativeSeconds = 0d;

		System.out.println(method.getShortName()+", N="+sampleCount+", trial="+(trial+1)
				+": production="+productionIterations+", terminal="+checkpoints[checkpoints.length-1]);
		for (long checkpoint : checkpoints) {
			long intervalIterations = checkpoint-previousIterations;
			Stopwatch watch = Stopwatch.createStarted();
			Result result = PointSetHillClimber.optimize(session, intervalIterations, random);
			watch.stop();
			double intervalSeconds = watch.elapsed().toNanos()*1e-9;
			cumulativeSeconds += intervalSeconds;
			acceptedCumulative += result.getAcceptedSwaps();
			values.add(new Checkpoint(checkpoint, result.getAcceptedSwaps(), acceptedCumulative,
					session.getCurrentValue(), intervalSeconds, cumulativeSeconds));
			System.out.printf(Locale.US, "\t%9d (%5.2fx): objective=%.8g, accepted=%d (%.3f%%), time=%.2fs%n",
					checkpoint, (double)checkpoint/productionIterations, session.getCurrentValue(),
					result.getAcceptedSwaps(), 100d*result.getAcceptedSwaps()/intervalIterations, intervalSeconds);
			previousIterations = checkpoint;
		}
		return new TrialResult(method, sampleCount, trial, productionIterations, initialObjective, List.copyOf(values));
	}

	private static long[] buildCheckpoints(long productionIterations) {
		TreeSet<Long> checkpoints = new TreeSet<>();
		for (int divisor : new int[] { 64, 32, 16, 8, 4, 2 })
			checkpoints.add(Math.max(1L, productionIterations/divisor));
		checkpoints.add(productionIterations);
		for (int multiplier=2; multiplier<=TERMINAL_ITERATION_MULTIPLIER; multiplier*=2)
			checkpoints.add(Math.multiplyExact(productionIterations, multiplier));
		return checkpoints.stream().mapToLong(Long::longValue).toArray();
	}

	private static void writeCheckpointCSV(File file, DimensionConfiguration configuration, int dimensions,
			List<TrialResult> results) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Configuration", "Dimensions", "Method", "Sample count", "Trial", "Iterations",
				"Production iterations", "Iteration multiple", "Objective", "Objective / initial",
				"Accepted in interval", "Interval acceptance rate", "Accepted cumulative",
				"Cumulative acceptance rate", "Interval seconds", "Cumulative seconds",
				"Further improvement through terminal (%)");
		for (TrialResult result : results) {
			double terminalObjective = result.terminalCheckpoint().objective();
			for (Checkpoint checkpoint : result.checkpoints()) {
				long previousIterations = checkpoint.iterations() == 0L ? 0L
						: previousIterations(result.checkpoints(), checkpoint);
				long intervalIterations = checkpoint.iterations()-previousIterations;
				double intervalAcceptance = intervalIterations == 0L ? 0d
						: (double)checkpoint.acceptedInterval()/intervalIterations;
				double cumulativeAcceptance = checkpoint.iterations() == 0L ? 0d
						: (double)checkpoint.acceptedCumulative()/checkpoint.iterations();
				double furtherImprovement = percentFurtherImprovement(checkpoint.objective(), terminalObjective);
				csv.addLine(configuration.prefix, Integer.toString(dimensions), result.method().getShortName(),
						Integer.toString(result.sampleCount()), Integer.toString(result.trial()+1),
						Long.toString(checkpoint.iterations()), Long.toString(result.productionIterations()),
						Double.toString((double)checkpoint.iterations()/result.productionIterations()),
						Double.toString(checkpoint.objective()),
						Double.toString(checkpoint.objective()/result.initialObjective()),
						Long.toString(checkpoint.acceptedInterval()), Double.toString(intervalAcceptance),
						Long.toString(checkpoint.acceptedCumulative()), Double.toString(cumulativeAcceptance),
						Double.toString(checkpoint.intervalSeconds()), Double.toString(checkpoint.cumulativeSeconds()),
						Double.toString(furtherImprovement));
			}
		}
		csv.writeToFile(file);
	}

	private static long previousIterations(List<Checkpoint> checkpoints, Checkpoint target) {
		long previous = 0L;
		for (Checkpoint checkpoint : checkpoints) {
			if (checkpoint == target)
				return previous;
			previous = checkpoint.iterations();
		}
		throw new IllegalStateException("Checkpoint not found");
	}

	private static void writeProductionSummaryCSV(File file, List<TrialResult> results) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Method", "Sample count", "Trials", "Production iterations", "Terminal iterations",
				"Mean production objective", "Mean terminal objective",
				"Mean further improvement (%)", "Median further improvement (%)",
				"Minimum further improvement (%)", "Maximum further improvement (%)",
				"Mean post-production acceptance rate", "Mean production seconds", "Mean terminal seconds");
		for (SamplingMethod method : METHODS) {
			for (int sampleCount : SAMPLE_COUNTS) {
				List<TrialResult> matching = results.stream().filter(result -> result.method() == method
						&& result.sampleCount() == sampleCount).toList();
				if (matching.isEmpty())
					continue;
				double[] productionObjectives = new double[matching.size()];
				double[] terminalObjectives = new double[matching.size()];
				double[] furtherImprovements = new double[matching.size()];
				double[] postProductionAcceptance = new double[matching.size()];
				double[] productionSeconds = new double[matching.size()];
				double[] terminalSeconds = new double[matching.size()];
				for (int i=0; i<matching.size(); i++) {
					TrialResult result = matching.get(i);
					Checkpoint production = result.productionCheckpoint();
					Checkpoint terminal = result.terminalCheckpoint();
					productionObjectives[i] = production.objective();
					terminalObjectives[i] = terminal.objective();
					furtherImprovements[i] = percentFurtherImprovement(production.objective(), terminal.objective());
					postProductionAcceptance[i] = (double)(terminal.acceptedCumulative()-production.acceptedCumulative())
							/(terminal.iterations()-production.iterations());
					productionSeconds[i] = production.cumulativeSeconds();
					terminalSeconds[i] = terminal.cumulativeSeconds();
				}
				TrialResult first = matching.getFirst();
				csv.addLine(method.getShortName(), Integer.toString(sampleCount), Integer.toString(matching.size()),
						Long.toString(first.productionIterations()), Long.toString(first.terminalCheckpoint().iterations()),
						Double.toString(StatUtils.mean(productionObjectives)),
						Double.toString(StatUtils.mean(terminalObjectives)),
						Double.toString(StatUtils.mean(furtherImprovements)),
						Double.toString(StatUtils.percentile(furtherImprovements, 50d)),
						Double.toString(StatUtils.min(furtherImprovements)),
						Double.toString(StatUtils.max(furtherImprovements)),
						Double.toString(StatUtils.mean(postProductionAcceptance)),
						Double.toString(StatUtils.mean(productionSeconds)),
						Double.toString(StatUtils.mean(terminalSeconds)));
			}
		}
		csv.writeToFile(file);
	}

	private static void printProductionSummary(List<TrialResult> results) {
		System.out.println("\n========== PRODUCTION CUTOFF SUMMARY ==========");
		for (SamplingMethod method : METHODS) {
			System.out.println(method.getName());
			for (int sampleCount : SAMPLE_COUNTS) {
				double[] improvements = results.stream().filter(result -> result.method() == method
						&& result.sampleCount() == sampleCount)
						.mapToDouble(result -> percentFurtherImprovement(result.productionCheckpoint().objective(),
								result.terminalCheckpoint().objective())).toArray();
				System.out.printf(Locale.US,
						"\tN=%5d: continuing from 1x to %dx improves objective by median %.3f%%, range [%.3f%%, %.3f%%]%n",
						sampleCount, TERMINAL_ITERATION_MULTIPLIER,
						StatUtils.percentile(improvements, 50d), StatUtils.min(improvements), StatUtils.max(improvements));
			}
		}
	}

	private static double percentFurtherImprovement(double objective, double terminalObjective) {
		return objective == 0d ? 0d : 100d*(objective-terminalObjective)/objective;
	}
}
