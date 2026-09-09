package scratch.kevin.sampling;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.numbers.core.Precision;
import org.opensha.commons.data.sampling.ArrayPointSet;
import org.opensha.commons.data.sampling.CategoricalSamplingDimension;
import org.opensha.commons.data.sampling.ContinuousSamplingDimension;
import org.opensha.commons.data.sampling.DimensionedPointSet;
import org.opensha.commons.data.sampling.PermutedPointSet;
import org.opensha.commons.data.sampling.PointSet;
import org.opensha.commons.data.sampling.SamplingDimension;
import org.opensha.commons.data.sampling.generator.*;
import org.opensha.commons.data.sampling.optimization.IncrementalPointSetScorer;
import org.opensha.commons.data.sampling.optimization.PointSetHillClimber;
import org.opensha.commons.data.sampling.optimization.PointSetOptimizationResult;
import org.opensha.commons.data.sampling.optimization.QuantizedIncrementalPointSetScorer;
import org.opensha.commons.data.sampling.scoring.ProjectionDiscrepancyScore;
import org.opensha.commons.data.sampling.scoring.ProjectionDiscrepancyScorer;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

public class InitialSamplingTests {

	public static void main(String[] args) {
		int numSamples = powerOfTwo(12);
//		int numSamples = 5000;
//		int numDimensions = 10;
		int numDimensions = 30;
		int numCategorical = 5;
		int maxCategoriesPerDimension = 5;
		int numScoringDimensions = 3;
		int numIterations = 100000;
		
		Random r = new Random(123456789l); // repeatable
//		Random r = new Random(); // different each time
		
		List<SamplingDimension> dimensions = buildDimensions(numDimensions, numCategorical, maxCategoriesPerDimension, r);
		
//		PointSetGenerator generator = new MonteCarloPointSetGenerator(r);
//		PointSetGenerator generator = new LatinHypercubePointSetGenerator(r);
//		PointSetGenerator generator = new SobolPointSetGenerator();
//		PointSetGenerator generator = new SobolPointSetGenerator(1);
		PointSetGenerator generator = new OwenScrambledSobolPointSetGenerator(r);
		
		System.out.println("\nGenerating "+numSamples+" "+numDimensions+"-D samples with generator: "+generator);
		
		Stopwatch watch = Stopwatch.createStarted();
		PointSet samples = generator.generate(numSamples, numDimensions);
		watch.stop();
		System.out.println("Done in "+timeStr(watch));
		int[] debugPoints = {0, 1, 2, 3, 4, numSamples-1};
		for (int debugPoint : debugPoints) {
			System.out.print("Point "+debugPoint+":\t[");
			for (int d=0; d<numDimensions; d++) {
				if (d > 0)
					System.out.print(", ");
				System.out.print((float)samples.get(debugPoint, d));
			}
			System.out.println("]");
		}

		ProjectionDiscrepancyScorer exactScorer = ProjectionDiscrepancyScorer.exact(16);
		
		System.out.println("\nScoring continuous case to order "+numScoringDimensions);
		watch.reset().start();
		ProjectionDiscrepancyScore score = exactScorer.score(samples, numScoringDimensions);
		watch.stop();
		System.out.println("Done in "+timeStr(watch));
		System.out.println("Continuous score:\t"+score);
//		System.exit(0);

//		System.out.println("\nRe-scoring continuous case using quantized scorer");
//		ProjectionDiscrepancyScorer quantizedScorer = ProjectionDiscrepancyScorer.quantized(100);
//		watch.reset().start();
//		ProjectionDiscrepancyScore quantizedScore = quantizedScorer.score(samples, numScoringDimensions);
//		watch.stop();
//		System.out.println("Done in "+timeStr(watch));
//		System.out.println("Continuous score:\t"+quantizedScore);
		
		// now make some categorical
		DimensionedPointSet dimensioned = new DimensionedPointSet(samples, dimensions);
		System.out.println("\nScoring dimensioned set");
		watch.reset().start();
		ProjectionDiscrepancyScore dimensionedScore = exactScorer.score(dimensioned, numScoringDimensions);
		watch.stop();
		System.out.println("Done in "+timeStr(watch));
		System.out.println("Dimensioned score:\t"+dimensionedScore);
		
		System.out.println("\nImproving pairwise with "+numIterations+" hill-climbing iterations");
		watch.reset().start();
		PermutedPointSet permuted = PermutedPointSet.independentDimensions(dimensioned);
		IncrementalPointSetScorer incrementalScorer = new QuantizedIncrementalPointSetScorer(permuted, 100);
		PointSetOptimizationResult result = PointSetHillClimber.optimize(incrementalScorer, numIterations, r);
		watch.stop();
		System.out.println("Done in "+timeStr(watch));
		System.out.println("Optimization result:\t"+result);
		
		System.out.println("\nScoring optimized version");
		watch.reset().start();
		ProjectionDiscrepancyScore optimizedScore = exactScorer.score(permuted, numScoringDimensions);
		watch.stop();
		System.out.println("Done in "+timeStr(watch));
		System.out.println("Optimized score:\t"+optimizedScore);
		
		// these tests can be uncommented if we need to check JVM or PointSet implementation performance again
//		System.out.println("\nConverting optimized version to an ArrayPointSet");
//		watch.reset().start();
//		DimensionedPointSet materialized = new DimensionedPointSet(new ArrayPointSet(permuted), dimensions);
//		watch.stop();
//		System.out.println("Done in "+timeStr(watch));
//		System.out.println("Optimized score:\t"+optimizedScore);
//		
//		System.out.println("\nRe-scoring ArrayPointSet view of optimized version");
//		watch.reset().start();
//		ProjectionDiscrepancyScore materializedScore = scorer.score(materialized, numScoringDimensions);
//		watch.stop();
//		System.out.println("Done in "+timeStr(watch));
//		System.out.println("Optimized score:\t"+materializedScore);
//		
//		System.out.println("\nRe-scoring the initial continuos case (JVM test)");
//		watch.reset().start();
//		ProjectionDiscrepancyScore score2 = scorer.score(samples, numScoringDimensions);
//		watch.stop();
//		System.out.println("Done in "+timeStr(watch));
//		System.out.println("Continuous re-score:\t"+score2);
	}
	
	static List<SamplingDimension> buildDimensions(int numDimensions, int numCategorical, int maxCategoriesPerDimension, Random r) {
		System.out.println("Building categories");
		List<SamplingDimension> dimensions = new ArrayList<>(numDimensions);
		Preconditions.checkState(numCategorical <= numDimensions);
		Preconditions.checkState(maxCategoriesPerDimension >= 2);
		// fill with sequential initially
		for (int i=0; i<numDimensions; i++)
			dimensions.add(ContinuousSamplingDimension.INSTANCE);
		int myNumCategorical = 0;
		while (myNumCategorical < numCategorical) {
			int index = r.nextInt(numDimensions);
			while (dimensions.get(index) != ContinuousSamplingDimension.INSTANCE)
				index = r.nextInt(numDimensions);
			int numCategories = maxCategoriesPerDimension > 2 ? 2 + r.nextInt(maxCategoriesPerDimension-2) : 2;
			boolean even = r.nextBoolean();
			double[] weights = new double[numCategories];
			if (even) {
				Arrays.fill(weights, 1d/numCategories);
			} else {
				for (int i=0; i<numCategories; i++)
					weights[i] = r.nextDouble();
				// snap to 5% bounds
				for (int i=0; i<numCategories; i++)
					weights[i] = Math.max(0.05, Math.round(weights[i]*20)/20d);
				double sum = StatUtils.sum(weights);
				if (!Precision.equals(sum, 1d, 1e-5)) {
					if (sum < 0.999) {
						// simple, add weight
						double remainder = 1d - sum;
						weights[r.nextInt(numCategories)] += remainder;
					} else {
						// need to remove weight
						double remainder = sum - 1d;
						while (!Precision.equals(sum, 1d, 1e-5)) {
							Preconditions.checkState(remainder > 0.04999);
							int removeIndex = r.nextInt(numCategories);
							if (weights[removeIndex] > 0.0999) {
								weights[removeIndex] -= 0.05d;
								sum = StatUtils.sum(weights);
								remainder = sum - 1d;
							}
						}
					}
				}
			}
			System.out.println("Replacing index "+index+" with categorical weights: "+Arrays.toString(weights));
			dimensions.set(index, CategoricalSamplingDimension.forWeights(weights));
			myNumCategorical++;
		}
		return dimensions;
	}
	
	private static int powerOfTwo(int n) {
		Preconditions.checkState(n <= 30);
		return 1 << n;
	}
	
	private static final DecimalFormat timeDF = new DecimalFormat("0.0");
	private static String timeStr(Stopwatch watch) {
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		if (secs < 90d)
			return timeDF.format(secs)+" s";
		double mins = secs/60d;
		if (mins < 90d)
			return timeDF.format(mins)+" m";
		double hours = mins / 60d;
		return timeDF.format(hours)+" h";
	}

}
