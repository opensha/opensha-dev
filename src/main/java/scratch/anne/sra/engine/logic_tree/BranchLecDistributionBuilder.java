package scratch.anne.sra.engine.logic_tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Turns a branch-weighted sample list for ANY single metric (a
 * likelihood-level's loss-at-that-level across branches, or the EAL
 * distribution across branches, or any future metric with the same
 * shape) into a summarized LevelDistribution: weighted mean, weighted
 * percentiles, and missing-branch bookkeeping.
 *
 * <p>Deliberately metric-agnostic -- "requestedLevel" is just a label
 * (a likelihood value, or a sentinel like EAL_LEVEL for the EAL row);
 * nothing here cares what the underlying values represent, only that
 * each BranchLevelSample carries (branchIndex, branchWeight, value).</p>
 */
public final class BranchLecDistributionBuilder {

    /** Sentinel "level" used to key the EAL distribution alongside likelihood-level distributions. */
    public static final double EAL_LEVEL = Double.NaN;

    public record BranchLevelSample(int branchIndex, double branchWeight, double loss) {}

    /** One metric's full weighted distribution: mean, percentiles, and missing-branch info. */
    public record LevelDistribution(double requestedLevel,
                                     List<BranchLevelSample> samples,
                                     double mean,
                                     Map<Double, Double> percentiles,
                                     int missingCount,
                                     double missingFraction) {}

    private BranchLecDistributionBuilder() {}

    /**
     * Summarizes ONE metric's samples into a LevelDistribution -- the
     * single building block used for both EAL and every likelihood
     * level. Callers gather samples (e.g. via streaming the
     * contributions Parquet, or reading branch-eals.csv) and hand them
     * here one metric at a time.
     */
    public static LevelDistribution summarize(double requestedLevel,
                                                List<BranchLevelSample> samples,
                                                double[] percentilesOfInterest) {

        int missingCount = 0;
        for (BranchLevelSample s : samples)
            if (Double.isNaN(s.loss()))
                missingCount++;
        double missingFraction = samples.isEmpty() ? 0d : (double) missingCount / samples.size();

        double mean = weightedMean(samples);
        Map<Double, Double> percentiles = weightedPercentiles(samples, percentilesOfInterest);

        return new LevelDistribution(requestedLevel, samples, mean, percentiles, missingCount, missingFraction);
    }

    /**
     * Convenience: summarizes multiple levels' worth of already-gathered
     * samples in one call (e.g. after a single streaming pass has
     * collected samplesByLevel for every likelihood level).
     */
    public static Map<Double, LevelDistribution> summarizeAll(Map<Double, List<BranchLevelSample>> samplesByLevel,
                                                                double[] percentilesOfInterest) {
        Map<Double, LevelDistribution> distributions = new java.util.LinkedHashMap<>();
        for (Map.Entry<Double, List<BranchLevelSample>> entry : samplesByLevel.entrySet()) {
            LevelDistribution dist = summarize(entry.getKey(), entry.getValue(), percentilesOfInterest);
            distributions.put(entry.getKey(), dist);

            System.out.printf("Level %s: %d/%d branches (%.1f%%) never reached this exceedance rate%n",
                    entry.getKey(), dist.missingCount(), entry.getValue().size(), dist.missingFraction() * 100d);
        }
        return distributions;
    }

    private static double weightedMean(List<BranchLevelSample> samples) {
        double weightedSum = 0d;
        double totalWeight = 0d;
        for (BranchLevelSample s : samples) {
            if (Double.isNaN(s.loss()))
                continue;
            weightedSum += s.loss() * s.branchWeight();
            totalWeight += s.branchWeight();
        }
        return totalWeight > 0d ? weightedSum / totalWeight : Double.NaN;
    }

    /**
     * Weighted percentiles over branch values for one metric, ignoring
     * NaN values (missing/unreached) and renormalizing weights over the
     * remaining branches.
     */
    private static Map<Double, Double> weightedPercentiles(List<BranchLevelSample> samples, double[] percentilesOfInterest) {

        List<BranchLevelSample> valid = new ArrayList<>();
        for (BranchLevelSample s : samples)
            if (!Double.isNaN(s.loss()))
                valid.add(s);

        valid.sort((a, b) -> Double.compare(a.loss(), b.loss()));

        double totalWeight = 0d;
        for (BranchLevelSample s : valid)
            totalWeight += s.branchWeight();

        Map<Double, Double> result = new java.util.LinkedHashMap<>();

        if (valid.isEmpty() || totalWeight <= 0d) {
            for (double p : percentilesOfInterest)
                result.put(p, Double.NaN);
            return result;
        }

        for (double p : percentilesOfInterest) {
            double targetFraction = p / 100d;
            double running = 0d;
            double selected = valid.get(valid.size() - 1).loss();
            for (BranchLevelSample s : valid) {
                running += s.branchWeight();
                if (running / totalWeight >= targetFraction) {
                    selected = s.loss();
                    break;
                }
            }
            result.put(p, selected);
        }

        return result;
    }
}