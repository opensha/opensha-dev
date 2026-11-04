package scratch.anne.sra.domain.risk;

import java.util.ArrayList;
import java.util.List;

import scratch.anne.sra.domain.hazard.RuptureKey;

/**
 * A finished loss-exceedance curve: parallel arrays of loss threshold
 * and cumulativeRate, ordered loss increasing / cumulativeRate
 * decreasing. Produced by DiscreteRankAccumulation or
 * LogNormalCOVAccumulation; consumed by anything that wants to plot it,
 * export it wholesale, or sample it at fixed rate levels.
 *
 * <p>Optionally carries a RuptureKey at each point -- populated by
 * DiscreteRankAccumulation (one specific rupture sits at each rank);
 * null for LogNormalCOVAccumulation (every rupture smears across every
 * threshold, so no single rupture is "responsible" for any one point).</p>
 */
public final class LossExceedanceCurve {

    private final double[] losses;      // loss values, ascending
    private final double[] cumulativeLikelihood;  // same order as losses, descending
    private final RuptureKey[] keys;        // nullable array, or null if strategy has no discrete key

    public LossExceedanceCurve(double[] losses, double[] cumulativeLikelihood, RuptureKey[] keys) {
        if (losses.length != cumulativeLikelihood.length)
            throw new IllegalArgumentException("losses/cumulativeRate length mismatch");
        if (keys != null && keys.length != losses.length)
            throw new IllegalArgumentException("keys length mismatch");
        this.losses = losses;
        this.cumulativeLikelihood = cumulativeLikelihood;
        this.keys = keys;
    }

    public int size() {
        return losses.length;
    }

    public double lossAt(int i) {
        return losses[i];
    }

    public double cumulativeLikelihoodAt(int i) {
        return cumulativeLikelihood[i];
    }

    /** May be null if this curve's strategy doesn't identify a discrete rupture per point. */
    public RuptureKey keyAt(int i) {
        return keys == null ? null : keys[i];
    }

    public record LossAtLikelihood(double requestedLikelihood, double loss, double cumulativeLikelihood, RuptureKey key) {}

    /**
     * For each requested exceedance-likelihood, finds the bracketing
     * point: the lowest-loss point on this curve whose cumulativeRate
     * is still >= that level.
     */
    public List<LossAtLikelihood> sampleAtLikelihoods(double[] rateLevels) {
        List<LossAtLikelihood> results = new ArrayList<>(rateLevels.length);

        for (double level : rateLevels) {
            int bracketIndex = -1;
            for (int i = 0; i < losses.length; i++) {
                if (cumulativeLikelihood[i] >= level)
                    bracketIndex = i;
                else
                    break;
            }

            if (bracketIndex < 0) {
                results.add(new LossAtLikelihood(level, Double.NaN, 0d, null));
            } else {
                results.add(new LossAtLikelihood(level, losses[bracketIndex], cumulativeLikelihood[bracketIndex], keyAt(bracketIndex)));
            }
        }

        return results;
    }

    /**
     * The cumulative likelihood of this step-function curve evaluated at
     * an arbitrary loss, not just at one of this curve's own rank points:
     * the cumulativeLikelihood carried by the nearest point whose loss is
     * {@code <= targetLoss} (0 if targetLoss is below this curve's
     * smallest loss) -- the correct value of the step function at any
     * point, not an approximation.
     */
    public double cumulativeLikelihoodAtOrBelow(double targetLoss) {
        int lo = 0, hi = size() - 1, result = -1;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (lossAt(mid) <= targetLoss) {
                result = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return result < 0 ? 0d : cumulativeLikelihoodAt(result);
    }
}
