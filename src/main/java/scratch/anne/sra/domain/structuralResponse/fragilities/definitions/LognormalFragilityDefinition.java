package scratch.anne.sra.domain.structuralResponse.fragilities.definitions;

import scratch.anne.sra.util.NumericUtil;

/**
 * A fragility curve defined parametrically by a lognormal distribution
 * (median, beta), as opposed to explicit points like
 * {@link DiscreteFragilityDefinition}.
 */
public class LognormalFragilityDefinition implements FragilityDefinition {

    private final double median;
    private final double beta;

    public LognormalFragilityDefinition(double median, double beta) {

        if (Double.isNaN(median) || Double.isInfinite(median) || median <= 0.0)
            throw new IllegalArgumentException(
                    "Median must be a positive finite value");

        if (Double.isNaN(beta) || Double.isInfinite(beta) || beta <= 0.0)
            throw new IllegalArgumentException(
                    "Beta must be a positive finite value");

        this.median = median;
        this.beta = beta;
    }

    // -------- Getters --------
    public double getMedian() { return median; }
    public double getBeta() { return beta; }

    // -------- Probability computation --------

    /**
     * Compute lognormal probability for a single IM value.
     * Returns P(Damage ≥ IM) using lognormal CDF.
     */
    public double getProbability(double imValue) {
        if (imValue <= 0.0) return 0.0;
        double lnRatio = Math.log(imValue / median);
        return NumericUtil.standardNormalCDF(lnRatio / beta);
    }

    /**
     * Compute lognormal probabilities for an array of IM values.
     * Returns array of P(Damage ≥ IM).
     */
    public double[] getProbability(double[] imValues) {
        double[] probs = new double[imValues.length];
        for (int i = 0; i < imValues.length; i++) {
            probs[i] = getProbability(imValues[i]);
        }
        return probs;
    }

    @Override
    public String toString() {
        return "LognormalFragilityDefinition{median=" + median + ", beta=" + beta + "}";
    }

    public String toVerboseString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Lognormal Fragility Definition\n");
        sb.append(String.format("%-10s : %.6f%n", "Median", median));
        sb.append(String.format("%-10s : %.6f%n", "Beta", beta));
        return sb.toString();
    }
}