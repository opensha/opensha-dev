package scratch.anne.sra.domain.structural_response.fragilities.definitions;

/**
 * A fragility curve defined by explicit (IM level, probability) points,
 * as opposed to a parametric form like {@link LognormalFragilityDefinition}.
 */
public class DiscreteFragilityDefinition implements FragilityDefinition {

    private final double[] imLevels;
    private final double[] probability;

    public DiscreteFragilityDefinition(double[] imLevels, double[] probability) {

        if (imLevels == null || probability == null)
            throw new IllegalArgumentException("IM levels and probability arrays cannot be null");

        if (imLevels.length == 0)
            throw new IllegalArgumentException("Arrays cannot be empty");

        if (imLevels.length != probability.length)
            throw new IllegalArgumentException(
                    "IM levels and probability arrays must have equal length");

        this.imLevels = imLevels.clone();
        this.probability = probability.clone();

        validate();
    }

    private void validate() {

        // IM levels must be strictly increasing
        for (int i = 1; i < imLevels.length; i++) {
            if (imLevels[i] <= imLevels[i - 1])
                throw new IllegalStateException(
                        "IM levels must be strictly increasing");
        }

        // probabilities must lie in [0,1]
        for (double p : probability) {
            if (p < 0.0 || p > 1.0)
                throw new IllegalStateException(
                        "Probabilities must be between 0 and 1");
        }

        // fragility exceedance curves should be non-decreasing
        for (int i = 1; i < probability.length; i++) {
            if (probability[i] < probability[i - 1])
                throw new IllegalStateException(
                        "Fragility probabilities must be non-decreasing with IM");
        }
    }

    // -------- Getters --------

    public double[] getImLevels() {
        return imLevels.clone();
    }

    public double[] getProbability() {
        return probability.clone();
    }

    public int size() {
        return imLevels.length;
    }

    @Override
    public String toString() {
        return "DiscreteFragilityDefinition{points=" + imLevels.length + "}";
    }

    public String toVerboseString() {

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%-12s %-12s%n", "IML", "ProbExceed"));

        for (int i = 0; i < imLevels.length; i++) {
            sb.append(String.format(
                    "%-12.6f %-12.6f%n",
                    imLevels[i],
                    probability[i]));
        }

        return sb.toString();
    }
}
