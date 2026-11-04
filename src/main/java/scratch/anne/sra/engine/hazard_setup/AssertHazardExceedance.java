package scratch.anne.sra.engine.hazard_setup;

import scratch.anne.sra.util.enums.HazardMetric;
import scratch.anne.sra.util.enums.HazardMetricType;

/**
 * Validation-only guard requiring an exceedance hazard metric (rate or
 * probability of exceedance, as opposed to occurrence) for
 * risk-convolution contexts.
 *
 * <p>Shared by {@code PortfolioRiskConvolutionCalculator},
 * {@code PortfolioPerRuptureRiskConvolutionCalculator}, and
 * {@code RuptureProbabilityCurveCollection.fromERF} -- previously each
 * had their own copy of a switch enumerating
 * {@code PROBABILITY_EXCEEDANCE}/{@code RATE_EXCEEDANCE} individually,
 * when the actual invariant is simply
 * {@link HazardMetric#getHazardMetricType()} {@code == EXCEEDANCE}
 * (the same check {@code RuptureLossResultsCollection}/
 * {@code RuptureProbabilityCurveCollection}'s own constructors already
 * make directly, without a switch). Two of the three call sites use
 * this purely as a validation guard (rate vs. probability is handled
 * entirely via the rupture likelihood, not the curve itself, at those
 * two sites); the one call site with a real transform
 * ({@code RATE_EXCEEDANCE -> calc.getAnnualizedRates(...)} in
 * {@code PortfolioRiskConvolutionCalculator}) keeps that transform
 * inline, since it is genuinely distinct logic, not duplicated.</p>
 */
public final class AssertHazardExceedance {

    private AssertHazardExceedance() {}

    /**
     * @throws IllegalArgumentException if {@code hazardMetric} isn't an
     *         exceedance metric (rate or probability of exceedance)
     */
    public static void require(HazardMetric hazardMetric) {
        if (hazardMetric.getHazardMetricType() != HazardMetricType.EXCEEDANCE) {
            throw new IllegalArgumentException(
                    "Hazard metric not supported for risk convolution (must be exceedance): " + hazardMetric);
        }
    }
}
