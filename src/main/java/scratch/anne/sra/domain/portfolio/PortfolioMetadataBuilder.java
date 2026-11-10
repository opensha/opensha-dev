package scratch.anne.sra.domain.portfolio;

import java.util.HashMap;
import java.util.Map;

import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio;

/**
 * Builds the flat {@code key -> value} reproducibility metadata map for
 * a portfolio -- which portfolio CSV produced this run's results, and
 * (when available) that portfolio's own total asset value (needed to
 * reproduce {@code AssetRiskForRuptureWriter}'s own threshold decisions
 * -- {@code fAbsolute * totalPortfolioValue}, see its
 * {@code writeThresholdValue} -- when re-reading an export later).
 *
 * <p>A pure function of the portfolio object itself, mirroring
 * {@code HazardParametersMetadataBuilder.build(HazardParameters)}: the
 * source file path is read from the portfolio's own
 * {@link PortfolioGetters#getMetadata()} (populated by
 * {@code PortfolioReader.readCSV}'s {@code portfolioSource} entry), not
 * supplied separately by the caller -- so any caller holding a
 * portfolio gets its provenance for free, the same way any caller
 * holding a {@code HazardParameters} does, rather than only the one
 * call site that happened to have the raw file path in scope.</p>
 *
 * <p>Deliberately separate from {@code HazardParametersMetadataBuilder}:
 * {@link scratch.anne.sra.domain.hazard.HazardParameters} has
 * no concept of a portfolio or vulnerability file at all (it describes
 * the ERF/hazard setup, not the asset side of a risk calculation), so
 * this metadata doesn't belong folded into that builder -- it's merged
 * in alongside it by callers that have both (e.g.
 * {@code RuptureResultsExporter}), not owned by either.</p>
 */
public final class PortfolioMetadataBuilder {

    public static final String PORTFOLIO_PREFIX = "portfolio.";

    private PortfolioMetadataBuilder() {}

    /**
     * @param portfolio any portfolio; its {@code portfolioSource} metadata
     *                  entry (full path, as written by
     *                  {@code PortfolioReader.readCSV}) is used if present,
     *                  else an empty string
     */
    public static Map<String, String> build(PortfolioGetters<?> portfolio) {
        Map<String, String> meta = new HashMap<>();
        meta.put(PORTFOLIO_PREFIX + "fileName", portfolio.getMetadata().getOrDefault("portfolioSource", ""));
        return meta;
    }

    /**
     * Same as {@link #build(PortfolioGetters)}, plus the portfolio's own
     * total asset value -- only meaningful for an {@code EXPECTED_LOSS}
     * portfolio ({@link RiskConvolutionPortfolio#getTotalPortfolioValue()}
     * throws for a {@code FAILURE_PROBABILITY}/fragility portfolio, which
     * has no concept of asset value). The {@code totalValue} entry is
     * simply omitted for any other risk metric type, rather than
     * throwing, so this overload is always safe to call on any
     * {@link RiskConvolutionPortfolio}.
     *
     * @param portfolio the risk-convolution-ready portfolio to read both
     *                   source path and (if applicable) total value from
     */
    public static Map<String, String> build(RiskConvolutionPortfolio portfolio) {
        Map<String, String> meta = build((PortfolioGetters<?>) portfolio);
        if (portfolio.getRiskMetricType() == RiskConvolutionAsset.RiskMetricType.EXPECTED_LOSS) {
            meta.put(PORTFOLIO_PREFIX + "totalValue", String.valueOf(portfolio.getTotalPortfolioValue()));
        }
        return meta;
    }
}
