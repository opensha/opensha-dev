package scratch.anne.sra.domain.portfolio;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/**
 * Builds the flat {@code key -> value} reproducibility metadata map for
 * a portfolio -- which portfolio CSV produced this run's results, and
 * that portfolio's own total asset value (needed to reproduce
 * {@code AssetRiskForRuptureWriter}'s own threshold decisions --
 * {@code fAbsolute * totalPortfolioValue}, see its
 * {@code writeThresholdValue} -- when re-reading an export later).
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
     * @param portfolioCsv        the portfolio CSV this run's assets were read from
     * @param totalPortfolioValue the portfolio's own total asset value
     */
    public static Map<String, String> build(Path portfolioCsv, double totalPortfolioValue) {
        Map<String, String> meta = new HashMap<>();
        meta.put(PORTFOLIO_PREFIX + "fileName", portfolioCsv == null ? "" : portfolioCsv.getFileName().toString());
        meta.put(PORTFOLIO_PREFIX + "totalValue", String.valueOf(totalPortfolioValue));
        return meta;
    }
}
