package scratch.anne.sra.domain.hazard;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/**
 * Builds the standard flat {@code key -> value} reproducibility metadata
 * map for a {@link HazardParameters} instance -- erfName, erfDuration,
 * hazardMetric, gmmName, solutionFile, plus every {@code erfConfig.*}/
 * {@code filterConfig.*} entry.
 *
 * <p>Shared by every writer that needs to record "exactly which ERF/
 * hazard setup produced this file" -- originally only
 * {@code RuptureResultsExporter} (the {@code _eal-per-rup.csv} metadata
 * block), extracted here so
 * {@code LogicTreeBranchRuptureContributionExporter}'s own Parquet
 * file-level metadata (see its javadoc) can record the IDENTICAL
 * HazardParameters provenance -- the same {@link HazardParameters}
 * instance reconstructed via {@code TrueMeanSetupFromRuptureResults.rebuild}
 * from that same {@code _eal-per-rup.csv} -- without duplicating this
 * field list by hand in a second place.</p>
 */
public final class HazardParametersMetadataBuilder {

    public static final String ERF_CONFIG_PREFIX = "erfConfig.";
    public static final String FILTER_CONFIG_PREFIX = "filterConfig.";

    private HazardParametersMetadataBuilder() {}

    /**
     * @param hp hazard parameters to flatten
     */
    public static Map<String, String> build(HazardParameters hp) {

        Map<String, String> meta = new HashMap<>();
        meta.put("erfName",     hp.getErfName());
        meta.put("erfDuration", String.valueOf(hp.getErfDuration()));
        meta.put("hazardMetric", hp.getHazardMetric().name());
        meta.put("gmmName",      hp.getGmmName());

        // solutionFile: null means the self-constructing route was used
        // (erfName built its own solution internally). Recorded
        // explicitly as an empty string rather than omitted, so the
        // absence itself is unambiguous on read-back rather than looking
        // like a missing/forgotten field.
        Path solutionFile = hp.getSolutionFile();
        meta.put("solutionFile", solutionFile == null ? "" : solutionFile.toString());

        // erfConfig/filterConfig: prefixed so these keys can't collide
        // with each other or with the fixed fields above.
        Map<String, String> erfConfig = hp.getErfConfig();
        if (erfConfig != null) {
            for (Map.Entry<String, String> e : erfConfig.entrySet())
                meta.put(ERF_CONFIG_PREFIX + e.getKey(), e.getValue());
        }

        Map<String, String> filterConfig = hp.getFilterConfig();
        if (filterConfig != null) {
            for (Map.Entry<String, String> e : filterConfig.entrySet())
                meta.put(FILTER_CONFIG_PREFIX + e.getKey(), e.getValue());
        }

        return meta;
    }
}
