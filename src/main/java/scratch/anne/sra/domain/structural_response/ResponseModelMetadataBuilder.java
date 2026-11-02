package scratch.anne.sra.domain.structural_response;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/**
 * Builds the flat {@code key -> value} reproducibility metadata map for
 * the response/vulnerability model library used by a run -- which
 * vulnerability library JSON file produced the response models this
 * run's assets were convolved against.
 *
 * <p>Deliberately separate from {@code HazardParametersMetadataBuilder}
 * (describes the ERF/hazard setup only) and
 * {@code PortfolioMetadataBuilder} (describes the asset side) -- this
 * is its own, third, independent piece of provenance, merged in
 * alongside the other two by callers that have all three (e.g.
 * {@code RuptureResultsExporter}).</p>
 */
public final class ResponseModelMetadataBuilder {

    public static final String RESPONSE_MODEL_PREFIX = "responseModel.";

    private ResponseModelMetadataBuilder() {}

    /**
     * @param vulnLibraryJson the vulnerability library JSON this run's response models were read from
     */
    public static Map<String, String> build(Path vulnLibraryJson) {
        Map<String, String> meta = new HashMap<>();
        meta.put(RESPONSE_MODEL_PREFIX + "fileName",
                vulnLibraryJson == null ? "" : vulnLibraryJson.getFileName().toString());
        return meta;
    }
}
