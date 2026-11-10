package scratch.anne.sra.domain.structuralResponse;

import java.util.HashMap;
import java.util.Map;

/**
 * Builds the flat {@code key -> value} reproducibility metadata map for
 * the response/vulnerability model library used by a run -- which
 * vulnerability library JSON file produced the response models this
 * run's assets were convolved against.
 *
 * <p>A pure function of the library object itself, mirroring
 * {@code HazardParametersMetadataBuilder.build(HazardParameters)}: the
 * source file path is read from the library's own
 * {@link LibraryGetters#getMetadata()} (populated by
 * {@code VulnerabilityLibraryReader}/{@code FragilityLibraryReader}'s
 * {@code librarySource} entry), not supplied separately by the caller --
 * so any caller holding a library gets its provenance for free, rather
 * than only the one call site that happened to have the raw file path
 * in scope.</p>
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
     * @param library any response model library; its {@code librarySource}
     *                metadata entry (full path, as written by
     *                {@code VulnerabilityLibraryReader}/{@code FragilityLibraryReader})
     *                is used if present, else an empty string
     */
    public static Map<String, String> build(LibraryGetters<?> library) {
        Map<String, String> meta = new HashMap<>();
        meta.put(RESPONSE_MODEL_PREFIX + "fileName", library.getMetadata().getOrDefault("librarySource", ""));
        return meta;
    }
}
