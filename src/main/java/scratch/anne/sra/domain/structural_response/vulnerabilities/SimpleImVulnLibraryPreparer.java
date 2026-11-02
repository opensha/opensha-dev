package scratch.anne.sra.domain.structural_response.vulnerabilities;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponse;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.util.ImValueTransformer;

/**
 * Preparer that converts a {@link ResponseModelLibrary} of {@link VulnerabilityModel}
 * directly into a {@link SimpleImResponseLibrary} of {@link SimpleImResponse}.
 *
 * <p>Supports optional filtering by name and optional IM transformation.
 */
public class SimpleImVulnLibraryPreparer {

    private SimpleImVulnLibraryPreparer() {
        // static utility class
    }

    /**
     * Convert all models in the library with a no-op transformer.
     */
    public static SimpleImResponseLibrary prepare(
            ResponseModelLibrary<VulnerabilityModel> modelLibrary) {
        return prepare(modelLibrary, null, new ImValueTransformer.NoImTransformation());
    }

    /**
     * Convert only the models whose names are in {@code filterNames} with a no-op transformer.
     */
    public static SimpleImResponseLibrary prepare(
            ResponseModelLibrary<VulnerabilityModel> modelLibrary,
            Set<String> filterNames) {
        return prepare(modelLibrary, filterNames, new ImValueTransformer.NoImTransformation());
    }

    /**
     * Core preparer that converts models into {@link SimpleImResponse} using a custom transformer.
     *
     * @param modelLibrary source vulnerability library (with metadata)
     * @param filterNames  optional subset of model names (null = all)
     * @param transformer  IM value transformer
     * @return SimpleImResponseLibrary of prepared ExpectedVulnerability objects
     */
    public static SimpleImResponseLibrary prepare(
            ResponseModelLibrary<VulnerabilityModel> modelLibrary,
            Set<String> filterNames,
            ImValueTransformer transformer) {

        // Convert filtered models
        List<SimpleImResponse> responses = modelLibrary.getModels().stream()
                .filter(m -> filterNames == null || filterNames.contains(m.getName()))
                .map(m -> convert(m, transformer))
                .collect(Collectors.toList());

        // Wrap into ExpectedVulnLibrary using metadata from original library
        return SimpleImResponseLibrary.of(responses, modelLibrary.getMetadata());
    }

    // ------------------------------------------------------------------------
    // Conversion helper
    // ------------------------------------------------------------------------

    /**
     * Convert a single VulnerabilityModel into an ExpectedVulnerability
     */
    private static SimpleImResponse convert(
            VulnerabilityModel model,
            ImValueTransformer transformer) {

        double[] imValues = model.getImValues();
        double[] meanDR   = model.getMeanDamageRatio();

        ImValueTransformer.Result result = transformer.transform(imValues, meanDR);

        return new SimpleImResponse(
                model.getName(),
                model.getImt(),
                model.getPeriod(),
                result.imValues,
                result.respValues
        );
    }
}