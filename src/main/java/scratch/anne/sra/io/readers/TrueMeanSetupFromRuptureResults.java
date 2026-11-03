package scratch.anne.sra.io.readers;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.PortfolioMetadataBuilder;
import scratch.anne.sra.domain.structural_response.ResponseModelMetadataBuilder;
import scratch.anne.sra.util.enums.HazardMetric;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Reconstructs the full TrueMeanERFDemo setup (everything through
 * updateForecast() and the grid-source lookup table) from a
 * RuptureResultsExporter CSV, instead of hardcoding the construction
 * sequence directly in a main() method.
 *
 * <p>This reads the same "# HazardParameters" key,value block that
 * RuptureResultsExporter.exportCsv(...) writes (erfName, erfDuration,
 * hazardMetric, gmmName, solutionFile, erfConfig.*, filterConfig.*, plus
 * any portfolio.* /responseModel.* provenance the caller supplied --
 * see PortfolioMetadataBuilder/ResponseModelMetadataBuilder), parses it
 * back into a HazardParameters, and drives ErfBuilder from that -- so
 * the reconstructed ERF is built through the exact same code path as
 * the original run, rather than a hand-copied duplicate of it. The full
 * raw block is also exposed via {@link TrueMeanSetup#rawMetadata()} (and
 * {@link TrueMeanSetup#portfolioAndResponseModelMetadata()} for just the
 * non-HazardParameters fields), so downstream pipeline stages (e.g.
 * BranchContributionsStage) can recover portfolio/response-model
 * provenance without any new file I/O of their own.</p>
 */
public final class TrueMeanSetupFromRuptureResults {

    private TrueMeanSetupFromRuptureResults() {}

    /** Everything the demo's branch-mapping loop needs, bundled together. */
    public record TrueMeanSetup(
            HazardParameters hazardParameters,
            FaultSystemSolution trueMeanSol,
            BaseFaultSystemSolutionERF trueMeanERF,
            GridSourceList meanGridList,
            Table<Integer, TectonicRegionType, Integer> meanLocToSourceTable,
            float gridMinMag,
            Map<String, String> rawMetadata
    ) {

        /**
         * Convenience accessor for the {@code portfolio.*}/
         * {@code responseModel.*} subset of {@link #rawMetadata()} --
         * see {@code PortfolioMetadataBuilder}/
         * {@code ResponseModelMetadataBuilder}, which both produce
         * their own prefixed keys into the same exported CSV's
         * metadata block this record was reconstructed from. Keys are
         * returned with their prefix intact (e.g. {@code portfolio.fileName}),
         * so callers merging this into another metadata map (e.g.
         * {@code BranchContributionsStage}'s own Parquet footer) don't
         * need to know the prefix constants themselves.
         */
        public Map<String, String> portfolioAndResponseModelMetadata() {
            Map<String, String> out = new java.util.HashMap<>();
            for (Map.Entry<String, String> e : rawMetadata.entrySet()) {
                String key = e.getKey();
                if (key.startsWith(PortfolioMetadataBuilder.PORTFOLIO_PREFIX)
                        || key.startsWith(ResponseModelMetadataBuilder.RESPONSE_MODEL_PREFIX)) {
                    out.put(key, e.getValue());
                }
            }
            return out;
        }
    }

    /**
     * Parses the "# HazardParameters" block at the top of an exported
     * RuptureResultsExporter CSV file (stops at the first blank line,
     * matching exportCsv's own layout) into a plain key,value map.
     */
    private static Map<String, String> parseHazardParametersBlock(Path exportedCsv) throws IOException {
        Map<String, String> raw = new LinkedHashMap<>();

        try (BufferedReader br = Files.newBufferedReader(exportedCsv)) {
            String line = br.readLine(); // "# HazardParameters" header line
            if (line == null || !line.startsWith("#")) {
                throw new IllegalArgumentException(
                        "Expected '# HazardParameters' header as first line of " + exportedCsv);
            }

            while ((line = br.readLine()) != null) {
                if (line.isBlank())
                    break; // exportCsv writes a blank line after the HazardParameters block

                String[] tokens = line.split(",", 2);
                if (tokens.length == 2) {
                    raw.put(tokens[0].trim(), tokens[1].trim());
                }
            }
        }

        return raw;
    }

    /**
     * Reconstructs HazardParameters from an exported CSV's metadata block.
     *
     * <p>Inverts RuptureResultsExporter.buildMetadata(...): un-prefixes
     * "erfConfig." and "filterConfig." keys back into their own maps, and
     * treats an empty solutionFile value as null (self-constructing
     * route), matching how buildMetadata wrote it out.</p>
     */
    public static HazardParameters parseHazardParameters(Path exportedCsv) throws IOException {

        Map<String, String> raw = parseHazardParametersBlock(exportedCsv);

        String erfName = require(raw, "erfName", exportedCsv);
        double erfDuration = Double.parseDouble(require(raw, "erfDuration", exportedCsv));
        HazardMetric hazardMetric =
                HazardMetric.valueOf(require(raw, "hazardMetric", exportedCsv));
        String gmmName = require(raw, "gmmName", exportedCsv);

        String solutionFileStr = raw.get("solutionFile");
        Path solutionFile = (solutionFileStr == null || solutionFileStr.isBlank())
                ? null
                : Path.of(solutionFileStr);

        Map<String, String> erfConfig = new HashMap<>();
        Map<String, String> filterConfig = new HashMap<>();

        for (Map.Entry<String, String> e : raw.entrySet()) {
            String key = e.getKey();
            if (key.startsWith("erfConfig.")) {
                erfConfig.put(key.substring("erfConfig.".length()), e.getValue());
            } else if (key.startsWith("filterConfig.")) {
                filterConfig.put(key.substring("filterConfig.".length()), e.getValue());
            }
            // fixed fields (erfName, erfDuration, hazardMetric, gmmName,
            // solutionFile, totalExpectedLoss) are already handled above
            // or are informational-only and not needed to rebuild the ERF.
        }

        return new HazardParameters(
                erfName,
                erfDuration,
                hazardMetric,
                gmmName,
                filterConfig,
                erfConfig,
                solutionFile
        );
    }

    private static String require(Map<String, String> raw, String key, Path source) {
        String v = raw.get(key);
        if (v == null || v.isBlank()) {
            throw new IllegalArgumentException(
                    "Missing required field '" + key + "' in exported HazardParameters block: " + source);
        }
        return v;
    }

    /**
     * Full reconstruction: parses HazardParameters from the exported CSV,
     * builds the ERF via ErfBuilder (throwing if the file's erfName isn't
     * a BaseFaultSystemSolutionERF, since the grid-source lookup table
     * below requires that specific type), and rebuilds the
     * meanGridList/meanLocToSourceTable exactly as TrueMeanERFDemo does.
     *
     * <p>Equivalent to TrueMeanERFDemo's setup through updateForecast()
     * (its line ~112), but driven by exported metadata instead of a
     * hardcoded construction sequence.</p>
     */
    public static TrueMeanSetup rebuild(Path exportedCsv) throws IOException {

        HazardParameters hazardParameters = parseHazardParameters(exportedCsv);
        Map<String, String> rawMetadata = parseHazardParametersBlock(exportedCsv);

        if (hazardParameters.getSolutionFile() == null) {
            throw new IllegalArgumentException(
                    "Exported HazardParameters has no solutionFile -- this reconstruction path "
                    + "is only valid for runs built via the solution-file route (see ErfBuilder). "
                    + "Source: " + exportedCsv);
        }

        var erf = scratch.anne.sra.domain.hazard.ErfBuilder.build(hazardParameters);

        if (!(erf instanceof BaseFaultSystemSolutionERF trueMeanERF)) {
            throw new IllegalStateException(
                    "Reconstructed ERF is not a BaseFaultSystemSolutionERF (got "
                    + erf.getClass() + "); cannot recover a FaultSystemSolution reference from it.");
        }

        // ErfBuilder doesn't hand back the FaultSystemSolution it loaded
        // internally, but BaseFaultSystemSolutionERF exposes it back out
        // via getSolution(), same solution instance that was passed to
        // setSolution(...) inside ErfBuilder.
        FaultSystemSolution trueMeanSol = trueMeanERF.getSolution();

        // --- Mirrors TrueMeanERFDemo lines ~44-51 ---
        // The pre-filter mag cutoff applied to meanGridList is assumed to
        // match griddedPointSourceMagCutoff (the same cutoff ErfBuilder
        // applies to the ERF's own GriddedSeismicitySettings), rather than
        // being tracked as a separate reproducibility field. Falls back to
        // 5f (the demo's original hardcoded value) if that key wasn't
        // present in erfConfig, e.g. for exports predating this field.
        double gridMinMagDouble = hazardParameters.getErfConfig()
                .entrySet().stream()
                .filter(e -> e.getKey().equals(
                        scratch.anne.sra.domain.hazard.ErfBuilder.GRIDDED_POINT_SOURCE_MAG_CUTOFF_KEY))
                .map(e -> Double.parseDouble(e.getValue()))
                .findFirst()
                .orElse(5d);
        float gridMinMag = (float) gridMinMagDouble;

        GridSourceList meanGridList = trueMeanSol.requireModule(GridSourceList.class);
        meanGridList = meanGridList.getAboveMinMag(gridMinMag);
        trueMeanSol.setGridSourceProvider(meanGridList);

        Table<Integer, TectonicRegionType, Integer> meanLocToSourceTable = HashBasedTable.create();
        for (int s = 0; s < meanGridList.getNumSources(); s++) {
            int locIndex = meanGridList.getLocationIndexForSource(s);
            TectonicRegionType trt = meanGridList.tectonicRegionTypeForSourceIndex(s);
            meanLocToSourceTable.put(locIndex, trt, s);
        }

        return new TrueMeanSetup(
                hazardParameters,
                trueMeanSol,
                trueMeanERF,
                meanGridList,
                meanLocToSourceTable,
                gridMinMag,
                rawMetadata
        );
    }
}