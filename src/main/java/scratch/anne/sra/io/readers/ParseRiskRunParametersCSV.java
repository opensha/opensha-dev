package scratch.anne.sra.io.readers;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.HazardParameters.HazardMetric;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.util.StringUtil;

public final class ParseRiskRunParametersCSV {

    private ParseRiskRunParametersCSV() {}

    private static final Set<String> RESERVED_KEYS = Set.of(
            "portfolio_csv",
            "vuln_library_json",
            "fragility_library_json",
            "file_tag",
            "log_im_step",
            "integration_method",
            "erf_class",
            "erf_duration",
            "gmm",
            "hazard_metric",
            "probability_targets"
    );

    /** Within an erf_config file: path to an external FaultSystemSolution to load. */
    public static final String SOLUTION_FILE_KEY = "solutionFile";

    /**
     * Within an erf_config file: the ERF class to wrap the loaded solution
     * with. Required if solutionFile is present. Deliberately independent
     * from input.csv's erf_class -- a solution-file ERF (e.g.
     * BaseFaultSystemSolutionERF) and a self-constructing ERF (e.g.
     * NSHM23_BranchAveragedERF) are architecturally unrelated and not
     * interchangeable, so the erf_config file must be fully self-contained.
     */
    public static final String ERF_CLASS_KEY = "erfClass";

    /**
     * Loads a generic two-column (key,value) CSV into a map.
     * Reused for both the main run config and any erf_config file.
     */
    public static Map<String, String> loadConfig(Path csvPath)
            throws IOException {

        Map<String, String> map = new HashMap<>();

        try (BufferedReader br = Files.newBufferedReader(csvPath)) {

            br.readLine(); // header

            String line;
            while ((line = br.readLine()) != null) {

                String[] tokens = line.split(",", 2);

                if (tokens.length == 2) {
                    map.put(
                            tokens[0].trim(),
                            tokens[1].trim()
                    );
                }
            }
        }

        return map;
    }

    /**
     * Holds the parsed contents of an erf_config file, when input.csv's
     * erf_class points at one rather than naming a class directly.
     *
     * <p>Two mutually exclusive construction routes:</p>
     * <ul>
     *   <li>{@code usesSolutionFile() == false}: input.csv's erf_class is
     *       an actual class name (no erf_config file involved at all).
     *       The named ERF is instantiated via its own zero-arg
     *       constructor and is expected to build/load its own solution
     *       internally (e.g. NSHM23_BranchAveragedERF).</li>
     *   <li>{@code usesSolutionFile() == true}: input.csv's erf_class
     *       pointed at an erf_config file, which supplies BOTH
     *       solutionFile and its own erfClass (e.g.
     *       BaseFaultSystemSolutionERF). That erfClass governs
     *       construction instead of input.csv's erf_class value.</li>
     * </ul>
     *
     * <p>{@code params} holds all remaining erf_config entries, applied
     * generically via {@code Parameter.setValue(...)} (see ErfBuilder)
     * after construction and before {@code updateForecast()} (e.g.
     * useRupMFDs, includeBackground, griddedSurfaceType).</p>
     */
    public record ErfConfig(Path solutionFile, String erfClass, Map<String, String> params) {

        public static final ErfConfig EMPTY = new ErfConfig(null, null, Map.of());

        public boolean usesSolutionFile() {
            return solutionFile != null;
        }
    }

    /**
     * Determines whether input.csv's erf_class value should be treated as
     * a path to an erf_config file, rather than a Java class name.
     *
     * <p>Rule: if it ends in ".csv" (case-insensitive), it's a path.
     * Class names never have a file extension, so this is unambiguous
     * and doesn't require the file to exist yet at check time.</p>
     */
    public static boolean isErfConfigPath(String erfClassValue) {
        return erfClassValue != null
                && erfClassValue.trim().toLowerCase().endsWith(".csv");
    }

    /**
     * Resolves an erf_class value that isErfConfigPath(...) against the
     * run folder: absolute paths are used as-is, relative paths (e.g.
     * just "erf_config.csv") are resolved relative to runFolder, so the
     * file can sit right next to input.csv without needing a hardcoded
     * fixed filename.
     */
    public static Path resolveErfConfigPath(Path runFolder, String erfClassValue) {
        Path p = Path.of(erfClassValue.trim());
        Path resolved = p.isAbsolute() ? p : runFolder.resolve(p);
        return resolved.normalize();
    }

    /**
     * Loads the erf_config, if input.csv's erf_class points at one.
     *
     * <p>If erf_class is an actual class name (not a .csv path), returns
     * {@link ErfConfig#EMPTY} -- the self-constructing route is used with
     * no extra parameters, exactly as if this feature didn't exist.</p>
     *
     * @throws IllegalArgumentException if erf_class looks like a path but
     *         the file doesn't exist, or if the referenced file specifies
     *         solutionFile without an accompanying erfClass.
     */
    public static ErfConfig loadErfConfig(Path runFolder, Map<String, String> config) throws IOException {

        String erfClassValue = config.get("erf_class");

        if (!isErfConfigPath(erfClassValue))
            return ErfConfig.EMPTY;

        Path erfConfigPath = resolveErfConfigPath(runFolder, erfClassValue);

        if (!Files.exists(erfConfigPath)) {
            throw new IllegalArgumentException(
                    "input.csv's erf_class ('" + erfClassValue + "') looks like a path to an "
                    + "erf_config file (ends in .csv) but no file was found at: "
                    + erfConfigPath.toAbsolutePath());
        }

        Map<String, String> raw = loadConfig(erfConfigPath);

        String solutionFileStr = raw.remove(SOLUTION_FILE_KEY);
        String erfClass = raw.remove(ERF_CLASS_KEY);

        Path solutionFile = (solutionFileStr == null || solutionFileStr.isBlank())
                ? null
                : Path.of(solutionFileStr);

        if (solutionFile != null && (erfClass == null || erfClass.isBlank())) {
            throw new IllegalArgumentException(
                    erfConfigPath + " specifies '" + SOLUTION_FILE_KEY + "' but is missing '"
                    + ERF_CLASS_KEY + "' -- the solution-file route requires its own ERF "
                    + "class (e.g. a BaseFaultSystemSolutionERF-family class), independent "
                    + "of input.csv's erf_class."
            );
        }

        return new ErfConfig(solutionFile, erfClass, Map.copyOf(raw));
    }

    /**
     * Extracts all non-reserved entries from the main run config, for
     * SourceFilterManager configuration only.
     *
     * <p>Note: ERF construction parameters live entirely in the separate
     * erf_config file, if any (see {@link #loadErfConfig}), not in this
     * file, so there is no key collision to worry about here.</p>
     *
     * Examples:
     * max_distance=200
     * min_mag=5.0
     */
    public static Map<String, String> extractSourceFilterConfig(
            Map<String, String> config) {

        Map<String, String> filterConfig = new HashMap<>();

        for (var e : config.entrySet()) {
            if (!RESERVED_KEYS.contains(e.getKey())) {
                filterConfig.put(e.getKey(), e.getValue());
            }
        }

        return filterConfig;
    }

    /**
     * Builds HazardParameters, resolving the ERF construction route.
     *
     * <p>Routing rule:</p>
     * <ul>
     *   <li>erfConfig.usesSolutionFile() == true: erfConfig's own erfClass
     *       governs construction; input.csv's erf_class (the path to the
     *       erf_config file) is not itself a class name and is not used
     *       for construction.</li>
     *   <li>otherwise: input.csv's erf_class is a real class name and
     *       governs construction directly, as before.</li>
     * </ul>
     *
     * <p>The resolved class name is stored as HazardParameters.erfName
     * either way, so downstream code (and export metadata) never needs
     * to know which route was taken -- getErfName() always reflects the
     * class actually instantiated.</p>
     */
    public static HazardParameters buildHazardParameters(
            Map<String, String> config,
            ErfConfig erfConfig) {

        String gmm = config.get("gmm");
        if (gmm == null || gmm.isBlank()) {
            throw new IllegalArgumentException("Missing required key: gmm");
        }

        String effectiveErfClass;
        if (erfConfig.usesSolutionFile()) {
            // already validated non-blank in loadErfConfig
            effectiveErfClass = erfConfig.erfClass();
        } else {
            String inputCsvErfClass = config.get("erf_class");
            if (inputCsvErfClass == null || inputCsvErfClass.isBlank()) {
                throw new IllegalArgumentException("Missing required key: erf_class");
            }
            if (isErfConfigPath(inputCsvErfClass)) {
                // Defensive: loadErfConfig should have already turned this
                // into a usesSolutionFile()==true ErfConfig if the file
                // was valid. Reaching here means the file existed but had
                // no solutionFile in it -- ambiguous, so fail loudly
                // rather than silently trying to instantiate a ".csv" as
                // a class name.
                throw new IllegalArgumentException(
                        "input.csv's erf_class ('" + inputCsvErfClass + "') points at an "
                        + "erf_config file, but that file has no '" + SOLUTION_FILE_KEY
                        + "' entry, so there is no ERF class to construct. Either add "
                        + SOLUTION_FILE_KEY + "/" + ERF_CLASS_KEY + " to that file, or set "
                        + "input.csv's erf_class to an actual class name.");
            }
            effectiveErfClass = inputCsvErfClass;
        }

        double erfDuration =
                StringUtil.parseDoubleOrDefault(
                        config.get("erf_duration"),
                        1.0
                );

        HazardMetric hazardMetric =
                parseHazardMetric(
                        config.get("hazard_metric")
                );

        Map<String, String> filterConfig =
                extractSourceFilterConfig(config);

        return new HazardParameters(
                effectiveErfClass,
                erfDuration,
                hazardMetric,
                gmm,
                filterConfig,
                erfConfig.params(),
                erfConfig.solutionFile()
        );
    }

    public static HazardMetric parseHazardMetric(String s) {

        if (s == null || s.isBlank()) {
            return HazardMetric.PROBABILITY_EXCEEDANCE;
        }

        String key = s.trim().toUpperCase();

        switch (key) {

            case "PROB":
            case "PROBABILITY":
            case "PROBABILITY_EXCEEDANCE":
                return HazardMetric.PROBABILITY_EXCEEDANCE;

            case "RATE":
            case "RATE_EXCEEDANCE":
                return HazardMetric.RATE_EXCEEDANCE;

            default:
                return HazardMetric.valueOf(key);
        }
    }

    public static RiskConvolution.IntegrationMethod parseIntegrationMethod(String s) {
        if (s == null || s.isBlank()) {
            return RiskConvolution.IntegrationMethod.CLOSED_FORM;
        }

        return RiskConvolution.IntegrationMethod.valueOf(s.trim().toUpperCase());
    }
}