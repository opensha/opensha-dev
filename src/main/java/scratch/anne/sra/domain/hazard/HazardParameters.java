package scratch.anne.sra.domain.hazard;

import java.nio.file.Path;
import java.util.Map;
import java.util.Objects;

import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;

/**
 * Immutable metadata describing how a hazard field was generated.
 */
public final class HazardParameters {

	public enum HazardMetricType {
	    EXCEEDANCE,
	    OCCURRENCE
	}

	public enum HazardMetric {
	    RATE_EXCEEDANCE(HazardMetricType.EXCEEDANCE),
	    PROBABILITY_EXCEEDANCE(HazardMetricType.EXCEEDANCE),
	    RATE_OCCURENCE(HazardMetricType.OCCURRENCE),
	    PROBABILITY_OCCURENCE(HazardMetricType.OCCURRENCE);

	    private final HazardMetricType metricType;

	    HazardMetric(HazardMetricType metricType) {
	        this.metricType = metricType;
	    }

	    public HazardMetricType getHazardMetricType() {
	        return metricType;
	    }
	}

    private final String erfName;
    private final double erfDuration;
    private final HazardMetric hazardMetric;
    private final String gmmName;

    private final Map<String, String> filterConfig;

    /**
     * Extra ERF construction parameters (e.g. useRupMFDs, includeBackground,
     * griddedSurfaceType, griddedPointSourceMagCutoff), applied generically
     * via Parameter.setValueFromString(String) after the ERF is instantiated
     * and before updateForecast() is called.
     *
     * <p>Empty by default -- an empty map means "build the ERF with no
     * customization," i.e. whatever defaults the ERF class itself provides.
     * ERF customization is strictly opt-in.</p>
     */
    private final Map<String, String> erfConfig;

    /**
     * Path to an external FaultSystemSolution to load and inject into the
     * ERF via setSolution(...), instead of letting the ERF class build its
     * own solution internally.
     *
     * <p>Null (the default) means: instantiate erfName via its own zero-arg
     * constructor and expect it to build/load its own solution (e.g.
     * NSHM23_BranchAveragedERF). Non-null means: erfName is expected to be
     * a BaseFaultSystemSolutionERF-family class that accepts an externally
     * loaded solution.</p>
     *
     * <p>This field's presence is what drives the ERF-construction routing
     * decision in the calculator; HazardParameters itself stays agnostic
     * about *why* -- it just records which route was taken.</p>
     */
    private final Path solutionFile;

    // ------------------------------------------------------------------
    // Convenience constructors (narrowest to widest), all delegating to
    // the full constructor below. Existing call sites with 4 or 5 args
    // keep compiling unchanged, defaulting new fields to empty/null.
    // ------------------------------------------------------------------

    /** No filter config, no erf config, no solution file. */
    public HazardParameters(
            String erfName,
            double erfDuration,
            HazardMetric hazardMetric,
            String gmmName
    ) {
        this(erfName, erfDuration, hazardMetric, gmmName, null, null, null);
    }

    /** Filter config only; no erf config, no solution file. */
    public HazardParameters(
            String erfName,
            double erfDuration,
            HazardMetric hazardMetric,
            String gmmName,
            Map<String, String> filterConfig
    ) {
        this(erfName, erfDuration, hazardMetric, gmmName, filterConfig, null, null);
    }

    /** Filter config + erf config; no solution file. */
    public HazardParameters(
            String erfName,
            double erfDuration,
            HazardMetric hazardMetric,
            String gmmName,
            Map<String, String> filterConfig,
            Map<String, String> erfConfig
    ) {
        this(erfName, erfDuration, hazardMetric, gmmName, filterConfig, erfConfig, null);
    }

    /** Full constructor. */
    public HazardParameters(
            String erfName,
            double erfDuration,
            HazardMetric hazardMetric,
            String gmmName,
            Map<String,String> filterConfig,
            Map<String,String> erfConfig,
            Path solutionFile
    ) {
        this.erfName = Objects.requireNonNull(erfName);
        this.erfDuration = erfDuration;
        this.hazardMetric = Objects.requireNonNull(hazardMetric);
        this.gmmName = Objects.requireNonNull(gmmName);

        this.filterConfig =
                filterConfig == null
                ? Map.of()
                : Map.copyOf(filterConfig);

        this.erfConfig =
                erfConfig == null
                ? Map.of()
                : Map.copyOf(erfConfig);

        this.solutionFile = solutionFile; // null is a valid, expected state

        // -------------------------------
        // Validate ERF class exists
        // -------------------------------
        // Note: this validates that erfName is *a* valid AbstractERF class,
        // regardless of which construction route (solutionFile or not) is
        // used -- BaseFaultSystemSolutionERF-family classes are themselves
        // AbstractERF subclasses, so this check remains correct either way.
        try {
            Class<?> erfClass = Class.forName(erfName);

            if (!org.opensha.sha.earthquake.AbstractERF.class.isAssignableFrom(erfClass)) {
                throw new IllegalArgumentException(
                        "ERF class must extend AbstractERF: " + erfName
                );
            }

        } catch (ClassNotFoundException e) {
            throw new IllegalArgumentException("ERF class not found: " + erfName, e);
        }

        // -------------------------------
        // Validate GMM exists
        // -------------------------------
        try {
            AttenRelRef.valueOf(gmmName.toUpperCase());
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("GMM not recognized: " + gmmName, e);
        }
    }

    public String getErfName() {
        return erfName;
    }

    public double getErfDuration() {
        return erfDuration;
    }

    public HazardMetric getHazardMetric() {
        return hazardMetric;
    }

    public HazardMetricType getHazardMetricType() {
        return hazardMetric.getHazardMetricType();
    }

    public String getGmmName() {
        return gmmName;
    }

    public Map<String,String> getFilterConfig() {
        return filterConfig;
    }

    public Map<String,String> getErfConfig() {
        return erfConfig;
    }

    /** Null means: erfName self-constructs via its own zero-arg constructor. */
    public Path getSolutionFile() {
        return solutionFile;
    }

    @Override
    public String toString() {
        return "HazardParameters{" +
                "erf='" + erfName + '\'' +
                ", duration=" + erfDuration +
                ", hazard metric=" + hazardMetric +
                ", gmm='" + gmmName + '\'' +
                ", solutionFile=" + solutionFile +
                '}';
    }

    public SourceFilterManager buildSourceManager() {

        return new BuildSourceFilterManager()
                .fromConfig(filterConfig)
                .build();
    }
}