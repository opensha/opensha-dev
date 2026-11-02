package scratch.anne.sra.domain.hazard;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

import org.opensha.commons.param.Parameter;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;

/**
 * Builds and initializes an AbstractERF from HazardParameters, handling
 * both construction routes:
 *
 * <ul>
 *   <li><b>Self-constructing route</b> (solutionFile == null): erfName is
 *       instantiated via its own zero-arg constructor and is expected to
 *       build/load its own solution internally (e.g.
 *       NSHM23_BranchAveragedERF).</li>
 *   <li><b>Solution-file route</b> (solutionFile != null): a
 *       FaultSystemSolution is loaded from disk and injected via
 *       setSolution(...) into an instance of erfName, which must be a
 *       BaseFaultSystemSolutionERF-family class.</li>
 * </ul>
 *
 * <p>After construction, both routes converge:</p>
 * <ol>
 *   <li>GriddedSeismicitySettings (surface type, finite rupture settings,
 *       mag cutoff) is applied if any of its dedicated keys are present --
 *       this is NOT a simple Parameter, so it can't go through the generic
 *       coercion loop; it's built directly via its own fluent builder.</li>
 *   <li>All remaining erfConfig entries are applied via the parameter's
 *       existing value type (see {@link #coerce}).</li>
 *   <li>Forecast duration is set and updateForecast() is called.</li>
 * </ol>
 */
public final class ErfBuilder {

    // ------------------------------------------------------------------
    // Dedicated erf_config.csv keys for GriddedSeismicitySettings and its
    // nested GriddedFiniteRuptureSettings. Pulled out of erfConfig before
    // the generic Parameter loop runs, same pattern as solutionFile/
    // erfClass are pulled out of erf_config.csv in ParseRiskRunParametersCSV.
    // ------------------------------------------------------------------
    public static final String GRIDDED_SURFACE_TYPE_KEY = "griddedSurfaceType";
    public static final String GRIDDED_POINT_SOURCE_MAG_CUTOFF_KEY = "griddedPointSourceMagCutoff";

    public static final String GRIDDED_FINITE_NUM_SURFACES_KEY = "griddedFiniteNumSurfaces";
    public static final String GRIDDED_FINITE_STRIKE_KEY = "griddedFiniteStrike";
    public static final String GRIDDED_FINITE_SAMPLE_ALONG_STRIKE_KEY = "griddedFiniteSampleAlongStrike";
    public static final String GRIDDED_FINITE_SAMPLE_DOWN_DIP_KEY = "griddedFiniteSampleDownDip";

    private ErfBuilder() {}

    /**
     * Builds and fully initializes an ERF from the given HazardParameters.
     *
     * @throws IllegalArgumentException if the ERF class cannot be
     *         instantiated, does not match the required type for the
     *         chosen route, the solution file cannot be loaded, or any
     *         erfConfig parameter cannot be resolved/parsed/applied.
     */
    public static AbstractERF build(HazardParameters hazardParameters) {

        AbstractERF erf = (hazardParameters.getSolutionFile() != null)
                ? buildFromSolutionFile(hazardParameters)
                : buildSelfConstructing(hazardParameters);

        // Work on a mutable copy: GriddedSeismicitySettings keys are pulled
        // out here so the generic loop below never sees them.
        Map<String, String> remaining = new LinkedHashMap<>(hazardParameters.getErfConfig());

        applyGriddedSeismicitySettingsIfPresent(erf, remaining, hazardParameters.getErfName());

        applyErfConfig(erf, remaining, hazardParameters.getErfName());

        erf.getTimeSpan().setDuration(hazardParameters.getErfDuration());
        erf.updateForecast();

        return erf;
    }

    // ------------------------------------------------------------------
    // Self-constructing route
    // ------------------------------------------------------------------

    private static AbstractERF buildSelfConstructing(HazardParameters hazardParameters) {
        try {
            return (AbstractERF) Class.forName(hazardParameters.getErfName())
                    .getDeclaredConstructor()
                    .newInstance();
        } catch (ClassNotFoundException | NoSuchMethodException | IllegalAccessException
                 | InstantiationException | InvocationTargetException e) {
            throw new IllegalArgumentException(
                    "Failed to instantiate ERF class: " + hazardParameters.getErfName(), e);
        }
    }

    // ------------------------------------------------------------------
    // Solution-file route
    // ------------------------------------------------------------------

    private static AbstractERF buildFromSolutionFile(HazardParameters hazardParameters) {

        Path solutionFile = hazardParameters.getSolutionFile();

        FaultSystemSolution solution;
        try {
            solution = FaultSystemSolution.load(solutionFile.toFile());
        } catch (IOException e) {
            throw new IllegalArgumentException(
                    "Failed to load FaultSystemSolution from: " + solutionFile, e);
        }

        Object erfInstance;
        try {
            erfInstance = Class.forName(hazardParameters.getErfName())
                    .getDeclaredConstructor()
                    .newInstance();
        } catch (ClassNotFoundException | NoSuchMethodException | IllegalAccessException
                 | InstantiationException | InvocationTargetException e) {
            throw new IllegalArgumentException(
                    "Failed to instantiate ERF class: " + hazardParameters.getErfName(), e);
        }

        if (!(erfInstance instanceof BaseFaultSystemSolutionERF fssErf)) {
            throw new IllegalArgumentException(
                    "erf_config.csv's erfClass must be a BaseFaultSystemSolutionERF when "
                    + "solutionFile is set, but got: " + erfInstance.getClass());
        }

        fssErf.setSolution(solution);
        return fssErf;
    }

    // ------------------------------------------------------------------
    // GriddedSeismicitySettings -- not a Parameter, built via its own
    // fluent builder from up to 6 dedicated keys.
    // ------------------------------------------------------------------

    /**
     * If any of the gridded-settings keys are present in erfConfig, builds
     * a GriddedSeismicitySettings starting from GriddedSeismicitySettings
     * .DEFAULT and applies only the fields that were specified, then calls
     * setGriddedSeismicitySettings(...) on the ERF.
     *
     * <p>If none of the keys are present, the ERF's own default
     * GriddedSeismicitySettings is left untouched.</p>
     *
     * <p>Consumes (removes) its keys from erfConfig so the generic
     * Parameter loop that runs afterward never sees them.</p>
     */
    private static void applyGriddedSeismicitySettingsIfPresent(
            AbstractERF erf, Map<String, String> erfConfig, String erfName) {

        String surfaceTypeStr = erfConfig.remove(GRIDDED_SURFACE_TYPE_KEY);
        String magCutoffStr = erfConfig.remove(GRIDDED_POINT_SOURCE_MAG_CUTOFF_KEY);
        String numSurfacesStr = erfConfig.remove(GRIDDED_FINITE_NUM_SURFACES_KEY);
        String strikeStr = erfConfig.remove(GRIDDED_FINITE_STRIKE_KEY);
        String sampleAlongStrikeStr = erfConfig.remove(GRIDDED_FINITE_SAMPLE_ALONG_STRIKE_KEY);
        String sampleDownDipStr = erfConfig.remove(GRIDDED_FINITE_SAMPLE_DOWN_DIP_KEY);

        boolean anySpecified = surfaceTypeStr != null || magCutoffStr != null
                || numSurfacesStr != null || strikeStr != null
                || sampleAlongStrikeStr != null || sampleDownDipStr != null;

        if (!anySpecified)
            return; // leave the ERF's own default GriddedSeismicitySettings untouched

        if (!(erf instanceof BaseFaultSystemSolutionERF fssErf)) {
            throw new IllegalArgumentException(
                    "Gridded seismicity settings were specified in erf_config.csv, but "
                    + erfName + " does not support GriddedSeismicitySettings "
                    + "(not a BaseFaultSystemSolutionERF).");
        }

        GriddedSeismicitySettings settings = GriddedSeismicitySettings.DEFAULT;

        if (surfaceTypeStr != null) {
            try {
                settings = settings.forSurfaceType(BackgroundRupType.valueOf(surfaceTypeStr.trim()));
            } catch (Exception ex) {
                throw new IllegalArgumentException(
                        "Invalid " + GRIDDED_SURFACE_TYPE_KEY + " value: " + surfaceTypeStr, ex);
            }
        }

        if (magCutoffStr != null) {
            try {
                settings = settings.forPointSourceMagCutoff(Double.parseDouble(magCutoffStr.trim()));
            } catch (Exception ex) {
                throw new IllegalArgumentException(
                        "Invalid " + GRIDDED_POINT_SOURCE_MAG_CUTOFF_KEY + " value: " + magCutoffStr, ex);
            }
        }

        boolean anyFiniteSpecified = numSurfacesStr != null || strikeStr != null
                || sampleAlongStrikeStr != null || sampleDownDipStr != null;

        if (anyFiniteSpecified) {
            // Start from DEFAULT_SINGLE (numSurfaces=1, strike=null/random,
            // sampleAlongStrike=false, sampleDownDip=false) and override
            // only the fields that were specified.
            GriddedFiniteRuptureSettings finiteSettings = GriddedFiniteRuptureSettings.DEFAULT_SINGLE;

            if (numSurfacesStr != null) {
                try {
                    finiteSettings = finiteSettings.forNumSurfaces(Integer.parseInt(numSurfacesStr.trim()));
                } catch (Exception ex) {
                    throw new IllegalArgumentException(
                            "Invalid " + GRIDDED_FINITE_NUM_SURFACES_KEY + " value: " + numSurfacesStr, ex);
                }
            }

            if (strikeStr != null) {
                // Explicit "null"/"random" string means: sample randomly
                // (matches the field's own null-means-random semantics).
                if (strikeStr.equalsIgnoreCase("null") || strikeStr.equalsIgnoreCase("random")) {
                    finiteSettings = finiteSettings.forRandStrike();
                } else {
                    try {
                        finiteSettings = finiteSettings.forStrike(Double.parseDouble(strikeStr.trim()));
                    } catch (Exception ex) {
                        throw new IllegalArgumentException(
                                "Invalid " + GRIDDED_FINITE_STRIKE_KEY + " value: " + strikeStr, ex);
                    }
                }
            }

            if (sampleAlongStrikeStr != null) {
                finiteSettings = finiteSettings.forSampleAlongStrike(
                        Boolean.parseBoolean(sampleAlongStrikeStr.trim()));
            }

            if (sampleDownDipStr != null) {
                finiteSettings = finiteSettings.forSampleDownDip(
                        Boolean.parseBoolean(sampleDownDipStr.trim()));
            }

            settings = settings.forFiniteRuptureSettings(finiteSettings);
        }

        fssErf.setGriddedSeismicitySettings(settings);
    }

    // ------------------------------------------------------------------
    // Shared: generic parameter application (everything that IS a
    // simple Parameter<E> -- e.g. Use Rupture MFDs, Background
    // Seismicity, Fault Grid Spacing, Aseismicity Reduces Area).
    // ------------------------------------------------------------------

    /**
     * Applies each remaining erfConfig entry to the ERF via the
     * parameter's existing value type.
     *
     * <p>Parameter&lt;E&gt; has no generic string-setting method (no
     * setValueFromString on the interface) -- only the type-parameterized
     * setValue(E). To stay generic across parameter types without a
     * hardcoded per-name switch, this inspects the RUNTIME TYPE of the
     * parameter's current value (its default, e.g. Boolean.TRUE or an enum
     * constant like IncludeBackgroundOption.INCLUDE) and parses the CSV
     * string into that same type before calling setValue(...) via an
     * unchecked raw-type call.
     *
     * <p>This is inherently a best-effort coercion: it only knows how to
     * parse Boolean, enum, and common Number types. Anything else throws,
     * with the parameter name and offending value included so the failure
     * is traceable back to a specific erf_config.csv entry. Composite/
     * non-Parameter settings (like GriddedSeismicitySettings) must be
     * pulled out and handled separately before this runs -- see
     * {@link #applyGriddedSeismicitySettingsIfPresent}.</p>
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    private static void applyErfConfig(AbstractERF erf, Map<String, String> erfConfig, String erfName) {
        for (Map.Entry<String, String> e : erfConfig.entrySet()) {
            String paramName = e.getKey();
            String rawValue = e.getValue();

            Parameter param;
            try {
                param = erf.getParameter(paramName);
            } catch (Exception ex) {
                throw new IllegalArgumentException(
                        "ERF parameter not found: '" + paramName + "' on " + erfName, ex);
            }

            Object typedValue;
            try {
                typedValue = coerce(param, rawValue);
            } catch (Exception ex) {
                throw new IllegalArgumentException(
                        "Failed to parse value for ERF parameter '" + paramName + "'="
                        + rawValue + " on " + erfName, ex);
            }

            try {
                param.setValue(typedValue);
            } catch (Exception ex) {
                throw new IllegalArgumentException(
                        "Failed to apply ERF parameter '" + paramName + "'=" + rawValue
                        + " (parsed as " + typedValue + ") to " + erfName, ex);
            }
        }
    }

    /**
     * Parses rawValue into whatever type the parameter's current/default
     * value is, so setValue(...) receives an object of a type the
     * parameter's own constraint checking expects.
     */
    @SuppressWarnings("rawtypes")
    private static Object coerce(Parameter param, String rawValue) {

        Object current = param.getValue();
        if (current == null)
            current = param.getDefaultValue();

        if (current == null) {
            throw new IllegalStateException(
                    "Parameter '" + param.getName() + "' has no current or default value "
                    + "to infer a type from; cannot coerce string value: " + rawValue);
        }

        Class<?> targetType = current.getClass();

        if (targetType == Boolean.class) {
            return Boolean.parseBoolean(rawValue.trim());
        }
        if (targetType == Double.class) {
            return Double.parseDouble(rawValue.trim());
        }
        if (targetType == Integer.class) {
            return Integer.parseInt(rawValue.trim());
        }
        if (targetType == Long.class) {
            return Long.parseLong(rawValue.trim());
        }
        if (targetType == String.class) {
            return rawValue;
        }
        if (targetType.isEnum()) {
            @SuppressWarnings("unchecked")
            Object enumValue = Enum.valueOf((Class<? extends Enum>) targetType, rawValue.trim());
            return enumValue;
        }

        throw new IllegalArgumentException(
                "Don't know how to coerce a String into parameter type "
                + targetType.getName() + " for parameter '" + param.getName() + "'. "
                + "Supported types: Boolean, Double, Integer, Long, String, enum. "
                + "Composite types (e.g. GriddedSeismicitySettings) must be handled as a "
                + "dedicated special case, not through generic coercion.");
    }
}