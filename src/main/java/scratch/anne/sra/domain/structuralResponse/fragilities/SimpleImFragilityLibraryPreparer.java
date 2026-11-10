package scratch.anne.sra.domain.structuralResponse.fragilities;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.structuralResponse.*;
import scratch.anne.sra.domain.structuralResponse.fragilities.definitions.*;
import scratch.anne.sra.util.ImValueTransformer;
import scratch.anne.sra.util.NumericUtil;
import scratch.anne.sra.util.enums.LimitState;

/**
 * Converts a {@link ResponseModelLibrary} of {@link FragilityModel}
 * into a {@link SimpleImResponseLibrary}.
 *
 * <h2>Design Overview</h2>
 *
 * Preparation is separated into two conceptual stages:
 *
 * <ul>
 *   <li><b>{@link PrepareSpec}</b> — user intent specification
 *       (filtering, limit-state selection, IM grid policy)</li>
 *
 *   <li><b>{@link ImConfiguration}</b> — execution plan
 *       describing how IM arrays are transformed</li>
 * </ul>
 *
 * <p>The specification describes <i>what</i> should happen.
 * The configuration describes <i>how</i> to execute it efficiently.
 *
 * <p>This separation prevents repeated configuration logic during
 * per-fragility conversion and guarantees consistent IM treatment
 * across discrete and lognormal fragilities.
 */
public final class SimpleImFragilityLibraryPreparer {

    private SimpleImFragilityLibraryPreparer() {}

    /**
     * Prepare an IM-response library from a fragility library.
     *
     * @param library fragility models
     * @param spec preparation specification
     * @return discretized IM-response library
     */
    public static SimpleImResponseLibrary prepare(
            ResponseModelLibrary<FragilityModel> library,
            PrepareSpec spec) {

        // ---------------- filtering (SPEC responsibility) ----------------
        List<FragilityModel> selectedModels =
                spec.selectModels(library.getModels());

        // ---------------- resolve IM configuration ----------------
        ImConfiguration imConfig =
                spec.resolve();

        // ---------------- convert ----------------
        List<SimpleImResponse> responses =
                selectedModels.stream()
                        .map(m -> convert(m, spec, imConfig))
                        .collect(Collectors.toList());

        return SimpleImResponseLibrary.of(
                responses,
                library.getMetadata());
    }

    // =====================================================================
    // CONVERSION
    // =====================================================================

    /**
     * Converts a single fragility model into a discretized IM response.
     */
    private static SimpleImResponse convert(
            FragilityModel model,
            PrepareSpec spec,
            ImConfiguration imConfig) {

        LimitStateFragility ls =
                spec.selectLimitState(model);

        double[] imValues;
        double[] respValues;

        // =========================================================
        // DISCRETE FRAGILITY
        // =========================================================
        if (ls.getDefinition() instanceof DiscreteFragilityDefinition def) {

            var result = imConfig.transformer.transform(
                    def.getImLevels(),
                    def.getProbability()
            );

            imValues = result.imValues;
            respValues = result.respValues;
        }

        // =========================================================
        // LOGNORMAL FRAGILITY
        // =========================================================
        else if (ls.getDefinition() instanceof LognormalFragilityDefinition def) {

            imValues = imConfig.resolveImValues(def.getMedian());
            respValues = def.getProbability(imValues);
        }

        else {
            throw new IllegalStateException(
                    "Unknown FragilityDefinition type: "
                            + ls.getDefinition().getClass());
        }

        return new SimpleImResponse(
                model.getName(),
                model.getImt(),
                model.getPeriod(),
                imValues,
                respValues
        );
    }

    // =====================================================================
    // RESOLVED IM Configuration
    // =====================================================================

    /**
     * Execution plan derived from {@link PrepareSpec}.
     *
     * <p>This object contains only mechanics required during
     * per-fragility execution. It intentionally contains
     * <b>no filtering logic</b>.
     */
    static final class ImConfiguration {

        final NumericUtil.ImArrayParam imParam;
        final double[] customImValues;
        final ImValueTransformer transformer;

        ImConfiguration(
                NumericUtil.ImArrayParam imParam,
                double[] customImValues,
                ImValueTransformer transformer) {

            this.imParam = imParam;
            this.customImValues = customImValues;
            this.transformer = transformer;
        }

        /**
         * Resolve IM values for a lognormal fragility.
         */
        double[] resolveImValues(double median) {

            if (customImValues != null) {
                return customImValues;
            }

            return NumericUtil.distributeImAroundMedian(
                    median,
                    imParam
            );
        }
    }

    // =====================================================================
    // PREPARE SPEC
    // =====================================================================

    /**
     * Immutable configuration describing how fragility models
     * should be prepared.
     *
     * <p>This object represents user intent only.
     */
    public static final class PrepareSpec {

        private final Set<String> fragilityNames;
        private final LimitState selectedLS;
        private final NumericUtil.ImArrayParam imParam;
        private final double[] customImValues;
        private final ImValueTransformer transformer;

        private PrepareSpec(Builder b) {
            this.fragilityNames = b.fragilityNames;
            this.selectedLS = b.selectedLS;
            this.imParam = b.imParam;
            this.customImValues = b.customImValues;
            this.transformer = b.transformer;
        }

        // ---------------------------------------------------------
        // FILTERING
        // ---------------------------------------------------------

        /**
         * Applies name filtering to a collection of fragility models.
         */
        List<FragilityModel> selectModels(
                List<FragilityModel> models) {

            if (fragilityNames == null)
                return models;

            return models.stream()
                    .filter(m -> fragilityNames.contains(m.getName()))
                    .collect(Collectors.toList());
        }

        /**
         * Selects the limit state for a model.
         */
        LimitStateFragility selectLimitState(
                FragilityModel model) {

            if (selectedLS == null)
                return model.getPrimaryLimitState();

            return model.getLimitStateFragilities().stream()
                    .filter(ls -> ls.getLimitState() == selectedLS)
                    .findFirst()
                    .orElseThrow(() ->
                            new IllegalArgumentException(
                                    "LimitState not found: " + selectedLS));
        }

        // ---------------------------------------------------------
        // RESOLUTION
        // ---------------------------------------------------------

        /**
         * Resolves this specification into an executable im configuration
         */
        ImConfiguration resolve() {

            NumericUtil.ImArrayParam resolvedParam =
                    customImValues == null
                            ? (imParam != null ? imParam : NumericUtil.ImArrayParam.builder().build())
                            : null;

            ImValueTransformer resolvedTransformer =
                    customImValues != null
                            ? new ImValueTransformer.CustomImTransformation(customImValues)
                            : (transformer != null
                                ? transformer
                                : new ImValueTransformer.NoImTransformation());

            return new ImConfiguration(
                    resolvedParam,
                    customImValues,
                    resolvedTransformer);
        }

        // ---------------------------------------------------------
        // BUILDER
        // ---------------------------------------------------------

        public static Builder builder() {
            return new Builder();
        }

        public static final class Builder {

            private Set<String> fragilityNames;
            private LimitState selectedLS;
            private NumericUtil.ImArrayParam imParam;
            private double[] customImValues;
            private ImValueTransformer transformer;

            /** Restrict preparation to selected model names. */
            public Builder filterNames(Set<String> names) {
                this.fragilityNames = names;
                return this;
            }

            /** Select limit state using enum. */
            public Builder selectedLS(LimitState ls) {
                this.selectedLS = ls;
                return this;
            }
            
            /**
             * Select limit state using string (case-insensitive).
             *
             * @param lsName limit state name
             * @throws IllegalArgumentException if invalid
             */
            public Builder selectedLS(String lsName) {

                if (lsName == null) {
                    this.selectedLS = null;
                    return this;
                }

                try {
                    this.selectedLS =
                            LimitState.valueOf(lsName.trim().toUpperCase());
                } catch (IllegalArgumentException e) {
                    throw new IllegalArgumentException(
                            "Invalid LimitState: " + lsName, e);
                }

                return this;
            }

            /** Use generated IM grid. */
            public Builder imParam(NumericUtil.ImArrayParam param) {
                this.imParam = param;
                this.customImValues = null;
                return this;
            }

            /**
             * Use explicit IM values.
             *
             * <p>Automatically propagates to discrete fragilities
             * via {@link ImValueTransformer.CustomImTransformation}.
             */
            public Builder customImValues(double[] values) {
                this.customImValues = values.clone();
                this.imParam = null;
                return this;
            }

            /** Custom transformer override. */
            public Builder transformer(ImValueTransformer transformer) {
                this.transformer = transformer;
                return this;
            }

            public PrepareSpec build() {

                if (imParam != null && customImValues != null) {
                    throw new IllegalStateException(
                            "Specify either imParam OR customImValues");
                }

                return new PrepareSpec(this);
            }
        }
    }
}