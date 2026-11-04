package scratch.anne.sra.engine.hazard_setup;

import java.util.Comparator;

import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.constraint.impl.DoubleDiscreteConstraint;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;

import scratch.anne.sra.util.StringUtil.ImtPeriod;
import scratch.anne.sra.util.enums.IMT;

/**
 * Sets a GMM's intensity measure, including its period if the IMT is SA.
 *
 * <p>Shared by {@code PortfolioRiskConvolutionCalculator},
 * {@code PortfolioPerRuptureRiskConvolutionCalculator}, and
 * {@code RuptureProbabilityCurveCollection.fromERF}. Callers only need
 * to say "set this GMM to this IMT" -- how an SA period that isn't one
 * of the GMM's allowed discrete values gets resolved is entirely
 * {@link PeriodResolutionStrategy}'s concern (see
 * {@link PeriodResolutionStrategy#DEFAULT}), not this class's or its
 * callers'.</p>
 */
public final class GmmImtSetter {

    private GmmImtSetter() {}

    /**
     * Sets {@code gmm}'s intensity measure to {@code imt}. For SA,
     * also sets the period, resolved via {@link PeriodResolutionStrategy#DEFAULT}
     * if it isn't exactly one of the GMM's allowed values.
     */
    public static void setImt(ScalarIMR gmm, ImtPeriod imt) {
        gmm.setIntensityMeasure(imt.imt.name());

        if (imt.imt == IMT.SA) {
            PeriodResolutionStrategy.DEFAULT.resolve(gmm, imt.period);
        }
    }

    /**
     * How to resolve an SA period that isn't exactly one of a GMM's
     * allowed discrete values.
     *
     * <p>{@link #DEFAULT} is the single line to change to switch the
     * system-wide behavior -- every {@link GmmImtSetter#setImt} call
     * uses whatever strategy is configured here.</p>
     */
    public enum PeriodResolutionStrategy {

        /**
         * Snap to the nearest allowed period.
         *
         * TODO snap to nearest period is a workaround, not a real solution
         * find a better way to handle a requested period that isn't one of a
         * GMM's allowed discrete values (e.g. interpolating between the
         * two nearest periods' results, rather than silently
         * substituting a different period).
         */
        SNAP_TO_NEAREST {
            @Override
            void resolve(ScalarIMR gmm, double requestedPeriod) {
                @SuppressWarnings("unchecked")
                Parameter<Double> p = (Parameter<Double>) gmm.getParameter(PeriodParam.NAME);
                var constraint = (DoubleDiscreteConstraint) p.getConstraint();
                double allowed = constraint.getAllowedDoubles().stream()
                        .min(Comparator.comparing(d -> Math.abs(d - requestedPeriod)))
                        .orElseThrow();
                p.setValue(allowed);
            }
        };

        /** The strategy every {@link GmmImtSetter#setImt} call uses unless told otherwise. */
        public static final PeriodResolutionStrategy DEFAULT = SNAP_TO_NEAREST;

        abstract void resolve(ScalarIMR gmm, double requestedPeriod);
    }
}
