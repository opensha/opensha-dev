package scratch.anne.sra.util;

import java.util.Arrays;

public class NumericUtil {
	
    /**
     * Compute log-transformed values of a response array.
     * Useful for building hazard curves or IM keys.
     *
     * @param values input array
     * @return new array of log(values)
     */
    public static double[] computeLogValues(double[] values) {
        double[] logValues = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            logValues[i] = Math.log(values[i]);
        }
        return logValues;
    }

    /**
     * Compute exp-transformed values, the inverse of {@link #computeLogValues}.
     *
     * @param values input array of log values
     * @return new array of exp(values)
     */
    public static double[] computeExpValues(double[] values) {
        double[] expValues = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            expValues[i] = Math.exp(values[i]);
        }
        return expValues;
    }
    
    
    /** Standard normal CDF (Φ) 
     * 
     * Using Abramowitz & Stegun approximation
     * */
    public static double standardNormalCDF(double x) {
        double t = 1.0 / (1.0 + 0.2316419 * Math.abs(x));
        double d = 0.3989423 * Math.exp(-x * x / 2.0);
        double prob = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
        return (x >= 0.0) ? 1.0 - prob : prob;
    }    

    
    /**
     * Generate log-spaced IM grid using fixed log10 step.
     *
     * @param min positive lower bound
     * @param max positive upper bound
     * @param logStep spacing in log10 space
     */
    public static double[] distributeLogSpacedValues(double min, double max, double logStep) {

        if (min <= 0 || max <= 0)
            throw new IllegalArgumentException("IM bounds must be positive");

        if (max <= min)
            throw new IllegalArgumentException("max must be > min");

        if (logStep <= 0)
            throw new IllegalArgumentException("logStep must be positive");

        int nPts = (int) Math.ceil(
            Math.log10(max / min) / logStep
        );

        double[] logArray = new double[nPts + 1];

        for (int i = 0; i <= nPts; i++) {
        	logArray[i] = min * Math.pow(10.0, i * logStep);
        }

        logArray[nPts] = max; // enforce exact bound

        return logArray;
    }
    
    
    /**
     * Generate log-spaced IM grid based on fragility median.
     *
     * @param median median of fragility
     * @param logStep spacing in log10 space
     * @param spread multiplier for min/max points
     */
    public static double[] distributeImAroundMedian(double median, ImArrayParam param) {
    	
    	double min = median / param.spread;
    	double max = median * param.spread;
    	
    	return distributeLogSpacedValues(min, max, param.logStep);
    }
    
    /** builder for the im array params in distributeImAroundMedian */
    public static class ImArrayParam {

        // defaults
        public static final double DEFAULT_LOG_STEP = 0.025;
        public static final double DEFAULT_SPREAD = 10.0;

        public final double logStep;
        public final double spread;

        private ImArrayParam(Builder b) {
            this.logStep = b.logStep;
            this.spread = b.spread;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static class Builder {

            private double logStep = DEFAULT_LOG_STEP;
            private double spread = DEFAULT_SPREAD;

            public Builder logStep(double logStep) {
                if (logStep <= 0)
                    throw new IllegalArgumentException("logStep must be positive");
                this.logStep = logStep;
                return this;
            }

            public Builder spread(double spread) {
                if (spread <= 1)
                    throw new IllegalArgumentException("spread must be > 1");
                this.spread = spread;
                return this;
            }

            public ImArrayParam build() {
                return new ImArrayParam(this);
            }
        }
    }
    
    /**
     * Throws if {@code imLevels} is not strictly increasing.
     *
     * @param name label used in the error message (e.g. the model/array's name)
     * @param imLevels array to check
     */
    public static void checkMonotonicIM(String name, double[] imLevels) {
        double prev = Double.NEGATIVE_INFINITY;
        for (double im : imLevels) {
            if (im <= prev) {
                throw new IllegalArgumentException(
                        "IM levels for vulnerability '" + name + "' are not strictly increasing: " +
                                Arrays.toString(imLevels));
            }
            prev = im;
        }
    }

    public static double linearInterpolate(
            double[] x,
            double[] y,
            double xi) {

        int n = x.length;

        int idx = Arrays.binarySearch(x, xi);

        // exact match
        if (idx >= 0) {
            return y[idx];
        }

        int insert = -idx - 1;

        // ---------------- below range ----------------
        if (insert == 0) {
            return y[0];
        }

        // ---------------- above range ----------------
        if (insert >= n) {
            return y[n - 1];
        }

        // ---------------- linear interpolation ----------------
        double x0 = x[insert - 1];
        double x1 = x[insert];
        double y0 = y[insert - 1];
        double y1 = y[insert];

        return y0 + (y1 - y0) * (xi - x0) / (x1 - x0);
    }
}
