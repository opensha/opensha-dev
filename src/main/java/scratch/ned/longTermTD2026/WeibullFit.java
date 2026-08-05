package scratch.ned.longTermTD2026;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

/**
* Maximum-likelihood estimation of a two-parameter Weibull distribution.
*
* Method:
* - Shape parameter estimated by maximizing the Weibull log-likelihood
* with Brent's one-dimensional optimization algorithm.
* - Scale parameter computed from the closed-form MLE:
*
* lambda = ( (1/n) * sum(x_i^k) )^(1/k)
*
* References:
* Johnson, Kotz & Balakrishnan (1994),
* Continuous Univariate Distributions, Vol. 1.
*
* Optimization implemented using Apache Commons Math BrentOptimizer.
*/
public class WeibullFit {

    public static class WeibullParams {
        public final double shape; // k
        public final double scale; // lambda

        public WeibullParams(double shape, double scale) {
            this.shape = shape;
            this.scale = scale;
        }

        @Override
        public String toString() {
            return "shape k = " + shape + ", scale lambda = " + scale;
        }
    }

    public static WeibullParams fit(double[] samples) {
        if (samples == null || samples.length < 2)
            throw new IllegalArgumentException("Need at least two samples");

        for (double x : samples) {
            if (x <= 0.0 || !Double.isFinite(x))
                throw new IllegalArgumentException("Weibull samples must be positive and finite");
        }

        final int n = samples.length;

        // For any fixed shape k, the MLE scale is:
        // lambda = (mean(x_i^k))^(1/k)
        UnivariateFunction logLikelihoodForShape = new UnivariateFunction() {
            @Override
            public double value(double k) {
                if (k <= 0.0 || !Double.isFinite(k)) return Double.NEGATIVE_INFINITY;

                double sumLogX = 0.0;
                double sumXk = 0.0;

                for (double x : samples) {
                    sumLogX += Math.log(x);
                    sumXk += Math.pow(x, k);
                }

                double lambda = Math.pow(sumXk / n, 1.0 / k);
                double logLambda = Math.log(lambda);

                // Weibull log-likelihood:
                // sum[ log(k) - k log(lambda) + (k - 1) log(x_i) - (x_i/lambda)^k ]
                return n * Math.log(k)
                        - n * k * logLambda
                        + (k - 1.0) * sumLogX
                        - sumXk / Math.pow(lambda, k);
            }
        };

        BrentOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);

        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(1000),
                new UnivariateObjectiveFunction(logLikelihoodForShape),
                GoalType.MAXIMIZE,
                new SearchInterval(0.05, 20.0)
        );

        double k = result.getPoint();

        double sumXk = 0.0;
        for (double x : samples) {
            sumXk += Math.pow(x, k);
        }

        double lambda = Math.pow(sumXk / n, 1.0 / k);

        return new WeibullParams(k, lambda);
    }

    public static void main(String[] args) {
        double[] data = { 12.4, 9.8, 15.1, 11.7, 14.2, 10.9, 13.3 };

        WeibullParams p = fit(data);
        System.out.println(p);
    }
}