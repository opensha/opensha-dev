package scratch.anne.sra.util;

import java.util.Arrays;
import java.util.TreeSet;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Provides the standard set of exceedance-rate levels used for LEC
 * sampling: a dense head (0.1-9, covering very frequent events), the
 * CEA-style standard return periods (5-10,000 years), and a sparse tail
 * of very rare rates. Deduplicated and sorted descending.
 *
 * <p>Currently fixed; no inputs. If/when this needs to be configurable
 * (different head/tail ranges, different CEA periods), add parameters
 * to {@link #get()} or convert it to an instance with a builder --
 * the three-piece construction below is written so each piece can be
 * swapped out independently later without touching the others.</p>
 */
public final class StandardRateLikelihoods {

    private StandardRateLikelihoods() {}

    /**
     * @return exceedance-rate levels, deduplicated, sorted descending
     *         (most frequent/highest rate first)
     */
    public static double[] get() {

//        // dense head: 0.1, 0.2, ..., 0.9, 1, 2, ..., 9 -- i.e.
//        // np.arange(1,10,1) * 10^mag for mag in [-1, 0]
//        double[] denseHead = IntStream.rangeClosed(-1, 0)
//                .boxed()
//                .flatMap(mag -> IntStream.range(1, 10)
//                        .mapToObj(i -> i * Math.pow(10, mag)))
//                .mapToDouble(Double::doubleValue)
//                .toArray();
        
     // dense head: evenly spaced in log10-space from 0.1 to 9, so points
        // appear linearly spaced on a log-scale plot (constant ratio between
        // consecutive points, rather than constant difference)
        double[] denseHead = logspace(-1, Math.log10(9), 20);

        // CEA-style standard return periods (years)
        double[] ceaReturnPeriods = {5, 10, 25, 50, 100, 200, 250, 345, 350, 400, 500, 1000, 2000, 5000, 10000};

        // sparse tail: rates directly, no return-period round-trip needed
        double[] sparseTailRates = {3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7};

        double[] denseHeadRates = Arrays.stream(denseHead).map(t -> 1d / t).toArray();
        double[] ceaRates = Arrays.stream(ceaReturnPeriods).map(t -> 1d / t).toArray();

        TreeSet<Double> unique = DoubleStream.concat(
                DoubleStream.concat(Arrays.stream(denseHeadRates), Arrays.stream(ceaRates)),
                Arrays.stream(sparseTailRates)
        ).boxed().collect(java.util.stream.Collectors.toCollection(TreeSet::new));

        double[] sortedAscending = unique.stream().mapToDouble(Double::doubleValue).toArray();

        // reverse for descending order
        double[] rates = new double[sortedAscending.length];
        for (int i = 0; i < sortedAscending.length; i++)
            rates[i] = sortedAscending[sortedAscending.length - 1 - i];

        return rates;
    }
    
    /**
     * n points evenly spaced in log10-space between 10^startExp and 10^endExp
     * (inclusive), analogous to numpy's logspace.
     */
    private static double[] logspace(double startExp, double endExp, int n) {
        double[] result = new double[n];
        double step = (endExp - startExp) / (n - 1);
        for (int i = 0; i < n; i++)
            result[i] = Math.pow(10, startExp + i * step);
        return result;
    }
}
