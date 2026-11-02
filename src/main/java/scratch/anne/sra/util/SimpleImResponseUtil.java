package scratch.anne.sra.util;

public class SimpleImResponseUtil {

    /**
     * Compute midpoints of a response array.
     * <p>
     * mid[i] = 0.5 * (edges[i] + edges[i+1])
     * Last value is copied to preserve array length.
     *
     * @param edges input array (DRs, probabilities, etc.)
     * @return new array of midpoints (same length as input)
     */
    public static double[] computeMidpoints(double[] edges) {
        int n = edges.length;
        double[] mid = new double[n];
        for (int i = 0; i < n - 1; i++) {
            mid[i] = 0.5 * (edges[i] + edges[i + 1]);
        }
        mid[n - 1] = edges[n - 1];
        return mid;
    }
}
