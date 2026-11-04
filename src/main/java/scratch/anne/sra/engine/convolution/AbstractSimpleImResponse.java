package scratch.anne.sra.engine.convolution;

import scratch.anne.sra.util.NumericUtil;

/**
 * Abstract base class for simple IM-based response functions.
 * Handles storing edges, computing midpoints, and IM values.
 * Concrete classes (vulnerability, fragility) extend this.
 */
public abstract class AbstractSimpleImResponse implements SimpleImResponseInterface {

    protected final double[] imValues;
    protected final double[] logimValues;
    protected final double[] respEdges;
    protected final double[] respMid;

    /**
     * Constructs a response from IM values and response edges.
     *
     * @param imValues  Array of intensity measure values
     * @param respEdges Array of response values (DR or probability)
     */
    protected AbstractSimpleImResponse(double[] imValues, double[] respEdges) {
        if (imValues == null || respEdges == null || imValues.length != respEdges.length) {
            throw new IllegalArgumentException("IM and response edges must be non-null and the same length");
        }
        this.imValues = imValues.clone();
        this.logimValues = NumericUtil.computeLogValues(imValues);
        this.respEdges = respEdges.clone();
        this.respMid = computeMidpoints(respEdges);
    }

    @Override
    public double[] getRespEdges() { return respEdges.clone(); }

    @Override
    public double[] getRespMid() { return respMid.clone(); }

    @Override
    public double[] getImValues() { return imValues.clone(); }
    
    @Override
    public double[] getLogImValues() { return logimValues.clone(); }

    /**
     * Compute midpoints of a response array.
     * <p>
     * mid[i] = 0.5 * (edges[i] + edges[i+1])
     * Last value is copied to preserve array length.
     *
     * @param edges input array (DRs, probabilities, etc.)
     * @return new array of midpoints (same length as input)
     */
    private static double[] computeMidpoints(double[] edges) {
        int n = edges.length;
        double[] mid = new double[n];
        for (int i = 0; i < n - 1; i++) {
            mid[i] = 0.5 * (edges[i] + edges[i + 1]);
        }
        mid[n - 1] = edges[n - 1];
        return mid;
    }
}
