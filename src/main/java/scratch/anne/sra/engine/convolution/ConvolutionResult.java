package scratch.anne.sra.engine.convolution;

/**
 * Result of one {@link RiskConvolution#compute()} call: the overall
 * risk value plus its per-IM-bin disaggregated contribution, from which
 * a cumulative disaggregation can be derived on demand via
 * {@link #getDisaggCumulative()}.
 */
public class ConvolutionResult {

	private final double[] imls;
    private final double risk;
    private final double[] disaggContribution;

    private double[] disaggCumulative;

    public ConvolutionResult(double[] imls, double[] disaggContribution, double risk) {
    	this.imls = imls;
        this.risk = risk;
        this.disaggContribution = disaggContribution.clone();
        this.disaggCumulative = null;  // only computed if called
    }
    
    private static double[] buildCumulative(double[] contribution) {

        double[] cumulative = new double[contribution.length];

        double running = 0.0;

        for (double v : contribution) {
            running += v;
        }

        if (running == 0.0) {
            // avoid divide-by-zero; return flat zeros
            return cumulative;
        }

        double norm = 1.0 / running;

        double delta = 0.0;

        for (int i = 0; i < contribution.length; i++) {
            delta += contribution[i];
            cumulative[i] = delta * norm;
        }

        // enforce exact closure (numerical safety)
        cumulative[cumulative.length - 1] = 1.0;

        return cumulative;
    }
    
    public double[] getImls() {
    	return imls;
    }
    
    public double getRisk() {
    	return risk;
    }
    
    public double[] getDisaggContribution() {
        return disaggContribution;
    }

    public double[] getDisaggCumulative() {
        if (disaggCumulative == null) {
            disaggCumulative = buildCumulative(disaggContribution);
        }
        return disaggCumulative;
    }
}
