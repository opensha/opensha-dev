package scratch.anne.sra.engine.convolution;

import org.opensha.commons.data.function.DiscretizedFunc;

import scratch.anne.sra.util.HazardFunctionUtil;

/**
 * Simple risk convolution integral, for a given hazard and a simple response given hazard function (e.g., a vulnerability or fragility).
 * <p>
 * The output value is: <br>
 * <ul>
 * <li> for vulnerability: Expected loss, normalized by the asset value (0 to 1)
 * <li> for fragility: Risk (rate or probability) of damage state exceedance (e.g., collapse or 
 * </ul>
 * Timespan depends on the input hazard, typically annual. <br>
 * Accepts (and outputs) either rate or probability of exceedance. NOTE: risk convolution is formally defined in rate-space. <br>
 * Uses Riemann or Closed-Form integration.
 * <blockquote>
 * 		Close-Form solution: <br>
 * 		Porter, K.A.,C.R. Scawthorn, and J.L. Beck, 2006. Cost-effectiveness of stronger
 *			woodframe buildings. <br>
 *          Earthquake Spectra 22 (1), February 2006, 239-266, Equation #2 <br>
 *		    http://www.sparisk.com/pubs/Porter-2006-woodframe.pdf
 * </blockquote>
 */
public class RiskConvolution {

    public enum IntegrationMethod { RIEMANN, CLOSED_FORM }

    private final SimpleImResponseInterface response;  // prepared response (includes imEdges, respEdges, respMid, etc.)
    private double[] hazardValues;           // either rate or probability of exceedance, at response's imEdges
                                            // can also be null, if hazard will be supplied as .compute(hazard)
    private final IntegrationMethod method;

 // ===================== Constructors =====================

    /**
     * Constructor using a pre-built SimpleImResponse and hazard array.
     * Integration method defaults to RIEMANN if not specified.
     *
     * @param response    SimpleImResponse (including im bin edges, response at im edges, mid-bin response
     * @param hazardValue Hazard value (rate or probability of exceedance at vuln's imEdges)
     * @param method      Integration method (default: RIEMANN)
     */
    public RiskConvolution(SimpleImResponseInterface response, double[] hazard, IntegrationMethod method) {
        if (hazard != null && response.getImValues().length != hazard.length) {
            throw new IllegalArgumentException("response IM and hazard array lengths must match");
        }
        this.response = response;
        this.hazardValues = (hazard != null) ? hazard : null;
        this.method = (method != null) ? method : IntegrationMethod.RIEMANN;
    }

    /**
     * Convenience constructor using a pre-built SimpleImResponse and hazard array. <br>
     * Uses default RIEMANN integration method.
     *
     * @param response   SimpleImResponse
     * @param hazard     Hazard array
     */
    public RiskConvolution(SimpleImResponseInterface response, double[] hazard) {
        this(response, hazard, null);
    }

    /**
     * Constructor using a pre-built SimpleImResponse and hazard DiscretizedFunc.
     * Checks that the hazard x-values match the vulnerability's IM edges.
     * Default integration method is RIEMANN if not specified.
     *
     * @param response   SimpleImResponse
     * @param hazardFunc Hazard as DiscretizedFunc (x = IM, y = hazard value)
     * @param method     Integration method; if null, defaults to RIEMANN
     */
    public RiskConvolution(SimpleImResponseInterface response, DiscretizedFunc hazardFunc, IntegrationMethod method) {
        this.response = response;
        this.hazardValues = HazardFunctionUtil.convertHazFuncToArray(hazardFunc, response.getImValues());
        this.method = (method != null) ? method : IntegrationMethod.RIEMANN;
    }

    /**
     * Convenience constructor using a pre-built SimpleImResponse  and hazard DiscretizedFunc.
     * Defaults to RIEMANN integration method.
     *
     * @param response   SimpleImResponse 
     * @param hazardFunc Hazard as DiscretizedFunc
     */
    public RiskConvolution(SimpleImResponseInterface response, DiscretizedFunc hazardFunc) {
        this(response, hazardFunc, null);
    }

    /**
     * Constructor that builds an internal SimpleImResponse  from raw IM and response arrays,
     * along with a hazard array. Integration method defaults to RIEMANN if not specified.
     *
     * @param imls       Intensity measure array
     * @param fim      Response-as-a-function-of-im array
     * @param hazard   Hazard values corresponding to IMs
     * @param method   Integration method; if null, defaults to RIEMANN
     */
    public RiskConvolution(double[] imls, double[] fim, double[] hazard, IntegrationMethod method) {
        if (imls.length != fim.length || imls.length != hazard.length) {
            throw new IllegalArgumentException("IM, F(IM), and hazard arrays must all be the same length");
        }
        this.response = new SimpleImResponseFunction(imls, fim);
        this.hazardValues = hazard;
        this.method = (method != null) ? method : IntegrationMethod.RIEMANN;
    }

    /**
     * Convenience constructor that builds a vulnerability from IM/response arrays
     * and uses default RIEMANN integration method.
     *
     * @param imls     Intensity measure array
     * @param fim      Response-as-a-function-of-im array
     * @param hazard Hazard array
     */
    public RiskConvolution(double[] imls, double[] fim, double[] hazardValues) {
        this(imls, fim, hazardValues, null);
    }

    /**
     * Constructor using only the SimpleImResponse.
     * Hazard must be supplied later via compute().
     * Integration method defaults to RIEMANN if not specified.
     *
     * @param response   SimpleImResponse
     * @param method     Integration method; if null, defaults to RIEMANN
     */
    public RiskConvolution(SimpleImResponseInterface response, IntegrationMethod method) {
        this.response = response;
        this.hazardValues = null;
        this.method = (method != null) ? method : IntegrationMethod.RIEMANN;
    }
    
    /**
     * Constructor using only the SimpleImResponse.
     * Hazard must be supplied later via compute().
     * Integration method defaults to RIEMANN.
     *
     * @param response SimpleImResponse
     */
    public RiskConvolution(SimpleImResponseInterface response) {
        this.response = response;
        this.hazardValues = null;
        this.method = null;
    }
    
    /**
     * Constructor using only IM and F(IM) arrays.
     * Builds SimpleImResponse internally; hazard must be supplied later via compute().
     * Integration method defaults to RIEMANN if not specified.
     *
     * @param imls      Intensity measure array
     * @param fim     Response-as-a-function-of-im array
     * @param method  Integration method; if null, defaults to RIEMANN
     */
    public RiskConvolution(double[] imls, double[] fim, IntegrationMethod method) {
        this(new SimpleImResponseFunction(imls, fim), method);
    }

    /**
     * Constructor using only IM and F(IM) arrays.
     * Builds SimpleImResponse internally; hazard must be supplied later via compute().
     * Integration method defaults to RIEMANN.
     *
     * @param imls      Intensity measure array
     * @param fim     Response-as-a-function-of-im array
     */
    public RiskConvolution(double[] imls, double[] fim) {
        this.response = new SimpleImResponseFunction(imls, fim);
        this.hazardValues = null;
        this.method = null;
    }

 // -------------------- Compute Overloads --------------------

    /** Use stored hazard array */
    public ConvolutionResult compute() {
        if (hazardValues == null) {
            throw new IllegalStateException("Hazard array is not set. Provide hazard to compute.");
        }
        return computeInternal(hazardValues);
    }

    /** Provide hazard as double array */
    public ConvolutionResult compute(double[] hazard) {
        if (hazard.length != response.getImValues().length) {
            throw new IllegalArgumentException("Hazard array length does not match vulnerability IM length");
        }
        this.hazardValues = hazard;
        return computeInternal(hazardValues);
    }

    /** Provide hazard as DiscretizedFunc (checks IMs) */
    public ConvolutionResult compute(DiscretizedFunc hazardFunc) {
        this.hazardValues = HazardFunctionUtil.convertHazFuncToArray(hazardFunc, response.getImValues());
        return computeInternal(hazardValues);
    }
    
 // ---------------------- Internal shared compute ----------------------

    /**
     * Shared internal computation, chooses method based on this.method
     */
    private ConvolutionResult computeInternal(double[] hazardValues) {
        if (method == IntegrationMethod.RIEMANN) {
            return computeRiemann(hazardValues);
        } else {
            return computeClosedForm(hazardValues);
        }
    }
    
    // ---------------------- Core Riemann Integration ----------------------
    private ConvolutionResult computeRiemann(double[] hazardValues) {
        
    	final double[] hazard = hazardValues;
        final double[] respMid = response.getRespMid();
        final int n = respMid.length;
        
        double risk = 0.0;
        double[] contribution = new double[n];

        
        // interior bins
        double hazLeft = hazard[0];
        for (int i = 0; i < respMid.length - 1; i++) {
        	double hazRight = hazard[i + 1];
            double deltaHazard = hazLeft - hazRight;
        	hazLeft = hazRight;
            contribution[i] = respMid[i] * deltaHazard;
            risk += contribution[i];
        }

        // tail bin
        contribution[n - 1] = respMid[n - 1] * hazard[n - 1];
        risk += contribution[n - 1];

        return new ConvolutionResult(response.getImValues(), contribution, risk);
    }


    // ---------------------- Porter Closed-Form Integration ----------------------
    /**
     * Performance and numerical optimization note:
     *
     * This implementation uses the hazard ratio directly rather than explicitly computing
     * the intermediate gradient term G = log(h[i]/h[i-1]) / ΔIML.
     *
     * This is mathematically equivalent in the closed-form expression, since:
     * exp(G · ΔIML) = h[i] / h[i-1]
     *
     * Using the ratio form:
     *  - eliminates one logarithm and one exponential evaluation per iteration
     *  - reduces numerical overhead in the hot loop
     *  - improves instruction throughput in the JIT-compiled kernel
     *
     * This change is purely computationally equivalent and does not alter the
     * mathematical result of the convolution.
     */
    private ConvolutionResult computeClosedForm(double[] hazardValues) {
    	
    	final double[] hazard = hazardValues;
        final double[] imls = response.getImValues();
        final double[] respEdge = response.getRespEdges();

        final int n = imls.length;

        double risk = 0.0;
        double[] contribution = new double[n];

        double imlLeft = imls[0];
        double hazLeft = hazard[0];
        double respLeft = respEdge[0];
        for (int i = 0; i < n - 1; i++) {

            double imlRight = imls[i + 1];
            double hazRight = hazard[i + 1];
            double respRight = respEdge[i + 1];

            double imlDelta = imlRight - imlLeft;
            double invIml = 1.0 / imlDelta;

            double respDelta = (respRight - respLeft) * invIml;

            double ratio = hazRight / hazLeft;
            double delta;
            
            if (hazLeft == 0.0) {
            	// hazard is 0 so risk is 0
                delta = 0.0;
            } else if (Math.abs(ratio - 1.0) < 1e-10) {
            	// hazard is flat so the integral is a rectangular bin
                delta = hazLeft * (respRight - respLeft);
            } else {
            	// log-linear assumption
                double invG = imlDelta / Math.log(ratio);
                double term1 = respLeft * hazLeft * (1.0 - ratio);
                double term2 = respDelta * hazLeft * (ratio * (imlDelta - invG) + invG);
                delta = term1 - term2;
            }

            contribution[i + 1] = delta;
            risk += delta;

            imlLeft = imlRight;
            hazLeft = hazRight;
            respLeft = respRight;
        }

        return new ConvolutionResult(imls, contribution, risk);
    }
    

}