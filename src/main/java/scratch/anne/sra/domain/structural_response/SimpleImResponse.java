package scratch.anne.sra.domain.structural_response;

import scratch.anne.sra.engine.convolution.AbstractSimpleImResponse;
import scratch.anne.sra.util.StringUtil;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.enums.IMT;

/**
 * Represents a discrete response function prepared for convolution.
 *
 * <p>This class extends {@link AbstractSimpleImResponse}, which provides the core
 * response data arrays used for risk convolution:
 * <ul>
 *   <li>{@code imValues} – intensity measure values at bin edges</li>
 *   <li>{@code respEdges} – response (damage ratio) values at bin edges</li>
 *   <li>{@code respMid} – midpoint response values, automatically computed by the abstract class</li>
 * </ul>
 *
 * <p>Additional metadata stored in this subclass:</p>
 * <ul>
 *   <li>{@code name} – optional vulnerability name</li>
 *   <li>{@code imt} – intensity measure type (optional)</li>
 *   <li>{@code period} – spectral period, if applicable (optional)</li>
 *   <li>{@code imtString} – string representation of IMT + period</li>
 * </ul>
 *
 * <p>Use getters that safely clone the intensity measure and response values.</p>
 *
 * @see AbstractSimpleImResponse
 */
public class SimpleImResponse extends AbstractSimpleImResponse implements NamedResponseModel {

    // ------------------------------------------------------------------------
    // Metadata
    // ------------------------------------------------------------------------
    private final String name;
    private final IMT imt;
    private final Double period;
    private final String imtString;
    private final ImKey imKey;

    // ------------------------------------------------------------------------
    // Constructor: full control, specify IM domain
    // ------------------------------------------------------------------------
    public SimpleImResponse(String name, IMT imt, Double period,
                                 double[] imValues, double[] responseEdges,
                                 ImKey.ImDomain domain) {
        super(imValues, responseEdges);
        this.name = name;
        this.imt = imt;
        this.period = period;
        this.imtString = StringUtil.imtToString(imt, period);
        
        double[] keyValues;
        switch (domain) {
            case LOG_IM:
                keyValues = this.logimValues;
                break;
            case LINEAR_IM:
            	keyValues = this.imValues;
                break;
            default:
                throw new IllegalArgumentException("Unknown ImDomain: " + domain);
        }

        // create the key
        this.imKey = new ImKey(this.imtString, keyValues, domain);

    }

    // ------------------------------------------------------------------------
    // Convenience constructor: assume log-transformed IM values
    // ------------------------------------------------------------------------
    public SimpleImResponse(String name, IMT imt, Double period,
                                 double[] imValues, double[] responseEdges) {
        this(name, imt, period, imValues, responseEdges, ImKey.ImDomain.LOG_IM);
    }

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------
    public String getName() { return name; }
    public IMT getImt() { return imt; }
    public Double getPeriod() { return period; }
    public String getImtString() { return imtString; }

    public ImKey getImKey() { return imKey; }

    public double[] getImValues() { return imValues.clone(); }
    public double[] getLogImValues() { return logimValues.clone(); }
    public double[] getRespEdges() { return respEdges.clone(); }
    public double[] getRespMid() { return respMid.clone(); }

    @Override
    public String toString() {
        return "SimpleImResponse{" +
               "name='" + name + '\'' +
               ", imt=" + imtString +
               ", nIM=" + imValues.length +
               '}';
    }
    
    public String responseTable() {

        StringBuilder sb = new StringBuilder();

        sb.append("Simple IM Response\n");
        sb.append("------------------\n");
        sb.append("Name : ").append(name).append("\n");
        sb.append("IMT  : ").append(imtString).append("\n");
        sb.append("N    : ").append(imValues.length).append("\n\n");

        sb.append(String.format("%15s | %15s%n", "IM Value", "Response"));
        sb.append("-----------------+-----------------\n");

        for (int i = 0; i < imValues.length; i++) {
            sb.append(String.format(
                    "%15.6e | %15.6e%n",
                    imValues[i],
                    respEdges[i]
            ));
        }

        return sb.toString();
    }
}