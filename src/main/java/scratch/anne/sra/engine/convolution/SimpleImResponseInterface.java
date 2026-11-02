package scratch.anne.sra.engine.convolution;

import scratch.anne.sra.domain.structural_response.NamedResponseModel;

/**
 * Generic response function over intensity measures (IM), providing:
 * <ul>
 *   <li>Fragility: probability of exceeding a damage state</li>
 *   <li>Vulnerability: expected loss ratio / consequence</li>
 * </ul>
 *
 * <p>This interface extends {@link NamedResponseModel}, so any concrete
 * implementation can optionally have a name, but must expose IM values
 * and response arrays.
 */
public interface SimpleImResponseInterface extends NamedResponseModel {

    /** Intensity measure levels (IMs) that define the response bins */
    double[] getImValues();
    
    /** Log of IMs that define the response bins */
    double[] getLogImValues();

    /** Response edges: per-IM values (fragility: P[D>ds], vulnerability: expected loss ratio) */
    double[] getRespEdges();

    /** Midpoints between respEdges for integration */
    double[] getRespMid();
}