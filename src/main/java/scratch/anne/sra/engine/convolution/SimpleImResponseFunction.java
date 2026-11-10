package scratch.anne.sra.engine.convolution;

import scratch.anne.sra.domain.structuralResponse.NamedResponseModel;
import scratch.anne.sra.util.enums.IMT;

/**
 * Generic response function over IMs when name or IMT is irrelevant.
 *
 * <p>This class implements {@link SimpleImResponseInterface} but returns null
 * for all {@link NamedResponseModel} getters.
 */
public class SimpleImResponseFunction extends AbstractSimpleImResponse {

    /**
     * Constructs a generic response.
     *
     * @param imValues  Array of intensity measure values
     * @param respEdges Array of response values (DR or probability)
     */
    public SimpleImResponseFunction(double[] imValues, double[] respEdges) {
        super(imValues, respEdges);
    }

    // ----------------------------
    // NamedResponseModel implementations (always null)
    // ----------------------------

    @Override
    public String getName() {
        return null;
    }

    @Override
    public IMT getImt() {
        return null;
    }
    
    @Override
    public Double getPeriod() { 
    	return null; 
    }

    @Override
    public String getImtString() {
        return null;
    }
}
