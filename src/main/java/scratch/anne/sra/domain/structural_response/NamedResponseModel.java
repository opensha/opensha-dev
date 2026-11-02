package scratch.anne.sra.domain.structural_response;

import scratch.anne.sra.util.enums.IMT;

/**
 * Interface for response models (vulnerabilities, fragilities, etc.)
 * that have a unique name.
 */
public interface NamedResponseModel {

    /** Unique name */
    String getName();

    /** IMT object for this response */
    IMT getImt();

    /** Get spectral period (null if IMT != SA) */
    Double getPeriod();
    
    /** String representation of IMT (e.g., "SA(1.0)") */
    String getImtString();

}
