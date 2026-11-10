package scratch.anne.sra.domain.structuralResponse;

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

    /**
     * Throws unless {@code period} is present if and only if {@code imt} is SA.
     * Shared by every {@code NamedResponseModel} implementation so the
     * SA/period rule can't drift between them.
     */
    static void validateImtPeriod(IMT imt, Double period) {
        if (imt == IMT.SA && period == null)
            throw new IllegalStateException("SA requires a period");
        if (imt != IMT.SA && period != null)
            throw new IllegalStateException("Only SA can have a period");
    }

}
