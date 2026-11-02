package scratch.anne.sra.util;

import java.util.Arrays;
import java.util.Objects;

import scratch.anne.sra.util.StringUtil.ImtPeriod;

/**
 * Utility keys for grouping assets by location and intensity measure.
 * <p>
 * Both keys are immutable and can be used as map keys or set elements.
 */
public final class AssetKeys {

    private AssetKeys() {
        // Prevent instantiation
    }

    // -----------------------
    // SiteKey
    // -----------------------

    /** Immutable site key for grouping by location (latitude, longitude, vs30) */
    public static final class SiteKey {
        private final double latitude;
        private final double longitude;
        private final double vs30;

        public SiteKey(double latitude, double longitude, double vs30) {
            this.latitude = latitude;
            this.longitude = longitude;
            this.vs30 = vs30;
        }

        public double getLat() { return latitude; }
        public double getLon() { return longitude; }
        public double getVs30() { return vs30; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof SiteKey)) return false;
            SiteKey siteKey = (SiteKey) o;
            return Double.compare(siteKey.latitude, latitude) == 0 &&
                   Double.compare(siteKey.longitude, longitude) == 0 &&
                   Double.compare(siteKey.vs30, vs30) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(latitude, longitude, vs30);
        }

        @Override
        public String toString() {
            return String.format("SiteKey[lat=%.5f, lon=%.5f, vs30=%.1f]",
                    latitude, longitude, vs30);
        }
    }

    // -----------------------
    // ImKey
    // -----------------------

    /** Immutable key for grouping by intensity measure type (IMT) and log-transformed IM values */
    public static final class ImKey {
	    private final String imtString;
	    private final ImDomain domain;
	    private final double[] values;
	    private final int hash;
	
	    public ImKey(String imtString, double[] values, ImDomain domain) {
	        this.imtString = Objects.requireNonNull(imtString);
	        this.domain = Objects.requireNonNull(domain);
	        this.values = values.clone();
	
	        this.hash = Objects.hash(
	                imtString,
	                domain,
	                Arrays.hashCode(this.values)
	        );
	    }
	
	    public String getImtString() { return imtString; }
	    public ImDomain getDomain() { return domain; }
	    public double[] getValues() { return values.clone(); }
	    
	    public ImtPeriod getImtPeriod() { return StringUtil.stringToImtPeriod(imtString); }
	    
	    /**
	     * Returns the IM values as log(IM), converting from linear if needed.
	     */
	    public double[] getLogValues() {
	    	if (domain == ImDomain.LOG_IM ) {
	    		return values;
	    	}
	    	else if (domain == ImDomain.LINEAR_IM) {
	            return NumericUtil.computeLogValues(values);
	        } else {
	        	throw new IllegalArgumentException("Unknown ImDomain: " + domain);
	        }
	    }
	    
	    /**
	     * Returns the IM values as linear(IM), converting from log if needed.
	     */
	    public double[] getLinearValues() {
	    	if (domain == ImDomain.LINEAR_IM ) {
	    		return values;
	    	}
	    	else if (domain == ImDomain.LOG_IM) {
	            return NumericUtil.computeExpValues(values);
	        } else {
	        	throw new IllegalArgumentException("Unknown ImDomain: " + domain);
	        }
	    }	    
	
	    @Override
	    public boolean equals(Object o) {
	        if (this == o) return true;
	        if (!(o instanceof ImKey)) return false;

	        ImKey other = (ImKey) o;

	        return Objects.equals(imtString, other.imtString)
	                && Arrays.equals(values, other.values);
	    }

	    @Override
	    public int hashCode() {
	        return hash;
	    }
	
	    @Override
	    public String toString() {
	        return imtString + " [" + domain + ", n=" + values.length + "]";
	    }
	    
	    public static enum ImDomain {
	        LINEAR_IM,
	        LOG_IM,
	        CUSTOM
	    }
	}

}