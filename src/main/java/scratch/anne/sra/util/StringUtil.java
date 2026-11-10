package scratch.anne.sra.util;

import java.text.DecimalFormat;

import scratch.anne.sra.util.enums.IMT;

public class StringUtil {
	
    /**
     * Return string of the im type
     *
     * @param imtEnum imt type (PGA, PGV, SA)
     * @param period spectral period, if SA
     * @return im type as a string with spectral period (if applicable)
     */
	
	public static String imtToString(IMT imtEnum, Double period) {
	    if (imtEnum == null) {
	        throw new IllegalArgumentException("IMT cannot be null");
	    }
	    if (imtEnum.name().equalsIgnoreCase("SA") && period != null) {
	        return String.format("%s(%.2g)", imtEnum.name(), period);
	    } else {
	        return imtEnum.name(); // period is ignored
	    }
	}

	/** class to return both enum and period together */
    public static class ImtPeriod {
        public final IMT imt;
        public final Double period;

        public ImtPeriod(IMT imt, Double period) {
            this.imt = imt;
            this.period = period;
        }
    }

    /**
     * Converts a string like "PGA" or "SA(1.0)" into an IMT enum + optional period.
     * Prints any errors immediately and propagates them.
     */
    public static ImtPeriod stringToImtPeriod(String imtStr) {

        if (imtStr == null) {
            throw new IllegalArgumentException("IMT string is null");
        }

        imtStr = imtStr.trim();

        try {
            if (imtStr.contains("(") && imtStr.endsWith(")")) {

                int parenIndex = imtStr.indexOf('(');
                String enumPart = imtStr.substring(0, parenIndex);
                String periodPart = imtStr.substring(parenIndex + 1, imtStr.length() - 1);

                IMT imtEnum = IMT.valueOf(enumPart.toUpperCase());
                double period = Double.parseDouble(periodPart);

                return new ImtPeriod(imtEnum, period);

            } else {
                IMT imtEnum = IMT.valueOf(imtStr.toUpperCase());
                return new ImtPeriod(imtEnum, null);
            }

        } catch (Exception e) {
            throw new IllegalArgumentException(
                "Invalid IMT string: \"" + imtStr + "\"",
                e   // preserve cause
            );
        }
    }
    
    
    /** Canonical string form of a likelihood level, avoiding scientific-notation/formatting drift across runs. */
    public static String formatLikelihood(double level) {
        return new DecimalFormat("0.0###E0").format(level);
    }
    
    /** Canonical string form of a percentile, avoiding (int) truncation dropping values like 2.5 or 97.5. */
    public static String formatPercentile(double p) {
        return new DecimalFormat("0.#").format(p);
    }
}
