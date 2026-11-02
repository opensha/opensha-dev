package scratch.anne.sra.util;

import java.util.Arrays;
import java.util.List;

public class InputParsingUtil {
	
    public static double parseDoubleOrThrow(String s) {
        if (s == null || s.isEmpty()) {
            throw new IllegalArgumentException("Expected numeric value but found empty string");
        }
        return Double.parseDouble(s.trim());
    }

    
    public static double[] parseJsonNumericArray(Object obj) {
        if (obj == null) return new double[0];
        if (obj instanceof List<?>) {
            List<?> list = (List<?>) obj;
            return list.stream()
                       .mapToDouble(v -> {
                           if (v instanceof Number) return ((Number) v).doubleValue();
                           else return Double.parseDouble(v.toString());
                       })
                       .toArray();
        } else if (obj instanceof String) {
            String[] tokens = ((String) obj).split(";");
            return Arrays.stream(tokens)
                         .mapToDouble(Double::parseDouble)
                         .toArray();
        } else {
            throw new IllegalArgumentException("Cannot parse numeric array: " + obj);
        }
    }

    
    /**
     * Cleans a CSV field by trimming whitespace and removing surrounding quotes.
     *
     * @param s raw CSV field
     * @return cleaned value (never null)
     */

    public static String cleanCsvString(String s) {
        if (s == null) return "";
        s = s.trim();

        if (s.startsWith("\"") && s.endsWith("\"") && s.length() > 1) {
            s = s.substring(1, s.length() - 1);
        }

        return s;
    }


	public static double parseDoubleOrDefault(String s, double defaultVal) {
	    if (s == null || s.isBlank()) {
	        return defaultVal;
	    }
	    return Double.parseDouble(s.trim());
	}
}
