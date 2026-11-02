package scratch.anne.sra.util;

import org.opensha.commons.data.function.DiscretizedFunc;

public class HazardFunctionUtil {
	
    /**
     * Converts hazard function to array, ensuring IMs match.
     * Optionally allows hazardFunc X-values to be log(IMs).
     *
     * @param hazardFunc hazard curve (x = IM or log(IM), y = hazard)
     * @param vulnIMs    vulnerability linear IM values
     */
    public static double[] convertHazFuncToArray(
            DiscretizedFunc hazardFunc,
            double[] vulnIMs) {

        if (hazardFunc.size() != vulnIMs.length) {
            throw new IllegalArgumentException(
                "Hazard function size does not match vulnerability IM length");
        }

        int n = vulnIMs.length;
        double[] hazardValues = new double[n];

        // First pass: determine if X matches linear or log(IM)
        boolean matchesLinear = true;
        boolean matchesLog = true;

        for (int i = 0; i < n; i++) {
            double x = hazardFunc.getX(i);

            if (x != vulnIMs[i]) {
                matchesLinear = false;
            }

            if (x != Math.log(vulnIMs[i])) {
                matchesLog = false;
            }

            if (!matchesLinear && !matchesLog) {
                throw new IllegalArgumentException(
                    String.format(
                        "IM values do not match at index %d: vuln IM=%f, hazard IM=%f",
                        i, vulnIMs[i], x
                    )
                );
            }
        }

        // Second pass: retrieve Y-values (no transformation needed)
        for (int i = 0; i < n; i++) {
            hazardValues[i] = hazardFunc.getY(i);
        }

        return hazardValues;
    }

}
