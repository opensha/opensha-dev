package scratch.anne.sra.domain.asset.fragility;

import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.engine.convolution.ConvolutionResult;

/**
 * Derives multi-year survival probability for a precarious balanced
 * rock (PBR), from a {@link FailureProbabilityAsset}'s annual probability
 * of failure and the rock's age.
 *
 * <p>PBRs that have stood for a known age without toppling constrain the
 * annual failure probability: surviving {@code age} years at a constant
 * annual probability {@code p} implies survival probability
 * {@code (1-p)^age} ({@link #getProbabilityOfSurvival()}). Conversely,
 * {@link #getHazardAdjustment} answers the inverse question -- by what
 * factor would {@code p} need to scale for the probability of survival
 * over {@code age} years to equal a given target?</p>
 */
public class PBRSurvivalAsset extends RiskConvolutionAsset {

    private final FailureProbabilityAsset base;
    private final int age;
    private final double annualProbabilityOfFailure;

    public PBRSurvivalAsset(FailureProbabilityAsset base) {

        super(
            base.getAssetID(),
            base.getLatitude(),
            base.getLongitude(),
            base.getVs30(),
            base.getModelName(),
            base.getAdditionalFields(),
            AssetType.FRAGILITY,
            RiskMetricType.PBR_SURVIVAL
        );

        this.base = base;

        // require age
        String ageStr = base.getAdditionalFields().get("age");
        if (ageStr == null || ageStr.isBlank()) {
            throw new IllegalArgumentException(
                "PBRSurvivalAsset requires additional field 'age'");
        }
        this.age = parseKaToYears(ageStr);

        // require probability of failure
        double p = base.getProbabilityOfFailure();
        if (Double.isNaN(p)) {
            throw new IllegalStateException(
                "Probability of failure must be set before creating PBRSurvivalAsset");
        }
        this.annualProbabilityOfFailure = p;
    }


    public FailureProbabilityAsset getBase() {
        return base;
    }

    public double getAnnualProbabilityOfFailure() {
        return annualProbabilityOfFailure;
    }

    public int getAge() {return age;}

    public double getProbabilityOfSurvival() {
    	// (1 - annual_p)^T, using a more computationally stable form than Math.pow(1-p, T)
        return Math.exp(age * Math.log1p(-annualProbabilityOfFailure));
    }
    

    public double getHazardAdjustment(double probabilityTarget) {
    	// more computationally stable form than (1.0 - Math.pow(P1, 1.0 / T)) / P;
    	double x = Math.log(probabilityTarget) / age;
        return -Math.expm1(x) / annualProbabilityOfFailure;
    }
    


    
    @Override
    public void setRiskConvolutionResult(ConvolutionResult result) {
        throw new UnsupportedOperationException(
            "PBRSurvivalAsset does not support setting a risk convolution result");
    }

    @Override
    public ConvolutionResult getRiskConvolutionResult() {
        throw new UnsupportedOperationException(
            "PBRSurvivalAsset does not expose a risk convolution result");
    }
    

    public static int parseKaToYears(String kaValue) {

        if (kaValue == null || kaValue.isBlank()) {
            throw new IllegalArgumentException("Age (ka) is missing");
        }

        kaValue = kaValue.trim().replace("\"", "");

        double ka = Double.parseDouble(kaValue);

        return (int) Math.round(ka * 1000.0);
    }
    

}
