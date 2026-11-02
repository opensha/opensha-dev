package scratch.anne.sra.domain.asset.fragility;

import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.engine.convolution.ConvolutionResult;

public class FailureProbabilityAsset extends RiskConvolutionAsset{

    private final FragilityAsset asset;
    
    private ConvolutionResult riskResult;

    private double probabilityOfFailure = Double.NaN;

    public FailureProbabilityAsset(FragilityAsset asset) {

        super(
            asset.getAssetID(),
            asset.getLatitude(),
            asset.getLongitude(),
            asset.getVs30(),
            asset.getModelName(),
            asset.getAdditionalFields(),
            AssetType.FRAGILITY,
            RiskMetricType.FAILURE_PROBABILITY
        );

        this.asset = asset;
    }

    public FragilityAsset getAsset() {
        return asset;
    }
    
    public RiskMetricType getRiskMetricType() {
    	return RiskMetricType.FAILURE_PROBABILITY;
    }

    public double getProbabilityOfFailure() {
        return probabilityOfFailure;
    }

    /** wrap the loss setter for generic workflows */
	@Override
	public void setRiskConvolutionResult(ConvolutionResult result) {
		this.riskResult = result;
		this.probabilityOfFailure = result.getRisk();
		
	}

	/** wrap the loss getter for generic workflows */
	@Override
	public ConvolutionResult getRiskConvolutionResult() {
		return this.riskResult;
	}
	
}
