package scratch.anne.sra.engine.perRupture;

/** One asset's conditional risk (loss or exceedance probability) for one rupture, as written by {@link scratch.anne.sra.io.writers.AssetRiskForRuptureWriter}. */
public record AssetRiskForRupture(
	    int assetId,
	    int sourceId,
	    int ruptureId,
	    double conditionalRisk
	) {}
