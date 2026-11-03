package scratch.anne.sra.engine.per_rupture;

public record AssetRiskForRupture(
	    int assetId,
	    int sourceId,
	    int ruptureId,
	    double conditionalRisk
	) {}
