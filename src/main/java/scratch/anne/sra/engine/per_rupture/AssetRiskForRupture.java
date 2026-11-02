package scratch.anne.sra.engine.per_rupture;

public record AssetRiskForRupture(
	    String assetId,
	    int sourceId,
	    int ruptureId,
	    double conditionalRisk
	) {}
