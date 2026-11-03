package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.fragility.FailureProbabilityAsset;
import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.asset.vulnerability.ExpectedLossAsset;
import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.AssetHazardCurve;
import scratch.anne.sra.domain.hazard.HazardCurve;
import scratch.anne.sra.domain.hazard.HazardCurveCollection;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.domain.structural_response.SimpleImResponse;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.engine.per_rupture.RuptureLossResultsCollection;
import scratch.anne.sra.io.writers.FileFormat;
import scratch.anne.sra.io.writers.HazardCurvesExporter;
import scratch.anne.sra.io.writers.RuptureResultsExporter;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.PortfolioGroupingUtils;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

import java.nio.file.Path;

/**
 * Projection of a base Portfolio into a calculation-ready working set of RiskConvolutionAssets.
 *
 * <p>This class converts heterogeneous portfolio assets into a unified computational form and
 * builds deterministic lookup structures required for hazard/risk convolution workflows.</p>
 *
 * <p>The portfolio is immutable in structure after construction. All grouping maps are derived
 * from canonical per-asset key assignments.</p>
 */
public final class RiskConvolutionPortfolio implements PortfolioGetters<RiskConvolutionAsset> {

    // ---------------------------------------------------------------------
    // Core assets
    // ---------------------------------------------------------------------

    private final List<RiskConvolutionAsset> assets;
    private final Map<Integer, RiskConvolutionAsset> assetsById;
    private final AssetType assetType;
    private final RiskMetricType riskMetricType;

    // ---------------------------------------------------------------------
    // Canonical per-asset keys (single source of truth)
    // ---------------------------------------------------------------------

    private final Map<RiskConvolutionAsset, ImKey> imKeyByAsset;
    private final Map<RiskConvolutionAsset, SiteKey> siteKeyByAsset;

    // ---------------------------------------------------------------------
    // Derived indices
    // ---------------------------------------------------------------------

    private final Map<ImKey, List<RiskConvolutionAsset>> assetsByImKey;
    private final Map<SiteKey, List<RiskConvolutionAsset>> assetsBySite;
    private final Map<SiteKey, Map<ImKey, List<RiskConvolutionAsset>>> assetsBySiteAndImKey;
    private final Map<SiteKey, Set<ImKey>> imKeysBySite;

    private final Map<String, Map<String, List<RiskConvolutionAsset>>> assetsByAdditionalField;
    private final List<String> additionalFieldNames;

    private final Metadata metadata;
    
    private ConvolutionMode convolutionMode = ConvolutionMode.UNCOMPUTED;
    
	 // ---------------------------------------------------------------------
	 // Hazard curves (used when ConvolutionMode.FULL_HCURVE)
	 // ---------------------------------------------------------------------
	
    private HazardCurveCollection hazardCurves;
    
	 // ---------------------------------------------------------------------
	 // Rupture results (used when ConvolutionMode.PER_RUPTURE)
	 // ---------------------------------------------------------------------
	
   private RuptureLossResultsCollection ruptureResults;

    // ---------------------------------------------------------------------
    // Construction
    // ---------------------------------------------------------------------

    public RiskConvolutionPortfolio(
            Portfolio<? extends AbstractAsset> base,
            SimpleImResponseLibrary responseLibrary
    ) {
        Objects.requireNonNull(base);
        Objects.requireNonNull(responseLibrary);


        // Build IM lookup: modelName → IMKey
        Map<String, ImKey> responseToIm = new HashMap<>();
        for (SimpleImResponse r : responseLibrary.getModels()) {
            responseToIm.put(r.getName(), r.getImKey());
        }

        List<RiskConvolutionAsset> projected = new ArrayList<>(base.size());

        Map<RiskConvolutionAsset, ImKey> tmpIm = new LinkedHashMap<>();
        Map<RiskConvolutionAsset, SiteKey> tmpSite = new LinkedHashMap<>();

        for (AbstractAsset a : base.getAssets()) {

            RiskConvolutionAsset riskAsset = convertToRiskAsset(a);

            ImKey imKey = responseToIm.get(a.getModelName());
            if (imKey == null) {
                throw new IllegalArgumentException(
                        "No IMKey for model: " + a.getModelName()
                );
            }

            SiteKey siteKey = new SiteKey(
                    a.getLatitude(),
                    a.getLongitude(),
                    a.getVs30()
            );

            riskAsset.setImKey(imKey);

            projected.add(riskAsset);
            tmpIm.put(riskAsset, imKey);
            tmpSite.put(riskAsset, siteKey);
        }

        this.assets = List.copyOf(projected);
        
        validateHomogeneous(assets);
        this.assetType = assets.get(0).getAssetType();
        this.riskMetricType = assets.get(0).getRiskMetricType();
        
        this.imKeyByAsset = Collections.unmodifiableMap(tmpIm);
        this.siteKeyByAsset = Collections.unmodifiableMap(tmpSite);

        // -----------------------------------------------------------------
        // Index: asset ID
        // -----------------------------------------------------------------
        Map<Integer, RiskConvolutionAsset> idMap = new LinkedHashMap<>();
        for (RiskConvolutionAsset a : assets) {
            idMap.put(a.getAssetID(), a);
        }
        this.assetsById = Collections.unmodifiableMap(idMap);

        // -----------------------------------------------------------------
        // Index: IM
        // -----------------------------------------------------------------
        Map<ImKey, List<RiskConvolutionAsset>> imMap = new HashMap<>();
        for (Map.Entry<RiskConvolutionAsset, ImKey> e : imKeyByAsset.entrySet()) {
            imMap.computeIfAbsent(e.getValue(), k -> new ArrayList<>())
                 .add(e.getKey());
        }
        this.assetsByImKey = imMap.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> List.copyOf(e.getValue())
                ));

        // -----------------------------------------------------------------
        // Index: Site
        // -----------------------------------------------------------------
        Map<SiteKey, List<RiskConvolutionAsset>> siteMap = new HashMap<>();
        for (Map.Entry<RiskConvolutionAsset, SiteKey> e : siteKeyByAsset.entrySet()) {
            siteMap.computeIfAbsent(e.getValue(), k -> new ArrayList<>())
                   .add(e.getKey());
        }
        this.assetsBySite = siteMap.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> List.copyOf(e.getValue())
                ));

        // -----------------------------------------------------------------
        // Index: Site → IM
        // -----------------------------------------------------------------
        Map<SiteKey, Map<ImKey, List<RiskConvolutionAsset>>> siteIm = new HashMap<>();

        for (RiskConvolutionAsset a : assets) {

            SiteKey site = siteKeyByAsset.get(a);
            ImKey im = imKeyByAsset.get(a);

            siteIm
                    .computeIfAbsent(site, k -> new HashMap<>())
                    .computeIfAbsent(im, k -> new ArrayList<>())
                    .add(a);
        }

        this.assetsBySiteAndImKey = siteIm.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().entrySet().stream()
                                .collect(Collectors.toMap(
                                        Map.Entry::getKey,
                                        v -> List.copyOf(v.getValue())
                                ))
                ));

        // -----------------------------------------------------------------
        // Index: Site → IM keys
        // -----------------------------------------------------------------
        Map<SiteKey, Set<ImKey>> tmpKeysBySite = new HashMap<>();

        for (RiskConvolutionAsset a : assets) {
            SiteKey site = siteKeyByAsset.get(a);
            ImKey im = imKeyByAsset.get(a);

            tmpKeysBySite
                    .computeIfAbsent(site, k -> new LinkedHashSet<>())
                    .add(im);
        }

        this.imKeysBySite = tmpKeysBySite.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> Set.copyOf(e.getValue())
                ));

        // -----------------------------------------------------------------
        // Additional fields
        // -----------------------------------------------------------------
        this.assetsByAdditionalField =
                PortfolioGroupingUtils.invertNestedFields(
                        assets,
                        RiskConvolutionAsset::getAdditionalFields
                );

        this.additionalFieldNames =
                List.copyOf(assetsByAdditionalField.keySet());

        this.metadata = base.getMetadata();
    }

    // ---------------------------------------------------------------------
    // Factory conversion
    // ---------------------------------------------------------------------

    private static RiskConvolutionAsset convertToRiskAsset(AbstractAsset asset) {

        if (asset instanceof VulnerabilityAsset v)
            return new ExpectedLossAsset(v);

        if (asset instanceof FragilityAsset f)
            return new FailureProbabilityAsset(f);

        throw new IllegalArgumentException(
                "Unsupported asset type: " + asset.getClass()
        );
    }


    // ---------------------------------------------------------------------
    // PortfolioGetters
    // ---------------------------------------------------------------------

    @Override
    public List<RiskConvolutionAsset> getAssets() {
        return assets;
    }

    @Override
    public Set<Integer> getAssetIDs() {
        return assetsById.keySet();
    }
    
    @Override
    public AssetType getAssetType() {
    	return assetType;
    }
    
    public RiskMetricType getRiskMetricType() {
    	return riskMetricType;
    }

    @Override
    public RiskConvolutionAsset getAssetByID(int id) {
        return assetsById.get(id);
    }

    @Override
    public int size() {
        return assets.size();
    }

    @Override
    public Map<SiteKey, List<RiskConvolutionAsset>> getSiteMap() {
        return assetsBySite;
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return assetsBySite.keySet();
    }
    
    @Override
    public SiteKey getSiteKey(RiskConvolutionAsset asset) {
        return siteKeyByAsset.get(asset);
    }

    @Override
    public List<RiskConvolutionAsset> getAssetsBySite(SiteKey siteKey) {
        return assetsBySite.getOrDefault(siteKey, List.of());
    }

    public Set<ImKey> getImKeys() {
        return assetsByImKey.keySet();
    }

    public List<RiskConvolutionAsset> getAssetsByImKey(ImKey key) {
        return assetsByImKey.getOrDefault(key, List.of());
    }

    public List<RiskConvolutionAsset> getAssetsBySiteAndImKey(
            SiteKey site,
            ImKey imKey
    ) {
        return assetsBySiteAndImKey
                .getOrDefault(site, Map.of())
                .getOrDefault(imKey, List.of());
    }

    public Set<ImKey> getImKeysBySite(SiteKey siteKey) {
        return imKeysBySite.getOrDefault(siteKey, Set.of());
    }

    @Override
    public List<String> getAdditionalFieldNames() {
        return additionalFieldNames;
    }

    @Override
    public List<String> getAdditionalFieldValues(String field) {
        return List.copyOf(
                assetsByAdditionalField
                        .getOrDefault(field, Map.of())
                        .keySet()
        );
    }

    @Override
    public List<RiskConvolutionAsset> getAssetsByAdditionalField(
            String field,
            String value
    ) {
        return assetsByAdditionalField
                .getOrDefault(field, Map.of())
                .getOrDefault(value, List.of());
    }

    @Override
    public Map<String, List<RiskConvolutionAsset>> getAdditionalFieldMap(
            String field
    ) {
        return assetsByAdditionalField.getOrDefault(field, Map.of());
    }

    @Override
    public Metadata getMetadata() {
        return metadata;
    }

    @Override
    public Set<String> getResponseModelNames() {
        return assets.stream()
                .map(RiskConvolutionAsset::getModelName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }
    
    
    // ----- Convolution mode for portfolio state ------
    public enum ConvolutionMode {
	    UNCOMPUTED,
	    FULL_HCURVE,
	    PER_RUPTURE
	}
    
    public void setConvolutionMode(ConvolutionMode convolutionMode) {
        if (this.convolutionMode != ConvolutionMode.UNCOMPUTED && this.convolutionMode != convolutionMode) {
            throw new IllegalStateException("Portfolio already used in a different convolution mode");
        }
        this.convolutionMode = convolutionMode;
    }
    
    public void assertConvolutionExecutedAs(ConvolutionMode expected) {
        if (this.convolutionMode != expected) {
            throw new IllegalStateException(
                "Expected risk convolution mode " + expected + " but was " + convolutionMode
            );
        }
    }
    
    // ----- hazard storage helpers ---------   
    public void setHazardCurves(HazardCurveCollection curves) {
        this.hazardCurves = curves;
    }
    
    
    public HazardCurveCollection getHazardCurves() {
        return hazardCurves;
    }
    
    
    public void exportHazard(Path file) {
        new HazardCurvesExporter()
                .export(hazardCurves, file);
    }

    
    public AssetHazardCurve getHazardForAsset(int id) {

        RiskConvolutionAsset asset = assetsById.get(id);

        SiteKey siteKey = siteKeyByAsset.get(asset);
        if (siteKey == null) {
            throw new IllegalArgumentException(
                    "Unknown asset: " + asset.getAssetID());
        }

        ImKey imKey = imKeyByAsset.get(asset);
        if (imKey == null) {
            throw new IllegalStateException(
                    "No IMKey for asset: " + asset.getAssetID());
        }

        HazardCurve curve = hazardCurves.get(siteKey, imKey);

        if (curve == null) {
            throw new IllegalStateException(
                    "No hazard stored for asset " + asset.getAssetID());
        }

        return new AssetHazardCurve(
                siteKey,
                imKey,
                imKey.getLinearValues(),
                curve.getHazard()
        );
    }
    
    // ----- rupture results storage helpers ---------   
    public void setRuptureResults(RuptureLossResultsCollection ruptureResults) {
        this.ruptureResults = ruptureResults;
    }
    
    
    public RuptureLossResultsCollection getRuptureResults() {
        return ruptureResults;
    }
    
    
    public void exportRuptureResults(Path file) {
        exportRuptureResults(file, FileFormat.CSV);
    }

    public void exportRuptureResults(Path file, FileFormat format) {
        new RuptureResultsExporter()
                .export(ruptureResults, file, format);
    }
    
    private void validateHomogeneous(List<RiskConvolutionAsset> assets) {
        if (assets == null || assets.isEmpty()) {
            return;
        }

        RiskMetricType type = assets.get(0).getRiskMetricType();

        for (RiskConvolutionAsset asset : assets) {
            if (asset.getRiskMetricType() != type) {
                throw new IllegalArgumentException(
                    "Portfolio must contain homogeneous risk metric types. Expected " 
                    + type + " but found " + asset.getRiskMetricType()
                );
            }
        }
    }
    
    
    // ---------------------------------------------------------------------
    // Asset value queries (EXPECTED_LOSS portfolios only)
    // ---------------------------------------------------------------------
    
    /**
     * Total value of all assets in the portfolio.
     *
     * @throws IllegalStateException if this is not an EXPECTED_LOSS portfolio
     */
    public double getTotalPortfolioValue() {
        requireExpectedLoss("getTotalPortfolioValue()");
        return assets.stream()
                .mapToDouble(a -> ((ExpectedLossAsset) a).getValue())
                .sum();
    }

    /**
     * Total value of all assets at the given site, across all IMs.
     *
     * @throws IllegalStateException if this is not an EXPECTED_LOSS portfolio
     */
    public double getTotalSiteValue(SiteKey siteKey) {
        requireExpectedLoss("getTotalAssetValue(SiteKey)");
        return assetsBySite.getOrDefault(siteKey, List.of())
                .stream()
                .mapToDouble(a -> ((ExpectedLossAsset) a).getValue())
                .sum();
    }

    /**
     * Total value of assets at the given site sensitive to the given IM.
     *
     * @throws IllegalStateException if this is not an EXPECTED_LOSS portfolio
     */
    public double getTotalSiteImValue(SiteKey siteKey, ImKey imKey) {
        requireExpectedLoss("getTotalAssetValue(SiteKey, ImKey)");
        return assetsBySiteAndImKey
                .getOrDefault(siteKey, Map.of())
                .getOrDefault(imKey, List.of())
                .stream()
                .mapToDouble(a -> ((ExpectedLossAsset) a).getValue())
                .sum();
    }

    private void requireExpectedLoss(String method) {
        if (riskMetricType != RiskMetricType.EXPECTED_LOSS) {
            throw new IllegalStateException(
                method + " is only valid for EXPECTED_LOSS portfolios, but this portfolio is "
                + riskMetricType);
        }
    }
}