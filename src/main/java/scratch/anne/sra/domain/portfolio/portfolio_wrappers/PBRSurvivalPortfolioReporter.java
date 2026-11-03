package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.fragility.*;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Reporting view producing PBR survival results.
 *
 * Wraps FailureProbabilityAssets into PBRSurvivalAssets
 * and exposes CSV export utilities.
 */
public final class PBRSurvivalPortfolioReporter
        implements Iterable<PBRSurvivalAsset>,
                   PortfolioGetters<PBRSurvivalAsset> {

    // ---------------------------------------------------------------------
    // BASE PORTFOLIO
    // ---------------------------------------------------------------------

    private final RiskConvolutionPortfolio base;
    private final List<PBRSurvivalAsset> assets;
    private final List<Double> probabilityTargets;

    // ---------------------------------------------------------------------
    // CONSTRUCTOR
    // ---------------------------------------------------------------------

    public PBRSurvivalPortfolioReporter(
            RiskConvolutionPortfolio base,
            List<Double> probabilityTargets
    ) {

        this.base = Objects.requireNonNull(base);

        base.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);

        this.probabilityTargets = List.copyOf(probabilityTargets);

        this.assets = base.getAssets().stream()
                .map(a -> (FailureProbabilityAsset) a)
                .map(PBRSurvivalAsset::new)   // validation happens here
                .toList();
    }
    
    // -------- default constructor --------
    public PBRSurvivalPortfolioReporter(RiskConvolutionPortfolio base) {
        this(base, List.of(0.05));
    }

    // ---------------------------------------------------------------------
    // BASIC ACCESS
    // ---------------------------------------------------------------------

    @Override
    public List<PBRSurvivalAsset> getAssets() {
        return assets;
    }

    @Override
    public int size() {
        return assets.size();
    }

    @Override
    public Set<Integer> getAssetIDs() {
        return base.getAssetIDs();
    }
    
    @Override
    public AssetType getAssetType() {
    	// risk convolution portfolio already confirmed homogeneity
    	return assets.get(0).getAssetType();
    }
    
    public RiskMetricType getRiskMetricType() {
    	// risk convolution portfolio already confirmed homogeneity
    	return assets.get(0).getRiskMetricType();
    }

    @Override
    public PBRSurvivalAsset getAssetByID(int id) {
        return assets.stream()
                .filter(a -> a.getAssetID() == id)
                .findFirst()
                .orElse(null);
    }

    // ---------------------------------------------------------------------
    // SITE GROUPING
    // ---------------------------------------------------------------------

    @Override
    public Map<SiteKey, List<PBRSurvivalAsset>> getSiteMap() {
        return base.getSiteMap().entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> new PBRSurvivalAsset(
                                        (FailureProbabilityAsset) a))
                                .toList()
                ));
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return base.getSiteKeys();
    }
    
    @Override
    public SiteKey getSiteKey(PBRSurvivalAsset asset) {
        return base.getSiteKey(asset);
    }

    @Override
    public List<PBRSurvivalAsset> getAssetsBySite(SiteKey key) {
        return base.getAssetsBySite(key).stream()
                .map(a -> new PBRSurvivalAsset((FailureProbabilityAsset) a))
                .toList();
    }

    // ---------------------------------------------------------------------
    // IM GROUPING
    // ---------------------------------------------------------------------

    public Set<ImKey> getImKeys() {
        return base.getImKeys();
    }

    // ---------------------------------------------------------------------
    // ADDITIONAL FIELDS
    // ---------------------------------------------------------------------

    @Override
    public List<String> getAdditionalFieldNames() {
        return base.getAdditionalFieldNames();
    }

    @Override
    public List<String> getAdditionalFieldValues(String field) {
        return base.getAdditionalFieldValues(field);
    }
    
    @Override
    public List<PBRSurvivalAsset> getAssetsByAdditionalField(
            String field,
            String value
    ) {
        return base.getAssetsByAdditionalField(field, value).stream()
                .map(a -> (PBRSurvivalAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    @Override
    public Map<String, List<PBRSurvivalAsset>> getAdditionalFieldMap(String field) {
        return base.getAdditionalFieldMap(field).entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> (PBRSurvivalAsset) a)
                                .collect(java.util.stream.Collectors.toList())
                ));
    }

    // ---------------------------------------------------------------------
    // METADATA
    // ---------------------------------------------------------------------

    @Override
    public Metadata getMetadata() {
        return base.getMetadata();
    }

    // ---------------------------------------------------------------------
    // ITERATION
    // ---------------------------------------------------------------------

    @Override
    public Iterator<PBRSurvivalAsset> iterator() {
        return assets.iterator();
    }

    // ---------------------------------------------------------------------
    // CSV EXPORT
    // ---------------------------------------------------------------------

    public void writeCSV(Path file) throws IOException {

        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(file))) {

            // ----- header -----

            pw.print("AssetID,Latitude,Longitude,Vs30,Model");

            for (String f : getAdditionalFieldNames())
                pw.print("," + f);

            pw.print(",AnnualProbFailure,Survival");

            for (double pt : probabilityTargets)
                pw.printf(",HazAdj_P%.3f", pt);

            pw.println();

            // ----- rows -----

            for (PBRSurvivalAsset a : assets) {

                pw.printf(
                        "%d,%.6f,%.6f,%.0f,\"%s\"",
                        a.getAssetID(),
                        a.getLatitude(),
                        a.getLongitude(),
                        a.getVs30(),
                        a.getModelName()
                );

                Map<String,String> extra = a.getAdditionalFields();

                for (String f : getAdditionalFieldNames())
                    pw.printf(",\"%s\"", extra.getOrDefault(f, ""));

                pw.printf(",%.6e", a.getAnnualProbabilityOfFailure());
                pw.printf(",%.6e", a.getProbabilityOfSurvival());

                for (double pt : probabilityTargets)
                    pw.printf(",%.6e", a.getHazardAdjustment(pt));

                pw.println();
            }
        }
    }

    // ---------------------------------------------------------------------
    // RESPONSE MODELS
    // ---------------------------------------------------------------------

    @Override
    public Set<String> getResponseModelNames() {
        return base.getResponseModelNames();
    }
}
