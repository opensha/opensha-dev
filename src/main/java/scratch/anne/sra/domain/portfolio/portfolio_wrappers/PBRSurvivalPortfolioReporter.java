package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.fragility.FailureProbabilityAsset;
import scratch.anne.sra.domain.asset.fragility.PBRSurvivalAsset;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Reporting view producing PBR survival results.
 *
 * <p>Wraps each {@code FailureProbabilityAsset} into a
 * {@code PBRSurvivalAsset} exactly once (see
 * {@link RiskConvolutionPortfolioView}'s own javadoc on why this
 * matters), and exposes CSV export utilities.</p>
 */
public final class PBRSurvivalPortfolioReporter
        implements Iterable<PBRSurvivalAsset>, PortfolioGetters<PBRSurvivalAsset> {

    private final RiskConvolutionPortfolioView<PBRSurvivalAsset> view;
    private final List<Double> probabilityTargets;

    public PBRSurvivalPortfolioReporter(RiskConvolutionPortfolio base, List<Double> probabilityTargets) {
        // validation happens in PBRSurvivalAsset's own constructor
        this.view = new RiskConvolutionPortfolioView<>(base, a -> new PBRSurvivalAsset((FailureProbabilityAsset) a));
        this.probabilityTargets = List.copyOf(probabilityTargets);
    }

    // -------- default constructor --------
    public PBRSurvivalPortfolioReporter(RiskConvolutionPortfolio base) {
        this(base, List.of(0.05));
    }

    // ---------------------------------------------------------------------
    // PortfolioGetters (DELEGATED TO view)
    // ---------------------------------------------------------------------

    @Override public List<PBRSurvivalAsset> getAssets() { return view.getAssets(); }
    @Override public Set<Integer> getAssetIDs() { return view.getAssetIDs(); }
    @Override public PBRSurvivalAsset getAssetByID(int assetID) { return view.getAssetByID(assetID); }
    @Override public AssetType getAssetType() { return view.getAssetType(); }
    public RiskMetricType getRiskMetricType() { return view.getRiskMetricType(); }
    @Override public int size() { return view.size(); }
    @Override public Map<SiteKey, List<PBRSurvivalAsset>> getSiteMap() { return view.getSiteMap(); }
    @Override public Set<SiteKey> getSiteKeys() { return view.getSiteKeys(); }
    @Override public SiteKey getSiteKey(PBRSurvivalAsset asset) { return view.getSiteKey(asset); }
    @Override public List<PBRSurvivalAsset> getAssetsBySite(SiteKey siteKey) { return view.getAssetsBySite(siteKey); }
    public Set<ImKey> getImKeys() { return view.getImKeys(); }
    public List<PBRSurvivalAsset> getAssetsByImKey(ImKey key) { return view.getAssetsByImKey(key); }
    public List<PBRSurvivalAsset> getAssetsBySiteAndImKey(SiteKey site, ImKey imKey) { return view.getAssetsBySiteAndImKey(site, imKey); }
    public Set<ImKey> getImKeysBySite(SiteKey site) { return view.getImKeysBySite(site); }
    @Override public List<String> getAdditionalFieldNames() { return view.getAdditionalFieldNames(); }
    @Override public List<String> getAdditionalFieldValues(String field) { return view.getAdditionalFieldValues(field); }
    @Override public List<PBRSurvivalAsset> getAssetsByAdditionalField(String field, String value) { return view.getAssetsByAdditionalField(field, value); }
    @Override public Map<String, List<PBRSurvivalAsset>> getAdditionalFieldMap(String field) { return view.getAdditionalFieldMap(field); }
    @Override public Set<String> getResponseModelNames() { return view.getResponseModelNames(); }
    @Override public Metadata getMetadata() { return view.getMetadata(); }
    @Override public Iterator<PBRSurvivalAsset> iterator() { return view.getAssets().iterator(); }

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

            for (PBRSurvivalAsset a : view.getAssets()) {

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
}
