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
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Reporting view over a RiskConvolutionPortfolio whose assets are
 * already materialized FailureProbabilityAsset instances.
 *
 * <p>
 * This class does NOT transform or copy asset state. It only provides
 * export utilities, delegating the generic asset-projection plumbing
 * to {@link RiskConvolutionPortfolioView}.
 * </p>
 */
public final class ProbFailurePortfolioReporter
        implements Iterable<FailureProbabilityAsset>, PortfolioGetters<FailureProbabilityAsset> {

    private final RiskConvolutionPortfolioView<FailureProbabilityAsset> view;

    public ProbFailurePortfolioReporter(RiskConvolutionPortfolio base) {
        this.view = new RiskConvolutionPortfolioView<>(base, a -> (FailureProbabilityAsset) a);
    }

    // ---------------------------------------------------------------------
    // PortfolioGetters (DELEGATED TO view)
    // ---------------------------------------------------------------------

    @Override public List<FailureProbabilityAsset> getAssets() { return view.getAssets(); }
    @Override public Set<Integer> getAssetIDs() { return view.getAssetIDs(); }
    @Override public FailureProbabilityAsset getAssetByID(int assetID) { return view.getAssetByID(assetID); }
    @Override public AssetType getAssetType() { return view.getAssetType(); }
    public RiskMetricType getRiskMetricType() { return view.getRiskMetricType(); }
    @Override public int size() { return view.size(); }
    @Override public Map<SiteKey, List<FailureProbabilityAsset>> getSiteMap() { return view.getSiteMap(); }
    @Override public Set<SiteKey> getSiteKeys() { return view.getSiteKeys(); }
    @Override public SiteKey getSiteKey(FailureProbabilityAsset asset) { return view.getSiteKey(asset); }
    @Override public List<FailureProbabilityAsset> getAssetsBySite(SiteKey siteKey) { return view.getAssetsBySite(siteKey); }
    public Set<ImKey> getImKeys() { return view.getImKeys(); }
    public List<FailureProbabilityAsset> getAssetsByImKey(ImKey key) { return view.getAssetsByImKey(key); }
    public List<FailureProbabilityAsset> getAssetsBySiteAndImKey(SiteKey site, ImKey imKey) { return view.getAssetsBySiteAndImKey(site, imKey); }
    public Set<ImKey> getImKeysBySite(SiteKey site) { return view.getImKeysBySite(site); }
    @Override public List<String> getAdditionalFieldNames() { return view.getAdditionalFieldNames(); }
    @Override public List<String> getAdditionalFieldValues(String field) { return view.getAdditionalFieldValues(field); }
    @Override public List<FailureProbabilityAsset> getAssetsByAdditionalField(String field, String value) { return view.getAssetsByAdditionalField(field, value); }
    @Override public Map<String, List<FailureProbabilityAsset>> getAdditionalFieldMap(String field) { return view.getAdditionalFieldMap(field); }
    @Override public Set<String> getResponseModelNames() { return view.getResponseModelNames(); }
    @Override public Metadata getMetadata() { return view.getMetadata(); }
    @Override public Iterator<FailureProbabilityAsset> iterator() { return view.getAssets().iterator(); }

    // ---------------------------------------------------------------------
    // CSV EXPORT
    // ---------------------------------------------------------------------

    public void writeCSV(Path file) throws IOException {

        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(file))) {

            pw.print("AssetID,Latitude,Longitude,Vs30,Model");

            for (String f : getAdditionalFieldNames()) {
                pw.print("," + f);
            }

            pw.println(",ProbFailure");

            for (FailureProbabilityAsset a : view.getAssets()) {

                pw.printf(
                        "%d,%.6f,%.6f,%.0f,\"%s\"",
                        a.getAssetID(),
                        a.getLatitude(),
                        a.getLongitude(),
                        a.getVs30(),
                        a.getModelName()
                );

                Map<String, String> extra = a.getAdditionalFields();

                for (String f : getAdditionalFieldNames()) {
                    pw.printf(",\"%s\"", extra.getOrDefault(f, ""));
                }

                pw.printf(",%.6e%n", a.getProbabilityOfFailure());

            }
        }
    }

    public void writeGroupedCSV(Path file) throws IOException {

        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(file))) {

            pw.println(
                    "GroupField,GroupValue,GroupFieldNumAssets," +
                    "AvgProbFailure,Min,Max,StdDev"
            );

            for (String field : getAdditionalFieldNames()) {

                for (Map.Entry<String, List<FailureProbabilityAsset>> e :
                        getAdditionalFieldMap(field).entrySet()) {

                    List<FailureProbabilityAsset> group = e.getValue();

                    Map<String, Double> stats =
                            RiskConvolutionPortfolioView.summarize(group, FailureProbabilityAsset::getProbabilityOfFailure);

                    pw.printf(
                            "\"%s\",\"%s\",\"%d\",%.6f,%.6f,%.6f,%.6f%n",
                            field,
                            e.getKey(),
                            stats.get("n").intValue(),
                            stats.get("avg"),
                            stats.get("min"),
                            stats.get("max"),
                            stats.get("stdDev")
                    );
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // UTILITIES
    // ---------------------------------------------------------------------

    /** Prints a summary of portfolio probabilities of failure including per-asset detail. */
    public void printSummary() {
        System.out.println("----- ProbFailurePortfolio Summary -----");
        System.out.printf("Total assets: %d%n", view.size());
    }
}
