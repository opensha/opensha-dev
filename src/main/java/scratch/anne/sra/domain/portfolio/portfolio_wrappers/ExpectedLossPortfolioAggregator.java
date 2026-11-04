package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.vulnerability.ExpectedLossAsset;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Aggregation + reporting view over a RiskConvolutionPortfolio whose
 * assets are already materialized ExpectedLossAsset instances.
 *
 * <p>
 * This class does NOT transform or copy asset state. It only provides
 * aggregated views and export utilities, delegating the generic
 * asset-projection plumbing to {@link RiskConvolutionPortfolioView}.
 * </p>
 */
public final class ExpectedLossPortfolioAggregator
        implements Iterable<ExpectedLossAsset>, PortfolioGetters<ExpectedLossAsset> {

    private final RiskConvolutionPortfolioView<ExpectedLossAsset> view;

    public ExpectedLossPortfolioAggregator(RiskConvolutionPortfolio base) {
        this.view = new RiskConvolutionPortfolioView<>(base, a -> (ExpectedLossAsset) a);
    }

    // ---------------------------------------------------------------------
    // PortfolioGetters (DELEGATED TO view)
    // ---------------------------------------------------------------------

    @Override public List<ExpectedLossAsset> getAssets() { return view.getAssets(); }
    @Override public Set<Integer> getAssetIDs() { return view.getAssetIDs(); }
    @Override public ExpectedLossAsset getAssetByID(int assetID) { return view.getAssetByID(assetID); }
    @Override public AssetType getAssetType() { return view.getAssetType(); }
    public RiskMetricType getRiskMetricType() { return view.getRiskMetricType(); }
    @Override public int size() { return view.size(); }
    @Override public Map<SiteKey, List<ExpectedLossAsset>> getSiteMap() { return view.getSiteMap(); }
    @Override public Set<SiteKey> getSiteKeys() { return view.getSiteKeys(); }
    @Override public SiteKey getSiteKey(ExpectedLossAsset asset) { return view.getSiteKey(asset); }
    @Override public List<ExpectedLossAsset> getAssetsBySite(SiteKey siteKey) { return view.getAssetsBySite(siteKey); }
    public Set<ImKey> getImKeys() { return view.getImKeys(); }
    public List<ExpectedLossAsset> getAssetsByImKey(ImKey key) { return view.getAssetsByImKey(key); }
    public List<ExpectedLossAsset> getAssetsBySiteAndImKey(SiteKey site, ImKey imKey) { return view.getAssetsBySiteAndImKey(site, imKey); }
    public Set<ImKey> getImKeysBySite(SiteKey site) { return view.getImKeysBySite(site); }
    @Override public List<String> getAdditionalFieldNames() { return view.getAdditionalFieldNames(); }
    @Override public List<String> getAdditionalFieldValues(String field) { return view.getAdditionalFieldValues(field); }
    @Override public List<ExpectedLossAsset> getAssetsByAdditionalField(String field, String value) { return view.getAssetsByAdditionalField(field, value); }
    @Override public Map<String, List<ExpectedLossAsset>> getAdditionalFieldMap(String field) { return view.getAdditionalFieldMap(field); }
    @Override public Set<String> getResponseModelNames() { return view.getResponseModelNames(); }
    @Override public Metadata getMetadata() { return view.getMetadata(); }
    @Override public Iterator<ExpectedLossAsset> iterator() { return view.getAssets().iterator(); }

    // ---------------------------------------------------------------------
    // AGGREGATIONS
    // ---------------------------------------------------------------------

    public double getTotalExpectedLoss() {
        return view.getAssets().stream()
                .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                .sum();
    }

    public double getTotalAssetValue() {
        return view.getAssets().stream()
                .mapToDouble(ExpectedLossAsset::getValue)
                .sum();
    }

    public Map<SiteKey, Double> getExpectedLossBySite() {
        Map<SiteKey, Double> result = new LinkedHashMap<>();
        for (SiteKey s : getSiteKeys()) {
            result.put(
                    s,
                    getAssetsBySite(s).stream()
                            .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                            .sum()
            );
        }
        return result;
    }

    public Map<ImKey, Double> getExpectedLossByImKey() {
        Map<ImKey, Double> result = new LinkedHashMap<>();
        for (ImKey k : getImKeys()) {
            result.put(
                    k,
                    getAssetsByImKey(k).stream()
                            .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                            .sum()
            );
        }
        return result;
    }

    public Map<String, Double> getExpectedLossByAdditionalField(String field) {
        Map<String, Double> result = new LinkedHashMap<>();
        for (Map.Entry<String, List<ExpectedLossAsset>> e :
                getAdditionalFieldMap(field).entrySet()) {

            result.put(
                    e.getKey(),
                    e.getValue().stream()
                            .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                            .sum()
            );
        }
        return result;
    }

    // ---------------------------------------------------------------------
    // CSV EXPORT
    // ---------------------------------------------------------------------

    public void writeCSV(Path file) throws IOException {

        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(file))) {

            pw.print("AssetID,Latitude,Longitude,Vs30,Value,Model");

            for (String f : getAdditionalFieldNames()) {
                pw.print("," + f);
            }

            pw.println(",ExpectedLoss");

            for (ExpectedLossAsset a : view.getAssets()) {

                pw.printf(
                        "%d,%.6f,%.6f,%.0f,%.6f,\"%s\"",
                        a.getAssetID(),
                        a.getLatitude(),
                        a.getLongitude(),
                        a.getVs30(),
                        a.getValue(),
                        a.getModelName()
                );

                Map<String, String> extra = a.getAdditionalFields();

                for (String f : getAdditionalFieldNames()) {
                    pw.printf(",\"%s\"", extra.getOrDefault(f, ""));
                }

                pw.printf(",%.6e%n", a.getExpectedLoss());
            }
        }
    }

    public void writeAggregatedCSV(Path file) throws IOException {

        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(file))) {

            pw.println(
                    "GroupField,GroupValue,NumAssets,TotalValue,TotalExpectedLoss," +
                    "Avg,Min,Max,StdDev"
            );

            // ----- Portfolio-wide summary -----
            Map<String, Double> allStats = RiskConvolutionPortfolioView.summarize(getAssets(), ExpectedLossAsset::getExpectedLoss);

            double totalValue = getAssets().stream()
                    .mapToDouble(ExpectedLossAsset::getValue)
                    .sum();

            pw.printf(
                    "\"%s\",\"%s\",%d,%.6f,%.6e,%.6f,%.6f,%.6f,%.6f%n",
                    "ALL",
                    "ALL",
                    allStats.get("n").intValue(),
                    totalValue,
                    allStats.get("total"),
                    allStats.get("avg"),
                    allStats.get("min"),
                    allStats.get("max"),
                    allStats.get("stdDev")
            );

            // ----- Group summaries -----
            for (String field : getAdditionalFieldNames()) {

                for (Map.Entry<String, List<ExpectedLossAsset>> e :
                        getAdditionalFieldMap(field).entrySet()) {

                    List<ExpectedLossAsset> group = e.getValue();

                    Map<String, Double> stats = RiskConvolutionPortfolioView.summarize(group, ExpectedLossAsset::getExpectedLoss);

                    totalValue = group.stream()
                            .mapToDouble(ExpectedLossAsset::getValue)
                            .sum();

                    pw.printf(
                            "\"%s\",\"%s\",%d,%.6f,%.6e,%.6f,%.6f,%.6f,%.6f%n",
                            field,
                            e.getKey(),
                            stats.get("n").intValue(),
                            totalValue,
                            stats.get("total"),
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

    /** Prints a summary of portfolio losses. */
    public void printSummary() {
        System.out.println("----- ELossPortfolio Summary -----");
        System.out.printf("Total assets: %d%n", view.size());
        System.out.printf("Total expected loss: %.2e%n", getTotalExpectedLoss());
    }
}
