package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.vulnerability.ExpectedLossAsset;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Aggregation + reporting view over a RiskConvolutionPortfolio whose
 * assets are already materialized ExpectedLossAsset instances.
 *
 * <p>
 * This class does NOT transform or copy asset state. It only provides
 * aggregated views and export utilities.
 * </p>
 */
public final class ExpectedLossPortfolioAggregator
        implements Iterable<ExpectedLossAsset>,
                   PortfolioGetters<ExpectedLossAsset> {

    // ---------------------------------------------------------------------
    // BASE PORTFOLIO
    // ---------------------------------------------------------------------

    private final RiskConvolutionPortfolio base;

    private final List<ExpectedLossAsset> assets;

    // ---------------------------------------------------------------------
    // CONSTRUCTOR
    // ---------------------------------------------------------------------

    public ExpectedLossPortfolioAggregator(RiskConvolutionPortfolio base) {

        this.base = Objects.requireNonNull(base);
        
        base.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
        
        this.assets = List.copyOf(
                base.getAssets().stream()
                        .map(a -> (ExpectedLossAsset) a)
                        .toList()
        );

    }

    // ---------------------------------------------------------------------
    // BASIC ACCESS (DELEGATED)
    // ---------------------------------------------------------------------

    @Override
    public List<ExpectedLossAsset> getAssets() {
        return assets;
    }

    @Override
    public Set<String> getAssetIDs() {
        return base.getAssetIDs();
    }

    @Override
    public ExpectedLossAsset getAssetByID(String assetID) {
        return (ExpectedLossAsset) base.getAssetByID(assetID);
    }
    
    @Override
    public AssetType getAssetType() {
    	// risk convolution portofolio already confirm homogeneity
    	return assets.get(0).getAssetType();
    }
    
    public RiskMetricType getRiskMetricType() {
    	// risk convolution portofolio already confirm homogeneity
    	return assets.get(0).getRiskMetricType();
    }

    @Override
    public int size() {
        return base.size();
    }

    // ---------------------------------------------------------------------
    // SITE GROUPING (DELEGATED)
    // ---------------------------------------------------------------------

    @Override
    public Map<SiteKey, List<ExpectedLossAsset>> getSiteMap() {
        return base.getSiteMap().entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> (ExpectedLossAsset) a)
                                .toList()
                ));
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return base.getSiteKeys();
    }
    
    @Override
    public SiteKey getSiteKey(ExpectedLossAsset asset) {
        return base.getSiteKey(asset);
    }

    @Override
    public List<ExpectedLossAsset> getAssetsBySite(SiteKey siteKey) {
        return base.getAssetsBySite(siteKey).stream()
                .map(a -> (ExpectedLossAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    // ---------------------------------------------------------------------
    // IM GROUPING (DELEGATED)
    // ---------------------------------------------------------------------

    public Set<ImKey> getImKeys() {
        return base.getImKeys();
    }

    public List<ExpectedLossAsset> getAssetsByImKey(ImKey key) {
        return base.getAssetsByImKey(key).stream()
                .map(a -> (ExpectedLossAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    public List<ExpectedLossAsset> getAssetsBySiteAndImKey(
            SiteKey site,
            ImKey imKey
    ) {
        return base.getAssetsBySiteAndImKey(site, imKey).stream()
                .map(a -> (ExpectedLossAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    public Set<ImKey> getImKeysBySite(SiteKey site) {
        return base.getImKeysBySite(site);
    }

    // ---------------------------------------------------------------------
    // ADDITIONAL FIELDS (DELEGATED)
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
    public List<ExpectedLossAsset> getAssetsByAdditionalField(
            String field,
            String value
    ) {
        return base.getAssetsByAdditionalField(field, value).stream()
                .map(a -> (ExpectedLossAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    @Override
    public Map<String, List<ExpectedLossAsset>> getAdditionalFieldMap(String field) {
        return base.getAdditionalFieldMap(field).entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> (ExpectedLossAsset) a)
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
    public Iterator<ExpectedLossAsset> iterator() {
        return assets.iterator();
    }

    // ---------------------------------------------------------------------
    // AGGREGATIONS
    // ---------------------------------------------------------------------

    public double getTotalExpectedLoss() {
        return assets.stream()
                .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                .sum();
    }

    public double getTotalAssetValue() {
        return assets.stream()
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

            for (ExpectedLossAsset a : assets) {

                pw.printf(
                        "\"%s\",%.6f,%.6f,%.0f,%.6f,\"%s\"",
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
            Map<String, Double> allStats = summarize(getAssets());

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

                    Map<String, Double> stats =
                            summarize(group);

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
        System.out.printf("Total assets: %d%n", assets.size());
        System.out.printf("Total expected loss: %.2e%n", getTotalExpectedLoss());
    }

    private static Map<String, Double> summarize(List<ExpectedLossAsset> assets) {

        Map<String, Double> out = new LinkedHashMap<>();

        if (assets.isEmpty()) return out;

        double total = assets.stream()
                .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                .sum();

        double avg = total / assets.size();

        double min = assets.stream()
                .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                .min().orElse(0);

        double max = assets.stream()
                .mapToDouble(ExpectedLossAsset::getExpectedLoss)
                .max().orElse(0);

        double var = assets.stream()
                .mapToDouble(a -> Math.pow(a.getExpectedLoss() - avg, 2))
                .sum() / assets.size();

        out.put("total", total);
        out.put("avg", avg);
        out.put("min", min);
        out.put("max", max);
        out.put("stdDev", Math.sqrt(var));
        out.put("n", (double) assets.size());

        return out;
    }

	@Override
	public Set<String> getResponseModelNames() {
		return base.getResponseModelNames();
	}
}
