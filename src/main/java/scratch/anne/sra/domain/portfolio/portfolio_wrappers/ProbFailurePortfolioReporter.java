package scratch.anne.sra.domain.portfolio.portfolio_wrappers;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.fragility.FailureProbabilityAsset;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Reporting reporting view over a RiskConvolutionPortfolio whose
 * assets are already materialized FailureProbabilityAsset instances.
 *
 * <p>
 * This class does NOT transform or copy asset state. It only provides
 * export utilities.
 * </p>
 */
public final class ProbFailurePortfolioReporter
        implements Iterable<FailureProbabilityAsset>,
                   PortfolioGetters<FailureProbabilityAsset> {

    // ---------------------------------------------------------------------
    // BASE PORTFOLIO
    // ---------------------------------------------------------------------

    private final RiskConvolutionPortfolio base;

    private final List<FailureProbabilityAsset> assets;

    // ---------------------------------------------------------------------
    // CONSTRUCTOR
    // ---------------------------------------------------------------------

    public ProbFailurePortfolioReporter(RiskConvolutionPortfolio base) {

        this.base = Objects.requireNonNull(base);
        
        base.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
        
        this.assets = List.copyOf(
                base.getAssets().stream()
                        .map(a -> (FailureProbabilityAsset) a)
                        .toList()
        );

    }

    // ---------------------------------------------------------------------
    // BASIC ACCESS (DELEGATED)
    // ---------------------------------------------------------------------

    @Override
    public List<FailureProbabilityAsset> getAssets() {
        return assets;
    }

    @Override
    public Set<Integer> getAssetIDs() {
        return base.getAssetIDs();
    }

    @Override
    public FailureProbabilityAsset getAssetByID(int assetID) {
        return (FailureProbabilityAsset) base.getAssetByID(assetID);
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
    public int size() {
        return base.size();
    }

    // ---------------------------------------------------------------------
    // SITE GROUPING (DELEGATED)
    // ---------------------------------------------------------------------

    @Override
    public Map<SiteKey, List<FailureProbabilityAsset>> getSiteMap() {
        return base.getSiteMap().entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> (FailureProbabilityAsset) a)
                                .toList()
                ));
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return base.getSiteKeys();
    }
    
    @Override
    public SiteKey getSiteKey(FailureProbabilityAsset asset) {
        return base.getSiteKey(asset);
    }

    @Override
    public List<FailureProbabilityAsset> getAssetsBySite(SiteKey siteKey) {
        return base.getAssetsBySite(siteKey).stream()
                .map(a -> (FailureProbabilityAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    // ---------------------------------------------------------------------
    // IM GROUPING (DELEGATED)
    // ---------------------------------------------------------------------

    public Set<ImKey> getImKeys() {
        return base.getImKeys();
    }

    public List<FailureProbabilityAsset> getAssetsByImKey(ImKey key) {
        return base.getAssetsByImKey(key).stream()
                .map(a -> (FailureProbabilityAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    public List<FailureProbabilityAsset> getAssetsBySiteAndImKey(
            SiteKey site,
            ImKey imKey
    ) {
        return base.getAssetsBySiteAndImKey(site, imKey).stream()
                .map(a -> (FailureProbabilityAsset) a)
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
    public List<FailureProbabilityAsset> getAssetsByAdditionalField(
            String field,
            String value
    ) {
        return base.getAssetsByAdditionalField(field, value).stream()
                .map(a -> (FailureProbabilityAsset) a)
                .collect(java.util.stream.Collectors.toList());
    }

    @Override
    public Map<String, List<FailureProbabilityAsset>> getAdditionalFieldMap(String field) {
        return base.getAdditionalFieldMap(field).entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> e.getValue().stream()
                                .map(a -> (FailureProbabilityAsset) a)
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
    public Iterator<FailureProbabilityAsset> iterator() {
        return assets.iterator();
    }


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

            for (FailureProbabilityAsset a : assets) {

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
                            summarize(group);

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
        System.out.printf("Total assets: %d%n", assets.size());
    }

    private static Map<String, Double> summarize(List<FailureProbabilityAsset> assets) {

        Map<String, Double> out = new LinkedHashMap<>();

        if (assets.isEmpty()) return out;

        double total = assets.stream()
                .mapToDouble(FailureProbabilityAsset::getProbabilityOfFailure)
                .sum();

        double avg = total / assets.size();

        double min = assets.stream()
                .mapToDouble(FailureProbabilityAsset::getProbabilityOfFailure)
                .min().orElse(0);

        double max = assets.stream()
                .mapToDouble(FailureProbabilityAsset::getProbabilityOfFailure)
                .max().orElse(0);

        double var = assets.stream()
                .mapToDouble(a -> Math.pow(a.getProbabilityOfFailure() - avg, 2))
                .sum() / assets.size();

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

