package scratch.anne.sra.domain.portfolio.portfolioWrappers;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.portfolio.PortfolioGetters;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * A view over a {@link RiskConvolutionPortfolio} -- i.e. a portfolio
 * that has already been run through risk convolution in
 * {@code FULL_HCURVE} mode -- projected onto a specific
 * {@link RiskConvolutionAsset} subtype {@code T}.
 *
 * <p>Held as a field (composition), not extended, by each reporting
 * class that needs this projection -- {@link ExpectedLossPortfolioAggregator}
 * ({@link scratch.anne.sra.domain.asset.vulnerability.ExpectedLossAsset}),
 * {@link ProbFailurePortfolioReporter} ({@link scratch.anne.sra.domain.asset.fragility.FailureProbabilityAsset}),
 * {@link PBRSurvivalPortfolioReporter} ({@link scratch.anne.sra.domain.asset.fragility.PBRSurvivalAsset}).
 * None of those classes <em>is</em> a generic asset-projection view;
 * they each report a specific risk metric and merely use this to avoid
 * re-implementing the same delegation plumbing.</p>
 *
 * <p>Unlike a view that re-derives its projected assets from
 * {@code base} on every call, this class projects every asset through
 * {@code cast} exactly ONCE, in its constructor, and every delegated
 * getter below operates on that already-materialized {@link #assets}
 * list. This matters specifically for a {@code cast} that constructs a
 * new object per call (e.g. {@code PBRSurvivalPortfolioReporter}'s
 * {@code a -> new PBRSurvivalAsset(...)}) -- re-deriving on every call
 * would silently hand back a different {@code PBRSurvivalAsset}
 * instance each time for the "same" asset, which is never desired.</p>
 *
 * @param <T> the specific {@link RiskConvolutionAsset} subtype this view projects onto
 */
final class RiskConvolutionPortfolioView<T extends RiskConvolutionAsset> implements PortfolioGetters<T> {

    private final RiskConvolutionPortfolio base;
    private final List<T> assets;

    /** siteKey -> this view's assets at that site, derived once from {@link #assets}. */
    private final Map<SiteKey, List<T>> siteMap;
    /** assetId -> this view's asset with that id, derived once from {@link #assets}. */
    private final Map<Integer, T> assetsById;
    /** this view's asset (the projection, e.g. a PBRSurvivalAsset) -> its site key. */
    private final Map<T, SiteKey> siteKeyByAsset;

    /**
     * @param base the underlying risk-convolution portfolio; must already
     *             have been run in {@code FULL_HCURVE} mode
     * @param cast projects one of {@code base}'s assets into this view's
     *             asset type {@code T} -- called exactly once per asset,
     *             here in the constructor. A plain cast for views that
     *             don't transform the asset (e.g. {@code a -> (ExpectedLossAsset) a}),
     *             or a wrapping constructor for views that do (e.g.
     *             {@code a -> new PBRSurvivalAsset((FailureProbabilityAsset) a)})
     */
    RiskConvolutionPortfolioView(RiskConvolutionPortfolio base, Function<RiskConvolutionAsset, T> cast) {
        this.base = Objects.requireNonNull(base);
        Objects.requireNonNull(cast);

        base.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);

        // pair each ORIGINAL base asset with its projection, so site-key
        // lookups below use the original's identity -- not the
        // projection's, which for a wrapping cast (e.g. PBRSurvivalAsset)
        // is a freshly constructed object base.getSiteKey(...) has never
        // seen and would otherwise return null for it
        List<RiskConvolutionAsset> originals = base.getAssets();

        List<T> projected = new ArrayList<>(originals.size());
        Map<T, SiteKey> siteKeyByProjected = new LinkedHashMap<>();
        for (RiskConvolutionAsset original : originals) {
            T p = cast.apply(original);
            projected.add(p);
            siteKeyByProjected.put(p, base.getSiteKey(original));
        }
        this.assets = List.copyOf(projected);
        this.siteKeyByAsset = Map.copyOf(siteKeyByProjected);

        Map<Integer, T> idMap = new LinkedHashMap<>();
        for (T a : assets)
            idMap.put(a.getAssetID(), a);
        this.assetsById = Map.copyOf(idMap);

        this.siteMap = assets.stream()
                .collect(Collectors.groupingBy(
                        siteKeyByAsset::get,
                        LinkedHashMap::new,
                        Collectors.toUnmodifiableList()));
    }

    // ---------------------------------------------------------------------
    // BASIC ACCESS (DELEGATED, FROM ALREADY-MATERIALIZED assets)
    // ---------------------------------------------------------------------

    @Override
    public List<T> getAssets() {
        return assets;
    }

    @Override
    public Set<Integer> getAssetIDs() {
        return assetsById.keySet();
    }

    @Override
    public T getAssetByID(int assetID) {
        return assetsById.get(assetID);
    }

    @Override
    public AssetType getAssetType() {
        // base already confirmed homogeneity
        return assets.get(0).getAssetType();
    }

    public RiskMetricType getRiskMetricType() {
        // base already confirmed homogeneity
        return assets.get(0).getRiskMetricType();
    }

    @Override
    public int size() {
        return assets.size();
    }

    // ---------------------------------------------------------------------
    // SITE GROUPING (DELEGATED)
    // ---------------------------------------------------------------------

    @Override
    public Map<SiteKey, List<T>> getSiteMap() {
        return siteMap;
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return siteMap.keySet();
    }

    @Override
    public SiteKey getSiteKey(T asset) {
        return siteKeyByAsset.get(asset);
    }

    @Override
    public List<T> getAssetsBySite(SiteKey siteKey) {
        return siteMap.getOrDefault(siteKey, List.of());
    }

    // ---------------------------------------------------------------------
    // IM GROUPING (DELEGATED -- key sets only come from base; asset lists
    // are re-filtered from the already-materialized assets, not base)
    // ---------------------------------------------------------------------

    public Set<ImKey> getImKeys() {
        return base.getImKeys();
    }

    public List<T> getAssetsByImKey(ImKey key) {
        return assets.stream()
                .filter(a -> key.equals(a.getImKey()))
                .toList();
    }

    public List<T> getAssetsBySiteAndImKey(SiteKey site, ImKey imKey) {
        return getAssetsBySite(site).stream()
                .filter(a -> imKey.equals(a.getImKey()))
                .toList();
    }

    public Set<ImKey> getImKeysBySite(SiteKey site) {
        return base.getImKeysBySite(site);
    }

    // ---------------------------------------------------------------------
    // ADDITIONAL FIELDS (DELEGATED, FROM ALREADY-MATERIALIZED assets)
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
    public List<T> getAssetsByAdditionalField(String field, String value) {
        return assets.stream()
                .filter(a -> value.equals(a.getAdditionalFields().get(field)))
                .toList();
    }

    @Override
    public Map<String, List<T>> getAdditionalFieldMap(String field) {
        return assets.stream()
                .filter(a -> a.getAdditionalFields().containsKey(field))
                .collect(Collectors.groupingBy(
                        a -> a.getAdditionalFields().get(field),
                        LinkedHashMap::new,
                        Collectors.toUnmodifiableList()));
    }

    // ---------------------------------------------------------------------
    // RESPONSE MODELS
    // ---------------------------------------------------------------------

    @Override
    public Set<String> getResponseModelNames() {
        return base.getResponseModelNames();
    }

    // ---------------------------------------------------------------------
    // METADATA
    // ---------------------------------------------------------------------

    @Override
    public Metadata getMetadata() {
        return base.getMetadata();
    }

    // ---------------------------------------------------------------------
    // SHARED AGGREGATION
    // ---------------------------------------------------------------------

    /**
     * Summarizes {@code assets} by one metric: total, mean, min, max,
     * sample standard deviation (population formula, dividing by n),
     * and count. Returns an empty map for an empty list.
     *
     * @param metric extracts the value to summarize from each asset
     *               (e.g. {@code ExpectedLossAsset::getExpectedLoss})
     */
    static <T> Map<String, Double> summarize(List<T> assets, ToDoubleFunction<T> metric) {

        Map<String, Double> out = new LinkedHashMap<>();

        if (assets.isEmpty())
            return out;

        double total = assets.stream().mapToDouble(metric).sum();
        double avg = total / assets.size();
        double min = assets.stream().mapToDouble(metric).min().orElse(0);
        double max = assets.stream().mapToDouble(metric).max().orElse(0);
        double var = assets.stream()
                .mapToDouble(a -> Math.pow(metric.applyAsDouble(a) - avg, 2))
                .sum() / assets.size();

        out.put("total", total);
        out.put("avg", avg);
        out.put("min", min);
        out.put("max", max);
        out.put("stdDev", Math.sqrt(var));
        out.put("n", (double) assets.size());

        return out;
    }
}
