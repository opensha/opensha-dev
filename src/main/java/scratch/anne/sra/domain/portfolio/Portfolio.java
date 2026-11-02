package scratch.anne.sra.domain.portfolio;

import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.AbstractAsset;
import scratch.anne.sra.domain.asset.AbstractAsset.AssetType;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.PortfolioGroupingUtils;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Immutable canonical portfolio.
 *
 * <p>
 * Holds assets and deterministic structural groupings derived from
 * asset data.
 * </p>
 *
 * <p>
 * All ordering is defined by original asset list iteration order and frozen
 * at construction time.
 * </p>
 */
public class Portfolio<T extends AbstractAsset> implements PortfolioGetters<T> {

    // ------------------------------------------------------------------------
    // Core data
    // ------------------------------------------------------------------------

    private final List<T> assets;
    private final Map<String, T> assetMap;
    
    private final AssetType assetType;

    /** Site grouping (lat, lon, vs30) */
    private final Map<SiteKey, List<T>> siteMap;
    private final Map<T, SiteKey> siteKeyByAsset;

    /** Additional field grouping: field → value → assets */
    private final Map<String, Map<String, List<T>>> additionalFieldGroups;

    /** Ordered list of additional field names (derived deterministically) */
    private final List<String> additionalFieldNames;

    // ------------------------------------------------------------------------
    // Metadata
    // ------------------------------------------------------------------------

    private final Metadata metadata;

    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------

    public Portfolio(List<T> assets,
                     List<String> additionalFieldNames,
                     Metadata metadata) {

        if (assets == null || assets.isEmpty()) {
            throw new IllegalArgumentException("Assets must be non-null and non-empty");
        }

        validateHomogeneous(assets);
        this.assetType = assets.get(0).getAssetType();

        // -------------------------
        // Core immutable storage
        // -------------------------
        this.assets = List.copyOf(assets);

        this.assetMap = assets.stream()
                .collect(Collectors.toUnmodifiableMap(
                        T::getAssetID,
                        a -> a
                ));

        this.metadata = (metadata == null)
                ? new Metadata.Builder().build()
                : metadata;

        // -------------------------
        // Site grouping (stable order)
        // -------------------------
        this.siteMap = PortfolioGroupingUtils.groupBy(
                assets,
                a -> new SiteKey(a.getLatitude(), a.getLongitude(), a.getVs30())
        );

        // -------------------------
        // Additional field grouping (deterministic inversion)
        // -------------------------
        this.additionalFieldGroups =
                PortfolioGroupingUtils.invertNestedFields(
                        assets,
                        AbstractAsset::getAdditionalFields
                );
        
        Map<T, SiteKey> tmp = new LinkedHashMap<>();
        for (Map.Entry<SiteKey, List<T>> e : siteMap.entrySet()) {
            SiteKey site = e.getKey();
            for (T asset : e.getValue()) {
                tmp.put((T) asset, site);
            }
        }
        this.siteKeyByAsset = Collections.unmodifiableMap(tmp);

        // -------------------------
        // Field name ordering (derived from actual data)
        // -------------------------
        this.additionalFieldNames =
                List.copyOf(additionalFieldGroups.keySet());
    }

    // ------------------------------------------------------------------------
    // Convenience constructors
    // ------------------------------------------------------------------------

    public Portfolio(List<T> assets, List<String> additionalFieldNames) {
        this(assets, additionalFieldNames, null);
    }

    public Portfolio(List<T> assets, Metadata metadata) {
        this(assets, null, metadata);
    }

    public Portfolio(List<T> assets) {
        this(assets, null, null);
    }

    // ------------------------------------------------------------------------
    // Validation
    // ------------------------------------------------------------------------

    private void validateHomogeneous(List<T> assets) {
        if (assets == null || assets.isEmpty()) {
            return;
        }

        AssetType type = assets.get(0).getAssetType();

        for (T asset : assets) {
            if (asset.getAssetType() != type) {
                throw new IllegalArgumentException(
                    "Portfolio must contain homogeneous asset types. Expected " 
                    + type + " but found " + asset.getAssetType()
                );
            }
        }
    }

    // ------------------------------------------------------------------------
    // Core getters
    // ------------------------------------------------------------------------

    @Override
    public List<T> getAssets() {
        return assets;
    }

    @Override
    public Set<String> getAssetIDs() {
        return assetMap.keySet();
    }
    
    @Override
    public AssetType getAssetType() {
    	return assetType;
    }
    

    @Override
    public T getAssetByID(String assetID) {
        T asset = assetMap.get(assetID);
        if (asset == null) {
            throw new NoSuchElementException("No asset with ID: " + assetID);
        }
        return asset;
    }

    @Override
    public Set<String> getResponseModelNames() {
        return assets.stream()
                .map(AbstractAsset::getModelName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    @Override
    public int size() {
        return assets.size();
    }

    @Override
    public Metadata getMetadata() {
        return metadata;
    }

    // ------------------------------------------------------------------------
    // Site grouping
    // ------------------------------------------------------------------------

    @Override
    public Map<SiteKey, List<T>> getSiteMap() {
        return siteMap;
    }

    @Override
    public Set<SiteKey> getSiteKeys() {
        return siteMap.keySet();
    }

    @Override
    public List<T> getAssetsBySite(SiteKey key) {
        return siteMap.getOrDefault(key, List.of());
    }
    
    @Override
    public SiteKey getSiteKey(T asset) {
        return siteKeyByAsset.get(asset);
    }

    // ------------------------------------------------------------------------
    // Additional field access
    // ------------------------------------------------------------------------

    @Override
    public List<String> getAdditionalFieldNames() {
        return additionalFieldNames;
    }

    @Override
    public List<String> getAdditionalFieldValues(String field) {
        Map<String, List<T>> valueMap = additionalFieldGroups.get(field);
        return valueMap == null
                ? List.of()
                : List.copyOf(valueMap.keySet());
    }

    @Override
    public List<T> getAssetsByAdditionalField(String field, String value) {
        Map<String, List<T>> valueMap = additionalFieldGroups.get(field);
        if (valueMap == null) return List.of();
        return valueMap.getOrDefault(value, List.of());
    }

    @Override
    public Map<String, List<T>> getAdditionalFieldMap(String fieldName) {
        Map<String, List<T>> map = additionalFieldGroups.get(fieldName);
        return map != null ? map : Map.of();
    }

    // ------------------------------------------------------------------------
    // Debug
    // ------------------------------------------------------------------------

    @Override
    public String toString() {
        return String.format(
                "Portfolio[nAssets=%d, description=%s, additionalFields=%s]",
                assets.size(),
                metadata.get("description"),
                additionalFieldNames
        );
    }
}