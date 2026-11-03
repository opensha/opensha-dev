package scratch.anne.sra.domain.asset;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Base class for assets in a portfolio.
 * <p>
 * Immutable. Each asset has an identifier, location (lat/lon), vs30, responseModel
 * and a map of additional info fields from the CSV.
 */
public abstract class AbstractAsset {

    protected final int assetID;
    protected final AssetType assetType;
    protected final double lat;
    protected final double lon;
    protected final double vs30;
    protected final String responseModel;
    protected final Map<String, String> additionalFields;

    /**
     * @param assetID     Asset identifier
     * @param lat         Latitude
     * @param lon         Longitude
     * @param vs30        Site vs30 value
     * @param responseModel  Identifier for a vulnerability or fragility
     * @param additionalFields Map of additional CSV fields (nullable)
     */
    protected AbstractAsset(int assetID, double lat, double lon, double vs30, String responseModel, Map<String, String> additionalFields, AssetType assetType) {
        if (responseModel == null) throw new IllegalArgumentException("Response model must be non-null");
        this.assetID = assetID;
        this.lat = lat;
        this.lon = lon;
        this.vs30 = vs30;
        this.responseModel = responseModel;
        this.additionalFields = additionalFields != null ? Collections.unmodifiableMap(new LinkedHashMap<>(additionalFields)) : Collections.emptyMap();
        
        this.assetType = Objects.requireNonNull(assetType);
    }

    public int getAssetID() { return assetID; }
    public double getLatitude() { return lat; }
    public double getLongitude() { return lon; }
    public double getVs30() { return vs30; }
    public String getModelName() { return responseModel; }
    public Map<String, String> getAdditionalFields() { return additionalFields; }
    
    public String getAdditionalFieldValue(String name) {
        if (!additionalFields.containsKey(name)) {
            throw new IllegalArgumentException("Missing information: " + name);
        }
        return additionalFields.get(name);
    }
    
    public AssetType getAssetType() {
    	return assetType;
    }

    @Override
    public String toString() {
        return String.format("Asset[id=%d, lat=%.5f, lon=%.5f, vs30=%.1f, responseModel=%s]", assetID, lat, lon, vs30, responseModel);
    }
    
    public enum AssetType {
    	VULNERABILITY,
        FRAGILITY,
    }
}
