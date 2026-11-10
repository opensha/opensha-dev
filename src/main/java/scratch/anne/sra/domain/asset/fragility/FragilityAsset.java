package scratch.anne.sra.domain.asset.fragility;

import java.util.*;

import scratch.anne.sra.domain.asset.AbstractAsset;

/**
 * Asset pointing to a fragility for the response model.
 */
public class FragilityAsset extends AbstractAsset {

    /**
     * Full constructor for FragilityAsset.
     *
     * @param assetID     unique identifier
     * @param lat         latitude
     * @param lon         longitude
     * @param vs30        vs30 at the site
     * @param responseModel   name of response model
     * @param extraFields any extra grouping fields
     */
    public FragilityAsset(int assetID, double lat, double lon, double vs30,
                              String responseModel, Map<String, String> extraFields) {
        super(assetID, lat, lon, vs30, responseModel, extraFields, AssetType.FRAGILITY);
    }
   
    
    /**
     * Static factory method to construct a FragilityAsset from a CSV row.
     *
     * @param row      CSV row values
     * @param colIndex mapping from column name to index
     * @return new FragilityAsset instance
     * @throws IllegalArgumentException if assetID is not a valid integer
     */
    public static FragilityAsset fromCSV(String[] row, Map<String,Integer> colIndex) {

        String rawId = row[colIndex.get("assetID")];
        int id;
        try {
            id = Integer.parseInt(rawId.trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                    "assetID must be an integer, found: \"" + rawId + "\"", e);
        }
        double lat = Double.parseDouble(row[colIndex.get("lat")]);
        double lon = Double.parseDouble(row[colIndex.get("lon")]);
        double vs30 = Double.parseDouble(row[colIndex.get("vs30")]);
        String responseModel = row[colIndex.get("fragilityModel")];

        // Required fields for vulnerability assets
        Set<String> required = Set.of("assetID", "lat", "lon", "vs30", "fragilityModel");

        Map<String,String> extra =
                scratch.anne.sra.io.readers.PortfolioReader.extractAdditionalFields(
                        row, colIndex, required
                );

        return new FragilityAsset(id, lat, lon, vs30, responseModel, extra);
    }

    @Override
    public String toString() {
        return String.format("FragilityAsset[id=%d, model=%s]", assetID, responseModel);
    }
}
