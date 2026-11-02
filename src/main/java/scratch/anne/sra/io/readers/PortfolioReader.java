package scratch.anne.sra.io.readers;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

import scratch.anne.sra.domain.asset.AbstractAsset;
import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.portfolio.*;
import scratch.anne.sra.util.InputParsingUtil;
import scratch.anne.sra.util.Metadata;

/**
 * Generic portfolio CSV reader for vulnerability or fragility assets.
 *
 * <p>This reader works for any asset type extending {@link AbstractAsset} that
 * provides a static {@code fromCSV(String[], Map<String,Integer>)} method.
 *
 * <p>The required columns differ depending on the asset type:
 * <ul>
 *   <li>VulnerabilityAsset: assetID, lat, lon, vs30, value, vulnModel</li>
 *   <li>FragilityAsset: assetID, lat, lon, vs30, fragilityModel</li>
 * </ul>
 *
 * <p>Any additional columns are collected as additional fields and preserved in the portfolio.
 */
public class PortfolioReader {

    /**
     * Reads a CSV into a typed Portfolio of assets with default metadata.
     *
     * @param <T>        the type of asset (VulnerabilityAsset or FragilityAsset)
     * @param file       path to the CSV file
     * @param assetClass the asset class to read
     * @return Portfolio containing all assets of type T
     * @throws IOException              if reading the file fails
     * @throws IllegalArgumentException if required columns are missing
     */
	public static <T extends AbstractAsset> Portfolio<T> readCSV(Path file, Class<T> assetClass)
	        throws IOException {
	    // default metadata constructed from CSV path
	    Metadata metadata = new Metadata.Builder()
	            .set("portfolioSource", file.toString())
	            .set("description", "Portfolio loaded from CSV")
	            .set("creationInfo", "Generated " + java.time.Instant.now())
	            .build();

        return readCSV(file, assetClass, metadata);
    }

    /**
     * Reads a CSV into a typed Portfolio of assets with provided metadata.
     *
     * @param <T>        the type of asset (VulnerabilityAsset or FragilityAsset)
     * @param file       path to the CSV file
     * @param assetClass the asset class to read
     * @param metadata   PortfolioMetadata object to attach
     * @return Portfolio containing all assets of type T
     * @throws IOException              if reading the file fails
     * @throws IllegalArgumentException if required columns are missing
     */
    public static <T extends AbstractAsset> Portfolio<T> readCSV(
            Path file, Class<T> assetClass, Metadata metadata)
            throws IOException {

        List<String> lines = Files.readAllLines(file);
        if (lines.isEmpty()) throw new IOException("CSV file is empty: " + file);

        // Remove BOM if present
        String firstLine = lines.get(0);
        if (firstLine.startsWith("\uFEFF")) firstLine = firstLine.substring(1);

        // Parse header
        String[] header = firstLine.split(",", -1);
        Map<String,Integer> colIndex = new LinkedHashMap<>();
        for (int i = 0; i < header.length; i++) colIndex.put(header[i].trim(), i);

        // Determine required columns based on asset type
        List<String> requiredCols;
        if (VulnerabilityAsset.class.isAssignableFrom(assetClass)) {
            requiredCols = Arrays.asList("assetID", "lat", "lon", "vs30", "value", "vulnModel");
        } else if (FragilityAsset.class.isAssignableFrom(assetClass)) {
            requiredCols = Arrays.asList("assetID", "lat", "lon", "vs30", "fragilityModel");
        } else {
            throw new IllegalArgumentException("Unsupported asset class: " + assetClass);
        }

        // Verify required columns exist
        for (String col : requiredCols) {
            if (!colIndex.containsKey(col)) {
                throw new IllegalArgumentException("Missing required column: " + col);
            }
        }

        // Identify additional fields
        List<String> additionalFieldNames = new ArrayList<>();
        for (String col : header) if (!requiredCols.contains(col)) additionalFieldNames.add(col);

        // Read assets
        List<T> assets = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.isEmpty()) continue;

            String[] row = Arrays.stream(line.split(",", -1))
                    .map(InputParsingUtil::cleanCsvString)
                    .toArray(String[]::new);

            try {
                @SuppressWarnings("unchecked")
                T asset = (T) assetClass
                        .getMethod("fromCSV", String[].class, Map.class)
                        .invoke(null, row, colIndex);
                assets.add(asset);
            } catch (Exception e) {
                throw new RuntimeException("Failed to create asset at line " + (i + 1), e);
            }
        }

        // Construct portfolio with metadata
        return new Portfolio<>(assets, additionalFieldNames, metadata);
    }

    // ------------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------------

    public static Map<String, String> extractAdditionalFields(
            String[] row,
            Map<String, Integer> colIndex,
            Set<String> requiredFields) {

        Map<String, String> extra = new LinkedHashMap<>();

        for (var entry : colIndex.entrySet()) {
            String col = entry.getKey();
            if (!requiredFields.contains(col)) {
                int idx = entry.getValue();
                String value = (idx < row.length) ? InputParsingUtil.cleanCsvString(row[idx]) : "";
                extra.put(col, value);
            }
        }
        return extra;
    }
}