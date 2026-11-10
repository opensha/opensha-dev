package scratch.anne.sra.examples.io;

import java.nio.file.Path;

import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.portfolio.*;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.util.AssetKeys;

/**
 * Demonstrates reading a portfolio CSV into a typed {@link Portfolio}
 * for both asset kinds this codebase supports -- {@link VulnerabilityAsset}
 * and {@link FragilityAsset} -- via {@link PortfolioReader#readCSV}, then
 * inspecting the result: metadata, size, additional fields, site
 * groupings, and a sample of individual assets.
 *
 */
public class PortfolioReaderExample {

    public static void main(String[] args) {
        try {
            // --- Example CSV files ---
        	Path resourceFolder = Path.of("C:/Users/ahulsey/git/opensha-dev/src/main/resources/scratch/anne/sra/");
            Path vulnPortfolioCSV = resourceFolder.resolve("P366_Porter-vulns/p366-portfolio_Porter-vuln_SHORT.csv");
            Path fragPortfolioCSV = resourceFolder.resolve("mini/p366-portfolio_placeholder-fragility_SHORT.csv");
            

            // --- Read a portfolio of VulnerabilityAssets ---
            Portfolio<VulnerabilityAsset> vulnPortfolio =
                    PortfolioReader.readCSV(vulnPortfolioCSV, VulnerabilityAsset.class);

            System.out.println("\nVuln Portfolio Loaded:");
            System.out.println("Metadata: " + vulnPortfolio.getMetadata());
            System.out.println("Total assets: " + vulnPortfolio.size());
            System.out.println("Additional fields: " + vulnPortfolio.getAdditionalFieldNames());

            // --- Summary of site groupings ---
            System.out.println("\nSite grouping summary:");
            for (AssetKeys.SiteKey key : vulnPortfolio.getSiteKeys()) {
                int count = vulnPortfolio.getAssetsBySite(key).size();
                System.out.printf("Site %s → %d assets%n", key, count);
            }

            // --- Summary of additional field groupings ---
            System.out.println("\nAdditional field groupings:");
            for (String field : vulnPortfolio.getAdditionalFieldNames()) {
                System.out.printf("Field '%s' values:%n", field);
                for (String value : vulnPortfolio.getAdditionalFieldValues(field)) {
                    int count = vulnPortfolio.getAssetsByAdditionalField(field, value).size();
                    System.out.printf("  Value '%s' → %d assets%n", value, count);
                }
            }

            // --- Optional: print first few assets for inspection ---
            System.out.println("\nSample assets:");
            vulnPortfolio.getAssets().stream().limit(5).forEach(asset -> {
                System.out.printf("AssetID=%s, Lat=%.3f, Lon=%.3f, Vs30=%.1f, Fragility=%s%n",
                        asset.getAssetID(),
                        asset.getLatitude(),
                        asset.getLongitude(),
                        asset.getVs30(),
                        asset.getModelName());
            });

         // --- Read a portfolio of FragilityAssets ---
            Portfolio<FragilityAsset> fragPortfolio =
                    PortfolioReader.readCSV(fragPortfolioCSV, FragilityAsset.class);

            System.out.println("\nFragility Portfolio Loaded:");
            System.out.println("Metadata: " + fragPortfolio.getMetadata());
            System.out.println("Total assets: " + fragPortfolio.size());
            System.out.println("Additional fields: " + fragPortfolio.getAdditionalFieldNames());

            // --- Summary of site groupings ---
            System.out.println("\nSite grouping summary:");
            for (AssetKeys.SiteKey key : fragPortfolio.getSiteKeys()) {
                int count = fragPortfolio.getAssetsBySite(key).size();
                System.out.printf("Site %s → %d assets%n", key, count);
            }

            // --- Summary of additional field groupings ---
            System.out.println("\nAdditional field groupings:");
            for (String field : fragPortfolio.getAdditionalFieldNames()) {
                System.out.printf("Field '%s' values:%n", field);
                for (String value : fragPortfolio.getAdditionalFieldValues(field)) {
                    int count = fragPortfolio.getAssetsByAdditionalField(field, value).size();
                    System.out.printf("  Value '%s' → %d assets%n", value, count);
                }
            }

            // --- Optional: print first few assets for inspection ---
            System.out.println("\nSample assets:");
            fragPortfolio.getAssets().stream().limit(5).forEach(asset -> {
                System.out.printf("AssetID=%s, Lat=%.3f, Lon=%.3f, Vs30=%.1f, Fragility=%s%n",
                        asset.getAssetID(),
                        asset.getLatitude(),
                        asset.getLongitude(),
                        asset.getVs30(),
                        asset.getModelName());
            });

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
