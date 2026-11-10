package scratch.anne.sra.examples.risk_convolution;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;

import scratch.anne.sra.domain.asset.AbstractAsset;
import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.engine.per_rupture.RuptureHazardCurveCollection;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.io.writers.RuptureProbabilityCurvesExporter;
import scratch.anne.sra.util.InputParsingUtil;
import scratch.anne.sra.util.NumericUtil;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.ImKey.ImDomain;


/**
 * For up to {@code nAssets} assets in a run folder's portfolio, builds
 * the shared ERF/GMM once and computes each asset's per-rupture
 * conditional exceedance-probability curves via
 * {@link RuptureHazardCurveCollection#fromERF}, writing one CSV
 * per asset via {@link RuptureProbabilityCurvesExporter}.
 *
 * <p><b>Choice points -- edit before running:</b></p>
 * <ul>
 *   <li>{@code baseFolder}/{@code inputFolder} select which run
 *       folder's {@code input.csv} (and portfolio it points to) drives
 *       the run.</li>
 *   <li>{@code nAssets} caps how many portfolio assets are processed.</li>
 *   <li>{@code imt}/{@code imMin}/{@code imMax}/{@code logImStep}
 *       define the IM grid each asset's curves are evaluated at.</li>
 * </ul>
 */
public class PerRuptureHazardCurvesExample {

    public static Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\my_scratch\\conversion to Java project\\BERM_test-outputs\\_vb\\csv_inputs\\tests\\disagg");
    public static Path inputFolder = Path.of("");
    
    static int nAssets = 1;
    
    static String imt = "SA(1.0)";
    static double imMin = 0.001;
    static double imMax = 10.0;
    static double logImStep = 0.02;

    public static void main(String[] args) throws Exception {

        Path runFolder = baseFolder.resolve(inputFolder);

        if (!Files.isDirectory(runFolder))
            throw new IllegalArgumentException("Run folder does not exist: " + runFolder);

        Path configCSV = runFolder.resolve("input.csv");

        if (!Files.exists(configCSV))
            throw new IllegalArgumentException("Missing input.csv in run folder: " + runFolder);

        System.out.println("Running rupture probability curves from folder:");
        System.out.println(runFolder);

        Map<String, String> config = ParseRiskRunParametersCSV.loadConfig(configCSV);

        // ------------------- 1. Resolve input files -------------------
        Path portfolioCSV = Path.of(config.get("portfolio_csv"));
        InputParsingUtil.verifyFileExists(portfolioCSV, "Portfolio CSV");

        // ------------------- 2. output paths -------------------
        String fileTag = config.getOrDefault("file_tag", runFolder.getFileName().toString());

        // ------------------- 3. hazard parameters -------------------
        
        // erf_config.csv is optional -- returns ErfConfig.EMPTY if absent,
        // meaning the ERF self-constructs from input.csv's erf_class with
        // no extra parameters, exactly as before this file existed.
        ParseRiskRunParametersCSV.ErfConfig erfConfig =
                ParseRiskRunParametersCSV.loadErfConfig(runFolder, config);
        
        HazardParameters hazardParameters =
                ParseRiskRunParametersCSV.buildHazardParameters(config, erfConfig);

        // ------------------- 4. IM grid -------------------
        double[] imGrid = NumericUtil.distributeLogSpacedValues(imMin, imMax, logImStep);
        ImKey imKey = new ImKey(imt, imGrid, ImDomain.LINEAR_IM);

        // ------------------- 5. Load portfolio -------------------
        System.out.println("Loading portfolio...");
        Portfolio<VulnerabilityAsset> portfolio =
                PortfolioReader.readCSV(portfolioCSV, VulnerabilityAsset.class);

        // ------------------- 6. build shared ERF once -------------------
        AbstractERF erf;
        try {
            erf = (AbstractERF)
                    Class.forName(hazardParameters.getErfName())
                            .getDeclaredConstructor()
                            .newInstance();
        } catch (Exception e) {
            throw new IllegalArgumentException(
                    "Failed to instantiate ERF class: " + hazardParameters.getErfName(), e);
        }
        erf.getTimeSpan().setDuration(hazardParameters.getErfDuration());
        erf.updateForecast();

        AttenRelRef gmmRef = AttenRelRef.valueOf(hazardParameters.getGmmName().toUpperCase());

        // ------------------- 7. loop over assets -------------------
        int nProcessed = 0;
        for (AbstractAsset asset : portfolio.getAssets()) {
        	
            if (nProcessed >= nAssets)
                break;

            System.out.println("Computing rupture exceedance curves for asset "
                    + asset.getAssetID() + " at "
                    + asset.getLatitude() + ", " + asset.getLongitude()
                    + " (vs30=" + asset.getVs30() + ")...");

            ScalarIMR gmm = gmmRef.get();
            gmm.setParamDefaults();

            RuptureHazardCurveCollection collection =
                    RuptureHazardCurveCollection.fromERF(
                            erf, gmm, hazardParameters, asset, imKey
                    );
            collection.printSummary();

            Path outputCSV = runFolder.resolve(
                    fileTag + "_site-" + nProcessed + "_rupture-probability-curves.csv");

            new RuptureProbabilityCurvesExporter().export(collection, outputCSV);

            System.out.println("Wrote " + outputCSV.toAbsolutePath());
            
            nProcessed++;
        }
    }

}