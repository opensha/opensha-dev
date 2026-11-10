package scratch.anne.sra.examples.riskConvolution;

import java.nio.file.*;
import java.util.*;

import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.HazardParametersMetadataBuilder;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.PortfolioMetadataBuilder;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.ExpectedLossPortfolioAggregator;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structuralResponse.ResponseModelLibrary;
import scratch.anne.sra.domain.structuralResponse.ResponseModelMetadataBuilder;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structuralResponse.vulnerabilities.SimpleImVulnLibraryPreparer;
import scratch.anne.sra.domain.structuralResponse.vulnerabilities.VulnerabilityModel;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.engine.portfolioWorkflow.PortfolioPerRuptureRiskConvolutionCalculator;
import scratch.anne.sra.engine.portfolioWorkflow.PortfolioRiskConvolutionCalculator;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.io.readers.VulnerabilityLibraryReader;
import scratch.anne.sra.io.writers.AssetRiskForRuptureWriter;
import scratch.anne.sra.util.ImValueTransformer;
import scratch.anne.sra.util.InputParsingUtil;

/**
 * Fully dynamic CSV-driven Expected Loss Portfolio runner with input
 * verification. Reads a run folder's {@code input.csv} (and optional
 * {@code erf_config.csv}), convolves the portfolio against hazard in
 * one or both convolution modes, and writes the corresponding outputs
 * (full per-asset, aggregated, and/or per-rupture CSVs; optionally a
 * hazard JSON and a thresholded per-asset-per-rupture Parquet export).
 *
 * <p><b>Choice points -- edit these before running:</b></p>
 * <ul>
 *   <li>{@code baseFolder}/{@code inputFolder} select which run
 *       folder's {@code input.csv} drives the run (several alternatives
 *       are commented out above the active one).</li>
 *   <li>{@code convolutionModes} selects which {@link ConvolutionMode}(s)
 *       to run -- {@code FULL_HCURVE}, {@code PER_RUPTURE}, or both
 *       (alternatives commented out above the active one).</li>
 *   <li>{@code writeHazard} toggles writing the hazard JSON
 *       ({@code FULL_HCURVE} mode only).</li>
 *   <li>{@code keepAssetLoss}/{@code fRelative}/{@code fAbsolute}
 *       control the optional thresholded per-asset-per-rupture Parquet
 *       writer ({@code PER_RUPTURE} mode only) -- see
 *       {@link AssetRiskForRuptureWriter#writeThresholdValue}.</li>
 * </ul>
 */
public class PortfolioExpectedLossExample {
		
//	public static Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\BERM\\results\\EAL\\rate\\riemann\\porter_vulns");
//	public static Path inputFolder = Path.of("PorterVulns_vs-360");
	
//	public static Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\my_scratch\\conversion to Java project\\BERM_test-outputs\\_vb\\csv_inputs");
//	public static Path inputFolder = Path.of("tests\\asset-per-rupture\\short_portfolio");
//	public static Path inputFolder = Path.of("tests\\asset-per-rupture\\sparse_portfolio");
	
	public static Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\BERM\\results\\EAL\\true_mean_solution");
	public static Path inputFolder = Path.of("RES1-W1-MC_Porter-vuln_vs-360");
	
//	public static List<ConvolutionMode> convolutionModes = List.of(ConvolutionMode.FULL_HCURVE,ConvolutionMode.PER_RUPTURE);
	public static List<ConvolutionMode> convolutionModes = List.of(ConvolutionMode.PER_RUPTURE);
//	public static List<ConvolutionMode> convolutionModes = List.of(ConvolutionMode.FULL_HCURVE);
	
	public static boolean writeHazard = true;
	
	public static boolean keepAssetLoss = false;
	public static double fRelative = 0d;
	public static double fAbsolute = 0d;
	

    public static void main(String[] args) throws Exception {

        Path runFolder = baseFolder.resolve(inputFolder);

        if (!Files.isDirectory(runFolder))
            throw new IllegalArgumentException("Run folder does not exist: " + runFolder);

        Path configCSV = runFolder.resolve("input.csv");

        if (!Files.exists(configCSV))
            throw new IllegalArgumentException("Missing input.csv in run folder: " + runFolder);

        System.out.println("Running portfolio from folder:");
        System.out.println(runFolder);
        
        Map<String, String> config =
                ParseRiskRunParametersCSV.loadConfig(configCSV);

        // erf_config.csv is optional -- returns ErfConfig.EMPTY if absent,
        // meaning the ERF self-constructs from input.csv's erf_class with
        // no extra parameters, exactly as before this file existed.
        ParseRiskRunParametersCSV.ErfConfig erfConfig =
                ParseRiskRunParametersCSV.loadErfConfig(runFolder, config);

        if (erfConfig.usesSolutionFile()) {
            System.out.println("erf_config.csv found with solutionFile set -- using solution-file "
                    + "route. input.csv's erf_class ('" + config.get("erf_class")
                    + "') will be IGNORED; erf_config.csv's erfClass ('" + erfConfig.erfClass()
                    + "') governs construction instead.");
        }

        // ------------------- 1. Resolve input files -------------------
        Path portfolioCSV = Path.of(config.get("portfolio_csv"));
        Path vulnLibraryJSON = Path.of(config.get("vuln_library_json"));

        InputParsingUtil.verifyFileExists(portfolioCSV, "Portfolio CSV");
        InputParsingUtil.verifyFileExists(vulnLibraryJSON, "Vulnerability JSON");

        // ------------------- 2. create output paths inside run folder -------------------
        String fileTag = config.getOrDefault("file_tag", runFolder.getFileName().toString());
        Path fullOutputCSV = runFolder.resolve(fileTag + "_eal-all-assets.csv");
        Path aggregatedOutputCSV = runFolder.resolve(fileTag + "_aggregated-eal.csv");
        Path perRuptureOutputCSV = runFolder.resolve(fileTag + "_eal-per-rup.csv");
        Path hazardJson = writeHazard
                ? runFolder.resolve(fileTag + "_hazard-list.json")
                : null;
        
        // ------------------- 3. prepare hazard and integration parameters ----------------
        HazardParameters hazardParameters =
                ParseRiskRunParametersCSV.buildHazardParameters(config, erfConfig);
        
        RiskConvolution.IntegrationMethod integrationMethod = ParseRiskRunParametersCSV.parseIntegrationMethod(config.get("integration_method"));
        double logImStep = InputParsingUtil.parseDoubleOrDefault(config.get("log_im_step"), Double.NaN);
        

        // ------------------- 4. Load portfolio -------------------
        System.out.println("Loading portfolio...");
        Portfolio<VulnerabilityAsset> basePortfolio =
                PortfolioReader.readCSV(portfolioCSV, VulnerabilityAsset.class);

        Set<String> portfolioVulns = basePortfolio.getResponseModelNames();

        // ------------------- 5. Load vulnerability library -------------------
        System.out.println("Loading full vulnerability library...");
        ResponseModelLibrary<VulnerabilityModel> vulnLib =
                VulnerabilityLibraryReader.readLibrary(vulnLibraryJSON);

        // ------------------- 6. Prepare Expected Vulnerability Library -------------------
        System.out.println("Preparing vulnerability library...");
        ImValueTransformer transformer = Double.isNaN(logImStep)
                ? new ImValueTransformer.NoImTransformation()
                : new ImValueTransformer.LogInterpTransformation(logImStep);
        SimpleImResponseLibrary expVulnLib =
                SimpleImVulnLibraryPreparer.prepare(vulnLib, portfolioVulns, transformer);

        // ------------------- 6b. asset writer, with full provenance -------------------
        Map<String, String> writerMetadata = new HashMap<>();
        writerMetadata.putAll(HazardParametersMetadataBuilder.build(hazardParameters));
        writerMetadata.putAll(PortfolioMetadataBuilder.build(basePortfolio));
        writerMetadata.putAll(ResponseModelMetadataBuilder.build(expVulnLib));

        AssetRiskForRuptureWriter assetWriter = keepAssetLoss
        		? new AssetRiskForRuptureWriter(runFolder, fileTag + "_asset-eal-per-rup", fRelative, fAbsolute, writerMetadata)
        		: null;

        for (ConvolutionMode convolutionMode : convolutionModes) {
	        if (convolutionMode == ConvolutionMode.FULL_HCURVE) {
	        	
	        	RiskConvolutionPortfolio riskConvolutionPortfolio = new RiskConvolutionPortfolio(basePortfolio, expVulnLib);
	        	
	        	PortfolioRiskConvolutionCalculator calculator =
	                    new PortfolioRiskConvolutionCalculator(
	                            riskConvolutionPortfolio,
	                            expVulnLib,
	                            hazardParameters,
	                            integrationMethod
	                    );
	            
	            riskConvolutionPortfolio = calculator.computeRisk();
	            riskConvolutionPortfolio.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
	            if (writeHazard) { riskConvolutionPortfolio.exportHazard(hazardJson); };
	            ExpectedLossPortfolioAggregator aggregator = new ExpectedLossPortfolioAggregator(riskConvolutionPortfolio);
	            aggregator.writeCSV(fullOutputCSV);
	            aggregator.writeAggregatedCSV(aggregatedOutputCSV);
	            
	            aggregator.printSummary();
	            System.out.println("Outputs written to run folder:");
	            System.out.println(fullOutputCSV);
	            System.out.println(aggregatedOutputCSV);
	            if (writeHazard) { System.out.println(hazardJson); };

	        }
	        else if (convolutionMode == ConvolutionMode.PER_RUPTURE) {
	        	
	        	RiskConvolutionPortfolio riskConvolutionPortfolio = new RiskConvolutionPortfolio(basePortfolio, expVulnLib);
	        	
	            PortfolioPerRuptureRiskConvolutionCalculator calculator =
	                    new PortfolioPerRuptureRiskConvolutionCalculator(
	                            riskConvolutionPortfolio,
	                            expVulnLib,
	                            hazardParameters,
	                            integrationMethod,
	                            assetWriter
	                    );
	            
	            riskConvolutionPortfolio = calculator.computeRisk();
	            riskConvolutionPortfolio.assertConvolutionExecutedAs(ConvolutionMode.PER_RUPTURE);
	            riskConvolutionPortfolio.exportRuptureResults(perRuptureOutputCSV);
	            
	            riskConvolutionPortfolio.getRuptureResults().printSummary();
	            
	            System.out.println("Outputs written to run folder:");
	            System.out.println(perRuptureOutputCSV);
        	
	        }
        }
    }
}