package scratch.anne.sra.examples.PBR;

import java.nio.file.*;
import java.util.*;

import scratch.anne.sra.domain.asset.fragility.FragilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.PBRSurvivalPortfolioReporter;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structural_response.fragilities.FragilityModel;
import scratch.anne.sra.domain.structural_response.fragilities.SimpleImFragilityLibraryPreparer;
import scratch.anne.sra.domain.structural_response.fragilities.SimpleImFragilityLibraryPreparer.PrepareSpec;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioRiskConvolutionCalculator;
import scratch.anne.sra.util.InputParsingUtil;
import scratch.anne.sra.io.readers.FragilityLibraryReader;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.util.IO;
import scratch.anne.sra.util.InputParsingUtil;

public class LJBPortfolioPFail {

    public static void main(String[] args) throws Exception {
    	
    	Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\precarious-balanced-rock\\results\\");
    	Path inputFolder = Path.of("LJBport_all-LJB1frag");
    	
    	boolean writeHazard = true;

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

        Map<String, String> filterConfig =
                ParseRiskRunParametersCSV.extractSourceFilterConfig(config);

        // ------------------- 1. Resolve input files -------------------
        Path portfolioCSV = Path.of(config.get("portfolio_csv"));
        Path fragLibraryJSON = Path.of(config.get("fragility_library_json"));

        IO.verifyFileExists(portfolioCSV, "Portfolio CSV");
        IO.verifyFileExists(fragLibraryJSON, "Fragility JSON");

        // ------------------- 2. create output paths inside run folder -------------------
        String fileTag = config.getOrDefault("file_tag", runFolder.getFileName().toString());
        Path outputCSV = runFolder.resolve(fileTag + ".csv");
        Path hazardJson = writeHazard
                ? runFolder.resolve(fileTag + "_hazard-list.json")
                : null;
        
        List<Double> probabilityTargets =
                InputParsingUtil.parseDoubleList(config.get("probability_targets"));
        if (probabilityTargets.isEmpty()) {
            throw new IllegalArgumentException(
                "Config missing 'probability_targets'");
        }
        
        
        // ------------------- 3. prepare hazard and integration parameters ----------------
        HazardParameters hazardParameters = new HazardParameters(
                config.get("erf_class"),
                InputParsingUtil.parseDoubleOrDefault(
                        config.get("erf_duration"),
                        1.0
                ),
                ParseRiskRunParametersCSV.parseHazardMetric(config.get("hazard_metric")),
                config.get("gmm"),
                filterConfig
        );
        
        RiskConvolution.IntegrationMethod integrationMethod = ParseRiskRunParametersCSV.parseIntegrationMethod(config.get("integration_method"));
        double logImStep = InputParsingUtil.parseDoubleOrDefault(config.get("log_im_step"), Double.NaN);

        // ------------------- 4. Load portfolio -------------------
        Portfolio<FragilityAsset> basePortfolio =
                PortfolioReader.readCSV(portfolioCSV, FragilityAsset.class);

        Set<String> portfolioFragilities = basePortfolio.getResponseModelNames();

        // ------------------- 5. Load fragility library -------------------
        ResponseModelLibrary<FragilityModel> baseFragilityLib =
                FragilityLibraryReader.readLibrary(fragLibraryJSON);

        // ------------------- 6. Prepare SimpleFragility Library -------------------
        PrepareSpec spec = PrepareSpec.builder()
        		.filterNames(portfolioFragilities).build();
        SimpleImResponseLibrary fragilityLib = SimpleImFragilityLibraryPreparer.prepare(baseFragilityLib, spec);

        // ------------------- 7. Wrap portfolio with IMKey mapping -------------------
        RiskConvolutionPortfolio riskConvolutionPortfolio = new RiskConvolutionPortfolio(basePortfolio, fragilityLib);

        // ------------------- 8. Create calculator -------------------       
        PortfolioRiskConvolutionCalculator calculator =
                new PortfolioRiskConvolutionCalculator(
                        riskConvolutionPortfolio,
                        fragilityLib,
                        hazardParameters,
                        integrationMethod
                );

        // ------------------- 9. Compute Probability of Failure -------------------
        riskConvolutionPortfolio = calculator.computeRisk();
        riskConvolutionPortfolio.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
        riskConvolutionPortfolio.exportHazard(hazardJson);
        
	    // ------------------- 10. Write PBR outputs to run folder -------------------
	    PBRSurvivalPortfolioReporter PBRreporter = new PBRSurvivalPortfolioReporter(riskConvolutionPortfolio, probabilityTargets);
	    
	    PBRreporter.writeCSV(outputCSV);

        System.out.println("Outputs written to run folder:");
        System.out.println(outputCSV);
        System.out.println(hazardJson);
        
    }
    
}