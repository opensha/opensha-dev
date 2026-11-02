package scratch.anne.sra.examples.workflow_tests;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structural_response.vulnerabilities.SimpleImVulnLibraryPreparer;
import scratch.anne.sra.domain.structural_response.vulnerabilities.VulnerabilityModel;
import scratch.anne.sra.engine.convolution.RiskConvolution.IntegrationMethod;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioPerRuptureRiskConvolutionCalculator;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.io.readers.VulnerabilityLibraryReader;
import scratch.anne.sra.io.writers.AssetRiskForRuptureWriter;
import scratch.anne.sra.util.IO;
import scratch.anne.sra.util.InputParsingUtil;


public class PortfolioELossAssetFilterTesting {
	
	public static void main(String[] args) throws Exception {
		
//		double[] fRelativeGrid = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 0};
//		double[] fAbsoluteGrid = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 0};
		
//		double[] fRelativeGrid = {1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1};
//		double[] fAbsoluteGrid = {1e-1, 1e-2, 1e-3, 1e-4, 2e-1, 3e-1, 4e-1, 5e-1};
//		double[] fAbsoluteGrid = {0};
		
//		double[] fRelativeGrid = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 0};
//		double[] fAbsoluteGrid = {1e-1, 1e-2, 1e-3, 1e-4, 2e-1, 3e-1, 4e-1, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 0};
		
//		double[] fRelativeGrid = {0};
//		double[] fAbsoluteGrid = {0};
		
		double[] fRelativeGrid = {1e-1, 1e-2};
		double[] fAbsoluteGrid = {1e-3, 1e-4, 1e-5};
		
		String fileTag = "";
		
		// ---------------------------------------------
		
		Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\my_scratch\\conversion to Java project\\BERM_test-outputs\\_vb\\csv_inputs");
		Path inputFolder = Path.of("tests\\asset-per-rupture\\asset-filter-testing\\distributed_portfolio");
	    Path runFolder = baseFolder.resolve(inputFolder);
	    
	    if (!Files.isDirectory(runFolder))
	        throw new IllegalArgumentException("Run folder does not exist: " + runFolder);
	
	    Path configCSV = runFolder.resolve("common_input.csv");
	    if (!Files.exists(configCSV))
	        throw new IllegalArgumentException("Missing input.csv in run folder: " + runFolder);

        
	    // ---------------------------------------------
	    Map<String, String> config =
	            ParseRiskRunParametersCSV.loadConfig(configCSV);
	    
        Map<String, String> filterConfig =
                ParseRiskRunParametersCSV.extractSourceFilterConfig(config);
	    
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
        
        IntegrationMethod integrationMethod = ParseRiskRunParametersCSV.parseIntegrationMethod(config.get("integration_method"));
	    
	    // ---------------------------------------------
	    
        Path portfolioCSV = Path.of(config.get("portfolio_csv"));
        Path vulnLibraryJSON = Path.of(config.get("vuln_library_json"));

        IO.verifyFileExists(portfolioCSV, "Portfolio CSV");
        IO.verifyFileExists(vulnLibraryJSON, "Vulnerability JSON");
        
        Portfolio<VulnerabilityAsset> basePortfolio =
                PortfolioReader.readCSV(portfolioCSV, VulnerabilityAsset.class);

        ResponseModelLibrary<VulnerabilityModel> vulnLib =
                VulnerabilityLibraryReader.readLibrary(vulnLibraryJSON);

        SimpleImResponseLibrary expVulnLib =
                SimpleImVulnLibraryPreparer.prepare(vulnLib, basePortfolio.getResponseModelNames());
	    
	    // --------------------------------------------------
        
        for (double fRel : fRelativeGrid) {
            for (double fAbs : fAbsoluteGrid) {
            	String fileName = String.format("%s_fRel-%s_fAbs-%s", fileTag, sciTag(fRel), sciTag(fAbs));
            	Path outputCSV = runFolder.resolve(fileName + ".csv");
            	
            	if (!Files.exists(outputCSV)) {
	            	AssetRiskForRuptureWriter assetWriter = new AssetRiskForRuptureWriter(runFolder, fileName, fRel, fAbs);
	            	
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
		            riskConvolutionPortfolio.exportRuptureResults(outputCSV);
            	}
            }
        }
	}
	
	static String sciTag(double v) {
	    if (v == 0.0) return "0";
	    int exp = (int) Math.floor(Math.log10(v));
	    double coeff = v / Math.pow(10, exp);
	    long coeffRounded = Math.round(coeff);
	    return coeffRounded + "e" + exp;
	}
	

}
