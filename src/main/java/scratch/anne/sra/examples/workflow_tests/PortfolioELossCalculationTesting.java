package scratch.anne.sra.examples.workflow_tests;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.util.enums.HazardMetric;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.ExpectedLossPortfolioAggregator;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structural_response.vulnerabilities.SimpleImVulnLibraryPreparer;
import scratch.anne.sra.domain.structural_response.vulnerabilities.VulnerabilityModel;
import scratch.anne.sra.engine.convolution.RiskConvolution.IntegrationMethod;
import scratch.anne.sra.engine.per_rupture.RuptureLossResultsCollection;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioPerRuptureRiskConvolutionCalculator;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioRiskConvolutionCalculator;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.io.readers.VulnerabilityLibraryReader;
import scratch.anne.sra.util.InputParsingUtil;


/**
 * Sensitivity check: runs a portfolio's expected-loss calculation across
 * every combination of ERF duration, hazard metric, integration method,
 * and convolution mode, appending each combination's total expected
 * loss to {@code sensitivity_summary.csv} in the run folder.
 *
 * <p>Column headers in that CSV are generic ({@code Variable1}..
 * {@code Variable4}, {@code totalValue}) with headroom for more
 * variables than the 4 currently swept, though none are currently
 * planned.</p>
 *
 * <p><b>Choice point -- edit before running:</b> {@code erfDurations}/
 * {@code hazardMetrics}/{@code integrationMethods}/{@code convolutionModes}
 * select which values of each dimension are swept (alternatives
 * commented out above the active ones); {@code baseFolder}/
 * {@code inputFolder} select which run folder's {@code common_input.csv}
 * supplies the portfolio/vulnerability library/ERF class/GMM.</p>
 */
public class PortfolioELossCalculationTesting {
	
	public static void main(String[] args) throws Exception {
		
//		List<Integer>             erfDurations = List.of(1, 25);
//		List<HazardMetric>        hazardMetrics = List.of(HazardMetric.PROBABILITY_EXCEEDANCE, HazardMetric.RATE_EXCEEDANCE);
//		List<IntegrationMethod>   integrationMethods = List.of(IntegrationMethod.CLOSED_FORM, IntegrationMethod.RIEMANN);
//		List<ConvolutionMode>     convolutionModes = List.of(ConvolutionMode.FULL_HCURVE, ConvolutionMode.PER_RUPTURE);
		
		List<Integer>             erfDurations = List.of(1);
		List<HazardMetric>        hazardMetrics = List.of(HazardMetric.RATE_EXCEEDANCE);
		List<IntegrationMethod>   integrationMethods = List.of(IntegrationMethod.RIEMANN);
		List<ConvolutionMode>     convolutionModes = List.of(ConvolutionMode.PER_RUPTURE);
//		List<ConvolutionMode>     convolutionModes = List.of(ConvolutionMode.FULL_HCURVE);
		
		// ---------------------------------------------
		
		Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\my_scratch\\conversion to Java project\\BERM_test-outputs\\_vb\\csv_inputs");
		Path inputFolder = Path.of("tests\\analysis_flag_sensitivity");
	    Path runFolder = baseFolder.resolve(inputFolder);
	    
	    if (!Files.isDirectory(runFolder))
	        throw new IllegalArgumentException("Run folder does not exist: " + runFolder);
	
	    Path configCSV = runFolder.resolve("common_input.csv");
	    if (!Files.exists(configCSV))
	        throw new IllegalArgumentException("Missing input.csv in run folder: " + runFolder);
	    
	    Map<String, String> config =
	            ParseRiskRunParametersCSV.loadConfig(configCSV);
	    
	    String gmmName = config.get("gmm");
	    
	    // ---------------------------------------------
	    
        Path portfolioCSV = Path.of(config.get("portfolio_csv"));
        Path vulnLibraryJSON = Path.of(config.get("vuln_library_json"));

        InputParsingUtil.verifyFileExists(portfolioCSV, "Portfolio CSV");
        InputParsingUtil.verifyFileExists(vulnLibraryJSON, "Vulnerability JSON");
        
        Portfolio<VulnerabilityAsset> basePortfolio =
                PortfolioReader.readCSV(portfolioCSV, VulnerabilityAsset.class);

        ResponseModelLibrary<VulnerabilityModel> vulnLib =
                VulnerabilityLibraryReader.readLibrary(vulnLibraryJSON);

        SimpleImResponseLibrary expVulnLib =
                SimpleImVulnLibraryPreparer.prepare(vulnLib, basePortfolio.getResponseModelNames());
	    
	    // --------------------------------------------------
        
        Path summaryFile = runFolder.resolve("sensitivity_summary.csv");
        boolean fileExists = Files.exists(summaryFile);

        try (BufferedWriter w = Files.newBufferedWriter(
                summaryFile,
                java.nio.file.StandardOpenOption.CREATE,
                java.nio.file.StandardOpenOption.APPEND)) {

            if (!fileExists) {
                writeHeader(w);
                w.flush();
            }
            
            for (int erfDuration : erfDurations) {
            	for (HazardMetric hazardMetric : hazardMetrics) {
            		for (IntegrationMethod integrationMethod: integrationMethods) {
            			for (ConvolutionMode convolutionMode: convolutionModes) {
            				
            				List<String> variableNames = List.of(
            				        String.valueOf(erfDuration)+"-yr",
            				        String.valueOf(hazardMetric),
            				        String.valueOf(integrationMethod),
            				        String.valueOf(convolutionMode)
            				        );
            				Path outputCSV = runFolder.resolve(String.join("_", variableNames) + ".csv");
            				System.out.println("Running: " + String.join("_", variableNames) + "\n");
            				
            				RiskConvolutionPortfolio riskConvolutionPortfolio = new RiskConvolutionPortfolio(basePortfolio, expVulnLib);
            				
            		        HazardParameters hazardParameters = new HazardParameters(
            		                config.get("erf_class"),
            		                erfDuration,
            		                hazardMetric,
            		                gmmName
            		        );
            		        
            		        if (convolutionMode == ConvolutionMode.FULL_HCURVE) {
            		        	PortfolioRiskConvolutionCalculator calculator =
            		                    new PortfolioRiskConvolutionCalculator(
            		                            riskConvolutionPortfolio,
            		                            expVulnLib,
            		                            hazardParameters,
            		                            integrationMethod
            		                    );
            		            
            		            riskConvolutionPortfolio = calculator.computeRisk();
            		            riskConvolutionPortfolio.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
            		            ExpectedLossPortfolioAggregator aggregator = new ExpectedLossPortfolioAggregator(riskConvolutionPortfolio);
            		            aggregator.writeAggregatedCSV(outputCSV);
            		            
            		            double totalValue = aggregator.getTotalExpectedLoss();
            		            
            		            System.out.printf(
            		            	    "%n%n%s%n\tTotalLoss: %.2g%n%n%n",
            		            	    String.join("_", variableNames),
            		            	    totalValue
            		            	);
            		            
            		            writeRow(w, variableNames, totalValue);
            		            w.flush();
            		        }
            		        else if (convolutionMode == ConvolutionMode.PER_RUPTURE) {           		        	
            		            PortfolioPerRuptureRiskConvolutionCalculator calculator =
            		                    new PortfolioPerRuptureRiskConvolutionCalculator(
            		                            riskConvolutionPortfolio,
            		                            expVulnLib,
            		                            hazardParameters,
            		                            integrationMethod
            		                    );
            		            
            		            riskConvolutionPortfolio = calculator.computeRisk();
            		            riskConvolutionPortfolio.assertConvolutionExecutedAs(ConvolutionMode.PER_RUPTURE);
            		            riskConvolutionPortfolio.exportRuptureResults(outputCSV);
            		            
            		            RuptureLossResultsCollection ruptureResults = riskConvolutionPortfolio.getRuptureResults();
            		            double totalValue = ruptureResults.getTotalExpectedLoss();
            		            
            		            String topContributors = ruptureResults.getTopContributors(5).stream()
            		            	    .map(r -> String.format(
            		            	            "\t\tSource: %d, Rup: %d, loss=%.2g, relative=%.2g",
            		            	            r.key().sourceId(),
            		            	            r.key().ruptureId(),
            		            	            r.lossContribution(),
            		            	    		r.relativeContribution()))
            		            	    .collect(Collectors.joining("\n"));
            		            

								System.out.printf(
								    "%n%n%s%n\tTotalLoss: %.2f%n\tTopContributors:%n%s%n%n%n",
								    String.join("_", variableNames),
								    totalValue,
								    topContributors
								);
            		            
            		            writeRow(w, variableNames, totalValue);
            		            w.flush();
            		        }
            			}
            		}
            	}
            	
            }
        }
	}
	
    private static final int NUM_COLUMNS = 10;
    private static final String VALUE_COL = "totalValue"; 
	
	private static void writeHeader(BufferedWriter w) throws IOException {
	    String[] header = new String[NUM_COLUMNS];

	    for (int i = 0; i < NUM_COLUMNS - 1; i++) {
	        header[i] = "Variable" + (i + 1);
	    }

	    header[NUM_COLUMNS - 1] = VALUE_COL;

	    w.write(String.join(",", header));
	    w.newLine();
	}
	
	private static void writeRow(
	        BufferedWriter w,
	        List<String> variableNames,
	        Object value) throws IOException {

	    String[] row = new String[NUM_COLUMNS];

	    int i = 0;

	    for (; i < variableNames.size() && i < NUM_COLUMNS - 1; i++) {
	        row[i] = variableNames.get(i);
	    }

	    while (i < NUM_COLUMNS - 1) {
	        row[i++] = "";
	    }

	    row[NUM_COLUMNS - 1] = String.valueOf(value);

	    w.write(String.join(",", row));
	    w.newLine();
	}
}
