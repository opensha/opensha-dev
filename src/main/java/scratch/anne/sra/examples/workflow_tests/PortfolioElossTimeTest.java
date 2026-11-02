package scratch.anne.sra.examples.workflow_tests;

import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

import scratch.anne.sra.domain.asset.vulnerability.VulnerabilityAsset;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.Portfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.ExpectedLossPortfolioAggregator;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structural_response.vulnerabilities.SimpleImVulnLibraryPreparer;
import scratch.anne.sra.domain.structural_response.vulnerabilities.VulnerabilityModel;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.engine.convolution.RiskConvolution.IntegrationMethod;
import scratch.anne.sra.engine.per_rupture.RuptureLossResultsCollection;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioPerRuptureRiskConvolutionCalculator;
import scratch.anne.sra.engine.portfolio_workflow.PortfolioRiskConvolutionCalculator;
import scratch.anne.sra.io.readers.ParseRiskRunParametersCSV;
import scratch.anne.sra.io.readers.PortfolioReader;
import scratch.anne.sra.io.readers.VulnerabilityLibraryReader;
import scratch.anne.sra.io.writers.AssetRiskForRuptureWriter;
import scratch.anne.sra.io.writers.FileFormat;
import scratch.anne.sra.util.IO;
import scratch.anne.sra.util.ImValueTransformer;
import scratch.anne.sra.util.InputParsingUtil;

public class PortfolioElossTimeTest {

    public static Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\openSRA\\software architecture\\my_scratch\\conversion to Java project\\BERM_test-outputs\\_vb\\csv_inputs");
    public static Path inputFolder = Path.of("tests\\time-tests\\short_portfolio");

    public static boolean writeHazard = false;
    public static double fRelative = 0d;
    public static double fAbsolute = 0d;

    public static void main(String[] args) throws Exception {

        Path runFolder = baseFolder.resolve(inputFolder);

        if (!Files.isDirectory(runFolder))
            throw new IllegalArgumentException("Run folder does not exist: " + runFolder);

        Path configCSV = runFolder.resolve("input.csv");
        if (!Files.exists(configCSV))
            throw new IllegalArgumentException("Missing input.csv in run folder: " + runFolder);

        Map<String, String> config = ParseRiskRunParametersCSV.loadConfig(configCSV);
        Map<String, String> filterConfig = ParseRiskRunParametersCSV.extractSourceFilterConfig(config);

        Path portfolioCSV     = Path.of(config.get("portfolio_csv"));
        Path vulnLibraryJSON  = Path.of(config.get("vuln_library_json"));
        IO.verifyFileExists(portfolioCSV,    "Portfolio CSV");
        IO.verifyFileExists(vulnLibraryJSON, "Vulnerability JSON");

        String fileTag = config.getOrDefault("file_tag", runFolder.getFileName().toString());

        HazardParameters hazardParameters = new HazardParameters(
                config.get("erf_class"),
                InputParsingUtil.parseDoubleOrDefault(config.get("erf_duration"), 1.0),
                ParseRiskRunParametersCSV.parseHazardMetric(config.get("hazard_metric")),
                config.get("gmm"),
                filterConfig
        );

        Portfolio<VulnerabilityAsset> basePortfolio =
                PortfolioReader.readCSV(portfolioCSV, VulnerabilityAsset.class);

        ResponseModelLibrary<VulnerabilityModel> vulnLib =
                VulnerabilityLibraryReader.readLibrary(vulnLibraryJSON);

        SimpleImResponseLibrary expVulnLib =
                SimpleImVulnLibraryPreparer.prepare(
                        vulnLib,
                        basePortfolio.getResponseModelNames(),
                        new ImValueTransformer.NoImTransformation()
                );

        int nTimes = 10;

        // summary CSV — append if exists, write header only if new
        Path summaryFile = runFolder.resolve(fileTag + "_method_comparison.csv");
        boolean fileExists = Files.exists(summaryFile);

        try (var w = Files.newBufferedWriter(summaryFile,
                StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {

            if (!fileExists) {
                w.write("run,convolutionMode,integrationMethod,keepAssetLoss,totalExpectedLoss,elapsedSeconds");
                w.newLine();
            }

            for (int run = 1; run <= nTimes; run++) {
            System.out.printf("%n===== Run %d / %d =====%n", run, nTimes);

            // --- FULL_HCURVE x {CLOSED_FORM, RIEMANN} ---
            for (IntegrationMethod method : List.of(IntegrationMethod.CLOSED_FORM, IntegrationMethod.RIEMANN)) {

                RiskConvolutionPortfolio portfolio = new RiskConvolutionPortfolio(basePortfolio, expVulnLib);

                long t0 = System.nanoTime();
                PortfolioRiskConvolutionCalculator calc = new PortfolioRiskConvolutionCalculator(
                        portfolio, expVulnLib, hazardParameters, method);
                portfolio = calc.computeRisk();
                double elapsed = (System.nanoTime() - t0) / 1e9;

                portfolio.assertConvolutionExecutedAs(ConvolutionMode.FULL_HCURVE);
                ExpectedLossPortfolioAggregator agg = new ExpectedLossPortfolioAggregator(portfolio);
                double totalLoss = agg.getTotalExpectedLoss();

                Path outputCSV = runFolder.resolve(
                        fileTag + "_full-hcurve_" + method + ".csv");
                agg.writeAggregatedCSV(outputCSV);

                System.out.printf("FULL_HCURVE / %-12s  loss=%.4g  time=%.1fs%n",
                        method, totalLoss, elapsed);

                w.write(String.join(",",
                        String.valueOf(run),
                        "FULL_HCURVE", method.name(), "N/A",
                        String.valueOf(totalLoss),
                        String.format("%.2f", elapsed)));
                w.newLine();
                w.flush();
            }

            // --- PER_RUPTURE x {CLOSED_FORM, RIEMANN} x {no asset loss, with asset loss} ---
            for (IntegrationMethod method : List.of(IntegrationMethod.CLOSED_FORM, IntegrationMethod.RIEMANN)) {
                for (boolean keepAssetLoss : List.of(false, true)) {

                    AssetRiskForRuptureWriter assetWriter = keepAssetLoss
                            ? new AssetRiskForRuptureWriter(
                                    runFolder,
                                    fileTag + "_asset-risk-per-rup_" + method,
                                    fRelative, fAbsolute)
                            : null;

                    RiskConvolutionPortfolio portfolio = new RiskConvolutionPortfolio(basePortfolio, expVulnLib);

                    long t0 = System.nanoTime();
                    PortfolioPerRuptureRiskConvolutionCalculator calc =
                            new PortfolioPerRuptureRiskConvolutionCalculator(
                                    portfolio, expVulnLib, hazardParameters, method, assetWriter);
                    portfolio = calc.computeRisk();
                    double elapsed = (System.nanoTime() - t0) / 1e9;

                    portfolio.assertConvolutionExecutedAs(ConvolutionMode.PER_RUPTURE);

                    Path outputCSV = runFolder.resolve(
                            fileTag + "_per-rup_" + method
                            + (keepAssetLoss ? "_with-asset-loss" : "") + ".csv");
                    portfolio.exportRuptureResults(outputCSV);

                    double totalLoss = portfolio.getRuptureResults().getTotalExpectedLoss();

                    System.out.printf("PER_RUPTURE  / %-12s  assetLoss=%-5s  loss=%.4g  time=%.1fs%n",
                            method, keepAssetLoss, totalLoss, elapsed);

                    w.write(String.join(",",
                            String.valueOf(run),
                            "PER_RUPTURE", method.name(), String.valueOf(keepAssetLoss),
                            String.valueOf(totalLoss),
                            String.format("%.2f", elapsed)));
                    w.newLine();
                    w.flush();
                }
            }
            } // end run loop
        }

        System.out.println("\nSummary written to: " + summaryFile);
    }
}