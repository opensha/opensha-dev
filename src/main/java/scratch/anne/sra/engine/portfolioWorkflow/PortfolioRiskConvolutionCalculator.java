package scratch.anne.sra.engine.portfolioWorkflow;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;

import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.domain.hazard.ErfBuilder;
import scratch.anne.sra.domain.hazard.HazardCurve;
import scratch.anne.sra.domain.hazard.HazardCurveCollection;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolioWrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponseLibrary;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.engine.hazardSetup.AssertHazardExceedance;
import scratch.anne.sra.engine.hazardSetup.GmmImtSetter;
import scratch.anne.sra.engine.hazardSetup.SiteBuilder;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;
import scratch.anne.sra.util.StringUtil.ImtPeriod;

/**
 * Portfolio workflow:
 * <ul>
 *   <li>Parallelization is performed at SiteKey level</li>
 *   <li>Loop over all siteKeys and imKeys to minimize hazard calculations</li>
 *   <li>InnerLoop over all relevant assets for risk calcs</li>
 * </ul>

 */
public class PortfolioRiskConvolutionCalculator {

    private final RiskConvolutionPortfolio portfolio;
    private final SimpleImResponseLibrary responseLib;
    private final HazardParameters hazardParameters;
    private final RiskConvolution.IntegrationMethod integrationMethod;
    
    private final AbstractERF erf;
    private final AttenRelRef gmmRef;
    
    private final SourceFilterManager filterManager;

    HazardCurveCollection hazardCurves;


    /** full constructor */
    public PortfolioRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters,
            RiskConvolution.IntegrationMethod integrationMethod
    ) {
        this.portfolio = portfolio;
        this.responseLib = responseLib;
        this.hazardParameters = hazardParameters;
        this.integrationMethod = integrationMethod;

        AssertHazardExceedance.require(hazardParameters.getHazardMetric());

    	this.erf = ErfBuilder.build(hazardParameters);

    	this.filterManager = hazardParameters.buildSourceManager();
      
        this.gmmRef = AttenRelRef.valueOf(hazardParameters.getGmmName().toUpperCase());
    }

    /** default integration method */
    public PortfolioRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters
    ) {
        this(portfolio, responseLib, hazardParameters,
                RiskConvolution.DEFAULT_METHOD);
    }

    /** hazard and risk calculation loops */
    public RiskConvolutionPortfolio computeRisk() {
    	
    	hazardCurves = new HazardCurveCollection(hazardParameters);

        ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
        ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();

        ScalarIMR baseGmm = gmmRef.get();

        List<SiteKey> siteKeys = new ArrayList<>(portfolio.getSiteKeys());

        int totalSiteImKeys = countTotalSiteImKeys();
        AtomicInteger counter = new AtomicInteger();

        long startTime = System.nanoTime();
        
        System.out.printf(
                "%nRunning hazard curve loop for %d Site & Im Keys %n",
                totalSiteImKeys
        );

        List<Site> sites = new ArrayList<>();

        for (SiteKey siteKey : siteKeys) {
            sites.add(SiteBuilder.fromSiteKey(siteKey, baseGmm));
        }

        List<CompletableFuture<Void>> futures = new ArrayList<>();

        for (int i = 0; i < sites.size(); i++) {

            final Site site = sites.get(i);
            final SiteKey siteKey = siteKeys.get(i);

            futures.add(CompletableFuture.runAsync(() -> {

                ScalarIMR gmm = null;
                HazardCurveCalculator calc = null;

                try {

                    synchronized (gmmDeque) {
                        gmm = gmmDeque.isEmpty() ? gmmRef.get() : gmmDeque.pop();
                    }

                    synchronized (calcDeque) {
                    	calc = calcDeque.isEmpty()
                    			// filter erf when doing the hazard calculation
                    		    ? new HazardCurveCalculator(filterManager)
                    		    : calcDeque.pop();
                    }

                    gmm.setSite(site);

                    for (ImKey imKey : portfolio.getImKeysBySite(siteKey)) {

                        ImtPeriod imt = imKey.getImtPeriod();
                        GmmImtSetter.setImt(gmm, imt);

                        DiscretizedFunc hazFunc =
                                new ArbitrarilyDiscretizedFunc();

                        for (double x : imKey.getLogValues())
                            hazFunc.set(x, 0d);
                       
                        
                        hazFunc = calc.getHazardCurve(hazFunc, site, gmm, erf);

                        if (hazardParameters.getHazardMetric() == scratch.anne.sra.util.enums.HazardMetric.RATE_EXCEEDANCE) {
                            // average rate over full ERF duration (not annualized)
                            hazFunc = calc.getAnnualizedRates(hazFunc, 1.0);
                        }
                        // else PROBABILITY_EXCEEDANCE: already correct
                        
                        double[] hazard = new double[hazFunc.size()];
                        for (int j = 0; j < hazard.length; j++)
                            hazard[j] = hazFunc.getY(j);

                        hazardCurves.put(
                                siteKey,
                                imKey,
                                new HazardCurve(hazard)
                        );
                        
                        for (RiskConvolutionAsset asset :
                            portfolio.getAssetsBySiteAndImKey(siteKey, imKey)) {

	                        RiskConvolution riskCalc = new RiskConvolution(
	                                responseLib.getByName(asset.getModelName()),
	                                hazard,
	                                integrationMethod
	                        );
	
	                        asset.setRiskConvolutionResult(riskCalc.compute());
                        }

                        // ----- PROGRESS -----
                        int done = counter.incrementAndGet();

                        if (done % 100 == 0 || done == totalSiteImKeys) {

                            long now = System.nanoTime();

                            double elapsed = (now - startTime) / 1e9;
                            double rate = done / elapsed;
                            double etaSeconds = (totalSiteImKeys - done) / rate;

                            LocalTime currentTime = LocalTime.now();
                            LocalTime etaTime = currentTime.plusSeconds((long) etaSeconds);

                            DateTimeFormatter fmt =
                                    DateTimeFormatter.ofPattern("h:mm:ss a");

                            System.out.printf(
                                    "Hazard %d / %d (%.1f%%) | %.1f sites/s | Now %s | ETA %s (%.1f hr remaining)%n",
                                    done,
                                    totalSiteImKeys,
                                    100.0 * done / totalSiteImKeys,
                                    rate,
                                    currentTime.format(fmt),
                                    etaTime.format(fmt),
                                    etaSeconds / 60.0 / 60.0
                            );
                        }
                    }

                } finally {
                    if (gmm != null) synchronized (gmmDeque) { gmmDeque.push(gmm); }
                    if (calc != null) synchronized (calcDeque) { calcDeque.push(calc); }
                }
            }));
        }

        futures.forEach(CompletableFuture::join);

        hazardCurves.freeze();
        portfolio.setHazardCurves(hazardCurves);
        
        portfolio.setConvolutionMode(ConvolutionMode.FULL_HCURVE);
       
        
        double elapsedSeconds = (System.nanoTime() - startTime) / 1e9;

        if (elapsedSeconds >= 3600) {
            System.out.printf("Completed in %.1f hours%n",
                    elapsedSeconds / 3600.0);
        } else if (elapsedSeconds >= 60) {
            System.out.printf("Completed in %.1f minutes%n",
                    elapsedSeconds / 60.0);
        } else {
            System.out.printf("Completed in %.1f seconds%n",
                    elapsedSeconds);
        }
        
        
        
        return portfolio;
    }
    
    private int countTotalSiteImKeys() {

        int total = 0;

        for (SiteKey siteKey : portfolio.getSiteKeys()) {
            total += portfolio.getImKeysBySite(siteKey).size();
        }

        return total;
    }

}