package scratch.anne.sra.engine.portfolio_workflow;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.sha.calc.sourceFilters.SourceFilter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilterUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy.CacheTypes;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;

import scratch.anne.sra.domain.asset.RiskConvolutionAsset;
import scratch.anne.sra.domain.asset.RiskConvolutionAsset.RiskMetricType;
import scratch.anne.sra.domain.asset.vulnerability.ExpectedLossAsset;
import scratch.anne.sra.domain.hazard.ErfBuilder;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio;
import scratch.anne.sra.domain.portfolio.portfolio_wrappers.RiskConvolutionPortfolio.ConvolutionMode;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.engine.convolution.RiskConvolution;
import scratch.anne.sra.engine.hazard_setup.GmmImtSetter;
import scratch.anne.sra.engine.hazard_setup.AssertHazardExceedance;
import scratch.anne.sra.engine.hazard_setup.SiteBuilder;
import scratch.anne.sra.engine.per_rupture.RuptureLossResultsCollection;
import scratch.anne.sra.io.writers.AssetRiskForRuptureWriter;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;
import scratch.anne.sra.util.StringUtil.ImtPeriod;

/**
 * Portfolio workflow:
 * <ul>
 *   <li>Parallelization is performed at SiteKey level</li>
 *   <li>Loop over all siteKeys and imKeys to minimize hazard calculations</li>
 *   <li>Loop over each rupture to accumulate risk results
 *   <li>InnerLoop over all relevant assets for risk calcs</li>
 * </ul>

 */
public class PortfolioPerRuptureRiskConvolutionCalculator {

    private final RiskConvolutionPortfolio portfolio;
    private final SimpleImResponseLibrary responseLib;
    private final HazardParameters hazardParameters;
    private final RiskConvolution.IntegrationMethod integrationMethod;
    
    private final RiskMetricType riskMetricType;
    
    private final AbstractERF erf;
    private final AttenRelRef gmmRef;
    
    private final SourceFilterManager filterManager;
    private final List<SourceFilter> sourceFilters;

    RuptureLossResultsCollection ruptureResults;
    private final AssetRiskForRuptureWriter assetWriter;

    /** full constructor */
    public PortfolioPerRuptureRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters,
            RiskConvolution.IntegrationMethod integrationMethod,
            AssetRiskForRuptureWriter assetWriter
    ) {
    	this.portfolio = portfolio;
    	this.riskMetricType = portfolio.getRiskMetricType();
    	
    	// if not calculating expected loss, the only save results are asset level
    	if (riskMetricType != RiskMetricType.EXPECTED_LOSS) {
    		throw new IllegalStateException("Per-rupture probability of failure assessment save asset-level results");
    	}
        
    	this.responseLib = responseLib;
    	this.hazardParameters = hazardParameters;
    	this.integrationMethod = integrationMethod;

    	AssertHazardExceedance.require(hazardParameters.getHazardMetric());

    	SurfaceCachingPolicy.force(CacheTypes.THREAD_LOCAL);

    	this.erf = ErfBuilder.build(hazardParameters);

    	this.filterManager = hazardParameters.buildSourceManager();
    	this.sourceFilters = filterManager.getEnabledFilters();
      
        this.gmmRef = AttenRelRef.valueOf(hazardParameters.getGmmName().toUpperCase());
        
        this.assetWriter = assetWriter;

    }

    /** specify integration method but no asset writer */
    public PortfolioPerRuptureRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters,
            RiskConvolution.IntegrationMethod integrationMethod
    ) {
        this(portfolio, responseLib, hazardParameters,
                integrationMethod,
                null);
    }
    
    /** default integration method with asset writer */
    public PortfolioPerRuptureRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters,
            AssetRiskForRuptureWriter assetWriter
    ) {
        this(portfolio, responseLib, hazardParameters,
                RiskConvolution.DEFAULT_METHOD,
                assetWriter);
    }
    
    /** default integration method and no asset writer */
    public PortfolioPerRuptureRiskConvolutionCalculator(
            RiskConvolutionPortfolio portfolio,
            SimpleImResponseLibrary responseLib,
            HazardParameters hazardParameters
    ) {
        this(portfolio, responseLib, hazardParameters,
                RiskConvolution.DEFAULT_METHOD,
                null);
    }

    /** hazard and risk calculation loops */
    public RiskConvolutionPortfolio computeRisk() {
    	
        System.out.printf(
                "%nInitializing rupture loss results collection...%n"
        );    	
    	ruptureResults = RuptureLossResultsCollection.fromERF(erf, hazardParameters);    	

        ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
        ScalarIMR baseGmm = gmmRef.get();
        List<SiteKey> siteKeys = new ArrayList<>(portfolio.getSiteKeys());

        int totalSiteImKeys = countTotalSiteImKeys();
        int counterIncrement = 10;
        AtomicInteger counter = new AtomicInteger();
        long startTime = System.nanoTime();
        
        System.out.printf(
                "%nRunning per-rupture loop for %d Site & Im Keys %n",
                totalSiteImKeys
        );
        
        
        final boolean thresholdedAssetWriter = assetWriter != null
                && riskMetricType == RiskMetricType.EXPECTED_LOSS;
        final double totalPortfolioValue = thresholdedAssetWriter
                ? portfolio.getTotalPortfolioValue()
                : 0d;        

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

                try {
                    synchronized (gmmDeque) {
                        gmm = gmmDeque.isEmpty() ? gmmRef.get() : gmmDeque.pop();
                    }

                    gmm.setSite(site);

                    for (ImKey imKey : portfolio.getImKeysBySite(siteKey)) {

                        ImtPeriod imt = imKey.getImtPeriod();
                        GmmImtSetter.setImt(gmm, imt);

                        DiscretizedFunc hazFunc =
                                new ArbitrarilyDiscretizedFunc();

                        for (double x : imKey.getLogValues())
                            hazFunc.set(x, 0d);
                        
                        
                        for (int sourceID = 0; sourceID < erf.getNumSources(); sourceID++) {
                            ProbEqkSource source = erf.getSource(sourceID);
                            if (SourceFilterUtils.canSkipSource(sourceFilters, source, site)) {
                                // source is outside filter limits
                                continue;
                            }          
                            
                            for (int ruptureID = 0; ruptureID < source.getNumRuptures(); ruptureID++) {
                                ProbEqkRupture rupture = source.getRupture(ruptureID);
	                            if (SourceFilterUtils.canSkipRupture(sourceFilters, rupture, site)) {
	                                // rupture is outside filter limits
	                                continue;
	                            }
	
	                            gmm.getExceedProbabilities(rupture, hazFunc);
	                            // rate vs probability is addressed outside, in 
		                        // the rupture likelihood, not the curve itself

		                        double[] hazard = new double[hazFunc.size()];
		                        for (int j = 0; j < hazard.length; j++)
		                            hazard[j] = hazFunc.getY(j);
			                        
		                        final double siteImValue = thresholdedAssetWriter
		                                ? portfolio.getTotalSiteImValue(siteKey, imKey)
		                                : 0d;
		                        Map<Integer, Double> assetLosses = thresholdedAssetWriter ? new HashMap<>() : null;
	                        	
		                        for (RiskConvolutionAsset asset :
		                            portfolio.getAssetsBySiteAndImKey(siteKey, imKey)) {
		                        	
			                        RiskConvolution riskCalc = new RiskConvolution(
			                                responseLib.getByName(asset.getModelName()),
			                                hazard,
			                                integrationMethod
			                        );
			                        
                                   if (riskMetricType == RiskMetricType.EXPECTED_LOSS) {
                                        ExpectedLossAsset vulnAsset = (ExpectedLossAsset) asset;
                                        double assetLossForRupture = vulnAsset.getValue() * riskCalc.compute().getRisk();
                                        ruptureResults.accumulateRuptureLoss(new RuptureKey(sourceID, ruptureID), assetLossForRupture);
                                        if (thresholdedAssetWriter) {
                                            assetLosses.put(asset.getAssetID(), assetLossForRupture);
                                        }
                                    } else {
                                        // not expected loss — always write, no threshold
                                        assetWriter.write(asset.getAssetID(), sourceID, ruptureID, riskCalc.compute().getRisk());
                                    }
                                }

		                        if (thresholdedAssetWriter) {
		                            double siteImLoss = assetLosses.values().stream()
		                                    .mapToDouble(Double::doubleValue).sum();
		                            double threshold = assetWriter.writeThresholdValue(siteImValue, totalPortfolioValue);
		                            if (siteImLoss >= threshold) {
		                                for (Map.Entry<Integer, Double> entry : assetLosses.entrySet()) {
		                                    assetWriter.write(entry.getKey(), sourceID, ruptureID, entry.getValue());
		                                }
		                            }
		                        }
	                        }
                        }

                        // ----- PROGRESS -----
                        int done = counter.incrementAndGet();

                        if (done % counterIncrement == 0 || done == totalSiteImKeys) {

                            long now = System.nanoTime();

                            double elapsed = (now - startTime) / 1e9;
                            double rate = done / elapsed;
                            double etaSeconds = (totalSiteImKeys - done) / rate;

                            LocalTime currentTime = LocalTime.now();
                            LocalTime etaTime = currentTime.plusSeconds((long) etaSeconds);

                            DateTimeFormatter fmt =
                                    DateTimeFormatter.ofPattern("h:mm:ss a");

                            System.out.printf(
                                    "Site&IM %d / %d (%.1f%%) | %.2g sites/s | Now %s | ETA %s (%.1f hr remaining)%n",
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

                }
            }));
        }

        futures.forEach(CompletableFuture::join);

        ruptureResults.freeze();
        portfolio.setRuptureResults(ruptureResults);
        
        portfolio.setConvolutionMode(ConvolutionMode.PER_RUPTURE);      
        
        double elapsedSeconds = (System.nanoTime() - startTime) / 1e9;

        if (elapsedSeconds >= 3600) {
	            System.out.printf("Per-rupture expected loss assessment completed in %.1f hours%n",
	                    elapsedSeconds / 3600.0);
	        } else if (elapsedSeconds >= 60) {
	            System.out.printf("Per-rupture expected loss assessment completed in %.1f minutes%n",
	                    elapsedSeconds / 60.0);
	        } else {
	            System.out.printf("Per-rupture expected loss assessment completed in %.1f seconds%n",
	                    elapsedSeconds);
	        }
        
        
        if (assetWriter != null) {
        	assetWriter.close();
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
