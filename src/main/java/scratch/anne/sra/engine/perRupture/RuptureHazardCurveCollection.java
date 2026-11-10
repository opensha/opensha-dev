package scratch.anne.sra.engine.perRupture;

import java.lang.reflect.InvocationTargetException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.SourceFilter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilterUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;

import scratch.anne.sra.domain.asset.AbstractAsset;
import scratch.anne.sra.domain.hazard.HazardCurve;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.engine.hazardSetup.GmmImtSetter;
import scratch.anne.sra.engine.hazardSetup.SiteBuilder;
import scratch.anne.sra.util.enums.HazardMetric.HazardMetricType;
import scratch.anne.sra.util.enums.HazardMetric;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;
import scratch.anne.sra.util.StringUtil.ImtPeriod;

import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.PointSurface;

/**
 * Collects, per rupture in an ERF, the likelihood (rate or probability of
 * exceedance, per {@link HazardParameters}), magnitude, distance to a single
 * asset's site, and a conditional exceedance-probability curve (as a
 * {@link HazardCurve}) evaluated at a fixed IM grid via a GMM.
 *
 * <p>Single-asset / single-site / single-IMT. Modeled on the per-rupture
 * hazard loop in {@code PortfolioPerRuptureRiskConvolutionCalculator},
 * including source filtering and IMT/period setup, but without any
 * vulnerability/loss convolution — this produces hazard-per-rupture data
 * only.</p>
 *
 * <p>The {@link SiteKey} and {@link ImKey} used to build the collection are
 * stored as metadata (like {@link HazardParameters}), since they are fixed
 * for the whole collection rather than varying per rupture.</p>
 *
 * <p>Parallel in structure to {@link RuptureResultsCollection}: same
 * {@link RuptureKey}, same frozen-after-build discipline. There is no
 * accumulation phase here (unlike loss), so {@link #freeze()} computes two
 * comparison hazard curves rather than loss contributions.</p>
 */
public final class RuptureHazardCurveCollection {

    private final HazardParameters hazardParameters;
    private final SiteKey siteKey;
    private final ImKey imKey;

    private final Map<RuptureKey, Double> ruptureLikelihoods;
    private final Map<RuptureKey, Double> ruptureMagnitudes;
    private final Map<RuptureKey, Double> ruptureDistances;
    private final Map<RuptureKey, HazardCurve> ruptureCurves;

    // retained for the direct (HazardCurveCalculator) comparison computation in freeze();
    // not otherwise part of this collection's public state
    private final AbstractERF erf;
    private final ScalarIMR gmm;
    private final Site site;
    private final SourceFilterManager filterManager;

    private boolean frozen;
    private HazardCurve aggregateHazardCurve;
    private HazardCurve directHazardCurve;

    public RuptureHazardCurveCollection(
            HazardParameters hazardParameters,
            SiteKey siteKey,
            ImKey imKey,
            Map<RuptureKey, Double> ruptureLikelihoods,
            Map<RuptureKey, Double> ruptureMagnitudes,
            Map<RuptureKey, Double> ruptureDistances,
            Map<RuptureKey, HazardCurve> ruptureCurves,
            AbstractERF erf,
            ScalarIMR gmm,
            Site site,
            SourceFilterManager filterManager
    ) {
        this.hazardParameters = hazardParameters;

        if (hazardParameters.getHazardMetricType() != HazardMetricType.EXCEEDANCE) {
            throw new IllegalArgumentException(
                "RuptureProbabilityCurveCollection requires EXCEEDANCE hazard metric type, but got: "
                + hazardParameters.getHazardMetricType()
            );
        }

        if (hazardParameters.getHazardMetric() != HazardMetric.RATE_EXCEEDANCE) {
            throw new IllegalArgumentException(
                "RuptureProbabilityCurveCollection currently only supports RATE_EXCEEDANCE, but got: "
                + hazardParameters.getHazardMetric()
            );
        }

        this.siteKey = siteKey;
        this.imKey = imKey;

        this.ruptureLikelihoods = ruptureLikelihoods;
        this.ruptureMagnitudes = ruptureMagnitudes;
        this.ruptureDistances = ruptureDistances;
        this.ruptureCurves = ruptureCurves;

        this.erf = erf;
        this.gmm = gmm;
        this.site = site;
        this.filterManager = filterManager;
    }


    public double getLikelihood(RuptureKey key) {
        return ruptureLikelihoods.get(key);
    }

    public double getMagnitude(RuptureKey key) {
        return ruptureMagnitudes.get(key);
    }

    public double getDistance(RuptureKey key) {
        return ruptureDistances.get(key);
    }

    public HazardCurve getCurve(RuptureKey key) {
        return ruptureCurves.get(key);
    }

    public int getRuptureCount() {
        return ruptureLikelihoods.size();
    }

    /** Hazard generation metadata (ERF, GMM, duration, metric, filters) for this collection. */
    public HazardParameters getHazardParameters() {
        return hazardParameters;
    }

    /** The single site (lat/lon/vs30) this collection was evaluated at. */
    public SiteKey getSiteKey() {
        return siteKey;
    }

    /** The single IMT + IM grid this collection was evaluated at. */
    public ImKey getImKey() {
        return imKey;
    }


    /**
     * Immutable result for export (CSV / Parquet).
     */
    public record RuptureExceedanceResult(
            RuptureKey key,
            double likelihood,
            double magnitude,
            double distance,
            HazardCurve curve
    ) {}

    public RuptureExceedanceResult get(RuptureKey key) {
        return new RuptureExceedanceResult(
                key,
                ruptureLikelihoods.get(key),
                ruptureMagnitudes.get(key),
                ruptureDistances.get(key),
                ruptureCurves.get(key)
        );
    }

    public Iterable<RuptureExceedanceResult> results() {
        return ruptureLikelihoods.keySet().stream()
                .map(this::get)
                .toList();
    }


    /**
     * Freezes the collection and computes two hazard curves for comparison:
     * <ul>
     *   <li>the aggregate curve — for each IML, the rate-weighted sum of
     *       per-rupture conditional exceedance probabilities across all
     *       ruptures already collected in this object (see
     *       {@link #getAggregateHazardCurve()})</li>
     *   <li>the direct curve — computed independently via
     *       {@link HazardCurveCalculator#getHazardCurve}, not by looping over
     *       ruptures ourselves (see {@link #getDirectHazardCurve()})</li>
     * </ul>
     * The two are expected to agree closely; this is primarily a correctness
     * check on the manual per-rupture aggregation.
     *
     * <p>Only valid when likelihood is a rate — summing
     * probabilities directly is not a valid aggregation.</p>
     */
    public void freeze() {
        if (frozen)
            return;

        // ----- aggregate curve: rate-weighted sum of stored per-rupture curves -----
        int numImls = imKey.getValues().length;
        double[] aggregate = new double[numImls];

        for (RuptureKey key : ruptureLikelihoods.keySet()) {
            double likelihood = ruptureLikelihoods.get(key);
            double[] curve = ruptureCurves.get(key).getHazard();

            for (int i = 0; i < numImls; i++)
                aggregate[i] += likelihood * curve[i];
        }

        aggregateHazardCurve = new HazardCurve(aggregate);

        // ----- direct curve: computed independently via HazardCurveCalculator -----
        DiscretizedFunc hazFunc = new ArbitrarilyDiscretizedFunc();
        for (double x : imKey.getLogValues())
            hazFunc.set(x, 0d);

        HazardCurveCalculator calc = new HazardCurveCalculator(filterManager);
        hazFunc = calc.getHazardCurve(hazFunc, site, gmm, erf);
        // average rate over full ERF duration (not annualized)
        hazFunc = calc.getAnnualizedRates(hazFunc, 1.0);

        double[] direct = new double[hazFunc.size()];
        for (int i = 0; i < direct.length; i++)
            direct[i] = hazFunc.getY(i);

        directHazardCurve = new HazardCurve(direct);

        frozen = true;
    }

    /**
     * The aggregate hazard curve (rate-weighted sum of per-rupture curves,
     * one value per IML), computed by {@link #freeze()}.
     */
    public HazardCurve getAggregateHazardCurve() {
        if (!frozen)
            throw new IllegalStateException("Aggregate hazard curve not yet computed; call freeze() first");
        return aggregateHazardCurve;
    }

    /**
     * The hazard curve computed directly by {@link HazardCurveCalculator}
     * (not by looping over ruptures ourselves), for comparison against
     * {@link #getAggregateHazardCurve()}. The two should agree closely; any
     * material difference suggests a bug in the manual per-rupture
     * aggregation, the source filtering, or an inconsistency in how the GMM
     * was configured between the two paths.
     */
    public HazardCurve getDirectHazardCurve() {
        if (!frozen)
            throw new IllegalStateException("Direct hazard curve not yet computed; call freeze() first");
        return directHazardCurve;
    }


    /**
     * Convenience overload: builds both the ERF and the GMM from
     * {@code hazardParameters} internally (mirroring the ERF-instantiation
     * logic in {@code PortfolioPerRuptureRiskConvolutionCalculator}'s
     * constructor — reflection on {@code getErfName()}, duration set,
     * {@code updateForecast()} called), rather than requiring the caller to
     * build and pass one in.
     *
     * <p>Intended for single-asset, non-portfolio callers where there's no
     * benefit to sharing one ERF instance across multiple calls. A portfolio
     * wrapper that reuses one ERF across many assets/threads should instead
     * build the ERF once and call one of the {@code fromERF} overloads.</p>
     */
    public static RuptureHazardCurveCollection fromHazardParameters(
            HazardParameters hazardParameters,
            AbstractAsset asset,
            ImKey imKey
    ) {
        AbstractERF erf;
        try {
            erf = (AbstractERF)
                    Class.forName(hazardParameters.getErfName())
                            .getDeclaredConstructor()
                            .newInstance();
        } catch (ClassNotFoundException |
                 NoSuchMethodException |
                 IllegalAccessException |
                 InstantiationException |
                 InvocationTargetException e) {

            throw new IllegalArgumentException(
                    "Failed to instantiate ERF class: " + hazardParameters.getErfName(), e
            );
        }

        erf.getTimeSpan().setDuration(hazardParameters.getErfDuration());
        erf.updateForecast();

        return fromERF(erf, hazardParameters, asset, imKey);
    }

    /**
     * Convenience overload: builds a GMM from {@code hazardParameters.getGmmName()}
     * internally (with default params), rather than requiring the caller to
     * supply one.
     *
     * <p>Useful when the caller already has a shared/reusable ERF (e.g. built
     * once for a portfolio) but doesn't need to manage GMM lifecycle
     * themselves for this call — for example a single-threaded loop over a
     * handful of assets. A portfolio wrapper doing per-thread GMM pooling
     * (matching the deque pattern in the real per-rupture workflow) should
     * build its own {@code ScalarIMR} per thread and call the full
     * {@code fromERF(erf, gmm, hazardParameters, asset, imKey)} overload
     * instead.</p>
     */
    public static RuptureHazardCurveCollection fromERF(
            AbstractERF erf,
            HazardParameters hazardParameters,
            AbstractAsset asset,
            ImKey imKey
    ) {
        AttenRelRef gmmRef = AttenRelRef.valueOf(hazardParameters.getGmmName().toUpperCase());
        ScalarIMR gmm = gmmRef.get();
        gmm.setParamDefaults();

        return fromERF(erf, gmm, hazardParameters, asset, imKey);
    }

    /**
     * Builds a rupture probability curve collection for a single asset.
     *
     * <p>{@code erf} is expected to already be built and forecast-updated
     * (duration set, {@code updateForecast()} called) — it is not mutated
     * here and may be shared across concurrent calls from a portfolio
     * wrapper. {@code gmm} is expected to be caller-owned for this call (its
     * site/IMT/rupture state will be mutated here) — do not share a single
     * {@code ScalarIMR} instance across concurrent calls.</p>
     *
     * <p>{@code hazardParameters} is used for its hazard metric (rate vs.
     * probability, for likelihood only — the curve itself is always the
     * GMM's conditional exceedance probability) and its source filters, and
     * is stored on the resulting collection as metadata. It is expected to
     * be consistent with how {@code erf} and {@code gmm} were actually built,
     * but that consistency is the caller's responsibility — this method does
     * not re-derive or validate erf/gmm identity against it.</p>
     *
     * @param erf              pre-built, forecast-updated ERF; safe to share read-only
     * @param gmm              caller-owned GMM for this call; will be mutated (site, IMT, rupture)
     * @param hazardParameters must have an EXCEEDANCE metric type; supplies
     *                         hazard metric and source filters, stored as metadata
     * @param asset            supplies lat/lon/vs30 for the site; stored as a {@link SiteKey}
     * @param imKey            IMT (+ period, if SA) and IM grid to evaluate at; stored as metadata
     */
    public static RuptureHazardCurveCollection fromERF(
            AbstractERF erf,
            ScalarIMR gmm,
            HazardParameters hazardParameters,
            AbstractAsset asset,
            ImKey imKey
    ) {

        SourceFilterManager filterManager = hazardParameters.buildSourceManager();
        List<SourceFilter> sourceFilters = filterManager.getEnabledFilters();

        // ----------------- Site setup -----------------
        SiteKey siteKey = new SiteKey(asset.getLatitude(), asset.getLongitude(), asset.getVs30());

        Site site = SiteBuilder.fromSiteKey(siteKey, gmm);

        gmm.setSite(site);

        ImtPeriod imt = imKey.getImtPeriod();
        GmmImtSetter.setImt(gmm, imt);

        DiscretizedFunc hazFunc = new ArbitrarilyDiscretizedFunc();
        for (double x : imKey.getLogValues())
            hazFunc.set(x, 0d);

        // ----------------- Per-rupture loop -----------------
        Map<RuptureKey, Double> ruptureLikelihoods = new LinkedHashMap<>();
        Map<RuptureKey, Double> ruptureMagnitudes = new LinkedHashMap<>();
        Map<RuptureKey, Double> ruptureDistances = new LinkedHashMap<>();
        Map<RuptureKey, HazardCurve> ruptureCurves = new LinkedHashMap<>();

        for (int sourceID = 0; sourceID < erf.getNumSources(); sourceID++) {

            ProbEqkSource source = erf.getSource(sourceID);
            if (SourceFilterUtils.canSkipSource(sourceFilters, source, site))
                continue;

            int numRuptures = erf.getNumRuptures(sourceID);

            for (int ruptureID = 0; ruptureID < numRuptures; ruptureID++) {

                ProbEqkRupture rupture = source.getRupture(ruptureID);

                if (SourceFilterUtils.canSkipRupture(sourceFilters, rupture, site))
                    continue;

                // likelihood representation depends on the hazard metric (rate vs probability);
                // the curve itself (below) is always the GMM's conditional exceedance probability
                double likelihood;

                switch (hazardParameters.getHazardMetric()) {

                    case RATE_EXCEEDANCE:
                        // average rate over full ERF duration (not annualized)
                        likelihood = rupture.getMeanAnnualRate(1.0);
                        break;

                    case PROBABILITY_EXCEEDANCE:
                        likelihood = rupture.getProbability();
                        break;

                    default:
                        throw new IllegalStateException(
                                "Unsupported hazard metric: "
                                        + hazardParameters.getHazardMetric()
                        );
                }

                double magnitude = rupture.getMag();               
                // getDistanceRup() throws for PointSurface.DistanceCorrectable instances (e.g.
                // NSHM-style point-source approximations). Note: DistanceCorrectable's own
                // getUncorrectedDistances(Location) is broken in this version of OpenSHA — it
                // calls super.getDistances(...), which internally makes a virtual call back to
                // getDistanceRup(...) and hits the overridden (always-throwing) version rather
                // than the base PointSurface implementation. Working around this by going
                // directly to the wrapped uncorrected surface via getUncorrectedSurface().
                
                // TODO ask for opensha fix of the getDistanceRup() method requiring my workaround 
                RuptureSurface surface = rupture.getRuptureSurface();
                double distance;
                if (surface instanceof PointSurface.DistanceCorrectable) {
                    PointSurface uncorrected =
                            ((PointSurface.DistanceCorrectable) surface).getUncorrectedSurface();
                    distance = uncorrected.getDistanceRup(site.getLocation());
                } else {
                    distance = surface.getDistanceRup(site.getLocation());
                }

                gmm.getExceedProbabilities(rupture, hazFunc);

                // the curve is always the conditional exceedance probability
                // regardless of hazard metric; rate vs. probability is
                // handled entirely in the likelihood above.

                double[] hazardVals = new double[hazFunc.size()];
                for (int i = 0; i < hazardVals.length; i++)
                    hazardVals[i] = hazFunc.getY(i);

                RuptureKey key = new RuptureKey(sourceID, ruptureID);

                ruptureLikelihoods.put(key, likelihood);
                ruptureMagnitudes.put(key, magnitude);
                ruptureDistances.put(key, distance);
                ruptureCurves.put(key, new HazardCurve(hazardVals));
            }
        }

        RuptureHazardCurveCollection collection = new RuptureHazardCurveCollection(
                hazardParameters,
                siteKey,
                imKey,
                ruptureLikelihoods,
                ruptureMagnitudes,
                ruptureDistances,
                ruptureCurves,
                erf,
                gmm,
                site,
                filterManager
        );

        collection.freeze();

        return collection;
    }


    // ---------------------------------------------------------------------
    // UTILITIES
    // ---------------------------------------------------------------------

    /** Prints a summary of the collection. */
    public void printSummary() {
        System.out.println("----- RuptureProbabilityCurveCollection Summary -----");
        System.out.printf("Total ruptures: %d%n", getRuptureCount());
        System.out.println(hazardParameters);
        System.out.println(siteKey);
        System.out.println(imKey);
        if (frozen) {
            System.out.println("Aggregate hazard curve: " + java.util.Arrays.toString(aggregateHazardCurve.getHazard()));
            System.out.println("Direct hazard curve:    " + java.util.Arrays.toString(directHazardCurve.getHazard()));
        }
    }

}