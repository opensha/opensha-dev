package scratch.anne.sra.engine.per_rupture;

import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;

import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.domain.hazard.HazardParameters.HazardMetricType;

/**
 * Per-rupture expected-loss results for an entire ERF.
 *
 * <p>Internally indexed by a dense {@code int flatIndex} (0..numRuptures-1)
 * rather than by {@code RuptureKey} objects, backed by plain primitive
 * arrays ({@code double[]} for likelihood/magnitude, {@code int[]} for
 * fltSysRupIndex) plus an array-backed {@link RuptureLossAccumulator} --
 * NOT a {@code Map<RuptureKey, ...>} per field as in earlier versions of
 * this class.</p>
 *
 * <p>Rationale: for an ERF with a large total rupture count (e.g. a
 * gridded+fault-system true-mean solution with 1M+ ruptures), eagerly
 * populating one {@code Map<RuptureKey, Double/Integer>} entry per
 * rupture, across several parallel maps, costs far more than the raw
 * data: every {@code RuptureKey} is its own object (duplicated across
 * each map that stores it), every {@code Double}/{@code Integer} is
 * boxed, and every hash-map entry carries its own node overhead (hash,
 * next-pointer, key/value references) on top of that. A
 * {@code ConcurrentHashMap<RuptureKey, DoubleAdder>} accumulator under
 * heavy parallel write contention (dozens of threads) compounds this
 * further. Flat arrays eliminate essentially all of that: a
 * {@code double[]} of size N is N*8 bytes, full stop, with zero
 * per-entry object/hash overhead -- the practical difference between
 * fitting comfortably in memory and an OutOfMemoryError on a
 * million-plus-rupture ERF.</p>
 *
 * <p>The {@code RuptureKey}-based public API is unchanged -- every
 * public method still takes/returns {@code RuptureKey}, so no caller of
 * this class needs to change. The {@code RuptureKey -> flatIndex}
 * translation happens once, internally, via {@link #flatIndexOf(RuptureKey)}
 * (a direct two-level array lookup, not a hash lookup -- sourceID and
 * ruptureID are already small, dense, ERF-assigned integers). "flatIndex"
 * is deliberately named to avoid any visual/semantic confusion with
 * {@code fltSysRupIndex} (a fault-system rupture index, an entirely
 * different, source-scoped concept) -- flatIndex is purely an internal
 * storage slot and is never exposed outside this class.</p>
 */
public final class RuptureLossResultsCollection {

    private final HazardParameters hazardParameters;

    /**
     * sourceID -> ruptureID -> dense flatIndex. Built once, up front,
     * single-threaded, in {@link #fromERF}; read-only afterward, so
     * concurrent lookups from the hot loop are safe with no locking.
     */
    private final int[][] flatIndexBySourceAndRupture;

    private final RuptureKey[] keyByFlatIndex;
    private final double[] likelihoods;
    private final double[] magnitudes;
    /**
     * Fault-system rupture index per flatIndex, or {@link #NO_FLT_SYS_RUP_INDEX}
     * when not applicable (gridded/background sources, or an ERF that
     * isn't a {@link BaseFaultSystemSolutionERF} at all). A primitive
     * {@code int[]} has no "absent" state the way a boxed {@code Integer}
     * / map-entry-absence does, so this sentinel plays that role instead
     * -- never a real fault-system rupture index, which is always {@code >= 0}.
     */
    private static final int NO_FLT_SYS_RUP_INDEX = -1;
    private final int[] fltSysRupIndices;

    private final RuptureLossAccumulator accumulator;

    private double totalExpectedLoss;
    private double[] lossContributions;
    private double[] relativeContributions;

    private boolean frozen;

    private RuptureLossResultsCollection(
            HazardParameters hazardParameters,
            int[][] flatIndexBySourceAndRupture,
            RuptureKey[] keyByFlatIndex,
            double[] likelihoods,
            double[] magnitudes,
            int[] fltSysRupIndices
    ) {
        this.hazardParameters = hazardParameters;

        if (hazardParameters.getHazardMetricType() != HazardMetricType.EXCEEDANCE) {
            throw new IllegalArgumentException(
                "RuptureResultsCollection requires EXCEEDANCE hazard metric type, but got: "
                + hazardParameters.getHazardMetricType()
            );
        }

        this.flatIndexBySourceAndRupture = flatIndexBySourceAndRupture;
        this.keyByFlatIndex = keyByFlatIndex;
        this.likelihoods = likelihoods;
        this.magnitudes = magnitudes;
        this.fltSysRupIndices = fltSysRupIndices;

        this.accumulator = new RuptureLossAccumulator(keyByFlatIndex.length);

        this.totalExpectedLoss = 0.0;
        this.lossContributions = new double[keyByFlatIndex.length];
        this.relativeContributions = new double[keyByFlatIndex.length];
    }


    public void accumulateRuptureLoss(RuptureKey key, double loss) {
        if (frozen)
            throw new IllegalStateException(
                "RuptureResultsCollection is frozen");
        accumulator.accumulateRuptureLoss(flatIndexOf(key), loss);
    }


    public double getLikelihood(RuptureKey key) {
        return likelihoods[flatIndexOf(key)];
    }


    public double getMagnitude(RuptureKey key) {
        return magnitudes[flatIndexOf(key)];
    }


    /**
     * Returns the fault-system rupture index for this rupture's source, or
     * {@code null} if
     *    a) the rupture is from gridded seismicity, or
     *    b) these results were not built from a {@link BaseFaultSystemSolutionERF}.
     */
    public Integer getFltSysRupIndex(RuptureKey key) {
        int value = fltSysRupIndices[flatIndexOf(key)];
        return value == NO_FLT_SYS_RUP_INDEX ? null : value;
    }


    public double getRuptureLoss(RuptureKey key) {
        return accumulator.getLoss(flatIndexOf(key));
    }


    public double getRuptureContribution(RuptureKey key) {

        if (!frozen)
            throw new IllegalStateException(
                "Rupture contributions not yet computed");

        return lossContributions[flatIndexOf(key)];
    }


    public double getTotalExpectedLoss() {
        if (!frozen)
            throw new IllegalStateException(
                "Total expected loss not yet computed");
        return totalExpectedLoss;
    }

    public int getRuptureCount() {
        return keyByFlatIndex.length;
    }


    public HazardParameters getHazardParameters() {
        return hazardParameters;
    }


    public RuptureResult getTopContributor() {

        if (!frozen)
            throw new IllegalStateException(
                    "Rupture contributions not yet computed");

        int bestIndex = -1;
        double bestValue = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < lossContributions.length; i++) {
            if (lossContributions[i] > bestValue) {
                bestValue = lossContributions[i];
                bestIndex = i;
            }
        }
        if (bestIndex < 0)
            throw new java.util.NoSuchElementException("No ruptures in this collection");

        return get(keyByFlatIndex[bestIndex]);
    }

    public RuptureKey getTopContributorKey() {

        return getTopContributor().key();
    }


    public List<RuptureResult> getTopContributors(int n) {

        if (!frozen)
            throw new IllegalStateException(
                    "Rupture contributions not yet computed");

        Integer[] order = new Integer[lossContributions.length];
        for (int i = 0; i < order.length; i++)
            order[i] = i;
        java.util.Arrays.sort(order, (a, b) -> Double.compare(lossContributions[b], lossContributions[a]));

        List<RuptureResult> results = new ArrayList<>(Math.min(n, order.length));
        for (int i = 0; i < Math.min(n, order.length); i++)
            results.add(get(keyByFlatIndex[order[i]]));
        return results;
    }


    public List<RuptureKey> getTopContributorsKeys(int n) {
        return getTopContributors(n).stream()
                .map(RuptureResult::key)
                .toList();
    }


    /**
     * Immutable result for export (CSV / Parquet).
     */
    public record RuptureResult(
            RuptureKey key,
            double likelihood,
            double magnitude,
            Integer fltSysRupIndex,
            double ruptureLoss,
            double lossContribution,
            double relativeContribution
    ) {}

    public RuptureResult get(RuptureKey key) {
        int flatIndex = flatIndexOf(key);
        int fltSysRupIndex = fltSysRupIndices[flatIndex];
        return new RuptureResult(
                key,
                likelihoods[flatIndex],
                magnitudes[flatIndex],
                fltSysRupIndex == NO_FLT_SYS_RUP_INDEX ? null : fltSysRupIndex,
                accumulator.getLoss(flatIndex),
                lossContributions[flatIndex],
                relativeContributions[flatIndex]
        );
    }

    public Iterable<RuptureResult> results() {
        List<RuptureResult> all = new ArrayList<>(keyByFlatIndex.length);
        for (RuptureKey key : keyByFlatIndex)
            all.add(get(key));
        return all;
    }


    /**
     * Builds a rupture results collection from an already-initialized ERF.
     *
     * <p>The supplied ERF must have already had
     * {@code updateForecast()} called.</p>
     *
     * <p>If {@code erf} is a {@link BaseFaultSystemSolutionERF}, the
     * fault-system rupture index for each rupture's source is also captured
     * (via {@link BaseFaultSystemSolutionERF#getFltSysRupIndexForSource(int)});
     * otherwise every rupture's fault-system rupture index is absent
     * ({@code null} via {@link #getFltSysRupIndex}).</p>
     */
    public static RuptureLossResultsCollection fromERF(
            AbstractERF erf,
            HazardParameters hazardParameters
    ) {

        BaseFaultSystemSolutionERF fssErf =
                (erf instanceof BaseFaultSystemSolutionERF b) ? b : null;
        // fault-system sources occupy indices [0, numFaultSystemSources) when background
        // seismicity (IncludeBackgroundOption.INCLUDE or EXCLUDE) is appended after them;
        // under ONLY, every source is a background source and none of them map back to a
        // fault-system rupture index
        int numFaultSystemSources = (fssErf != null) ? fssErf.getNumFaultSystemSources() : 0;

        int numSources = erf.getNumSources();

        // ------------------------------------------------------------
        // Pass 1: count total ruptures across all sources, so every
        // array below can be allocated exactly once, at its final size
        // -- no resizing, no transient double-capacity copies.
        // ------------------------------------------------------------
        int totalRuptures = 0;
        int[] numRupturesBySource = new int[numSources];
        for (int sourceID = 0; sourceID < numSources; sourceID++) {
            int numRuptures = erf.getNumRuptures(sourceID);
            numRupturesBySource[sourceID] = numRuptures;
            totalRuptures += numRuptures;
        }

        int[][] flatIndexBySourceAndRupture = new int[numSources][];
        RuptureKey[] keyByFlatIndex = new RuptureKey[totalRuptures];
        double[] likelihoods = new double[totalRuptures];
        double[] magnitudes = new double[totalRuptures];
        int[] fltSysRupIndices = new int[totalRuptures];

        // ------------------------------------------------------------
        // Pass 2: populate every array at the dense flatIndex assigned
        // by iterating sources/ruptures in the same fixed order used
        // to size them above.
        // ------------------------------------------------------------
        int flatIndex = 0;
        for (int sourceID = 0; sourceID < numSources; sourceID++) {

            int numRuptures = numRupturesBySource[sourceID];
            flatIndexBySourceAndRupture[sourceID] = new int[numRuptures];

            // fault-system rupture index is a property of the source, not of
            // the individual rupture, so compute it once per source; only
            // valid for sources that are actually fault-system sources
            int fltSysRupIndexForSource = (fssErf != null && sourceID < numFaultSystemSources)
                    ? fssErf.getFltSysRupIndexForSource(sourceID)
                    : NO_FLT_SYS_RUP_INDEX;

            for (int ruptureID = 0; ruptureID < numRuptures; ruptureID++) {

                ProbEqkRupture rupture = erf.getSource(sourceID).getRupture(ruptureID);

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

                flatIndexBySourceAndRupture[sourceID][ruptureID] = flatIndex;
                keyByFlatIndex[flatIndex] = new RuptureKey(sourceID, ruptureID);
                likelihoods[flatIndex] = likelihood;
                magnitudes[flatIndex] = rupture.getMag();
                fltSysRupIndices[flatIndex] = fltSysRupIndexForSource;

                flatIndex++;
            }
        }

        return new RuptureLossResultsCollection(
                hazardParameters,
                flatIndexBySourceAndRupture,
                keyByFlatIndex,
                likelihoods,
                magnitudes,
                fltSysRupIndices
        );
    }

    public void freeze() {
        if (frozen)
            return;

        accumulator.freeze();

        frozen = true;

        totalExpectedLoss = computeLossOverAllRuptures();
    }

    private double computeLossOverAllRuptures() {

        if (!frozen)
            throw new IllegalStateException(
                "Total loss can only be calculated after all ruptures are accumulated and frozen");

        // for both rate and probability, total loss is the sum of (likelihood * conditional loss)
        // the difference is handled in the rupture likelihood, not in the combination equation
        //TODO confirm whether prob_exc is appropriate for per-rupture

        double total = 0.0;

        // first pass: absolute contributions
        for (int i = 0; i < keyByFlatIndex.length; i++) {
            double contribution = likelihoods[i] * accumulator.getLoss(i);
            lossContributions[i] = contribution;
            total += contribution;
        }

        // second pass: relative contributions
        if (total > 0.0) {
            for (int i = 0; i < keyByFlatIndex.length; i++)
                relativeContributions[i] = lossContributions[i] / total;
        } else {
            // edge case: no loss -- arrays already default to 0.0
        }

        return total;
    }


    /**
     * Translates a RuptureKey into its dense flatIndex via a direct
     * two-level array lookup (NOT a hash lookup) -- sourceId/ruptureId
     * are already small, dense, ERF-assigned integers.
     */
    private int flatIndexOf(RuptureKey key) {
        int sourceId = key.sourceId();
        if (sourceId < 0 || sourceId >= flatIndexBySourceAndRupture.length)
            throw new IllegalArgumentException("Unknown rupture key: " + key);
        int[] ruptureIds = flatIndexBySourceAndRupture[sourceId];
        int ruptureId = key.ruptureId();
        if (ruptureId < 0 || ruptureId >= ruptureIds.length)
            throw new IllegalArgumentException("Unknown rupture key: " + key);
        return ruptureIds[ruptureId];
    }


    // ---------------------------------------------------------------------
    // UTILITIES
    // ---------------------------------------------------------------------

    /** Prints a summary of portfolio losses. */
    public void printSummary() {
        System.out.println("----- ELossPortfolio Summary -----");
        System.out.printf("Total ruptures: %d%n", getRuptureCount());
        System.out.printf("Total expected loss: %.2e%n", getTotalExpectedLoss());
    }

}
