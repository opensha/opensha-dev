package scratch.anne.sra.engine.logic_tree;

import java.util.Map;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.engine.per_rupture.RuptureLossResultsCollection;

/**
 * Recomputes expected annualized loss (EAL) for a SINGLE logic-tree branch,
 * given per-rupture rates from that branch and a FIXED, precomputed
 * per-rupture loss lookup (from a frozen RuptureLossResultsCollection built
 * against the true mean ERF).
 *
 * <p>Per-rupture losses come from an expensive GMM/vulnerability
 * convolution and are the same for every branch (confirmed: only the
 * RATES vary per branch, not the losses). This class does none of that
 * convolution work -- it's a cheap dot product:</p>
 *
 * <pre>
 *   branchEAL = sum over ruptures [ branchRate(rupture) * fixedLoss(rupture) ]
 * </pre>
 *
 * <p><b>Thread safety / parallelism:</b> each branch should get its OWN
 * instance of this class. All mutable state (the running sum) is local to
 * the instance -- nothing is shared or mutated across branches. The only
 * shared object is {@code fixedLossLookup}, which is read-only after
 * construction, so many BranchEALCalculator instances can safely run
 * concurrently against the same lookup, one per thread/branch.</p>
 *
 * <p>This deliberately does NOT reuse RuptureLossResultsCollection's
 * accumulate/freeze lifecycle -- that class is single-use (frozen once,
 * never reset) and was built for the one-time expensive convolution, not
 * for cheap repeated recomputation across potentially millions of
 * branches. Reusing it here would require adding an unfreeze/reset path
 * that isn't otherwise needed and would undermine its single-use
 * invariant (accumulate-then-freeze) for its original purpose.</p>
 */
public final class BranchEALCalculator {

    /** Read-only, shared across all branches. Built once from the true mean collection. */
    private final Map<RuptureKey, Double> fixedLossLookup;

    /** Local to this branch only. */
    private double runningTotal = 0d;

    public BranchEALCalculator(Map<RuptureKey, Double> fixedLossLookup) {
        this.fixedLossLookup = fixedLossLookup;
    }

    /**
     * Convenience factory: builds the fixed loss lookup directly from a
     * frozen RuptureLossResultsCollection (the true-mean collection,
     * already run through the expensive GMM/vulnerability convolution).
     *
     * @throws IllegalStateException if the collection is not yet frozen
     */
    public static Map<RuptureKey, Double> buildFixedLossLookup(RuptureLossResultsCollection trueMeanCollection) {
        Map<RuptureKey, Double> lookup = new java.util.HashMap<>();
        for (RuptureLossResultsCollection.RuptureResult r : trueMeanCollection.results()) {
            lookup.put(r.key(), r.ruptureLoss());
        }
        return lookup;
    }

    /**
     * Accumulates one fault-based rupture's contribution to this branch's
     * EAL, given the true-mean source/rupture index it maps to (see
     * TrueMeanERFDemo's fault-mapping loop) and this branch's rate for
     * that rupture.
     *
     * <p>If the given key has no entry in the fixed loss lookup (e.g. the
     * rupture had zero rate / zero loss in the true-mean run and was never
     * recorded), it contributes zero rather than throwing -- a rupture
     * with rate > 0 on this branch but no recorded loss is a real
     * possibility if the true-mean run's mapped rupture happened to have
     * zero loss, not necessarily a bug.</p>
     */
    public void accumulateFault(int trueMeanSourceIndex, int trueMeanRupIndexInSource, double branchRate) {
        accumulate(trueMeanSourceIndex, trueMeanRupIndexInSource, branchRate);
    }

    /**
     * Accumulates one gridded rupture's contribution to this branch's EAL,
     * given the true-mean source/rupture index it maps to (see
     * TrueMeanERFDemo's gridded-mapping loop) and this branch's rate for
     * that rupture.
     */
    public void accumulateGridded(int trueMeanSourceIndex, int trueMeanRupIndexInSource, double branchRate) {
        accumulate(trueMeanSourceIndex, trueMeanRupIndexInSource, branchRate);
    }

    private void accumulate(int trueMeanSourceIndex, int trueMeanRupIndexInSource, double branchRate) {
        if (branchRate == 0d)
            return; // no contribution possible; skip the lookup entirely

        RuptureKey key = new RuptureKey(trueMeanSourceIndex, trueMeanRupIndexInSource);
        Double fixedLoss = fixedLossLookup.get(key);
        if (fixedLoss == null)
            return; // no recorded loss for this rupture in the true-mean run; contributes zero

        runningTotal += branchRate * fixedLoss;
    }

    /** The accumulated EAL for this branch, once all ruptures have been added. */
    public double getEAL() {
        return runningTotal;
    }
}
