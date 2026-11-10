package scratch.anne.sra.engine.perRupture;

import java.util.concurrent.atomic.DoubleAdder;

/**
 * Thread-safe accumulator for rupture-level losses, indexed by a dense
 * {@code flatIndex} (0..numRuptures-1) rather than a {@code RuptureKey}
 * object.
 *
 * <p>Backed by a plain {@code DoubleAdder[]}, pre-populated once (single-
 * threaded, in the constructor) with one {@code DoubleAdder} per
 * flatIndex. This avoids the per-entry overhead of a
 * {@code ConcurrentHashMap<RuptureKey, DoubleAdder>} -- no hash table,
 * no {@code RuptureKey} object allocations, no boxed entry nodes -- while
 * keeping identical thread-safety: concurrent calls to
 * {@link #accumulateRuptureLoss(int, double)} for the same flatIndex are
 * still backed by the same {@code DoubleAdder} instance's own internal
 * striping, exactly as before. The dense flatIndex itself is assigned
 * once, up front, single-threaded, by whoever builds the index mapping
 * (see {@link RuptureLossResultsCollection}), so there is no risk of two
 * different ruptures colliding on the same slot. "flatIndex" is named to
 * avoid confusion with the unrelated, source-scoped fltSysRupIndex
 * concept also tracked by {@link RuptureLossResultsCollection}.</p>
 */
public final class RuptureLossAccumulator {

    private boolean frozen = false;

    private final DoubleAdder[] lossByRupture;

    /**
     * Pre-populates an empty DoubleAdder for every flatIndex in
     * {@code [0, numRuptures)}.
     */
    public RuptureLossAccumulator(int numRuptures) {
        this.lossByRupture = new DoubleAdder[numRuptures];
        for (int i = 0; i < numRuptures; i++)
            lossByRupture[i] = new DoubleAdder();
    }

    /**
     * Adds loss contribution for a rupture.
     *
     * <p>This is the ONLY mutation method intended for hot-loop usage.</p>
     */
    public void accumulateRuptureLoss(int flatIndex, double loss) {
        if (frozen)
            throw new IllegalStateException(
                "RuptureLossAccumulator is frozen");

        checkIndex(flatIndex);
        lossByRupture[flatIndex].add(loss);
    }

    /**
     * Gets accumulated loss for a rupture.
     * Returns 0 if no loss was recorded.
     */
    public double getLoss(int flatIndex) {
        checkIndex(flatIndex);
        return lossByRupture[flatIndex].sum();
    }

    /**
     * Optional: number of ruptures tracked by this accumulator.
     */
    public int size() {
        return lossByRupture.length;
    }

    public void freeze() {
        frozen = true;
    }

    private void checkIndex(int flatIndex) {
        if (flatIndex < 0 || flatIndex >= lossByRupture.length)
            throw new IllegalArgumentException("Unknown flatIndex: " + flatIndex);
    }
}
