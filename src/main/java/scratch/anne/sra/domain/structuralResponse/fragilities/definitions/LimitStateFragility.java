package scratch.anne.sra.domain.structuralResponse.fragilities.definitions;

import scratch.anne.sra.util.enums.LimitState;

/**
 * Pairs a {@link LimitState} (e.g. SLIGHT, COLLAPSE) with the
 * {@link FragilityDefinition} describing that limit state's exceedance
 * curve. One {@link FragilityModel} holds a list of these, one per
 * limit state it defines.
 */
public class LimitStateFragility {

    private final LimitState limitState;
    private final FragilityDefinition definition;

    public LimitStateFragility(LimitState limitState,
                                FragilityDefinition definition) {

        if (limitState == null)
            throw new IllegalArgumentException("Limit state required");

        if (definition == null)
            throw new IllegalArgumentException("Definition cannot be null");

        this.limitState = limitState;
        this.definition = definition;
    }

    public LimitState getLimitState() {
        return limitState;
    }

    public FragilityDefinition getDefinition() {
        return definition;
    }

    public boolean isOrdered() {
        return limitState.isOrdered();
    }

    public boolean isCompeting() {
        return limitState.isCompeting();
    }
}