package scratch.anne.sra.domain.structural_response.fragilities.definitions;

import scratch.anne.sra.util.enums.LimitState;

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