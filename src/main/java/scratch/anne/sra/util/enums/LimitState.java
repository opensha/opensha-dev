package scratch.anne.sra.util.enums;

/**
 * Enum defining fragility limit states, including the relationship between them.
 * Ordered limit states occur in sequence. Competing are mutually exclusive.
 */
public enum LimitState {

    SLIGHT(LimitStateType.ORDERED),
    MODERATE(LimitStateType.ORDERED),
    EXTENSIVE(LimitStateType.ORDERED),
    COMPLETE(LimitStateType.ORDERED),
    
    COLLAPSE(LimitStateType.COMPETING), // separate failure mode

    TOPPLE(LimitStateType.COMPETING); // separate failure mode

    private final LimitStateType limitStateType;

    LimitState(LimitStateType limitStateType) {
        this.limitStateType = limitStateType;
    }

    public LimitStateType getLimitStateType() {
        return limitStateType;
    }

    public boolean isOrdered() {
        return limitStateType == LimitStateType.ORDERED;
    }

    public boolean isCompeting() {
        return limitStateType == LimitStateType.COMPETING;
    }
    
    /**
     * Determines how the limit state participates in fragility calculations.
     */
    public enum LimitStateType {
        ORDERED,     // Part of a cumulative sequence (P(DS ≥ i))
        COMPETING    // Independent failure mode, mutually exclusive
    }
}