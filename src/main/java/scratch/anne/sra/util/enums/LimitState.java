package scratch.anne.sra.util.enums;

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
}