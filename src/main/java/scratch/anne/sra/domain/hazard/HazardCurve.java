package scratch.anne.sra.domain.hazard;

public final class HazardCurve {

    private final double[] hazard;

    public HazardCurve(double[] hazard) {
        this.hazard = hazard;
    }

    public double[] getHazard() {
        return hazard;
    }
}
