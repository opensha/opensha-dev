package scratch.anne.sra.domain.hazard;

/**
 * A bare hazard curve: rate or probability of exceedance at each of a
 * response model's IM values (see {@link scratch.anne.sra.util.AssetKeys.ImKey}).
 * Carries no site/IM identity of its own.
 *
 * <p>For a hazard curve tagged with the {@code SiteKey}/{@code ImKey} it
 * belongs to, see {@link AssetHazardCurve}.</p>
 */
public final class HazardCurve {

    private final double[] hazard;

    public HazardCurve(double[] hazard) {
        this.hazard = hazard;
    }

    public double[] getHazard() {
        return hazard;
    }
}
