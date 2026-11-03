package scratch.anne.sra.domain.hazard;

import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * A {@link HazardCurve} tagged with the site and intensity measure it
 * belongs to, plus the linear IM values the hazard was evaluated at.
 * Returned when exporting or inspecting a single asset's hazard, as
 * opposed to the bare per-(site, IM) storage in {@link HazardCurveCollection}.
 */
public final class AssetHazardCurve {
    public final SiteKey siteKey;
    public final ImKey imKey;
    public final double[] imls;
    public final double[] hazard;

    public AssetHazardCurve(SiteKey siteKey,
                           ImKey imKey,
                           double[] imls,
                           double[] hazard) {
        this.siteKey = siteKey;
        this.imKey = imKey;
        this.imls = imls;
        this.hazard = hazard;
    }
}
