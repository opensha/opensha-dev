package scratch.anne.sra.domain.hazard;

import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

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
