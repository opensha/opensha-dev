package scratch.anne.sra.engine.hazardSetup;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.imr.ScalarIMR;

import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Builds an OpenSHA {@link Site} from a {@link SiteKey}, cloning the
 * GMM's own site parameters and setting Vs30.
 *
 * <p>Shared by {@code PortfolioRiskConvolutionCalculator},
 * {@code PortfolioPerRuptureRiskConvolutionCalculator}, and
 * {@code RuptureProbabilityCurveCollection.fromERF} -- previously
 * duplicated verbatim in all three.</p>
 */
public final class SiteBuilder {

    private SiteBuilder() {}

    /**
     * @param siteKey lat/lon/vs30 to build the site from
     * @param gmm     GMM whose site parameters are cloned onto the new site
     */
    public static Site fromSiteKey(SiteKey siteKey, ScalarIMR gmm) {
        Site site = new Site(new Location(siteKey.getLat(), siteKey.getLon()));

        for (Parameter<?> p : gmm.getSiteParams())
            site.addParameter((Parameter<?>) p.clone());

        @SuppressWarnings("unchecked")
        Parameter<Double> vs30 = (Parameter<Double>) site.getParameter("Vs30");
        vs30.setValue(siteKey.getVs30());

        return site;
    }
}
