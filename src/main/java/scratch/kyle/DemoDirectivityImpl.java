package scratch.kyle;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;
//import org.opensha.commons.param.ParameterList;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.mod.AbstractAttenRelMod;

import com.google.common.base.Preconditions;

public class DemoDirectivityImpl extends AbstractAttenRelMod {

	@Override
	public String getShortName() {
		return "DemoDirectivity";
	}

	@Override
	public String getName() {
		return "Demo Direcivity Model";
	}

	@Override
	public void setIMRParams(ScalarIMR imr) {
		// do nothing, this would be used to override parameter settings in the original GMPE if needed
	}

	@Override
	public double getModMean(ScalarIMR imr) {
		EqkRupture rup = imr.getEqkRupture();
		Site site = imr.getSite();
		RuptureSurface surf = rup.getRuptureSurface();
		Location hypo = rup.getHypocenterLocation();
		Preconditions.checkNotNull(hypo, "No hypocenter");
		
		Location siteLoc = site.getLocation();
		
		// here's rx for the site location relative to the rupture trace:
		double rx = surf.getDistanceX(siteLoc);
		// here's ry for the site location relative to the rupture trace
		double ry = GriddedSurfaceUtils.getDistanceY(surf.getUpperEdge(), siteLoc);
		
		// in ln space
		double origMean = imr.getMean();
		
		// TODO actually modify it
		double mean = origMean;
		// super simple fake directivity model, as a placeholder
		double distToSurf = surf.getDistanceX(siteLoc);
		double length = surf.getAveLength();
		double distToHypo = LocationUtils.linearDistanceFast(hypo, siteLoc);
		if (distToHypo > distToSurf)
			mean *= ((distToHypo-distToSurf)/(0.5*length));
		
		// also in ln space
		return mean;
	}

	@Override
	public double getModStdDev(ScalarIMR imr) {
		double origStdDev = imr.getStdDev();
		
		// TODO modify it?
		return origStdDev;
	}

	@Override
	public ParameterList getModParams() {
		// can add adjustable parameters if needed, let me know
		return null;
	}

}
