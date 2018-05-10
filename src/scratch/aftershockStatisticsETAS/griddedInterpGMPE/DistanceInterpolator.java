package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.ScalarIMR;

import com.google.common.base.Preconditions;

public class DistanceInterpolator extends AbstractGMPEInterpolation.LogSpacedDouble {

	public DistanceInterpolator(boolean includeZero, double minNonZeroDist, double maxDist, int numBins) {
		super("Horizontal Distance", minNonZeroDist, maxDist, numBins, includeZero);
	}

	@Override
	public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {
		// set site location for current rupture
		double dist = getValue(index);
		EqkRupture rup = gmpe.getEqkRupture();
		Preconditions.checkState(rup.getRuptureSurface() instanceof PointSurface);
		Location rupLoc = ((PointSurface)rup.getRuptureSurface()).getLocation();
		Location siteLoc = LocationUtils.location(rupLoc, 0, dist);
		gmpe.setSiteLocation(siteLoc);
	}

//	@Override
//	public Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index) {
//		return ruptures;
//	}

	@Override
	public Double detectCurrentVal(ScalarIMR gmpe, Site site) {
		EqkRupture rup = gmpe.getEqkRupture();
		Preconditions.checkState(rup.getRuptureSurface() instanceof PointSurface);
		Location rupLoc = ((PointSurface)rup.getRuptureSurface()).getLocation();
		return LocationUtils.horzDistanceFast(gmpe.getSite().getLocation(), rupLoc);
	}

}
