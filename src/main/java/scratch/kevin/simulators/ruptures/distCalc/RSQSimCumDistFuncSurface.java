package scratch.kevin.simulators.ruptures.distCalc;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;

import scratch.kevin.simulators.ruptures.distCalc.SimRuptureDistCalcUtils.DistanceType;
import scratch.kevin.simulators.ruptures.distCalc.SimRuptureDistCalcUtils.LocationElementDistanceCache;
import scratch.kevin.simulators.ruptures.distCalc.SimRuptureDistCalcUtils.Scalar;

public class RSQSimCumDistFuncSurface extends CompoundSurface {
	
	private RSQSimEvent event;
	private Scalar scalar;
	private double threshold;

	public RSQSimCumDistFuncSurface(RSQSimEvent event, Scalar scalar, double threshold,
			List<? extends FaultSection> subSects) {
		super(buildSurfs(subSects));
		this.event = event;
		this.scalar = scalar;
		this.threshold = threshold;
	}
	
	private static List<RuptureSurface> buildSurfs(List<? extends FaultSection> subSects) {
		double gridSpacing = 1d;

		List<RuptureSurface> rupSurfs = new ArrayList<>();
		for (FaultSection sect : subSects)
			rupSurfs.add(sect.getFaultSurface(gridSpacing, false, false));
		return rupSurfs;
	}

	@Override
	public CompoundSurfaceDistances calcDistances(Location loc) {
		CompoundSurfaceDistances dists = super.calcDistances(loc);
		
		LocationElementDistanceCache siteLocDistCache = SimRuptureDistCalcUtils.buildSiteLocDistCache(loc);
		
		DiscretizedFunc rJBFunc = SimRuptureDistCalcUtils.calcDistScalarFunc(event, loc, siteLocDistCache,
				DistanceType.R_JB, scalar);
		double distanceJB = Double.NaN;
		double targetVal = threshold * rJBFunc.getY(rJBFunc.size()-1);
		for (Point2D pt : rJBFunc) {
			if (pt.getY() >= targetVal) {
				distanceJB = pt.getX();
				break;
			}
		}
		
		DiscretizedFunc rRupFunc = SimRuptureDistCalcUtils.calcDistScalarFunc(event, loc, siteLocDistCache,
				DistanceType.R_RUP, scalar);
		double distanceRup = Double.NaN;
		targetVal = threshold * rRupFunc.getY(rRupFunc.size()-1);
		for (Point2D pt : rRupFunc) {
			if (pt.getY() >= targetVal) {
				distanceRup = pt.getX();
				break;
			}
		}
		
		dists = new CompoundSurfaceDistances(distanceRup, distanceJB, dists.getDistanceSeis(), dists.distXIndex);
		
		return dists;
	}

}
