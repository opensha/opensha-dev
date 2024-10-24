package scratch.kevin.pointSources;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm.DistanceCorrection2013;
import org.opensha.sha.faultSurface.FiniteApproxPointSurface;
import org.opensha.sha.faultSurface.PointSurface;

public class SupersampledApproxPointSurfaceNshm extends PointSurface {

	private FiniteApproxPointSurface ptSurf;
	private Location centerLoc;
	private Region cell;
	private GriddedRegion supersampledCell;
	private double mag;
	
	private DistanceCorrection2013 distCorr = new DistanceCorrection2013();
	
	// should we use the mean or the median?
	private static final boolean MEDIAN_DIST = false;
	// if true, calculate corrected to each supersampled location and then average
	// if false, calculate raw to each supersampled location, average, then correct that
	private static final boolean AVERAGE_CORRECTED = true;

	public SupersampledApproxPointSurfaceNshm(FiniteApproxPointSurface ptSurf, double mag, double gridSpacing, double superSampleGridSpacing) {
		super(ptSurf.getLocation());
		this.ptSurf = ptSurf;
		this.mag = mag;
		this.centerLoc = ptSurf.getLocation();
		cell = new Region(new Location(centerLoc.lat-0.5*gridSpacing, centerLoc.lon-0.5*gridSpacing),
				new Location(centerLoc.lat+0.5*gridSpacing, centerLoc.lon+0.5*gridSpacing));
		supersampledCell = new GriddedRegion(cell, superSampleGridSpacing, new Location(0.5*superSampleGridSpacing, 0.5*superSampleGridSpacing));
	}

	@Override
	public double getAveDip() {
		return ptSurf.getAveDip();
	}

	@Override
	public double getDepth() {
		return ptSurf.getDepth();
	}

	@Override
	public double getAveRupTopDepth() {
		return ptSurf.getAveRupTopDepth();
	}

	@Override
	public double getDistanceRup(Location siteLoc) {
		return ptSurf.getDistanceRup(getDistanceJB(siteLoc));
	}

	@Override
	public double getDistanceJB(Location siteLoc) {
		double[] dists = new double[supersampledCell.getNodeCount()];
		for (int i=0; i<dists.length; i++)
			dists[i] = LocationUtils.horzDistanceFast(supersampledCell.getLocation(i), siteLoc);
		if (AVERAGE_CORRECTED)
			for (int i=0; i<dists.length; i++)
				dists[i] = distCorr.getCorrectedDistanceJB(mag, ptSurf, dists[i]);
		double avg;
		if (MEDIAN_DIST)
			avg = DataUtils.median(dists);
		else
			avg = StatUtils.mean(dists);
		if (!AVERAGE_CORRECTED)
			// correct the average
			avg = distCorr.getCorrectedDistanceJB(mag, ptSurf, avg);
		return avg;
	}

	@Override
	public double getDistanceSeis(Location siteLoc) {
		return ptSurf.getDistanceSeis(siteLoc);
	}

	@Override
	public double getQuickDistance(Location siteLoc) {
		return ptSurf.getQuickDistance(siteLoc);
	}

	@Override
	public double getDistanceX(Location siteLoc) {
		// this isn't quite right, but they're just using the sign anyway
		return ptSurf.getDistanceX(siteLoc);
	}

	@Override
	public double getAveWidth() {
		return ptSurf.getAveWidth();
	}

}
