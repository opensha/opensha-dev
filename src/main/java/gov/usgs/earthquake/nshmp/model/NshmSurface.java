package gov.usgs.earthquake.nshmp.model;

import java.util.ListIterator;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.CacheEnabledSurface;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;

/**
 * Rupture surface implementation for USGS NSHMs. Most methods throw an
 * UnsupportedOperationException except those required for hazard calculations
 * with current GMMs (dip, width, rRup, rJb, rX, zTor)
 *
 * @author U.S. Geological Survey
 */
public class NshmSurface implements CacheEnabledSurface {

	private final gov.usgs.earthquake.nshmp.fault.surface.RuptureSurface delegate;

	// distance metrics for reference site; this should
	// work for single threaded calculations
	private Location location;
	private Distance distance;

	public NshmSurface(gov.usgs.earthquake.nshmp.fault.surface.RuptureSurface delegate) {
		this.delegate = delegate;
	}

	// return nshmp-haz rupture centroid as OpenSHA location for
	// use in computing min distance to a fault system subsection
	Location centroid() {
		return NshmUtil.toOpenShaLocation(delegate.centroid());
	}

	// OpenSHA RupureSurface interface methods

	// @formatter:off
	@Override public double getAveDip() { return delegate.dip(); }
	@Override public double getAveWidth() { return delegate.width(); }
	@Override public double getAveRupTopDepth() { return delegate.depth(); }
	@Override public double getArea() { return delegate.area(); }

	@Override
	public synchronized double getDistanceRup(Location location) {
		if (location != this.location) {
			setDistances(location);
		}
		return distance.rRup;
	}

	@Override
	public synchronized double getDistanceJB(Location location) {
		if (location != this.location) {
			setDistances(location);
		}
		return distance.rJB;
	}

	@Override
	public synchronized double getDistanceX(Location location) {
		if (location != this.location) {
			setDistances(location);
		}
		return distance.rX;
	}

	@Override
	public synchronized double getDistanceSeis(Location location) {
		// distanceSeis isn't used by any modern GMM so we're just returning rRup
		if (location != this.location) {
			setDistances(location);
		}
		return distance.rRup;
	}

	private void setDistances(Location location) {
		this.distance = delegate.distanceTo(NshmUtil.fromOpenShaLocation(location));
		this.location = location;
	}

	// Needed by Compound surface initialization

	@Override
	public Location getFirstLocOnUpperEdge() {
		// this will only be asked for by OpenSHA CompoundSurface
		return NshmUtil.toOpenShaLocation(
				((gov.usgs.earthquake.nshmp.fault.surface.GriddedSurface) delegate)
				.getFirstLocOnUpperEdge());
	}

	@Override
	public Location getLastLocOnUpperEdge() {
		// this will only be asked for by OpenSHA CompoundSurface
		return NshmUtil.toOpenShaLocation(
				((gov.usgs.earthquake.nshmp.fault.surface.GriddedSurface) delegate)
				.getLastLocOnUpperEdge());
	}
	
	// Caching

	@Override
	public SurfaceDistances calcDistances(Location location) {
		Distance distance = delegate.distanceTo(NshmUtil.fromOpenShaLocation(location));
		return new SurfaceDistances(distance.rRup, distance.rJB, distance.rRup);
	}

	@Override
	public double calcQuickDistance(Location location) {
		return LocationUtils.horzDistanceFast(centroid(), location);
	}

	@Override
	public double calcDistanceX(Location location) {
		Distance distance = delegate.distanceTo(NshmUtil.fromOpenShaLocation(location));
		return distance.rX;
	}

	@Override
	public synchronized void clearCache() {
		this.location = null;
		this.distance = null;
	}

	// Unnecessary methods for hazard calculations

	@Override public ListIterator<Location> getLocationsIterator() { throw new UnsupportedOperationException(); }
	@Override public LocationList getEvenlyDiscritizedPerimeter() { throw new UnsupportedOperationException(); }
	@Override public LocationList getPerimeter() { throw new UnsupportedOperationException(); }
	@Override public boolean isPointSurface() { throw new UnsupportedOperationException(); }
	@Override public double getAveStrike() { throw new UnsupportedOperationException(); }
	@Override public double getAveLength() { throw new UnsupportedOperationException(); }
	@Override public double getAreaInsideRegion(Region region) { throw new UnsupportedOperationException(); }
	@Override public LocationList getEvenlyDiscritizedListOfLocsOnSurface() { throw new UnsupportedOperationException(); }
	@Override public FaultTrace getEvenlyDiscritizedUpperEdge() { throw new UnsupportedOperationException(); }
	@Override public LocationList getEvenlyDiscritizedLowerEdge() { throw new UnsupportedOperationException(); }
	@Override public double getAveGridSpacing() { throw new UnsupportedOperationException(); }
	@Override public double getQuickDistance(Location siteLoc) { throw new UnsupportedOperationException(); }
	@Override public double getAveDipDirection() { throw new UnsupportedOperationException(); }
	@Override public FaultTrace getUpperEdge() { throw new UnsupportedOperationException(); }
	@Override public double getFractionOfSurfaceInRegion(Region region) { throw new UnsupportedOperationException(); }
	@Override public String getInfo() { throw new UnsupportedOperationException(); }
	@Override public double getMinDistance(RuptureSurface surface) { throw new UnsupportedOperationException(); }
	@Override public RuptureSurface getMoved(LocationVector v) { throw new UnsupportedOperationException(); }
	@Override public RuptureSurface copyShallow() { throw new UnsupportedOperationException(); }
}
