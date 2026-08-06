package org.opensha.nshmp.shaded.model;

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
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrection;
import org.opensha.sha.util.TectonicRegionType;

import org.opensha.nshmp.shaded.fault.surface.NshmpDefaultGriddedSurface;
import org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface;
import org.opensha.nshmp.shaded.model.NshmpGridSource.PointSurface;
import org.opensha.nshmp.shaded.model.NshmpGridSourceFinite.FiniteSurface;

/**
 * NshmpRupture surface implementation for USGS NSHMs. Most methods throw an
 * UnsupportedOperationException except those required for hazard calculations
 * with current GMMs (dip, width, rRup, rJb, rX, zTor)
 *
 * @author U.S. Geological Survey
 */
public class NshmSurface implements CacheEnabledSurface {

	private final org.opensha.nshmp.shaded.fault.surface.NshmpRuptureSurface delegate;

	// distance metrics for reference site; this should
	// work for single threaded calculations
	private Location location;
	private NshmpDistance distance;

	public NshmSurface(org.opensha.nshmp.shaded.fault.surface.NshmpRuptureSurface delegate) {
		this.delegate = delegate;
	}
	
	/**
	 * This creates a point surface that will work with existing OpenSHA point-source optimizations
	 * @param delegate
	 * @return
	 */
	public static org.opensha.sha.faultSurface.PointSurface buildPointSurface(
			gov.usgs.earthquake.nshmp.fault.surface.RuptureSurface delegate) {
		// this is the point surface
		double len = 0d;
		try {
			len = delegate.length();
		} catch (Exception e) {}
		org.opensha.sha.faultSurface.PointSurface surf = new org.opensha.sha.faultSurface.PointSurface(
				NshmUtil.toOpenShaLocation(delegate.centroid()), delegate.dip(), delegate.depth(),
				delegate.depth() + delegate.width()*Math.sin(Math.toRadians(delegate.dip())), len);
		return new org.opensha.sha.faultSurface.PointSurface.DistanceCorrecting(
				surf, new DelegatePointSourceCorrection(delegate), null, Double.NaN);
	}
	
	/**
	 * Delegate point source correction that passes through to NSHMP-haz. This is required for point source
	 * optimizations to work with wrapped point sources. As part of that, equals/hashCode have to be constant
	 * for all delegate corrections.
	 */
	private static class DelegatePointSourceCorrection implements PointSourceDistanceCorrection.Single {
		
		private gov.usgs.earthquake.nshmp.fault.surface.RuptureSurface delegate;

		private DelegatePointSourceCorrection(gov.usgs.earthquake.nshmp.fault.surface.RuptureSurface delegate) {
			this.delegate = delegate;
		}

		@Override
		public SurfaceDistances getCorrectedDistance(Location location, org.opensha.sha.faultSurface.PointSurface surf,
				TectonicRegionType trt, double mag, double horzDist) {
			Distance distance = delegate.distanceTo(NshmUtil.fromOpenShaLocation(location));
			return new SurfaceDistances.Precomputed(location, distance.rRup, distance.rJB, distance.rX);
		}
		
		private static final int hashCode = DelegatePointSourceCorrection.class.hashCode();

		@Override
		public int hashCode() {
			return hashCode;
		}

		@Override
		public boolean equals(Object obj) {
//			return obj instanceof DelegatePointSourceCorrection;
			if (!(obj instanceof DelegatePointSourceCorrection))
				return false;
			DelegatePointSourceCorrection other = (DelegatePointSourceCorrection)obj;
//			System.out.println("Delegate classes:\t"+delegate.getClass()+"\t"+other.delegate.getClass());
			return other.delegate.getClass().equals(delegate.getClass());
		}
		
	}

	// return nshmp-haz rupture centroid as OpenSHA location for
	// use in computing min distance to a fault system subsection
	public Location centroid() {
		return NshmUtil.toOpenShaLocation(delegate.centroid());
	}

	// OpenSHA RupureSurface interface methods

	// @formatter:off
	@Override public double getAveDip() { return delegate.dip(); }
	@Override public double getAveWidth() { return delegate.width(); }
	@Override public double getAveHorizontalWidth() { return delegate.width()*Math.cos(Math.toRadians(delegate.dip())); }
	@Override public double getArea() { return delegate.area(); }


	@Override public double getAveRupTopDepth() {
		if (delegate instanceof NshmpDefaultGriddedSurface) {
			return ((NshmpDefaultGriddedSurface) delegate).get(0, 0).depth;
		}
		return delegate.depth();
	}

	@Override public double getAveRupBottomDepth() {
		if (delegate instanceof NshmpDefaultGriddedSurface) {
			return ((NshmpDefaultGriddedSurface) delegate).get(((NshmpDefaultGriddedSurface) delegate).getNumRows()-1, 0).depth;
		}
		return delegate.depth() + delegate.width()*Math.sin(Math.toRadians(delegate.dip()));
	}

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
	public synchronized SurfaceDistances getDistances(Location location) {
		if (location != this.location) {
			setDistances(location);
		}
		return new SurfaceDistances.Precomputed(location, distance.rRup, distance.rJB, distance.rX);
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
				((org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface) delegate)
				.getFirstLocOnUpperEdge());
	}

	@Override
	public Location getLastLocOnUpperEdge() {
		// this will only be asked for by OpenSHA CompoundSurface
		return NshmUtil.toOpenShaLocation(
				((org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface) delegate)
				.getLastLocOnUpperEdge());
	}

	@Override
	public Location getFirstLocOnLowerEdge() {
		// this will only be asked for by OpenSHA CompoundSurface
		org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface gridDelegate = (org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface) delegate;
		return NshmUtil.toOpenShaLocation(gridDelegate.getLocation(gridDelegate.getNumRows()-1, 0));
	}

	@Override
	public Location getLastLocOnLowerEdge() {
		// this will only be asked for by OpenSHA CompoundSurface
		org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface gridDelegate = (org.opensha.nshmp.shaded.fault.surface.NshmpGriddedSurface) delegate;
		return NshmUtil.toOpenShaLocation(gridDelegate.getLocation(gridDelegate.getNumRows()-1, gridDelegate.getNumCols()-1));
	}

	// Caching

	@Override
	public SurfaceDistances calcDistances(Location location) {
		NshmpDistance distance = delegate.distanceTo(NshmUtil.fromOpenShaLocation(location));
		return new SurfaceDistances.Precomputed(location, distance.rRup, distance.rJB, distance.rX);
	}

	@Override
	public double calcQuickDistance(Location location) {
		return LocationUtils.horzDistanceFast(centroid(), location);
	}

	@Override
	public synchronized void clearCache() {
		this.location = null;
		this.distance = null;
	}

	// ERF Calculations

	@Override
	public LocationList getEvenlyDiscritizedListOfLocsOnSurface() {
		if (delegate instanceof FiniteSurface) {
			return LocationList.of(NshmUtil.toOpenShaLocation(((FiniteSurface) delegate).loc));
		} else if (delegate instanceof PointSurface) {
			return LocationList.of(NshmUtil.toOpenShaLocation(((PointSurface) delegate).loc));
		}
		return NshmUtil.toOpenShaLocationList(
				((NshmpGriddedSurface) delegate).getEvenlyDiscritizedListOfLocsOnSurface());
	}

	@Override
	public double getAveLength() {
		return delegate.length();
	}

	@Override public double getQuickDistance(Location siteLoc) {
		return calcQuickDistance(siteLoc);
	}

	// Unnecessary methods for hazard calculations

	@Override public ListIterator<Location> getLocationsIterator() { throw new UnsupportedOperationException(); }
	@Override public LocationList getEvenlyDiscritizedPerimeter() { throw new UnsupportedOperationException(); }
	@Override public LocationList getPerimeter() { throw new UnsupportedOperationException(); }
	@Override public boolean isPointSurface() { throw new UnsupportedOperationException(); }
	@Override public double getAveStrike() { throw new UnsupportedOperationException(); }
	@Override public double getAreaInsideRegion(Region region) { throw new UnsupportedOperationException(); }
	@Override public FaultTrace getEvenlyDiscritizedUpperEdge() { throw new UnsupportedOperationException(); }
	@Override public LocationList getEvenlyDiscritizedLowerEdge() { throw new UnsupportedOperationException(); }
	@Override public double getAveGridSpacing() { throw new UnsupportedOperationException(); }
	@Override public double getAveDipDirection() { throw new UnsupportedOperationException(); }
	@Override public FaultTrace getUpperEdge() { throw new UnsupportedOperationException(); }
	@Override public double getFractionOfSurfaceInRegion(Region region) { throw new UnsupportedOperationException(); }
	@Override public String getInfo() { throw new UnsupportedOperationException(); }
	@Override public double getMinDistance(RuptureSurface surface) { throw new UnsupportedOperationException(); }
	@Override public RuptureSurface getMoved(LocationVector v) { throw new UnsupportedOperationException(); }
	@Override public RuptureSurface copyShallow() { throw new UnsupportedOperationException(); }
}
