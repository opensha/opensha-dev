package scratch.kevin.simulators.ruptures.rotation;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.simulators.srf.SRF_PointData;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;

public class GPRotatedRupture {
	
	public final int eventID;
	public final BBP_SourceFile src;
	public final List<SRF_PointData> srf;
	public final Location centroid;
	
	public GPRotatedRupture(int eventID, BBP_SourceFile src, List<SRF_PointData> srf) {
		this.eventID = eventID;
		this.src = src;
		srf = new ArrayList<>(srf);
		for (int i=0; i<srf.size(); i++) {
			SRF_PointData point = srf.get(i);
			Location loc = point.getLocation();
			if (loc.getLongitude() > 180) {
				loc = new Location(loc.getLatitude(), loc.getLongitude()-360);
				srf.set(i, point.translated(loc));
			}
			Preconditions.checkState(loc.getLongitude() <= 180);
		}
		this.centroid = calcCentroid(srf);
		this.srf = srf;
	}
	
	private Location calcCentroid(List<SRF_PointData> srf) {
		List<Double> weights = new ArrayList<>();
		List<Double> lats = new ArrayList<>();
		List<Double> lons = new ArrayList<>();
		
		for (int i=0; i<srf.size(); i++) {
			SRF_PointData point = srf.get(i);
			weights.add(FaultMomentCalc.getMoment(point.getArea(), point.getTotalSlip()));
			Location loc = point.getLocation();
			lats.add(loc.getLatitude());
			lons.add(loc.getLongitude());
		}
		double lat = FaultUtils.getScaledAngleAverage(weights, lats);
		double lon = FaultUtils.getScaledAngleAverage(weights, lons);
		while (lon > 180)
			lon -= 360;
		while (lat > 90)
			lat -= 360;
		return new Location(lat, lon);
	}
	
	public GPRotatedRupture getTransRotated(Location rotOrigin, double rotationAz,
			LocationVector transVector, boolean transFirst) {
		
		Preconditions.checkArgument((rotOrigin != null && rotationAz != 0) || transVector != null,
				"Must do translation or rotation (or both)");
		
		List<SRF_PointData> newSRF = new ArrayList<>();
		
		for (int i=0; i<srf.size(); i++) {
			SRF_PointData point = srf.get(i);
			
			SRF_PointData newPoint = null;
			if (transFirst && transVector != null)
				newPoint = translate(point, transVector);
			
			if (rotOrigin != null && rotationAz != 0d)
				newPoint = rotate(newPoint == null ? point : newPoint, rotOrigin, rotationAz);
			
			if (!transFirst && transVector != null)
				newPoint = translate(newPoint == null ? point : newPoint, transVector);
			newSRF.add(newPoint);
		}
		
		BBP_PlanarSurface origSurf = src.getSurface();
		
		FocalMechanism rotMech = origSurf.getFocalMechanism();
		if (rotationAz != 0d)
			rotMech = rotateMech(origSurf.getFocalMechanism(), rotationAz);
		
		Location topCenter = origSurf.getTopCenter();
		if (transFirst && transVector != null)
			topCenter = LocationUtils.location(topCenter, transVector);
		
		if (rotOrigin != null && rotationAz != 0d)
			topCenter = rotateLoc(topCenter, rotOrigin, rotationAz);
		
		if (!transFirst && transVector != null)
			topCenter = LocationUtils.location(topCenter, transVector);
		
		BBP_PlanarSurface rotSurf = new BBP_PlanarSurface(topCenter, origSurf.getLength(), origSurf.getWidth(), rotMech);
		
		BBP_SourceFile rotSrc = new BBP_SourceFile(rotSurf, src.getMag(), src.getHypoAlongStrike(),
				src.getHypoDownDip(), src.getdWid(), src.getdLen(), src.getCornerFreq(), src.getSeed());
		
		return new GPRotatedRupture(eventID, rotSrc, newSRF);
	}
	
	private static FocalMechanism rotateMech(FocalMechanism mech, double azimuth) {
		double newStrike = mech.getStrike() + azimuth;
		while (newStrike >= 360)
			newStrike -= 360;
		while (newStrike < 0)
			newStrike += 360;
		return new FocalMechanism(newStrike, mech.getDip(), mech.getRake());
	}
	
	/*
	 * returns cloned rotated element around given origin. azimuth in decimal degrees
	 */
	private static SRF_PointData rotate(SRF_PointData point, Location origin, double azimuth) {
		Location rotLoc = rotateLoc(point.getLocation(), origin, azimuth);
		
		FocalMechanism mech = point.getFocalMech();
		FocalMechanism rotMech = null;
		if (mech != null)
			rotMech = rotateMech(mech, azimuth);
		
		return point.translated(rotLoc, rotMech);
	}
	
	private static Location rotateLoc(Location loc, Location origin, double azimuth) {
		LocationVector vector = LocationUtils.vector(origin, loc);
		vector.set(vector.getAzimuth()+azimuth, vector.getHorzDistance(), vector.getVertDistance());
		return LocationUtils.location(origin, vector);
	}
	
	/*
	 * returns cloned translated element.
	 */
	private static SRF_PointData translate(SRF_PointData point, LocationVector vector) {
		Preconditions.checkState(vector.getVertDistance() == 0d, "Vertical should always be zero: %s", vector);
		Location transLoc = LocationUtils.location(point.getLocation(), vector);
		return point.translated(transLoc);
	}
}
