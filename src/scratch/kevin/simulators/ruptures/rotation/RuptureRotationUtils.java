package scratch.kevin.simulators.ruptures.rotation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.utm.UTM;
import org.opensha.commons.geo.utm.WGS84;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.RectangularElement;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.TriangularElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RuptureRotationUtils {
	
	public static RSQSimEvent getRotated(RSQSimEvent event, Location origin, double rotationAz) {
		return getTransRotated(event, origin, rotationAz, null, false);
	}
	
	public static RSQSimEvent getTranslated(RSQSimEvent event, LocationVector vector) {
		return getTransRotated(event, null, 0d, vector, true);
	}
	
	public static RSQSimEvent getTransRotated(RSQSimEvent event, Location rotOrigin, double rotationAz,
			LocationVector transVector, boolean transFirst) {
		List<RSQSimEventRecord> newRecords = new ArrayList<>();
		
		Preconditions.checkArgument((rotOrigin != null && rotationAz != 0) || transVector != null,
				"Must do translation or rotation (or both)");
		
		for (EventRecord record : event) {
			List<SimulatorElement> elems = record.getElements();
			List<SimulatorElement> newElems = new ArrayList<>();
			
			for (int i=0; i<elems.size(); i++) {
				SimulatorElement elem = elems.get(i);
				
				SimulatorElement newElem = null;
				if (transFirst && transVector != null)
					newElem = translate(elem, transVector);
				
				if (rotOrigin != null && rotationAz != 0d)
					newElem = rotate(newElem == null ? elem : newElem, rotOrigin, rotationAz);
				
				if (!transFirst && transVector != null)
					newElem = translate(newElem == null ? elem : newElem, transVector);
				newElems.add(newElem);
			}
			
			newRecords.add(new RelocatedRecord((RSQSimEventRecord)record, newElems));
		}
		
		RSQSimEvent newEvent = new RSQSimEvent(newRecords);
		newEvent.setNextEventTime(event.getNextEventTime());
		return newEvent;
	}
	
	public static RSQSimEvent getMirroredNS(RSQSimEvent event, double latitude) {
		List<RSQSimEventRecord> newRecords = new ArrayList<>();
		
		for (EventRecord record : event) {
			List<SimulatorElement> elems = record.getElements();
			List<SimulatorElement> newElems = new ArrayList<>();
			
			for (int i=0; i<elems.size(); i++) {
				SimulatorElement elem = elems.get(i);
				
				SimulatorElement newElem = mirrorNS(elem, latitude);
				newElems.add(newElem);
			}
			
			newRecords.add(new RelocatedRecord((RSQSimEventRecord)record, newElems));
		}
		
		RSQSimEvent newEvent = new RSQSimEvent(newRecords);
		newEvent.setNextEventTime(event.getNextEventTime());
		return newEvent;
	}
	
	private static class RelocatedRecord extends RSQSimEventRecord {

		private List<SimulatorElement> newElements;

		public RelocatedRecord(RSQSimEventRecord origRec, List<SimulatorElement> newElements) {
			super(null);
			this.newElements = newElements;
			
			int[] origElemIDs = origRec.getElementIDs();
			Preconditions.checkState(newElements.size() == origElemIDs.length);
			double[] origElemSlips = origRec.getElementSlips();
			Preconditions.checkState(origElemIDs.length == origElemSlips.length);
			double[] times = origRec.getElementTimeFirstSlips();
			for (int i=0; i<origElemIDs.length; i++) {
				if (times == null)
					addSlip(origElemIDs[i], origElemSlips[i]);
				else
					addSlip(origElemIDs[i], origElemSlips[i], times[i]);
			}
			
			Preconditions.checkState(origRec.hasElementSlipsAndIDs());
			Preconditions.checkState(origElemIDs.length == origRec.getElementSlips().length, "Bad array lengths. orig ID len=%s, orig slip len=%s",
					origElemIDs.length, origRec.getElementSlips().length);
			
			setArea(origRec.getArea());
			setDuration(origRec.getDuration());
			Preconditions.checkState(hasElementSlipsAndIDs());
			if (origRec.getFirstPatchToSlip() >= 0)
				setFirstPatchToSlip(origRec.getFirstPatchToSlip());
			setID(origRec.getID());
			setLength(origRec.getLength());
			setMagnitude(origRec.getMagnitude());
			setMoment(origRec.getMoment());
			setSectionID(origRec.getSectionID());
			setTime(origRec.getTime());
		}

		@Override
		public List<SimulatorElement> getElements() {
			return newElements;
		}
		
	}
	
	/*
	 * returns cloned rotated element around given origin. azimuth in decimal degrees
	 */
	private static SimulatorElement rotate(SimulatorElement elem, Location origin, double azimuth) {
		Vertex[] verts = elem.getVertices();
		Vertex[] rotVerts = new Vertex[verts.length];
		
		for (int i=0; i<verts.length; i++) {
			LocationVector vector = LocationUtils.vector(origin, verts[i]);
			vector.set(vector.getAzimuth()+azimuth, vector.getHorzDistance(), vector.getVertDistance());
			Location rotLoc = LocationUtils.location(origin, vector);
			rotVerts[i] = new Vertex(rotLoc, verts[i].getID(), verts[i].getDAS(), verts[i].getTraceFlag());
		}
		
		FocalMechanism mech = elem.getFocalMechanism();
		FocalMechanism rotMech = null;
		if (mech != null) {
			double newStrike = mech.getStrike() + azimuth;
			while (newStrike >= 360)
				newStrike -= 360;
			while (newStrike < 0)
				newStrike += 360;
			rotMech = new FocalMechanism(newStrike, mech.getDip(), mech.getRake());
		}
		
		if (elem instanceof TriangularElement)
			return new TriangularElement(elem.getID(), rotVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(),
					elem.getNumAlongStrike(), elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), rotMech);
		else if (elem instanceof RectangularElement)
			return new RectangularElement(elem.getID(), rotVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(), elem.getNumAlongStrike(),
					elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), rotMech, ((RectangularElement)elem).isPerfect());
		throw new IllegalStateException("Only supports triangular and rectangular elements");
	}
	
	/*
	 * returns cloned mirrored element around given origin
	 */
	private static SimulatorElement mirrorNS(SimulatorElement elem, double latitude) {
		Vertex[] verts = elem.getVertices();
		Vertex[] rotVerts = new Vertex[verts.length];
		
		for (int i=0; i<verts.length; i++) {
			Location ptAtMirror = new Location(latitude, verts[i].getLongitude(), verts[i].getDepth());
			LocationVector vector = LocationUtils.vector(verts[i], ptAtMirror);
			Location rotLoc = LocationUtils.location(ptAtMirror, vector);
			rotVerts[i] = new Vertex(rotLoc, verts[i].getID(), verts[i].getDAS(), verts[i].getTraceFlag());
		}
		
		FocalMechanism mech = elem.getFocalMechanism();
		FocalMechanism rotMech = null;
		if (mech != null) {
			// mirror the strike around 0
			double newStrike = -mech.getStrike();
			while (newStrike >= 360)
				newStrike -= 360;
			while (newStrike < 0)
				newStrike += 360;
			rotMech = new FocalMechanism(newStrike, mech.getDip(), mech.getRake());
		}
		
		if (elem instanceof TriangularElement)
			return new TriangularElement(elem.getID(), rotVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(),
					elem.getNumAlongStrike(), elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), rotMech);
		else if (elem instanceof RectangularElement)
			return new RectangularElement(elem.getID(), rotVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(), elem.getNumAlongStrike(),
					elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), rotMech, ((RectangularElement)elem).isPerfect());
		throw new IllegalStateException("Only supports triangular and rectangular elements");
	}
	
	/*
	 * returns cloned translated element.
	 */
	private static SimulatorElement translate(SimulatorElement elem, LocationVector vector) {
		Preconditions.checkState(vector.getVertDistance() == 0d, "Vertical should always be zero: %s", vector);
		Vertex[] verts = elem.getVertices();
		Vertex[] transVerts = new Vertex[verts.length];
		
		for (int i=0; i<verts.length; i++) {
			Location transLoc = LocationUtils.location(verts[i], vector);
			transVerts[i] = new Vertex(transLoc, verts[i].getID(), verts[i].getDAS(), verts[i].getTraceFlag());
		}
		
		if (elem instanceof TriangularElement)
			return new TriangularElement(elem.getID(), transVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(),
					elem.getNumAlongStrike(), elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), elem.getFocalMechanism());
		else if (elem instanceof RectangularElement)
			return new RectangularElement(elem.getID(), transVerts, elem.getSectionName(), elem.getFaultID(), elem.getSectionID(), elem.getNumAlongStrike(),
					elem.getNumDownDip(), elem.getSlipRate(), elem.getAseisFactor(), elem.getFocalMechanism(), ((RectangularElement)elem).isPerfect());
		throw new IllegalStateException("Only supports triangular and rectangular elements");
	}
	
	private static final boolean centroid_utm = true; 
	public static Location calcRuptureCentroid(RSQSimEvent event) {
		List<SimulatorElement> elems = event.getAllElements();
		double[] slips = event.getAllElementSlips();
		Preconditions.checkState(!elems.isEmpty());
		
		double lat, lon;
		
		if (centroid_utm) {
			double totWeight = 0d;
			double northing = 0d;
			double easting = 0d;
			
			int zone = -1;
			char letter = 'Z';
			for (int i=0; i<elems.size(); i++) {
				SimulatorElement elem = elems.get(i);
				double weight = FaultMomentCalc.getMoment(elem.getArea(), slips[i]);
				totWeight += weight;
				Location loc = elem.getCenterLocation();
				WGS84 wgs = new WGS84(loc.getLatitude(), loc.getLongitude());
				UTM utm;
				if (i == 0) {
					utm = new UTM(wgs);
					zone = utm.getZone();
					letter = utm.getLetter();
				} else {
					utm = new UTM(wgs, zone, letter);
				}
				northing += weight*utm.getNorthing();
				easting += weight*utm.getEasting();
			}
			northing /= totWeight;
			easting /= totWeight;
			
			UTM utm = new UTM(zone, letter, easting, northing);
			WGS84 wgs = new WGS84(utm);
			lat = wgs.getLatitude();
			lon = wgs.getLongitude();
		} else {
			List<Double> weights = new ArrayList<>();
			List<Double> lats = new ArrayList<>();
			List<Double> lons = new ArrayList<>();
			
			for (int i=0; i<elems.size(); i++) {
				SimulatorElement elem = elems.get(i);
				weights.add(FaultMomentCalc.getMoment(elem.getArea(), slips[i]));
				Location loc = elem.getCenterLocation();
				lats.add(loc.getLatitude());
				lons.add(loc.getLongitude());
			}
			lat = FaultUtils.getScaledAngleAverage(weights, lats);
			lon = FaultUtils.getScaledAngleAverage(weights, lons);
		}
		
		while (lon > 180)
			lon -= 360;
		while (lat > 90)
			lat -= 360;
		
		return new Location(lat, lon);
	}
	
	protected static final boolean D = false;
	
	public static RSQSimEvent getInitiallyOriented(RSQSimCatalog catalog, RSQSimEvent rupture,
			Location centroid) {
		List<FaultSectionPrefData> allSubSects = catalog.getU3SubSects();
		int offset;
		try {
			offset = RSQSimUtils.getSubSectIndexOffset(catalog.getElements(), catalog.getU3SubSects());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		// rotate it such that average strike is 0
		List<Double> strikes = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		ArrayList<SimulatorElement> rupElems = rupture.getAllElements();
		double[] elemSlips = rupture.getAllElementSlips();
		for (int i=0; i<rupElems.size(); i++) {
			SimulatorElement elem = rupElems.get(i);
			FaultSectionPrefData sect = allSubSects.get(elem.getSectionID()-offset);
			double elemStrike = elem.getFocalMechanism().getStrike();
			
			// check to see if it's flipped ~180 from the section strike (happens often for SS faults)
			// and correct if necessary
			double sectStrike = sect.getFaultTrace().getAveStrike();
			double strikeDiff = angleDiff(sectStrike, elemStrike);
			if (strikeDiff > 120)
				elemStrike += 180;
			
			strikes.add(elemStrike);
			weights.add(FaultMomentCalc.getMoment(elem.getArea(), elemSlips[i]));
		}
		double aveStrike = FaultUtils.getScaledAngleAverage(weights, strikes);
		if (D) System.out.println("Average strike: "+aveStrike);
		RSQSimEvent rotated = RuptureRotationUtils.getRotated(rupture, centroid, -aveStrike);
		
		if (RSQSimRotatedRupVariabilityConfig.HYPO_NORTH) {
			// now make sure the hypocenter is on the North side of the centroid
			Location hypocenter = RSQSimUtils.getHypocenter(rotated);
			if (hypocenter.getLatitude() < centroid.getLatitude()) {
				if (D) System.out.println("Mirroring");
				// flip the rupture horizontally. don't spin it, as that would mess up
				// Aki & Richards convention, mirror it
				rotated = RuptureRotationUtils.getMirroredNS(rotated, centroid.getLatitude());
			}
		}
		
		return rotated;
	}
	
	protected static double angleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		return Math.abs(angleDiff);
	}
	
	public static QuadSurface getIdealizedQuadSurfaceRepresentation(RSQSimEvent oriented, Location centroid) {
		double centroidDepth = Double.NaN;
		double centroidDist = Double.POSITIVE_INFINITY;
		
		double minLat = Double.POSITIVE_INFINITY;
		double maxLat = Double.NEGATIVE_INFINITY;
		List<Double> dips = new ArrayList<>();
		double minDepth = Double.POSITIVE_INFINITY;
		double maxDepth = Double.NEGATIVE_INFINITY;
		for (SimulatorElement elem : oriented.getAllElements()) {
			for (Location loc : elem.getVertices()) {
				double lat = loc.getLatitude();
				minLat = Double.min(minLat, lat);
				maxLat = Double.max(maxLat, lat);
				minDepth = Double.min(minDepth, loc.getDepth());
				maxDepth = Double.max(maxDepth, loc.getDepth());
				
				double cDist = LocationUtils.horzDistanceFast(loc, centroid);
				if (cDist < centroidDist) {
					centroidDepth = loc.getDepth();
					centroidDist = cDist;
				}
			}
			dips.add(elem.getFocalMechanism().getDip());
		}
		
		Location southernExtent = new Location(minLat, centroid.getLongitude());
		Location northernExtent = new Location(maxLat, centroid.getLongitude());
		
		double aveDip = StatUtils.mean(Doubles.toArray(dips));
		boolean dipping = aveDip < 80;
		
		if (dipping) {
			// not strike-slip, figure out trace
			double horzOffset = (centroidDepth-minDepth)/Math.tan(Math.toRadians(aveDip));
//			System.out.println("aveDip="+aveDip+"centroidDepth="+centroidDepth+", minDepth="+minDepth+", horzOffset="+horzOffset);
			double west = 1.5*Math.PI;
			// move west to be the surface projection of trace
			southernExtent = LocationUtils.location(southernExtent, west, horzOffset);
			northernExtent = LocationUtils.location(northernExtent, west, horzOffset);
		}
		FaultTrace tr = new FaultTrace(null);
		tr.add(new Location(southernExtent.getLatitude(), southernExtent.getLongitude(), minDepth));
		tr.add(new Location(northernExtent.getLatitude(), northernExtent.getLongitude(), minDepth));
		double width = (maxDepth-minDepth)/Math.sin(Math.toRadians(aveDip));
		return new QuadSurface(tr, aveDip, width);
	}

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		
		File outputDir = new File(catalog.getCatalogDir(), "rotation_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int numToTest = 5;
		int numRots = 6;
		double minMag = 7d;
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).load();
		Collections.shuffle(events, new Random(events.size()));
		events = events.subList(0, Integer.min(numToTest, events.size()));
		
		double locRectWidth = 0.01;
		
		Location origin = new Location(34, -118);
		double targetDist = 100d;
		
		// centroid rotations
		for (int i=0; i<events.size(); i++) {
			System.out.println("Centroid rotation "+i);
			RSQSimEvent event = events.get(i);
			Location centroid = calcRuptureCentroid(event);
			
			// centroid annotation
			List<XYAnnotation> anns = new ArrayList<>();
			anns.add(getLocationAnn(locRectWidth, centroid, Color.BLUE));
			
			List<SimulatorElement> plotElems = new ArrayList<>();
			double rotAngle = 360d/(double)numRots;
			for (int j=0; j<numRots; j++) {
				double angle = rotAngle * (j+1);
				RSQSimEvent rotEvent = getRotated(event, centroid, angle);
				plotElems.addAll(rotEvent.getAllElements());
			}
			RupturePlotGenerator.writeMapPlot(plotElems, event, null, outputDir, "rot_centroid_test_"+i, null, null, null, null, null, null, anns);
		}
		
		// rotate about origin at fixed distance
		
		boolean rJB = false;
		
		for (int i=0; i<events.size(); i++) {
			System.out.println("Origin rotation "+i);
			RSQSimEvent event = events.get(i);
			
			Location closest = null;
			double minDist = Double.POSITIVE_INFINITY;
			for (SimulatorElement elem : event.getAllElements()) {
				for (Vertex v : elem.getVertices()) {
					double elemDist = LocationUtils.horzDistanceFast(origin, v);
					if (elemDist < minDist) {
						minDist = elemDist;
						closest = v;
					}
				}
			}
			
			// first translate rupture to match expected distance
			System.out.println("\tOrig min dist: "+minDist);
			LocationVector rupToOrigin = LocationUtils.vector(closest, origin);
			LocationVector transVector = new LocationVector(rupToOrigin.getAzimuth(), rupToOrigin.getHorzDistance()-targetDist, 0d);
			event = getTranslated(event, transVector);
			System.out.println("\tTrans min dist: "+calcMinDist(origin, event, rJB));
			
			// origin annotation
			List<XYAnnotation> anns = new ArrayList<>();
			anns.add(getLocationAnn(locRectWidth*10, origin, Color.BLUE));
			
			// now rotate
			List<SimulatorElement> plotElems = new ArrayList<>();
			double rotAngle = 360d/(double)numRots;
			for (int j=0; j<numRots; j++) {
				double angle = rotAngle * (j+1);
				RSQSimEvent rotEvent = getRotated(event, origin, angle);
				System.out.println("\trot "+j+" min dist: "+calcMinDist(origin, rotEvent, rJB));
				plotElems.addAll(rotEvent.getAllElements());
			}
			
			// add tiny annotations at the extremes to force it to plot everything
			MinMaxAveTracker latTrack = new MinMaxAveTracker();
			MinMaxAveTracker lonTrack = new MinMaxAveTracker();
			for (SimulatorElement elem : plotElems) {
				for (Vertex v : elem.getVertices()) {
					latTrack.addValue(v.getLatitude());
					lonTrack.addValue(v.getLongitude());
				}
			}
			Location[] rectangle = new Location[4];
			rectangle[0] = new Location(latTrack.getMax(), lonTrack.getMax());
			rectangle[1] = new Location(latTrack.getMax(), lonTrack.getMin());
			rectangle[2] = new Location(latTrack.getMin(), lonTrack.getMax());
			rectangle[3] = new Location(latTrack.getMin(), lonTrack.getMin());
			for (Location loc : rectangle)
				anns.add(getLocationAnn(1e-10, loc, Color.WHITE));
			
			RupturePlotGenerator.writeMapPlot(plotElems, event, null, outputDir, "rot_origin_dist_test_"+i, null, null, null, null, null, null, anns);
		}
	}
	
	static double calcMinDist(Location loc, RSQSimEvent event, boolean rJB) {
		double minDist = Double.POSITIVE_INFINITY;
		for (SimulatorElement elem : event.getAllElements()) {
			for (Vertex v : elem.getVertices()) {
				double elemDist;
				if (rJB)
					elemDist = LocationUtils.horzDistanceFast(loc, v);
				else
					elemDist = LocationUtils.linearDistanceFast(loc, v);
				if (elemDist < minDist) {
					minDist = elemDist;
				}
			}
		}
		return minDist;
	}
	
	static XYPolygonAnnotation getLocationAnn(double locRectWidth, Location loc, Color c) {
		double[] poly = new double[10];
		double lat = loc.getLatitude();
		double lon = loc.getLongitude();
		double ux = lon+0.5*locRectWidth;
		double lx = lon-0.5*locRectWidth;
		double uy = lat+0.5*locRectWidth;
		double ly = lat-0.5*locRectWidth;
		poly[0] = ux;
		poly[1] = uy;
		poly[2] = lx;
		poly[3] = uy;
		poly[4] = lx;
		poly[5] = ly;
		poly[6] = ux;
		poly[7] = ly;
		poly[8] = ux;
		poly[9] = uy;
		return new XYPolygonAnnotation(poly, null, null, c);
	}

}
