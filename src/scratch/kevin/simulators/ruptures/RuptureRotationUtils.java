package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.RectangularElement;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.TriangularElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;

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
		
		return new RSQSimEvent(newRecords);
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
		
		return new RSQSimEvent(newRecords);
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
			setNextSlipTimes(origRec.getNextSlipTimes());
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
	
	static Location calcRuptureCentroid(RSQSimEvent event) {
		List<SimulatorElement> elems = event.getAllElements();
		Preconditions.checkState(!elems.isEmpty());
		
		List<Double> lats = new ArrayList<>();
		List<Double> lons = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		for (SimulatorElement elem : elems) {
			weights.add(elem.getArea());
			Location loc = elem.getCenterLocation();
			lats.add(loc.getLatitude());
			lons.add(loc.getLongitude());
		}
		
		double lat = FaultUtils.getScaledAngleAverage(weights, lats);
		if (lat > 90)
			lat -= 360;
		double lon = FaultUtils.getScaledAngleAverage(weights, lons);
		if (lon > 180)
			lon -= 360;
		return new Location(lat, lon);
		
		// TODO the below isn't working right
//		List<LocationVector> vectors = new ArrayList<>();
//		List<Double> weights = new ArrayList<>();
//		double totWeight = 0d;
//		
//		vectors.add(new LocationVector(0d, 0d, 0d));
//		double weight0 = elems.get(0).getArea();
//		weights.add(weight0);
//		totWeight += weight0;
//		
//		Location loc0 = elems.get(0).getCenterLocation();
//		loc0 = new Location(loc0.getLatitude(), loc0.getLongitude()); // strip depth
//		for (int i=1; i<elems.size(); i++) {
//			SimulatorElement elem = elems.get(i);
//			Location loc1 = elem.getCenterLocation();
//			loc1 = new Location(loc1.getLatitude(), loc1.getLongitude()); // strip depth
//			LocationVector vector = LocationUtils.vector(loc0, loc1);
//			vectors.add(vector);
//			double weight = elem.getArea();
//			weights.add(weight);
//			totWeight += weight;
//		}
//		
//		List<Double> angles = new ArrayList<>();
//		double hDistSum = 0d;
//		for (int i=0; i<vectors.size(); i++) {
//			LocationVector vector = vectors.get(i);
//			angles.add(vector.getAzimuth());
//			hDistSum += vector.getHorzDistance()*weights.get(i);
//		}
//		hDistSum /= totWeight;
//		for (int i=0; i<weights.size(); i++)
//			weights.set(i, weights.get(0)/totWeight);
//		double az = FaultUtils.getScaledAngleAverage(weights, angles);
//		
//		LocationVector centroidVector = new LocationVector(az, hDistSum, 0d);
//		return LocationUtils.location(loc0, centroidVector);
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
			System.out.println("\tTrans min dist: "+calcMinDist(origin, event));
			
			// origin annotation
			List<XYAnnotation> anns = new ArrayList<>();
			anns.add(getLocationAnn(locRectWidth*10, origin, Color.BLUE));
			
			// now rotate
			List<SimulatorElement> plotElems = new ArrayList<>();
			double rotAngle = 360d/(double)numRots;
			for (int j=0; j<numRots; j++) {
				double angle = rotAngle * (j+1);
				RSQSimEvent rotEvent = getRotated(event, origin, angle);
				System.out.println("\trot "+j+" min dist: "+calcMinDist(origin, rotEvent));
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
	
	static double calcMinDist(Location loc, RSQSimEvent event) {
		double minDist = Double.POSITIVE_INFINITY;
		for (SimulatorElement elem : event.getAllElements()) {
			for (Vertex v : elem.getVertices()) {
				double elemDist = LocationUtils.horzDistanceFast(loc, v);
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
