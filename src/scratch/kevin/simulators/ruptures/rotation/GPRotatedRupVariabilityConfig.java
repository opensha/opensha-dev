package scratch.kevin.simulators.ruptures.rotation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.io.FileUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Somerville_2006_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.simulators.RectangularElement;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class GPRotatedRupVariabilityConfig extends RotatedRupVariabilityConfig<GPRotatedRupture> {
	
	private List<GPRotatedRupture> ruptures;
	
	public GPRotatedRupVariabilityConfig(List<Site> sites, List<GPRotatedRupture> ruptures,
			double[] distances, int numSourceAz, int numSiteToSourceAz) {
		this(ruptures, buildRotations(sites, getIDs(ruptures), distances, numSourceAz, numSiteToSourceAz));
	}
	
	private static List<Integer> getIDs(List<GPRotatedRupture> ruptures) {
		List<Integer> ids = new ArrayList<>();
		for (int i=0; i<ruptures.size(); i++)
			ids.add(i);
		return ids;
	}

	public GPRotatedRupVariabilityConfig(List<GPRotatedRupture> ruptures, List<RotationSpec> rotations) {
		super(rotations);
		this.ruptures = ruptures;
	}

	public List<GPRotatedRupture> getRuptureList() {
		return ruptures;
	}
	
	@Override
	public boolean hasRuptures() {
		return true;
	}

	private static final double trans_p_diff_thresh = 0.5;
	private static final double trans_abs_diff_thresh = 0.2;
	private static final int min_translations = 2;
	private static final int max_translations = 100;
	
	@Override
	protected GPRotatedRupture loadRotatedRupture(RotationSpec rotation) {
		GPRotatedRupture rupture = getInitiallyOrientedRupture(rotation.eventID);
		
//		Location siteLoc = key.site
//		Location southernPoint
		
		if (rotation.sourceAz != null) {
			// first rotate the rupture around its centroid
			RotationSpec centroidRotSpec = new RotationSpec(-1, null, rotation.eventID, null, rotation.sourceAz, null);
			GPRotatedRupture rotated = rotationCache.getIfPresent(centroidRotSpec);
			if (rotated == null) {
				// not yet cached
				if (D) System.out.println("Rotating about centroid to az="+rotation.sourceAz);
				rotated = rupture.getTransRotated(rupture.centroid, rotation.sourceAz, null, false);
				rotationCache.put(centroidRotSpec, rotated);
			}
			rupture = rotated;
		}
		Preconditions.checkNotNull(rupture);
		
		if (rotation.distance != null) {
			// now translate it to the supplied distance
			RotationSpec transSpec = new RotationSpec(-1, rotation.site, rotation.eventID, rotation.distance, rotation.sourceAz, null);
			GPRotatedRupture translated = rotationCache.getIfPresent(transSpec);
			if (translated == null) {
				// not yet cached, have to do it
				if (D) System.out.println("Translating to distance: "+rotation.distance);
				
				// first move the centroid to the desired position
				if (D) System.out.println("Centroid: "+rupture.centroid);
				Location targetCentroidLoc = LocationUtils.location(rotation.site.getLocation(), 0d, rotation.distance);
				if (D) System.out.println("Target centroid: "+targetCentroidLoc);
				LocationVector initialVector = LocationUtils.vector(rupture.centroid, targetCentroidLoc);
				if (D) System.out.println("Initial translation: "+initialVector);
				translated = rupture.getTransRotated(null, 0d, initialVector, true);
				
				// now the centroid is the rJB away from the site, but we want the actual rupture to be that distance away
				
				// locate the point that is at the southernmost rupture latitude
				double minLat = Double.POSITIVE_INFINITY;
				for (Location loc : translated.src.getSurface().getRectangle())
					minLat = Double.min(minLat, loc.getLatitude());
				Location southOfCentroidLoc = new Location(minLat, translated.centroid.getLongitude());
				double distSouthCentroid = LocationUtils.horzDistanceFast(translated.centroid, southOfCentroidLoc);
				
				if (D) System.out.println("Current south-of-centroid is "+distSouthCentroid+" km south at: "+southOfCentroidLoc);
				
				// move the rupture North that amount. now the southernmost point will be on the the line of latitude
				// that is rJB away from the site (may form a triangle though with that line and the site)
				if (D) System.out.println("Translating north "+distSouthCentroid+" km");
				LocationVector southCentroidVector = new LocationVector(0d, distSouthCentroid, 0d);
				translated = translated.getTransRotated(null, 0d, southCentroidVector, true);
				
				// now adjust as necessary to account for nonplanar and buried ruptures
				
				int numTrans = 0;
				double minDist = Double.NaN, pDiff = Double.NaN, absDiff = Double.NaN;
				double angleDiff = Double.NaN;
				double rupAngle = Double.NaN;
				double transDist = Double.NaN;
				double origTransDist = Double.NaN;
				LocationVector siteToRup = null;
				LocationVector transVector = null;
				Location closest = null;
				while (true) {
					if (D) System.out.println("Translate loop "+numTrans);
					closest = null;
					minDist = Double.POSITIVE_INFINITY;
					LocationList surf = translated.src.getSurface().getQuadSurface().getEvenlyDiscritizedListOfLocsOnSurface();
					for (Location loc : surf) {
						double elemDist;
						if (BBP_PartBValidationConfig.DIST_JB)
							elemDist = LocationUtils.horzDistanceFast(rotation.site.getLocation(), loc);
						else
							elemDist = LocationUtils.linearDistanceFast(rotation.site.getLocation(), loc);
						if (elemDist < minDist) {
							minDist = elemDist;
							closest = loc;
						}
					}
					
					pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
					absDiff = Math.abs(minDist - rotation.distance);
					if (D) System.out.println("Closest is "+minDist+" away: "+closest);
					if (numTrans >= min_translations && (pDiff < trans_p_diff_thresh || absDiff < trans_abs_diff_thresh)
							|| numTrans == max_translations)
						break;
					
					siteToRup = LocationUtils.vector(rotation.site.getLocation(), closest);
					if (D) System.out.println("Vector from site to rupture: "+siteToRup);
					if (BBP_PartBValidationConfig.DIST_JB) {
						origTransDist = siteToRup.getHorzDistance()-rotation.distance;
					} else {
						// find rJB for the desired rRup
						double zClose = closest.getDepth();
						double targetRjb = Math.sqrt(rotation.distance*rotation.distance - zClose*zClose);
						if (D) System.out.println("Target rJB: "+targetRjb+" for clozest with z="+zClose);
						origTransDist = siteToRup.getHorzDistance()-targetRjb;
					}
					if (D) System.out.println("Orig trans dist: "+origTransDist);
					// only move north/south
					rupAngle = siteToRup.getAzimuth();
					angleDiff = RuptureRotationUtils.angleDiff(rupAngle, 0d);
					if (D) System.out.println("Angle diff: "+angleDiff);
					// cap it at 45 degrees as we can get stuck in a loop otherwise
					transDist = origTransDist*Math.cos(Math.toRadians(Math.min(angleDiff, 45)));
					if (D) System.out.println("Trans dist: "+transDist);
					// positive transDist means we are too far North, so we move south
					transVector = new LocationVector(180d, transDist, 0d);
//					translated = RuptureRotationUtils.getTranslated(translated, transVector);
					translated = translated.getTransRotated(null, 0d, transVector, true);
					numTrans++;
				}
				
				if (D) System.out.println("Done with loop with dist: "+minDist);
				
				Preconditions.checkState(pDiff < trans_p_diff_thresh || absDiff < trans_abs_diff_thresh,
						"Translation didn't work after %s rounds for event %s! target: %s, actual: %s"
						+ "\n\tangle: %s, angleDiff: %s, origTransDist: %s, transDist: %s"
						+ "\n\tBefore last translation, siteToRup: %s"
						+ "\n\tLast translation: %s"
						+ "\n\tClosest: %s",
						numTrans, rupture.eventID, rotation.distance, minDist, rupAngle, angleDiff, origTransDist, transDist,
						siteToRup, transVector, closest);
				rotationCache.put(transSpec, translated);
			} else {
//				double minDist = RuptureRotationUtils.calcMinDist(rotation.site.getLocation(), translated,
//						BBP_PartBValidationConfig.DIST_JB);
				Location siteLoc = rotation.site.getLocation();
				QuadSurface surf = translated.src.getSurface().getQuadSurface();
				double minDist = BBP_PartBValidationConfig.DIST_JB ? surf.getDistanceJB(siteLoc) : surf.getDistanceRup(siteLoc);
				double pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
				double absDiff = Math.abs(minDist - rotation.distance);
				Preconditions.checkState(pDiff < trans_p_diff_thresh || absDiff < trans_abs_diff_thresh,
						"Cached translation is wrong! target: %s, actual: %s",
						rotation.distance, minDist);
			}
			rupture = translated;
		}
		
		if (rotation.siteToSourceAz != null)
			// rotate it around the site
			rupture = rupture.getTransRotated(rotation.site.getLocation(), rotation.siteToSourceAz, null, true);
		rotationCache.put(rotation, rupture);
		return rupture;
	}

	@Override
	protected GPRotatedRupture loadInitialOrientationRupture(Integer eventID) {
		return ruptures.get(eventID);
	}

	@Override
	protected Location loadCentroid(GPRotatedRupture rupture) {
		return rupture.centroid;
	}

	@Override
	public void plotRotations(File outputDir, String prefix, List<RotationSpec> rotations, boolean highlightCentroid)
			throws IOException {
		// origin annotation
		List<XYAnnotation> anns = new ArrayList<>();
		HashSet<Site> sites = new HashSet<>();

		System.out.println("Plotting map of "+rotations.size()+" rotations");

		// now rotate
		// just going to use fake elements here to make it easy
		List<SimulatorElement> plotElems = new ArrayList<>();
		GPRotatedRupture first = null;
		for (RotationSpec rotation : rotations) {
			Site site = rotation.site;
			if (!sites.contains(site)) {
				anns.add(RuptureRotationUtils.getLocationAnn(0.02, site.getLocation(), Color.BLUE));
				sites.add(site);
			}
			GPRotatedRupture rotated = getRotatedRupture(rotation);
			if (first == null)
				first = rotated;
			Vertex[] vertices = new Vertex[4];
			Location[] corners = rotated.src.getSurface().getRectangle();
			for (int i=0; i<vertices.length; i++)
				vertices[i] = new Vertex(corners[i]);
			SimulatorElement elem = new RectangularElement(-1, vertices, null, -1, -1, 0, 0, 0d, 0d, null, false);
			plotElems.add(elem);
		}

		if (highlightCentroid) {
			anns.add(RuptureRotationUtils.getLocationAnn(0.01, first.centroid, Color.GREEN));
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
			anns.add(RuptureRotationUtils.getLocationAnn(1e-10, loc, Color.WHITE));

		RupturePlotGenerator.OTHER_ELEM_COLOR = new Color(100, 100, 100);
		RupturePlotGenerator.writeMapPlot(plotElems, null, null, outputDir, prefix, null, null, null, null, null, null, anns);
	}
	
	public static void writeRuptures(File rupDir, List<GPRotatedRupture> ruptures) throws IOException {
		for (GPRotatedRupture rupture : ruptures) {
			File srcFile = new File(rupDir, "rup_"+rupture.eventID+".src");
			File srfFile = new File(rupDir, "rup_"+rupture.eventID+".srf");
			
			rupture.src.writeToFile(srcFile);
			SRF_PointData.writeSRF(srfFile, rupture.srf, 2d);
		}
	}
	
	public static List<GPRotatedRupture> readRuptures(File rupDir, int num) throws IOException {
		List<GPRotatedRupture> rups = new ArrayList<>();
		for (int i=0; i<num; i++) {
			File srcFile = new File(rupDir, "rup_"+i+".src");
			File srfFile = new File(rupDir, "rup_"+i+".srf");
			
			BBP_SourceFile src = BBP_SourceFile.readFile(srcFile);
			List<SRF_PointData> srf = SRF_PointData.readSRF(srfFile);
			
			rups.add(new GPRotatedRupture(i, src, srf));
		}
		return rups;
	}
	
	public static GPRotatedRupVariabilityConfig loadCSV(File csvFile, File rupDir, List<Site> sites)
			throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		List<RotationSpec> rotations = new ArrayList<>(csv.getNumRows()-1);
		Map<String, Site> sitesMap = new HashMap<>();
		if (sites != null)
			for (Site site : sites)
				sitesMap.put(site.getName(), site);
		
		int maxEventID = 0;
		
		for (int row=1; row<csv.getNumRows(); row++) {
			int index = Integer.parseInt(csv.get(row, 0));
			Preconditions.checkState(index == row-1);
			String siteName = csv.get(row, 1);
			Site site = sitesMap.get(siteName);
			if (site == null) {
				Location loc = new Location(Double.parseDouble(csv.get(row, 2)), Double.parseDouble(csv.get(row, 3)));
				site = new Site(loc, siteName);
				sitesMap.put(siteName, site);
			}
			int eventID = Integer.parseInt(csv.get(row, 4));
			maxEventID = Integer.max(maxEventID, eventID);
			float distance = Float.parseFloat(csv.get(row, 5));
			float sourceAz = Float.parseFloat(csv.get(row, 6));
			float siteToSourceAz = Float.parseFloat(csv.get(row, 7));
			rotations.add(new RotationSpec(index, site, eventID, distance, sourceAz, siteToSourceAz));
		}
		
		List<GPRotatedRupture> ruptures = readRuptures(rupDir, maxEventID+1);
		
		return new GPRotatedRupVariabilityConfig(ruptures, rotations);
	}
	
	public static List<GPRotatedRupture> buildRuptures(Scenario scenario, int num, double elemArea) throws IOException {
		List<Integer> indexes = new ArrayList<>();
		for (int i=0; i<num; i++)
			indexes.add(i);
		return buildRuptures(scenario, indexes, elemArea);
	}
	
	public static List<GPRotatedRupture> buildRuptures(Scenario scenario, List<Integer> indexes, double elemArea) throws IOException {
		Location topCenter = new Location(34d, -118d);
		
		double mag = scenario.getMagnitude();
		double rake;
		switch (scenario.getFaultStyle()) {
		case NORMAL:
			rake = -90;
			break;
		case REVERSE:
			rake = 90;
			break;
		case STRIKE_SLIP:
			rake = 180;
			break;

		default:
			throw new IllegalStateException();
		}
		WC1994_MagLengthRelationship wc = new WC1994_MagLengthRelationship();
		double length = wc.getMedianLength(mag, rake);
		System.out.println("WC 1994 Length for M="+(float)mag+": "+(float)length+" km");
		Somerville_2006_MagAreaRel magArea = new Somerville_2006_MagAreaRel();
		double area = magArea.getMedianArea(mag);
		System.out.println("Somerville 2006 Area for M="+(float)mag+": "+(float)area+" km");
		double width = area/length;
		System.out.println("Calculated wdith for M="+(float)mag+": "+(float)width+" km");
		
		double sqrtElemArea = Math.sqrt(elemArea);
		
		FocalMechanism mech = new FocalMechanism(0d, scenario.getDip(), rake);
		
		Random r = new Random((long)area*(1l + indexes.get(0)));

		List<GPRotatedRupture> rups = new ArrayList<>();
		
		File tmpDir = Files.createTempDir();
		
		for (int index : indexes) {
			BBP_PlanarSurface surface = new BBP_PlanarSurface(topCenter, length, width, mech);
			double hypoAlongStrike = length*(r.nextDouble()-0.5); // 0 is center
			double hypoDownDip = width*r.nextDouble(); // TODO
			BBP_SourceFile src = new BBP_SourceFile(surface, mag, hypoAlongStrike, hypoDownDip,
					sqrtElemArea, sqrtElemArea, 0.15, 1000+index);
			
			File subDir = new File(tmpDir, "rup_"+index);
			subDir.mkdir();
			File srcFile = new File(subDir, "rup.src");
			src.writeToFile(srcFile);
			
			BBP_Wrapper wrapper = new BBP_Wrapper(RSQSimBBP_Config.VM, RSQSimBBP_Config.METHOD, srcFile,
					null, null, null, subDir);
			wrapper.setSRFGenOnly(true);
			wrapper.run();
			
			File srfFile = null;
			for (File file : subDir.listFiles())
				if (file.getName().endsWith(".srf"))
					srfFile = file;
			Preconditions.checkNotNull(srfFile, "Couldn't file SRF file in "+subDir.getAbsolutePath());
			List<SRF_PointData> srf = SRF_PointData.readSRF(srfFile);
			
			rups.add(new GPRotatedRupture(index, src, srf));
		}
		
		FileUtils.deleteDirectory(tmpDir);
		
		return rups;
	}
	
	public static void main(String[] args) throws IOException {
		List<GPRotatedRupture> ruptures = buildRuptures(Scenario.M7p2_VERT_SS_SURFACE, 1, 1d);
		
		Site site = new Site(new Location(34, -118));
		List<Site> sites = new ArrayList<>();
		sites.add(site);
		
		double[] distances = { 20d };
		int numSourceAz = 4;
		int numSiteToSourceAz = 4;
		
		GPRotatedRupVariabilityConfig config = new GPRotatedRupVariabilityConfig(
				sites, ruptures, distances, numSourceAz, numSiteToSourceAz);
		
		for (RotationSpec rot : config.getRotations())
			config.getRotatedRupture(rot);
	}

}
