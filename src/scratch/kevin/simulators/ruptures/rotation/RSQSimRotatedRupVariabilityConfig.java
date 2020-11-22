package scratch.kevin.simulators.ruptures.rotation;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class RSQSimRotatedRupVariabilityConfig extends RotatedRupVariabilityConfig<RSQSimEvent> {
	
	private RSQSimCatalog catalog;
	private Map<Integer, RSQSimEvent> idToOrigMap;
	
	public RSQSimRotatedRupVariabilityConfig(RSQSimCatalog catalog, List<Site> sites, List<RSQSimEvent> ruptures,
			double[] distances, int numSourceAz, int numSiteToSourceAz) {
		this(catalog, ruptures, buildRotations(sites, getIDs(ruptures), distances, numSourceAz, numSiteToSourceAz));
	}
	
	private static List<Integer> getIDs(List<RSQSimEvent> events) {
		List<Integer> ids = new ArrayList<>();
		for (RSQSimEvent e : events)
			ids.add(e.getID());
		return ids;
	}
	
	@SuppressWarnings("unchecked")
	public RSQSimRotatedRupVariabilityConfig(RSQSimCatalog catalog, Collection<RSQSimEvent> ruptures, List<RotationSpec> rotations) {
		super(rotations);
		this.catalog = catalog;
		if (ruptures != null)
			setRuptures(ruptures);
	}
	
	@Override
	public synchronized boolean hasRuptures() {
		return idToOrigMap != null;
	}
	
	public synchronized void setRuptures(Collection<RSQSimEvent> ruptures) {
		idToOrigMap = new HashMap<>();
		for (RSQSimEvent rupture : ruptures)
			idToOrigMap.put(rupture.getID(), rupture);
	}
	
	@Override
	protected Location loadCentroid(RSQSimEvent rupture) {
		return RuptureRotationUtils.calcRuptureCentroid(rupture);
	}
	
	@Override
	protected RSQSimEvent loadInitialOrientationRupture(Integer eventID) {
		RSQSimEvent rupture = idToOrigMap.get(eventID);
		if (D) System.out.println("Initial orientation for "+rupture.getID());
		Location centroid = getCentroid(rupture);
		if (D) System.out.println("Initial centroid: "+centroid);
		
		return RuptureRotationUtils.getInitiallyOriented(catalog, rupture, centroid);
	}
	
	private static final double trans_p_diff_thresh = 0.5;
	private static final double trans_abs_diff_thresh = 0.2;
	private static final int min_translations = 2;
	private static final int max_translations = 100;
	
	@Override
	protected RSQSimEvent loadRotatedRupture(RotationSpec rotation) {
		RSQSimEvent rupture = getInitiallyOrientedRupture(rotation.eventID);
		Preconditions.checkNotNull(rupture);
		
		Location centroid = getCentroid(rupture);
		if (D) System.out.println("Rotating for "+rupture.getID()+" with centroid "+centroid);
		
		if (rotation.sourceAz != null) {
			// first rotate the rupture around its centroid
			RotationSpec centroidRotSpec = new RotationSpec(-1, null, rotation.eventID, null, rotation.sourceAz, null);
			RSQSimEvent rotated = rotationCache.getIfPresent(centroidRotSpec);
			if (rotated == null) {
				// not yet cached
				if (D) System.out.println("Rotating about centroid to az="+rotation.sourceAz);
				rotated = RuptureRotationUtils.getRotated(rupture, centroid, rotation.sourceAz);
				rotationCache.put(centroidRotSpec, rotated);
			}
			rupture = rotated;
		}
		Preconditions.checkNotNull(rupture);
		
		if (rotation.distance != null) {
			// now translate it to the supplied distance
			RotationSpec transSpec = new RotationSpec(-1, rotation.site, rotation.eventID, rotation.distance, rotation.sourceAz, null);
			RSQSimEvent translated = rotationCache.getIfPresent(transSpec);
			if (translated == null) {
				// not yet cached, have to do it
				if (D) System.out.println("Translating to distance: "+rotation.distance);
				
				// first move the centroid to the desired position
				if (D) System.out.println("Centroid: "+centroid);
				Location targetCentroidLoc = LocationUtils.location(rotation.site.getLocation(), 0d, rotation.distance);
				if (D) System.out.println("Target centroid: "+targetCentroidLoc);
				LocationVector initialVector = LocationUtils.vector(centroid, targetCentroidLoc);
				if (D) System.out.println("Initial translation: "+initialVector);
				translated = RuptureRotationUtils.getTranslated(rupture, initialVector);
				
				centroid = targetCentroidLoc;
				// now the centroid is the rJB away from the site, but we want the actual rupture to be that distance away
				
				// locate the point that is at the southernmost rupture latitude
				double minLat = Double.POSITIVE_INFINITY;
				for (SimulatorElement elem : translated.getAllElements())
					for (Vertex v : elem.getVertices())
						minLat = Double.min(minLat, v.getLatitude());
				Location southOfCentroidLoc = new Location(minLat, centroid.getLongitude());
				double distSouthCentroid = LocationUtils.horzDistanceFast(centroid, southOfCentroidLoc);
				
				if (D) System.out.println("Current south-of-centroid is "+distSouthCentroid+" km south at: "+southOfCentroidLoc);
				
				// move the rupture North that amount. now the southernmost point will be on the the line of latitude
				// that is rJB away from the site (may form a triangle though with that line and the site)
				if (D) System.out.println("Translating north "+distSouthCentroid+" km");
				LocationVector southCentroidVector = new LocationVector(0d, distSouthCentroid, 0d);
				translated = RuptureRotationUtils.getTranslated(translated, southCentroidVector);
				
//				if (rupture.getID() == 86330 && rotation.sourceAz == null && rotation.siteToSourceAz == null) {
//					System.out.println("=====DEBUG=====");
//					System.out.println("Site: "+rotation.site.getLocation());
//					System.out.println("Target location: "+targetCentroidLoc);
//					System.out.println("Original vector: "+initialVector.getAzimuth()+", "+initialVector.getHorzDistance());
//					Location newCentroid = RuptureRotationUtils.calcRuptureCentroid(translated);
//					System.out.println("New centroid: "+newCentroid);
////					initialVector = LocationUtils.vector(newCentroid, targetLoc);
////					System.out.println("New vector: "+initialVector.getAzimuth()+", "+initialVector.getHorzDistance());
//					LocationVector siteToCentroid = LocationUtils.vector(rotation.site.getLocation(), newCentroid);
//					System.out.println("Site to centroid: "+siteToCentroid.getAzimuth()+", "+siteToCentroid.getHorzDistance());
//					System.out.println("Min dist: "+RuptureRotationUtils.calcMinDist(rotation.site.getLocation(), translated));
//				}
				
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
					for (SimulatorElement elem : translated.getAllElements()) {
						for (Vertex v : elem.getVertices()) {
							double elemDist;
							if (BBP_PartBValidationConfig.DIST_JB)
								elemDist = LocationUtils.horzDistanceFast(rotation.site.getLocation(), v);
							else
								elemDist = LocationUtils.linearDistanceFast(rotation.site.getLocation(), v);
							if (elemDist < minDist) {
								minDist = elemDist;
								closest = v;
							}
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
					translated = RuptureRotationUtils.getTranslated(translated, transVector);
					numTrans++;
				}
				
				if (D) System.out.println("Done with loop with dist: "+minDist);
				
				if (numTrans == 20) {
					System.out.println("DEBUGGIN A FAIL!");
					Location newCentroid = RuptureRotationUtils.calcRuptureCentroid(translated);
					System.out.println("\tCentroid should be: "+LocationUtils.location(centroid, southCentroidVector));
					System.out.println("\tCentroid is: "+newCentroid);
					System.out.println("\tVector to centroid: "+LocationUtils.vector(rotation.site.getLocation(), newCentroid));
					System.out.println("\tClosest: "+closest);
					System.out.println("\tVector to closest: "+LocationUtils.vector(rotation.site.getLocation(), closest));
				}
				
				Preconditions.checkState(pDiff < trans_p_diff_thresh || absDiff < trans_abs_diff_thresh,
						"Translation didn't work after %s rounds for event %s! target: %s, actual: %s"
						+ "\n\tangle: %s, angleDiff: %s, origTransDist: %s, transDist: %s"
						+ "\n\tBefore last translation, siteToRup: %s"
						+ "\n\tLast translation: %s"
						+ "\n\tClosest: %s",
						numTrans, rupture.getID(), rotation.distance, minDist, rupAngle, angleDiff, origTransDist, transDist,
						siteToRup, transVector, closest);
				rotationCache.put(transSpec, translated);
			} else {
				double minDist = RuptureRotationUtils.calcMinDist(rotation.site.getLocation(), translated,
						BBP_PartBValidationConfig.DIST_JB);
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
			rupture = RuptureRotationUtils.getRotated(rupture, rotation.site.getLocation(), rotation.siteToSourceAz);
		rotationCache.put(rotation, rupture);
		
		return rupture;
	}
	
	public static RSQSimRotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile) throws IOException {
		return loadCSV(catalog, csvFile, null, null);
	}
	
	public static RSQSimRotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile, List<RSQSimEvent> events) throws IOException {
		return loadCSV(catalog, csvFile, events, null);
	}
	
	public static RSQSimRotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile, List<RSQSimEvent> events, List<Site> sites)
			throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		List<RotationSpec> rotations = new ArrayList<>(csv.getNumRows()-1);
		Map<String, Site> sitesMap = new HashMap<>();
		if (sites != null)
			for (Site site : sites)
				sitesMap.put(site.getName(), site);
		
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
			float distance = Float.parseFloat(csv.get(row, 5));
			float sourceAz = Float.parseFloat(csv.get(row, 6));
			float siteToSourceAz = Float.parseFloat(csv.get(row, 7));
			rotations.add(new RotationSpec(index, site, eventID, distance, sourceAz, siteToSourceAz));
		}
		
		return new RSQSimRotatedRupVariabilityConfig(catalog, events, rotations);
	}
	
	public void plotRotations(File outputDir, String prefix, List<RotationSpec> rotations, boolean highlightCentroid)
			throws IOException {
		// origin annotation
		List<XYAnnotation> anns = new ArrayList<>();
		HashSet<Site> sites = new HashSet<>();

		System.out.println("Plotting map of "+rotations.size()+" rotations");

		// now rotate
		List<SimulatorElement> plotElems = new ArrayList<>();
		RSQSimEvent first = null;
		for (RotationSpec rotation : rotations) {
			Site site = rotation.site;
			if (!sites.contains(site)) {
				anns.add(RuptureRotationUtils.getLocationRectAnn(0.02, site.getLocation(), Color.BLUE));
				sites.add(site);
			}
			RSQSimEvent rotated = getRotatedRupture(rotation);
			Location hypo = RSQSimUtils.getHypocenter(rotated);
			if (first == null) {
				first = rotated;
				
				// draw a circle at the given distance
				double distKM = rotation.distance.doubleValue()*0.99;
				Location siteLoc = site.getLocation();
				Location topLoc = LocationUtils.location(siteLoc, 0d, distKM);
				Location botLoc = LocationUtils.location(siteLoc, Math.PI, distKM);
				Location rightLoc = LocationUtils.location(siteLoc, 0.5*Math.PI, distKM);
				Location leftLoc = LocationUtils.location(siteLoc, 1.5*Math.PI, distKM);
				double circleHeight = topLoc.getLatitude() - botLoc.getLatitude();
				double circleWidth = rightLoc.getLongitude() - leftLoc.getLongitude();
				// docs say says "upper left", but seems to be lower left corner of framing rectangle
				double x = siteLoc.getLongitude() - 0.5*circleWidth;
				double y = siteLoc.getLatitude() - 0.5*circleHeight;
				Shape shape = new Ellipse2D.Double(x, y, circleWidth, circleHeight);
				Stroke stroke = new BasicStroke(2f, BasicStroke.CAP_BUTT,
						BasicStroke.JOIN_BEVEL,0,new float[] {9},0);
				Color circleColor = new Color(0, 0, 0, 60);
				XYShapeAnnotation distCircle =  new XYShapeAnnotation(shape, stroke, circleColor, null);
				anns.add(distCircle);
				Location distLineLoc = LocationUtils.location(siteLoc, Math.PI/6d, distKM);
				XYLineAnnotation distLine = new XYLineAnnotation(
						siteLoc.getLongitude(), siteLoc.getLatitude(),
						distLineLoc.getLongitude(), distLineLoc.getLatitude(), stroke, circleColor);
				anns.add(distLine);
				Location distTextLoc = LocationUtils.location(siteLoc, Math.PI/6d, distKM*0.5);
				XYTextAnnotation distLabel = new XYTextAnnotation(" "+rotation.distance.intValue()+" km",
						distTextLoc.getLongitude(), distTextLoc.getLatitude());
				distLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
				distLabel.setTextAnchor(TextAnchor.TOP_LEFT);
				anns.add(distLabel);
				
//				anns.add(RuptureRotationUtils.getLocationCircleAnn(0.015, hypo, new Color(200, 0, 0)));
				if (highlightCentroid)
					anns.add(RuptureRotationUtils.getLocationCircleAnn(
							0.015, RuptureRotationUtils.calcRuptureCentroid(rotated), new Color(0, 200, 0)));
			} else {
				XYPolygonAnnotation rectHypoPoly = new XYPolygonAnnotation(
						RupturePlotGenerator.star(hypo.getLongitude(), hypo.getLatitude(), 0.0075), new BasicStroke(1f), Color.BLACK, new Color(200, 0, 0, 127));
				anns.add(rectHypoPoly);
//				anns.add(RuptureRotationUtils.getLocationCircleAnn(0.008, hypo, new Color(255, 120, 120)));
				if (highlightCentroid)
					anns.add(RuptureRotationUtils.getLocationCircleAnn(
							0.008, RuptureRotationUtils.calcRuptureCentroid(rotated), new Color(0, 200, 0, 127)));
			}
			plotElems.addAll(rotated.getAllElements());
		}
		
//		// clear out the element time of first slips so that it won't plot the hypocenter star
//		for (EventRecord rec : first)
//			rec.setElementTimeFirstSlips(null);
		
//		if (highlightCentroid) {
//			anns.add(RuptureRotationUtils.getLocationRectAnn(0.01, RuptureRotationUtils.calcRuptureCentroid(first), Color.GREEN));
//		}

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
			anns.add(RuptureRotationUtils.getLocationRectAnn(1e-10, loc, Color.WHITE));

		RupturePlotGenerator.OTHER_ELEM_COLOR = new Color(100, 100, 100);
		RupturePlotGenerator.HYPO_RADIUS = 0.015;
		RupturePlotGenerator.HYPO_COLOR = new Color(200, 0, 0);
		RupturePlotGenerator.writeMapPlot(plotElems, first, null, outputDir, prefix, null, null, null, null, null, null, anns);
	}
	
	public RSQSimRotatedRupVariabilityConfig forSites(List<Site> sites) {
		List<RotationSpec> masterRotations = new ArrayList<>();
		for (Site site : sites) {
			List<RotationSpec> siteRotations = getRotationsForQuantities(Quantity.SITE, site);
			Preconditions.checkNotNull(siteRotations);
			Preconditions.checkState(!siteRotations.isEmpty());
			List<RotationSpec> modRotations = new ArrayList<>();
			for (RotationSpec rot : siteRotations)
				modRotations.add(new RotationSpec(rot.index, site, rot.eventID, rot.distance, rot.sourceAz, rot.siteToSourceAz));
			masterRotations.addAll(modRotations);
		}
		return forRotationSubset(masterRotations);
	}
	
	public RSQSimRotatedRupVariabilityConfig forRotationSubset(List<RotationSpec> rotations) {
		return new RSQSimRotatedRupVariabilityConfig(catalog, idToOrigMap == null ? null : idToOrigMap.values(), rotations);
	}
	
	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		int skipYears = 5000;
		int maxRuptures = 10;
		boolean buildAllRuptures = false;
		boolean plotExamples = true;
		
		// to debug an event
////		RSQSimEvent debugEvent = catalog.loader().byID(970307);
////		RSQSimEvent debugEvent = catalog.loader().byID(1992428);
////		RSQSimEvent debugEvent = catalog.loader().byID(7122655);
//		RSQSimEvent debugEvent = catalog.loader().byID(5802150);
//		List<RSQSimEvent> debugEvents = new ArrayList<>();
//		debugEvents.add(debugEvent);
//		Site debugSite = new Site(new Location(34.0192, -118.286));
//		List<Site> debugSites = new ArrayList<>();
//		debugSites.add(debugSite);
//		RotatedRupVariabilityConfig debugConfig = new RotatedRupVariabilityConfig(catalog, debugSites, debugEvents,
//				new double[] {20}, 36, 1);
//		for (RotationSpec rotation : debugConfig.getRotations())
//			debugConfig.getRotatedRupture(rotation);
//		File debugOut = new File("/tmp/event_"+debugEvent.getID());
//		Preconditions.checkState(debugOut.exists() || debugOut.mkdir());
//		debugConfig.plotRotations(debugOut, "rotation_test", debugConfig.getRotations(), true);
//		System.exit(0);
		
//		File bbpDir = new File("/data/kevin/bbp/parallel/2019_01_17-rundir2585-rotatedRups-m6p6_vert_ss_surface-50.0km"
//				+ "-36srcAz-4siteSrcAz-100rups-skipYears5000-noHF-csLASites");
		File bbpDir = null;
		
//		int[] debugIDs = {14400, 14401, 14402, 14401};
		int[] debugIDs = null;
		
//		Scenario scenario = Scenario.M6p6_VERT_SS_SURFACE;
		Scenario scenario = Scenario.M6p6_REVERSE;
		
		System.out.println("Loading ruptures for scenario");
		List<RSQSimEvent> ruptures = scenario.getMatches(catalog, skipYears);
		System.out.println("Loaded "+ruptures.size()+" ruptures");
		if (ruptures.size() > maxRuptures) {
			ruptures = ruptures.subList(0, maxRuptures);
			System.out.println("Trimmed to "+ruptures.size()+" ruptures");
		}
		
		RSQSimRotatedRupVariabilityConfig config;
		
		if (bbpDir == null) {
			List<BBP_Site> bbpSites = RSQSimBBP_Config.getCyberShakeInitialLASites();
			
			double[] distances = BBP_PartBValidationConfig.OFFICIAL_DISTANCES;
			int numSourceAz = 10;
			int numSiteToSourceAz = 10;
			
			List<Site> sites = new ArrayList<>();
			for (BBP_Site bbpSite : bbpSites)
				sites.add(bbpSite.buildGMPE_Site(null));
			
			config = new RSQSimRotatedRupVariabilityConfig(catalog,
					sites, ruptures, distances, numSourceAz, numSiteToSourceAz);
		} else {
			File csvFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			config = RSQSimRotatedRupVariabilityConfig.loadCSV(catalog, csvFile, ruptures);
		}
		
		List<RotationSpec> rotations = config.getRotations();
		
		System.out.println("Have "+rotations.size()+" rotations");
//		System.out.println("First 100:");
//		for (int i=0; i<100 && i<rotations.size(); i++)
//			System.out.println("\t"+rotations.get(i));
		
		if (buildAllRuptures) {
			System.out.println("Rotating ruptures");
			for (int i=0; i<rotations.size(); i++) {
				if (i % 1000 == 0) System.out.println("\tBuilding Rotated Rupture "+i+"/"+rotations.size()
					+"\t(cache size: "+config.rotationCache.size()+")");
				
				config.getRotatedRupture(rotations.get(i));
			}
			System.out.println("Done rotating ruptures");
		}
		
		if (plotExamples) {
			File mainDir = new File(catalog.getCatalogDir(), "rotation_tests");
			Preconditions.checkState(mainDir.exists() || mainDir.mkdir());
			
			File plotDir = new File(mainDir, "maps_"+scenario.getPrefix());
			Preconditions.checkState(plotDir.exists() || plotDir.mkdir());
			
			Site site = config.getRotations().get(0).site;
			float distance = config.getRotations().get(0).distance;
			int index = 0;
			for (RSQSimEvent rupture : ruptures) {
				System.out.println("Plotting "+index+" for event "+rupture.getID());
				config.plotRotations(plotDir, "path_rotation_"+index, config.getRotationsForQuantities(
						Quantity.SITE, site, Quantity.EVENT_ID, rupture.getID(), Quantity.DISTANCE, distance, Quantity.SOURCE_AZIMUTH, 0f), true);
				config.plotRotations(plotDir, "centroid_rotation_"+index, config.getRotationsForQuantities(
						Quantity.SITE, site, Quantity.EVENT_ID, rupture.getID(), Quantity.DISTANCE, distance, Quantity.SITE_TO_SOURTH_AZIMUTH, 0f), true);
				index++;
			}
		}
		
		if (debugIDs != null) {
			File debugDir = new File(catalog.getCatalogDir(), "rotation_tests");
			Preconditions.checkState(debugDir.exists() || debugDir.mkdir());
			
			for (int index : debugIDs) {
				RotationSpec rotation = rotations.get(index);
				System.out.println("Debugging: "+rotation.getPrefix());
				RSQSimEvent rupture = config.getRotatedRupture(rotation);
				
				List<BBP_Site> sites = new ArrayList<>();
				Site site = rotation.site;
				sites.add(new BBP_Site(site.getName(), site.getLocation(), RSQSimBBP_Config.VM.getVs30(),
						RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ));
				
				double dist = RuptureRotationUtils.calcMinDist(site.getLocation(), rupture, BBP_PartBValidationConfig.DIST_JB);
				System.out.println("Distance: "+dist);
				
				File outputDir = new File(debugDir, rotation.getPrefix());
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				File srcFile = new File(outputDir, "event_"+rupture.getID()+".src");
				
				BBP_PlanarSurface surface;
				if (RSQSimBBP_Config.U3_SURFACES)
					surface = RSQSimBBP_Config.planarEquivalentU3Surface(catalog, rupture, RSQSimBBP_Config.ADJ_WIDTH_MATCH_AREA);
				else
					surface = RSQSimBBP_Config.estimateBBP_PlanarSurface(rupture);
				BBP_SourceFile bbpSource = RSQSimBBP_Config.buildBBP_Source(rupture, surface, RSQSimBBP_Config.DEFAULT_SEED);
				bbpSource.writeToFile(srcFile);
				
				File srfFile = new File(outputDir, "event_"+rupture.getID()+".srf");
				RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(rupture);
				System.out.println("Generating SRF for dt="+(float)RSQSimBBP_Config.SRF_DT+", "+RSQSimBBP_Config.SRF_INTERP_MODE);
				List<SRF_PointData> srf = RSQSimSRFGenerator.buildSRF(func, rupture.getAllElements(), RSQSimBBP_Config.SRF_DT, RSQSimBBP_Config.SRF_INTERP_MODE);
				SRF_PointData.writeSRF(srfFile, srf, RSQSimBBP_Config.SRF_VERSION);
				
				File sitesFile = new File(outputDir, "sites.stl");
				BBP_Site.writeToFile(sitesFile, sites);
				
				BBP_Wrapper bbpWrap = new BBP_Wrapper(RSQSimBBP_Config.VM, RSQSimBBP_Config.METHOD, srcFile, null, srfFile, sitesFile, outputDir);
				bbpWrap.setDoHF(RSQSimBBP_Config.DO_HF);
				bbpWrap.run();
				
				List<XYAnnotation> anns = new ArrayList<>();
				anns.add(RuptureRotationUtils.getLocationRectAnn(0.02, site.getLocation(), Color.BLUE));
				RupturePlotGenerator.writeMapPlot(null, rupture, func, outputDir, "event_"+rupture.getID(), null, null, null, null, null, null, anns);
			}
		}
	}

}
