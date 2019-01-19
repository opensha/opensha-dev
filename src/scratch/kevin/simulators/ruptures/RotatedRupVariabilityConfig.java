package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ExecutionException;

import org.jfree.chart.annotations.XYAnnotation;
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
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVariabilityConfig {

	private RSQSimCatalog catalog;
	private Map<Integer, RSQSimEvent> idToOrigMap;
	private List<RotationSpec> rotations;
	private Map<Site, Table<Integer, Float, List<RotationSpec>>> siteRotations;
	
	// caches
	private static final int max_cache_size = 5000;
	private LoadingCache<RotationSpec, RSQSimEvent> rotationCache;
	private LoadingCache<Integer, RSQSimEvent> initialOrientationCache;
	private Map<RSQSimEvent, Location> centroidCache;
	
	public RotatedRupVariabilityConfig(RSQSimCatalog catalog, List<Site> sites, List<RSQSimEvent> ruptures,
			double[] distances, int numSourceAz, int numSiteToSourceAz) {
		this(catalog, ruptures, buildRotations(sites, ruptures, distances, numSourceAz, numSiteToSourceAz));
	}
	
	public RotatedRupVariabilityConfig(RSQSimCatalog catalog, List<RSQSimEvent> ruptures, List<RotationSpec> rotations) {
		this.catalog = catalog;
		if (ruptures != null) {
			idToOrigMap = new HashMap<>();
			for (RSQSimEvent rupture : ruptures)
				idToOrigMap.put(rupture.getID(), rupture);
			
			rotationCache = CacheBuilder.newBuilder().maximumSize(max_cache_size).build(
					new CacheLoader<RotationSpec, RSQSimEvent>() {

						@Override
						public RSQSimEvent load(RotationSpec key) throws Exception {
							// TODO Auto-generated method stub
							return RotatedRupVariabilityConfig.this.loadRupture(key);
						}
				
			});
			initialOrientationCache = CacheBuilder.newBuilder().build(
					new CacheLoader<Integer, RSQSimEvent>() {

						@Override
						public RSQSimEvent load(Integer key) throws Exception {
							return getInitialOrientation(idToOrigMap.get(key));
						}
				
			});
			centroidCache = new HashMap<>();
		}
		this.rotations = rotations;
		siteRotations = new HashMap<>();
		for (RotationSpec rotation : rotations) {
			Table<Integer, Float, List<RotationSpec>> rotationsTable = siteRotations.get(rotation.site);
			if (rotationsTable == null) {
				rotationsTable = HashBasedTable.create();
				siteRotations.put(rotation.site, rotationsTable);
			}
			List<RotationSpec> subRots = rotationsTable.get(rotation.eventID, rotation.distance);
			if (subRots == null) {
				subRots = new ArrayList<>();
				rotationsTable.put(rotation.eventID, rotation.distance, subRots);
			}
			subRots.add(rotation);
		}
	}
	
	public static class RotationSpec {
		public final int index;
		public final Site site;
		public final int eventID;
		public final Float distance;
		public final Float sourceAz;
		public final Float siteToSourceAz;
		
		public RotationSpec(int index, Site site, int eventID, Float distance, Float sourceAz, Float siteToSourceAz) {
			this.index = index;
			this.site = site;
			this.eventID = eventID;
			this.distance = distance;
			if (sourceAz != null && sourceAz == 0d)
				sourceAz = null;
			this.sourceAz = sourceAz;
			if (siteToSourceAz != null && siteToSourceAz == 0f)
				siteToSourceAz = null;
			this.siteToSourceAz = siteToSourceAz;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((distance == null) ? 0 : distance.hashCode());
			result = prime * result + eventID;
			result = prime * result + ((site == null) ? 0 : site.hashCode());
			result = prime * result + ((siteToSourceAz == null) ? 0 : siteToSourceAz.hashCode());
			result = prime * result + ((sourceAz == null) ? 0 : sourceAz.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			RotationSpec other = (RotationSpec) obj;
			if (distance == null) {
				if (other.distance != null)
					return false;
			} else if (!distance.equals(other.distance))
				return false;
			if (eventID != other.eventID)
				return false;
			if (site == null) {
				if (other.site != null)
					return false;
			} else if (!site.equals(other.site))
				return false;
			if (siteToSourceAz == null) {
				if (other.siteToSourceAz != null)
					return false;
			} else if (!siteToSourceAz.equals(other.siteToSourceAz))
				return false;
			if (sourceAz == null) {
				if (other.sourceAz != null)
					return false;
			} else if (!sourceAz.equals(other.sourceAz))
				return false;
			return true;
		}

		@Override
		public String toString() {
			return "RotationSpec [index="+index+"\tsiteLoc="+site+"\teventID="+eventID+"\tdistance="
					+distance+"\tsourceAz="+sourceAz+"\tsiteToSourceAz="+siteToSourceAz+"]";
		}
		
		private String siteName() {
			if (site == null)
				return "null";
			else if (site.getName() != null)
				return site.getName();
			return (float)site.getLocation().getLatitude()+"_"+(float)site.getLocation().getLongitude();
		}
		
		public String getPrefix() {
			return "i"+index+"_"+siteName()+"_event"+eventID+"_dist"+floatStr(distance)
				+"_srcAz"+floatStr(sourceAz)+"_siteSrcAz"+floatStr(siteToSourceAz);
		}
	}
	
	private static String doubleStr(Double d) {
		if (d == null)
			return "0.0";
		return floatStr(d.floatValue());
	}
	
	private static String floatStr(Float f) {
		if (f == null)
			return "0.0";
		return f.toString();
	}
	
	private static List<RotationSpec> buildRotations(List<Site> sites, List<RSQSimEvent> ruptures, double[] distances,
			int numSourceAz, int numSiteToSourceAz) {
		List<RotationSpec> rotations = new ArrayList<>();
		int index = 0;
		double sourceDeltaAz = 360/(double)numSourceAz;
		double siteToSourceDeltaAz = 360/(double)numSiteToSourceAz;
		for (Site site : sites) {
			for (double distance : distances) {
				for (RSQSimEvent rupture : ruptures) {
					int eventID = rupture.getID();
					for (int nSrc=0; nSrc<numSourceAz; nSrc++) {
						double sourceAz = nSrc*sourceDeltaAz;
						for (int nSite=0; nSite<numSiteToSourceAz; nSite++) {
							double siteToSourceAz = nSite*siteToSourceDeltaAz;
							rotations.add(new RotationSpec(index++, site, eventID, (float)distance, (float)sourceAz, (float)siteToSourceAz));
						}
					}
				}
			}
		}
		return rotations;
	}
	
	public List<RotationSpec> getRotations() {
		return rotations;
	}
	
	public List<RotationSpec> getRotationsForEvent(Site site, int eventID, float distance) {
		Table<Integer, Float, List<RotationSpec>> rotationsTable = siteRotations.get(site);
		Preconditions.checkState(rotationsTable != null, "No rotations for site %s", site.getName());
		return rotationsTable.get(eventID, distance);
	}
	
	public List<RotationSpec> getSiteToSourceRotations(Site site, int eventID, float distance, Float sourceAzimuth) {
		List<RotationSpec> eventRots = getRotationsForEvent(site, eventID, distance);
		if (eventRots == null)
			return null;
		List<RotationSpec> ret = new ArrayList<>();
		if (sourceAzimuth == 0f)
			sourceAzimuth = null;
		
		for (RotationSpec rotation : eventRots) {
			if (Objects.equals(sourceAzimuth, rotation.sourceAz))
				ret.add(rotation);
		}
		return ret;
	}
	
	public List<RotationSpec> getSourceRotations(Site site, int eventID, float distance, Float siteToSourceAzimuth) {
		List<RotationSpec> eventRots = getRotationsForEvent(site, eventID, distance);
		if (eventRots == null)
			return null;
		List<RotationSpec> ret = new ArrayList<>();
		if (siteToSourceAzimuth == 0f)
			siteToSourceAzimuth = null;
		
		for (RotationSpec rotation : eventRots) {
			if (Objects.equals(siteToSourceAzimuth, rotation.siteToSourceAz))
				ret.add(rotation);
		}
		return ret;
	}
	
	private RSQSimEvent getInitialOrientation(RSQSimEvent rupture) {
		Location centroid = centroidCache.get(rupture);
		if (centroid == null) {
			centroid = RuptureRotationUtils.calcRuptureCentroid(rupture);
			centroidCache.put(rupture, centroid);
		}
		
		List<FaultSectionPrefData> sects = catalog.getSubSectsForRupture(rupture, 0.2);
		
		// rotate it such that average strike is 0
		List<Double> strikes = new ArrayList<>();
//		for (SimulatorElement elem : rupture.getAllElements())
//			strikes.add(elem.getFocalMechanism().getStrike());
		for (FaultSectionPrefData sect : sects)
			strikes.add(sect.getFaultTrace().getAveStrike());
		double aveStrike = FaultUtils.getAngleAverage(strikes);
		return RuptureRotationUtils.getRotated(rupture, centroid, -aveStrike);
	}
	
	private RSQSimEvent loadRupture(RotationSpec rotation) {
		RSQSimEvent rupture;
		try {
			rupture = initialOrientationCache.get(rotation.eventID);
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		Preconditions.checkNotNull(rupture);
		
		Location centroid = centroidCache.get(rupture);
		if (centroid == null) {
			centroid = RuptureRotationUtils.calcRuptureCentroid(rupture);
			centroidCache.put(rupture, centroid);
		}
		
		if (rotation.sourceAz != null) {
			// first rotate the rupture around its centroid
			RotationSpec centroidRotSpec = new RotationSpec(-1, null, rotation.eventID, null, rotation.sourceAz, null);
			RSQSimEvent rotated = rotationCache.getIfPresent(centroidRotSpec);
			if (rotated == null) {
				// not yet cached
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
				
				// first move the centroid to the desired position, then adjust as needed
				Location targetLoc = LocationUtils.location(rotation.site.getLocation(), 0d, rotation.distance);
				LocationVector initialVector = LocationUtils.vector(centroid, targetLoc);
				translated = RuptureRotationUtils.getTranslated(rupture, initialVector);
				
				int numTrans = 0;
				double minDist = Double.NaN, pDiff = Double.NaN, absDiff = Double.NaN;
				while (numTrans < 5) { // max 5 translations
					Location closest = null;
					minDist = Double.POSITIVE_INFINITY;
					for (SimulatorElement elem : translated.getAllElements()) {
						for (Vertex v : elem.getVertices()) {
							double elemDist = LocationUtils.horzDistanceFast(rotation.site.getLocation(), v);
							if (elemDist < minDist) {
								minDist = elemDist;
								closest = v;
							}
						}
					}
					
					pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
					absDiff = Math.abs(minDist - rotation.distance);
					if (numTrans > 0 && pDiff < 0.5 || absDiff < 0.5)
						break;
					
					LocationVector rupToOrigin = LocationUtils.vector(closest, rotation.site.getLocation());
					LocationVector transVector = new LocationVector(rupToOrigin.getAzimuth(),
							rupToOrigin.getHorzDistance()-rotation.distance, 0d);
					translated = RuptureRotationUtils.getTranslated(translated, transVector);
					numTrans++;
				}
				
				Preconditions.checkState(pDiff < 0.5 || absDiff < 0.5,
						"Translation didn't work after %s rounds! target: %s, actual: %s, orig: %s",
						numTrans, rotation.distance, minDist);
				rotationCache.put(transSpec, translated);
			} else {
				double minDist = RuptureRotationUtils.calcMinDist(rotation.site.getLocation(), translated);
				double pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
				double absDiff = Math.abs(minDist - rotation.distance);
				Preconditions.checkState(pDiff < 0.5 || absDiff < 0.5,
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
	
	public synchronized RSQSimEvent getRotatedRupture(RotationSpec rotation) {
		Preconditions.checkNotNull(idToOrigMap);
		try {
			return rotationCache.get(rotation);
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	public void writeCSV(File csvFile) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("Index", "Site Name", "Latitude", "Longitude", "Event ID", "Distance",
				"Source Rotation Azimuth", "Site-To-Source Rotation Azimuth");
		for (RotationSpec rotation : rotations) {
			Preconditions.checkState(rotation.site.getName() != null && !rotation.site.getName().isEmpty(),
					"Must supply site names for CSV representation");
			csv.addLine(rotation.index+"", rotation.site.getName(), doubleStr(rotation.site.getLocation().getLatitude()),
					doubleStr(rotation.site.getLocation().getLongitude()), rotation.eventID+"", floatStr(rotation.distance),
					floatStr(rotation.sourceAz), floatStr(rotation.siteToSourceAz));
		}
		csv.writeToFile(csvFile);
	}
	
	public static RotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile) throws IOException {
		return loadCSV(catalog, csvFile, null, null);
	}
	
	public static RotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile, List<RSQSimEvent> events) throws IOException {
		return loadCSV(catalog, csvFile, events, null);
	}
	
	public static RotatedRupVariabilityConfig loadCSV(RSQSimCatalog catalog, File csvFile, List<RSQSimEvent> events, List<Site> sites)
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
		
		return new RotatedRupVariabilityConfig(catalog, events, rotations);
	}
	
	public void plotRotations(File outputDir, String prefix, List<RotationSpec> rotations) throws IOException {
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
				anns.add(RuptureRotationUtils.getLocationAnn(0.02, site.getLocation(), Color.BLUE));
				sites.add(site);
			}
			RSQSimEvent rotated = getRotatedRupture(rotation);
			if (first == null)
				first = rotated;
			plotElems.addAll(rotated.getAllElements());
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
		RupturePlotGenerator.writeMapPlot(plotElems, first, null, outputDir, prefix, null, null, null, null, null, null, anns);
	}
	
	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		int skipYears = 5000;
		int maxRuptures = 100;
		boolean buildAllRuptures = false;
		boolean plotExamples = true;
		
//		File bbpDir = new File("/data/kevin/bbp/parallel/2019_01_17-rundir2585-rotatedRups-m6p6_vert_ss_surface-50.0km"
//				+ "-36srcAz-4siteSrcAz-100rups-skipYears5000-noHF-csLASites");
		File bbpDir = null;
		
//		int[] debugIDs = {14400, 14401, 14402, 14401};
		int[] debugIDs = null;
		
		Scenario scenario = Scenario.M6p6_VERT_SS_SURFACE;
		
		System.out.println("Loading ruptures for scenario");
		List<RSQSimEvent> ruptures = BBP_PartBValidationConfig.Scenario.M6p6_VERT_SS_SURFACE.getMatches(catalog, skipYears);
		System.out.println("Loaded "+ruptures.size()+" ruptures");
		if (ruptures.size() > maxRuptures) {
			ruptures = ruptures.subList(0, maxRuptures);
			System.out.println("Trimmed to "+ruptures.size()+" ruptures");
		}
		
		RotatedRupVariabilityConfig config;
		
		if (bbpDir == null) {
			List<BBP_Site> bbpSites = RSQSimBBP_Config.getCyberShakeInitialLASites();
			
			double[] distances = BBP_PartBValidationConfig.DISTANCES;
			int numSourceAz = 10;
			int numSiteToSourceAz = 10;
			
			List<Site> sites = new ArrayList<>();
			for (BBP_Site bbpSite : bbpSites)
				sites.add(bbpSite.buildGMPE_Site(VelocityModel.LA_BASIN));
			
			config = new RotatedRupVariabilityConfig(catalog,
					sites, ruptures, distances, numSourceAz, numSiteToSourceAz);
		} else {
			File csvFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			config = RotatedRupVariabilityConfig.loadCSV(catalog, csvFile, ruptures);
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
				config.plotRotations(plotDir, "path_rotation_"+index, config.getSiteToSourceRotations(site, rupture.getID(), distance, 0f));
				config.plotRotations(plotDir, "centroid_rotation_"+index, config.getSourceRotations(site, rupture.getID(), distance, 0f));
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
				
				double dist = RuptureRotationUtils.calcMinDist(site.getLocation(), rupture);
				System.out.println("Distance: "+dist);
				
				File outputDir = new File(debugDir, rotation.getPrefix());
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				File srcFile = new File(outputDir, "event_"+rupture.getID()+".src");
				
				BBP_PlanarSurface surface;
				if (RSQSimBBP_Config.U3_SURFACES)
					surface = RSQSimBBP_Config.planarEquivalentU3Surface(catalog, rupture, RSQSimBBP_Config.MIN_SUB_SECT_FRACT, RSQSimBBP_Config.ADJ_WIDTH_MATCH_AREA);
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
				anns.add(RuptureRotationUtils.getLocationAnn(0.02, site.getLocation(), Color.BLUE));
				RupturePlotGenerator.writeMapPlot(null, rupture, func, outputDir, "event_"+rupture.getID(), null, null, null, null, null, null, anns);
			}
		}
	}

}
