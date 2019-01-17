package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RotatedRupVariabilityConfig {

	private Map<Integer, RSQSimEvent> idToOrigMap;
	private List<RotationSpec> rotations;
	
	// caches
	private static final int max_cache_size = 5000;
	private LoadingCache<RotationSpec, RSQSimEvent> rotationCache;
	private Map<RSQSimEvent, Location> centroidCache;
	
	public RotatedRupVariabilityConfig(List<Site> sites, List<RSQSimEvent> ruptures,
			double[] distances, int numSourceAz, int numSiteToSourceAz) {
		this(ruptures, buildRotations(sites, ruptures, distances, numSourceAz, numSiteToSourceAz));
	}
	
	public RotatedRupVariabilityConfig(List<RSQSimEvent> ruptures, List<RotationSpec> rotations) {
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
			centroidCache = new HashMap<>();
		}
		this.rotations = rotations;
	}
	
	public static class RotationSpec {
		public final int index;
		public final Site site;
		public final int eventID;
		public final Double distance;
		public final Double sourceAz;
		public final Double siteToSourceAz;
		
		public RotationSpec(int index, Site site, int eventID, Double distance, Double sourceAz, Double siteToSourceAz) {
			this.index = index;
			this.site = site;
			this.eventID = eventID;
			this.distance = distance;
			if (sourceAz != null && sourceAz == 0d)
				sourceAz = null;
			this.sourceAz = sourceAz;
			if (siteToSourceAz != null && siteToSourceAz == 0d)
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
		
		private String doubleStr(Double d) {
			if (d == null)
				return "0.0";
			return d.floatValue()+"";
		}
		
		private String siteName() {
			if (site == null)
				return "null";
			else if (site.getName() != null)
				return site.getName();
			return (float)site.getLocation().getLatitude()+"_"+(float)site.getLocation().getLongitude();
		}
		
		public String getPrefix() {
			return "i"+index+"_"+siteName()+"_event"+eventID+"_dist"+doubleStr(distance)
				+"_srcAz"+doubleStr(sourceAz)+"_siteSrcAz"+doubleStr(siteToSourceAz);
		}
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
							rotations.add(new RotationSpec(index++, site, eventID, distance, sourceAz, siteToSourceAz));;
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
	
	private RSQSimEvent loadRupture(RotationSpec rotation) {
		RSQSimEvent rupture = idToOrigMap.get(rotation.eventID);
		Preconditions.checkNotNull(rupture);
		
		if (rotation.sourceAz != null) {
			// first rotate the rupture around it's centroid
			RotationSpec centroidRotSpec = new RotationSpec(-1, null, rotation.eventID, null, rotation.sourceAz, null);
			RSQSimEvent rotated = rotationCache.getIfPresent(centroidRotSpec);
			if (rotated == null) {
				// not yet cached
				Location centroid = centroidCache.get(rupture);
				if (centroid == null) {
					centroid = RuptureRotationUtils.calcRuptureCentroid(rupture);
					centroidCache.put(rupture, centroid);
				}
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
				translated = rupture;
				
				double origDist = -1;
				
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
					if (origDist < 0)
						origDist = minDist;
					
					pDiff = DataUtils.getPercentDiff(minDist, rotation.distance);
					absDiff = Math.abs(minDist - rotation.distance);
					if (pDiff < 0.5 || absDiff < 0.5)
						break;
					
					LocationVector rupToOrigin = LocationUtils.vector(closest, rotation.site.getLocation());
					LocationVector transVector = new LocationVector(rupToOrigin.getAzimuth(),
							rupToOrigin.getHorzDistance()-rotation.distance, 0d);
					translated = RuptureRotationUtils.getTranslated(translated, transVector);
					numTrans++;
				}
				
				Preconditions.checkState(pDiff < 0.5 || absDiff < 0.5,
						"Translation didn't work after %s rounds! target: %s, actual: %s, orig: %s",
						numTrans, rotation.distance, minDist, origDist);
				rotationCache.put(transSpec, translated);
				rupture = translated;
			}
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
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		int skipYears = 5000;
		
		List<BBP_Site> bbpSites = RSQSimBBP_Config.getCyberShakeInitialLASites();
		
		double[] distances = BBP_PartBValidationConfig.DISTANCES;
		int numSourceAz = 36;
		int numSiteToSourceAz = 36;
		int maxRuptures = 100;
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site bbpSite : bbpSites)
			sites.add(bbpSite.buildGMPE_Site(VelocityModel.LA_BASIN));
		
		System.out.println("Loading ruptures for scenario");
		List<RSQSimEvent> ruptures = BBP_PartBValidationConfig.Scenario.M6p6_VERT_SS_SURFACE.getMatches(catalog, skipYears);
		System.out.println("Loaded "+ruptures.size()+" ruptures");
		if (ruptures.size() > maxRuptures) {
			ruptures = ruptures.subList(0, maxRuptures);
			System.out.println("Trimmed to "+ruptures.size()+" ruptures");
		}
		
		RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(
				sites, ruptures, distances, numSourceAz, numSiteToSourceAz);
		
		List<RotationSpec> rotations = config.getRotations();
		
		System.out.println("Generated "+rotations.size()+" rotations");
		System.out.println("First 100:");
		for (int i=0; i<100 && i<rotations.size(); i++)
			System.out.println("\t"+rotations.get(i));
		
		System.out.println("Rotating ruptures");
		for (int i=0; i<rotations.size(); i++) {
			if (i % 1000 == 0) System.out.println("\tBuilding Rotated Rupture "+i+"/"+rotations.size()
				+"\t(cache size: "+config.rotationCache.size()+")");
			
			config.getRotatedRupture(rotations.get(i));
		}
		System.out.println("Done rotating ruptures");
	}

}
