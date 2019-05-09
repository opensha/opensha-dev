package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ExecutionException;

import org.jfree.chart.annotations.XYAnnotation;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.NamedComparator;
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
import org.opensha.sha.simulators.utils.RSQSimUtils;
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
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVariabilityConfig {

	private RSQSimCatalog catalog;
	private Map<Integer, RSQSimEvent> idToOrigMap;
	private List<RotationSpec> rotations;
	private Map<Quantity, List<?>> quantitiesMap;
	private Table<Quantity, Object, List<RotationSpec>> quantityRotationsCache;
	
	// caches
	private static final int max_cache_size = 5000;
	private LoadingCache<RotationSpec, RSQSimEvent> rotationCache;
	private LoadingCache<Integer, RSQSimEvent> initialOrientationCache;
	private Map<RSQSimEvent, Location> centroidCache;
	
	public enum Quantity {
		SITE("Site", "Unique site locations. If 3-d, each will have unique velocity profiles.", new NamedComparator()),
		EVENT_ID("Rupture", "Unique (but similar in faulting style and magnitude) ruptures which match the "
				+ "given scenario."),
		DISTANCE(BBP_PartBValidationConfig.DIST_JB ? "Joyner-Boore Distance" : "Distance",
				BBP_PartBValidationConfig.DIST_JB ? "Shortest horizontal distance between the site and the surface projection of the rupture."
						: "3-dimensional distance between the site and the rupture surface."),
		SOURCE_AZIMUTH("Rupture Strike", "Rupture strike conforming to the Aki & Richards (1980) convention, where dipping "
				+ "faults dip to the right of the rupture. If path rotation is also performed, this azimuth is relative to the path."),
		SITE_TO_SOURTH_AZIMUTH("Path", "Path from the site to the centroid of the rupture, in azimuthal degrees (0 is North)");
		
		private String name;
		private String description;
		private Comparator<?> comparator;
		
		private Quantity(String name, String description) {
			this(name, description, null);
		}
		
		private Quantity(String name, String description, Comparator<?> comparator) {
			this.name = name;
			this.description = description;
			this.comparator = comparator;
		}
		
		@Override
		public String toString() {
			return name;
		}
		
		public String getName() {
			return name;
		}
		
		public String getDescription() {
			return description;
		}
		
		boolean zeroStoredAsNull() {
			return this == SOURCE_AZIMUTH || this == SITE_TO_SOURTH_AZIMUTH;
		}
	}
	
	public RotatedRupVariabilityConfig(RSQSimCatalog catalog, List<Site> sites, List<RSQSimEvent> ruptures,
			double[] distances, int numSourceAz, int numSiteToSourceAz) {
		this(catalog, ruptures, buildRotations(sites, ruptures, distances, numSourceAz, numSiteToSourceAz));
	}
	
	@SuppressWarnings("unchecked")
	public RotatedRupVariabilityConfig(RSQSimCatalog catalog, Collection<RSQSimEvent> ruptures, List<RotationSpec> rotations) {
		this.catalog = catalog;
		if (ruptures != null)
			setRuptures(ruptures);
		this.rotations = rotations;
		// build quantity lists
		Map<Quantity, HashSet<Object>> quantitySetMap = new HashMap<>();
		for (Quantity quantity : Quantity.values())
			quantitySetMap.put(quantity, new HashSet<>());
		for (RotationSpec rotation : rotations) {
			for (Quantity quantity : Quantity.values()) {
				Object value = rotation.getValue(quantity);
				Preconditions.checkNotNull(value, "Null quantity");
				quantitySetMap.get(quantity).add(value);
			}
		}
		
		quantitiesMap = new HashMap<>();
		for (Quantity quantity : Quantity.values()) {
			HashSet<Object> values = quantitySetMap.get(quantity);
//			System.out.println(quantity.name+": "+values.size()+" unique values");
			List<?> list = new ArrayList<>(values);
			if (quantity.comparator == null) {
				List<Comparable<? super Object>> compList = new ArrayList<>();
				for (Object val : list) {
					Preconditions.checkState(val instanceof Comparable<?>);
					compList.add((Comparable<? super Object>) val);
				}
				Collections.sort(compList);
				list = compList;
			} else {
				list.sort((Comparator<? super Object>) quantity.comparator);
			}
			quantitiesMap.put(quantity, list);
		}

		quantityRotationsCache = HashBasedTable.create();
	}
	
	public synchronized boolean hasRuptures() {
		return idToOrigMap != null;
	}
	
	public synchronized void setRuptures(Collection<RSQSimEvent> ruptures) {
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
	
	@SuppressWarnings("unchecked")
	public <T> List<T> getValues(Class<T> type, Quantity quantity) {
		return (List<T>) quantitiesMap.get(quantity);
	}
	
	public Map<Quantity, List<?>> getQuantitiesMap() {
		return quantitiesMap;
	}
	
	public static class RotationSpec {
		public final int index;
		public final Site site;
		public final int eventID;
		public final Float distance;
		public final Float sourceAz;
		public final Float siteToSourceAz;
		private final Map<Quantity, Object> quantities;
		
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
			
			quantities = new HashMap<>();
			if (site != null)
				quantities.put(Quantity.SITE, site);
			if (eventID >= 0)
				quantities.put(Quantity.EVENT_ID, eventID);
			if (distance != null)
				quantities.put(Quantity.DISTANCE, distance);
			if (sourceAz != null)
				quantities.put(Quantity.SOURCE_AZIMUTH, sourceAz);
			if (siteToSourceAz != null)
				quantities.put(Quantity.SITE_TO_SOURTH_AZIMUTH, siteToSourceAz);
		}
		
		public boolean hasQuantity(Quantity quantity, Object value) {
			if (quantity.zeroStoredAsNull() && Objects.equals(value, 0f))
				value = null;
			if (quantity == Quantity.SITE) {
				Site mySite = (Site)quantities.get(quantity);
				if (value == null)
					return mySite == null;
				else if (mySite == null)
					return false;
				return locEqual(mySite.getLocation(), ((Site)value).getLocation());
			}
			return Objects.equals(value, quantities.get(quantity));
		}
		
		private boolean locEqual(Location loc1, Location loc2) {
			float lat1 = (float)loc1.getLatitude();
			float lon1 = (float)loc1.getLongitude();
			float lat2 = (float)loc2.getLatitude();
			float lon2 = (float)loc2.getLongitude();
			return LocationUtils.areSimilar(loc1, loc2) || (lat1 == lat2 && lon1 == lon2);
		}
		
		public Object getValue(Quantity quantity) {
			Object value = quantities.get(quantity);
			if (value == null && quantity.zeroStoredAsNull())
				return 0f;
			return value;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((distance == null) ? 0 : distance.hashCode());
			result = prime * result + eventID;
			if (site != null) {
				Float lat = (float)site.getLocation().getLatitude();
				Float lon = (float)site.getLocation().getLongitude();
				result = prime * result + lat.hashCode();
				result = prime * result + lon.hashCode();
			}
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
			} else if (!locEqual(site.getLocation(), other.site.getLocation()))
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
	
	private List<RotationSpec> getCachedRotations(Quantity quantity, Object value) {
		List<RotationSpec> ret = quantityRotationsCache.get(quantity, value);
		if (ret == null) {
			ret = new ArrayList<>();
			for (RotationSpec rotation : getRotations()) {
				if (rotation.hasQuantity(quantity, value))
					ret.add(rotation);
			}
			Preconditions.checkNotNull(quantity, "Null quantity?");
			Preconditions.checkNotNull(value, "Null value for quantity %s", quantity);
			quantityRotationsCache.put(quantity, value, ret);
		}
		return ret;
	}
	
	public List<RotationSpec> getRotationsForQuantities(Quantity q1, Object v1) {
		return getRotationsForQuantities(new Quantity[] {q1}, new Object[] {v1});
	}
	
	public List<RotationSpec> getRotationsForQuantities(Quantity q1, Object v1, Quantity q2, Object v2) {
		return getRotationsForQuantities(new Quantity[] {q1, q2}, new Object[] {v1, v2});
	}
	
	public List<RotationSpec> getRotationsForQuantities(Quantity q1, Object v1, Quantity q2, Object v2, Quantity q3, Object v3) {
		return getRotationsForQuantities(new Quantity[] {q1, q2, q3}, new Object[] {v1, v2, v3});
	}
	
	public List<RotationSpec> getRotationsForQuantities(Quantity q1, Object v1, Quantity q2, Object v2, Quantity q3, Object v3,
			Quantity q4, Object v4) {
		return getRotationsForQuantities(new Quantity[] {q1, q2, q3, q4}, new Object[] {v1, v2, v3, v4});
	}
	
	public List<RotationSpec> getRotationsForQuantities(Quantity[] quantities, Object[] values) {
		Preconditions.checkArgument(quantities.length == values.length);
		if (quantities.length == 0)
			return getRotations();
		
		List<RotationSpec> rotations = getCachedRotations(quantities[0], values[0]);
		if (quantities.length == 1)
			return rotations;
		
		return getRotationsForQuantities(rotations, Arrays.copyOfRange(quantities, 1, quantities.length),
				Arrays.copyOfRange(values, 1, values.length));
	}
	
	public static List<RotationSpec> getRotationsForQuantities(List<RotationSpec> rotations, Quantity[] quantities, Object[] values) {
		Preconditions.checkArgument(quantities.length == values.length);
		if (quantities.length == 0)
			return rotations;
		
		List<RotationSpec> ret = new ArrayList<>();
		rotationLoop:
		for (RotationSpec rotation : rotations) {
			for (int i=0; i<quantities.length; i++)
				if (!rotation.hasQuantity(quantities[i], values[i]))
					continue rotationLoop;
			ret.add(rotation);
		}
		return ret;
	}
	
	private static final boolean HYPO_NORTH = false;
	
	private static double angleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		return Math.abs(angleDiff);
	}
	
	private static final boolean D = false;
	
	private RSQSimEvent getInitialOrientation(RSQSimEvent rupture) {
		if (D) System.out.println("Initial orientation for "+rupture.getID());
		Location centroid = centroidCache.get(rupture);
		if (centroid == null) {
			centroid = RuptureRotationUtils.calcRuptureCentroid(rupture);
			centroidCache.putIfAbsent(rupture, centroid);
		}
		if (D) System.out.println("Initial centroid: "+centroid);
		
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
		
		if (HYPO_NORTH) {
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
	
	private static final double trans_p_diff_thresh = 0.5;
	private static final double trans_abs_diff_thresh = 0.2;
	private static final int min_translations = 2;
	private static final int max_translations = 100;
	
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
			centroidCache.putIfAbsent(rupture, centroid);
		}
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
					angleDiff = angleDiff(rupAngle, 0d);
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
				anns.add(RuptureRotationUtils.getLocationAnn(0.02, site.getLocation(), Color.BLUE));
				sites.add(site);
			}
			RSQSimEvent rotated = getRotatedRupture(rotation);
			if (first == null)
				first = rotated;
			plotElems.addAll(rotated.getAllElements());
		}
		
		if (highlightCentroid) {
			anns.add(RuptureRotationUtils.getLocationAnn(0.01, RuptureRotationUtils.calcRuptureCentroid(first), Color.GREEN));
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
	
	public RotatedRupVariabilityConfig forSites(List<Site> sites) {
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
	
	public RotatedRupVariabilityConfig forRotationSubset(List<RotationSpec> rotations) {
		return new RotatedRupVariabilityConfig(catalog, idToOrigMap == null ? null : idToOrigMap.values(), rotations);
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
		
		RotatedRupVariabilityConfig config;
		
		if (bbpDir == null) {
			List<BBP_Site> bbpSites = RSQSimBBP_Config.getCyberShakeInitialLASites();
			
			double[] distances = BBP_PartBValidationConfig.OFFICIAL_DISTANCES;
			int numSourceAz = 10;
			int numSiteToSourceAz = 10;
			
			List<Site> sites = new ArrayList<>();
			for (BBP_Site bbpSite : bbpSites)
				sites.add(bbpSite.buildGMPE_Site());
			
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
