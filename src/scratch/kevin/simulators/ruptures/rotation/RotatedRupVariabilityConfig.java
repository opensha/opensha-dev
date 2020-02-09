package scratch.kevin.simulators.ruptures.rotation;

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
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.NamedComparator;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;

public abstract class RotatedRupVariabilityConfig<E> {

	private List<RotationSpec> rotations;
	private Map<Quantity, List<?>> quantitiesMap;
	private Table<Quantity, Object, List<RotationSpec>> quantityRotationsCache;
	
	// caches
	private static final int max_cache_size = 5000;
	protected LoadingCache<RotationSpec, E> rotationCache;
	private LoadingCache<Integer, E> initialOrientationCache;
	private Map<E, Location> centroidCache;
	
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
	
	@SuppressWarnings("unchecked")
	public RotatedRupVariabilityConfig(List<RotationSpec> rotations) {
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
		
		rotationCache = CacheBuilder.newBuilder().maximumSize(max_cache_size).build(
				new CacheLoader<RotationSpec, E>() {

					@Override
					public E load(RotationSpec key) throws Exception {
						// TODO Auto-generated method stub
						return loadRotatedRupture(key);
					}
			
		});
		initialOrientationCache = CacheBuilder.newBuilder().build(
				new CacheLoader<Integer, E>() {

					@Override
					public E load(Integer key) throws Exception {
						return loadInitialOrientationRupture(key);
					}
			
		});
		centroidCache = new HashMap<>();
	}
	
	public abstract boolean hasRuptures();

	protected abstract E loadRotatedRupture(RotationSpec key);
	protected abstract E loadInitialOrientationRupture(Integer eventID);
	protected abstract Location loadCentroid(E rupture);
	
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
			if (value instanceof Collection<?>) {
				for (Object subValue : ((Collection<?>)value))
					if (hasQuantity(quantity, subValue))
						return true;
				return false;
			}
//			if (quantity.zeroStoredAsNull() && Objects.equals(value, 0f))
//				value = null;
//			if (quantity == Quantity.SITE) {
//				Site mySite = (Site)quantities.get(quantity);
//				if (value == null)
//					return mySite == null;
//				else if (mySite == null)
//					return false;
//				return locEqual(mySite.getLocation(), ((Site)value).getLocation());
//			}
//			return Objects.equals(value, quantities.get(quantity));
			Object myVal = quantities.get(quantity);
			if (myVal == null && quantity.zeroStoredAsNull())
				return value == null || (Float)value == 0f;
			if (quantity == Quantity.SITE) {
				Site mySite = (Site)myVal;
				if (value == null)
					return mySite == null;
				else if (mySite == null)
					return false;
				return locEqual(mySite.getLocation(), ((Site)value).getLocation());
			}
			if (value instanceof Integer)
				return ((Integer)value).intValue() == (Integer)myVal;
			if (value instanceof Float)
				return ((Float)value).floatValue() == (Float)myVal;
			return Objects.equals(value, myVal);
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
	
	protected static List<RotationSpec> buildRotations(List<Site> sites, List<Integer> eventIDs, double[] distances,
			int numSourceAz, int numSiteToSourceAz) {
		List<RotationSpec> rotations = new ArrayList<>();
		int index = 0;
		double sourceDeltaAz = 360/(double)numSourceAz;
		double siteToSourceDeltaAz = 360/(double)numSiteToSourceAz;
		for (Site site : sites) {
			for (double distance : distances) {
				for (int eventID : eventIDs) {
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
		if (value instanceof Collection<?>) {
			List<RotationSpec> ret = new ArrayList<>();
			for (Object subValue : (Collection<?>)value)
				ret.addAll(getCachedRotations(quantity, subValue));
			return ret;
		}
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
		Predicate<RotationSpec> predicate = new Predicate<RotatedRupVariabilityConfig.RotationSpec>() {
			
			@Override
			public boolean test(RotationSpec rotation) {
				for (int i=0; i<quantities.length; i++)
					if (!rotation.hasQuantity(quantities[i], values[i]))
						return false;
				return true;
			}
		};
		
//		return rotations.parallelStream().filter(predicate).collect(Collectors.toList());
		return rotations.stream().filter(predicate).collect(Collectors.toList());
		
//		List<RotationSpec> ret = new ArrayList<>();
//		rotationLoop:
//		for (RotationSpec rotation : rotations) {
//			for (int i=0; i<quantities.length; i++)
//				if (!rotation.hasQuantity(quantities[i], values[i]))
//					continue rotationLoop;
//			ret.add(rotation);
//		}
//		return ret;
	}
	
	protected static final boolean HYPO_NORTH = false;
	
	protected static final boolean D = false;
	
	public Location getCentroid(E rupture) {
		Location centroid = centroidCache.get(rupture);
		if (centroid == null) {
			centroid = loadCentroid(rupture);
			centroidCache.putIfAbsent(rupture, centroid);
		}
		return centroid;
	}
	
	public E getInitiallyOrientedRupture(Integer eventID) {
		try {
			return initialOrientationCache.get(eventID);
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	public synchronized E getRotatedRupture(RotationSpec rotation) {
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
	
	public abstract void plotRotations(File outputDir, String prefix, List<RotationSpec> rotations, boolean highlightCentroid)
			throws IOException;

}
