package scratch.aftershockStatistics;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;

import scratch.aftershockStatistics.OAFTectonicRegime;

import java.util.Map;
import java.util.HashMap;

// OAFRegion is a region for defining a tectonic regime.
// It includes a Region which defines a polygonal area on the earth's surface,
// and minimum and maximum depths which define a depth range.
// It also includs the tectonic regime.

public class OAFRegion {

	// region - The region.

	private Region region;

	// min_depth - The minimum depth, in km (depth is positive down).

	private double min_depth;

	// max_depth - The maximum depth, in km (depth is positive down).

	private double max_depth;

	// MIN_DEPTH_UNBOUNDED - Value to use for minimum depth when there is no bound.

	public static final double MIN_DEPTH_UNBOUNDED = -1.0e10;

	// MAX_DEPTH_UNBOUNDED - Value to use for maximum depth when there is no bound.

	public static final double MAX_DEPTH_UNBOUNDED = 1.0e10;

	// regime - The tectonic regime.

	private OAFTectonicRegime regime;

	// get_regime - Gets the tectonic regime.

	public OAFTectonicRegime get_regime () {
		return regime;
	}

	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "OAFRegion\n"
			+ "\tRegime: " + get_regime() + "\n"
			+ "\tMinLat: " + region.getMinLat() + "\n"
			+ "\tMinLon: " + region.getMinLon() + "\n"
			+ "\tMaxLat: " + region.getMaxLat() + "\n"
			+ "\tMaxLon: " + region.getMaxLon() + "\n"
			+ "\tMinDepth: " + min_depth + "\n"
			+ "\tMaxDepth: " + max_depth;
		return str;
	}

	// Constructor saves the regime, region, and depth range.

	public OAFRegion (OAFTectonicRegime the_regime, Region the_region, double the_min_depth, double the_max_depth) {
		regime = the_regime;
		region = the_region;
		min_depth = the_min_depth;
		max_depth = the_max_depth;
	}

	// contains - Determine whether the given location is inside the region.
	// See notes for Region.contains.

	public boolean contains(Location loc) {
		return loc.getDepth() >= min_depth
			&& loc.getDepth() <= max_depth
			&& region.contains(loc);
	}

}
