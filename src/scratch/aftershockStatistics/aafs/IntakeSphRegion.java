package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;

// IntakeSphRegion is a region for defining earthquake intake into the AAFS server.
// It includes a SphRegion which defines a region of the earth's surface,
// and minimum magnitudes.
// It also includes the region name.
//
// Because magnitudes are often inaccurate when an event is first reported,
// there are two minimum magnitudes.  The intake_mag is applied when the first
// report is received, and selects events that will be considered for processing.
// The min_mag is applied when the first forecast is generated; the event is
// dropped if its magnitude is not at least min_mag at that time.

public class IntakeSphRegion {

	// name - The intake region name.  Cannot be null.

	private String name;

	// region - The region.

	private SphRegion region;

	// min_mag - The minimum magnitude required for an event to be processed by AAFS.

	private double min_mag;

	// intake_mag - The minimum magnitude required for an event to be considered.

	private double intake_mag;

	// Get the name.

	public String get_name() {
		return name;
	}

	// Get the region.

	public SphRegion get_region() {
		return region;
	}

	// Get the minimum magnitude.

	public double get_min_mag() {
		return min_mag;
	}

	// Get the intake magnitude.

	public double get_intake_mag() {
		return intake_mag;
	}

	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "IntakeSphRegion\n"
			+ "\tName: " + name + "\n"
			+ "\tMinLat: " + region.getMinLat() + "\n"
			+ "\tMinLon: " + region.getMinLon() + "\n"
			+ "\tMaxLat: " + region.getMaxLat() + "\n"
			+ "\tMaxLon: " + region.getMaxLon() + "\n"
			+ "\tMinMag: " + min_mag + "\n"
			+ "\tIntakeMag: " + intake_mag;
		return str;
	}

	// Constructor saves the name, region, and magnitudes.

	public IntakeSphRegion (String the_name, SphRegion the_region, double the_min_mag, double the_intake_mag) {
		name = the_name;
		region = the_region;
		min_mag = the_min_mag;
		intake_mag = the_intake_mag;
	}

	// contains - Determine whether the given location is inside the region,
	// and satisfies the minimum magnitude criterion.

	public boolean contains (SphLatLon loc, double mag) {
		return mag >= min_mag
			&& region.contains(loc);
	}

	// contains_intake - Determine whether the given location is inside the region,
	// and satisfies the intake magnitude criterion.

	public boolean contains_intake (SphLatLon loc, double mag) {
		return mag >= intake_mag
			&& region.contains(loc);
	}

}
