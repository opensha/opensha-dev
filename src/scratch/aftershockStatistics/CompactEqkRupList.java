package scratch.aftershockStatistics;

import java.util.AbstractList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;

import org.apache.commons.math3.distribution.UniformRealDistribution;

/**
 * Compact representation of earthquake rupture sequence.
 * Author: Michael Barall 03/29/2018.
 *
 * CompactEqkRupList is designed to store earthquake
 * sequences in a compact form that requires minimum storage.
 *
 * All the information about an earthquake is compressed into two 64-bit integers.
 * One contains the latitude, longitude, and depth to a resolution of 10 meters.
 * The other contains the magnitude in units of 0.001 and the time in milliseconds.
 *
 * The allowed data ranges and representations are:
 *
 * Latitude - Values from -90.0000 to 90.0000, in units of 0.0001 degree,
 * encoded in 21 bits as integers from -900000 to 900000 in excess 2^20 representation.
 *
 * Longitude - Values from -360.0000 to 360.0000, in units of 0.0001 degree,
 * encoded in 23 bits as integers from -3600000 to 3600000 in excess 2^22 representation.
 *
 * Depth - Values from -100.00 to 1000.00, in units of 0.01 kilometers,
 * encoded in 17 bits as integers from -10000 to 100000 in excess 2^12 representation.
 *
 * Magnitude - Values from -12.000 to 12.000, in units of 0.001 magnitude,
 * encoded in 15 bits as integers from -12000 to 12000 in excess 2^14 representation.
 *
 * Time - Values from -30000000000000 to 30000000000000, in units of 1 millisecond,
 * encoded in 46 bits in excess 2^45 representation.  The time range is approximately
 * +/- 951 years.  Negative times are allowed so that times relative to an arbitrary
 * origin are representable.
 *
 * Latitude, longitude, and depth are packed into a 64-bit integer as follows:
 * Bits 0-20 = latitude
 * Bits 21-43 = longitude
 * Bits 44-60 = depth
 * Bits 61-63 = zero
 *
 * Magnitude and time are packed into a 64-bit integer as follows:
 * Bits 0-14 = magnitude
 * Bits 15-60 = time
 * Bits 61-63 = zero
 *
 * The encoding scheme is such that neither 64-bit integer can be zero.
 * So, the value zero can be used to indicate no data.
 *
 * CompactEqkRupList is used to hold information about an earthquake sequence,
 * in compressed representation.  The data is stored in two arrays of long.
 */
public class CompactEqkRupList extends AbstractList<ObsEqkRupture> {

	//----- Encoding -----

	// Latitude encoding.

	public static final double LAT_MIN_VALUE = -90.0000;
	public static final double LAT_MAX_VALUE =  90.0000;
	public static final double LAT_SCALE =  10000.0;
	public static final long LAT_MASK = 0x1FFFFFL;
	public static final long LAT_SHIFT = 0L;
	public static final long LAT_OFFSET =  0x100000L;
	public static final long LAT_MIN_ENC = 0x100000L - 900000L;
	public static final long LAT_MAX_ENC = 0x100000L + 900000L;

	// Longitude encoding.

	public static final double LON_MIN_VALUE = -360.0000;
	public static final double LON_MAX_VALUE =  360.0000;
	public static final double LON_SCALE =  10000.0;
	public static final long LON_MASK = 0x7FFFFFL;
	public static final long LON_SHIFT = 21L;
	public static final long LON_OFFSET =  0x400000L;
	public static final long LON_MIN_ENC = 0x400000L - 3600000L;
	public static final long LON_MAX_ENC = 0x400000L + 3600000L;

	// Depth encoding.

	public static final double DEPTH_MIN_VALUE = -100.00;
	public static final double DEPTH_MAX_VALUE = 1000.00;
	public static final double DEPTH_SCALE =  100.0;
	public static final long DEPTH_MASK = 0x1FFFFL;
	public static final long DEPTH_SHIFT = 44L;
	public static final long DEPTH_OFFSET =  0x4000L;
	public static final long DEPTH_MIN_ENC = 0x4000L -  10000L;
	public static final long DEPTH_MAX_ENC = 0x4000L + 100000L;

	// Magnitude encoding.

	public static final double MAG_MIN_VALUE = -12.000;
	public static final double MAG_MAX_VALUE =  12.000;
	public static final double MAG_SCALE =  1000.0;
	public static final long MAG_MASK = 0x7FFFL;
	public static final long MAG_SHIFT = 0L;
	public static final long MAG_OFFSET =  0x4000L;
	public static final long MAG_MIN_ENC = 0x4000L - 12000L;
	public static final long MAG_MAX_ENC = 0x4000L + 12000L;

	// Time encoding.

	public static final long TIME_MIN_VALUE = -30000000000000L;
	public static final long TIME_MAX_VALUE =  30000000000000L;
	public static final long TIME_RESOLUTION = 1L;
	public static final long TIME_MASK = 0x3FFFFFFFFFFFL;
	public static final long TIME_SHIFT = 15L;
	public static final long TIME_OFFSET =  0x200000000000L;
	public static final long TIME_MIN_ENC = 0x200000000000L - 30000000000000L;
	public static final long TIME_MAX_ENC = 0x200000000000L + 30000000000000L;




	//----- Sequence storage -----

	// eqk_count - Number of earthquakes.

	private int eqk_count;

	// capacity - The allocated size of the arrays.

	private int capacity;

	// lat_lon_depth_list - Compressed latitude, longitude, and depth for each earthquake.
	// An entry may be zero if there is no location data for the particular earthquake.

	private long[] lat_lon_depth_list;

	// mag_time_list - Compressed magnitude and time for each earthquake.
	// An entry may be zero if there is no time and magnitude data for the particular earthquake.

	private long[] mag_time_list;




	//----- Getters -----

	// get_eqk_count - Number of earthquakes.

	public int get_eqk_count () {
		return eqk_count;
	}

	// get_capacity - The allocated size of the arrays.

	public int get_capacity () {
		return capacity;
	}

	// get_lat_lon_depth_list - Compressed latitude, longitude, and depth for each earthquake.

	public long[] get_lat_lon_depth_list () {
		return lat_lon_depth_list;
	}

	// get_mag_time_list - Compressed magnitude and time for each earthquake.

	public long[] get_mag_time_list () {
		return mag_time_list;
	}

	// as_ObsEqkRupList - Convert to ObsEqkRupList and return it.

	public ObsEqkRupList as_ObsEqkRupList () {
		ObsEqkRupList rups = new ObsEqkRupList();
		for (int index = 0; index < eqk_count; ++index) {
			rups.add (extract_rupture (null, lat_lon_depth_list[index], mag_time_list[index]));
		}
		return rups;
	}




	//----- Static compression functions -----

	// combine_lat_lon_depth - Combine latitude, longitude, and depth.

	public static long combine_lat_lon_depth (double lat, double lon, double depth) {

		long lat_enc = Math.round (lat * LAT_SCALE) + LAT_OFFSET;
		if (lat_enc < LAT_MIN_ENC || lat_enc > LAT_MAX_ENC) {
			throw new IllegalArgumentException("CompactEqkRupList.combine_lat_lon_depth: Latitude out-of-range: " + lat);
		}

		long lon_enc = Math.round (lon * LON_SCALE) + LON_OFFSET;
		if (lon_enc < LON_MIN_ENC || lon_enc > LON_MAX_ENC) {
			throw new IllegalArgumentException("CompactEqkRupList.combine_lat_lon_depth: Longitude out-of-range: " + lon);
		}

		long depth_enc = Math.round (depth * DEPTH_SCALE) + DEPTH_OFFSET;
		if (depth_enc < DEPTH_MIN_ENC || depth_enc > DEPTH_MAX_ENC) {
			throw new IllegalArgumentException("CompactEqkRupList.combine_lat_lon_depth: Depth out-of-range: " + depth);
		}

		return (depth_enc << DEPTH_SHIFT) | (lon_enc << LON_SHIFT) | lat_enc;
	}

	// combine_lat_lon_depth - Combine latitude, longitude, and depth.

	public static long combine_lat_lon_depth (Location loc) {

		if (loc == null) {
			return 0L;
		}

		return combine_lat_lon_depth (loc.getLatitude(), loc.getLongitude(), loc.getDepth());
	}

	// combine_lat_lon_depth - Combine latitude, longitude, and depth.

	public static long combine_lat_lon_depth (ObsEqkRupture rup) {
		return combine_lat_lon_depth (rup.getHypocenterLocation());
	}

	// extract_lat - Extract latitude.

	public static double extract_lat (long lat_lon_depth) {

		if (lat_lon_depth == 0L) {
			return 0.0;
		}

		return ( (double)((lat_lon_depth & LAT_MASK) - LAT_OFFSET) ) / LAT_SCALE;
	}

	// extract_lon - Extract longitude.

	public static double extract_lon (long lat_lon_depth) {

		if (lat_lon_depth == 0L) {
			return 0.0;
		}

		return ( (double)(((lat_lon_depth >> LON_SHIFT) & LON_MASK) - LON_OFFSET) ) / LON_SCALE;
	}

	// extract_depth - Extract depth.

	public static double extract_depth (long lat_lon_depth) {

		if (lat_lon_depth == 0L) {
			return 0.0;
		}

		return ( (double)(((lat_lon_depth >> DEPTH_SHIFT) & DEPTH_MASK) - DEPTH_OFFSET) ) / DEPTH_SCALE;
	}

	// extract_location - Extract location (returns null if lat_lon_depth is zero).

	public static Location extract_location (long lat_lon_depth) {

		if (lat_lon_depth == 0L) {
			return null;
		}

		double lat = ( (double)((lat_lon_depth & LAT_MASK) - LAT_OFFSET) ) / LAT_SCALE;
		double lon = ( (double)(((lat_lon_depth >> LON_SHIFT) & LON_MASK) - LON_OFFSET) ) / LON_SCALE;
		double depth = ( (double)(((lat_lon_depth >> DEPTH_SHIFT) & DEPTH_MASK) - DEPTH_OFFSET) ) / DEPTH_SCALE;

		return new Location(lat, lon, depth);
	}

	// combine_mag_time - Combine magnitude and time.

	public static long combine_mag_time (double mag, long time) {

		long mag_enc = Math.round (mag * MAG_SCALE) + MAG_OFFSET;
		if (mag_enc < MAG_MIN_ENC || mag_enc > MAG_MAX_ENC) {
			throw new IllegalArgumentException("CompactEqkRupList.combine_mag_time: Magnitude out-of-range: " + mag);
		}

		long time_enc = time + TIME_OFFSET;
		if (time_enc < TIME_MIN_ENC || time_enc > TIME_MAX_ENC) {
			throw new IllegalArgumentException("CompactEqkRupList.combine_mag_time: Time out-of-range: " + time);
		}

		return (time_enc << TIME_SHIFT) | mag_enc;
	}

	// combine_mag_time - Combine magnitude and time.

	public static long combine_mag_time (ObsEqkRupture rup) {
		return combine_mag_time (rup.getMag(), rup.getOriginTime());
	}

	// extract_mag - Extract magnitude.

	public static double extract_mag (long mag_time) {

		if (mag_time == 0L) {
			return 0.0;
		}

		return ( (double)((mag_time & MAG_MASK) - MAG_OFFSET) ) / MAG_SCALE;
	}

	// extract_time - Extract time.

	public static long extract_time (long mag_time) {

		if (mag_time == 0L) {
			return 0L;
		}

		return ((mag_time >> TIME_SHIFT) & TIME_MASK) - TIME_OFFSET;
	}

	// extract_rupture - Extract all values into a rupture.

	public static ObsEqkRupture extract_rupture (String eventId, long lat_lon_depth, long mag_time) {
		return new ObsEqkRupture(eventId, extract_time (mag_time), 
						extract_location (lat_lon_depth), extract_mag (mag_time));
	}




	//----- Implementation of AbstractList -----

	// size - Returns the number of elements in this collection.

	@Override
	public int size () {
		return eqk_count;
	}

	// get - Returns the element at the specified position in this list.

	@Override
	public ObsEqkRupture get (int index) {
		if (!( index >= 0 && index < eqk_count )) {
			throw new IndexOutOfBoundsException("CompactEqkRupList.get: Invalid index: index = " + index + ", size = " + eqk_count);
		}
		return extract_rupture (null, lat_lon_depth_list[index], mag_time_list[index]);
	}

	// set - Replaces the element at the specified position in this list with the specified element (optional operation).
	// Returns the element previously at the specified position.

	@Override
	public ObsEqkRupture set (int index, ObsEqkRupture element) {
		ObsEqkRupture previous = get (index);
		lat_lon_depth_list[index] = combine_lat_lon_depth (element);
		mag_time_list[index] = combine_mag_time (element);
		return previous;
	}

	// set_only - Replaces the element at the specified position in this list with the specified element (optional operation).
	// This version (not part of AbstractList) does not bother to compute the previous element.

	public void set_only (int index, ObsEqkRupture element) {
		if (!( index >= 0 && index < eqk_count )) {
			throw new IndexOutOfBoundsException("CompactEqkRupList.set_only: Invalid index: index = " + index + ", size = " + eqk_count);
		}
		lat_lon_depth_list[index] = combine_lat_lon_depth (element);
		mag_time_list[index] = combine_mag_time (element);
		return;
	}

	// add - Inserts the specified element at the specified position in this list (optional operation). 
	// Shifts the element currently at that position (if any) and any subsequent elements to the right (adds one to their indices).

	@Override
	public void add (int index, ObsEqkRupture element) {
		++modCount;
		insert_elements (index, 1);
		set_only (index, element);
		return;
	}

	// add - Appends the specified element to the end of this list (optional operation).
	// Always returns true.

	@Override
	public boolean add (ObsEqkRupture element) {
		++modCount;
		if (eqk_count == capacity) {
			increase_capacity (1);
		}
		++eqk_count;
		set_only (eqk_count - 1, element);
		return true;
	}

	// removeRange - Removes from this list all of the elements whose index is between fromIndex, inclusive, and toIndex, exclusive.
	// Shifts any succeeding elements to the left (reduces their index).

	@Override
	protected void removeRange (int fromIndex, int toIndex) {
		if (fromIndex != toIndex) {
			++modCount;
		}
		remove_elements (fromIndex, toIndex - fromIndex);
		return;
	}

	// remove - Removes the element at the specified position in this list (optional operation).
	// Shifts any subsequent elements to the left (subtracts one from their indices).
	// Returns the element that was removed from the list. 

	@Override
	public ObsEqkRupture remove (int index) {
		++modCount;
		ObsEqkRupture previous = get (index);
		remove_elements (index, 1);
		return previous;
	}

	// remove_only - Removes the element at the specified position in this list (optional operation).
	// Shifts any subsequent elements to the left (subtracts one from their indices).
	// This version (not part of AbstractList) does not bother to compute the previous element.

	public void remove_only (int index) {
		++modCount;
		remove_elements (index, 1);
		return;
	}

	// addAll - Inserts all of the elements in the specified collection into this list at the specified position (optional operation).
	// Shifts the element currently at that position (if any) and any subsequent elements to the right (increases their indices).
	// The new elements will appear in this list in the order that they are returned by the specified collection's iterator.
	// Returns true if this list changed as a result of the call.

	@Override
	public boolean addAll (int index, Collection<? extends ObsEqkRupture> c) {

		// Get number of elements to add

		int len = c.size();

		// Nothing to do if collection is empty

		if (len == 0) {
			return false;
		}
		++modCount;

		// Create the space

		insert_elements (index, len);

		// Iterate over collection, inserting each element

		int n = 0;
		for (ObsEqkRupture rup : c) {
			if (n >= len) {
				throw new RuntimeException("CompactEqkRupList.addAll: Collection size larger than expected");
			}
			set_only (index + n, rup);
			++n;
		}
		if (n != len) {
			throw new RuntimeException("CompactEqkRupList.addAll: Collection size smaller than expected");
		}

		return true;
	}




	//----- Storage management -----

	// increase_capacity - Add more capacity to the list.

	private void increase_capacity (int amount) {

		// Check for maximum capacity

		if (Integer.MAX_VALUE - capacity < amount) {
			throw new RuntimeException("CompactEqkRupList.increase_capacity: Exceeded maximum capacity");
		}

		// Double capacity, but at least 100, and not more than Integer.MAX_VALUE

		int new_capacity = Math.max (100, capacity + Math.min (Math.max (capacity, amount), Integer.MAX_VALUE - capacity));

		// Reallocate the arrays, with zero padding

		lat_lon_depth_list = Arrays.copyOf (lat_lon_depth_list, new_capacity);
		mag_time_list = Arrays.copyOf (mag_time_list, new_capacity);
		capacity = new_capacity;
		return;
	}

	// alloc_storage - Initial allocation of storage.

	private void alloc_storage (int initial_capacity) {

		// Initial size of arrays is initial_capacity

		lat_lon_depth_list = new long[initial_capacity];
		mag_time_list = new long[initial_capacity];
		capacity = initial_capacity;
		eqk_count = 0;

		// Zero-initialize the arrays (strictly this is unnecessary, but I don't like to rely on implicit initialization)

		Arrays.fill (lat_lon_depth_list, 0L);
		Arrays.fill (mag_time_list, 0L);
		return;
	}

	// insert_elements - Insert elements in the arrays, beginning at the given index, for the given length.

	private void insert_elements (int index, int len) {
		if (!( index >= 0 && len >= 0 && index <= eqk_count )) {
			throw new IndexOutOfBoundsException("CompactEqkRupList.insert_elements: Invalid index or length: index = " + index + ", len = " + len + ", size = " + eqk_count);
		}

		// Nothing to do if length is zero

		if (len != 0) {

			// Increase capacity if needed

			if (capacity - eqk_count < len) {
				increase_capacity (len - (capacity - eqk_count));
			}

			// If the space we need is internal, shift the array contents

			if (index < eqk_count) {
				System.arraycopy (lat_lon_depth_list, index, lat_lon_depth_list, index + len, eqk_count - index);
				System.arraycopy (mag_time_list, index, mag_time_list, index + len, eqk_count - index);

				// Zero-fill the new elements (if index == eqk_count then they're already zero,
				// and there's no need to fill past the original end of the array because it's already zero)

				Arrays.fill (lat_lon_depth_list, index, Math.min (eqk_count, index + len), 0L);
				Arrays.fill (mag_time_list, index, Math.min (eqk_count, index + len), 0L);
			}

			// Adjust the earthquake count

			eqk_count += len;
		}
		return;
	}

	// remove_elements - Remove elements in the arrays, beginning at the given index, for the given length.

	private void remove_elements (int index, int len) {
		if (!( index >= 0 && len >= 0 && index + len <= eqk_count )) {
			throw new IndexOutOfBoundsException("CompactEqkRupList.remove_elements: Invalid index or length: index = " + index + ", len = " + len + ", size = " + eqk_count);
		}

		// Nothing to do if length is zero

		if (len != 0) {

			// If the elements being remove are internal, shift the array contents

			if (index + len < eqk_count) {
				System.arraycopy (lat_lon_depth_list, index + len, lat_lon_depth_list, index, eqk_count - (index + len));
				System.arraycopy (mag_time_list, index + len, mag_time_list, index, eqk_count - (index + len));
			}

			// Zero-fill the vacated space

			Arrays.fill (lat_lon_depth_list, eqk_count - len, eqk_count, 0L);
			Arrays.fill (mag_time_list, eqk_count - len, eqk_count, 0L);

			// Adjust the earthquake count

			eqk_count -= len;
		}
		return;
	}




	//----- Construction -----

	// Create an empty list with default capacity.

	public CompactEqkRupList () {
		alloc_storage (100);
	}

	// Create an empty list with specified initial capacity.

	public CompactEqkRupList (int initial_capacity) {
		alloc_storage (initial_capacity);
	}

	// Create a list that is a copy of another list.
	// Note: The resulting arrays have length exactly equal to c.size().

	public CompactEqkRupList (Collection<? extends ObsEqkRupture> c) {
		alloc_storage (c.size());
		addAll (0, c);
	}

	// Create a list with the given arrays.
	// Note: The supplied arrays are stored in this object (not copied).
	// The arrays must have the same length, and that length becomes the capacity.
	// If longer than eqk_count, the arrays must be zero-padded if the list is
	// considered modifiable (not necessary if the list is treated as read-only).
	// This is intended primarily for reading lists stored in the OAF database.

	public CompactEqkRupList (int eqk_count, long[] lat_lon_depth_list, long[] mag_time_list) {
		if (!( eqk_count >= 0
			&& lat_lon_depth_list.length == mag_time_list.length
			&& eqk_count <= lat_lon_depth_list.length )) {
			throw new IndexOutOfBoundsException("CompactEqkRupList.CompactEqkRupList: Inconsistent lengths");
		}

		this.eqk_count = eqk_count;
		this.capacity = lat_lon_depth_list.length;
		this.lat_lon_depth_list = lat_lon_depth_list;
		this.mag_time_list = mag_time_list;
	}




	//----- Testing -----




	// test_make_random_rupture - Make a random earthquake rupture.

	private static ObsEqkRupture test_make_random_rupture (UniformRealDistribution rangen) {

		// Bounds

		double lat_min = -90.0;
		double lat_max = 90.0;

		double lon_min = -180.0;
		double lon_max = 360.0;

		double depth_min = -5.0;
		double depth_max = 700.0;

		double mag_min = MAG_MIN_VALUE;
		double mag_max = MAG_MAX_VALUE;

		long time_min = TIME_MIN_VALUE;
		long time_max = TIME_MAX_VALUE;

		// Random values

		double lat = rangen.sample() * (lat_max - lat_min) + lat_min;
		double lon = rangen.sample() * (lon_max - lon_min) + lon_min;
		double depth = rangen.sample() * (depth_max - depth_min) + depth_min;
		double mag = rangen.sample() * (mag_max - mag_min) + mag_min;
		long time = Math.round(rangen.sample() * ((double)time_max - (double)time_min)) + time_min;

		// One time in 100, use a null Location

		if (rangen.sample() < 0.01) {
			return new ObsEqkRupture (null, time, null, mag);
		}

		// Return the rupture
	
		return new ObsEqkRupture (null, time, new Location (lat, lon, depth), mag);
	}




	// test_compare_ruptures - Compare two ruptures, within tolerance, return true if match.

	private static boolean test_compare_ruptures (ObsEqkRupture rup1, ObsEqkRupture rup2, boolean f_verbose) {
	
		// Tolerances

		double lat_tol = 1.0/LAT_SCALE;

		double lon_tol = 1.0/LON_SCALE;

		double depth_tol = 1.0/DEPTH_SCALE;

		double mag_tol = 1.0/MAG_SCALE;

		long time_tol = TIME_RESOLUTION;

		// Values for rup1

		Location r_loc = rup1.getHypocenterLocation();
		double r_lat = ((r_loc == null) ? -999.0 : r_loc.getLatitude());
		double r_lon = ((r_loc == null) ? -999.0 : r_loc.getLongitude());
		double r_depth = ((r_loc == null) ? -999.0 : r_loc.getDepth());
		double r_mag = rup1.getMag();
		long r_time = rup1.getOriginTime();

		// Values for rup2

		Location s_loc = rup2.getHypocenterLocation();
		double s_lat = ((s_loc == null) ? -999.0 : s_loc.getLatitude());
		double s_lon = ((s_loc == null) ? -999.0 : s_loc.getLongitude());
		double s_depth = ((s_loc == null) ? -999.0 : s_loc.getDepth());
		double s_mag = rup2.getMag();
		long s_time = rup2.getOriginTime();

		// Comparison

		if (!( Math.abs (r_lat - s_lat) < lat_tol
			&& Math.abs (r_lon - s_lon) < lon_tol
			&& Math.abs (r_depth - s_depth) < depth_tol
			&& Math.abs (r_mag - s_mag) < mag_tol
			&& Math.abs (r_time - s_time) < time_tol )) {

			if (f_verbose) {
				System.out.println ("Mismatch:\n" +
					"lat: " + r_lat + " - " + s_lat + " = " + (r_lat - s_lat) + "\n" +
					"lon: " + r_lon + " - " + s_lon + " = " + (r_lon - s_lon) + "\n" +
					"depth: " + r_depth + " - " + s_depth + " = " + (r_depth - s_depth) + "\n" +
					"mag: " + r_mag + " - " + s_mag + " = " + (r_mag - s_mag) + "\n" +
					"time: " + r_time + " - " + s_time + " = " + (r_time - s_time)
				);
			}

			return false;
		}
	
		return true;
	}




	// test_compare_sequences - Compare two rupture sequences, within tolerance, return true if match.

	private static boolean test_compare_sequences (List<ObsEqkRupture> seq1, List<ObsEqkRupture> seq2, boolean f_verbose) {

		// Size

		int n1 = seq1.size();
		int n2 = seq2.size();

		if (n1 != n2) {
			if (f_verbose) {
				System.out.println ("Size mismatch: " + n1 + " versus " + n2);
			}
			return false;
		}

		// Compare contents

		int errors = 0;
		for (int index = 0; index < n1; ++index) {
			if (!( test_compare_ruptures (seq1.get(index), seq2.get(index), f_verbose && errors < 10) )) {
				++errors;
			}
		}

		if (f_verbose) {
			System.out.println ("Error count: " + errors);
		}

		if (errors != 0) {
			return false;
		}

		return true;
	}




	// test_random_index - Return a random integer from 0 to n-1.

	private static int test_random_index (UniformRealDistribution rangen, int n) {
		int index = (int)Math.floor (rangen.sample() * (double)n);
		if (index < 0) {
			index = 0;
		}
		if (index >= n) {
			index = n-1;
		}
		return index;
	}




	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("CompactEqkRupList : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  num_events
		// First, construct a ObsEqkRupList of size num_events, with random contents.
		// Then, convert it to CompactEqkRupList.
		// Finally, iterate through the CompactEqkRupList and test that the original
		// random contents is recovered, within tolerances.
		// Note: It's OK to use num_events values up to 1 million, but larger values
		// may crash the JVM due to large memory usage.

		if (args[0].equalsIgnoreCase ("test1")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("CompactEqkRupList : Invalid 'test1' subcommand");
				return;
			}
			int num_events = Integer.parseInt(args[1]);

			// Bounds and tolerances

			double lat_min = -90.0;
			double lat_max = 90.0;
			double lat_tol = 1.0/CompactEqkRupList.LAT_SCALE;

			double lon_min = -180.0;
			double lon_max = 360.0;
			double lon_tol = 1.0/CompactEqkRupList.LON_SCALE;

			double depth_min = -5.0;
			double depth_max = 700.0;
			double depth_tol = 1.0/CompactEqkRupList.DEPTH_SCALE;

			double mag_min = CompactEqkRupList.MAG_MIN_VALUE;
			double mag_max = CompactEqkRupList.MAG_MAX_VALUE;
			double mag_tol = 1.0/CompactEqkRupList.MAG_SCALE;

			long time_min = CompactEqkRupList.TIME_MIN_VALUE;
			long time_max = CompactEqkRupList.TIME_MAX_VALUE;
			long time_tol = CompactEqkRupList.TIME_RESOLUTION;

			System.out.println (
				"lat_min = " + lat_min + "\n" +
				"lat_max = " + lat_max + "\n" +
				"lat_tol = " + lat_tol + "\n" +
				"lon_min = " + lon_min + "\n" +
				"lon_max = " + lon_max + "\n" +
				"lon_tol = " + lon_tol + "\n" +
				"depth_min = " + depth_min + "\n" +
				"depth_max = " + depth_max + "\n" +
				"depth_tol = " + depth_tol + "\n" +
				"mag_min = " + mag_min + "\n" +
				"mag_max = " + mag_max + "\n" +
				"mag_tol = " + mag_tol + "\n" +
				"time_min = " + time_min + "\n" +
				"time_max = " + time_max + "\n" +
				"time_tol = " + time_tol + "\n" +
				"num_events = " + num_events
			);

			// Random number generator

			UniformRealDistribution rangen = new UniformRealDistribution();

			// Generate random values

			System.out.println ("Generating random data ...");

			int n;

			double[] lat = new double[num_events];
			for (n = 0; n < num_events; ++n) {
				lat[n] = rangen.sample() * (lat_max - lat_min) + lat_min;
			}

			double[] lon = new double[num_events];
			for (n = 0; n < num_events; ++n) {
				lon[n] = rangen.sample() * (lon_max - lon_min) + lon_min;
			}

			double[] depth = new double[num_events];
			for (n = 0; n < num_events; ++n) {
				depth[n] = rangen.sample() * (depth_max - depth_min) + depth_min;
			}

			double[] mag = new double[num_events];
			for (n = 0; n < num_events; ++n) {
				mag[n] = rangen.sample() * (mag_max - mag_min) + mag_min;
			}

			long[] time = new long[num_events];
			for (n = 0; n < num_events; ++n) {
				time[n] = Math.round(rangen.sample() * ((double)time_max - (double)time_min)) + time_min;
			}

			// Create the rupture list

			System.out.println ("Creating rupture list ...");

			ObsEqkRupList rup_list = new ObsEqkRupList();

			for (n = 0; n < num_events; ++n) {
				rup_list.add (new ObsEqkRupture (
					null,
					time[n],
					new Location (lat[n], lon[n], depth[n]),
					mag[n]
				));
			}

			// Compress the rupture list

			System.out.println ("Compressing rupture list ...");

			CompactEqkRupList compact_list = new CompactEqkRupList (rup_list);

			// Compare results

			System.out.println ("Comparing results ...");

			n = 0;
			int errors = 0;

			for (ObsEqkRupture rup : compact_list) {
				double r_lat = rup.getHypocenterLocation().getLatitude();
				double r_lon = rup.getHypocenterLocation().getLongitude();
				double r_depth = rup.getHypocenterLocation().getDepth();
				double r_mag = rup.getMag();
				long r_time = rup.getOriginTime();

				if (!( Math.abs (r_lat - lat[n]) < lat_tol
					&& Math.abs (r_lon - lon[n]) < lon_tol
					&& Math.abs (r_depth - depth[n]) < depth_tol
					&& Math.abs (r_mag - mag[n]) < mag_tol
					&& Math.abs (r_time - time[n]) < time_tol )) {

					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatch:\n" +
							"lat: " + r_lat + " - " + lat[n] + " = " + (r_lat - lat[n]) + "\n" +
							"lon: " + r_lon + " - " + lon[n] + " = " + (r_lon - lon[n]) + "\n" +
							"depth: " + r_depth + " - " + depth[n] + " = " + (r_depth - depth[n]) + "\n" +
							"mag: " + r_mag + " - " + mag[n] + " = " + (r_mag - mag[n]) + "\n" +
							"time: " + r_time + " - " + time[n] + " = " + (r_time - time[n])
						);
					}
				}

				++n;
			}

			System.out.println ("Error count: " + errors);

			return;
		}




		// Subcommand : Test #2
		// Command format:
		//  test1  num_events
		// Construct ObsEqkRupList and CompactEqkRupList sequences, and perform the
		// same series of add, remove, get, and set actions on each.  Testing is
		// done by comparing the two sequences.  Random ruptures are used.
		// Note: The maximum sequence size is num_events*4.  Setting num_events to 500
		// will force the storage to be reallocated 4 times.  Execution time is quadratic
		// in num_events.  (The addAll test increases list size to about num_events*10.)

		if (args[0].equalsIgnoreCase ("test2")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("CompactEqkRupList : Invalid 'test2' subcommand");
				return;
			}
			int num_events = Integer.parseInt(args[1]);

			// Random number generator

			UniformRealDistribution rangen = new UniformRealDistribution();

			// Quick test of the comparitor.

			ObsEqkRupList seq1 = new ObsEqkRupList();
			CompactEqkRupList seq2 = new CompactEqkRupList();

			int n;
			int index;

			ObsEqkRupture rup;
			ObsEqkRupture rup1;
			ObsEqkRupture rup2;

			rup = test_make_random_rupture (rangen);
			seq1.add (rup);
			seq2.add (rup);

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			rup = test_make_random_rupture (rangen);
			seq1.add (rup);
			rup = test_make_random_rupture (rangen);
			seq2.add (rup);

			if (test_compare_sequences (seq1, seq2, false)) {
				System.out.println ("Comparitor failed to detect difference.");
				return;
			}

			// Clear function

			System.out.println ("Testing clear ...");

			seq1.clear();
			seq2.clear();

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Add internal

			System.out.println ("Testing add internal ...");

			for (n = 0; n < num_events; ++n) {
				index = test_random_index (rangen, seq1.size() + 1);
				rup = test_make_random_rupture (rangen);
				seq1.add (index, rup);
				seq2.add (index, rup);
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Add at end

			System.out.println ("Testing add at end ...");

			for (n = 0; n < num_events * 3; ++n) {
				rup = test_make_random_rupture (rangen);
				seq1.add (rup);
				seq2.add (rup);
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Set

			System.out.println ("Testing set ...");

			for (n = 0; n < num_events; ++n) {
				if (n == 0) {
					index = 0;
				} else if (n == 1) {
					index = seq1.size() - 1;
				} else {
					index = test_random_index (rangen, seq1.size());
				}
				rup = test_make_random_rupture (rangen);
				rup1 = seq1.set (index, rup);
				rup2 = seq2.set (index, rup);

				if (!( test_compare_ruptures (rup1, rup2, true) )) {
					return;
				}
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Remove

			System.out.println ("Testing remove ...");

			for (n = 0; n < num_events; ++n) {
				if (n == 0) {
					index = 0;
				} else if (n == 1) {
					index = seq1.size() - 1;
				} else {
					index = test_random_index (rangen, seq1.size());
				}
				rup1 = seq1.remove (index);
				rup2 = seq2.remove (index);

				if (!( test_compare_ruptures (rup1, rup2, true) )) {
					return;
				}
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Add internal

			System.out.println ("Testing add internal ...");

			for (n = 0; n < num_events; ++n) {
				if (n == 0) {
					index = 0;
				} else if (n == 1) {
					index = seq1.size();
				} else {
					index = test_random_index (rangen, seq1.size() + 1);
				}
				rup = test_make_random_rupture (rangen);
				seq1.add (index, rup);
				seq2.add (index, rup);
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Add all

			System.out.println ("Testing add all ...");

			for (n = 0; n < num_events; ++n) {
				if (n == 0) {
					index = 0;
				} else if (n == 1) {
					index = seq1.size() - 1;
				} else {
					index = test_random_index (rangen, seq1.size());
				}
				ArrayList<ObsEqkRupture> c = new ArrayList<ObsEqkRupture>();
				int count = test_random_index (rangen, 13);
				for (int i = 0; i < count; ++i) {
					c.add (test_make_random_rupture (rangen));
				}
				seq1.addAll (index, c);
				seq2.addAll (index, c);
			}

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Clear

			System.out.println ("Testing clear ...");

			seq1.clear();
			seq2.clear();

			if (!( test_compare_sequences (seq1, seq2, true) )) {
				return;
			}

			// Done

			System.out.println ("All test completed successfully");

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("CompactEqkRupList : Unrecognized subcommand : " + args[0]);
		return;

	}





}
