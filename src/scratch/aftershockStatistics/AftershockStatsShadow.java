package scratch.aftershockStatistics;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;

import static org.opensha.commons.geo.GeoTools.TO_DEG;
import static org.opensha.commons.geo.GeoTools.TO_RAD;

import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;
import scratch.aftershockStatistics.util.SphRegionCircle;
import scratch.aftershockStatistics.util.SimpleUtils;



/**
 * Calculate aftershock event shadowing.
 * Author: Michael Barall 07/26/2018.
 *
 * It is generally not desireable to generate aftershock forecasts for an
 * earthquake that is an aftershock or foreshock of another, larger earthquake.
 * Such an earthquake is said to be shadowed by the larger earthquake.
 *
 * This class contains the machinery to determine if an earthquake is shadowed.
 */
public class AftershockStatsShadow {




	// This class holds a candidate shadowing event.

	public static class CandidateShadow {
	
		// The rupture.

		public ObsEqkRupture rupture;

		// The centroid radius, in kilometers.
		// Events that contribute to the centroid must lie within this distance of the hypocenter.

		public double centroid_radius;

		// The minimum magnitude for events that contribute to the centroid.

		public double centroid_min_mag;

		// The lower and upper limits of origin time for events that contribute to the centroid.
		// These are times in milliseconds since the epoch.

		public long centroid_time_lo;
		public long centroid_time_hi;

		// The The sample radius, in kilometers.
		// An event is considered to be an aftershock if it lies within this distance of the centroid.

		public double sample_radius;

		// Accumulators for vectors that define the centroid.

		public double x;
		public double y;
		public double z;

		// Computed centroid

		public Location candidate_centroid;

		// Candidate parameters: event ID, origin time, magnitude, hypocenter.

		String candidate_event_id;
		long candidate_time;
		double candidate_mag;
		Location candidate_hypo;

		// Constructor.

		public CandidateShadow (
			ObsEqkRupture rupture,
			double centroid_radius,
			double centroid_min_mag,
			long centroid_time_lo,
			long centroid_time_hi,
			double sample_radius) {
			
			this.rupture = rupture;
			this.centroid_radius = centroid_radius;
			this.centroid_min_mag = centroid_min_mag;
			this.centroid_time_lo = centroid_time_lo;
			this.centroid_time_hi = centroid_time_hi;
			this.sample_radius = sample_radius;

			this.x = 0.0;
			this.y = 0.0;
			this.z = 0.0;

			this.candidate_centroid = null;

			this.candidate_event_id = rupture.getEventId();
			this.candidate_time = rupture.getOriginTime();
			this.candidate_mag = rupture.getMag();
			this.candidate_hypo = rupture.getHypocenterLocation();

			// Accumulate the unit vector for the rupture itself (see AftershockStatsCalc.getSphCentroid)

			double lat = this.candidate_hypo.getLatRad();
			double lon = this.candidate_hypo.getLonRad();

			this.x += (Math.cos(lat) * Math.cos(lon));
			this.y += (Math.cos(lat) * Math.sin(lon));
			this.z += Math.sin(lat);
		}

		// Call Comcat to obtain a list of possible aftershocks, and accumulate their unit vectors.
		// An exception from this function likely indicates a problem with Comcat.

		public void accum_from_comcat (ComcatAccessor accessor, long system_time_now) {

			// If there is a nonempty time and space region in which we need to look for aftershocks ...

			if (centroid_time_lo < system_time_now
				&& centroid_time_lo < centroid_time_hi
				&& centroid_radius > 0.0) {

				// Construct a circle around the rupture hypocenter with the centroid radius

				SphLatLon candidate_sph_hypo = new SphLatLon (candidate_hypo);
				SphRegion centroid_region = SphRegion.makeCircle (candidate_sph_hypo, centroid_radius);

				// Get a list of possible aftershocks by calling Comcat
				// Aftershocks must lie in the centroid region, within the centroid times,
				// and have magnitude at least equal to the centroid magnitude

				double min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;
				double max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;

				boolean wrapLon = false;
				int limit_per_call = 0;
				int max_calls = 0;

				ObsEqkRupList aftershocks = accessor.fetchEventList (candidate_event_id,
							centroid_time_lo, centroid_time_hi,
							min_depth, max_depth,
							centroid_region, wrapLon,
							centroid_min_mag, limit_per_call, max_calls);

				System.out.println ("AftershockStatsShadow.accum_from_comcat: Found " + aftershocks.size() + " aftershocks within " + String.format ("%.3f", centroid_radius) + " km of candidate event " + candidate_event_id);
		
				// For each aftershock ...

				for (ObsEqkRupture aftershock : aftershocks) {

					// Get the aftershock parameters

					Location aftershock_hypo = aftershock.getHypocenterLocation();

					// Add unit vector for centroid calculation (see AftershockStatsCalc.getSphCentroid)

					double lat = aftershock_hypo.getLatRad();
					double lon = aftershock_hypo.getLonRad();

					x += (Math.cos(lat) * Math.cos(lon));
					y += (Math.cos(lat) * Math.sin(lon));
					z += Math.sin(lat);
				}
			}

			return;
		}

		// Get the centroid, as determined by the accumulators (see AftershockStatsCalc.getSphCentroid).
		// Also save the result into candidate_centroid.

		public Location get_centroid () {

			double lat;
			double lon;

			// If the vector is very small, just return the candidate location

			if (x*x + y*y + z*z < 1.0e-4) {
				lat = candidate_hypo.getLatitude();
				lon = candidate_hypo.getLongitude();
				if (lon > 180.0) {
					lon -= 360.0;
				}
			}

			// Otherwise, convert rectangular to spherical coordinates

			else {
				lat = Math.atan2 (z, Math.hypot(x, y)) * TO_DEG;
				lon = Math.atan2 (y, x) * TO_DEG;
			}

			// Make sure the angles are in range, since they were converted from radians

			if (lat > 90.0) {lat = 90.0;}
			if (lat < -90.0) {lat = -90.0;}
			if (lon > 180.0) {lon = 180.0;}
			if (lon < -180.0) {lon = -180.0;}

			// Centroid

			candidate_centroid = new Location (lat, lon);
			return candidate_centroid;
		}
	}




	// Default parameter values for find_shadow.

	public static final long YEAR_IN_MILLIS = 31536000000L;		// 1 year = 365 days

	public static final double DEF_SEARCH_RADIUS = 2000.0;		// default search radius = 2000 km

	public static final double DEF_CENTROID_MAG_FLOOR = 2.5;	// default centroid magnitude floor

	public static final double DEF_LARGE_MAG = 8.0;				// default large magnitude

	


	// Determine if a given mainshock is shadowed.
	// Parameters:
	//  mainshock = The mainshock to check for shadowing.
	//  time_now = The time at which the check is to be done.  Typically time_now is close
	//    to the current time.  But it can also be in the past, to determine if the
	//    mainshock was shadowed at a past time.
	//  search_radius = The radius of a circle surrounding the mainshock, in km.  The
	//    program searches within the circle to find events larger than the mainshock that
	//    might be shadowing the mainshock.  A recommended value is 2000 km, which is likely
	//    large enough to find the largest possible earthquake (assuming that the centroid
	//    radius multiplier in the magnitude-of-completeness parameters is equal to 1.0).
	//    The value must be positive.
	//  search_time_lo, search_time_hi = Time interval, expressed in milliseconds since
	//    the epoch.  The program searches for possible shadowing events that lie between
	//    the two times.  Typically these are chosen to bracket the mainshock origin time,
	//    and indicate how far apart in time events must be so they don't shadow each other.
	//    Typical values are 1 year before and after the mainshock origin time.  It is
	//    permitted for search_time_hi to be larger than time_now;  the effect is as if
	//    search_time_hi were set equal to time_now.  The value of search_time_lo must be
	//    strictly less than search_time_hi, time_now, and the current time as reported by
	//    System.currentTimeMillis().
	//  centroid_rel_time_lo, centroid_rel_time_hi = Relative time interval, expressed in
	//    milliseconds.  For each candidate shadowing event, the program determines a time
	//    interval by adding these values to the event origin time;  aftershocks occurring
	//    during that time interval are used to compute the centroid (except that
	//    aftershocks occuring after time_now are not considered).  Typical values are 0
	//    and 1 year.  It is required that centroid_rel_time_lo be non-negative, and that
	//    centroid_rel_time_hi be greater than centroid_rel_time_lo.
	//  centroid_mag_floor = Minimum magnitude for aftershocks that are used to compute
	//    the centroid.  The program also considers that centroid magnitude in the
	//    magnitude-of-completeness parameters, and uses the larger of that centroid
	//    magnitude and centroid_mag_floor.  Typical values are 2.5 (suitable for California)
	//    to 3.5 (suitable world-wide).  Note that it is not necessary to vary this with
	//    the mainshock location, if the appropriate variation is set in the magnitude-of-
	//    completeness parameters.  Be cautious about using values smaller than 3.0, as
	//    the amount of data to process increases very rapidly with reduced magnitude.
	//  large_mag = Minimum magnitude for a candidate shadowing event to be considered large.
	//    To perform the centroid algorithm, a separate call to Comcat is made for each large
	//    candidate.  All small candidates are grouped together, and the possible aftershocks
	//    for all of them are retrieved in a single call to Comcat.  The latter call to
	//    to Comcat may cover a region whose radius is as large as three times the Wells and
	//    Coppersmith radius of large_mag (assuming the centroid radius multiplier in the
	//    magnitude-of-completeness parameters is 1.0).  A typical value is 8.0, which gives
	//    a W&C radius of 200 km, and so limits the combined call to a radius of 600 km.
	//    This mechanism avoids the possiblity that the combined call might attempt to
	//    retrieve all small earthquakes over a very large area.  Set to 10.0 to disable.
	//  separation = A 2-element array that is used to return the separation between
	//    the mainshock and the shadowing event.  If the mainshock is shadowed, then
	//    separation[0] receives the separation in kilometers, and separation[1] receives
	//    the separation in days (positive means the mainshock occurs after the shadowing
	//    event).  Can be null if separation is not required.
	// Returns:
	// If the mainshock is not shadowed, then the return value is null.
	// If the mainshock is shadowed, then the return value is the shadowing earthquake.
	//   If there are multiple shadowing earthquakes, then the one with largest magnitude
	//   is returned.  If magnitudes are tied, then the earliest earthquake is returned.
	// An exception from this function likely means a Comcat failure.

	public static ObsEqkRupture find_shadow (ObsEqkRupture mainshock, long time_now,
					double search_radius, long search_time_lo, long search_time_hi,
					long centroid_rel_time_lo, long centroid_rel_time_hi,
					double centroid_mag_floor, double large_mag, double[] separation) {

		// Parameter validation

		if (!( mainshock != null )) {
			throw new IllegalArgumentException ("AftershockStatsShadow.find_shadow: No mainshock supplied");
		}

		long system_time_now = System.currentTimeMillis();

		if (!( search_time_lo < system_time_now
			&& search_time_lo < search_time_hi
			&& search_time_lo < time_now )) {
			throw new IllegalArgumentException ("AftershockStatsShadow.find_shadow: Invalid search times"
				+ ": search_time_lo = " + search_time_lo
				+ ", search_time_hi = " + search_time_hi
				+ ", time_now = " + time_now
				+ ", system_time_now = " + system_time_now
				);
		}

		if (!( centroid_rel_time_lo >= 0L
			&& centroid_rel_time_lo < centroid_rel_time_hi )) {
			throw new IllegalArgumentException ("AftershockStatsShadow.find_shadow: Invalid centroid relative times"
				+ ": centroid_rel_time_lo = " + centroid_rel_time_lo
				+ ", centroid_rel_time_hi = " + centroid_rel_time_hi
			);
		}

		if (!( search_radius > 0.0 )) {
			throw new IllegalArgumentException ("AftershockStatsShadow.find_shadow: Invalid search radius"
				+ ": search_radius = " + search_radius
			);
		}

		if (separation != null) {
			if (separation.length < 2) {
				throw new IllegalArgumentException ("AftershockStatsShadow.find_shadow: Separation array is too short");
			}
		}

		// List of candidate shadowing earthquakes

		ArrayList<CandidateShadow> candidates = new ArrayList<CandidateShadow>();

		// List of candidate shadowing earthquakes that are combined into a single aftershock call to Comcat

		ArrayList<CandidateShadow> combined_candidates = new ArrayList<CandidateShadow>();

		// A Comcat accessor to use

		ComcatAccessor accessor = new ComcatAccessor();

		// Fetch object for magnitude of completeness parameters

		MagCompPage_ParametersFetch mag_comp_fetch = new MagCompPage_ParametersFetch();

		// Wells and Coppersmith relation

		WC1994_MagLengthRelationship wcMagLen = new WC1994_MagLengthRelationship();

		// Get the mainshock parameters

		String mainshock_event_id = mainshock.getEventId();
		long mainshock_time = mainshock.getOriginTime();
		double mainshock_mag = mainshock.getMag();
		Location mainshock_hypo = mainshock.getHypocenterLocation();

		// Construct a circle around the mainshock with the search radius

		SphLatLon mainshock_sph_hypo = new SphLatLon (mainshock_hypo);
		SphRegion search_region = SphRegion.makeCircle (mainshock_sph_hypo, search_radius);

		// Get a list of potential candidates by calling Comcat
		// Potentials must lie in the search region, within the search times,
		// have magnitude at least equal to the mainshock, and not be the mainshock

		double min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;
		double max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;

		boolean wrapLon = false;
		int limit_per_call = 0;
		int max_calls = 0;

		ObsEqkRupList potentials = accessor.fetchEventList (mainshock_event_id,
					search_time_lo, Math.min (search_time_hi, time_now),
					min_depth, max_depth,
					search_region, wrapLon,
					mainshock_mag, limit_per_call, max_calls);

		System.out.println ("AftershockStatsShadow.find_shadow: Found " + potentials.size() + " potential shadowing events for mainshock " + mainshock_event_id);

		// Combined values for centroid radius, minimum magnitude, start time, and stop time

		double combined_centroid_radius = 0.0;
		double combined_centroid_min_mag = 10.0;
		long combined_centroid_time_lo = Long.MAX_VALUE;
		long combined_centroid_time_hi = Long.MIN_VALUE;

		// Examine the potentials and accumulate the candidates

		for (ObsEqkRupture potential : potentials) {

			// Get the potential parameters

			String potential_event_id = potential.getEventId();
			long potential_time = potential.getOriginTime();
			double potential_mag = potential.getMag();
			Location potential_hypo = potential.getHypocenterLocation();

			// The accessor should have checked that event ID is different and the potential
			// magnitude and time are in range, but repeat the check anyway;
			// magnitude ties are broken by considering the earlier origin time

			if ( (!(potential_event_id.equals (mainshock_event_id)))
				&& potential_time <= Math.min (search_time_hi, time_now)
				&& potential_time >= search_time_lo
				&& ( potential_mag > mainshock_mag
					|| (potential_mag == mainshock_mag && potential_time < mainshock_time) ) ) {

				// Get the magnitude of completeness parameters

				MagCompPage_Parameters mag_comp_params = mag_comp_fetch.get (potential_hypo);

				// Get the Wells and Coppersmith radius

				double wc_radius = wcMagLen.getMedianLength (potential_mag);

				// Get the centroid radius

				double centroid_radius = wc_radius * mag_comp_params.get_radiusCentroid();

				// Get the sample radius

				double sample_radius = wc_radius * mag_comp_params.get_radiusSample();

				// Get the distance to the mainshock

				double dist = LocationUtils.horzDistance (mainshock_hypo, potential_hypo);

				// If the potential is close enough to the mainshock so it could possibly shadow it ...

				if (sample_radius + centroid_radius > dist) {

					// Get the centroid minimum magnitude

					double centroid_min_mag = Math.max (mag_comp_params.get_magCentroid(), centroid_mag_floor);

					// Get the centroid start and end times

					long centroid_time_lo = potential_time + centroid_rel_time_lo;
					long centroid_time_hi = Math.min (potential_time + centroid_rel_time_hi, time_now);

					// Create the candidate

					CandidateShadow candidate = new CandidateShadow (
													potential,
													centroid_radius,
													centroid_min_mag,
													centroid_time_lo,
													centroid_time_hi,
													sample_radius);

					// Add to the list of candidates

					candidates.add (candidate);

					// If the candidate is large ...

					if (potential_mag >= large_mag) {

						// Accumulate its aftershocks now, in a separate call to Comcat

						candidate.accum_from_comcat (accessor, system_time_now);
					}

					// Otherwise the candidate is small ...

					else {

						// If the candidate can accept aftershocks ...

						if (centroid_time_lo < centroid_time_hi) {

							// Adjust the combined values to include what this candidate needs

							combined_centroid_radius = Math.max (combined_centroid_radius, dist + centroid_radius);
							combined_centroid_min_mag = Math.min (combined_centroid_min_mag, centroid_min_mag);
							combined_centroid_time_lo = Math.min (combined_centroid_time_lo, centroid_time_lo);
							combined_centroid_time_hi = Math.max (combined_centroid_time_hi, centroid_time_hi);

							// Place on the list of combined candidates

							combined_candidates.add (candidate);
						}
					}
				}
			}
		}

		System.out.println ("AftershockStatsShadow.find_shadow: Found " + candidates.size() + " candidate shadowing events for mainshock " + mainshock_event_id);

		// If there is a nonempty time and space region in which we need to look for aftershocks ...

		if (combined_centroid_time_lo < system_time_now
			&& combined_centroid_time_lo < combined_centroid_time_hi
			&& combined_centroid_radius > 0.0) {

			// Construct a circle around the mainshock with the combined centroid radius

			SphRegion centroid_region = SphRegion.makeCircle (mainshock_sph_hypo, combined_centroid_radius);

			// Get a list of possible aftershocks by calling Comcat
			// Aftershocks must lie in the combined centroid region, within the combined centroid times,
			// and have magnitude at least equal to the combined centroid magnitude

			ObsEqkRupList aftershocks = accessor.fetchEventList (null,
						combined_centroid_time_lo, combined_centroid_time_hi,
						min_depth, max_depth,
						centroid_region, wrapLon,
						combined_centroid_min_mag, limit_per_call, max_calls);

			System.out.println ("AftershockStatsShadow.find_shadow: Found " + aftershocks.size() + " possible aftershocks within " + String.format ("%.3f", combined_centroid_radius) + " km of mainshock " + mainshock_event_id);
		
			// For each aftershock ...

			for (ObsEqkRupture aftershock : aftershocks) {

				// Get the aftershock parameters

				String aftershock_event_id = aftershock.getEventId();
				long aftershock_time = aftershock.getOriginTime();
				double aftershock_mag = aftershock.getMag();
				Location aftershock_hypo = aftershock.getHypocenterLocation();

				// Get unit vector for centroid calculation (see AftershockStatsCalc.getSphCentroid)

				double lat = aftershock_hypo.getLatRad();
				double lon = aftershock_hypo.getLonRad();

				double x = (Math.cos(lat) * Math.cos(lon));
				double y = (Math.cos(lat) * Math.sin(lon));
				double z = Math.sin(lat);

				// For each candidate ...

				for (CandidateShadow candidate : combined_candidates) {

					// This aftershock contributes to the candidate centroid if it lies in the
					// time interval, within the radius, has sufficient magnitude, and is not the candidate

					if (   aftershock_time >= candidate.centroid_time_lo
						&& aftershock_time <= candidate.centroid_time_hi
						&& aftershock_mag >= candidate.centroid_min_mag
						&& LocationUtils.horzDistance (aftershock_hypo, candidate.candidate_hypo) <= candidate.centroid_radius
						&& (!( aftershock_event_id.equals (candidate.candidate_event_id) )) ) {

						// Add to the centroid accumulators

						candidate.x += x;
						candidate.y += y;
						candidate.z += z;
					}
				}
			}
		}

		// This is the best candidate so far

		CandidateShadow best_candidate = null;

		// For each candidate ...

		for (CandidateShadow candidate : candidates) {

			// Get the centroid

			Location centroid = candidate.get_centroid();

			// If the centroid is within the sample radius of the mainshock, then it's shadowing the mainshock

			if (LocationUtils.horzDistance (mainshock_hypo, centroid) <= candidate.sample_radius) {
			
				// If this is a new best candidate, save it

				if (best_candidate == null
					|| candidate.candidate_mag > best_candidate.candidate_mag
					|| (candidate.candidate_mag == best_candidate.candidate_mag && candidate.candidate_time < best_candidate.candidate_time) ) {
					
					best_candidate = candidate;
				}
			}
		}

		// If no candidates shadow the mainshock, then return null

		if (best_candidate == null) {
			System.out.println ("AftershockStatsShadow.find_shadow: Mainshock " + mainshock_event_id + " is not shadowed");
			return null;
		}

		// Return the largest-magnitude event that shadows the mainshock

		double best_distance = LocationUtils.horzDistance (mainshock_hypo, best_candidate.candidate_hypo);
		double best_time_offset = ((double)(mainshock_time - best_candidate.candidate_time))/ComcatAccessor.day_millis;

		if (separation != null) {
			separation[0] = best_distance;
			separation[1] = best_time_offset;
		}

		System.out.println ("AftershockStatsShadow.find_shadow: Mainshock " + mainshock_event_id + " is shadowed by event " + best_candidate.candidate_event_id);
		System.out.println (String.format ("AftershockStatsShadow.find_shadow: Mainshock magnitude = %.2f, shadowing event magnitude = %.2f", mainshock_mag, best_candidate.candidate_mag));
		System.out.println (String.format ("AftershockStatsShadow.find_shadow: Distance = %.3f km, time offset = %.3f days", best_distance, best_time_offset));

		return best_candidate.rupture;
	}




	// Constructor.

	public AftershockStatsShadow () {
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("AftershockStatsShadow : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  event_id  days
		// Fetch information for an event, and display it.
		// Then, test if the event is shadowed.  The time_now is set to the specified
		//  number of days after the event, and other parameters are set to defaults.
		// Then, display the shadowing event if one is found.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Two additional arguments

			if (args.length != 3) {
				System.err.println ("AftershockStatsShadow : Invalid 'test1' subcommand");
				return;
			}

			try {

				String event_id = args[1];
				double days = Double.parseDouble (args[2]);

				// Say hello

				System.out.println ("Fetching event: " + event_id);

				// Create the accessor

				ComcatAccessor accessor = new ComcatAccessor();

				// Get the rupture

				ObsEqkRupture rup = accessor.fetchEvent (event_id, false, true);

				// Display its information

				if (rup == null) {
					System.out.println ("Null return from fetchEvent");
					System.out.println ("http_status = " + accessor.get_http_status_code());
					return;
				}

				String rup_event_id = rup.getEventId();
				long rup_time = rup.getOriginTime();
				double rup_mag = rup.getMag();
				Location rup_hypo = rup.getHypocenterLocation();
				double rup_lat = rup_hypo.getLatitude();
				double rup_lon = rup_hypo.getLongitude();
				double rup_depth = rup_hypo.getDepth();

				System.out.println ("rup_event_id = " + rup_event_id);
				System.out.println ("rup_time = " + rup_time + " (" + SimpleUtils.time_to_string(rup_time) + ")");
				System.out.println ("rup_mag = " + rup_mag);
				System.out.println ("rup_lat = " + rup_lat);
				System.out.println ("rup_lon = " + rup_lon);
				System.out.println ("rup_depth = " + rup_depth);

				// Get find shadow parameters

				long time_now = rup_time + (long)(days*ComcatAccessor.day_millis);
				double search_radius = DEF_SEARCH_RADIUS;
				long search_time_lo = rup_time - YEAR_IN_MILLIS;
				long search_time_hi = rup_time + YEAR_IN_MILLIS;
				long centroid_rel_time_lo = 0L;
				long centroid_rel_time_hi = YEAR_IN_MILLIS;
				double centroid_mag_floor = DEF_CENTROID_MAG_FLOOR;
				double large_mag = DEF_LARGE_MAG;
				double[] separation = new double[2];

				System.out.println ("");
				System.out.println ("find_shadow parameters:");
				System.out.println ("time_now = " + time_now + " (" + SimpleUtils.time_to_string(time_now) + ")");
				System.out.println ("search_radius = " + search_radius);
				System.out.println ("search_time_lo = " + search_time_lo + " (" + SimpleUtils.time_to_string(search_time_lo) + ")");
				System.out.println ("search_time_hi = " + search_time_hi + " (" + SimpleUtils.time_to_string(search_time_hi) + ")");
				System.out.println ("centroid_rel_time_lo = " + centroid_rel_time_lo);
				System.out.println ("centroid_rel_time_hi = " + centroid_rel_time_hi);
				System.out.println ("centroid_mag_floor = " + centroid_mag_floor);
				System.out.println ("large_mag = " + large_mag);

				// Run find_shadow

				System.out.println ("");
				System.out.println ("Finding shadow:");

				ObsEqkRupture shadow = find_shadow (rup, time_now,
					search_radius, search_time_lo, search_time_hi,
					centroid_rel_time_lo, centroid_rel_time_hi,
					centroid_mag_floor, large_mag, separation);

				// Display results

				System.out.println ("");

				if (shadow == null) {
					System.out.println ("Event is not shadowed");
				} else {
					System.out.println ("Event is shadowed by:");

					String shadow_event_id = shadow.getEventId();
					long shadow_time = shadow.getOriginTime();
					double shadow_mag = shadow.getMag();
					Location shadow_hypo = shadow.getHypocenterLocation();
					double shadow_lat = shadow_hypo.getLatitude();
					double shadow_lon = shadow_hypo.getLongitude();
					double shadow_depth = shadow_hypo.getDepth();

					System.out.println ("shadow_event_id = " + shadow_event_id);
					System.out.println ("shadow_time = " + shadow_time + " (" + SimpleUtils.time_to_string(shadow_time) + ")");
					System.out.println ("shadow_mag = " + shadow_mag);
					System.out.println ("shadow_lat = " + shadow_lat);
					System.out.println ("shadow_lon = " + shadow_lon);
					System.out.println ("shadow_depth = " + shadow_depth);

					System.out.println ("separation_km = " + separation[0]);
					System.out.println ("separation_days = " + separation[1]);
				}

            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("AftershockStatsShadow : Unrecognized subcommand : " + args[0]);
		return;

	}

}
