package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

import java.time.Duration;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;

import scratch.aftershockStatistics.OAFParameterSet;

/**
 * Configuration file for AAFS server actions.
 * Author: Michael Barall 04/29/2018.
 *
 * All fields are public, since this is just a buffer for reading and writing files.
 *
 * JSON file format:
 *
 *	"ActionConfigFile" = Integer giving file version number, currently 24001.
 *	"forecast_min_gap" = String giving minimum allowed gap between forcasts, in java.time.Duration format.
 *	"forecast_max_delay" = String giving maximum allowed delay in reporting a forecast to PDL, in java.time.Duration format.
 *	"comcat_clock_skew" = Assumed maximum difference between our clock and ComCat clock, in java.time.Duration format.
 *	"comcat_origin_skew" = Assumed maximum change in mainshock origin time, in java.time.Duration format.
 *	"comcat_retry_min_gap" = String giving minimum allowed gap between ComCat retries, in java.time.Duration format.
 *  "seq_spec_min_lag" = String giving minimum time lag at which sequence-specific forecasts can be generated, in java.time.Duration format.
 *  "advisory_dur_week" = String giving minimum time lag at which one-week advisories can be generated, in java.time.Duration format.
 *  "advisory_dur_month" = String giving minimum time lag at which one-month advisories can be generated, in java.time.Duration format.
 *  "advisory_dur_year" = String giving minimum time lag at which one-year advisories can be generated, in java.time.Duration format.
 *	"forecast_lags" = [ Array giving a list of time lags at which forecasts are generated, in increasing order.
 *		element = String giving time lag since mainshock, in java.time.Duration format.
 *	]
 *	"comcat_retry_lags" = [ Array giving a list of time lags at which ComCat forecast operations are retried, in increasing order.
 *		element = String giving time lag since first attempt, in java.time.Duration format.
 *	]
 *	"comcat_intake_lags" = [ Array giving a list of time lags at which ComCat intake operations are retried, in increasing order.
 *		element = String giving time lag since first attempt, in java.time.Duration format.
 *	]
 *	"pdl_report_retry_lags" = [ Array giving a list of time lags at which PDL report attempts are retried, in increasing order.
 *		element = String giving time lag since first attempt, in java.time.Duration format.
 *	]
 *	"pdl_intake_regions" = [ Array giving a list of intake regions, in the order they are checked.
 *		element = IntakeSphRegion JSON representation, giving region and magnitude thresholds.
 *	]
 */
public class ActionConfigFile {

	//----- Parameter values -----

	// Minimum gap between forecasts, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long forecast_min_gap;

	// Maximum delay in reporting a forecast to PDL, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long forecast_max_delay;

	// Assumed maximum difference between our clock and ComCat's clock, in milliseconds.
	// (Specifically, if an earthquake occurs at time T then it should be visible in
	// ComCat by the time our clock reads T + comcat_clock_skew.)
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long comcat_clock_skew;

	// Assumed maximum change in mainshock origin time, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long comcat_origin_skew;

	// Minimum gap between ComCat retries, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long comcat_retry_min_gap;

	// Minimum time after an earthquake at which sequence-specific forecasts can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long seq_spec_min_lag;

	// Minimum time after an earthquake at which one-week advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long advisory_dur_week;

	// Minimum time after an earthquake at which one-month advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long advisory_dur_month;

	// Minimum time after an earthquake at which one-year advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long advisory_dur_year;

	// Time lags at which forecasts are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the mainshock.  Must have at least 1 element.
	// The difference between successive elements must be at least the minimum gap.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> forecast_lags;

	// Time lags at which ComCat forecast retries are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the initial attempt.
	// The difference between successive elements must be at least the minimum gap.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> comcat_retry_lags;

	// Time lags at which ComCat intake retries are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the initial attempt.
	// This is used when looking for initial appearance in ComCat of an event reported by PDL.
	// The difference between successive elements must be at least the minimum gap.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> comcat_intake_lags;

	// Time lags at which PDL report retries are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the initial attempt.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> pdl_report_retry_lags;

	// Regions for filtering events received from PDL.

	public ArrayList<IntakeSphRegion> pdl_intake_regions;


	//----- Construction -----

	// Default constructor.

	public ActionConfigFile () {
		clear();
	}

	// Clear the contents.

	public void clear () {
		forecast_min_gap = 0L;
		forecast_max_delay = 0L;
		comcat_clock_skew = 0L;
		comcat_origin_skew = 0L;
		comcat_retry_min_gap = 0L;
		seq_spec_min_lag = 0L;
		advisory_dur_week = 0L;
		advisory_dur_month = 0L;
		advisory_dur_year = 0L;
		forecast_lags = new ArrayList<Long>();
		comcat_retry_lags = new ArrayList<Long>();
		comcat_intake_lags = new ArrayList<Long>();
		pdl_report_retry_lags = new ArrayList<Long>();
		pdl_intake_regions = new ArrayList<IntakeSphRegion>();
		return;
	}

	// Check that a lag time is valid: a whole number of seconds, between 1 and 10^9 seconds.

	private static final long MIN_LAG = 1000L;
	private static final long MAX_LAG = 1000000000000L;
	private static final long UNIT_LAG = 1000L;

	private boolean is_valid_lag (long lag) {
		return (lag >= MIN_LAG && lag <= MAX_LAG && lag % UNIT_LAG == 0L);
	}

	private boolean is_valid_lag (long lag, long min_lag) {
		return (lag >= min_lag && lag <= MAX_LAG && lag % UNIT_LAG == 0L);
	}

	// Check that values are valid, throw an exception if not.

	public void check_invariant () {

		if (!( is_valid_lag(forecast_min_gap) )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_min_gap: " + forecast_min_gap);
		}

		if (!( is_valid_lag(forecast_max_delay) )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_max_delay: " + forecast_max_delay);
		}

		if (!( is_valid_lag(comcat_clock_skew) )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_clock_skew: " + comcat_clock_skew);
		}

		if (!( is_valid_lag(comcat_origin_skew) )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_origin_skew: " + comcat_origin_skew);
		}

		if (!( is_valid_lag(comcat_retry_min_gap) )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_retry_min_gap: " + comcat_retry_min_gap);
		}

		if (!( is_valid_lag(seq_spec_min_lag) )) {
				throw new RuntimeException("ActionConfigFile: Invalid seq_spec_min_lag: " + seq_spec_min_lag);
		}

		if (!( is_valid_lag(advisory_dur_week) )) {
				throw new RuntimeException("ActionConfigFile: Invalid advisory_dur_week: " + advisory_dur_week);
		}

		if (!( is_valid_lag(advisory_dur_month) )) {
				throw new RuntimeException("ActionConfigFile: Invalid advisory_dur_month: " + advisory_dur_month);
		}

		if (!( is_valid_lag(advisory_dur_year) )) {
				throw new RuntimeException("ActionConfigFile: Invalid advisory_dur_year: " + advisory_dur_year);
		}

		int n = forecast_lags.size();
		long min_lag = MIN_LAG;
		
		if (!( n > 0 )) {
				throw new RuntimeException("ActionConfigFile: Empty forecast_lags list");
		}

		for (int i = 0; i < n; ++i) {
			long forecast_lag = forecast_lags.get(i);
			if (!( is_valid_lag(forecast_lag, min_lag) )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_lag: " + forecast_lag + ", index = " + i);
			}
			min_lag = forecast_lag + forecast_min_gap;
		}

		n = comcat_retry_lags.size();
		min_lag = MIN_LAG;

		for (int i = 0; i < n; ++i) {
			long comcat_retry_lag = comcat_retry_lags.get(i);
			if (!( is_valid_lag(comcat_retry_lag, min_lag) )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_retry_lag: " + comcat_retry_lag + ", index = " + i);
			}
			min_lag = comcat_retry_lag + comcat_retry_min_gap;
		}

		n = comcat_intake_lags.size();
		min_lag = MIN_LAG;

		for (int i = 0; i < n; ++i) {
			long comcat_intake_lag = comcat_intake_lags.get(i);
			if (!( is_valid_lag(comcat_intake_lag, min_lag) )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_intake_lag: " + comcat_intake_lag + ", index = " + i);
			}
			min_lag = comcat_intake_lag + comcat_retry_min_gap;
		}

		n = pdl_report_retry_lags.size();
		min_lag = MIN_LAG;

		for (int i = 0; i < n; ++i) {
			long pdl_report_retry_lag = pdl_report_retry_lags.get(i);
			if (!( is_valid_lag(pdl_report_retry_lag, min_lag) )) {
				throw new RuntimeException("ActionConfigFile: Invalid pdl_report_retry_lag: " + pdl_report_retry_lag + ", index = " + i);
			}
			min_lag = pdl_report_retry_lag + comcat_retry_min_gap;
		}

		return;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		result.append ("ActionConfigFile:" + "\n");

		result.append ("forecast_min_gap = " + Duration.ofMillis(forecast_min_gap).toString() + "\n");
		result.append ("forecast_max_delay = " + Duration.ofMillis(forecast_max_delay).toString() + "\n");
		result.append ("comcat_clock_skew = " + Duration.ofMillis(comcat_clock_skew).toString() + "\n");
		result.append ("comcat_origin_skew = " + Duration.ofMillis(comcat_origin_skew).toString() + "\n");
		result.append ("comcat_retry_min_gap = " + Duration.ofMillis(comcat_retry_min_gap).toString() + "\n");
		result.append ("seq_spec_min_lag = " + Duration.ofMillis(seq_spec_min_lag).toString() + "\n");
		result.append ("advisory_dur_week = " + Duration.ofMillis(advisory_dur_week).toString() + "\n");
		result.append ("advisory_dur_month = " + Duration.ofMillis(advisory_dur_month).toString() + "\n");
		result.append ("advisory_dur_year = " + Duration.ofMillis(advisory_dur_year).toString() + "\n");

		result.append ("forecast_lags = [" + "\n");
		for (int i = 0; i < forecast_lags.size(); ++i) {
			long forecast_lag = forecast_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(forecast_lag).toString() + "\n");
		}
		result.append ("]" + "\n");

		result.append ("comcat_retry_lags = [" + "\n");
		for (int i = 0; i < comcat_retry_lags.size(); ++i) {
			long comcat_retry_lag = comcat_retry_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(comcat_retry_lag).toString() + "\n");
		}
		result.append ("]" + "\n");

		result.append ("comcat_intake_lags = [" + "\n");
		for (int i = 0; i < comcat_intake_lags.size(); ++i) {
			long comcat_intake_lag = comcat_intake_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(comcat_intake_lag).toString() + "\n");
		}
		result.append ("]" + "\n");

		result.append ("pdl_report_retry_lags = [" + "\n");
		for (int i = 0; i < pdl_report_retry_lags.size(); ++i) {
			long pdl_report_retry_lag = pdl_report_retry_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(pdl_report_retry_lag).toString() + "\n");
		}
		result.append ("]" + "\n");

		result.append ("pdl_intake_regions = [" + "\n");
		for (int i = 0; i < pdl_intake_regions.size(); ++i) {
			result.append (pdl_intake_regions.get(i).toString() + "\n");
		}
		result.append ("]" + "\n");

		return result.toString();
	}


	//----- Service functions -----

	// Get the first element of forecast_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_forecast_lag (long min_lag) {

		// Binary search

		int index = Collections.binarySearch (forecast_lags, new Long(min_lag));

		// If not found, convert to index of next larger element

		if (index < 0) {
			index = -(index + 1);
		}

		// If past end of list, then return -1

		if (index >= forecast_lags.size()) {
			return -1L;
		}

		// Return the lag value from the list

		return forecast_lags.get(index).longValue();
	}

	// Get the first element of comcat_retry_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_comcat_retry_lag (long min_lag) {

		// Binary search

		int index = Collections.binarySearch (comcat_retry_lags, new Long(min_lag));

		// If not found, convert to index of next larger element

		if (index < 0) {
			index = -(index + 1);
		}

		// If past end of list, then return -1

		if (index >= comcat_retry_lags.size()) {
			return -1L;
		}

		// Return the lag value from the list

		return comcat_retry_lags.get(index).longValue();
	}

	// Get the first element of comcat_intake_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_comcat_intake_lag (long min_lag) {

		// Binary search

		int index = Collections.binarySearch (comcat_intake_lags, new Long(min_lag));

		// If not found, convert to index of next larger element

		if (index < 0) {
			index = -(index + 1);
		}

		// If past end of list, then return -1

		if (index >= comcat_intake_lags.size()) {
			return -1L;
		}

		// Return the lag value from the list

		return comcat_intake_lags.get(index).longValue();
	}

	// Get the first element of pdl_report_retry_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_pdl_report_retry_lag (long min_lag) {

		// Binary search

		int index = Collections.binarySearch (pdl_report_retry_lags, new Long(min_lag));

		// If not found, convert to index of next larger element

		if (index < 0) {
			index = -(index + 1);
		}

		// If past end of list, then return -1

		if (index >= pdl_report_retry_lags.size()) {
			return -1L;
		}

		// Return the lag value from the list

		return pdl_report_retry_lags.get(index).longValue();
	}

	// Get the pdl intake region that satisfies the min_mag criterion.
	// If found, the region is returned.
	// If not found, null is returned.

	public IntakeSphRegion get_pdl_intake_region_for_min_mag (double lat, double lon, double mag) {

		// Construct the SphLatLon object

		double the_lon = lon;
		while (the_lon < -180.0) {
			the_lon += 360.0; 
		}
		while (the_lon > 180.0) {
			the_lon -= 360.0; 
		}
		SphLatLon loc = new SphLatLon (lat, the_lon);

		// Search the list of regions

		for (IntakeSphRegion intake_region : pdl_intake_regions) {
			if (intake_region.contains (loc, mag)) {
				return intake_region;
			}
		}
		return null;
	} 

	// Get the pdl intake region that satisfies the intake_mag criterion.
	// If found, the region is returned.
	// If not found, null is returned.

	public IntakeSphRegion get_pdl_intake_region_for_intake_mag (double lat, double lon, double mag) {

		// Construct the SphLatLon object

		double the_lon = lon;
		while (the_lon < -180.0) {
			the_lon += 360.0; 
		}
		while (the_lon > 180.0) {
			the_lon -= 360.0; 
		}
		SphLatLon loc = new SphLatLon (lat, the_lon);

		// Search the list of regions

		for (IntakeSphRegion intake_region : pdl_intake_regions) {
			if (intake_region.contains_intake (loc, mag)) {
				return intake_region;
			}
		}
		return null;
	} 




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 24001;

	private static final String M_VERSION_NAME = "ActionConfigFile";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 24000;
	protected static final int MARSHAL_ACTION_CFG = 24001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_ACTION_CFG;
	}

	// Marshal a duration.

	public static void marshal_duration (MarshalWriter writer, String name, long duration) {
		writer.marshalString (name, Duration.ofMillis(duration).toString());
		return;
	}

	// Unmarshal a duration.

	public static long unmarshal_duration (MarshalReader reader, String name) {
		return Duration.parse(reader.unmarshalString(name)).toMillis();
	}

	// Marshal a duration list.

	public static void marshal_duration_list (MarshalWriter writer, String name, List<Long> durations) {
		int n = durations.size();
		writer.marshalArrayBegin (name, n);
		for (Long duration : durations) {
			marshal_duration (writer, null, duration.longValue());
		}
		writer.marshalArrayEnd ();
		return;
	}

	// Unmarshal a duration list.

	public static ArrayList<Long> unmarshal_duration_list (MarshalReader reader, String name) {
		ArrayList<Long> duration_list = new ArrayList<Long>();
		int n = reader.unmarshalArrayBegin (name);
		for (int i = 0; i < n; ++i) {
			duration_list.add (new Long (unmarshal_duration (reader, null)));
		}
		reader.unmarshalArrayEnd ();
		return duration_list;
	}

	// Marshal an intake region.

	public static void marshal_intake_region (MarshalWriter writer, String name, IntakeSphRegion intake_region) {

		// Begin the JSON object

		writer.marshalMapBegin (name);

		// Write the name

		writer.marshalString ("name", intake_region.get_name());

		// Write the magnitudes
				
		writer.marshalDouble ("min_mag", intake_region.get_min_mag());
		writer.marshalDouble ("intake_mag", intake_region.get_intake_mag());

		// Write the spherical region

		SphRegion.marshal_poly (writer, "region", intake_region.get_region()) ;

		// End the JSON object

		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal an intake region.

	public static IntakeSphRegion unmarshal_intake_region (MarshalReader reader, String name) {

		// Begin the JSON object

		reader.unmarshalMapBegin (name);

		// Get the name

		String region_name = reader.unmarshalString ("name");

		// Get the magnitudes
				
		double min_mag = reader.unmarshalDouble ("min_mag");
		double intake_mag = reader.unmarshalDouble ("intake_mag");

		// Get the spherical region

		SphRegion sph_region = SphRegion.unmarshal_poly (reader, "region") ;

		if (sph_region == null) {
			throw new RuntimeException("ActionConfigFile: No spherical region specified");
		}

		// End the JSON object

		reader.unmarshalMapEnd ();

		// Form the region

		return new IntakeSphRegion (region_name, sph_region, min_mag, intake_mag);
	}

	// Marshal an intake region list.

	public static void marshal_intake_region_list (MarshalWriter writer, String name, List<IntakeSphRegion> intake_regions) {
		int n = intake_regions.size();
		writer.marshalArrayBegin (name, n);
		for (IntakeSphRegion intake_region : intake_regions) {
			marshal_intake_region (writer, null, intake_region);
		}
		writer.marshalArrayEnd ();
		return;
	}

	// Unmarshal an intake region list.

	public static ArrayList<IntakeSphRegion> unmarshal_intake_region_list (MarshalReader reader, String name) {
		ArrayList<IntakeSphRegion> intake_region_list = new ArrayList<IntakeSphRegion>();
		int n = reader.unmarshalArrayBegin (name);
		for (int i = 0; i < n; ++i) {
			intake_region_list.add (unmarshal_intake_region (reader, null));
		}
		reader.unmarshalArrayEnd ();
		return intake_region_list;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Error check

		check_invariant();

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		marshal_duration           (writer, "forecast_min_gap"     , forecast_min_gap     );
		marshal_duration           (writer, "forecast_max_delay"   , forecast_max_delay   );
		marshal_duration           (writer, "comcat_clock_skew"    , comcat_clock_skew    );
		marshal_duration           (writer, "comcat_origin_skew"   , comcat_origin_skew   );
		marshal_duration           (writer, "comcat_retry_min_gap" , comcat_retry_min_gap );
		marshal_duration           (writer, "seq_spec_min_lag"     , seq_spec_min_lag     );
		marshal_duration           (writer, "advisory_dur_week"    , advisory_dur_week    );
		marshal_duration           (writer, "advisory_dur_month"   , advisory_dur_month   );
		marshal_duration           (writer, "advisory_dur_year"    , advisory_dur_year    );
		marshal_duration_list      (writer, "forecast_lags"        , forecast_lags        );
		marshal_duration_list      (writer, "comcat_retry_lags"    , comcat_retry_lags    );
		marshal_duration_list      (writer, "comcat_intake_lags"   , comcat_intake_lags   );
		marshal_duration_list      (writer, "pdl_report_retry_lags", pdl_report_retry_lags);
		marshal_intake_region_list (writer, "pdl_intake_regions"   , pdl_intake_regions   );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		forecast_min_gap      = unmarshal_duration           (reader, "forecast_min_gap"     );
		forecast_max_delay    = unmarshal_duration           (reader, "forecast_max_delay"   );
		comcat_clock_skew     = unmarshal_duration           (reader, "comcat_clock_skew"    );
		comcat_origin_skew    = unmarshal_duration           (reader, "comcat_origin_skew"   );
		comcat_retry_min_gap  = unmarshal_duration           (reader, "comcat_retry_min_gap" );
		seq_spec_min_lag      = unmarshal_duration           (reader, "seq_spec_min_lag"     );
		advisory_dur_week     = unmarshal_duration           (reader, "advisory_dur_week"    );
		advisory_dur_month    = unmarshal_duration           (reader, "advisory_dur_month"   );
		advisory_dur_year     = unmarshal_duration           (reader, "advisory_dur_year"    );
		forecast_lags         = unmarshal_duration_list      (reader, "forecast_lags"        );
		comcat_retry_lags     = unmarshal_duration_list      (reader, "comcat_retry_lags"    );
		comcat_intake_lags    = unmarshal_duration_list      (reader, "comcat_intake_lags"   );
		pdl_report_retry_lags = unmarshal_duration_list      (reader, "pdl_report_retry_lags");
		pdl_intake_regions    = unmarshal_intake_region_list (reader, "pdl_intake_regions"   );

		// Error check

		check_invariant();

		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public ActionConfigFile unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, ActionConfigFile obj) {

		writer.marshalMapBegin (name);

		if (obj == null) {
			writer.marshalInt (M_TYPE_NAME, MARSHAL_NULL);
		} else {
			writer.marshalInt (M_TYPE_NAME, obj.get_marshal_type());
			obj.do_marshal (writer);
		}

		writer.marshalMapEnd ();

		return;
	}

	// Unmarshal object, polymorphic.

	public static ActionConfigFile unmarshal_poly (MarshalReader reader, String name) {
		ActionConfigFile result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("ActionConfigFile.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_ACTION_CFG:
			result = new ActionConfigFile();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

	// Unmarshal object from a configuration file.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.

	public static ActionConfigFile unmarshal_config (String filename, Class<?> requester) {
		MarshalReader reader = OAFParameterSet.load_file_as_json (filename,requester);
		return (new ActionConfigFile()).unmarshal (reader, null);
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ActionConfigFile : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1
		// Unmarshal from the configuration file, and display it.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ActionConfigFile : Invalid 'test1' subcommand");
				return;
			}

			// Read the configuration file

			ActionConfigFile action_cfg = unmarshal_config ("ActionConfig.json", ActionConfig.class);

			// Display it

			System.out.println (action_cfg.toString());

			return;
		}

		// Subcommand : Test #2
		// Command format:
		//  test2
		// Unmarshal from the configuration file, and display it.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled results.

		if (args[0].equalsIgnoreCase ("test2")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ActionConfigFile : Invalid 'test2' subcommand");
				return;
			}

			// Read the configuration file

			ActionConfigFile action_cfg = unmarshal_config ("ActionConfig.json", ActionConfig.class);

			// Display it

			System.out.println (action_cfg.toString());

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ActionConfigFile.marshal_poly (store, null, action_cfg);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			action_cfg = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			action_cfg = ActionConfigFile.unmarshal_poly (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (action_cfg.toString());

			return;
		}

		// Subcommand : Test #3
		// Command format:
		//  test3
		// Unmarshal from the configuration file, and display it.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled results.
		// This differs from test #2 only in that it uses the non-static marshal methods.

		if (args[0].equalsIgnoreCase ("test3")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ActionConfigFile : Invalid 'test3' subcommand");
				return;
			}

			// Read the configuration file

			ActionConfigFile action_cfg = unmarshal_config ("ActionConfig.json", ActionConfig.class);

			// Display it

			System.out.println (action_cfg.toString());

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			action_cfg.marshal (store, null);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			action_cfg = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			action_cfg = (new ActionConfigFile()).unmarshal (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (action_cfg.toString());

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ActionConfigFile : Unrecognized subcommand : " + args[0]);
		return;

	}

}
