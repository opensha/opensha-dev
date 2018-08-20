package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.ArrayList;

import java.time.Duration;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;

import scratch.aftershockStatistics.OAFParameterSet;

/**
 * Configuration for AAFS server actions.
 * Author: Michael Barall 04/29/2018.
 *
 * To use, create an object of this class, and then call its methods to obtain configuration parameters.
 *
 * Parameters come from a configuration file, in the format of ActionConfigFile.
 */
public class ActionConfig {

	//----- Parameter set -----

	// Cached parameter set.

	private static ActionConfigFile cached_param_set = null;

	// Parameter set.

	private ActionConfigFile param_set;

	// Get the parameter set.

	private static synchronized ActionConfigFile get_param_set () {

		// If we have a cached parameter set, return it

		if (cached_param_set != null) {
			return cached_param_set;
		}

		// Working data

		ActionConfigFile wk_param_set = null;

		// Any error reading the parameters aborts the program

		try {

			// Read the configuation file

			wk_param_set = ActionConfigFile.unmarshal_config ("ActionConfig.json", ActionConfig.class);

		} catch (Exception e) {
			e.printStackTrace();
            System.err.println("ActionConfig: Error loading parameter file ActionConfig.json, unable to continue");
            System.exit(0);
			//throw new RuntimeException("ActionConfig: Error loading parameter file ActionConfig.json", e);
		}

		// Save the parameter set

		cached_param_set = wk_param_set;
		return cached_param_set;
	}

	// unload_data - Remove the cached data from memory.
	// The data will be reloaded the next time one of these objects is created.
	// Any existing objects will continue to use the old data.
	// This makes it possible to load new parameter values without restarting the program.

	public static synchronized void unload_data () {
		cached_param_set = null;
		return;
	}


	//----- Construction -----

	// Default constructor.

	public ActionConfig () {
		param_set = get_param_set ();
	}

	// Display our contents

	@Override
	public String toString() {
		return "ActionConfig:\n" + param_set.toString();
	}


	//----- Parameter access -----

	// Get minimum gap between forecasts, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_forecast_min_gap () {
		return param_set.forecast_min_gap;
	}

	// Get maximum delay in reporting a forecast to PDL, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_forecast_max_delay () {
		return param_set.forecast_max_delay;
	}

	// Get assumed maximum difference between our clock and ComCat's clock, in milliseconds.
	// (Specifically, if an earthquake occurs at time T then it should be visible in
	// ComCat by the time our clock reads T + comcat_clock_skew.)
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_comcat_clock_skew () {
		return param_set.comcat_clock_skew;
	}

	// Get assumed maximum change in mainshock origin time, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_comcat_origin_skew () {
		return param_set.comcat_origin_skew;
	}

	// Get minimum gap between ComCat retries, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_comcat_retry_min_gap () {
		return param_set.comcat_retry_min_gap;
	}

	// Get minimum ComCat retry lag for missing events, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_comcat_retry_missing () {
		return param_set.comcat_retry_missing;
	}

	// Get minimum time after an earthquake at which sequence-specific forecasts can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_seq_spec_min_lag () {
		return param_set.seq_spec_min_lag;
	}

	// Get minimum time after an earthquake at which one-week advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_advisory_dur_week () {
		return param_set.advisory_dur_week;
	}

	// Get minimum time after an earthquake at which one-month advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_advisory_dur_month () {
		return param_set.advisory_dur_month;
	}

	// Get minimum time after an earthquake at which one-year advisories can be generated, in milliseconds.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_advisory_dur_year () {
		return param_set.advisory_dur_year;
	}

	// Get default value of the maximum forecast lag.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_def_max_forecast_lag () {
		return param_set.def_max_forecast_lag;
	}

	// Get forecast lag at which a timeline not passing the intake filter can be withdrawn.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_withdraw_forecast_lag () {
		return param_set.withdraw_forecast_lag;
	}

	// Get option selecting how to handle stale forecasts.
	// (A forecast is stale if another forecast could be issued immediately.)

	public int get_stale_forecast_option () {
		return param_set.stale_forecast_option;
	}

	// Get flag, indicating if stale forecasts should be skipped.

	public boolean get_skip_stale_forecasts () {
		return param_set.stale_forecast_option == ActionConfigFile.SFOPT_SKIP;
	}

	// Get flag, indicating if stale forecasts should be omitted.

	public boolean get_omit_stale_forecasts () {
		return param_set.stale_forecast_option == ActionConfigFile.SFOPT_OMIT;
	}

	// Get search radius to search for shadowing events, in km.

	public double get_shadow_search_radius () {
		return param_set.shadow_search_radius;
	}

	// Get amount of time to look back from the mainshock, to search for shadowing events.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_shadow_lookback_time () {
		return param_set.shadow_lookback_time;
	}

	// Get minimum magnitude to use for computing centroids, when searching for shadowing events.

	public double get_shadow_centroid_mag () {
		return param_set.shadow_centroid_mag;
	}

	// Get minimum magnitude for a candidate shadowing event to be considered large.

	public double get_shadow_large_mag () {
		return param_set.shadow_large_mag;
	}

	// Get period for the short polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_short_period () {
		return param_set.poll_short_period;
	}

	// Get lookback time for the short polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_short_lookback () {
		return param_set.poll_short_lookback;
	}

	// Get time gap between intake actions for the short polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_short_intake_gap () {
		return param_set.poll_short_intake_gap;
	}

	// Get period for the long polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_long_period () {
		return param_set.poll_long_period;
	}

	// Get lookback time for the long polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_long_lookback () {
		return param_set.poll_long_lookback;
	}

	// Get time gap between intake actions for the long polling cycle.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_poll_long_intake_gap () {
		return param_set.poll_long_intake_gap;
	}

	// Get maximum allowed age for PDL intake.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_pdl_intake_max_age () {
		return param_set.pdl_intake_max_age;
	}

	// Get maximum allowed time in future for PDL intake.
	// Must be a whole number of seconds, between 1 and 10^9 seconds.

	public long get_pdl_intake_max_future () {
		return param_set.pdl_intake_max_future;
	}

	// Get default value of injectable text for PDL JSON files, or "" for none.

	public String get_def_injectable_text () {
		return param_set.def_injectable_text;
	}

	// Get the first element of forecast_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_forecast_lag (long min_lag) {
		return param_set.get_next_forecast_lag (min_lag);
	}

	// Get the first element of forecast_lags that is >= the supplied min_lag and <= the supplied max_lag.
	// The return is -1 if there is no element in the given range.
	// If max_lag <= 0, then def_max_forecast_lag is used as the upper bound.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_forecast_lag (long min_lag, long max_lag) {
		return param_set.get_next_forecast_lag (min_lag, max_lag);
	}

	// Get the first element of comcat_retry_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_comcat_retry_lag (long min_lag) {
		return param_set.get_next_comcat_retry_lag (min_lag);
	}

	// Get the first element of comcat_intake_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_comcat_intake_lag (long min_lag) {
		return param_set.get_next_comcat_intake_lag (min_lag);
	}

	// Get the first element of pdl_report_retry_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_pdl_report_retry_lag (long min_lag) {
		return param_set.get_next_pdl_report_retry_lag (min_lag);
	}

	// Get the pdl intake region that satisfies the min_mag criterion.
	// If found, the region is returned.
	// If not found, null is returned.

	public IntakeSphRegion get_pdl_intake_region_for_min_mag (double lat, double lon, double mag) {
		return param_set.get_pdl_intake_region_for_min_mag (lat, lon, mag);
	}

	// Get the pdl intake region that satisfies the intake_mag criterion.
	// If found, the region is returned.
	// If not found, null is returned.

	public IntakeSphRegion get_pdl_intake_region_for_intake_mag (double lat, double lon, double mag) {
		return param_set.get_pdl_intake_region_for_intake_mag (lat, lon, mag);
	}

	// Get the minimum magnitude for the min_mag criterion in any intake region.
	// The result is 10.0 if there are no intake regions.

	public double get_pdl_intake_region_min_min_mag () {
		return param_set.get_pdl_intake_region_min_min_mag ();
	}

	// Get the minimum magnitude for the intake_mag criterion in any intake region.
	// The result is 10.0 if there are no intake regions.

	public double get_pdl_intake_region_min_intake_mag () {
		return param_set.get_pdl_intake_region_min_intake_mag ();
	}


	//----- Service functions -----

	private static final long UNIT_LAG = 1000L;

	// Convert a lag value so it can be stored in an int.
	// Note: This should be applied only to timing parameters returned by functions
	// in this class, which are guaranteed to be multiples of UNIT_LAG, and (after
	// division by UNIT_LAG) representable as an int.

	public int lag_to_int (long lag) {
		return (int)(lag / UNIT_LAG);
	}

	// Convert an int to a lag value.

	public long int_to_lag (int k) {
		return ((long)k) * UNIT_LAG;
	}

	// Given a lag value, round it down to a multiple of the lag unit (which is 1 second).
	// Note: Assumes lag >= 0L.

	public long floor_unit_lag (long lag) {
		return lag - (lag % UNIT_LAG);
	}

	// Same, except make the return value greater than prior_lag.

	public long floor_unit_lag (long lag, long prior_lag) {
		long x = Math.max (lag, prior_lag + UNIT_LAG);
		return x - (x % UNIT_LAG);
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ActionConfig : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1
		// Create an object, and display the parameters.
		// Then read out the time lag lists.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ActionConfig : Invalid 'test1' subcommand");
				return;
			}

			// Create a configuration object

			ActionConfig action_config = new ActionConfig();

			// Display it

			System.out.println (action_config.toString());

			// Display some calculated values

			System.out.println ("");
			System.out.println ("skip_stale_forecasts = " + action_config.get_skip_stale_forecasts());
			System.out.println ("omit_stale_forecasts = " + action_config.get_omit_stale_forecasts());
			System.out.println ("pdl_intake_region_min_min_mag = " + action_config.get_pdl_intake_region_min_min_mag());
			System.out.println ("pdl_intake_region_min_intake_mag = " + action_config.get_pdl_intake_region_min_intake_mag());

			// Display list of forecast time lags

			System.out.println ("");

			long min_lag = 0L;
			for (;;) {
				long forecast_lag = action_config.get_next_forecast_lag (min_lag);
				if (forecast_lag < 0L) {
					break;
				}
				System.out.println (Duration.ofMillis(forecast_lag).toString() + "  " + forecast_lag);
				min_lag = forecast_lag + action_config.get_forecast_min_gap ();
			}

			// Display list of forecast time lags, to default limit only

			System.out.println ("");

			min_lag = 0L;
			for (;;) {
				long forecast_lag = action_config.get_next_forecast_lag (min_lag, 0L);
				if (forecast_lag < 0L) {
					break;
				}
				System.out.println (Duration.ofMillis(forecast_lag).toString() + "  " + forecast_lag);
				min_lag = forecast_lag + action_config.get_forecast_min_gap ();
			}

			// Display list of ComCat retry time lags

			System.out.println ("");

			min_lag = 0L;
			for (;;) {
				long comcat_retry_lag = action_config.get_next_comcat_retry_lag (min_lag);
				if (comcat_retry_lag < 0L) {
					break;
				}
				System.out.println (Duration.ofMillis(comcat_retry_lag).toString() + "  " + comcat_retry_lag);
				min_lag = comcat_retry_lag + action_config.get_comcat_retry_min_gap ();
			}

			// Display list of ComCat intake time lags

			System.out.println ("");

			min_lag = 0L;
			for (;;) {
				long comcat_intake_lag = action_config.get_next_comcat_intake_lag (min_lag);
				if (comcat_intake_lag < 0L) {
					break;
				}
				System.out.println (Duration.ofMillis(comcat_intake_lag).toString() + "  " + comcat_intake_lag);
				min_lag = comcat_intake_lag + action_config.get_comcat_retry_min_gap ();
			}

			// Display list of PDL report retry time lags

			System.out.println ("");

			min_lag = 0L;
			for (;;) {
				long pdl_report_retry_lag = action_config.get_next_pdl_report_retry_lag (min_lag);
				if (pdl_report_retry_lag < 0L) {
					break;
				}
				System.out.println (Duration.ofMillis(pdl_report_retry_lag).toString() + "  " + pdl_report_retry_lag);
				min_lag = pdl_report_retry_lag + action_config.get_comcat_retry_min_gap ();
			}

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ActionConfig : Unrecognized subcommand : " + args[0]);
		return;

	}

}
