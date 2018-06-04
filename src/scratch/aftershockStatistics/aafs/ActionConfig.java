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

	// Get the first element of forecast_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_forecast_lag (long min_lag) {
		return param_set.get_next_forecast_lag (min_lag);
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
