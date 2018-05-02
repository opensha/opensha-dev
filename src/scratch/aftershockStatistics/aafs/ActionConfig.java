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

	// Get minimum gap between forecasts, in milliseconds.  Must be positive.

	public long get_forecast_min_gap () {
		return param_set.forecast_min_gap;
	}

	// Get the first element of forecast_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_forecast_lag (long min_lag) {
		return param_set.get_next_forecast_lag (min_lag);
	}

	// Get maximum delay in issuing the final forecast, in milliseconds.  Must be positive.

	public long get_forecast_max_delay () {
		return param_set.forecast_max_delay;
	}

	// Get assumed maximum difference between our clock and ComCat's clock, in milliseconds.

	public long get_comcat_clock_skew () {
		return param_set.comcat_clock_skew;
	}

	// Get minimum gap between ComCat retries, in milliseconds.  Must be positive.

	public long get_comcat_retry_min_gap () {
		return param_set.comcat_retry_min_gap;
	}

	// Get the first element of comcat_retry_lags that is >= the supplied min_lag.
	// The return is -1 if the supplied min_lag is greater than all elements.
	// If a value is found, it is guaranteed to be a whole number of seconds, from 1 to 10^9 seconds.

	public long get_next_comcat_retry_lag (long min_lag) {
		return param_set.get_next_comcat_retry_lag (min_lag);
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

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ActionConfig : Unrecognized subcommand : " + args[0]);
		return;

	}

}
