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
 *	"forecast_lags" = [ Array giving a list of time lags at which forecasts are generated, in increasing order.
 *		element = String giving time lag since mainshock, in java.time.Duration format.
 *	]
 *	"forecast_max_delay" = String giving maximum allowed delay for the last forcast, in java.time.Duration format.
 *	"comcat_clock_skew" = Assumed maximum difference between our clock and ComCat clock, in java.time.Duration format.
 *	"comcat_retry_min_gap" = String giving minimum allowed gap between ComCat retries, in java.time.Duration format.
 *	"comcat_retry_lags" = [ Array giving a list of time lags at which ComCat operations are retried, in increasing order.
 *		element = String giving time lag since first attempt, in java.time.Duration format.
 *	]
 */
public class ActionConfigFile {

	//----- Parameter values -----

	// Minimum gap between forecasts, in milliseconds.  Must be positive.

	public long forecast_min_gap;

	// Time lags at which forecasts are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the mainshock.  Must have at least 1 element.
	// The difference between successive elements must be at least the minimum gap.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> forecast_lags;

	// Maximum delay in issuing the final forecast, in milliseconds.  Must be positive.

	public long forecast_max_delay;

	// Assumed maximum difference between our clock and ComCat's clock, in milliseconds.

	public long comcat_clock_skew;

	// Minimum gap between ComCat retries, in milliseconds.  Must be positive.

	public long comcat_retry_min_gap;

	// Time lags at which ComCat retries are generated, in milliseconds.  Must be in increasing order.
	// This is time lag since the initial attempt.
	// The difference between successive elements must be at least the minimum gap.
	// Each element must be a whole number of seconds, between 1 and 10^9 seconds.

	public ArrayList<Long> comcat_retry_lags;


	//----- Construction -----

	// Default constructor.

	public ActionConfigFile () {
		clear();
	}

	// Clear the contents.

	public void clear () {
		forecast_min_gap = 0L;
		forecast_lags = new ArrayList<Long>();
		forecast_max_delay = 0L;
		comcat_clock_skew = 0L;
		comcat_retry_min_gap = 0L;
		comcat_retry_lags = new ArrayList<Long>();
		return;
	}

	// Check that values are valid, throw an exception if not.

	public void check_invariant () {

		if (!( forecast_min_gap > 0L )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_min_gap: " + forecast_min_gap);
		}

		int n = forecast_lags.size();
		
		if (!( n > 0 )) {
				throw new RuntimeException("ActionConfigFile: Empty forecast_lags list");
		}

		long min_lag = 1000L;

		for (int i = 0; i < n; ++i) {
			long forecast_lag = forecast_lags.get(i);
			if (!( forecast_lag >= min_lag && forecast_lag % 1000L == 0L && forecast_lag <= 1000000000000L )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_lag: " + forecast_lag + ", index = " + i);
			}
			min_lag = forecast_lag + forecast_min_gap;
		}

		if (!( forecast_max_delay > 0L )) {
				throw new RuntimeException("ActionConfigFile: Invalid forecast_max_delay: " + forecast_max_delay);
		}

		if (!( comcat_clock_skew >= 0L )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_clock_skew: " + comcat_clock_skew);
		}

		if (!( comcat_retry_min_gap > 0L )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_retry_min_gap: " + comcat_retry_min_gap);
		}

		n = comcat_retry_lags.size();
		min_lag = 1000L;

		for (int i = 0; i < n; ++i) {
			long comcat_retry_lag = comcat_retry_lags.get(i);
			if (!( comcat_retry_lag >= min_lag && comcat_retry_lag % 1000L == 0L && comcat_retry_lag <= 1000000000000L )) {
				throw new RuntimeException("ActionConfigFile: Invalid comcat_retry_lag: " + comcat_retry_lag + ", index = " + i);
			}
			min_lag = comcat_retry_lag + comcat_retry_min_gap;
		}
	
		return;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		result.append ("ActionConfigFile:" + "\n");
		result.append ("forecast_min_gap = " + Duration.ofMillis(forecast_min_gap).toString() + "\n");
		result.append ("forecast_lags = [" + "\n");
		for (int i = 0; i < forecast_lags.size(); ++i) {
			long forecast_lag = forecast_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(forecast_lag).toString() + "\n");
		}
		result.append ("]" + "\n");
		result.append ("forecast_max_delay = " + Duration.ofMillis(forecast_max_delay).toString() + "\n");
		result.append ("comcat_clock_skew = " + Duration.ofMillis(comcat_clock_skew).toString() + "\n");
		result.append ("comcat_retry_min_gap = " + Duration.ofMillis(comcat_retry_min_gap).toString() + "\n");
		result.append ("comcat_retry_lags = [" + "\n");
		for (int i = 0; i < comcat_retry_lags.size(); ++i) {
			long comcat_retry_lag = comcat_retry_lags.get(i);
			result.append ("  " + i + ":  " + Duration.ofMillis(comcat_retry_lag).toString() + "\n");
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

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Error check

		check_invariant();

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		marshal_duration      (writer, "forecast_min_gap"    , forecast_min_gap    );
		marshal_duration_list (writer, "forecast_lags"       , forecast_lags       );
		marshal_duration      (writer, "forecast_max_delay"  , forecast_max_delay  );
		marshal_duration      (writer, "comcat_clock_skew"   , comcat_clock_skew   );
		marshal_duration      (writer, "comcat_retry_min_gap", comcat_retry_min_gap);
		marshal_duration_list (writer, "comcat_retry_lags"   , comcat_retry_lags   );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		forecast_min_gap     = unmarshal_duration      (reader, "forecast_min_gap"    );
		forecast_lags        = unmarshal_duration_list (reader, "forecast_lags"       );
		forecast_max_delay   = unmarshal_duration      (reader, "forecast_max_delay"  );
		comcat_clock_skew    = unmarshal_duration      (reader, "comcat_clock_skew"   );
		comcat_retry_min_gap = unmarshal_duration      (reader, "comcat_retry_min_gap");
		comcat_retry_lags    = unmarshal_duration_list (reader, "comcat_retry_lags"   );

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
