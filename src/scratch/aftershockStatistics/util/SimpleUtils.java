package scratch.aftershockStatistics.util;

import java.io.StringWriter;
import java.io.PrintWriter;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;


/**
 * Class to hold some simple utility functions.
 * Author: Michael Barall 05/29/2018.
 *
 * All functions in this class are static.
 */
public class SimpleUtils {




	// Get a stack trace as a string.

	public static String getStackTraceAsString (Throwable e) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw, true);
		e.printStackTrace(pw);
		return sw.getBuffer().toString();
	}




	// Convert a time (in milliseconds after the epoch) to a human-readable string.

	public static String time_to_string (long the_time) {
		SimpleDateFormat fmt = new SimpleDateFormat ("yyyy-MM-dd HH:mm:ss z");
		fmt.setTimeZone (TimeZone.getTimeZone ("UTC"));
		return fmt.format (new Date (the_time));
	}




	// Convert a time (in milliseconds after the epoch) to a human-readable string.
	// This version does not have the "UTC" suffix (but the time is still UTC).

	public static String time_to_string_no_z (long the_time) {
		SimpleDateFormat fmt = new SimpleDateFormat ("yyyy-MM-dd HH:mm:ss");
		fmt.setTimeZone (TimeZone.getTimeZone ("UTC"));
		return fmt.format (new Date (the_time));
	}




	// Given a time (in milliseconds after the epoch), produce a string which
	// is its numerical value followed by the human-readable form in parentheses.

	public static String time_raw_and_string (long the_time) {
		return the_time + " (" + time_to_string(the_time) + ")";
	}




	// Given information about an event, produce a one-line summary.

	public static String event_info_one_line (long event_time,
			double event_mag, double event_lat, double event_lon, double event_depth) {

		StringBuilder result = new StringBuilder();

		result.append ("time = " + time_to_string (event_time));

		if (event_mag > -99.0 && event_mag < 99.0) {
			result.append (String.format(", mag = %.3f", event_mag));
		} else {
			result.append (", mag = " + event_mag);
		}

		if (event_lat > -999.0 && event_lat < 999.0) {
			result.append (String.format(", lat = %.5f", event_lat));
		} else {
			result.append (", lat = " + event_lat);
		}

		if (event_lon > -999.0 && event_lon < 999.0) {
			result.append (String.format(", lon = %.5f", event_lon));
		} else {
			result.append (", lon = " + event_lon);
		}

		if (event_depth > -9999.0 && event_depth < 9999.0) {
			result.append (String.format(", depth = %.3f", event_depth));
		} else {
			result.append (", depth = " + event_depth);
		}

		return result.toString();
	}




	// Given information about an event, produce a one-line summary.
	// This produces an event ID, followed by info in parentheses

	public static String event_id_and_info_one_line (String event_id, long event_time,
			double event_mag, double event_lat, double event_lon, double event_depth) {

		return event_id
				+ " ("
				+ event_info_one_line (event_time, event_mag, event_lat, event_lon, event_depth)
				+ ")";
	}




}
