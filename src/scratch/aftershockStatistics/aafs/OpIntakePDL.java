package scratch.aftershockStatistics.aafs;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Date;

import java.text.SimpleDateFormat;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Operation payload for intake of an event as a sync event.
 * Author: Michael Barall 05/20/2018.
 */
public class OpIntakePDL extends DBPayload {

	//----- Constants and variables -----

	// Floating point value to indicate no data.

	public static final double NO_DATA = 1.0e10;

	// Threshold value for test if there is data, values >= NO_DATA_THRESHOLD indicate no data.

	public static final double NO_DATA_THRESHOLD = 0.9e10;

	// PDL command line parameter names

	public static final String PDL_NAME_STATUS = "status";
	public static final String PDL_NAME_ACTION = "action";
	public static final String PDL_NAME_TYPE = "type";
	public static final String PDL_NAME_EVENT_ID = "preferred-eventid";
	public static final String PDL_NAME_MAGNITUDE = "preferred-magnitude";
	public static final String PDL_NAME_LATITUDE = "preferred-latitude";
	public static final String PDL_NAME_LONGITUDE = "preferred-longitude";
	public static final String PDL_NAME_DEPTH = "preferred-depth";
	public static final String PDL_NAME_EVENT_TIME = "preferred-eventtime";

	// PDL status values

	public static final String PDL_STATUS_UPDATE = "UPDATE";
	public static final String PDL_STATUS_DELETE = "DELETE";

	// PDL action values

	public static final String PDL_ACTION_EVENT_ADDED = "EVENT_ADDED";
	public static final String PDL_ACTION_EVENT_UPDATED = "EVENT_UPDATED";
	public static final String PDL_ACTION_EVENT_DELETED = "EVENT_DELETED";
	public static final String PDL_ACTION_EVENT_ARCHIVED = "EVENT_ARCHIVED";
	public static final String PDL_ACTION_PRODUCT_ADDED = "PRODUCT_ADDED";
	public static final String PDL_ACTION_PRODUCT_UPDATED = "PRODUCT_UPDATED";
	public static final String PDL_ACTION_PRODUCT_DELETED = "PRODUCT_DELETED";
	public static final String PDL_ACTION_PRODUCT_ARCHIVED = "PRODUCT_ARCHIVED";

	// PDL type values

	public static final String PDL_TYPE_ORIGIN = "origin";

	// Arguments appearing on the PDL command line.

	public String[] pdl_args;

	// PDL status parameter, or "" if none.

	public String pdl_status;

	// PDL action parameter, or "" if none.

	public String pdl_action;

	// PDL type parameter, or "" if none.

	public String pdl_type;

	// Event ID as received from PDL.

	public String event_id;

	// Mainshock time as received from PDL, or 0L if none, in milliseconds since the epoch.

	public long mainshock_time;

	// Mainshock magnitude as received from PDL, or NO_DATA if not available.

	public double mainshock_mag;

	// Mainshock latitude, in degrees, from -90 to +90, as received from PDL, or NO_DATA if not available.

	public double mainshock_lat;

	// Mainshock longitude, in degrees, from -180 to +180,  as received from PDL, or NO_DATA if not available.

	public double mainshock_lon;

	// Mainshock depth, in kilometers, positive underground,  as received from PDL, or NO_DATA if not available.

	public double mainshock_depth;

	// Parameters supplied by the analyst, or null if none.

	public AnalystOptions analyst_options;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpIntakePDL () {}


	// Set up the contents, for no analyst data.

	public void setup (String[] args, int argix_lo, int argix_hi) {
		parse_pdl_command_line (args, argix_lo, argix_hi);
		analyst_options = null;
		return;
	}


	// Set up the contents, with analyst data

	public void setup (String[] args, int argix_lo, int argix_hi, 
						AnalystOptions the_analyst_options) {
		parse_pdl_command_line (args, argix_lo, argix_hi);
		analyst_options = the_analyst_options;
		return;
	}


	// Parse the PDL command line.

	public void parse_pdl_command_line (String[] args, int argix_lo, int argix_hi) {

		// Set all PDL parameters to empty values

		pdl_status = "";
		pdl_action = "";
		pdl_type = "";
		event_id = "";
		mainshock_time = 0L;
		mainshock_mag = NO_DATA;
		mainshock_lat = NO_DATA;
		mainshock_lon = NO_DATA;
		mainshock_depth = NO_DATA;

		// List of arguments

		ArrayList<String> pdl_arg_list = new ArrayList<String>();

		// Loop to process arguments

		for (int ix = argix_lo; ix < argix_hi; ++ix) {

			// A quoted argument on the command line may be split among multiple arguments here.
			// We reconstruct them by checking if an argument begins with a quote.
			// If so, we concatenate arguments until we see the final quote.

			String param = args[ix].trim();
			if (param.isEmpty()) {
				continue;
			}

			if (param.startsWith ("\"")) {
				param = param.substring (1);		// strip the leading quote
				while (ix + 1 < argix_hi && !(param.endsWith ("\""))) {
					++ix;
					param = param + " " + (args[ix].trim());	// append a space and the next argument
				}
				if (param.endsWith ("\"")) {
					param = param.substring (0, param.length() - 1);	// strip the trailing quote
				}
			}

			// Save the parameter

			pdl_arg_list.add (param);

			// Try to split parameter into name and value

			if (param.startsWith ("--")) {
				int eqix = param.indexOf ("=", 2);
				if (eqix >= 0) {
					String name = param.substring (2, eqix);	// name is part before the equal sign
					while (name.endsWith (" ")) {
						name = name.substring (0, name.length() - 1);	// trim spaces before the equal sign
					}
					String value = param.substring (eqix + 1);	// value is part after the equal sign
					while (value.startsWith (" ")) {
						value = value.substring (1);		// trim spaces after the equal sign
					}

					// Switch over the parameters we want to recognize

					switch (name.toLowerCase()) {

					case PDL_NAME_STATUS:
						pdl_status = value;
						break;

					case PDL_NAME_ACTION:
						pdl_action = value;
						break;

					case PDL_NAME_TYPE:
						pdl_type = value;
						break;

					case PDL_NAME_EVENT_ID:
						event_id = value;
						break;

					case PDL_NAME_MAGNITUDE:
						try {
							mainshock_mag = Double.parseDouble (value);
						} catch (Exception e) {
							mainshock_mag = NO_DATA;
						}
						break;

					case PDL_NAME_LATITUDE:
						try {
							mainshock_lat = Double.parseDouble (value);
						} catch (Exception e) {
							mainshock_lat = NO_DATA;
						}
						break;

					case PDL_NAME_LONGITUDE:
						try {
							mainshock_lon = Double.parseDouble (value);
						} catch (Exception e) {
							mainshock_lon = NO_DATA;
						}
						break;

					case PDL_NAME_DEPTH:
						try {
							mainshock_depth = Double.parseDouble (value);
						} catch (Exception e) {
							mainshock_depth = NO_DATA;
						}
						break;

					case PDL_NAME_EVENT_TIME:
						try {
							SimpleDateFormat fmt = new SimpleDateFormat ("yyyy-MM-dd'T'HH:mm:ss.SSSXXX");
							Date date = fmt.parse (value);
							mainshock_time = date.getTime();
						} catch (Exception e) {
							mainshock_time = 0L;
						}
						break;
					}
				}
			}
		}

		// Get the parameters into an array

		int n = pdl_arg_list.size();
		pdl_args = new String[n];

		for (int i = 0; i < n; ++i) {
			pdl_args[i] = pdl_arg_list.get(i);
		}

		return;
	}


	// Return the effective analyst parameters, or null if none.

	public ForecastParameters get_eff_analyst_params () {
		if (analyst_options != null) {
			return analyst_options.analyst_params;
		}
		return null;
	}


	// Return true if we have an event id.

	public boolean has_event_id () {
		return event_id.length() > 0;
	}


	// Return true if we have latitude, longitude, depth, and magnitude.

	public boolean has_lat_lon_depth_mag () {
		return mainshock_mag < NO_DATA_THRESHOLD
			&& mainshock_lat < NO_DATA_THRESHOLD
			&& mainshock_lon < NO_DATA_THRESHOLD
			&& mainshock_depth < NO_DATA_THRESHOLD;
	}


	// Return true if we have mainshock time.

	public boolean has_event_time () {
		return mainshock_time > 0L;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 31001;

	private static final String M_VERSION_NAME = "OpIntakePDL";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalStringArray               ("pdl_args"           , pdl_args           );
		writer.marshalString                    ("pdl_status"         , pdl_status         );
		writer.marshalString                    ("pdl_action"         , pdl_action         );
		writer.marshalString                    ("pdl_type"           , pdl_type           );
		writer.marshalString                    ("event_id"           , event_id           );
		writer.marshalLong                      ("mainshock_time"     , mainshock_time     );
		writer.marshalDouble                    ("mainshock_mag"      , mainshock_mag      );
		writer.marshalDouble                    ("mainshock_lat"      , mainshock_lat      );
		writer.marshalDouble                    ("mainshock_lon"      , mainshock_lon      );
		writer.marshalDouble                    ("mainshock_depth"    , mainshock_depth    );

		AnalystOptions.marshal_poly     (writer, "analyst_options"    , analyst_options    );

		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Superclass

		super.do_umarshal (reader);

		// Contents

		pdl_args            = reader.unmarshalStringArray               ("pdl_args"           );
		pdl_status          = reader.unmarshalString                    ("pdl_status"         );
		pdl_action          = reader.unmarshalString                    ("pdl_action"         );
		pdl_type            = reader.unmarshalString                    ("pdl_type"           );
		event_id            = reader.unmarshalString                    ("event_id"           );
		mainshock_time      = reader.unmarshalLong                      ("mainshock_time"     );
		mainshock_mag       = reader.unmarshalDouble                    ("mainshock_mag"      );
		mainshock_lat       = reader.unmarshalDouble                    ("mainshock_lat"      );
		mainshock_lon       = reader.unmarshalDouble                    ("mainshock_lon"      );
		mainshock_depth     = reader.unmarshalDouble                    ("mainshock_depth"    );

		analyst_options     = AnalystOptions.unmarshal_poly     (reader, "analyst_options"    );

		return;
	}

	// Marshal object.

	@Override
	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	@Override
	public OpIntakePDL unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a pending task.

	@Override
	public OpIntakePDL unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

}
