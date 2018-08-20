package scratch.aftershockStatistics.aafs;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import scratch.aftershockStatistics.util.EventNotFoundException;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;

import scratch.aftershockStatistics.AftershockStatsCalc;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.ComcatException;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

/**
 * Mainshock properties a forecast.
 * Author: Michael Barall 07/29/2018.
 *
 * All fields are public, since there is little benefit to having lots of getters and setters.
 *
 * This code was originally part of ForecastParameters and was split out to a separate class.
 */
public class ForecastMainshock {

	//----- Root parameters -----

	// Mainshock event id, used for sending query to Comcat.

	public String query_event_id = "";


	//----- Mainshock parameters -----

	// Mainshock parameter available flag.

	public boolean mainshock_avail = false;

	// Mainshock event id, as received from Comcat.

	public String mainshock_event_id = null;

	// Mainshock network.

	public String mainshock_network = null;

	// Mainshock code.

	public String mainshock_code = null;

	// Mainshock id list;

	public String[] mainshock_id_list = null;

	// Mainshock time, in milliseconds since the epoch.

	public long mainshock_time = 0L;

	// Mainshock magnitude.

	public double mainshock_mag = 0.0;

	// Mainshock latitude, in degrees, from -90 to +90.

	public double mainshock_lat = 0.0;

	// Mainshock longitude, in degrees, from -180 to +180.

	public double mainshock_lon = 0.0;

	// Mainshock depth, in kilometers, positive underground.

	public double mainshock_depth = 0.0;

	// Set mainshock parameters to default.

	public void set_default_mainshock_params () {
		mainshock_event_id = null;
		mainshock_network = null;
		mainshock_code = null;
		mainshock_id_list = null;
		mainshock_time = 0L;
		mainshock_mag = 0.0;
		mainshock_lat = 0.0;
		mainshock_lon = 0.0;
		mainshock_depth = 0.0;
		return;
	}

	// Set mainshock parameters from rupture information.

	public void set_eqk_rupture (ObsEqkRupture rup) {

		mainshock_event_id = rup.getEventId();
		if (mainshock_event_id == null || mainshock_event_id.isEmpty()) {
			throw new EventNotFoundException ("ForecastMainshock.set_eqk_rupture: Comcat did not return event id: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		}

		Map<String, String> eimap = ComcatAccessor.extendedInfoToMap (rup, ComcatAccessor.EITMOPT_NULL_TO_EMPTY);
		mainshock_network = eimap.get (ComcatAccessor.PARAM_NAME_NETWORK);
		if (mainshock_network == null || mainshock_network.isEmpty()) {
			throw new EventNotFoundException ("ForecastMainshock.set_eqk_rupture: Comcat did not return event network: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		}
		mainshock_code = eimap.get (ComcatAccessor.PARAM_NAME_CODE);
		if (mainshock_code == null || mainshock_code.isEmpty()) {
			throw new EventNotFoundException ("ForecastMainshock.set_eqk_rupture: Comcat did not return event code: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		}

		String comcat_idlist = eimap.get (ComcatAccessor.PARAM_NAME_IDLIST);
		//if (comcat_idlist == null || comcat_idlist.isEmpty()) {
		//	throw new EventNotFoundException ("ForecastMainshock.set_eqk_rupture: Comcat did not return event id list: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		//}
		List<String> idlist = ComcatAccessor.idsToList (comcat_idlist, mainshock_event_id);
		if (idlist.isEmpty()) {
			throw new EventNotFoundException ("ForecastMainshock.set_eqk_rupture: Unable to construct event id list: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		}
		mainshock_id_list = idlist.toArray (new String[0]);

		mainshock_time = rup.getOriginTime();
		mainshock_mag = rup.getMag();
		Location hypo = rup.getHypocenterLocation();
		mainshock_lat = hypo.getLatitude();
		mainshock_lon = hypo.getLongitude();
		mainshock_depth = hypo.getDepth();

		if (mainshock_lon > 180.0) {
			mainshock_lon -= 360.0;
		}
		if (mainshock_lon < -180.0) {
			mainshock_lon = 180.0;
		}

		if (mainshock_lat > 90.0) {
			mainshock_lat = 90.0;
		}
		else if (mainshock_lat < -90.0) {
			mainshock_lat = -90.0;
		}

		return;
	}

	// Get location from mainshock parameters.

	public Location get_eqk_location () {
		return new Location (mainshock_lat, mainshock_lon, mainshock_depth);
	}

	// Get spherical location from mainshock parameters.

	public SphLatLon get_sph_eqk_location () {
		return new SphLatLon (mainshock_lat, mainshock_lon);
	}

	// Get rupture information from mainshock parameters.

	public ObsEqkRupture get_eqk_rupture () {
		return new ObsEqkRupture (mainshock_event_id, mainshock_time, get_eqk_location(), mainshock_mag);
	}

	// Fetch mainshock parameters.
	// Note: Must set query_event_id first.

	public void fetch_mainshock_params (ForecastMainshock prior_params) {

		// Fetch parameters from Comcat

		ObsEqkRupture rup;

		try {
			ComcatAccessor accessor = new ComcatAccessor();
			rup = accessor.fetchEvent(query_event_id, false, true);		// request extended info
		} catch (Exception e) {
			throw new ComcatException ("ForecastMainshock.fetch_mainshock_params: Comcat exception: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id), e);
		}

		if (rup == null) {
			throw new EventNotFoundException ("ForecastMainshock.fetch_mainshock_params: Comcat did not find event: query_event_id = " + ((query_event_id == null) ? "null" : query_event_id));
		}

		set_eqk_rupture (rup);
		mainshock_avail = true;
		return;
	}


	//----- Transient parameters -----

	// The timeline ID associated with this query, or null if none.
	// Note: This parameter is not marshaled/unmarshaled.

	public String timeline_id = null;

	// Set transient parameters to default.

	public void set_default_transient_params () {
		timeline_id = null;
		return;
	}


	//----- Construction -----

	// Default constructor.

	public ForecastMainshock () {}

	// Set up the mainshock parameters, after setting everything else to default.
	// This version always throws an exception if the event is not successfully fetched.

	public void setup_mainshock_only (String the_query_event_id) {
		setup_all_default (the_query_event_id);

		fetch_mainshock_params (null);
	
		return;
	}

	// Set up the mainshock parameters, after setting everything else to default.
	// This version does not throw EventNotFoundException.
	// If the event is not found, the function returns with mainshock_avail set to false.
	// The return value is mainshock_avail.

	public boolean setup_mainshock_poll (String the_query_event_id) {
		setup_all_default (the_query_event_id);

		try {
			fetch_mainshock_params (null);
		} catch (EventNotFoundException e) {
			mainshock_avail = false;
		}
	
		return mainshock_avail;
	}

	// Set everything to default.

	public void setup_all_default (String the_query_event_id) {
		query_event_id = the_query_event_id;

		mainshock_avail = false;
		set_default_mainshock_params();

		set_default_transient_params();
	
		return;
	}

	// Copy from another object.

	public void copy_from (ForecastMainshock other) {
		mainshock_avail = other.mainshock_avail;
		mainshock_event_id = other.mainshock_event_id;
		mainshock_network = other.mainshock_network;
		mainshock_code = other.mainshock_code;
		mainshock_id_list = ((other.mainshock_id_list == null) ? other.mainshock_id_list : other.mainshock_id_list.clone());
		mainshock_time = other.mainshock_time;
		mainshock_mag = other.mainshock_mag;
		mainshock_lat = other.mainshock_lat;
		mainshock_lon = other.mainshock_lon;
		mainshock_depth = other.mainshock_depth;
		timeline_id = other.timeline_id;
		return;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();

		result.append ("ForecastMainshock:" + "\n");

		result.append ("query_event_id = " + query_event_id + "\n");

		result.append ("mainshock_avail = " + mainshock_avail + "\n");
		if (mainshock_avail) {
			result.append ("mainshock_event_id = " + mainshock_event_id + "\n");
			result.append ("mainshock_network = " + mainshock_network + "\n");
			result.append ("mainshock_code = " + mainshock_code + "\n");
			result.append ("mainshock_id_list = " + Arrays.toString (mainshock_id_list) + "\n");
			result.append ("mainshock_time = " + mainshock_time + "\n");
			result.append ("mainshock_mag = " + mainshock_mag + "\n");
			result.append ("mainshock_lat = " + mainshock_lat + "\n");
			result.append ("mainshock_lon = " + mainshock_lon + "\n");
			result.append ("mainshock_depth = " + mainshock_depth + "\n");
		}

		if (timeline_id != null) {
			result.append ("timeline_id = " + timeline_id + "\n");
		}

		return result.toString();
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 39001;

	private static final String M_VERSION_NAME = "ForecastMainshock";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 39000;
	protected static final int MARSHAL_FCAST_MAIN = 39001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_FCAST_MAIN;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalString ("query_event_id" , query_event_id );

		writer.marshalBoolean ("mainshock_avail"     , mainshock_avail     );
		if (mainshock_avail) {
			writer.marshalString      ("mainshock_event_id", mainshock_event_id);
			writer.marshalString      ("mainshock_network" , mainshock_network );
			writer.marshalString      ("mainshock_code"    , mainshock_code    );
			writer.marshalStringArray ("mainshock_id_list" , mainshock_id_list );
			writer.marshalLong        ("mainshock_time"    , mainshock_time    );
			writer.marshalDouble      ("mainshock_mag"     , mainshock_mag     );
			writer.marshalDouble      ("mainshock_lat"     , mainshock_lat     );
			writer.marshalDouble      ("mainshock_lon"     , mainshock_lon     );
			writer.marshalDouble      ("mainshock_depth"   , mainshock_depth   );
		}
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		query_event_id  = reader.unmarshalString ("query_event_id" );

		mainshock_avail      = reader.unmarshalBoolean ("mainshock_avail");
		if (mainshock_avail) {
			mainshock_event_id = reader.unmarshalString      ("mainshock_event_id");
			mainshock_network  = reader.unmarshalString      ("mainshock_network" );
			mainshock_code     = reader.unmarshalString      ("mainshock_code"    );
			mainshock_id_list  = reader.unmarshalStringArray ("mainshock_id_list" );
			mainshock_time     = reader.unmarshalLong        ("mainshock_time"    );
			mainshock_mag      = reader.unmarshalDouble      ("mainshock_mag"     );
			mainshock_lat      = reader.unmarshalDouble      ("mainshock_lat"     );
			mainshock_lon      = reader.unmarshalDouble      ("mainshock_lon"     );
			mainshock_depth    = reader.unmarshalDouble      ("mainshock_depth"   );
		} else {
			set_default_mainshock_params();
		}

		set_default_transient_params();

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

	public ForecastMainshock unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, ForecastMainshock obj) {

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

	public static ForecastMainshock unmarshal_poly (MarshalReader reader, String name) {
		ForecastMainshock result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("ForecastMainshock.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_FCAST_MAIN:
			result = new ForecastMainshock();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ForecastMainshock : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  query_event_id
		// Get info for the event, and display it.

		if (args[0].equalsIgnoreCase ("test1")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastMainshock : Invalid 'test1' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_query_event_id);

			// Display it

			System.out.println ("");
			System.out.println (fcmain.toString());

			return;
		}




		// Subcommand : Test #2
		// Command format:
		//  test2  query_event_id
		// Get info for the event, and display it.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled info.

		if (args[0].equalsIgnoreCase ("test2")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastMainshock : Invalid 'test2' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_query_event_id);

			// Display it

			System.out.println ("");
			System.out.println (fcmain.toString());

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ForecastMainshock.marshal_poly (store, null, fcmain);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			fcmain = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			fcmain = ForecastMainshock.unmarshal_poly (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (fcmain.toString());

			return;
		}




		// Subcommand : Test #3
		// Command format:
		//  test3  query_event_id
		// Get info for the event, and display it.
		// This version does not throw an exception if the event is not found.

		if (args[0].equalsIgnoreCase ("test3")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastMainshock : Invalid 'test3' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_poll (the_query_event_id);

			// Display it

			System.out.println ("");
			System.out.println (fcmain.toString());

			return;
		}




		// Subcommand : Test #4
		// Command format:
		//  test4  query_event_id
		// Get info for the event, and display it.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled info.
		// This version does not throw an exception if the event is not found.

		if (args[0].equalsIgnoreCase ("test4")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastMainshock : Invalid 'test4' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_poll (the_query_event_id);

			// Display it

			System.out.println ("");
			System.out.println (fcmain.toString());

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ForecastMainshock.marshal_poly (store, null, fcmain);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			fcmain = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			fcmain = ForecastMainshock.unmarshal_poly (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (fcmain.toString());

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("ForecastMainshock : Unrecognized subcommand : " + args[0]);
		return;

	}

}
