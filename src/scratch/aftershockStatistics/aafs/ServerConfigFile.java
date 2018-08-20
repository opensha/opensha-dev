package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;
import scratch.aftershockStatistics.util.TimeSplitOutputStream;

import scratch.aftershockStatistics.OAFParameterSet;

import scratch.aftershockStatistics.pdl.PDLSenderConfig;

/**
 * Configuration file for AAFS server.
 * Author: Michael Barall 06/03/2018.
 *
 * All fields are public, since this is just a buffer for reading and writing files.
 *
 * JSON file format:
 *
 *	"ServerConfigFile" = Integer giving file version number, currently 34001.
 *	"db_host" = String giving database host name or IP address.
 *	"db_port" = Integer giving database port number.
 *	"db_name" = String giving database name.  Used for both database access and user authentication.
 *	"db_user" = String giving database user name.  This name provides read/write access.
 *	"db_password" = String giving database password.
 *	"activemq_host" = String giving ActiveMQ host name or IP address.
 *	"activemq_port" = Integer giving ActiveMQ port number.
 *	"activemq_user" = String giving ActiveMQ user name.
 *	"activemq_password" = String giving ActiveMQ password.
 *	"log_con_aafs" = String giving pattern for AAFS console log filenames, in the format of SimpleDateFormat, or "" if none.
 *	"log_con_intake" = String giving pattern for intake console log filenames, in the format of SimpleDateFormat, or "" if none.
 *	"log_con_control" = String giving pattern for control console log filenames, in the format of SimpleDateFormat, or "" if none.
 *	"log_summary" = String giving pattern for summary log filenames, in the format of SimpleDateFormat, or "" if none.
 *	"comcat_url" = String giving Comcat URL.
 *	"pdl_enable" = Integer giving PDL enable option: 0 = none, 1 = development, 2 = production.
 *	"pdl_key_filename" = String giving PDL signing key filename, can be empty string for none.
 *	"pdl_dev_senders" = [ Array giving a list of PDL sender configurations for development PDL, in priority order.
 *		element = { Structure giving PDL server configuration.
 *			"host" = String giving PDL sender host name or IP address.
 *			"port" = Integer giving PDL sender port number.
 *			"connectTimeout" = Integer giving PDL sender connection timeout, in milliseconds.
 *		}
 *	]
 *	"pdl_prod_senders" = [ Array giving a list of PDL sender configurations for production PDL, in priority order.
 *		element = { Structure giving PDL server configuration.
 *			"host" = String giving PDL sender host name or IP address.
 *			"port" = Integer giving PDL sender port number.
 *			"connectTimeout" = Integer giving PDL sender connection timeout, in milliseconds.
 *		}
 *	]
 */
public class ServerConfigFile {

	//----- Parameter values -----

	// Database host name or IP address.

	public String db_host;

	// Database port number.

	public int db_port;

	// Database name.  Used for both database access and user authentication.

	public String db_name;

	// Database user name.  This name provides read/write access.

	public String db_user;

	// Database password.

	public String db_password;

	// ActiveMQ host name or IP address.

	public String activemq_host;

	// ActiveMQ port number.

	public int activemq_port;

	// ActiveMQ user name.

	public String activemq_user;

	// ActiveMQ password.

	public String activemq_password;

	// Pattern for AAFS console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String log_con_aafs;

	// Pattern for intake console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String log_con_intake;

	// Pattern for control console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String log_con_control;

	// Pattern for summary log filenames, in the format of SimpleDateFormat, or "" if none.

	public String log_summary;

	// Comcat URL.

	public String comcat_url;

	// PDL enable option.

	public static final int PDLOPT_MIN = 0;
	public static final int PDLOPT_NONE = 0;		// No PDL access
	public static final int PDLOPT_DEV = 1;			// PDL development server
	public static final int PDLOPT_PROD = 2;		// PDL production server
	public static final int PDLOPT_MAX = 2;

	public static final int PDLOPT_UNSPECIFIED = -1;	// PDL access is unspecified

	public int pdl_enable;

	// PDL signing key filename, can be empty string for none.  Cannot be null.

	public String pdl_key_filename;

	// List of PDL development senders.

	public ArrayList<PDLSenderConfig> pdl_dev_senders;

	// List of PDL production senders.

	public ArrayList<PDLSenderConfig> pdl_prod_senders;


	//----- Construction -----

	// Default constructor.

	public ServerConfigFile () {
		clear();
	}

	// Clear the contents.

	public void clear () {
		db_host = "";
		db_port = 0;
		db_name = "";
		db_user = "";
		db_password = "";
		activemq_host = "";
		activemq_port = 0;
		activemq_user = "";
		activemq_password = "";
		log_con_aafs = "";
		log_con_intake = "";
		log_con_control = "";
		log_summary = "";
		comcat_url = "";
		pdl_enable = PDLOPT_NONE;
		pdl_key_filename = "";
		pdl_dev_senders = new ArrayList<PDLSenderConfig>();
		pdl_prod_senders = new ArrayList<PDLSenderConfig>();
		return;
	}

	// Check that values are valid, throw an exception if not.

	public void check_invariant () {

		if (!( db_host != null && db_host.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid db_host: " + ((db_host == null) ? "<null>" : db_host));
		}

		if (!( db_port >= 1024 && db_port <= 65535 )) {
			throw new RuntimeException("ServerConfigFile: Invalid db_port: " + db_port);
		}

		if (!( db_name != null && db_name.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid db_name: " + ((db_name == null) ? "<null>" : db_name));
		}

		if (!( db_user != null && db_user.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid db_user: " + ((db_user == null) ? "<null>" : db_user));
		}

		if (!( db_password != null && db_password.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid db_password: " + ((db_password == null) ? "<null>" : db_password));
		}

		if (!( activemq_host != null && activemq_host.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid activemq_host: " + ((activemq_host == null) ? "<null>" : activemq_host));
		}

		if (!( activemq_port >= 1024 && activemq_port <= 65535 )) {
			throw new RuntimeException("ServerConfigFile: Invalid activemq_port: " + activemq_port);
		}

		if (!( activemq_user != null && activemq_user.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid activemq_user: " + ((activemq_user == null) ? "<null>" : activemq_user));
		}

		if (!( activemq_password != null && activemq_password.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid activemq_password: " + ((activemq_password == null) ? "<null>" : activemq_password));
		}

		if (!( log_con_aafs != null && (log_con_aafs.isEmpty() || TimeSplitOutputStream.is_valid_pattern(log_con_aafs)) )) {
			throw new RuntimeException("ServerConfigFile: Invalid log_con_aafs: " + ((log_con_aafs == null) ? "<null>" : log_con_aafs));
		}

		if (!( log_con_intake != null && (log_con_intake.isEmpty() || TimeSplitOutputStream.is_valid_pattern(log_con_intake)) )) {
			throw new RuntimeException("ServerConfigFile: Invalid log_con_intake: " + ((log_con_intake == null) ? "<null>" : log_con_intake));
		}

		if (!( log_con_control != null && (log_con_control.isEmpty() || TimeSplitOutputStream.is_valid_pattern(log_con_control)) )) {
			throw new RuntimeException("ServerConfigFile: Invalid log_con_control: " + ((log_con_control == null) ? "<null>" : log_con_control));
		}

		if (!( log_summary != null && (log_summary.isEmpty() || TimeSplitOutputStream.is_valid_pattern(log_summary)) )) {
			throw new RuntimeException("ServerConfigFile: Invalid log_summary: " + ((log_summary == null) ? "<null>" : log_summary));
		}

		if (!( comcat_url != null && comcat_url.trim().length() > 0 )) {
			throw new RuntimeException("ServerConfigFile: Invalid comcat_url: " + ((comcat_url == null) ? "<null>" : comcat_url));
		}

		if (!( pdl_enable >= PDLOPT_MIN && pdl_enable <= PDLOPT_MAX )) {
			throw new RuntimeException("ServerConfigFile: Invalid pdl_enable: " + pdl_enable);
		}

		if (!( pdl_key_filename != null )) {
			throw new RuntimeException("ServerConfigFile: Invalid pdl_key_filename: " + ((pdl_key_filename == null) ? "<null>" : pdl_key_filename));
		}

		if (!( pdl_dev_senders != null )) {
			throw new RuntimeException("ServerConfigFile: pdl_dev_senders list is null");
		}

		if ( pdl_enable == PDLOPT_DEV && pdl_dev_senders.size() == 0 ) {
			throw new RuntimeException("ServerConfigFile: pdl_dev_senders is empty, but pdl_enable = " + pdl_enable);
		}

		if (!( pdl_prod_senders != null )) {
			throw new RuntimeException("ServerConfigFile: pdl_prod_senders list is null");
		}

		if ( pdl_enable == PDLOPT_PROD && pdl_prod_senders.size() == 0 ) {
			throw new RuntimeException("ServerConfigFile: pdl_prod_senders is empty, but pdl_enable = " + pdl_enable);
		}

		return;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		result.append ("ServerConfigFile:" + "\n");

		result.append ("db_host = " + ((db_host == null) ? "<null>" : db_host) + "\n");
		result.append ("db_port = " + db_port + "\n");
		result.append ("db_name = " + ((db_name == null) ? "<null>" : db_name) + "\n");
		result.append ("db_user = " + ((db_user == null) ? "<null>" : db_user) + "\n");
		result.append ("db_password = " + ((db_password == null) ? "<null>" : db_password) + "\n");
		result.append ("activemq_host = " + ((activemq_host == null) ? "<null>" : activemq_host) + "\n");
		result.append ("activemq_port = " + activemq_port + "\n");
		result.append ("activemq_user = " + ((activemq_user == null) ? "<null>" : activemq_user) + "\n");
		result.append ("activemq_password = " + ((activemq_password == null) ? "<null>" : activemq_password) + "\n");
		result.append ("log_con_aafs = " + ((log_con_aafs == null) ? "<null>" : log_con_aafs) + "\n");
		result.append ("log_con_intake = " + ((log_con_intake == null) ? "<null>" : log_con_intake) + "\n");
		result.append ("log_con_control = " + ((log_con_control == null) ? "<null>" : log_con_control) + "\n");
		result.append ("log_summary = " + ((log_summary == null) ? "<null>" : log_summary) + "\n");
		result.append ("comcat_url = " + ((comcat_url == null) ? "<null>" : comcat_url) + "\n");
		result.append ("pdl_enable = " + pdl_enable + "\n");
		result.append ("pdl_key_filename = " + ((pdl_key_filename == null) ? "<null>" : pdl_key_filename) + "\n");

		result.append ("pdl_dev_senders = [" + "\n");
		for (int i = 0; i < pdl_dev_senders.size(); ++i) {
			PDLSenderConfig pdl_sender = pdl_dev_senders.get(i);
			result.append ("  " + i + ":  " + pdl_sender.toString() + "\n");
		}
		result.append ("]" + "\n");

		result.append ("pdl_prod_senders = [" + "\n");
		for (int i = 0; i < pdl_prod_senders.size(); ++i) {
			PDLSenderConfig pdl_sender = pdl_prod_senders.get(i);
			result.append ("  " + i + ":  " + pdl_sender.toString() + "\n");
		}
		result.append ("]" + "\n");

		return result.toString();
	}


	//----- Service functions -----

	// Get the currently selected list of PDL senders.
	// This returns a copy of the list, so the original cannot be modified.

	public List<PDLSenderConfig> get_pdl_senders () {
		ArrayList<PDLSenderConfig> pdl_senders = new ArrayList<PDLSenderConfig>();

		switch (pdl_enable) {

		case PDLOPT_DEV:
			for (PDLSenderConfig pdl_sender : pdl_dev_senders) {
				pdl_senders.add (pdl_sender);
			}
			break;

		case PDLOPT_PROD:
			for (PDLSenderConfig pdl_sender : pdl_prod_senders) {
				pdl_senders.add (pdl_sender);
			}
			break;

		}

		return pdl_senders;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 34001;

	private static final String M_VERSION_NAME = "ServerConfigFile";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 34000;
	protected static final int MARSHAL_ACTION_CFG = 34001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_ACTION_CFG;
	}

	// Marshal a PDL sender configuration.

	public static void marshal_pdl_sender (MarshalWriter writer, String name, PDLSenderConfig pdl_sender) {
		String host = pdl_sender.get_host();
		int port = pdl_sender.get_port();
		int connectTimeout = pdl_sender.get_connectTimeout();

		writer.marshalMapBegin (name);
		writer.marshalString ("host", host);
		writer.marshalInt ("port", port);
		writer.marshalInt ("connectTimeout", connectTimeout);
		writer.marshalMapEnd ();

		return;
	}

	// Unmarshal a PDL sender configuration.

	public static PDLSenderConfig unmarshal_pdl_sender (MarshalReader reader, String name) {

		reader.unmarshalMapBegin (name);
		String host = reader.unmarshalString ("host");
		int port = reader.unmarshalInt ("port");
		int connectTimeout = reader.unmarshalInt ("connectTimeout");
		reader.unmarshalMapEnd ();

		return new PDLSenderConfig (host, port, connectTimeout);
	}

	// Marshal a PDL sender configuration list.

	public static void marshal_pdl_sender_list (MarshalWriter writer, String name, List<PDLSenderConfig> pdl_sender_list) {
		int n = pdl_sender_list.size();
		writer.marshalArrayBegin (name, n);
		for (PDLSenderConfig pdl_sender : pdl_sender_list) {
			marshal_pdl_sender (writer, null, pdl_sender);
		}
		writer.marshalArrayEnd ();
		return;
	}

	// Unmarshal a PDL sender configuration list.

	public static ArrayList<PDLSenderConfig> unmarshal_pdl_sender_list (MarshalReader reader, String name) {
		ArrayList<PDLSenderConfig> pdl_sender_list = new ArrayList<PDLSenderConfig>();
		int n = reader.unmarshalArrayBegin (name);
		for (int i = 0; i < n; ++i) {
			pdl_sender_list.add (unmarshal_pdl_sender (reader, null));
		}
		reader.unmarshalArrayEnd ();
		return pdl_sender_list;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Error check

		check_invariant();

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalString    (        "db_host"          , db_host          );
		writer.marshalInt       (        "db_port"          , db_port          );
		writer.marshalString    (        "db_name"          , db_name          );
		writer.marshalString    (        "db_user"          , db_user          );
		writer.marshalString    (        "db_password"      , db_password      );
		writer.marshalString    (        "activemq_host"    , activemq_host    );
		writer.marshalInt       (        "activemq_port"    , activemq_port    );
		writer.marshalString    (        "activemq_user"    , activemq_user    );
		writer.marshalString    (        "activemq_password", activemq_password);
		writer.marshalString    (        "log_con_aafs"     , log_con_aafs     );
		writer.marshalString    (        "log_con_intake"   , log_con_intake   );
		writer.marshalString    (        "log_con_control"  , log_con_control  );
		writer.marshalString    (        "log_summary"      , log_summary      );
		writer.marshalString    (        "comcat_url"       , comcat_url       );
		writer.marshalInt       (        "pdl_enable"       , pdl_enable       );
		writer.marshalString    (        "pdl_key_filename" , pdl_key_filename );
		marshal_pdl_sender_list (writer, "pdl_dev_senders"  , pdl_dev_senders  );
		marshal_pdl_sender_list (writer, "pdl_prod_senders" , pdl_prod_senders );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		db_host           = reader.unmarshalString    (        "db_host"          );
		db_port           = reader.unmarshalInt       (        "db_port"          );
		db_name           = reader.unmarshalString    (        "db_name"          );
		db_user           = reader.unmarshalString    (        "db_user"          );
		db_password       = reader.unmarshalString    (        "db_password"      );
		activemq_host     = reader.unmarshalString    (        "activemq_host"    );
		activemq_port     = reader.unmarshalInt       (        "activemq_port"    );
		activemq_user     = reader.unmarshalString    (        "activemq_user"    );
		activemq_password = reader.unmarshalString    (        "activemq_password");
		log_con_aafs      = reader.unmarshalString    (        "log_con_aafs"     );
		log_con_intake    = reader.unmarshalString    (        "log_con_intake"   );
		log_con_control   = reader.unmarshalString    (        "log_con_control"  );
		log_summary       = reader.unmarshalString    (        "log_summary"      );
		comcat_url        = reader.unmarshalString    (        "comcat_url"       );
		pdl_enable        = reader.unmarshalInt       (        "pdl_enable"       );
		pdl_key_filename  = reader.unmarshalString    (        "pdl_key_filename" );
		pdl_dev_senders   = unmarshal_pdl_sender_list (reader, "pdl_dev_senders"  );
		pdl_prod_senders  = unmarshal_pdl_sender_list (reader, "pdl_prod_senders" );

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

	public ServerConfigFile unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, ServerConfigFile obj) {

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

	public static ServerConfigFile unmarshal_poly (MarshalReader reader, String name) {
		ServerConfigFile result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("ServerConfigFile.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_ACTION_CFG:
			result = new ServerConfigFile();
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

	public static ServerConfigFile unmarshal_config (String filename, Class<?> requester) {
		MarshalReader reader = OAFParameterSet.load_file_as_json (filename,requester);
		return (new ServerConfigFile()).unmarshal (reader, null);
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ServerConfigFile : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1
		// Unmarshal from the configuration file, and display it.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ServerConfigFile : Invalid 'test1' subcommand");
				return;
			}

			// Read the configuration file

			ServerConfigFile action_cfg = unmarshal_config ("ServerConfig.json", ServerConfig.class);

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
				System.err.println ("ServerConfigFile : Invalid 'test2' subcommand");
				return;
			}

			// Read the configuration file

			ServerConfigFile action_cfg = unmarshal_config ("ServerConfig.json", ServerConfig.class);

			// Display it

			System.out.println (action_cfg.toString());

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ServerConfigFile.marshal_poly (store, null, action_cfg);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			action_cfg = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			action_cfg = ServerConfigFile.unmarshal_poly (retrieve, null);
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
				System.err.println ("ServerConfigFile : Invalid 'test3' subcommand");
				return;
			}

			// Read the configuration file

			ServerConfigFile action_cfg = unmarshal_config ("ServerConfig.json", ServerConfig.class);

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
			action_cfg = (new ServerConfigFile()).unmarshal (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (action_cfg.toString());

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ServerConfigFile : Unrecognized subcommand : " + args[0]);
		return;

	}

}
