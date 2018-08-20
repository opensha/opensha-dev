package scratch.aftershockStatistics.aafs;

import java.util.List;

import scratch.aftershockStatistics.OAFParameterSet;

import scratch.aftershockStatistics.pdl.PDLSenderConfig;

/**
 * Configuration for AAFS server actions.
 * Author: Michael Barall 06/03/2018.
 *
 * To use, create an object of this class, and then call its methods to obtain configuration parameters.
 *
 * Parameters come from a configuration file, in the format of ServerConfigFile.
 */
public class ServerConfig {

	//----- Parameter set -----

	// Cached parameter set.

	private static ServerConfigFile cached_param_set = null;

	// Parameter set.

	private ServerConfigFile param_set;

	// Get the parameter set.

	private static synchronized ServerConfigFile get_param_set () {

		// If we have a cached parameter set, return it

		if (cached_param_set != null) {
			return cached_param_set;
		}

		// Working data

		ServerConfigFile wk_param_set = null;

		// Any error reading the parameters aborts the program

		try {

			// Read the configuation file

			wk_param_set = ServerConfigFile.unmarshal_config ("ServerConfig.json", ServerConfig.class);

		} catch (Exception e) {
			e.printStackTrace();
            System.err.println("ServerConfig: Error loading parameter file ServerConfig.json, unable to continue");
            System.exit(0);
			//throw new RuntimeException("ServerConfig: Error loading parameter file ServerConfig.json", e);
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

	public ServerConfig () {
		param_set = get_param_set ();
	}

	// Display our contents

	@Override
	public String toString() {
		return "ServerConfig:\n" + param_set.toString();
	}


	//----- Parameter access -----

	// Database host name or IP address.

    public String getDb_host() {
        return param_set.db_host;
    }

	// Database port number.

    public int getDb_port() {
        return param_set.db_port;
    }

	// Database name.  Used for both database access and user authentication.

    public String getDb_name() {
        return param_set.db_name;
    }

	// Database user name.  This name provides read/write access.

    public String getDb_user() {
        return param_set.db_user;
    }

	// Database password.

    public String getDb_password() {
        return param_set.db_password;
    }

	// ActiveMQ host name or IP address.

    public String getActivemq_host() {
        return param_set.activemq_host;
    }

	// ActiveMQ port number.

    public String getActivemq_port() {
        return String.valueOf (param_set.activemq_port);
    }

	// ActiveMQ user name.

    public String getActivemq_user() {
        return param_set.activemq_user;
    }

	// ActiveMQ password.

    public String getActivemq_password() {
        return param_set.activemq_password;
    }

	// Pattern for AAFS console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String get_log_con_aafs() {
        return param_set.log_con_aafs;
    }

	// Pattern for intake console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String get_log_con_intake() {
        return param_set.log_con_intake;
    }

	// Pattern for control console log filenames, in the format of SimpleDateFormat, or "" if none.

	public String get_log_con_control() {
        return param_set.log_con_control;
    }

	// Pattern for summary log filenames, in the format of SimpleDateFormat, or "" if none.

	public String get_log_summary() {
        return param_set.log_summary;
    }

	// Comcat URL.

    public String get_comcat_url() {
        return param_set.comcat_url;
    }

	// PDL enable option.

    public int get_pdl_enable() {
        return param_set.pdl_enable;
    }

	// PDL signing key filename, can be empty string for none.

    public String get_pdl_key_filename() {
        return param_set.pdl_key_filename;
    }

	// Get the currently selected list of PDL senders.
	// This returns a copy of the list, so the original cannot be modified.

	public List<PDLSenderConfig> get_pdl_senders() {
		return param_set.get_pdl_senders();
	}


	//----- Parameter modification -----

	// Get the internal server configuration file.
	// Note: This is provided so that the GUI can adjust parameters based on
	// command-line parameters or user input.
	// Note: Calling unload_data will revert all parameters to the values in
	// the configuration file.

	public ServerConfigFile get_server_config_file () {
		return param_set;
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ServerConfig : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1
		// Create an object, and display the parameters in the underlying file.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ServerConfig : Invalid 'test1' subcommand");
				return;
			}

			// Create a configuration object

			ServerConfig server_config = new ServerConfig();

			// Display it

			System.out.println (server_config.toString());

			return;
		}




		// Subcommand : Test #2
		// Command format:
		//  test2
		// Create an object, and display the parameters fetched through the object.

		if (args[0].equalsIgnoreCase ("test2")) {

			// Zero additional argument

			if (args.length != 1) {
				System.err.println ("ServerConfig : Invalid 'test2' subcommand");
				return;
			}

			// Create a configuration object

			ServerConfig server_config = new ServerConfig();

			// Display it

			System.out.println("db_host = " + server_config.getDb_host());
			System.out.println("db_port = " + server_config.getDb_port());
			System.out.println("db_name = " + server_config.getDb_name());
			System.out.println("db_user = " + server_config.getDb_user());
			System.out.println("db_password = " + server_config.getDb_password());
			System.out.println("activemq_host = " + server_config.getActivemq_host());
			System.out.println("activemq_port = " + server_config.getActivemq_port());
			System.out.println("activemq_user = " + server_config.getActivemq_user());
			System.out.println("activemq_password = " + server_config.getActivemq_password());
			System.out.println("log_con_aafs = " + server_config.get_log_con_aafs());
			System.out.println("log_con_intake = " + server_config.get_log_con_intake());
			System.out.println("log_con_control = " + server_config.get_log_con_control());
			System.out.println("log_summary = " + server_config.get_log_summary());
			System.out.println("comcat_url = " + server_config.get_comcat_url());
			System.out.println("pdl_enable = " + server_config.get_pdl_enable());
			System.out.println("pdl_key_filename = " + server_config.get_pdl_key_filename());

			List<PDLSenderConfig> pdl_senders = server_config.get_pdl_senders();
			System.out.println("pdl_senders = [");
			for (int i = 0; i < pdl_senders.size(); ++i) {
				PDLSenderConfig pdl_sender = pdl_senders.get(i);
				System.out.println("  " + i + ":  " + pdl_sender.toString());
			}
			System.out.println("]");

			// Adjust PDL enable to development, and display senders

			server_config.get_server_config_file().pdl_enable = ServerConfigFile.PDLOPT_DEV;
			pdl_senders = server_config.get_pdl_senders();
			System.out.println("pdl_senders (DEV) = [");
			for (int i = 0; i < pdl_senders.size(); ++i) {
				PDLSenderConfig pdl_sender = pdl_senders.get(i);
				System.out.println("  " + i + ":  " + pdl_sender.toString());
			}
			System.out.println("]");

			// Adjust PDL enable to production, and display senders

			server_config.get_server_config_file().pdl_enable = ServerConfigFile.PDLOPT_PROD;
			pdl_senders = server_config.get_pdl_senders();
			System.out.println("pdl_senders (PROD) = [");
			for (int i = 0; i < pdl_senders.size(); ++i) {
				PDLSenderConfig pdl_sender = pdl_senders.get(i);
				System.out.println("  " + i + ":  " + pdl_sender.toString());
			}
			System.out.println("]");

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("ServerConfig : Unrecognized subcommand : " + args[0]);
		return;

	}

}
