package scratch.aftershockStatistics.aafs;

import java.util.Properties;

import scratch.aftershockStatistics.OAFParameterSet;

/**
 * Author: Michael Barall.
 * Imported from aftershockStatistics.cmu on 03/15/2018.
 */
public class ServerConfig {

	// Our property table.

	private static Properties prop_table = null;

	// Values of specific properties.

    private static String db_host = null;
    private static int db_port = 0;
    private static String db_name = null;
    private static String db_user = null;
    private static String db_password = null;
    private static String activemq_host = null;
    private static String activemq_port = null;
    private static String activemq_user = null;
    private static String activemq_password = null;

	// List of required properties.

	private static final String[] required_props = {
		"db_host",
		"db_port",
		"db_name",
		"db_user",
		"db_password",
		"activemq_host",
		"activemq_port",
		"activemq_user",
		"activemq_password"
	};

	// Load the property table.

	private static synchronized void load_prop() {

		// If properties already loaded, do nothing

		if (prop_table != null) {
			return;
		}

		// Working data

		Properties wk_prop_table = null;
		int wk_db_port = 0;

		// Any error reading the properties aborts the program

		try {

			// Read the property file

			wk_prop_table = OAFParameterSet.load_properties ("ServerConfig.txt", ServerConfig.class, null, required_props);

			// Convert port number

			try {
				wk_db_port = Integer.parseInt(wk_prop_table.getProperty("db_port"));
			} catch (NumberFormatException e) {
				throw new RuntimeException("ServerConfig: Malformed port number", e);
			}

		} catch (Exception e) {
			e.printStackTrace();
            System.err.println("ServerConfig: Error loading property file ServerConfig.txt, unable to continue");
            System.exit(0);
			//throw new RuntimeException("ServerConfig: Unable to load property file ServerConfig.txt", e);
		}

		// Save the properties

        db_host = wk_prop_table.getProperty("db_host");
        db_port = wk_db_port;
        db_name = wk_prop_table.getProperty("db_name");
        db_user = wk_prop_table.getProperty("db_user");
        db_password = wk_prop_table.getProperty("db_password");
        activemq_host = wk_prop_table.getProperty("activemq_host");
        activemq_port = wk_prop_table.getProperty("activemq_port");
        activemq_user = wk_prop_table.getProperty("activemq_user");
        activemq_password = wk_prop_table.getProperty("activemq_password");

		prop_table = wk_prop_table;
		return;
	}

	// Constructor loads the property table if needed.

    public ServerConfig() {
		load_prop();
    }

	// Property getter functions.

    public String getDb_host() {
        return db_host;
    }

    public int getDb_port() {
        return db_port;
    }

    public String getDb_name() {
        return db_name;
    }

    public String getDb_user() {
        return db_user;
    }

    public String getDb_password() {
        return db_password;
    }

    public String getActivemq_host() {
        return activemq_host;
    }

    public String getActivemq_port() {
        return activemq_port;
    }

    public String getActivemq_user() {
        return activemq_user;
    }

    public String getActivemq_password() {
        return activemq_password;
    }

	// Simple test loads the properties and prints them out.
	
	public static void main(String[] args) {
		ServerConfig config = new ServerConfig();

		System.out.println("db_host = " + config.getDb_host());
		System.out.println("db_port = " + config.getDb_port());
		System.out.println("db_name = " + config.getDb_name());
		System.out.println("db_user = " + config.getDb_user());
		System.out.println("db_password = " + config.getDb_password());
		System.out.println("activemq_host = " + config.getActivemq_host());
		System.out.println("activemq_port = " + config.getActivemq_port());
		System.out.println("activemq_user = " + config.getActivemq_user());
		System.out.println("activemq_password = " + config.getActivemq_password());

		return;
	}
}
