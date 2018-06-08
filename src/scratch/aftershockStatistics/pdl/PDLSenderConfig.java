package scratch.aftershockStatistics.pdl;


/**
 * Configuration settings for a PDL sender.
 * Author: Michael Barall 06/03/2018.
 */
public class PDLSenderConfig {

	// The host name.

	private String host;

	// The port number.

	private int port;

	// The connection timeout, in milliseconds (default value is 15000 or 15 seconds).
	// We allow values from 1 second to 300 seconds.

	private int connectTimeout;

	// Minimum and maximum port numbers.

	public static final int MIN_PORT = 1024;
	public static final int MAX_PORT = 65535;

	// Minimum and maximum connection timeouts, and default.

	public static final int MIN_CONNECT_TIMEOUT = 1000;
	public static final int MAX_CONNECT_TIMEOUT = 300000;
	public static final int DEF_CONNECT_TIMEOUT = 15000;

	// Constructor.

	public PDLSenderConfig (String host, int port, int connectTimeout) {
		if (!( host != null
			&& host.trim().length() > 0
			&& port >= MIN_PORT
			&& port <= MAX_PORT
			&& connectTimeout >= MIN_CONNECT_TIMEOUT
			&& connectTimeout <= MAX_CONNECT_TIMEOUT )) {

			throw new IllegalArgumentException (
				"PDLSenderConfig: Invalid PDL sender configuration: host = "
				+ ((host == null) ? "<null>" : host) + ", port = " + port + ", connectTimeout = " + connectTimeout);
		}
	
		this.host = host.trim();
		this.port = port;
		this.connectTimeout = connectTimeout;
	}

	// Getters.

	public String get_host () {
		return host;
	}

	public int get_port () {
		return port;
	}

	public int get_connectTimeout () {
		return connectTimeout;
	}

	// Display our contents.

	@Override
	public String toString() {
		return "host = " + ((host == null) ? "<null>" : host) + ", port = " + port + ", connectTimeout = " + connectTimeout; 
	}

}
