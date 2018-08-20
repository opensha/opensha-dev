package scratch.aftershockStatistics.aafs;

import java.nio.file.Files;
import java.nio.file.Paths;

import scratch.aftershockStatistics.pdl.PDLProductBuilderOaf;
import scratch.aftershockStatistics.pdl.PDLSender;
import gov.usgs.earthquake.product.Product;


/**
 * PDL command-line interface.
 * Author: Michael Barall 07/23/2018.
 */

// This class contains a static routine that interprets command-line options.
// Then it sets the PDL configuration, and/or sends a PDL product.
//
// To select the PDL destination, use one of these:
// --pdl=dryrun
// --pdl=dev
// --pdl=prod
// If the --pdl option is not specified, the default is taken from the pdl_default function parameter;
// if pdl_default is not specified then the default is taken from the server configuration file.
//
// To specify a PDL key file, use:
// --privateKey=PRIVATEKEYFILE
// If the --privateKey option is not specified, then products are sent unsigned.
//
// To delete a product, then in addition to the above you should also include all of these:
// --delete
// --code=PRODUCTCODE
// --eventsource=EVENTNETWORK
// --eventsourcecode=EVENTCODE
// The value of --code identifies the product that is to be deleted.  The value of --code is typically an event ID.
// The values of --eventsource and --eventsourcecode identify the event with which the product is associated;
// these determine which event page displays the product.
//
// If a JSON file exists on disk, then it can be sent as a product by including:
// --update=JSONFILENAME
// --code=PRODUCTCODE
// --eventsource=EVENTNETWORK
// --eventsourcecode=EVENTCODE
// The value of --code identifies the product that is to be sent.  The value of --code is typically an event ID.
// The product replaces any prior product that was sent with the same --code.
// The values of --eventsource and --eventsourcecode identify the event with which the product is associated;
// these determine which event page displays the product.
//
// The command line is considered to be "consumed" if a product or delete was sent to PDL,
// or if an error was reported to the user.  If the command line is consumed, then the caller
// should exit without taking further action.  If the command line is not consumed, then the
// caller should continue normally;  this permits the user to configure PDL access.

public class PDLCmd {


	// Parameter names

	public static final String PNAME_PDL = "--pdl";						// PDL access option
	public static final String PNAME_KEYFILE = "--privateKey";			// File containing private key
	public static final String PNAME_CODE = "--code";					// Product identifier
	public static final String PNAME_EVENT_NETWORK = "--eventsource";	// Network for the event
	public static final String PNAME_EVENT_CODE = "--eventsourcecode";	// Network's code for the event
	public static final String PNAME_UPDATE = "--update";				// Send product update
	public static final String PNAME_DELETE = "--delete";				// Delete product

	// Values for the PNAME_PDL parameter

	public static final String PVAL_DRYRUN = "dryrun";					// No PDL access (drt run)
	public static final String PVAL_DEV = "dev";						// PDL development server
	public static final String PVAL_PROD = "prod";						// PDL production server

	// String for splitting parameter into name and value

	public static final String PSPLIT = "=";




	// Execute PDL command line.
	// Parameters:
	//  args = Command line parameters.
	//  lo = Index of first argument to process.  Ignore args[0] thru args[lo-1].
	//  f_config = True to permit updating the ServerConfig without sending anything.
	//  f_send = True to permit sending to PDL.
	//  pdl_default = Default PDL enable option, as defined in ServerConfigFile.
	//                Can be PDLOPT_UNSPECIFIED to use the configured value.
	// Returns true if the command has been consumed.  This is always possible even if
	//  f_send is false (e.g., if an error was reported to the user).
	// Returns false if ServerConfig has been updated without sending anything.  This is
	//  possible only if f_config is true.

	public static boolean exec_pdl_cmd (String[] args, int lo, boolean f_config, boolean f_send, int pdl_default) {
	
		// Check arguments

		if (!( args != null
			&& (lo >= 0 && lo <= args.length)
			&& (f_config || f_send)
			&& ((pdl_default >= ServerConfigFile.PDLOPT_MIN && pdl_default <= ServerConfigFile.PDLOPT_MAX)
				|| pdl_default == ServerConfigFile.PDLOPT_UNSPECIFIED) )) {
			throw new IllegalArgumentException ("PDLCmdexec_pdl_cmd: Invalid arguments");
		}

		// Parameter values

		int pdl_enable = ServerConfigFile.PDLOPT_UNSPECIFIED;
		String keyfile = null;
		String code = null;
		String event_network = null;
		String event_code = null;
		String update = null;
		boolean delete = false;

		// Scan parameters

		for (int i = lo; i < args.length; ++i) {
		
			// Get the argument

			String arg = args[i].trim();

			// Split into name and value

			String name;
			String value;

			int splitix = arg.indexOf (PSPLIT);

			if (splitix < 0) {
				name = arg;
				value = null;
			} else {
				name = arg.substring (0, splitix).trim();
				value = arg.substring (splitix + PSPLIT.length()).trim();
			}

			// PDL access option

			if (name.equalsIgnoreCase (PNAME_PDL)) {
				if (pdl_enable != ServerConfigFile.PDLOPT_UNSPECIFIED) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				if (value.equalsIgnoreCase (PVAL_DRYRUN)) {
					pdl_enable = ServerConfigFile.PDLOPT_NONE;
				}
				else if (value.equalsIgnoreCase (PVAL_DEV)) {
					pdl_enable = ServerConfigFile.PDLOPT_DEV;
				}
				else if (value.equalsIgnoreCase (PVAL_PROD)) {
					pdl_enable = ServerConfigFile.PDLOPT_PROD;
				}
				else {
					System.out.println ("Invalid value in command-line option: " + arg);
					System.out.println ("Valid values are: " + PVAL_DRYRUN + ", " + PVAL_DEV + ", " + PVAL_PROD);
					return true;
				}
			}

			// Key file option

			else if (name.equalsIgnoreCase (PNAME_KEYFILE)) {
				if (keyfile != null) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				keyfile = value;
				if (!( Files.exists (Paths.get (keyfile)) )) {
					System.out.println ("Key file does not exist: " + keyfile);
					return true;
				}
			}

			// Product identifier option

			else if (name.equalsIgnoreCase (PNAME_CODE)) {
				if (!( f_send )) {
					System.out.println ("Unrecognized command-line option: " + arg);
					return true;
				}
				if (code != null) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				code = value;
			}

			// Event network option

			else if (name.equalsIgnoreCase (PNAME_EVENT_NETWORK)) {
				if (!( f_send )) {
					System.out.println ("Unrecognized command-line option: " + arg);
					return true;
				}
				if (event_network != null) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				event_network = value;
			}

			// Event code option

			else if (name.equalsIgnoreCase (PNAME_EVENT_CODE)) {
				if (!( f_send )) {
					System.out.println ("Unrecognized command-line option: " + arg);
					return true;
				}
				if (event_code != null) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				event_code = value;
			}

			// Update option

			else if (name.equalsIgnoreCase (PNAME_UPDATE)) {
				if (!( f_send )) {
					System.out.println ("Unrecognized command-line option: " + arg);
					return true;
				}
				if (update != null) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (delete) {
					System.out.println ("Command-line options cannot include both " + PNAME_UPDATE + " and " + PNAME_DELETE);
					return true;
				}
				if (value == null || value.isEmpty()) {
					System.out.println ("Missing value in command-line option: " + arg);
					return true;
				}
				update = value;
				if (!( Files.exists (Paths.get (update)) )) {
					System.out.println ("Product file does not exist: " + update);
					return true;
				}
			}

			// Delete option

			else if (name.equalsIgnoreCase (PNAME_DELETE)) {
				if (!( f_send )) {
					System.out.println ("Unrecognized command-line option: " + arg);
					return true;
				}
				if (delete) {
					System.out.println ("Duplicate command-line option: " + arg);
					return true;
				}
				if (update != null) {
					System.out.println ("Command-line options cannot include both " + PNAME_UPDATE + " and " + PNAME_DELETE);
					return true;
				}
				if (!( value == null || value.isEmpty() )) {
					System.out.println ("Command-line option " + PNAME_DELETE + " may not have a value: " + arg);
					return true;
				}
				delete = true;
			}

			// Unrecognized option

			else {
				System.out.println ("Unrecognized command-line option: " + arg);
				return true;
			}
		}

		// If no PDL access option, apply default

		if (pdl_enable == ServerConfigFile.PDLOPT_UNSPECIFIED) {
			pdl_enable = pdl_default;
		}

		// If config-only is not permitted, then we must have either update or delete

		if (!( f_config )) {
			if (!( delete || update != null )) {
				System.out.println ("Command-line options must include one of " + PNAME_UPDATE + " and " + PNAME_DELETE);
				return true;
			}
		}

		// Server configuration

		ServerConfig server_config = new ServerConfig();

		// If PDL access is specified (including by default), enter it into server configuration

		if (pdl_enable != ServerConfigFile.PDLOPT_UNSPECIFIED) {
			server_config.get_server_config_file().pdl_enable = pdl_enable;
		}

		// If key file is specified, enter it into server configuration

		if (keyfile != null) {
			server_config.get_server_config_file().pdl_key_filename = keyfile;
		}

		// Send update

		if (update != null) {

			// Check for required options

			if (!( code != null
				&& event_network != null
				&& event_code != null )) {
				System.out.println ("Cannot send PDL update because one or more of the following command-line options are missing:");
				System.out.println (PNAME_CODE + ", " + PNAME_EVENT_NETWORK + ", " + PNAME_EVENT_CODE);
				return true;
			}

			// Perform the send

			boolean send_ok = false;

			try {

				// Read in the file

				String jsonText = new String (Files.readAllBytes (Paths.get (update)));

				// Modification time, 0 means now

				long modifiedTime = 0L;

				// Review status, false means automatically generated

				boolean isReviewed = true;

				// Build the product

				Product product = PDLProductBuilderOaf.createProduct (code, event_network, event_code, isReviewed, jsonText, modifiedTime);

				// Sign the product

				PDLSender.signProduct (product);

				// Send the product, true means it is text

				PDLSender.sendProduct (product, true);

				// Success

				send_ok = true;
			}

			catch (Exception e) {
				System.out.println ("Exception occurred while attempting to send update to PDL.");
				e.printStackTrace();
			}

			// Inform user of result

			if (send_ok) {
				System.out.println ("PDL update was sent successfully.");
			} else {
				System.out.println ("PDL update was NOT sent successfully.");
			}
		
			return true;
		}

		// Send delete

		if (delete) {

			// Check for required options

			if (!( code != null
				&& event_network != null
				&& event_code != null )) {
				System.out.println ("Cannot send PDL delete because one or more of the following command-line options are missing:");
				System.out.println (PNAME_CODE + ", " + PNAME_EVENT_NETWORK + ", " + PNAME_EVENT_CODE);
				return true;
			}

			// Perform the send

			boolean send_ok = false;

			try {

				// Modification time, 0 means now

				long modifiedTime = 0L;

				// Review status, false means automatically generated

				boolean isReviewed = true;

				// Build the product

				Product product = PDLProductBuilderOaf.createDeletionProduct (code, event_network, event_code, isReviewed, modifiedTime);

				// Sign the product

				PDLSender.signProduct (product);

				// Send the product, true means it is text

				PDLSender.sendProduct (product, true);

				// Success

				send_ok = true;
			}

			catch (Exception e) {
				System.out.println ("Exception occurred while attempting to send delete to PDL.");
				e.printStackTrace();
			}

			// Inform user of result

			if (send_ok) {
				System.out.println ("PDL delete was sent successfully.");
			} else {
				System.out.println ("PDL delete was NOT sent successfully.");
			}
		
			return true;
		}

		// We only adjusted the configuration

		return false;
	}




	// Entry point.
	
	public static void main(String[] args) {

		// Execute PDL command

		int lo = 0;
		boolean f_config = false;
		boolean f_send = true;
		int pdl_default = ServerConfigFile.PDLOPT_UNSPECIFIED;

		exec_pdl_cmd (args, lo, f_config, f_send, pdl_default);
		return;
	}
}
