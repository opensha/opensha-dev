package scratch.aftershockStatistics.aafs;

import java.util.List;

import java.io.Closeable;

import scratch.aftershockStatistics.aafs.MongoDBUtil;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.AftershockStatsCalc;
import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.RJ_AftershockModel_SequenceSpecific;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.TimeSplitOutputStream;
import scratch.aftershockStatistics.util.ConsoleRedirector;


/**
 * AAFS server command-line interface.
 * Author: Michael Barall 05/23/2018.
 */
public class ServerCmd {




	// cmd_start - Run the server.

	public static void cmd_start(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerCmd : Invalid 'start' subcommand");
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_aafs(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

			// Summary log

			TimeSplitOutputStream sum_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_summary(), start_time);

		){

			try {

				// Say hello
			
				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);
				System.out.println (VersionInfo.get_title());
				System.out.println ("");
				System.out.println ("AAFS server is starting at " + SimpleUtils.time_to_string (start_time));

				// Get a task dispatcher

				TaskDispatcher dispatcher = new TaskDispatcher();

				// Install the log files

				dispatcher.set_console_log_tsop (con_tsop);
				dispatcher.set_summary_log_tsop (sum_tsop);

				// Run it

				dispatcher.run();

				// Display final status

				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);

				int dispatcher_state = dispatcher.get_dispatcher_state();
				if (dispatcher_state == TaskDispatcher.STATE_SHUTDOWN) {
					System.out.println ("AAFS server exited normally at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				} else {
					System.out.println ("Server exited abnormally at " + SimpleUtils.time_to_string (ServerClock.get_true_time()) + ", final state code = " + dispatcher_state);
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("AAFS server had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("AAFS server had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_stop - Stop the server.

	public static void cmd_stop(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerCmd : Invalid 'stop' subcommand");
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_control(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

		){

			try {

				// Say hello
			
				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);
				System.out.println ("Sending shutdown command to AAFS server at " + SimpleUtils.time_to_string (start_time));

				// Post the shutdown task

				boolean result = TaskDispatcher.post_shutdown ("ServerCmd");

				// Display result

				if (result) {
					System.out.println ("Shutdown command was sent to AAFS server at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
					System.out.println ("It takes about 30 seconds for the shutdown to be complete.");
				} else {
					System.out.println ("Unable to send shutdown command to AAFS server at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("Shutdown command had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("Shutdown command had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_pdl_intake - Intake an event from PDL.

	public static void cmd_pdl_intake(String[] args) {

		// At least one additional argument

		if (args.length < 2) {
			System.err.println ("ServerCmd : Invalid 'intake' subcommand");
			return;
		}

		OpIntakePDL payload = new OpIntakePDL();

		try {
			payload.setup (args, 1, args.length);
		}

		// If PDL arguments do not parse, just drop it silently
		
		catch (Exception e) {
			return;
		}

		// If no event id, just drop it silently

		if (!( payload.has_event_id() )) {
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_intake(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

		){

			try {

				// Get the timestamp

				String timestamp = SimpleUtils.time_to_string_no_z (start_time);

				// If we don't have location, magnitude, and time, drop it

				if (!( payload.has_lat_lon_depth_mag() && payload.has_event_time() )) {
					System.out.println (timestamp + " DROP " + payload.event_id + " (incomplete data)");
					return;
				}

				// Build event info

				String event_info = SimpleUtils.event_id_and_info_one_line (
										payload.event_id,
										payload.mainshock_time,
										payload.mainshock_mag,
										payload.mainshock_lat,
										payload.mainshock_lon,
										payload.mainshock_depth);

				// Test if event passes the intake filter, using the minimum magnitude criterion

				ActionConfig action_config = new ActionConfig();

				long the_time = ServerClock.get_time();

				long sched_time = the_time;

				IntakeSphRegion intake_region = action_config.get_pdl_intake_region_for_min_mag (
						payload.mainshock_lat, payload.mainshock_lon, payload.mainshock_mag);

				if (intake_region == null) {

					// Didn't pass, check using the intake magnitude criterion

					intake_region = action_config.get_pdl_intake_region_for_intake_mag (
							payload.mainshock_lat, payload.mainshock_lon, payload.mainshock_mag);

					// Schedule task for projected time of first forecast, ignoring origin skew

					sched_time = payload.mainshock_time
									+ action_config.get_next_forecast_lag(0L)
									+ action_config.get_comcat_clock_skew();

					// If we didn't pass, or if the scheduled time has already passed, then drop event

					if (intake_region == null || sched_time < the_time) {
						System.out.println (timestamp + " INTAKE-FILTER " + event_info);
						return;
					}
				}

				// Test if event is in acceptable age range

				long the_min_time = start_time - action_config.get_pdl_intake_max_age();
				long the_max_time = start_time + action_config.get_pdl_intake_max_future();

				if (!( payload.mainshock_time >= the_min_time && payload.mainshock_time <= the_max_time )) {
				
					// Out of acceptable age range, drop event

					System.out.println (timestamp + " AGE-RANGE " + event_info);
					return;
				}

				// Post the task

				String event_id = payload.event_id;

				int opcode = TaskDispatcher.OPCODE_INTAKE_PDL;
				int stage = 0;

				boolean result = TaskDispatcher.post_task (event_id, sched_time, the_time, "ServerCmd", opcode, stage, payload.marshal_task());

				// Display result

				if (result) {
					System.out.println (timestamp + " ACCEPT " + event_info);
				} else {
					System.out.println (timestamp + " SUBMISSION-FAILURE " + event_info);
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("PDL intake had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("PDL intake had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_add_event - Tell the server to watch an event.

	public static void cmd_add_event(String[] args) {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerCmd : Invalid 'cmd_add_event' subcommand");
			return;
		}

		String event_id = args[1];

		OpIntakeSync payload = new OpIntakeSync();
		payload.setup ();

		// Say hello

		System.out.println ("Add event, event_id = " + event_id);

		// Post the task

		int opcode = TaskDispatcher.OPCODE_INTAKE_SYNC;
		int stage = 0;

		long the_time = ServerClock.get_time();

		boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerCmd", opcode, stage, payload.marshal_task());

		// Display result

		if (result) {
			System.out.println ("Event was sent to AAFS server.");
		} else {
			System.out.println ("Unable to send event to AAFS server.");
		}

		return;
	}




	// cmd_start_with_pdl - Run the server, and accept PDL options on the command line.

	public static void cmd_start_with_pdl(String[] args) {

		// Any number of additional arguments

		if (args.length < 1) {
			System.err.println ("ServerCmd : Invalid 'start_with_pdl' subcommand");
			return;
		}

		// Read PDL options

		int lo = 1;
		boolean f_config = true;
		boolean f_send = false;
		int pdl_default = ServerConfigFile.PDLOPT_UNSPECIFIED;

		if (PDLCmd.exec_pdl_cmd (args, lo, f_config, f_send, pdl_default)) {
			System.out.println ("AAFS server not started due to error in PDL options.");
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_aafs(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

			// Summary log

			TimeSplitOutputStream sum_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_summary(), start_time);

		){

			try {

				// Say hello
			
				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);
				System.out.println (VersionInfo.get_title());
				System.out.println ("");
				System.out.println ("AAFS server is starting at " + SimpleUtils.time_to_string (start_time));

				// Get a task dispatcher

				TaskDispatcher dispatcher = new TaskDispatcher();

				// Install the log files

				dispatcher.set_console_log_tsop (con_tsop);
				dispatcher.set_summary_log_tsop (sum_tsop);

				// Run it

				dispatcher.run();

				// Display final status

				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);

				int dispatcher_state = dispatcher.get_dispatcher_state();
				if (dispatcher_state == TaskDispatcher.STATE_SHUTDOWN) {
					System.out.println ("AAFS server exited normally at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				} else {
					System.out.println ("Server exited abnormally at " + SimpleUtils.time_to_string (ServerClock.get_true_time()) + ", final state code = " + dispatcher_state);
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("AAFS server had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("AAFS server had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_start_comcat_poll - Tell the server to start polling Comcat.

	public static void cmd_start_comcat_poll(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerCmd : Invalid 'start_comcat_poll' subcommand");
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_control(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

		){

			try {

				String event_id = ServerComponent.EVID_POLL;

				OpPollComcatStart payload = new OpPollComcatStart();
				payload.setup ();

				// Say hello

				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);
				System.out.println ("Sending command to start polling Comcat at " + SimpleUtils.time_to_string (start_time));

				// Post the task

				int opcode = TaskDispatcher.OPCODE_POLL_COMCAT_START;
				int stage = 0;

				long the_time = ServerClock.get_time();

				boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerCmd", opcode, stage, payload.marshal_task());

				// Display result

				if (result) {
					System.out.println ("Command to start polling Comcat was sent to AAFS server at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
					System.out.println ("It takes about 30 seconds for the command to take effect.");
				} else {
					System.out.println ("Unable to send AAFS server command to start polling Comcat at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("Command to start polling had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("Command to start polling had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_stop_comcat_poll - Tell the server to stop polling Comcat.

	public static void cmd_stop_comcat_poll(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerCmd : Invalid 'stop_comcat_poll' subcommand");
			return;
		}

		// Get current time

		long start_time = ServerClock.get_true_time();

		try (

			// Console redirection and log

			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (
										(new ServerConfig()).get_log_con_control(), start_time);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));

		){

			try {

				String event_id = ServerComponent.EVID_POLL;

				OpPollComcatStop payload = new OpPollComcatStop();
				payload.setup ();

				// Say hello

				System.out.println (ServerComponent.LOG_SEPARATOR_LINE);
				System.out.println ("Sending command to stop polling Comcat at " + SimpleUtils.time_to_string (start_time));

				// Post the task

				int opcode = TaskDispatcher.OPCODE_POLL_COMCAT_STOP;
				int stage = 0;

				long the_time = ServerClock.get_time();

				boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerCmd", opcode, stage, payload.marshal_task());

				// Display result

				if (result) {
					System.out.println ("Command to stop polling Comcat was sent to AAFS server at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
					System.out.println ("It takes about 30 seconds for the command to take effect.");
				} else {
					System.out.println ("Unable to send AAFS server command to stop polling Comcat at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				}

			}

			// Report any uncaught exceptions

			catch (Exception e) {
				System.out.println ("Command to stop polling had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
				e.printStackTrace();
			}
		}

		// Report any uncaught exceptions

		catch (Exception e) {
			System.out.println ("Command to stop polling had an uncaught exception at " + SimpleUtils.time_to_string (ServerClock.get_true_time()));
            e.printStackTrace();
		}

		return;
	}




	// cmd_show_version - Print version info to the standard output.

	public static void cmd_show_version(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerCmd : Invalid 'show_version' subcommand");
			return;
		}

		// Say hello

		System.out.println (VersionInfo.get_title());

		return;
	}




	// Entry point.
	
	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ServerCmd : Missing subcommand");
			return;
		}

		switch (args[0].toLowerCase()) {

		// Subcommand : start
		// Command format:
		//  start
		// Run the server.

		case "start":
			try {
				cmd_start(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : stop
		// Command format:
		//  stop
		// Stop the server.

		case "stop":
			try {
				cmd_stop(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : pdl_intake
		// Command format:
		//  pdl_intake  arg...
		// Post a PDL intake command for the given command line.
		//
		// This command is normally invoked from a bash script, which in turn is
		// invoked from a PDL indexer listener.  It expects to receive the command
		// line parameters from the indexer listener.  The command in the script
		// file should look like this:
		//
		// java -cp jar-file-list scratch.aftershockStatistics.aafs.ServerCmd pdl_intake "$@"
		//
		// Notice that $@ appears inside quotes.  The quotes are necessary for bash
		// to correctly pass the script parameters to Java.

		case "pdl_intake":
			try {
				cmd_pdl_intake(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : add_event
		// Command format:
		//  add_event  event_id
		// Post a sync intake command for the given event.

		case "add_event":
			try {
				cmd_add_event(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : start_with_pdl
		// Command format:
		//  start_with_pdl  [pdl_option...]
		// Run the server, after parsing PDL options.

		case "start_with_pdl":
			try {
				cmd_start_with_pdl(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : start_comcat_poll
		// Command format:
		//  start_comcat_poll 
		// Post a command to start polling Comcat.

		case "start_comcat_poll":
			try {
				cmd_start_comcat_poll(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : stop_comcat_poll
		// Command format:
		//  stop_comcat_poll
		// Post a command to stop polling Comcat.

		case "stop_comcat_poll":
			try {
				cmd_stop_comcat_poll(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		// Subcommand : show_version
		// Command format:
		//  show_version
		// Print version info to the standard output.

		case "show_version":
			try {
				cmd_show_version(args);
            } catch (Exception e) {
                e.printStackTrace();
			}
			return;

		}

		// Unrecognized subcommand.

		System.err.println ("ServerCmd : Unrecognized subcommand : " + args[0]);
		return;
	}
}
