package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.Arrays;

import java.io.Reader;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;

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

import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;

import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.TimeSplitOutputStream;
import scratch.aftershockStatistics.util.ConsoleRedirector;

import gov.usgs.earthquake.product.Product;
import scratch.aftershockStatistics.pdl.PDLProductBuilderOaf;
import scratch.aftershockStatistics.pdl.PDLSender;


/**
 * Holds a set of tests for the AAFS server code.
 * Author: Michael Barall 03/16/2018.
 */
public class ServerTest {




	// Test #1 - Add a few elements to the task pending queue.

	public static void test1(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test1' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
		
			String event_id;
			long sched_time;
			long submit_time;
			String submit_id;
			int opcode;
			int stage;
			MarshalWriter details;
		
			event_id = "Event_2";
			sched_time = 20100L;
			submit_time = 20000L;
			submit_id = "Submitter_2";
			opcode = 102;
			stage = 2;
			details = PendingTask.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_2");
			details.marshalLong (null, 21010L);
			details.marshalLong (null, 21020L);
			details.marshalDouble (null, 21030.0);
			details.marshalDouble (null, 21040.0);
			details.marshalArrayEnd ();
			PendingTask.submit_task (event_id, sched_time, submit_time,
				submit_id, opcode, stage, details);
		
			event_id = "Event_4";
			sched_time = 40100L;
			submit_time = 40000L;
			submit_id = "Submitter_4_no_details";
			opcode = 104;
			stage = 4;
			details = null;
			PendingTask.submit_task (event_id, sched_time, submit_time,
				submit_id, opcode, stage, details);
		
			event_id = "Event_1";
			sched_time = 10100L;
			submit_time = 10000L;
			submit_id = "Submitter_1";
			opcode = 101;
			stage = 1;
			details = PendingTask.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_1");
			details.marshalLong (null, 11010L);
			details.marshalLong (null, 11020L);
			details.marshalDouble (null, 11030.0);
			details.marshalDouble (null, 11040.0);
			details.marshalArrayEnd ();
			PendingTask.submit_task (event_id, sched_time, submit_time,
				submit_id, opcode, stage, details);
		
			event_id = "Event_5";
			sched_time = 50100L;
			submit_time = 50000L;
			submit_id = "Submitter_5";
			opcode = 105;
			stage = 5;
			details = PendingTask.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_5");
			details.marshalLong (null, 51010L);
			details.marshalLong (null, 51020L);
			details.marshalDouble (null, 51030.0);
			details.marshalDouble (null, 51040.0);
			details.marshalArrayEnd ();
			PendingTask.submit_task (event_id, sched_time, submit_time,
				submit_id, opcode, stage, details);
		
			event_id = "Event_3";
			sched_time = 30100L;
			submit_time = 30000L;
			submit_id = "Submitter_3";
			opcode = 103;
			stage = 3;
			details = PendingTask.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_3");
			details.marshalLong (null, 31010L);
			details.marshalLong (null, 31020L);
			details.marshalDouble (null, 31030.0);
			details.marshalDouble (null, 31040.0);
			details.marshalArrayEnd ();
			PendingTask.submit_task (event_id, sched_time, submit_time,
				submit_id, opcode, stage, details);

		}

		return;
	}




	// Test #2 - Display the pending task queue, unsorted.

	public static void test2(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test2' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of pending tasks

			List<PendingTask> tasks = PendingTask.get_all_tasks_unsorted();

			// Display them

			for (PendingTask task : tasks) {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #3 - Display the pending task queue, sorted by execution time.

	public static void test3(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test3' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of pending tasks

			List<PendingTask> tasks = PendingTask.get_all_tasks();

			// Display them

			for (PendingTask task : tasks) {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #4 - Display the first task in the pending task queue, according to execution time.

	public static void test4(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test4' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the task

			PendingTask task = PendingTask.get_first_task();

			// Display it

			if (task == null) {
				System.out.println ("null");
			} else {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #5 - Display the first task in the pending task queue, before cutoff time, according to execution time.

	public static void test5(String[] args) {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test5' subcommand");
			return;
		}

		long cutoff_time = Long.parseLong(args[1]);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the task

			PendingTask task = PendingTask.get_first_ready_task (cutoff_time);

			// Display it

			if (task == null) {
				System.out.println ("null");
			} else {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #6 - Activate the first document before the cutoff time, and display the retrieved document.

	public static void test6(String[] args) {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test6' subcommand");
			return;
		}

		long cutoff_time = Long.parseLong(args[1]);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the task

			PendingTask task = PendingTask.activate_first_ready_task (cutoff_time);

			// Display it

			if (task == null) {
				System.out.println ("null");
			} else {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #7 - Activate the first document before the cutoff time, and stage it.

	public static void test7(String[] args) {

		// Three or four additional arguments

		if (args.length != 4 && args.length != 5) {
			System.err.println ("ServerTest : Invalid 'test7' subcommand");
			return;
		}

		long cutoff_time = Long.parseLong(args[1]);
		long exec_time = Long.parseLong(args[2]);
		int stage = Integer.parseInt(args[3]);
		String event_id = null;
		if (args.length >= 5) {
			event_id = args[4];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the task

			PendingTask task = PendingTask.activate_first_ready_task (cutoff_time);

			// Stage it

			if (task != null) {
				PendingTask.stage_task (task, exec_time, stage, event_id);
			}

			// Display it

			if (task == null) {
				System.out.println ("null");
			} else {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #8 - Activate the first document before the cutoff time, and delete it.

	public static void test8(String[] args) {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test8' subcommand");
			return;
		}

		long cutoff_time = Long.parseLong(args[1]);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the task

			PendingTask task = PendingTask.activate_first_ready_task (cutoff_time);

			// Delete it

			if (task != null) {
				PendingTask.delete_task (task);
			}

			// Display it

			if (task == null) {
				System.out.println ("null");
			} else {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #9 - Run task dispatcher.

	public static void test9(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test9' subcommand");
			return;
		}

		// Get a task dispatcher

		TaskDispatcher dispatcher = new TaskDispatcher();

		// Run it

		dispatcher.run();

		// Display final status

		System.out.println ("Dispatcher final state: " + dispatcher.get_dispatcher_state());

		return;
	}




	// Test #10 - Post a shutdown task.

	public static void test10(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test10' subcommand");
			return;
		}

		// Get a task dispatcher

		TaskDispatcher dispatcher = new TaskDispatcher();

		// Post the shutdown task

		boolean result = dispatcher.post_shutdown ("ServerTest");

		// Display result

		System.out.println ("Post shutdown result: " + result);

		return;
	}




	// Test #11 - Scan the pending task queue, sorted, and write a log entry for each.

	public static void test11(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test11' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of pending tasks

			List<PendingTask> tasks = PendingTask.get_all_tasks();

			// Write the log entries

			for (PendingTask task : tasks) {
				LogEntry.submit_log_entry (task, task.get_sched_time() + 100L, task.get_opcode() + 1000, "Result_for_" + task.get_opcode());
			}

		}

		return;
	}




	// Test #12 - Scan the pending task queue, sorted, and search the log for each.

	public static void test12(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test12' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of pending tasks

			List<PendingTask> tasks = PendingTask.get_all_tasks();

			// Search for the log entries, display task and matching log entry

			for (PendingTask task : tasks) {
				System.out.println (task.toString());
				LogEntry entry = LogEntry.get_log_entry_for_key (task.get_record_key());
				if (entry == null) {
					System.out.println ("LogEntry: null");
				} else {
					System.out.println (entry.toString());
				}
			}

		}

		return;
	}




	// Test #13 - Search the log for log time and/or event id.

	public static void test13(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test13' subcommand");
			return;
		}

		long log_time_lo = Long.parseLong(args[1]);
		long log_time_hi = Long.parseLong(args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching log entries

			List<LogEntry> entries = LogEntry.get_log_entry_range (log_time_lo, log_time_hi, event_id);

			// Display them

			for (LogEntry entry : entries) {
				System.out.println (entry.toString());
			}

		}

		return;
	}




	// Test #14 - Search the log for log time and/or event id, and delete the entries.

	public static void test14(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test14' subcommand");
			return;
		}

		long log_time_lo = Long.parseLong(args[1]);
		long log_time_hi = Long.parseLong(args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching log entries

			List<LogEntry> entries = LogEntry.get_log_entry_range (log_time_lo, log_time_hi, event_id);

			// Display them, and delete

			for (LogEntry entry : entries) {
				System.out.println (entry.toString());
				LogEntry.delete_log_entry (entry);
			}

		}

		return;
	}




	// Test #15 - Post a task with given event id, opcode, stage, and details.

	public static void test15(String[] args) {

		// Three or four additional arguments

		if (args.length != 4 && args.length != 5) {
			System.err.println ("ServerTest : Invalid 'test15' subcommand");
			return;
		}

		String event_id = args[1];
		if (event_id.equalsIgnoreCase ("-")) {
			event_id = "";
		}
		int opcode = Integer.parseInt(args[2]);
		int stage = Integer.parseInt(args[3]);
		MarshalWriter details = PendingTask.begin_details();
		if (args.length == 5) {
			details.marshalMapBegin (null);
			details.marshalString ("value", args[4]);
			details.marshalMapEnd ();
		}

		// Post the task

		long the_time = ServerClock.get_time();

		TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, details);

		return;
	}




	// Test #16 - Execute the next task.

	public static void test16(String[] args) {

		// Zero or one additional arguments

		if (args.length < 1 || args.length > 2) {
			System.err.println ("ServerTest : Invalid 'test16' subcommand");
			return;
		}

		boolean f_adjust_time = false;
		if (args.length >= 2) {
			f_adjust_time = Boolean.parseBoolean (args[1]);
		}

		// Get a task dispatcher

		TaskDispatcher dispatcher = new TaskDispatcher();

		// Run one task

		boolean f_verbose = true;

		dispatcher.run_next_task (f_verbose, f_adjust_time);

		return;
	}




	// Test #17 - Write a catalog snapshot.

	public static void test17(String[] args) {

		// Three additional arguments

		if (args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test17' subcommand");
			return;
		}

		double start_time_days = Double.parseDouble (args[1]);
		double end_time_days   = Double.parseDouble (args[2]);
		String event_id = args[3];

		long start_time = Math.round(start_time_days * 86400000L);
		long end_time   = Math.round(end_time_days   * 86400000L);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Create a simulated aftershock sequence, using the method of RJ_AftershockModel_SequenceSpecific
			// (Start and end times should be in the range of 0 to 30 days)
			
			double a = -1.67;
			double b = 0.91;
			double c = 0.05;
			double p = 1.08;
			double magMain = 7.5;
			double magCat = 2.5;
			double capG = 1.25;
			double capH = 0.75;

			ObsEqkRupList aftershockList = AftershockStatsCalc.simAftershockSequence(a, b, magMain, magCat, capG, capH, p, c, start_time_days, end_time_days);

			CompactEqkRupList rupture_list = new CompactEqkRupList (aftershockList);

			// Write the rupture sequence

			CatalogSnapshot entry_in = CatalogSnapshot.submit_catalog_shapshot (null, event_id, start_time, end_time, rupture_list);

			System.out.println (entry_in.toString());

			// Search for it

			CatalogSnapshot entry_out = CatalogSnapshot.get_catalog_shapshot_for_key (entry_in.get_record_key());

			System.out.println (entry_out.toString());

			// Use the retrieved rupture sequence to make a sequence-specific model
		
			double min_a = -2.0;
			double max_a = -1.0;
			int num_a = 101;

			double min_p = 0.9; 
			double max_p = 1.2; 
			int num_p = 31;
		
			double min_c=0.05;
			double max_c=0.05;
			int num_c=1;

			ObsEqkRupture mainShock = new ObsEqkRupture("0", 0L, null, magMain);

			CompactEqkRupList rupture_list_out = entry_out.get_rupture_list();

			// Make the model, it will output some information

			RJ_AftershockModel_SequenceSpecific seqModel =
				new RJ_AftershockModel_SequenceSpecific(mainShock, rupture_list_out,
			 								magCat, capG, capH,
											b, start_time_days, end_time_days,
											min_a, max_a, num_a, 
											min_p, max_p, num_p, 
											min_c, max_c, num_c);

		}

		return;
	}




	// Test #18 - Search the catalog snapshots for end time and/or event id.

	public static void test18(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test18' subcommand");
			return;
		}

		double end_time_lo_days = Double.parseDouble (args[1]);
		double end_time_hi_days = Double.parseDouble (args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		long end_time_lo = Math.round(end_time_lo_days * 86400000L);
		long end_time_hi = Math.round(end_time_hi_days * 86400000L);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching catalog snapshots

			List<CatalogSnapshot> entries = CatalogSnapshot.get_catalog_snapshot_range (end_time_lo, end_time_hi, event_id);

			// Display them

			for (CatalogSnapshot entry : entries) {
				System.out.println (entry.toString());
			}

		}

		return;
	}




	// Test #19 - Search the catalog snapshots for end time and/or event id, and delete the matching entries.

	public static void test19(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test19' subcommand");
			return;
		}

		double end_time_lo_days = Double.parseDouble (args[1]);
		double end_time_hi_days = Double.parseDouble (args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		long end_time_lo = Math.round(end_time_lo_days * 86400000L);
		long end_time_hi = Math.round(end_time_hi_days * 86400000L);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching catalog snapshots

			List<CatalogSnapshot> entries = CatalogSnapshot.get_catalog_snapshot_range (end_time_lo, end_time_hi, event_id);

			// Display them, and delete

			for (CatalogSnapshot entry : entries) {
				System.out.println (entry.toString());
				CatalogSnapshot.delete_catalog_snapshot (entry);
			}

		}

		return;
	}




	// Test #20 - Add a few elements to the timeline.

	public static void test20(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test20' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
		
			String event_id;
			String[] comcat_ids;
			long action_time;
			int actcode;
			MarshalWriter details;
		
			event_id = "Event_2";
			comcat_ids = new String[]{"ccid_21", "ccid_22", "ccid_23"};
			action_time = 20102L;
			actcode = 102;
			details = TimelineEntry.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_2");
			details.marshalLong (null, 21010L);
			details.marshalLong (null, 21020L);
			details.marshalDouble (null, 21030.0);
			details.marshalDouble (null, 21040.0);
			details.marshalArrayEnd ();
			TimelineEntry.submit_timeline_entry (null, action_time, event_id,
				comcat_ids, actcode, details);
		
			event_id = "Event_4";
			comcat_ids = new String[]{"ccid_41", "ccid_42", "ccid_43", "ccid_23"};
			action_time = 40104L;
			actcode = 104;
			details = null;
			TimelineEntry.submit_timeline_entry (null, action_time, event_id,
				comcat_ids, actcode, details);
		
			event_id = "Event_1";
			comcat_ids = new String[]{"ccid_11", "ccid_12", "ccid_13"};
			action_time = 10101L;
			actcode = 101;
			details = TimelineEntry.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_1");
			details.marshalLong (null, 11010L);
			details.marshalLong (null, 11020L);
			details.marshalDouble (null, 11030.0);
			details.marshalDouble (null, 11040.0);
			details.marshalArrayEnd ();
			TimelineEntry.submit_timeline_entry (null, action_time, event_id,
				comcat_ids, actcode, details);
		
			event_id = "Event_5";
			comcat_ids = new String[]{"ccid_51", "ccid_52", "ccid_53", "ccid_23", "ccid_13"};
			action_time = 50105L;
			actcode = 105;
			details = TimelineEntry.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_5");
			details.marshalLong (null, 51010L);
			details.marshalLong (null, 51020L);
			details.marshalDouble (null, 51030.0);
			details.marshalDouble (null, 51040.0);
			details.marshalArrayEnd ();
			TimelineEntry.submit_timeline_entry (null, action_time, event_id,
				comcat_ids, actcode, details);
		
			event_id = "Event_3";
			comcat_ids = new String[]{"ccid_31", "ccid_32"};
			action_time = 30103L;
			actcode = 103;
			details = TimelineEntry.begin_details();
			details.marshalArrayBegin (null, 5);
			details.marshalString (null, "Details_3");
			details.marshalLong (null, 31010L);
			details.marshalLong (null, 31020L);
			details.marshalDouble (null, 31030.0);
			details.marshalDouble (null, 31040.0);
			details.marshalArrayEnd ();
			TimelineEntry.submit_timeline_entry (null, action_time, event_id,
				comcat_ids, actcode, details);

		}

		return;
	}




	// Test #21 - Search the timeline for action time and/or event id; using list.

	public static void test21(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test21' subcommand");
			return;
		}

		long action_time_lo = Long.parseLong(args[1]);
		long action_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] action_time_div_rem = null;
		if (div_rem > 0L) {
			action_time_div_rem = new long[2];
			action_time_div_rem[0] = div_rem / 1000L;
			action_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching timeline entries

			List<TimelineEntry> entries = TimelineEntry.get_timeline_entry_range (action_time_lo, action_time_hi, event_id, comcat_ids, action_time_div_rem);

			// Display them

			for (TimelineEntry entry : entries) {
				System.out.println (entry.toString());
			}

		}

		return;
	}




	// Test #22 - Search the timeline for action time and/or event id; using iterator.

	public static void test22(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test22' subcommand");
			return;
		}

		long action_time_lo = Long.parseLong(args[1]);
		long action_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] action_time_div_rem = null;
		if (div_rem > 0L) {
			action_time_div_rem = new long[2];
			action_time_div_rem[0] = div_rem / 1000L;
			action_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over matching timeline entries

				RecordIterator<TimelineEntry> entries = TimelineEntry.fetch_timeline_entry_range (action_time_lo, action_time_hi, event_id, comcat_ids, action_time_div_rem);
			){

				// Display them

				for (TimelineEntry entry : entries) {
					System.out.println (entry.toString());
				}
			}
		}

		return;
	}




	// Test #23 - Search the timeline for action time and/or event id; and re-fetch the entries.

	public static void test23(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test23' subcommand");
			return;
		}

		long action_time_lo = Long.parseLong(args[1]);
		long action_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] action_time_div_rem = null;
		if (div_rem > 0L) {
			action_time_div_rem = new long[2];
			action_time_div_rem[0] = div_rem / 1000L;
			action_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching timeline entries

			List<TimelineEntry> entries = TimelineEntry.get_timeline_entry_range (action_time_lo, action_time_hi, event_id, comcat_ids, action_time_div_rem);

			// Display them, and re-fetch

			for (TimelineEntry entry : entries) {
				System.out.println (entry.toString());
				TimelineEntry refetch = TimelineEntry.get_timeline_entry_for_key (entry.get_record_key());
				System.out.println (refetch.toString());
			}
		}

		return;
	}




	// Test #24 - Search the timeline for action time and/or event id; and delete the entries.

	public static void test24(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test24' subcommand");
			return;
		}

		long action_time_lo = Long.parseLong(args[1]);
		long action_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] action_time_div_rem = null;
		if (div_rem > 0L) {
			action_time_div_rem = new long[2];
			action_time_div_rem[0] = div_rem / 1000L;
			action_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching timeline entries

			List<TimelineEntry> entries = TimelineEntry.get_timeline_entry_range (action_time_lo, action_time_hi, event_id, comcat_ids, action_time_div_rem);

			// Display them, and delete

			for (TimelineEntry entry : entries) {
				System.out.println (entry.toString());
				TimelineEntry.delete_timeline_entry (entry);
			}
		}

		return;
	}




	// Test #25 - Display the pending task queue, sorted by execution time, using iterator.

	public static void test25(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test25' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over pending tasks

				RecordIterator<PendingTask> tasks = PendingTask.fetch_all_tasks();
			){

				// Display them

				for (PendingTask task : tasks) {
					System.out.println (task.toString());
				}
			}
		}

		return;
	}




	// Test #26 - Search the log for log time and/or event id, using iterator.

	public static void test26(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test26' subcommand");
			return;
		}

		long log_time_lo = Long.parseLong(args[1]);
		long log_time_hi = Long.parseLong(args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over matching log entries

				RecordIterator<LogEntry> entries = LogEntry.fetch_log_entry_range (log_time_lo, log_time_hi, event_id);
			){

				// Display them

				for (LogEntry entry : entries) {
					System.out.println (entry.toString());
				}
			}
		}

		return;
	}




	// Test #27 - Search the catalog snapshots for end time and/or event id, using iterator.

	public static void test27(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test27' subcommand");
			return;
		}

		double end_time_lo_days = Double.parseDouble (args[1]);
		double end_time_hi_days = Double.parseDouble (args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		long end_time_lo = Math.round(end_time_lo_days * 86400000L);
		long end_time_hi = Math.round(end_time_hi_days * 86400000L);

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over matching catalog snapshots

				RecordIterator<CatalogSnapshot> entries = CatalogSnapshot.fetch_catalog_snapshot_range (end_time_lo, end_time_hi, event_id);
			){

				// Display them

				for (CatalogSnapshot entry : entries) {
					System.out.println (entry.toString());
				}
			}
		}

		return;
	}




	// Test #28 - Post a console message task with given stage and message.

	public static void test28(String[] args) {

		// Two additional arguments

		if (args.length != 3) {
			System.err.println ("ServerTest : Invalid 'test28' subcommand");
			return;
		}

		int stage = Integer.parseInt(args[1]);

		OpConsoleMessage payload = new OpConsoleMessage();
		payload.message = args[2];

		String event_id = "";
		int opcode = TaskDispatcher.OPCODE_CON_MESSAGE;

		// Post the task

		long the_time = ServerClock.get_time();

		TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		return;
	}




	// Test #29 - Search the task queue for execution time and/or event id.

	public static void test29(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test29' subcommand");
			return;
		}

		long exec_time_lo = Long.parseLong(args[1]);
		long exec_time_hi = Long.parseLong(args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching tasks

			List<PendingTask> tasks = PendingTask.get_task_entry_range (exec_time_lo, exec_time_hi, event_id);

			// Display them

			for (PendingTask task : tasks) {
				System.out.println (task.toString());
			}

		}

		return;
	}




	// Test #30 - Search the task queue for execution time and/or event id, using iterator.

	public static void test30(String[] args) {

		// Two or three additional arguments

		if (args.length != 3 && args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test30' subcommand");
			return;
		}

		long exec_time_lo = Long.parseLong(args[1]);
		long exec_time_hi = Long.parseLong(args[2]);
		String event_id = null;
		if (args.length == 4) {
			event_id = args[3];
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over matching tasks

				RecordIterator<PendingTask> tasks = PendingTask.fetch_task_entry_range (exec_time_lo, exec_time_hi, event_id);
			){

				// Display them

				for (PendingTask task : tasks) {
					System.out.println (task.toString());
				}
			}
		}

		return;
	}




	// Test #31 - Search the timeline for action time and/or event id; get most recent.

	public static void test31(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test31' subcommand");
			return;
		}

		long action_time_lo = Long.parseLong(args[1]);
		long action_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] action_time_div_rem = null;
		if (div_rem > 0L) {
			action_time_div_rem = new long[2];
			action_time_div_rem[0] = div_rem / 1000L;
			action_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the most recent matching timeline entry

			TimelineEntry entry = TimelineEntry.get_recent_timeline_entry (action_time_lo, action_time_hi, event_id, comcat_ids, action_time_div_rem);

			// Display it

			if (entry == null) {
				System.out.println ("null");
			} else {
				System.out.println (entry.toString());
			}

		}

		return;
	}




	// Test #32 - Post a sync intake command for the given event.

	public static void test32(String[] args) {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test32' subcommand");
			return;
		}

		String event_id = args[1];

		OpIntakeSync payload = new OpIntakeSync();
		payload.setup ();

		int opcode = TaskDispatcher.OPCODE_INTAKE_SYNC;
		int stage = 0;

		// Post the task

		long the_time = ServerClock.get_time();

		TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		return;
	}




	// Test #33 - Parse a PDL intake command for the given command line.

	public static void test33(String[] args) {

		// At least one additional argument

		if (args.length < 2) {
			System.err.println ("ServerTest : Invalid 'test33' subcommand");
			return;
		}

		OpIntakePDL payload = new OpIntakePDL();

		payload.setup (args, 1, args.length);

		System.out.println ("PDL arguments:");
		for (String s : payload.pdl_args) {
			System.out.println (s);
		}

		System.out.println ("Parsed values:");
		System.out.println ("pdl_status = " + payload.pdl_status);
		System.out.println ("pdl_action = " + payload.pdl_action);
		System.out.println ("pdl_type = " + payload.pdl_type);
		System.out.println ("event_id = " + payload.event_id);
		System.out.println ("mainshock_time = " + payload.mainshock_time);
		System.out.println ("mainshock_mag = " + payload.mainshock_mag);
		System.out.println ("mainshock_lat = " + payload.mainshock_lat);
		System.out.println ("mainshock_lon = " + payload.mainshock_lon);
		System.out.println ("mainshock_depth = " + payload.mainshock_depth);

		return;
	}




	// Test #34 - Post a PDL intake command for the given command line.

	public static void test34(String[] args) {

		// At least one additional argument

		if (args.length < 2) {
			System.err.println ("ServerTest : Invalid 'test34' subcommand");
			return;
		}

		OpIntakePDL payload = new OpIntakePDL();

		payload.setup (args, 1, args.length);

		System.out.println ("PDL arguments:");
		for (String s : payload.pdl_args) {
			System.out.println (s);
		}

		System.out.println ("Parsed values:");
		System.out.println ("pdl_status = " + payload.pdl_status);
		System.out.println ("pdl_action = " + payload.pdl_action);
		System.out.println ("pdl_type = " + payload.pdl_type);
		System.out.println ("event_id = " + payload.event_id);
		System.out.println ("mainshock_time = " + payload.mainshock_time);
		System.out.println ("mainshock_mag = " + payload.mainshock_mag);
		System.out.println ("mainshock_lat = " + payload.mainshock_lat);
		System.out.println ("mainshock_lon = " + payload.mainshock_lon);
		System.out.println ("mainshock_depth = " + payload.mainshock_depth);

		String event_id = payload.event_id;

		int opcode = TaskDispatcher.OPCODE_INTAKE_PDL;
		int stage = 0;

		// Post the task

		long the_time = ServerClock.get_time();

		TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		return;
	}




	// Test #35 - Add a few elements to the alias families.

	public static void test35(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test35' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
		
			String timeline_id;
			String[] comcat_ids;
			String[] removed_ids;
			long family_time;
			AliasAssignment assignment;
			AliasAssignmentList assignments;
		

			timeline_id = "event_1";
			comcat_ids = new String[]{"ccid_11", "ccid_12", "ccid_13"};
			removed_ids = new String[0];
			family_time = 10101L;
			assignment = new AliasAssignment();
			assignments = new AliasAssignmentList();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			AliasFamily.submit_alias_family (null, family_time, assignments);
		

			timeline_id = "event_2";
			comcat_ids = new String[]{"ccid_21", "ccid_22", "ccid_23"};
			removed_ids = new String[]{"rmid_21", "rmid_22"};
			family_time = 20102L;
			assignment = new AliasAssignment();
			assignments = new AliasAssignmentList();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			AliasFamily.submit_alias_family (null, family_time, assignments);
		

			timeline_id = "event_3";
			comcat_ids = new String[]{"ccid_31"};
			removed_ids = new String[]{"rmid_31"};
			family_time = 30103L;
			assignment = new AliasAssignment();
			assignments = new AliasAssignmentList();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			AliasFamily.submit_alias_family (null, family_time, assignments);
		

			timeline_id = "event_1";
			comcat_ids = new String[]{"ccid_42", "ccid_11"};
			removed_ids = new String[]{"ccid_12", "ccid_13"};
			family_time = 40104L;
			assignment = new AliasAssignment();
			assignments = new AliasAssignmentList();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			timeline_id = "event_4";
			comcat_ids = new String[]{"ccid_41", "ccid_12"};
			removed_ids = new String[]{"ccid_13"};
			assignment = new AliasAssignment();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			AliasFamily.submit_alias_family (null, family_time, assignments);
		

			timeline_id = "event_2";
			comcat_ids = new String[]{"ccid_21"};
			removed_ids = new String[]{"rmid_21", "rmid_22", "ccid_23", "ccid_22"};
			family_time = 50105L;
			assignment = new AliasAssignment();
			assignments = new AliasAssignmentList();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			timeline_id = "event_3";
			comcat_ids = new String[]{"ccid_31", "ccid_23"};
			removed_ids = new String[]{"rmid_31"};
			assignment = new AliasAssignment();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			timeline_id = "event_5";
			comcat_ids = new String[]{"ccid_22"};
			removed_ids = new String[0];
			assignment = new AliasAssignment();
			assignment.set_timeline_id (timeline_id);
			assignment.set_comcat_ids_from_array (comcat_ids);
			for (String id : removed_ids) {
				assignment.add_removed_id (id);
			}
			assignments.add_assignment (assignment);
			AliasFamily.submit_alias_family (null, family_time, assignments);

		}

		return;
	}




	// Test #36 - Search the alias families for family time and/or event id; using list.

	public static void test36(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test36' subcommand");
			return;
		}

		long family_time_lo = Long.parseLong(args[1]);
		long family_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] family_time_div_rem = null;
		if (div_rem > 0L) {
			family_time_div_rem = new long[2];
			family_time_div_rem[0] = div_rem / 1000L;
			family_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching alias family entries

			List<AliasFamily> entries = AliasFamily.get_alias_family_range (family_time_lo, family_time_hi, event_id, comcat_ids, family_time_div_rem);

			// Display them

			for (AliasFamily entry : entries) {
				System.out.println (entry.toString());
				System.out.println (entry.get_assignments().toString());
			}

		}

		return;
	}




	// Test #37 - Search the alias families for family time and/or event id; using iterator.

	public static void test37(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test37' subcommand");
			return;
		}

		long family_time_lo = Long.parseLong(args[1]);
		long family_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] family_time_div_rem = null;
		if (div_rem > 0L) {
			family_time_div_rem = new long[2];
			family_time_div_rem[0] = div_rem / 1000L;
			family_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			try (

				// Get an iterator over matching alias family entries

				RecordIterator<AliasFamily> entries = AliasFamily.fetch_alias_family_range (family_time_lo, family_time_hi, event_id, comcat_ids, family_time_div_rem);
			){

				// Display them

				for (AliasFamily entry : entries) {
					System.out.println (entry.toString());
					System.out.println (entry.get_assignments().toString());
				}
			}
		}

		return;
	}




	// Test #38 - Search the alias families for family time and/or event id; and re-fetch the entries.

	public static void test38(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test38' subcommand");
			return;
		}

		long family_time_lo = Long.parseLong(args[1]);
		long family_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] family_time_div_rem = null;
		if (div_rem > 0L) {
			family_time_div_rem = new long[2];
			family_time_div_rem[0] = div_rem / 1000L;
			family_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching alias family entries

			List<AliasFamily> entries = AliasFamily.get_alias_family_range (family_time_lo, family_time_hi, event_id, comcat_ids, family_time_div_rem);

			// Display them, and re-fetch

			for (AliasFamily entry : entries) {
				System.out.println (entry.toString());
				System.out.println (entry.get_assignments().toString());
				AliasFamily refetch = AliasFamily.get_alias_family_for_key (entry.get_record_key());
				System.out.println (refetch.toString());
				System.out.println (refetch.get_assignments().toString());
			}
		}

		return;
	}




	// Test #39 - Search the alias families for family time and/or event id; and delete the entries.

	public static void test39(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test39' subcommand");
			return;
		}

		long family_time_lo = Long.parseLong(args[1]);
		long family_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] family_time_div_rem = null;
		if (div_rem > 0L) {
			family_time_div_rem = new long[2];
			family_time_div_rem[0] = div_rem / 1000L;
			family_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the list of matching alias family entries

			List<AliasFamily> entries = AliasFamily.get_alias_family_range (family_time_lo, family_time_hi, event_id, comcat_ids, family_time_div_rem);

			// Display them, and delete

			for (AliasFamily entry : entries) {
				System.out.println (entry.toString());
				System.out.println (entry.get_assignments().toString());
				AliasFamily.delete_alias_family (entry);
			}
		}

		return;
	}




	// Test #40 - Search the alias families for family time and/or event id; get most recent.

	public static void test40(String[] args) {

		// Three or more additional arguments

		if (args.length < 4) {
			System.err.println ("ServerTest : Invalid 'test40' subcommand");
			return;
		}

		long family_time_lo = Long.parseLong(args[1]);
		long family_time_hi = Long.parseLong(args[2]);
		long div_rem = Long.parseLong(args[3]);
		long[] family_time_div_rem = null;
		if (div_rem > 0L) {
			family_time_div_rem = new long[2];
			family_time_div_rem[0] = div_rem / 1000L;
			family_time_div_rem[1] = div_rem % 1000L;
		}
		String event_id = null;
		if (args.length >= 5) {
			if (!( args[4].equalsIgnoreCase("-") )) {
				event_id = args[4];
			}
		}
		String[] comcat_ids = null;
		if (args.length >= 6) {
			comcat_ids = Arrays.copyOfRange (args, 5, args.length);
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get the most recent matching alias family entry

			AliasFamily entry = AliasFamily.get_recent_alias_family (family_time_lo, family_time_hi, event_id, comcat_ids, family_time_div_rem);

			// Display it

			if (entry == null) {
				System.out.println ("null");
			} else {
				System.out.println (entry.toString());
				System.out.println (entry.get_assignments().toString());
			}

		}

		return;
	}




	// Test #41 - Delete a product from PDL-Development.

	public static void test41(String[] args) throws Exception {

		// Three additional arguments

		if (args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test41' subcommand");
			return;
		}

		String eventID = args[1];
		String eventNetwork = args[2];
		String eventCode = args[3];

		// Direct operation to PDL-Development

		ServerConfig server_config = new ServerConfig();
		server_config.get_server_config_file().pdl_enable = ServerConfigFile.PDLOPT_DEV;

		// Construct the deletion product

		boolean isReviewed = false;
		long modifiedTime = 0L;
		Product product = PDLProductBuilderOaf.createDeletionProduct (eventID, eventNetwork, eventCode, isReviewed, modifiedTime);

		// Send to PDL

		PDLSender.signProduct(product);
		PDLSender.sendProduct(product, true);

		return;
	}




	// Test #42 - Read an alias list from a file, and store it in the database.

	public static void test42(String[] args) throws Exception {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test42' subcommand");
			return;
		}

		String filename = args[1];

		// Read the alias list

		AliasAssignmentList assignments;
		long family_time;

		try (
			BufferedReader file_reader = new BufferedReader (new FileReader (filename));
		){
			MarshalImpJsonReader reader = new MarshalImpJsonReader (file_reader);
			
			assignments = (new AliasAssignmentList()).unmarshal (reader, null);

			family_time = assignments.get_family_time();
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Write the database entry
		
			AliasFamily.submit_alias_family (null, family_time, assignments);
		}

		return;
	}




	// Test #43 - Get current alias information for a timeline.

	public static void test43(String[] args) throws Exception {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test43' subcommand");
			return;
		}

		String timeline_id = args[1];

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get a task dispatcher and server group

			TaskDispatcher dispatcher = new TaskDispatcher();
			ServerGroup sg = dispatcher.get_server_group();

			// Set up task context

			dispatcher.setup_task_context();

			// Get timeline alias information

			ForecastMainshock fcmain = new ForecastMainshock();

			int rescode = sg.alias_sup.get_mainshock_for_timeline_id (timeline_id, fcmain);

			// Write the result code

			System.out.println ("result code = " + sg.get_rescode_as_string (rescode));

			// Display mainshock info if we got it

			if (rescode == ServerComponent.RESCODE_SUCCESS) {
				System.out.println (fcmain.toString());
			}
		}

		return;
	}




	// Test #44 - Get current alias information for an event.

	public static void test44(String[] args) throws Exception {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test44' subcommand");
			return;
		}

		String event_id = args[1];

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get a task dispatcher and server group

			TaskDispatcher dispatcher = new TaskDispatcher();
			ServerGroup sg = dispatcher.get_server_group();

			// Set up task context

			dispatcher.setup_task_context();

			// Get timeline alias information

			ForecastMainshock fcmain = new ForecastMainshock();

			int rescode = sg.alias_sup.get_mainshock_for_event_id (event_id, fcmain);

			// Write the result code

			System.out.println ("result code = " + sg.get_rescode_as_string (rescode));

			// Display mainshock info if we got it

			if (rescode == ServerComponent.RESCODE_SUCCESS || rescode == ServerComponent.RESCODE_ALIAS_NEW_EVENT) {
				System.out.println (fcmain.toString());
			}
		}

		return;
	}




	// Test #45 - Get current alias information for an event, create timeline if new event.

	public static void test45(String[] args) throws Exception {

		// One additional argument

		if (args.length != 2) {
			System.err.println ("ServerTest : Invalid 'test45' subcommand");
			return;
		}

		String event_id = args[1];

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get a task dispatcher and server group

			TaskDispatcher dispatcher = new TaskDispatcher();
			ServerGroup sg = dispatcher.get_server_group();

			// Set up task context

			dispatcher.setup_task_context();

			// Get timeline alias information

			ForecastMainshock fcmain = new ForecastMainshock();

			int rescode = sg.alias_sup.get_mainshock_for_event_id (event_id, fcmain);

			// Write the result code

			System.out.println ("result code = " + sg.get_rescode_as_string (rescode));

			// If it's a new event, create the timeline

			if (rescode == ServerComponent.RESCODE_ALIAS_NEW_EVENT) {
				sg.alias_sup.write_mainshock_to_new_timeline (fcmain);
			}

			// Display mainshock info if we got it

			if (rescode == ServerComponent.RESCODE_SUCCESS || rescode == ServerComponent.RESCODE_ALIAS_NEW_EVENT) {
				System.out.println (fcmain.toString());
			}
		}

		return;
	}




	// Test #46 - Test console redirection and time split output streams.

	public static void test46(String[] args) throws IOException {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test46' subcommand");
			return;
		}

		long day_millis = 86400000L;

		long day_1 = System.currentTimeMillis();
		long day_2 = day_1 + day_millis;
		long day_3 = day_2 + day_millis;
		long day_4 = day_3 + day_millis;
		long day_5 = day_4 + day_millis;
		long day_6 = day_5 + day_millis;
		long day_7 = day_6 + day_millis;
		long day_8 = day_7 + day_millis;

		String pattern_test = "'logtest/logs/'yyyy-MM-dd'-test.log'";
		String pattern_out = "'logtest/logs/'yyyy-MM-dd'-out.log'";
		String pattern_err = "'logtest/logs/'yyyy-MM-dd'-err.log'";

		System.out.println ("Start time = " + day_1 + " (" + SimpleUtils.time_to_string (day_1) + ")");

		// No redirection

		System.out.println ("con - out - line 1");
		System.err.println ("con - err - line 1");

		// Redirection and write within a single day, create directories

		try (
			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (pattern_test, day_1);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));
		){
			System.out.println ("day 1 - out - line 2");
			System.err.println ("day 1 - err - line 2");
			System.out.println ("day 1 - out - line 3");
			System.err.println ("day 1 - err - line 3");

			con_tsop.redirect (day_1);
		
			System.out.println ("day 1 - out - line 4");
			System.err.println ("day 1 - err - line 4");
			System.out.println ("day 1 - out - line 5");
			System.err.println ("day 1 - err - line 5");
		}

		// Redirection and write within two days, append to file for first day

		try (
			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (pattern_test, day_1);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));
		){
			System.out.println ("day 1 - out - line 6");
			System.err.println ("day 1 - err - line 6");
			System.out.println ("day 1 - out - line 7");
			System.err.println ("day 1 - err - line 7");

			con_tsop.redirect (day_2);
		
			System.out.println ("day 2 - out - line 8");
			System.err.println ("day 2 - err - line 8");
			System.out.println ("day 2 - out - line 9");
			System.err.println ("day 2 - err - line 9");
		}

		// Redirection and write within three days, test lazy open by not writing for two days

		try (
			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop (pattern_test, day_3);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));
		){
			con_tsop.redirect (day_4);
		
			System.out.println ("day 4 - out - line 10");
			System.err.println ("day 4 - err - line 10");

			con_tsop.redirect (day_4);
		
			System.out.println ("day 4 - out - line 11");
			System.err.println ("day 4 - err - line 11");

			con_tsop.redirect (day_5);
		}

		// Redirection and write within two days, with tee and separated files

		try (
			TimeSplitOutputStream con_tsop_out = TimeSplitOutputStream.make_tsop (pattern_out, day_6);
			TimeSplitOutputStream con_tsop_err = TimeSplitOutputStream.make_tsop (pattern_err, day_6);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop_out, con_tsop_err, false, true);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop_out,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop_err,
										ConsoleRedirector.get_new_err (con_red));
		){
			System.out.println ("day 6 - tee - out - line 12");
			System.err.println ("day 6 - tee - err - line 12");
			System.out.println ("day 6 - tee - out - line 13");
			System.err.println ("day 6 - tee - err - line 13");

			con_tsop_out.redirect (day_7);
			con_tsop_err.redirect (day_7);
		
			System.out.println ("day 7 - tee - out - line 14");
			System.err.println ("day 7 - tee - err - line 14");

			con_tsop_out.redirect (day_7);
			con_tsop_err.redirect (day_7);

			System.out.println ("day 7 - tee - out - line 15");
			System.err.println ("day 7 - tee - err - line 15");
		}

		// Redirection with empty pattern

		try (
			TimeSplitOutputStream con_tsop = TimeSplitOutputStream.make_tsop ("", day_8);
			ConsoleRedirector con_red = ConsoleRedirector.make_redirector (con_tsop, false, false);
			Closeable auto_out = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_out (con_red));
			Closeable auto_err = TimeSplitOutputStream.add_auto_upstream (con_tsop,
										ConsoleRedirector.get_new_err (con_red));
		){
			System.out.println ("day 8 - empty - out - line 16");
			System.err.println ("day 8 - empty - err - line 16");
			System.out.println ("day 8 - empty - out - line 17");
			System.err.println ("day 8 - empty - err - line 17");
		}

		// No redirection

		System.out.println ("con - out - line 18");
		System.err.println ("con - err - line 18");

		return;
	}




	// Test #47 - Delete all database tables, allowing a fresh start.

	public static void test47(String[] args) throws Exception {

		// Three additional arguments

		if (args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test47' subcommand");
			return;
		}

		if (!( args[1].equals ("delete")
			&& args[2].equals ("all")
			&& args[3].equals ("tables") )) {
			System.err.println ("ServerTest : Wrong confirmation for 'test47' subcommand");
			return;
		}

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){
			
			System.out.println ("ServerTest : Deleting all database tables");

			// Get the list of task entries and delete

			List<PendingTask> task_entries = PendingTask.get_task_entry_range (0L, 0L, null);
		
			System.out.println ("ServerTest : Deleting " + task_entries.size() + " task entries");

			for (PendingTask entry : task_entries) {
				PendingTask.delete_task (entry);
			}

			// Get the list of log entries and delete

			List<LogEntry> log_entries = LogEntry.get_log_entry_range (0L, 0L, null);
		
			System.out.println ("ServerTest : Deleting " + log_entries.size() + " log entries");

			for (LogEntry entry : log_entries) {
				LogEntry.delete_log_entry (entry);
			}

			// Get the list of catalog snapshots and delete

			List<CatalogSnapshot> cat_entries = CatalogSnapshot.get_catalog_snapshot_range (0L, 0L, null);
		
			System.out.println ("ServerTest : Deleting " + cat_entries.size() + " catalog snapshot entries");

			for (CatalogSnapshot entry : cat_entries) {
				CatalogSnapshot.delete_catalog_snapshot (entry);
			}

			// Get the list of timeline entries and delete

			List<TimelineEntry> tline_entries = TimelineEntry.get_timeline_entry_range (0L, 0L, null, null, null);
		
			System.out.println ("ServerTest : Deleting " + tline_entries.size() + " timeline entries");

			for (TimelineEntry entry : tline_entries) {
				TimelineEntry.delete_timeline_entry (entry);
			}

			// Get the list of alias family entries and delete

			List<AliasFamily> alfam_entries = AliasFamily.get_alias_family_range (0L, 0L, null, null, null);
		
			System.out.println ("ServerTest : Deleting " + alfam_entries.size() + " alias family entries");

			for (AliasFamily entry : alfam_entries) {
				AliasFamily.delete_alias_family (entry);
			}
			
			System.out.println ("ServerTest : Deleted all database tables");

		}

		return;
	}




	// Test #48 - Post a task to start polling Comcat.

	public static void test48(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test48' subcommand");
			return;
		}

		String event_id = ServerComponent.EVID_POLL;

		OpPollComcatStart payload = new OpPollComcatStart();
		payload.setup ();

		// Post the task

		int opcode = TaskDispatcher.OPCODE_POLL_COMCAT_START;
		int stage = 0;

		long the_time = ServerClock.get_time();

		boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		// Display result

		System.out.println ("Post poll Comcat start result: " + result);

		return;
	}




	// Test #49 - Post a task to stop polling Comcat.

	public static void test49(String[] args) {

		// No additional arguments

		if (args.length != 1) {
			System.err.println ("ServerTest : Invalid 'test49' subcommand");
			return;
		}

		String event_id = ServerComponent.EVID_POLL;

		OpPollComcatStop payload = new OpPollComcatStop();
		payload.setup ();

		// Post the task

		int opcode = TaskDispatcher.OPCODE_POLL_COMCAT_STOP;
		int stage = 0;

		long the_time = ServerClock.get_time();

		boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		// Display result

		System.out.println ("Post poll Comcat stop result: " + result);

		return;
	}




	// Test #50 - Post an analyst intervention task.

	public static void test50(String[] args) {

		// Three additional arguments

		if (args.length != 4) {
			System.err.println ("ServerTest : Invalid 'test50' subcommand");
			return;
		}

		String event_id = args[1];
		int state_change = Integer.parseInt(args[2]);
		boolean f_create_timeline = Boolean.parseBoolean (args[3]);

		// Set up the payload

		OpAnalystIntervene payload = new OpAnalystIntervene();
		payload.setup (state_change, f_create_timeline);

		// Post the task

		int opcode = TaskDispatcher.OPCODE_ANALYST_INTERVENE;
		int stage = 0;

		long the_time = ServerClock.get_time();

		boolean result = TaskDispatcher.post_task (event_id, the_time, the_time, "ServerTest", opcode, stage, payload.marshal_task());

		// Display result

		System.out.println ("Post analyst intervention result: " + result);

		return;
	}




	// Test dispatcher.
	
	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ServerTest : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1
		// Add a few items to the pending task queue.

		if (args[0].equalsIgnoreCase ("test1")) {

			try {
				test1(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #2
		// Command format:
		//  test2
		// Display the pending task queue, unsorted.

		if (args[0].equalsIgnoreCase ("test2")) {

			try {
				test2(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #3
		// Command format:
		//  test3
		// Display the pending task queue, sorted by execution time.

		if (args[0].equalsIgnoreCase ("test3")) {

			try {
				test3(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #4
		// Command format:
		//  test4
		// Display the first task in the pending task queue, according to execution time.

		if (args[0].equalsIgnoreCase ("test4")) {

			try {
				test4(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #5
		// Command format:
		//  test5  cutoff_time
		// Display the first task in the pending task queue, before cutoff time, according to execution time.

		if (args[0].equalsIgnoreCase ("test5")) {

			try {
				test5(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #6
		// Command format:
		//  test6  cutoff_time
		// Activate the first document before the cutoff time, and display the retrieved document.

		if (args[0].equalsIgnoreCase ("test6")) {

			try {
				test6(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #7
		// Command format:
		//  test7  cutoff_time  exec_time  stage  [event_id]
		// Activate the first document before the cutoff time, and stage it.

		if (args[0].equalsIgnoreCase ("test7")) {

			try {
				test7(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #8
		// Command format:
		//  test6  cutoff_time
		// Activate the first document before the cutoff time, and delete it.

		if (args[0].equalsIgnoreCase ("test8")) {

			try {
				test8(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #9
		// Command format:
		//  test9
		// Run task dispatcher.

		if (args[0].equalsIgnoreCase ("test9")) {

			try {
				test9(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #10
		// Command format:
		//  test10
		// Post a shutdown task.

		if (args[0].equalsIgnoreCase ("test10")) {

			try {
				test10(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #11
		// Command format:
		//  test11
		// Scan the pending task queue, sorted, and write a log entry for each.

		if (args[0].equalsIgnoreCase ("test11")) {

			try {
				test11(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #12
		// Command format:
		//  test12
		// Scan the pending task queue, sorted, and search the log for each.

		if (args[0].equalsIgnoreCase ("test12")) {

			try {
				test12(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #13
		// Command format:
		//  test13  log_time_lo  log_time_hi  [event_id]
		// Search the log for log time and/or event id.
		// Log times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test13")) {

			try {
				test13(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #14
		// Command format:
		//  test14  log_time_lo  log_time_hi  [event_id]
		// Search the log for log time and/or event id, and delete the matching entries.
		// Log times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test14")) {

			try {
				test14(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #15
		// Command format:
		//  test15  event_id  opcode  stage  [details]
		// Post a task with given event id, opcode, stage, and details.
		// Event id can be "-" for an empty string.

		if (args[0].equalsIgnoreCase ("test15")) {

			try {
				test15(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #16
		// Command format:
		//  test16  [f_adjust_time]
		// Execute the next task.
		// If f_adjust_time is "true" then adjust clock to be the execution time of the task.
		// If f_adjust_time is omitted then the default value is "false".

		if (args[0].equalsIgnoreCase ("test16")) {

			try {
				test16(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #17
		// Command format:
		//  test17  start_time_days  end_time_days  event_id
		// Write a catalog snapshot.

		if (args[0].equalsIgnoreCase ("test17")) {

			try {
				test17(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #18
		// Command format:
		//  test18  end_time_lo_days  end_time_hi_days  [event_id]
		// Search the catalog snapshots for end time and/or event id.
		// Times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test18")) {

			try {
				test18(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #19
		// Command format:
		//  test19  end_time_lo_days  end_time_hi_days  [event_id]
		// Search the catalog snapshots for end time and/or event id, and delete the matching entries.
		// Times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test19")) {

			try {
				test19(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #20
		// Command format:
		//  test20
		// Add a few elements to the timeline.

		if (args[0].equalsIgnoreCase ("test20")) {

			try {
				test20(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #21
		// Command format:
		//  test21  action_time_lo  action_time_hi  action_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; using list.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The action_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test21")) {

			try {
				test21(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #22
		// Command format:
		//  test22  action_time_lo  action_time_hi  action_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; using iterator.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The action_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test22")) {

			try {
				test22(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #23
		// Command format:
		//  test23  action_time_lo  action_time_hi  action_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; and re-fetch the entries.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The action_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test23")) {

			try {
				test23(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #24
		// Command format:
		//  test24  action_time_lo  action_time_hi  action_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; and delete the entries.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The action_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test24")) {

			try {
				test24(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #25
		// Command format:
		//  test25
		// Display the pending task queue, sorted by execution time, using iterator.

		if (args[0].equalsIgnoreCase ("test25")) {

			try {
				test25(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #26
		// Command format:
		//  test26  log_time_lo  log_time_hi  [event_id]
		// Search the log for log time and/or event id, using iterator.
		// Log times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test26")) {

			try {
				test26(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #27
		// Command format:
		//  test27  end_time_lo_days  end_time_hi_days  [event_id]
		// Search the catalog snapshots for end time and/or event id, using iterator.
		// Times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test27")) {

			try {
				test27(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #28
		// Command format:
		//  test28  stage  message
		// Post a console message task with given stage and message.

		if (args[0].equalsIgnoreCase ("test28")) {

			try {
				test28(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #29
		// Command format:
		//  test29  exec_time_lo  exec_time_hi  [event_id]
		// Search the task queue for execution time and/or event id.
		// Execution times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test29")) {

			try {
				test29(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #30
		// Command format:
		//  test30  exec_time_lo  exec_time_hi  [event_id]
		// Search the task queue for execution time and/or event id, using iterator.
		// Execution times can be 0 for no bound, event id can be omitted for no restriction.

		if (args[0].equalsIgnoreCase ("test30")) {

			try {
				test30(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #31
		// Command format:
		//  test31  action_time_lo  action_time_hi  action_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; get most recent.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The action_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test31")) {

			try {
				test31(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #32
		// Command format:
		//  test32  event_id
		// Post a sync intake command for the given event.

		if (args[0].equalsIgnoreCase ("test32")) {

			try {
				test32(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #33
		// Command format:
		//  test33  arg...
		// Parse a PDL intake command for the given command line.

		if (args[0].equalsIgnoreCase ("test33")) {

			try {
				test33(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #34
		// Command format:
		//  test34  arg...
		// Post a PDL intake command for the given command line.

		if (args[0].equalsIgnoreCase ("test34")) {

			try {
				test34(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #35
		// Command format:
		//  test35
		// Add a few elements to the alias families.

		if (args[0].equalsIgnoreCase ("test35")) {

			try {
				test35(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #36
		// Command format:
		//  test36  family_time_lo  family_time_hi  family_time_div_rem  [event_id]  [comcat_id]...
		// Search the alias families for family time and/or event id and/or comcat id; using list.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The family_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test36")) {

			try {
				test36(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #37
		// Command format:
		//  test37  family_time_lo  family_time_hi  family_time_div_rem  [event_id]  [comcat_id]...
		// Search the alias families for family time and/or event id and/or comcat id; using iterator.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The family_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test37")) {

			try {
				test37(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #38
		// Command format:
		//  test38  family_time_lo  family_time_hi  family_time_div_rem  [event_id]  [comcat_id]...
		// Search the alias families for family time and/or event id and/or comcat id; and re-fetch the entries.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The family_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test38")) {

			try {
				test38(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #39
		// Command format:
		//  test39  family_time_lo  family_time_hi  family_time_div_rem  [event_id]  [comcat_id]...
		// Search the alias families for family time and/or event id and/or comcat id; and delete the entries.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The family_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test39")) {

			try {
				test39(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #40
		// Command format:
		//  test40  family_time_lo  family_time_hi  family_time_div_rem  [event_id]  [comcat_id]...
		// Search the timeline for action time and/or event id and/or comcat id; get most recent.
		// Times can be 0 for no bound, event id can be omitted or equal to "-" for no restriction.
		// If any comcat_id are given, the entry must match at least one of them.
		// The family_time_div_rem can be 0 for no modulus test, otherwise it is divisor * 1000 + remainder.

		if (args[0].equalsIgnoreCase ("test40")) {

			try {
				test40(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #41
		// Command format:
		//  test41  eventID  eventNetwork  eventCode
		// Delete a product from PDL-Development.
		// Note this always uses PDL-Development regardless of the ServerConfig setting.

		if (args[0].equalsIgnoreCase ("test41")) {

			try {
				test41(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #42
		// Command format:
		//  test42  filename
		// Read an alias list from a file, and store it in the database.

		if (args[0].equalsIgnoreCase ("test42")) {

			try {
				test42(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #43
		// Command format:
		//  test43  timeline_id
		// Get current alias information for a timeline.

		if (args[0].equalsIgnoreCase ("test43")) {

			try {
				test43(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #44
		// Command format:
		//  test44  event_id
		// Get current alias information for an event.

		if (args[0].equalsIgnoreCase ("test44")) {

			try {
				test44(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #45
		// Command format:
		//  test45  event_id
		// Get current alias information for an event, create timeline if new event.

		if (args[0].equalsIgnoreCase ("test45")) {

			try {
				test45(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #46
		// Command format:
		//  test46
		// Test console redirection and time split output streams.

		if (args[0].equalsIgnoreCase ("test46")) {

			try {
				test46(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #47
		// Command format:
		//  test47 "delete" "all" "tables"
		// Delete all the database tables, allowing a fresh start.

		if (args[0].equalsIgnoreCase ("test47")) {

			try {
				test47(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #48
		// Command format:
		//  test48
		// Post a task to start polling Comcat.

		if (args[0].equalsIgnoreCase ("test48")) {

			try {
				test48(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #49
		// Command format:
		//  test49
		// Post a task to stop polling Comcat.

		if (args[0].equalsIgnoreCase ("test49")) {

			try {
				test49(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Subcommand : Test #50
		// Command format:
		//  test50  event_id  state_change  f_create_timeline
		// Post a task for analyst intervention.
		// The event_id can be either a Comcat ID or a timeline ID.
		// The state_change is an integer:
		//  1 = no change, 2 = start timeline, 3 = stop timeline, 4 = withdraw timeline.
		// The f_create_timeline is true to create the timeline if it doesn't exist.

		if (args[0].equalsIgnoreCase ("test50")) {

			try {
				test50(args);
            } catch (Exception e) {
                e.printStackTrace();
			}

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ServerTest : Unrecognized subcommand : " + args[0]);
		return;
	}
}
