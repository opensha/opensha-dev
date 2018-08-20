package scratch.aftershockStatistics.aafs;

import java.util.List;

import java.io.IOException;
import java.io.PrintStream;
import java.io.OutputStream;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.TimeSplitOutputStream;

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatException;
import scratch.aftershockStatistics.ComcatConflictException;
import scratch.aftershockStatistics.ComcatRemovedException;
import scratch.aftershockStatistics.ComcatAccessor;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

/**
 * Support functions for logging.
 * Author: Michael Barall 06/23/2018.
 */
public class LogSupport extends ServerComponent {


	//----- Logging state -----


	// The summary log output destination, or null if none.

	private PrintStream summary_log_out = null;




	//----- Destination control -----


	// Set the summary log destination.
	// Destination can be null if no summary log is desired.

	public void set_summary_log_out (OutputStream dest) {

		if (dest == null) {
			summary_log_out = null;
		} else {
			boolean autoFlush = true;
			summary_log_out = new PrintStream (dest, autoFlush);
		}

		return;
	}




	//----- Internal reporting functions -----


	// Report an action.
	// Parameters:
	//  name = Name of action, must be non-null and non-empty
	//  items = List of items describing action, any null or empty items are ignored.

	private void report_action (String name, String... items) {

		// If we have a summary log ...

		if (summary_log_out != null) {

			StringBuilder result = new StringBuilder();

			// Time stamp

			result.append (SimpleUtils.time_to_string_no_z (ServerClock.get_true_time()));

			// Name

			result.append (" " + name);

			// Items

			String sep = ": ";

			for (String item : items) {
				if (!( item == null || item.isEmpty() )) {
					result.append (sep + item);
					sep = ", ";
				}
			}

			// Write the line

			summary_log_out.println (result.toString());
		}

		return;
	}



	// Report an exception.
	// Parameters:
	//  e = Exception, must be non-null.
	// Note: You should call report_action first.

	private void report_exception (Throwable e) {

		// If we have a summary log ...

		if (summary_log_out != null) {

			// Write the stack trace

			summary_log_out.println (SimpleUtils.getStackTraceAsString (e));
		}

		return;
	}



	// Report information.
	// Parameters:
	//  info = Information, can be multiple lines.
	// Note: You should call report_action first.

	private void report_info (String info) {

		// If we have a summary log ...

		if (summary_log_out != null) {

			// Write the information

			if (!( info == null || info.isEmpty() )) {
				summary_log_out.println (info);
			}
		}

		return;
	}




	//----- Reporting functions -----




	// Report an uncaught dispatcher exception.
	// Note: The task can be null.

	public void report_dispatcher_exception (PendingTask task, Throwable e) {
		report_action ("DISPATCHER-EXCEPTION");
		if (task != null) {
			report_info ("Failing task: " + task.toString());
		}
		report_exception (e);
		return;
	}




	// Report invalid task exception.

	public void report_invalid_task_exception (PendingTask task, Exception e) {
		report_action ("INVALID-TASK-EXCEPTION");
		report_info ("Invalid task: " + task.toString());
		report_info ("Dump: " + task.dump_details());
		report_exception (e);
		return;
	}




	// Report task begin.

	public void report_task_begin (PendingTask task) {
		report_action ("TASK-BEGIN",
					get_opcode_as_string (task.get_opcode()),
					"event_id = " + task.get_event_id(),
					"stage = " + task.get_stage());
		return;
	}




	// Report task restart.

	public void report_task_restart (PendingTask task) {
		report_action ("TASK-RESTART",
					get_opcode_as_string (task.get_opcode()),
					"event_id = " + task.get_event_id(),
					"stage = " + task.get_stage());
		return;
	}




	// Report task end.

	public void report_task_end (PendingTask task, int rescode) {
		report_action ("TASK-END",
					get_opcode_as_string (task.get_opcode()),
					get_rescode_as_string (rescode));
		return;
	}




	// Report task delete.

	public void report_task_delete (PendingTask task, int rescode) {
		report_action ("TASK-DELETE",
					get_opcode_as_string (task.get_opcode()),
					get_rescode_as_string (rescode));
		return;
	}




	// Report task stage.
	// Note: event_id can be null if the event ID is not being changed.

	public void report_task_stage (PendingTask task, int rescode, String event_id, int stage, long exec_time) {
		report_action ("TASK-STAGE",
					get_opcode_as_string (task.get_opcode()),
					get_rescode_as_string (rescode),
					(event_id == null) ? ("event_id = " + task.get_event_id()) : ("new_event_id = " + event_id),
					"stage = " + stage,
					"exec_time = " + SimpleUtils.time_to_string (exec_time));
		return;
	}




	// Report catalog snapshot saved.

	public void report_catalog_saved (String event_id, int eqk_count) {
		report_action ("CATALOG-SAVED",
					event_id,
					"eqk_count = " + eqk_count);
		return;
	}




	// Report timeline appended.

	public void report_timeline_appended (TimelineStatus tstatus) {

		if (tstatus.is_first_entry() && tstatus.forecast_mainshock != null) {
			report_action ("INTAKE-EVENT",
				SimpleUtils.event_id_and_info_one_line (
					tstatus.forecast_mainshock.mainshock_event_id,
					tstatus.forecast_mainshock.mainshock_time,
					tstatus.forecast_mainshock.mainshock_mag,
					tstatus.forecast_mainshock.mainshock_lat,
					tstatus.forecast_mainshock.mainshock_lon,
					tstatus.forecast_mainshock.mainshock_depth)
			);
		}

		report_action ("TIMELINE-APPENDED",
					tstatus.event_id,
					tstatus.get_actcode_as_string (),
					tstatus.get_fc_status_as_string (),
					tstatus.get_pdl_status_as_string (),
					(tstatus.result_has_shadowing())
						? (tstatus.get_fc_result_as_string () + " (" + tstatus.shadowing_event_id + ")")
						: (tstatus.get_fc_result_as_string ()),
					"action_time = " + SimpleUtils.time_raw_and_string (tstatus.action_time));
		return;
	}




	// Report forecast request.

	public void report_forecast_request (String event_id, long sched_time) {
		report_action ("FORECAST-REQUEST",
					event_id,
					"sched_time = " + SimpleUtils.time_raw_and_string (sched_time));
		return;
	}




	// Report PDL report request.

	public void report_pdl_report_request (String event_id, long sched_time) {
		report_action ("PDL-REPORT-REQUEST",
					event_id,
					"sched_time = " + SimpleUtils.time_raw_and_string (sched_time));
		return;
	}




	// Report expire request.

	public void report_expire_request (String event_id) {
		report_action ("EXPIRE-REQUEST",
					event_id);
		return;
	}




	// Report dispatcher start.

	public void report_dispatcher_start () {
		report_info (LOG_SEPARATOR_LINE);
		report_action ("DISPATCHER-START");
		report_info (VersionInfo.get_title());
		return;
	}




	// Report dispatcher shutdown.

	public void report_dispatcher_shutdown () {
		report_action ("DISPATCHER-SHUTDOWN");
		return;
	}




	// Report dispatcher immediate shutdown.

	public void report_dispatcher_immediate_shutdown () {
		report_action ("DISPATCHER-IMMEDIATE-SHUTDOWN");
		return;
	}




	// Report dispatcher restart.

	public void report_dispatcher_restart () {
		report_action ("DISPATCHER-RESTART");
		return;
	}




	// Report timeline entry deleted.

	public void report_timeline_entry_deleted (String event_id) {
		report_action ("TIMELINE-ENTRY-DELETED",
					event_id);
		return;
	}




	// Report timeline task unwind.

	public void report_timeline_unwind (String event_id) {
		report_action ("TIMELINE-TASK-UNWIND",
					event_id);
		return;
	}




	// Report corrupt timeline exception.

	public void report_corrupt_timeline_exception (TimelineEntry tentry, PendingTask task, Exception e) {
		report_action ("CORRUPT-TIMELINE-EXCEPTION",
					task.get_event_id());
		report_info ("Timeline entry synopsis:\n" + tentry.toString());
		report_info ("Timeline entry details dump:\n" + tentry.dump_details());
		report_exception (e);
		return;
	}




	// Report timeline error appended.

	public void report_timeline_error_appended (String event_id) {
		report_action ("TIMELINE-ERROR-APPENDED",
					event_id);
		return;
	}




	// Report Comcat exception.
	// Note: event_id can be null;

	public void report_comcat_exception (String event_id, Exception e) {

		if (e instanceof ComcatRemovedException) {
			report_action ("COMCAT-REMOVED",
						event_id);
			String message = e.getMessage();
			if (!( message == null || message.isEmpty() )) {
				report_info (message);
			}
		}

		else {
			report_action ("COMCAT-EXCEPTION",
						event_id);
			report_exception (e);
		}

		return;
	}




	// Report event shadowed.

	public void report_event_shadowed (String event_id, String shadow_id,
					double event_mag, double shadow_mag, double distance, double interval) {
		report_action ("EVENT-SHADOWED",
					"event_id = " + event_id,
					"shadow_id = " + shadow_id,
					"event_mag = " + String.format("%.2f", event_mag),
					"shadow_mag = " + String.format("%.2f", shadow_mag),
					"distance = " + String.format("%.3f", distance) + " km",
					"interval = " + String.format("%.3f", interval) + " days");
		return;
	}




	// Report event is a foreshock.

	public void report_event_foreshock (String event_id, String cat_max_id,
					double event_mag, double cat_max_mag) {
		report_action ("EVENT-FORESHOCK",
					"event_id = " + event_id,
					"cat_max_id = " + cat_max_id,
					"event_mag = " + String.format("%.2f", event_mag),
					"cat_max_mag = " + String.format("%.2f", cat_max_mag));
		return;
	}




	// Report successful PDL send.

	public void report_pdl_send_ok (TimelineStatus tstatus) {
		report_action ("PDL-SEND-OK",
					"eventID = " + sg.alias_sup.timeline_id_to_pdl_code (tstatus.event_id),
					"eventNetwork = " + tstatus.forecast_mainshock.mainshock_network,
					"eventCode = " + tstatus.forecast_mainshock.mainshock_code);
		return;
	}




	// Report PDL send exception.

	public void report_pdl_send_exception (TimelineStatus tstatus, Exception e) {
		report_action ("PDL-SEND-EXCEPTION",
					"eventID = " + sg.alias_sup.timeline_id_to_pdl_code (tstatus.event_id),
					"eventNetwork = " + tstatus.forecast_mainshock.mainshock_network,
					"eventCode = " + tstatus.forecast_mainshock.mainshock_code);
		report_exception (e);
		return;
	}




	// Report Comcat poll done.

	public void report_comcat_poll_done (long poll_lookback, int count_no_timeline, int count_withdrawn_timeline) {
		double poll_lookback_days = ((double)poll_lookback)/ComcatAccessor.day_millis;
		report_action ("COMCAT-POLL-DONE",
					"poll_lookback = " + String.format ("%.3f", poll_lookback_days) + " days",
					"count_no_timeline = " + count_no_timeline,
					"count_withdrawn_timeline = " + count_withdrawn_timeline);
		return;
	}




	// Report alias family created.

	public void report_alias_family_created (long family_time, String info) {
		report_action ("ALIAS-CREATED",
					"family_time = " + SimpleUtils.time_raw_and_string (family_time));
		report_info (info);
		return;
	}




	// Report alias family updated.

	public void report_alias_family_updated (long family_time, String info) {
		report_action ("ALIAS-UPDATED",
					"family_time = " + SimpleUtils.time_raw_and_string (family_time));
		report_info (info);
		return;
	}




	//----- Construction -----


	// Default constructor.

	public LogSupport () {
		summary_log_out = null;
	}


	// Set up this component by linking to the server group.
	// A subclass may override this to perform additional setup operations.

	@Override
	public void setup (ServerGroup the_sg) {
		super.setup (the_sg);

		summary_log_out = null;

		return;
	}

}
