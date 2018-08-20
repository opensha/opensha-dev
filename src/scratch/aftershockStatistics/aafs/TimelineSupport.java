package scratch.aftershockStatistics.aafs;

import java.util.List;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;

import scratch.aftershockStatistics.ComcatException;
import scratch.aftershockStatistics.ComcatConflictException;
import scratch.aftershockStatistics.ComcatRemovedException;
import scratch.aftershockStatistics.CompactEqkRupList;

/**
 * Support functions for timelines.
 * Author: Michael Barall 06/23/2018.
 */
public class TimelineSupport extends ServerComponent {




	//----- Task execution subroutines : Timeline operations -----




	// Delete all waiting tasks with the given event id that are delayed timeline actions.
	// Note: The currently active task is not deleted, even if it would match.

	public void delete_delayed_timeline_tasks (String event_id) {
		sg.task_sup.delete_all_waiting_tasks (event_id, OPCODE_GEN_FORECAST, OPCODE_GEN_PDL_REPORT, OPCODE_GEN_EXPIRE);
		return;
	}




	// Get the most recent entry on the timeline, and extract its status structure.
	// If the task is restarted, this routine either completes the task or rolls back
	// the timeline state so the task can start over.
	// Return values are:
	//  RESCODE_TIMELINE_EXISTS = The timeline exists, and its state is in tstatus.
	//                            If restarted, anything that was previously done is undone.
	//  RESCODE_TIMELINE_NOT_FOUND = The timeline does not exist, and tstatus is indeterminate.
	//                               If restarted, anything that was previously done is undone.
	//  RESCODE_TASK_RETRY_SUCCESS = The task was restarted, and has now been completed successfully.
	//  RESCODE_TIMELINE_CORRUPT, RESCODE_TASK_CORRUPT (or anything else) = Error, operation cannot proceed.

	public int open_timeline (PendingTask task, TimelineStatus tstatus, DBPayload payload) {

		TimelineEntry tentry = null;
		CatalogSnapshot catsnap = null;

		// Get the payload

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// If this task is a delayed timeline operation ...

			switch (task.get_opcode()) {
			case OPCODE_GEN_FORECAST:
			case OPCODE_GEN_PDL_REPORT:
			case OPCODE_GEN_EXPIRE:

				// Get the most recent timeline entry for this event

				tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, task.get_event_id(), null, null);

				if (tentry != null) {

					// Get the status for this timeline entry

					try {
						tstatus.unmarshal_timeline (tentry);
					}

					// Invalid timeline entry

					catch (Exception e2) {

						// Display the error and log the task

						display_timeline_corrupt (tentry, task, e2);
						return RESCODE_TIMELINE_CORRUPT;
					}

					// Issue any new delayed command that is needed to replace the failed one

					next_auto_timeline (tstatus);
				}
				break;
			}

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		// If task is restarting ...

		if (task.is_restarted()) {

			// If a timeline entry was written for this task, get it

			tentry = TimelineEntry.get_timeline_entry_for_key (task.get_record_key());

			// If a catalog snapshot was written for this task, get it

			catsnap = CatalogSnapshot.get_catalog_shapshot_for_key (task.get_record_key());

			// If we got a timeline entry ...

			if (tentry != null) {

				// Get the status for this timeline entry

				try {
					tstatus.unmarshal_timeline (tentry);
				}

				// Invalid timeline entry

				catch (Exception e) {

					// Display the error and log the task

					display_timeline_corrupt (tentry, task, e);
					return RESCODE_TIMELINE_CORRUPT;
				}

				// If we have any required catalog snapshot, then we can complete the task ...

				if ( (!(tstatus.has_catalog_snapshot())) || catsnap != null ) {

					//--- Complete the pending task
				
					// Remove any delayed commands

					delete_delayed_timeline_tasks (tstatus.event_id);

					// Issue any new delayed command that is needed

					next_auto_timeline (tstatus);

					// Task is now completed

					return RESCODE_TASK_RETRY_SUCCESS;
				}

				// Unwind the task, by deleting the timeline entry

				sg.log_sup.report_timeline_entry_deleted (task.get_event_id());

				TimelineEntry.delete_timeline_entry (tentry);
			}

			//--- Undo the pending task, so it can start over from the beginning

			sg.log_sup.report_timeline_unwind (task.get_event_id());

			// Delete the catalog entry if we have it

			if (catsnap != null) {
				CatalogSnapshot.delete_catalog_snapshot (catsnap);
			}
				
			// Remove any delayed commands

			delete_delayed_timeline_tasks (tstatus.event_id);

			// Get the most recent timeline entry for this event

			tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, task.get_event_id(), null, null);

			if (tentry == null) {
				return RESCODE_TIMELINE_NOT_FOUND;
			}

			// Get the status for this timeline entry

			try {
				tstatus.unmarshal_timeline (tentry);
			}

			// Invalid timeline entry

			catch (Exception e) {

				// Display the error and log the task

				display_timeline_corrupt (tentry, task, e);
				return RESCODE_TIMELINE_CORRUPT;
			}

			// If the current task is not a delayed task, re-issue any delayed task for the current timeline entry

			switch (task.get_opcode()) {
			case OPCODE_GEN_FORECAST:
			case OPCODE_GEN_PDL_REPORT:
			case OPCODE_GEN_EXPIRE:
				break;

			default:
				next_auto_timeline (tstatus);
				break;
			}

			// Found timeline entry

			return RESCODE_TIMELINE_EXISTS;
		}

		//--- Prepare for the pending task

		// Get the most recent timeline entry for this event

		tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, task.get_event_id(), null, null);

		if (tentry == null) {
			return RESCODE_TIMELINE_NOT_FOUND;
		}

		// Get the status for this timeline entry

		try {
			tstatus.unmarshal_timeline (tentry);
		}

		// Invalid timeline entry

		catch (Exception e) {

			// Display the error and log the task

			display_timeline_corrupt (tentry, task, e);
			return RESCODE_TIMELINE_CORRUPT;
		}

		// Found timeline entry

		return RESCODE_TIMELINE_EXISTS;
	}




	// Display an exception caused by a corrupt timeline entry, and set the log remark.
	// Also, if the timeline entry is not an error entry, then write an error entry.

	public void display_timeline_corrupt (TimelineEntry tentry, PendingTask task, Exception e) {

		// Display messages
		
		sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry is corrupt:\n"
			+ "event_id = " + task.get_event_id() + "\n"
			+ "Timeline entry synopsis:\n" + tentry.toString() + "\n"
			+ "Timeline entry details dump:\n" + tentry.dump_details() + "\n"
			+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));

		sg.log_sup.report_corrupt_timeline_exception (tentry, task, e);

		// If the timeline entry is not marked as error, then add a new timeline entry

		if (tentry.get_actcode() != TimelineStatus.ACTCODE_ERROR) {

			// Delete any pending delayed tasks

			delete_delayed_timeline_tasks (tentry.get_event_id());

			// Write an error entry into the timeline

			String[] comcat_ids = new String[1];
			comcat_ids[0] = EVID_ERROR;

			TimelineEntry.submit_timeline_entry (
				task.get_record_key(),										// key
				Math.max(sg.task_disp.get_time(), tentry.get_action_time() + 1L),	// action_time
				tentry.get_event_id(),										// event_id
				comcat_ids,													// comcat_ids
				TimelineStatus.ACTCODE_ERROR,								// actcode
				null);														// details

			sg.log_sup.report_timeline_error_appended (tentry.get_event_id());
		}

		return;
	}




	// Process an exception caused by a ComCat failure.
	// Display a message, set the log remark, set the timeline to ComCat fail state, and write the timeline entry.

	private void process_timeline_comcat_fail (PendingTask task, TimelineStatus tstatus, Exception e) {

		// Display messages
			
		sg.log_sup.report_comcat_exception (task.get_event_id(), e);
		
		sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline stopped due to ComCat failure:\n"
			+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
			+ "event_id = " + task.get_event_id() + "\n"
			+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));

		// Set to ComCat fail state

		tstatus.set_state_comcat_fail (sg.task_disp.get_time());

		// Write timeline entry

		append_timeline (task, tstatus);
		return;
	}




	// Process an exception caused by an event removed from ComCat.
	// Display a message, set the log remark, set the timeline to withdrawn state, and write the timeline entry.

	private void process_timeline_comcat_removed (PendingTask task, TimelineStatus tstatus, Exception e) {

		// Display messages
			
		sg.log_sup.report_comcat_exception (task.get_event_id(), e);
		
		sg.task_disp.set_display_taskres_log ("TASK-INFO: Timeline stopped due to event deleted or merged in ComCat:\n"
			+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
			+ "event_id = " + task.get_event_id() + "\n"
			+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));

		// Withdraw the timeline

		tstatus.set_state_withdrawn (sg.task_disp.get_time(), null);

		// Write timeline entry

		append_timeline (task, tstatus);
		return;
	}




	// Process a retry caused by a ComCat failure.
	// Stage the task to retry the Comcat operation.
	// If retries are exhausted, then fail the operation.
	// Return values:
	//  RESCODE_STAGE_COMCAT_RETRY = The task is being staged for retry.
	//  RESCODE_TIMELINE_COMCAT_FAIL = Retries exhausted, the timeline is stopped in Comcat fail state.
	//  RESCODE_TIMELINE_EVENT_REMOVED = Retries exhausted, event deleted or merged in Comcat, timeline is in withdrawn state.

	public int process_timeline_comcat_retry (PendingTask task, TimelineStatus tstatus, Exception e) {

		// Get the next ComCat retry lag

		long min_lag = sg.task_disp.get_action_config().int_to_lag (task.get_stage()) + 1L;

		if (e instanceof ComcatRemovedException) {
			if (min_lag < sg.task_disp.get_action_config().get_comcat_retry_missing()) {
				min_lag = sg.task_disp.get_action_config().get_comcat_retry_missing();
			}
		}

		long next_comcat_retry_lag = sg.task_disp.get_action_config().get_next_comcat_retry_lag (min_lag);

		// If there is another retry, stage the task

		if (next_comcat_retry_lag >= 0L) {

			sg.log_sup.report_comcat_exception (task.get_event_id(), e);

			sg.task_disp.set_taskres_stage (task.get_sched_time() + next_comcat_retry_lag,
								sg.task_disp.get_action_config().lag_to_int (next_comcat_retry_lag));
			return RESCODE_STAGE_COMCAT_RETRY;
		}

		// Retries exhausted, process the error and log the task

		if (e instanceof ComcatRemovedException) {
			process_timeline_comcat_removed (task, tstatus, e);
			return RESCODE_TIMELINE_EVENT_REMOVED;
		}

		process_timeline_comcat_fail (task, tstatus, e);
		return RESCODE_TIMELINE_COMCAT_FAIL;
	}




	// Append an entry to the timeline.
	// The entry is tagged with the key from the current task.
	// Also write a catalog snapshot to the database if necessary.
	// Also issue the next delayed command for the timeline, if necessary.
	// Parameters:
	//  task = Current task.
	//  tstatus = Timeline status.
	//  last_pdl_lag = Lag for the last PDL report attempt, or -1L if none.  Defaults to -1L if omitted.
	// Note: Typically last_pdl_lag == 0L if the caller has already attempted a PDL report, -1L otherwise.

	public void append_timeline (PendingTask task, TimelineStatus tstatus) {

		append_timeline (task, tstatus, -1L);
		return;
	}

	public void append_timeline (PendingTask task, TimelineStatus tstatus, long last_pdl_lag) {

		// If the current task is not a delayed task, delete any delayed task for the current timeline entry

		switch (task.get_opcode()) {
		case OPCODE_GEN_FORECAST:
		case OPCODE_GEN_PDL_REPORT:
		case OPCODE_GEN_EXPIRE:
			break;

		default:
			//delete_delayed_timeline_tasks (task.get_event_id());
			delete_delayed_timeline_tasks (tstatus.event_id);
			break;
		}

		// If there is a catalog snapshot available, write it to the database

		CompactEqkRupList catalog_aftershocks = tstatus.get_catalog_snapshot();

		if (catalog_aftershocks != null) {

			// Write catalog snapshot to database

			CatalogSnapshot.submit_catalog_shapshot (
				task.get_record_key(),							// key
				tstatus.event_id,								// event_id
				tstatus.forecast_results.catalog_start_time,	// start_time
				tstatus.forecast_results.catalog_end_time,		// end_time
				catalog_aftershocks);							// rupture_list

			// Display message
		
			sg.task_disp.display_taskinfo ("TASK-INFO: Catalog snapshot saved:\n"
				+ "event_id = " + tstatus.event_id + "\n"
				+ "catalog_eqk_count = " + tstatus.forecast_results.catalog_eqk_count);

			sg.log_sup.report_catalog_saved (tstatus.event_id, tstatus.forecast_results.catalog_eqk_count);

		}

		// Write timeline entry

		TimelineEntry.submit_timeline_entry (
			task.get_record_key(),				// key
			tstatus.action_time,				// action_time
			tstatus.event_id,					// event_id
			tstatus.comcat_ids,					// comcat_ids
			tstatus.actcode,					// actcode
			tstatus.marshal_timeline());		// details

		// Display message
		
		sg.task_disp.display_taskinfo ("TASK-INFO: Timeline appended:\n"
			+ "event_id = " + tstatus.event_id + "\n"
			+ "actcode = " + tstatus.get_actcode_as_string () + "\n"
			+ "action_time = " + tstatus.action_time + "\n"
			+ "fc_status = " + tstatus.get_fc_status_as_string () + "\n"
			+ "pdl_status = " + tstatus.get_pdl_status_as_string () + "\n"
			+ "fc_result = " + tstatus.get_fc_result_as_string ());

		sg.log_sup.report_timeline_appended (tstatus);

		// Issue any new delayed command that is needed

		next_auto_timeline (tstatus, last_pdl_lag);

		return;
	}




	// Determine the next forecast lag for a timeline entry.
	// Parameters:
	//  tstatus = Timeline status.
	// Returns -1L if there is no next forecast lag.
	// If there is a next forecast lag, the return value is
	// positive and a multiple of the lag unit.
	// Note: Forecast lag is relative to the mainshock origin time.
	// Note: This does not consider the possibility of a PDL report.

	public long get_next_forecast_lag (TimelineStatus tstatus) {

		long next_forecast_lag = -1L;

		// If the timeline state requests a forecast ...

		if (tstatus.is_forecast_state()) {

			long min_lag;
			long min_extra_lag;

			// Minimum lag is 0L if this is the first forecast,
			// otherwise it is a minimum spacing after the last lag

			if (tstatus.last_forecast_lag < 0L) {
				min_lag = 0L;
				min_extra_lag = 0L;
			} else {
				min_lag = tstatus.last_forecast_lag + sg.task_disp.get_action_config().get_forecast_min_gap();
				min_extra_lag = tstatus.last_forecast_lag;
			}

			// Get next forecast lag from configured schedule

			next_forecast_lag = sg.task_disp.get_action_config().get_next_forecast_lag (min_lag, tstatus.analyst_options.max_forecast_lag);

			// If no next forecast lag on the configured schedule ...

			if (next_forecast_lag < 0L) {
			
				// Use the requested extra forecast lag, if any

				if (tstatus.analyst_options.extra_forecast_lag >= 0L) {

					// Make sure the value is a multiple of the lag unit, and greater than the last lag

					next_forecast_lag = 
						sg.task_disp.get_action_config().floor_unit_lag (tstatus.analyst_options.extra_forecast_lag, min_extra_lag);
				} else {
					next_forecast_lag = -1L;
				}
			}

			// Otherwise, we have a forecast lag from the schedule ...

			else {
			
				// If there is a requested extra forecast lag ...

				if (tstatus.analyst_options.extra_forecast_lag >= 0L) {

					// Use the smaller of the scheduled and extra lags

					next_forecast_lag = Math.min (next_forecast_lag,
						sg.task_disp.get_action_config().floor_unit_lag (tstatus.analyst_options.extra_forecast_lag, min_extra_lag));
				}
			}
		}

		// Return next lag

		return next_forecast_lag;
	}




	// Determine the next PDL report stage for a timeline entry.
	// Parameters:
	//  tstatus = Timeline status.
	//  next_forecast_lag = Scheduled time of the next forecast, or -1L if none.
	//  last_pdl_lag = Lag for the last PDL report attempt, or -1L if none.
	//	base_pdl_time = Start time for the PDL report sequence.
	// Returns -1L if there is no next PDL report lag.
	// If there is a next PDL report lag, the return value is
	// non-negative, a multiple of the lag unit, and greather than last_pdl_lag.
	// Note: PDL report lag is relative to the start time of the PDL report sequence.
	// Note: If last_pdl_lag == -1L and there is a next PDL report lag, then the return value is 0L.

	public long get_next_pdl_lag (TimelineStatus tstatus, long next_forecast_lag, long last_pdl_lag, long base_pdl_time) {

		long next_pdl_lag = -1L;

		// If the timeline state requests a PDL retry ...

		if (tstatus.is_pdl_retry_state()) {

			// Get next PDL lag from the schedule, must be > last_pdl_lag

			if (last_pdl_lag >= 0L) {
				next_pdl_lag = sg.task_disp.get_action_config().get_next_pdl_report_retry_lag (last_pdl_lag + 1L);
			} else {
				next_pdl_lag = 0L;
			}

			// If there is another lag in the schedule ...

			if (next_pdl_lag >= 0L) {

				// This will be the PDL time ceiling (the retry must occur before this time)

				long pdl_time_ceiling = Long.MAX_VALUE;

				// If there is a next forecast, limit to a time shortly before the forecast

				if (next_forecast_lag >= 0L) {
					pdl_time_ceiling = Math.min (pdl_time_ceiling,
						tstatus.forecast_mainshock.mainshock_time + next_forecast_lag - sg.task_disp.get_action_config().get_forecast_min_gap());
				}

				// If there is a previous forecast (should always be), limit to a maximum time after it
				// Note: Originally the "else" did not appear on the following line.  With "else" the
				// staleness test is applied only after the last forecast.  Wihtout "else" the staleness
				// test is applied to every forecast.

				else if (tstatus.last_forecast_lag >= 0L) {
					pdl_time_ceiling = Math.min (pdl_time_ceiling,
						tstatus.forecast_mainshock.mainshock_time + tstatus.last_forecast_lag + sg.task_disp.get_action_config().get_forecast_max_delay());
				}

				// Kill the PDL retry if it would not occur before the time ceiling
				// (The max below is the projected execution time of the PDL retry)

				if (Math.max (base_pdl_time + next_pdl_lag, sg.task_disp.get_time()) >= pdl_time_ceiling) {
					next_pdl_lag = -1L;
				}

			} else {
				next_pdl_lag = -1L;
			}
		}

		// Return next stage

		return next_pdl_lag;
	}




	// Post the next automatic database command for a timeline entry.
	// Parameters:
	//  tstatus = Timeline status.
	//  last_pdl_lag = Lag for the last PDL report attempt, or -1L if none.  Defaults to -1L if omitted.
	// Returns true if a task was posted.
	// Note: Typically last_pdl_lag == 0L if the caller has already attempted a PDL report, -1L otherwise.

	public boolean next_auto_timeline (TimelineStatus tstatus) {

		return next_auto_timeline (tstatus, -1L);
	}


	public boolean next_auto_timeline (TimelineStatus tstatus, long last_pdl_lag) {

		// Get the next forecast lag, or -1 if none

		long next_forecast_lag = get_next_forecast_lag (tstatus);

		// Get the next PDL report lag, or -1 if none, assuming the report begins now

		long next_pdl_lag = get_next_pdl_lag (tstatus, next_forecast_lag, last_pdl_lag, sg.task_disp.get_time());

		// If a PDL report is desired, submit the task

		if (next_pdl_lag >= 0L) {
		
			OpGeneratePDLReport pdl_payload = new OpGeneratePDLReport();
			pdl_payload.setup (tstatus.action_time, tstatus.last_forecast_lag, sg.task_disp.get_time());

			PendingTask.submit_task (
				tstatus.event_id,											// event id
				sg.task_disp.get_time() + next_pdl_lag,						// sched_time
				sg.task_disp.get_time(),									// submit_time
				SUBID_AAFS,													// submit_id
				OPCODE_GEN_PDL_REPORT,										// opcode
				sg.task_disp.get_action_config().lag_to_int (next_pdl_lag),	// stage
				pdl_payload.marshal_task());								// details

			sg.log_sup.report_pdl_report_request (tstatus.event_id, sg.task_disp.get_time() + next_pdl_lag);

			return true;
		}

		// If a forecast is desired, submit the task

		if (next_forecast_lag >= 0L) {
		
			OpGenerateForecast forecast_payload = new OpGenerateForecast();
			forecast_payload.setup (tstatus.action_time, tstatus.last_forecast_lag, next_forecast_lag);

			PendingTask.submit_task (
				tstatus.event_id,												// event id
				tstatus.forecast_mainshock.mainshock_time + next_forecast_lag 
				+ sg.task_disp.get_action_config().get_comcat_clock_skew()
				+ sg.task_disp.get_action_config().get_comcat_origin_skew(),	// sched_time
				sg.task_disp.get_time(),										// submit_time
				SUBID_AAFS,														// submit_id
				OPCODE_GEN_FORECAST,											// opcode
				0,																// stage
				forecast_payload.marshal_task());								// details

			sg.log_sup.report_forecast_request (tstatus.event_id,
				tstatus.forecast_mainshock.mainshock_time + next_forecast_lag 
				+ sg.task_disp.get_action_config().get_comcat_clock_skew()
				+ sg.task_disp.get_action_config().get_comcat_origin_skew());

			return true;
		}

		// If timeline state requests an action, submit an expire command

		if (tstatus.is_forecast_state() || tstatus.is_pdl_retry_state()) {
		
			OpGenerateExpire expire_payload = new OpGenerateExpire();
			expire_payload.setup (tstatus.action_time, tstatus.last_forecast_lag);

			PendingTask.submit_task (
				tstatus.event_id,										// event id
				sg.task_sup.get_prompt_exec_time(),						// sched_time
				sg.task_disp.get_time(),								// submit_time
				SUBID_AAFS,												// submit_id
				OPCODE_GEN_EXPIRE,										// opcode
				0,														// stage
				expire_payload.marshal_task());							// details

			sg.log_sup.report_expire_request (tstatus.event_id);

			return true;
		}

		// No command required

		return false;
	}




	// Set up Comcat retry for an intake command, while processing a Comcat exception.
	// Return values:
	//  RESCODE_STAGE_COMCAT_RETRY = The task is being staged for retry.
	//  RESCODE_INTAKE_COMCAT_FAIL = Retries exhausted, the task has failed.

	public int intake_setup_comcat_retry (PendingTask task, ComcatException e) {

		sg.log_sup.report_comcat_exception (task.get_event_id(), e);

		// For PDL intake only, delete any other PDL intake commands for this event, so we don't have multiple retries going on

		if (task.get_opcode() == OPCODE_INTAKE_PDL) {
			sg.task_sup.delete_all_waiting_tasks (task.get_event_id(), OPCODE_INTAKE_PDL);
		}

		// Get the next ComCat retry lag

		long next_comcat_intake_lag = sg.task_disp.get_action_config().get_next_comcat_intake_lag (
										sg.task_disp.get_action_config().int_to_lag (task.get_stage()) + 1L );

		// If there is another retry, stage the task

		if (next_comcat_intake_lag >= 0L) {
			sg.task_disp.set_taskres_stage (task.get_sched_time() + next_comcat_intake_lag,
								sg.task_disp.get_action_config().lag_to_int (next_comcat_intake_lag));
			return RESCODE_STAGE_COMCAT_RETRY;
		}

		// Retries exhausted, display the error and log the task
		
		sg.task_disp.set_display_taskres_log ("TASK-ERR: Intake failed due to ComCat failure:\n"
			+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
			+ "event_id = " + task.get_event_id() + "\n"
			+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));

		return RESCODE_INTAKE_COMCAT_FAIL;
	}




	// Convert event ID to timeline ID for an intake command.
	// Returns values:
	//  RESCODE_SUCCESS = The task already contains a timeline ID.
	//  RESCODE_STAGE_TIMELINE_ID or RESCODE_STAGE_COMCAT_RETRY = The task is being staged
	//    for retry, either to start over with a timeline ID in place of an event ID, or
	//    to retry a failed Comcat operation.
	//  RESCODE_INTAKE_COMCAT_FAIL = Comcat retries exhausted, the command has failed.
	//  RESCODE_ALIAS_EVENT_NOT_IN_COMCAT = The event ID is not in Comcat.

	public int intake_event_id_to_timeline_id (PendingTask task) {

		// If the task already contains a timeline ID, just return

		if (sg.alias_sup.is_timeline_id (task.get_event_id())) {
			return RESCODE_SUCCESS;
		}

		// Get mainshock parameters for an event ID

		ForecastMainshock fcmain = new ForecastMainshock();
		int retval;

		try {
			retval = sg.alias_sup.get_mainshock_for_event_id (task.get_event_id(), fcmain);
		}

		// Handle Comcat exception

		catch (ComcatException e) {
			return intake_setup_comcat_retry (task, e);
		}

		// If event not in Comcat, then return

		if (retval == RESCODE_ALIAS_EVENT_NOT_IN_COMCAT) {
			return RESCODE_ALIAS_EVENT_NOT_IN_COMCAT;
		}

		// If the event ID has not been seen before, create the alias timeline

		if (retval == RESCODE_ALIAS_NEW_EVENT) {
			sg.alias_sup.write_mainshock_to_new_timeline (fcmain);
		}

		// Stage the task, using the timeline ID in place of the event ID, for immediate execution

		sg.task_disp.set_taskres_stage (sg.task_sup.get_prompt_exec_time(),		// could use EXEC_TIME_MIN_WAITING
										task.get_stage(),
										fcmain.timeline_id);
		return RESCODE_STAGE_TIMELINE_ID;
	}




	//----- Construction -----


	// Default constructor.

	public TimelineSupport () {}

}
