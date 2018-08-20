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
 * Execute task: Shutdown.
 * Author: Michael Barall 06/25/2018.
 */
public class ExIntakePDL extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_intake_pdl (task);
	}




	// Intake an event for PDL.

	private int exec_intake_pdl (PendingTask task) {

		// Convert event ID to timeline ID if needed

		//int etres = sg.timeline_sup.intake_event_id_to_timeline_id (task);
		int etres = intake_pdl_event_id_to_timeline_id (task);
		if (etres != RESCODE_SUCCESS) {
			return etres;
		}

		//--- Get payload and timeline status

		OpIntakePDL payload = new OpIntakePDL();
		TimelineStatus tstatus = new TimelineStatus();

		int rescode = sg.timeline_sup.open_timeline (task, tstatus, payload);

		switch (rescode) {

		case RESCODE_TIMELINE_EXISTS:
			return RESCODE_DELETE_TIMELINE_EXISTS;		// Just delete, so that log is not flooded with PDL notifications

		case RESCODE_TIMELINE_NOT_FOUND:
			break;

		default:
			return rescode;
		}

		//--- Mainshock data

		// Get mainshock parameters

		ForecastMainshock fcmain = new ForecastMainshock();
		int retval;

		try {
			retval = sg.alias_sup.get_mainshock_for_timeline_id (task.get_event_id(), fcmain);
		}

		// An exception here triggers a ComCat retry

		catch (ComcatException e) {
			return sg.timeline_sup.intake_setup_comcat_retry (task, e);
		}

		// If timeline is not found or stopped, then drop the event

		if (retval != RESCODE_SUCCESS) {
			return RESCODE_DELETE_TIMELINE_NO_ALIAS;		// Just delete, so that log is not flooded with PDL notifications
		}

		//--- Intake check

		// Search intake regions, using the minimum magnitude criterion

		IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
			fcmain.mainshock_lat, fcmain.mainshock_lon, fcmain.mainshock_mag);

		if (intake_region == null) {

			// Didn't pass, check the original PDL values using the minimum magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
				payload.mainshock_lat, payload.mainshock_lon, payload.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with PDL notifications
			}

			// Now search intake regions, using the intake magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_intake_mag (
				fcmain.mainshock_lat, fcmain.mainshock_lon, fcmain.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with PDL notifications
			}
		}

		//--- Final steps

		// Set track state
			
		tstatus.set_state_track (
			sg.task_disp.get_time(),
			sg.task_disp.get_action_config(),
			task.get_event_id(),
			fcmain,
			TimelineStatus.FCORIG_PDL,
			TimelineStatus.FCSTAT_ACTIVE_INTAKE);

		// If the command contains analyst data, save it

		if (payload.analyst_options != null) {
			tstatus.set_analyst_data (payload.analyst_options);
		}

		// Write the new timeline entry

		sg.timeline_sup.append_timeline (task, tstatus);

		// Log the task

		return RESCODE_SUCCESS;
	}




	// Convert event ID to timeline ID for an intake PDL command.
	// Returns values:
	//  RESCODE_SUCCESS = The task already contains a timeline ID.
	//  RESCODE_STAGE_TIMELINE_ID or RESCODE_STAGE_COMCAT_RETRY = The task is being staged
	//    for retry, either to start over with a timeline ID in place of an event ID, or
	//    to retry a failed Comcat operation.
	//  RESCODE_INTAKE_COMCAT_FAIL = Comcat retries exhausted, the command has failed.
	//  RESCODE_DELETE_INTAKE_FILTERED = The task should be deleted without being logged.
	//  RESCODE_TASK_CORRUPT = The task payload is corrupted.
	// Note: This function is the same as TimelineSupport.intake_event_id_to_timeline_id,
	// except that it checks the intake filter and timeline existence.

	private int intake_pdl_event_id_to_timeline_id (PendingTask task) {

		// If the task already contains a timeline ID, just return

		if (sg.alias_sup.is_timeline_id (task.get_event_id())) {
			return RESCODE_SUCCESS;
		}

		//--- Get payload

		OpIntakePDL payload = new OpIntakePDL();

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		//  //--- Test for timeline existence
		//  //
		//  // Note: The test is moved below.  Putting it here would avoid an extra Comcat
		//  // query when the timeline exists.  But a PDL notification might coincide with
		//  // an alias reassignment, and so we do the Comcat query to check aliases.
		//  
		//  // Try to identify the timeline for this event
		//  
		//  String timeline_id = sg.alias_sup.get_timeline_id_for_primary_id (potential_event_id);
		//  
		//  // Get the corresponding timeline entry, or null if none
		//  
		//  TimelineEntry tentry = null;
		//  
		//  if (timeline_id != null) {
		//  	tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, timeline_id, null, null);
		//  }
		//  
		//  // If timeline found, then drop the event
		//  
		//  if (tentry != null) {
		//  	return RESCODE_DELETE_TIMELINE_EXISTS;		// Just delete, so that log is not flooded with PDL notifications
		//  }

		//--- Mainshock data

		// Get mainshock parameters for an event ID

		ForecastMainshock fcmain = new ForecastMainshock();
		int retval;

		try {
			retval = sg.alias_sup.get_mainshock_for_event_id (task.get_event_id(), fcmain);
		}

		// Handle Comcat exception

		catch (ComcatException e) {
			return sg.timeline_sup.intake_setup_comcat_retry (task, e);
		}

		// If event not in Comcat, then drop the event

		if (retval == RESCODE_ALIAS_EVENT_NOT_IN_COMCAT) {
			return RESCODE_DELETE_NOT_IN_COMCAT;		// Just delete, so that log is not flooded with notifications
		}

		//--- Intake check

		// Search intake regions, using the minimum magnitude criterion

		IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
			fcmain.mainshock_lat, fcmain.mainshock_lon, fcmain.mainshock_mag);

		if (intake_region == null) {

			// Didn't pass, check the original PDL values using the minimum magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
				payload.mainshock_lat, payload.mainshock_lon, payload.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with PDL notifications
			}

			// Now search intake regions, using the intake magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_intake_mag (
				fcmain.mainshock_lat, fcmain.mainshock_lon, fcmain.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with PDL notifications
			}
		}

		//--- Final steps (part 1)

		// If the event ID has not been seen before, create the alias timeline

		if (retval == RESCODE_ALIAS_NEW_EVENT) {
			sg.alias_sup.write_mainshock_to_new_timeline (fcmain);
		}

		//--- Test for timeline existence
		//
		// This test could be omitted, and timeline existence would be checked during
		// the staged task.  But doing it here avoids the Comcat query in the staged task
		// when the timeline exists.

		// Get the corresponding timeline entry, or null if none

		TimelineEntry tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, fcmain.timeline_id, null, null);

		// If timeline found, then drop the event

		if (tentry != null) {
			return RESCODE_DELETE_TIMELINE_EXISTS;		// Just delete, so that log is not flooded with PDL notifications
		}

		//--- Final steps (part 2)

		// Stage the task, using the timeline ID in place of the event ID, for immediate execution

		sg.task_disp.set_taskres_stage (sg.task_sup.get_prompt_exec_time(),		// could use EXEC_TIME_MIN_WAITING
										task.get_stage(),
										fcmain.timeline_id);
		return RESCODE_STAGE_TIMELINE_ID;
	}




	//----- Construction -----


	// Default constructor.

	public ExIntakePDL () {}

}
