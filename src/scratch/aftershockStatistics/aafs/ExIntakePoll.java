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

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.ComcatException;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

/**
 * Execute task: Intake event discovered by polling.
 * Author: Michael Barall 08/05/2018.
 */
public class ExIntakePoll extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_intake_poll (task);
	}




	// Intake an event for poll.

	private int exec_intake_poll (PendingTask task) {

		// Convert delayed command to non-delayed command if needed

		int dcres = intake_poll_delayed (task);
		if (dcres != RESCODE_SUCCESS) {
			return dcres;
		}

		// Convert event ID to timeline ID if needed

		//int etres = sg.timeline_sup.intake_event_id_to_timeline_id (task);
		int etres = intake_poll_event_id_to_timeline_id (task);
		if (etres != RESCODE_SUCCESS) {
			return etres;
		}

		//--- Get payload and timeline status

		OpIntakePoll payload = new OpIntakePoll();
		TimelineStatus tstatus = new TimelineStatus();

		int rescode = sg.timeline_sup.open_timeline (task, tstatus, payload);

		switch (rescode) {

		case RESCODE_TIMELINE_EXISTS:

			// If in a state that permits re-activation ...

			if (tstatus.can_intake_poll_start()) {

				//--- Mainshock data

				// Get mainshock parameters

				ForecastMainshock fcmain2 = new ForecastMainshock();
				int retval2;

				try {
					retval2 = sg.alias_sup.get_mainshock_for_timeline_id (task.get_event_id(), fcmain2);
				}

				// An exception here triggers a ComCat retry

				catch (ComcatException e) {
					return sg.timeline_sup.intake_setup_comcat_retry (task, e);
				}

				// If timeline is not found or stopped, then drop the event

				if (retval2 != RESCODE_SUCCESS) {
					return RESCODE_DELETE_TIMELINE_NO_ALIAS;	// Just delete, so that log is not flooded with poll notifications
				}

				//--- Intake check

				// Search intake regions, using the minimum magnitude criterion

				IntakeSphRegion intake_region2 = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
					fcmain2.mainshock_lat, fcmain2.mainshock_lon, fcmain2.mainshock_mag);

				// If none, then drop the event

				if (intake_region2 == null) {
					return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with poll notifications
				}

				//--- Final steps

				// Set status update
			
				tstatus.set_state_status_update (
					sg.task_disp.get_time(),
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

			return RESCODE_DELETE_TIMELINE_BAD_STATE;		// Just delete, so that log is not flooded with PDL notifications

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
			return RESCODE_DELETE_TIMELINE_NO_ALIAS;		// Just delete, so that log is not flooded with poll notifications
		}

		//--- Intake check

		// Search intake regions, using the minimum magnitude criterion

		IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
			fcmain.mainshock_lat, fcmain.mainshock_lon, fcmain.mainshock_mag);

		// If none, then drop the event

		if (intake_region == null) {
			return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with poll notifications
		}

		//--- Final steps

		// Set track state
			
		tstatus.set_state_track (
			sg.task_disp.get_time(),
			sg.task_disp.get_action_config(),
			task.get_event_id(),
			fcmain,
			TimelineStatus.FCORIG_POLL,
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




	// Convert event ID to timeline ID for an intake poll command.
	// Returns values:
	//  RESCODE_SUCCESS = The task already contains a timeline ID.
	//  RESCODE_STAGE_TIMELINE_ID or RESCODE_STAGE_COMCAT_RETRY = The task is being staged
	//    for retry, either to start over with a timeline ID in place of an event ID, or
	//    to retry a failed Comcat operation.
	//  RESCODE_INTAKE_COMCAT_FAIL = Comcat retries exhausted, the command has failed.
	//  RESCODE_DELETE_NOT_IN_COMCAT or RESCODE_DELETE_INTAKE_FILTERED = The task
	//    should be deleted without being logged.
	//  RESCODE_TASK_CORRUPT = The task payload is corrupted.
	// Note: This function is the same as TimelineSupport.intake_event_id_to_timeline_id,
	// except that it checks the intake filter.

	private int intake_poll_event_id_to_timeline_id (PendingTask task) {

		// If the task already contains a timeline ID, just return

		if (sg.alias_sup.is_timeline_id (task.get_event_id())) {
			return RESCODE_SUCCESS;
		}

		//--- Get payload

		OpIntakePoll payload = new OpIntakePoll();

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

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

		// If none, then drop the event

		if (intake_region == null) {
			return RESCODE_DELETE_INTAKE_FILTERED;		// Just delete, so that log is not flooded with notifications
		}

		//--- Final steps

		// If the event ID has not been seen before, create the alias timeline

		if (retval == RESCODE_ALIAS_NEW_EVENT) {
			sg.alias_sup.write_mainshock_to_new_timeline (fcmain);
		}

		// Delete any other intake poll command with the same timeline ID

		sg.task_sup.delete_all_waiting_tasks (fcmain.timeline_id, OPCODE_INTAKE_POLL);

		// Stage the task, using the timeline ID in place of the event ID, for immediate execution

		sg.task_disp.set_taskres_stage (sg.task_sup.get_prompt_exec_time(),		// could use EXEC_TIME_MIN_WAITING
										task.get_stage(),
										fcmain.timeline_id);
		return RESCODE_STAGE_TIMELINE_ID;
	}




	// Replace event ID with delayed event ID for a delayed command.
	// Returns values:
	//  RESCODE_SUCCESS = The task is not a delayed command.
	//  RESCODE_STAGE_EVENT_ID = The task is being staged for retry, to
	//    convert a delayed command into a non-delayed command.
	//  RESCODE_TASK_CORRUPT = The task payload is corrupted.

	private int intake_poll_delayed (PendingTask task) {

		// If the task is not a delayed command, just return

		if (!( task.get_event_id().equals (EVID_POLL) )) {
			return RESCODE_SUCCESS;
		}

		//--- Get payload

		OpIntakePoll payload = new OpIntakePoll();

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		//--- Final steps

		// Delete any other intake poll command with the same event ID

		sg.task_sup.delete_all_waiting_tasks (payload.delayed_event_id, OPCODE_INTAKE_POLL);

		// Stage the task, using the delayed event ID in place of the event ID, for immediate execution

		sg.task_disp.set_taskres_stage (sg.task_sup.get_prompt_exec_time(),		// could use EXEC_TIME_MIN_WAITING
										task.get_stage(),
										payload.delayed_event_id);
		return RESCODE_STAGE_EVENT_ID;
	}




	//----- Construction -----


	// Default constructor.

	public ExIntakePoll () {}

}
