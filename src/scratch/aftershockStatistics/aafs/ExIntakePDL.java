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

		int etres = sg.timeline_sup.intake_event_id_to_timeline_id (task);
		if (etres != RESCODE_SUCCESS) {
			return etres;
		}

		//--- Get payload and timeline status

		OpIntakePDL payload = new OpIntakePDL();
		TimelineStatus tstatus = new TimelineStatus();

		int rescode = sg.timeline_sup.open_timeline (task, tstatus, payload);

		switch (rescode) {

		case RESCODE_TIMELINE_EXISTS:
			return RESCODE_DELETE;		// Just delete, so that log is not flooded with PDL notifications

		case RESCODE_TIMELINE_NOT_FOUND:
			break;

		default:
			return rescode;
		}

		//--- Mainshock data

		// Get mainshock parameters

		ForecastParameters forecast_params = new ForecastParameters();

		try {
			sg.alias_sup.get_mainshock_for_timeline_id_ex (task.get_event_id(), forecast_params);
		}

		// An exception here triggers a ComCat retry

		catch (ComcatException e) {
			return sg.timeline_sup.intake_setup_comcat_retry (task, e);
		}

		//--- Intake check

		// Search intake regions, using the minimum magnitude criterion

		IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
			forecast_params.mainshock_lat, forecast_params.mainshock_lon, forecast_params.mainshock_mag);

		if (intake_region == null) {

			// Didn't pass, check the original PDL values using the minimum magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
				payload.mainshock_lat, payload.mainshock_lon, payload.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE;		// Just delete, so that log is not flooded with PDL notifications
			}

			// Now search intake regions, using the intake magnitude criterion

			intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_intake_mag (
				forecast_params.mainshock_lat, forecast_params.mainshock_lon, forecast_params.mainshock_mag);

			// If none, then drop the event

			if (intake_region == null) {
				return RESCODE_DELETE;		// Just delete, so that log is not flooded with PDL notifications
			}
		}

		//--- Final steps

		// Set track state
			
		tstatus.set_state_track (
			sg.task_disp.get_time(),
			sg.task_disp.get_action_config(),
			task.get_event_id(),
			forecast_params,
			TimelineStatus.FCORIG_PDL,
			TimelineStatus.FCSTAT_ACTIVE_INTAKE);

		// If the command contains analyst data, save it

		if (payload.f_has_analyst) {
			tstatus.set_analyst_data  (
				payload.analyst_id,
				payload.analyst_remark,
				payload.analyst_time,
				payload.analyst_params,
				payload.extra_forecast_lag);
		}

		// Write the new timeline entry

		sg.timeline_sup.append_timeline (task, tstatus);

		// Log the task

		return RESCODE_SUCCESS;
	}




	//----- Construction -----


	// Default constructor.

	public ExIntakePDL () {}

}
