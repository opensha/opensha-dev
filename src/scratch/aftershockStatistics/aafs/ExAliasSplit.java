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
 * Execute task: Notify alias timeline split.
 * Author: Michael Barall 07/19/2018.
 */
public class ExAliasSplit extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_alias_split (task);
	}




	// Notify alias timeline split.

	private int exec_alias_split (PendingTask task) {

		// Convert event ID to timeline ID if needed (should never happen for this operation)

		int etres = sg.timeline_sup.intake_event_id_to_timeline_id (task);
		if (etres != RESCODE_SUCCESS) {
			return etres;
		}

		//--- Get payload and timeline status

		OpAliasSplit payload = new OpAliasSplit();
		TimelineStatus tstatus = new TimelineStatus();

		int rescode = sg.timeline_sup.open_timeline (task, tstatus, payload);

		switch (rescode) {

		case RESCODE_TIMELINE_EXISTS:
			return RESCODE_TIMELINE_EXISTS;		// Log task as not referring to a new timeline

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
			return RESCODE_INTAKE_FILTERED;		// Log task as not passing intake filter
		}

		//--- Final steps

		// Set track state
			
		tstatus.set_state_track (
			sg.task_disp.get_time(),
			sg.task_disp.get_action_config(),
			task.get_event_id(),
			forecast_params,
			TimelineStatus.FCORIG_SPLIT,
			TimelineStatus.FCSTAT_ACTIVE_INTAKE);

		// Write the new timeline entry

		sg.timeline_sup.append_timeline (task, tstatus);

		// Log the task

		return RESCODE_SUCCESS;
	}




	//----- Construction -----


	// Default constructor.

	public ExAliasSplit () {}

}
