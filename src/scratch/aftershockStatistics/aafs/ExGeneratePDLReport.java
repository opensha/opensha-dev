package scratch.aftershockStatistics.aafs;

import java.util.List;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;

import scratch.aftershockStatistics.CompactEqkRupList;

/**
 * Execute task: Generate PDL report retry.
 * Author: Michael Barall 06/25/2018.
 */
public class ExGeneratePDLReport extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_gen_pdl_report (task);
	}




	// Generate PDL report retry.

	private int exec_gen_pdl_report (PendingTask task) {

		//--- Get payload and timeline status

		OpGeneratePDLReport payload = new OpGeneratePDLReport();
		TimelineStatus tstatus = new TimelineStatus();

		int rescode = sg.timeline_sup.open_timeline (task, tstatus, payload);

		switch (rescode) {

		case RESCODE_TIMELINE_EXISTS:
			break;

		case RESCODE_TIMELINE_NOT_FOUND:
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry not found:\n"
				+ "event_id = " + task.get_event_id());
			return rescode;

		default:
			return rescode;
		}

		//--- Timeline state check

		// Check that timeline is sending a PDL report

		if (!( tstatus.is_pdl_retry_state() )) {
		
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry is not sending a PDL report:\n"
				+ "event_id = " + task.get_event_id() + "\n"
				+ "tstatus.fc_status = " + tstatus.get_fc_status_as_string() + "\n"
				+ "tstatus.pdl_status = " + tstatus.get_pdl_status_as_string());

			sg.timeline_sup.next_auto_timeline (tstatus);
			return RESCODE_TIMELINE_NOT_PDL_PEND;
		}

		// Check state matches the command

		if (!( payload.action_time == tstatus.action_time
			&& payload.last_forecast_lag == tstatus.last_forecast_lag )) {
		
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry state does not match task:\n"
				+ "event_id = " + task.get_event_id() + "\n"
				+ "payload.action_time = " + payload.action_time + "\n"
				+ "tstatus.action_time = " + tstatus.action_time + "\n"
				+ "payload.last_forecast_lag = " + payload.last_forecast_lag + "\n"
				+ "tstatus.last_forecast_lag = " + tstatus.last_forecast_lag);

			sg.timeline_sup.next_auto_timeline (tstatus);
			return RESCODE_TIMELINE_TASK_MISMATCH;
		}

		//--- PDL report
			
		// Attempt to send the report

		try {
			sg.pdl_sup.send_pdl_report (tstatus);
		}

		// Exception here means PDL report did not succeed

		catch (Exception e) {

			sg.log_sup.report_pdl_send_exception (tstatus, e);

			// Get current PDL lag from the stage

			long new_next_pdl_lag = sg.task_disp.get_action_config().int_to_lag (task.get_stage());

			// Get the next forecast lag, or -1 if none

			long new_next_forecast_lag = sg.timeline_sup.get_next_forecast_lag (tstatus);

			// Get time of PDL retry

			new_next_pdl_lag = sg.timeline_sup.get_next_pdl_lag (tstatus, new_next_forecast_lag, new_next_pdl_lag, payload.base_pdl_time);

			// If there is another retry, stage the task

			if (new_next_pdl_lag >= 0L) {
				sg.task_disp.set_taskres_stage (payload.base_pdl_time + new_next_pdl_lag,
									sg.task_disp.get_action_config().lag_to_int (new_next_pdl_lag));

				return RESCODE_STAGE_PDL_RETRY;
			}

			// PDL report failed

			tstatus.set_state_pdl_update (sg.task_disp.get_time(), TimelineStatus.PDLSTAT_FAILURE);
		
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Unable to send forecast report to PDL:\n"
				+ "event_id = " + tstatus.event_id + "\n"
				+ "last_forecast_lag = " + tstatus.last_forecast_lag + "\n"
				+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));

			// Write the new timeline entry

			sg.timeline_sup.append_timeline (task, tstatus);

			// Log the task

			return RESCODE_TIMELINE_PDL_FAIL;
		}

		sg.log_sup.report_pdl_send_ok (tstatus);

		//--- Final steps

		// PDL report succeeded
			
		tstatus.set_state_pdl_update (sg.task_disp.get_time(), TimelineStatus.PDLSTAT_SUCCESS);

		// Write the new timeline entry

		sg.timeline_sup.append_timeline (task, tstatus);

		// Log the task

		return RESCODE_SUCCESS;
	}




	//----- Construction -----


	// Default constructor.

	public ExGeneratePDLReport () {}

}
