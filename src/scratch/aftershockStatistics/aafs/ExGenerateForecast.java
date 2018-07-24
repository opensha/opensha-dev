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
 * Execute task: Generate forecast.
 * Author: Michael Barall 06/25/2018.
 */
public class ExGenerateForecast extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_gen_forecast (task);
	}




	// Generate forecast.

	private int exec_gen_forecast (PendingTask task) {

		//--- Get payload and timeline status

		OpGenerateForecast payload = new OpGenerateForecast();
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

		// Get the expected forecast lag

		long next_forecast_lag = sg.timeline_sup.get_next_forecast_lag (tstatus);

		// Check that timeline is active

		if (next_forecast_lag < 0L) {
		
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry is not active:\n"
				+ "event_id = " + task.get_event_id() + "\n"
				+ "tstatus.fc_status = " + tstatus.get_fc_status_as_string());

			sg.timeline_sup.next_auto_timeline (tstatus);
			return RESCODE_TIMELINE_NOT_ACTIVE;
		}

		// Check state matches the command

		if (!( payload.action_time == tstatus.action_time
			&& payload.next_forecast_lag == next_forecast_lag
			&& payload.last_forecast_lag == tstatus.last_forecast_lag )) {
		
			sg.task_disp.set_display_taskres_log ("TASK-ERR: Timeline entry state does not match task:\n"
				+ "event_id = " + task.get_event_id() + "\n"
				+ "payload.action_time = " + payload.action_time + "\n"
				+ "tstatus.action_time = " + tstatus.action_time + "\n"
				+ "payload.next_forecast_lag = " + payload.next_forecast_lag + "\n"
				+ "next_forecast_lag = " + next_forecast_lag + "\n"
				+ "payload.last_forecast_lag = " + payload.last_forecast_lag + "\n"
				+ "tstatus.last_forecast_lag = " + tstatus.last_forecast_lag);

			sg.timeline_sup.next_auto_timeline (tstatus);
			return RESCODE_TIMELINE_TASK_MISMATCH;
		}

		//--- Forecast

		// Get mainshock parameters

		ForecastParameters forecast_params = new ForecastParameters();

		try {
			sg.alias_sup.get_mainshock_for_timeline_id_ex (task.get_event_id(), forecast_params);
		}

		// An exception here triggers a ComCat retry

		catch (ComcatException e) {
			return sg.timeline_sup.process_timeline_comcat_retry (task, tstatus, e);
		}

		// Check if it's too soon to do this forecast (might happen if mainshock origin time has changed)

		if (forecast_params.mainshock_time + next_forecast_lag
			+ sg.task_disp.get_action_config().get_comcat_clock_skew() > sg.task_disp.get_time()) {
		
			// Stage the task

			sg.task_disp.set_taskres_stage (forecast_params.mainshock_time + next_forecast_lag
								+ sg.task_disp.get_action_config().get_comcat_clock_skew()
								+ sg.task_disp.get_action_config().get_comcat_origin_skew(),
								task.get_stage());
			return RESCODE_STAGE;
		}

		// If intake need to be re-checked ...

		if (tstatus.is_intake_state()) {

			// Search intake regions

			IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
				forecast_params.mainshock_lat, forecast_params.mainshock_lon, forecast_params.mainshock_mag);

			// If none, then withdraw the timeline

			if (intake_region == null) {
		
				sg.task_disp.set_display_taskres_log ("TASK-INFO: Timeline entry withdrawn:\n"
					+ "event_id = " + task.get_event_id() + "\n"
					+ "forecast_params.mainshock_lat = " + forecast_params.mainshock_lat + "\n"
					+ "forecast_params.mainshock_lon = " + forecast_params.mainshock_lon + "\n"
					+ "forecast_params.mainshock_mag = " + forecast_params.mainshock_mag);

				tstatus.set_state_withdrawn (sg.task_disp.get_time(), forecast_params.mainshock_time);
				sg.timeline_sup.append_timeline (task, tstatus);

				return RESCODE_TIMELINE_WITHDRAW;
			}
		}

		// Fetch parameters, part 2 (model and search parameters), and calculate results

		ForecastResults forecast_results = new ForecastResults();

		try {
			forecast_params.fetch_all_2 (next_forecast_lag, tstatus.analyst_params);

			long advisory_lag;

			if (next_forecast_lag >= sg.task_disp.get_action_config().get_advisory_dur_year()) {
				advisory_lag = ForecastResults.ADVISORY_LAG_YEAR;
			} else if (next_forecast_lag >= sg.task_disp.get_action_config().get_advisory_dur_month()) {
				advisory_lag = ForecastResults.ADVISORY_LAG_MONTH;
			} else if (next_forecast_lag >= sg.task_disp.get_action_config().get_advisory_dur_week()) {
				advisory_lag = ForecastResults.ADVISORY_LAG_WEEK;
			} else {
				advisory_lag = ForecastResults.ADVISORY_LAG_DAY;
			}

			forecast_results.calc_all (
				forecast_params.mainshock_time + next_forecast_lag,
				advisory_lag,
				forecast_params,
				next_forecast_lag >= sg.task_disp.get_action_config().get_seq_spec_min_lag());
		}

		// An exception here triggers a ComCat retry

		catch (Exception e) {
			return sg.timeline_sup.process_timeline_comcat_retry (task, tstatus, e);
		}

		// Select report for PDL, if any

		forecast_results.pick_pdl_model();

		// If we have an earthquake catalog ...

		if (forecast_results.catalog_result_avail) {

			// If there is an aftershock with larger magnitude than the mainshock ...

			if (forecast_results.catalog_eqk_count > 0 && forecast_results.catalog_max_mag > forecast_params.mainshock_mag) {
			
				// Set timeline to foreshock state

				tstatus.set_state_foreshock (sg.task_disp.get_time(), forecast_results.catalog_max_event_id);

				// Display message

				sg.task_disp.set_display_taskres_log ("TASK-INFO: Foreshock detected:\n"
					+ "event_id = " + tstatus.event_id + "\n"
					+ "mainshock_mag = " + forecast_params.mainshock_mag + "\n"
					+ "catalog_max_event_id = " + forecast_results.catalog_max_event_id + "\n"
					+ "catalog_max_mag = " + forecast_results.catalog_max_mag);

				// Write the timeline entry

				sg.timeline_sup.append_timeline (task, tstatus);

				// Log the task

				return RESCODE_TIMELINE_FORESHOCK;
			}
		}

		// Insert forecast into timeline status

		tstatus.set_state_forecast (sg.task_disp.get_time(), sg.task_disp.get_action_config(), forecast_params, forecast_results);

		// Get the next forecast lag, or -1 if none

		long new_next_forecast_lag = sg.timeline_sup.get_next_forecast_lag (tstatus);

		// If no next forecast, mark the timeline expired

		if (new_next_forecast_lag < 0L) {
			tstatus.set_fc_status (TimelineStatus.FCSTAT_STOP_EXPIRED);
		}

		//--- PDL report

		// If a PDL report is requested based on state ...

		if (tstatus.is_pdl_retry_state()) {

			// Get the next PDL report lag, or -1 if none, assuming the report begins now

			long new_next_pdl_lag = sg.timeline_sup.get_next_pdl_lag (tstatus, new_next_forecast_lag, -1L, sg.task_disp.get_time());

			// Check for secondary status

			if (!( sg.pdl_sup.is_pdl_primary() )) {
				tstatus.set_pdl_status (TimelineStatus.PDLSTAT_SECONDARY);
			}

			// Otherwise, if no PDL report is requested based on timing, mark it bypassed

			else if (new_next_pdl_lag < 0L) {
				tstatus.set_pdl_status (TimelineStatus.PDLSTAT_BYPASSED);
			}

			// Otherwise, we need to send a PDL report ...

			else {
			
				// Attempt to send the report

				try {
					sg.pdl_sup.send_pdl_report (tstatus);
					tstatus.set_pdl_status (TimelineStatus.PDLSTAT_SUCCESS);
					new_next_pdl_lag = -1L;
				}

				// Exception here means PDL report did not succeed

				catch (Exception e) {

					// Get time of PDL retry

					tstatus.set_pdl_status (TimelineStatus.PDLSTAT_PENDING);	// in case it was changed in the try block
					new_next_pdl_lag = sg.timeline_sup.get_next_pdl_lag (tstatus, new_next_forecast_lag, 0L, sg.task_disp.get_time());

					// If no retry, report the failure now

					if (new_next_pdl_lag < 0L) {
						tstatus.set_pdl_status (TimelineStatus.PDLSTAT_FAILURE);
		
						sg.task_disp.display_taskinfo ("TASK-ERR: Unable to send forecast report to PDL:\n"
							+ "event_id = " + tstatus.event_id + "\n"
							+ "last_forecast_lag = " + tstatus.last_forecast_lag + "\n"
							+ "Stack trace:\n" + SimpleUtils.getStackTraceAsString(e));
					}
				}
			}
		}

		//--- Final steps

		// Write the new timeline entry

		sg.timeline_sup.append_timeline (task, tstatus, 0L);

		// Log the task

		return RESCODE_SUCCESS;
	}




	//----- Construction -----


	// Default constructor.

	public ExGenerateForecast () {}

}
