package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.Collections;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;
import scratch.aftershockStatistics.util.ObsEqkRupMaxTimeComparator;

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.ComcatException;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

/**
 * Execute task: Run Comcat poll.
 * Author: Michael Barall 08/05/2018.
 */
public class ExPollComcatRun extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_poll_comcat_run (task);
	}




	// Run Comcat poll.

	private int exec_poll_comcat_run (PendingTask task) {

		//--- Get payload

		OpPollComcatRun payload = new OpPollComcatRun();

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// Issue a new task to replace the corrupted one
			// Delay by a short poll period, to avoid going into a tight loop in case of repeating errors

			sg.poll_sup.kick_off_polling (sg.task_disp.get_action_config().get_poll_short_period());

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		//--- Polling operation

		// Find which poll to do

		int which_poll = sg.poll_sup.check_if_poll_time();

		// No poll, just stage the task

		if (which_poll == 0) {

			sg.task_disp.set_taskres_stage (sg.task_disp.get_time()
								+ sg.task_disp.get_action_config().get_poll_short_period(),
								task.get_stage());
			return RESCODE_STAGE_TOO_SOON;
		}

		// Get lookback depending on whether it is a short or long poll

		long poll_lookback;

		if (which_poll == 2) {
			poll_lookback = sg.task_disp.get_action_config().get_poll_long_lookback();
		} else {
			poll_lookback = sg.task_disp.get_action_config().get_poll_short_lookback();
		}

		// Create the accessor

		ComcatAccessor accessor = new ComcatAccessor();

		// Search the entire world, for minimum magnitude equal to the lowest in any intake region

		SphRegion search_region = SphRegion.makeWorld ();
		double min_mag = sg.task_disp.get_action_config().get_pdl_intake_region_min_min_mag();

		// Search within the lookback of the current time minus the clock skew

		long search_time_hi = sg.task_disp.get_time() - sg.task_disp.get_action_config().get_comcat_clock_skew();
		long search_time_lo = search_time_hi - poll_lookback;

		// This is the start of the search interval, for a short poll

		long search_time_short = search_time_hi - sg.task_disp.get_action_config().get_poll_short_lookback();

		// Call Comcat to get a list of potential events

		String exclude_id = null;

		double min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;
		double max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;

		boolean wrapLon = false;
		int limit_per_call = 0;
		int max_calls = 0;

		ObsEqkRupList potentials = null;

		try {
			potentials = accessor.fetchEventList (exclude_id,
					search_time_lo, search_time_hi,
					min_depth, max_depth,
					search_region, wrapLon,
					min_mag, limit_per_call, max_calls);
		}
		catch (Exception e) {
		
			// In case of Comcat failure, just delay to next polling time

			sg.task_disp.set_taskres_stage (sg.task_disp.get_time()
								+ sg.task_disp.get_action_config().get_poll_short_period(),
								task.get_stage());
			return RESCODE_STAGE_COMCAT_RETRY;
		}

		double poll_lookback_days = ((double)poll_lookback) / ((double)DURATION_DAY);
		System.out.println ("TASK-INFO: Comcat poll found " + potentials.size() + " potential events in " + String.format ("%.3f", poll_lookback_days) + " days");

		// Process potential events in temporal order, most recent first

		if (potentials.size() > 1) {
			Collections.sort (potentials, new ObsEqkRupMaxTimeComparator());
		}

		// Loop over potential events

		int count_no_timeline = 0;
		int count_withdrawn_timeline = 0;

		long total_delay = 0L;
		long short_delay = 0L;
		boolean f_seen_long = false;

		for (ObsEqkRupture potential : potentials) {

			// Get the potential parameters

			String potential_event_id = potential.getEventId();
			long potential_time = potential.getOriginTime();
			double potential_mag = potential.getMag();
			Location potential_hypo = potential.getHypocenterLocation();
			double potential_lat = potential_hypo.getLatitude();
			double potential_lon = potential_hypo.getLongitude();

			// Search intake regions, using the minimum magnitude criterion

			IntakeSphRegion intake_region = sg.task_disp.get_action_config().get_pdl_intake_region_for_min_mag (
				potential_lat, potential_lon, potential_mag);

			// If we passed the intake filter ...

			if (intake_region != null) {

				// Try to identify the timeline for this event

				String timeline_id = sg.alias_sup.get_timeline_id_for_primary_id (potential_event_id);

				// Get the corresponding timeline entry, or null if none

				TimelineEntry tentry = null;

				if (timeline_id != null) {
					tentry = TimelineEntry.get_recent_timeline_entry (0L, 0L, timeline_id, null, null);
				}

				// If no timeline found ...

				if (tentry == null) {

					// Submit a poll intake task for the event

					OpIntakePoll intake_payload = new OpIntakePoll();
					intake_payload.setup (potential_event_id);

					PendingTask.submit_task (
						EVID_POLL,									// event id
						sg.task_disp.get_time() + total_delay,		// sched_time
						sg.task_disp.get_time(),					// submit_time
						SUBID_AAFS,									// submit_id
						OPCODE_INTAKE_POLL,							// opcode
						0,											// stage
						intake_payload.marshal_task());				// details

					++count_no_timeline;

					// Adjust total and short delay time

					if (which_poll == 2 && potential_time < search_time_short) {
						f_seen_long = true;
					}
					if (f_seen_long) {
						total_delay += sg.task_disp.get_action_config().get_poll_long_intake_gap();
					} else {
						total_delay += sg.task_disp.get_action_config().get_poll_short_intake_gap();
						short_delay = total_delay;
					}
				}

				// Otherwise, found timeline entry ...

				else {

					// If the last action was one that could indicate a withdrawn timeline ...

					if (TimelineStatus.can_actcode_intake_poll_start (tentry.get_actcode())) {

						// Get the status for this timeline entry
					
						TimelineStatus tstatus = new TimelineStatus();
					
						try {
							tstatus.unmarshal_timeline (tentry);
						}
					
						// Invalid timeline entry
					
						catch (Exception e2) {
							tstatus = null;
						}
					
						// If we got the timeline status ...
					
						if (tstatus != null) {
					
							// If it is in a state that can be awakened by a poll ...

							if (tstatus.can_intake_poll_start()) {

								// Submit a poll intake task for the timeline

								OpIntakePoll intake_payload = new OpIntakePoll();
								intake_payload.setup (timeline_id);

								PendingTask.submit_task (
									EVID_POLL,									// event id
									sg.task_disp.get_time() + total_delay,		// sched_time
									sg.task_disp.get_time(),					// submit_time
									SUBID_AAFS,									// submit_id
									OPCODE_INTAKE_POLL,							// opcode
									0,											// stage
									intake_payload.marshal_task());				// details

								++count_withdrawn_timeline;

								// Adjust total and short delay time

								if (which_poll == 2 && potential_time < search_time_short) {
									f_seen_long = true;
								}
								if (f_seen_long) {
									total_delay += sg.task_disp.get_action_config().get_poll_long_intake_gap();
								} else {
									total_delay += sg.task_disp.get_action_config().get_poll_short_intake_gap();
									short_delay = total_delay;
								}
							}
						}
					}
				}
			}
		}

		System.out.println ("TASK-INFO: Comcat poll found " + count_no_timeline + " potential events with no corresponding timeline");
		System.out.println ("TASK-INFO: Comcat poll found " + count_withdrawn_timeline + " potential events with a withdrawn timeline");

		sg.log_sup.report_comcat_poll_done (poll_lookback, count_no_timeline, count_withdrawn_timeline);

		//--- Final steps

		// Update the poll timers

		sg.poll_sup.update_last_poll_time (which_poll, total_delay, short_delay);

		// Stage the task to repeat at the next poll time

		sg.task_disp.set_taskres_stage (sg.poll_sup.get_next_poll_time(),
							task.get_stage());
		return RESCODE_STAGE_REPEATING_TASK;

		//sg.task_disp.set_taskres_stage (sg.task_disp.get_time()
		//					+ sg.task_disp.get_action_config().get_poll_short_period(),
		//					task.get_stage());
		//return RESCODE_STAGE_REPEATING_TASK;
	}




	//----- Construction -----


	// Default constructor.

	public ExPollComcatRun () {}

}
