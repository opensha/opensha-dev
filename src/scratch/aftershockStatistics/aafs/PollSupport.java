package scratch.aftershockStatistics.aafs;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.EventNotFoundException;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.ComcatException;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

/**
 * Support functions for polling Comcat.
 * Author: Michael Barall 08/05/2018.
 */
public class PollSupport extends ServerComponent {




	//----- Timing variables -----

	// True if polling is currently enabled, false if not.

	private boolean f_polling_enabled;

	// True to allow the next poll to be long.

	private boolean f_long_poll_ok;

	// Time at which the next short poll can occur, in milliseconds since the epoch.

	private long next_short_poll_time;

	// Time at which the next long poll can occur, in milliseconds since the epoch.

	private long next_long_poll_time;




	//----- Polling subroutines -----




//	// Delete all waiting polling tasks.
//	// The currently active task is not deleted, even if it is a polling task.
//
//	public void delete_all_waiting_polling_tasks () {
//		sg.task_sup.delete_all_waiting_tasks (EVID_POLL, OPCODE_POLL_COMCAT_RUN, OPCODE_POLL_COMCAT_START, OPCODE_POLL_COMCAT_STOP);
//		return;
//	}




	// Delete all polling tasks (the poll Comcat run/start/stop tasks, and delayed poll intake tasks).
	// The currently active task is deleted, if it is a polling task.

	public void delete_all_existing_polling_tasks () {
		sg.task_sup.delete_all_tasks_for_event (EVID_POLL, OPCODE_INTAKE_POLL, OPCODE_POLL_COMCAT_RUN, OPCODE_POLL_COMCAT_START, OPCODE_POLL_COMCAT_STOP);
		return;
	}




	// Delete the poll Comcat run task, and all delayed poll intake tasks, if any.
	// The currently active task is not deleted, even if it would match.

	public void delete_waiting_poll_run_tasks () {
		sg.task_sup.delete_all_waiting_tasks (EVID_POLL, OPCODE_INTAKE_POLL, OPCODE_POLL_COMCAT_RUN);
		return;
	}




	// Set polling to the disabled state.
	// Also deletes the poll Comcat run task, if any.
	// Note: This is to be called from the execution function of a polling task.

	public void set_polling_disabled () {
		delete_waiting_poll_run_tasks ();

		f_polling_enabled = false;

		return;
	}




	// Set polling to the enabled state.
	// Also deletes the poll Comcat run task, if any, and issues a new one.
	// Note: This is to be called from the execution function of a polling task.
	// Note: If polling is already enabled, the cycle is reset.

	public void set_polling_enabled () {
		delete_waiting_poll_run_tasks ();

		f_polling_enabled = true;
		f_long_poll_ok = false;
		next_short_poll_time = sg.task_disp.get_time();
		next_long_poll_time = sg.task_disp.get_time();

		// Kick off polling

		kick_off_polling (0L);
		return;
	}




	// Initialize polling to the disabled state.
	// Also deletes all polling tasks.
	// Note: This is to be called from the task dispatcher, not during execution of a task.
	// This must be called before the first task is executed from the task queue.

	public void init_polling_disabled () {
		delete_all_existing_polling_tasks ();

		f_polling_enabled = false;
		f_long_poll_ok = false;
		next_short_poll_time = 0L;
		next_long_poll_time = 0L;

		return;
	}




	// Issue the poll Comcat run task.
	// The first poll happens no earlier than the given delay after the current time, in milliseconds.

	public void kick_off_polling (long delay) {

		OpPollComcatRun poll_run_payload = new OpPollComcatRun();
		poll_run_payload.setup ();

		PendingTask.submit_task (
			EVID_POLL,										// event id
			Math.max (sg.task_disp.get_time() + delay, get_next_poll_time()),	// sched_time
			sg.task_disp.get_time(),						// submit_time
			SUBID_AAFS,										// submit_id
			OPCODE_POLL_COMCAT_RUN,							// opcode
			0,												// stage
			poll_run_payload.marshal_task());				// details

		return;
	}




	// Check if it is time to perform a poll.
	// Return value indicates which poll:
	//  0 = No poll now.
	//  1 = Do short poll.
	//  2 = Do long poll.

	public int check_if_poll_time () {

		// If polling is enabled ...

		if (f_polling_enabled) {

			// Jitter allowance is 40% of the short poll time

			long jitter = sg.task_disp.get_action_config().get_poll_short_period() * 2L / 5L;

			// Effective time is current time plus jitter allowance

			long eff_time = sg.task_disp.get_time() + jitter;

			// If it's time for a short poll ...

			if (eff_time >= next_short_poll_time) {

				// If long poll is allowed ...

				if (f_long_poll_ok) {

					// If it's time for a long poll, return it

					if (eff_time >= next_long_poll_time) {
						return 2;
					}
				}

				// Otherwise, return short poll

				return 1;
			}
		}

		// No poll

		return 0;
	}




	// Update the last poll time.
	// Parameters:
	//  which_poll = 1 for short poll, 2 for long poll (must be return value from check_if_poll_time).
	//  total_delay = Delay time consumed by all delayed intake commands, in milliseconds.
	//  short_delay = Delay time consumed by delayed intake commands within range of a short poll, in milliseconds.
	// Note: For a short poll, short_delay should equal total_delay.

	public void update_last_poll_time (int which_poll, long total_delay, long short_delay) {

		// Jitter allowance is 60% of the short poll time
		// (Jitter allowances here and in check_if_poll_time should sum to 100%)

		long jitter = sg.task_disp.get_action_config().get_poll_short_period() * 3L / 5L;

		// End time of delay, requiring at least the jitter allowance after current time

		long delay_end_time = sg.task_disp.get_time() + Math.max (short_delay, jitter);

		// End time of period

		long period_end_time = next_short_poll_time
			+ (sg.task_disp.get_action_config().get_poll_short_period());

		// New time is the later of the two

		next_short_poll_time = Math.max (delay_end_time, period_end_time);

		// Allow long poll

		f_long_poll_ok = true;

		// If we just did a long poll ...

		if (which_poll == 2) {

			// Jitter allowance is 60% of the long poll time

			jitter = sg.task_disp.get_action_config().get_poll_long_period() * 3L / 5L;

			// End time of delay, requiring at least the jitter allowance after current time

			delay_end_time = sg.task_disp.get_time() + Math.max (total_delay, jitter);

			// End time of period

			period_end_time = next_long_poll_time
				+ (sg.task_disp.get_action_config().get_poll_long_period());

			// New time is the later of the two

			next_long_poll_time = Math.max (delay_end_time, period_end_time);

			// Cannot do two long polls in a row

			f_long_poll_ok = false;
		}

		return;
	}




	// Get the recommended execution time of the next poll Comcat run command.
	// The returned value is always at least the current time.

	public long get_next_poll_time () {

		// It's the time of the next short poll, but at least the current time

		return Math.max (next_short_poll_time, sg.task_disp.get_time());
	}




	//----- Construction -----


	// Default constructor.

	public PollSupport () {

		f_polling_enabled = false;
		f_long_poll_ok = false;
		next_short_poll_time = 0L;
		next_long_poll_time = 0L;

	}


	// Set up this component by linking to the server group.
	// A subclass may override this to perform additional setup operations.

	@Override
	public void setup (ServerGroup the_sg) {
		super.setup (the_sg);

		f_polling_enabled = false;
		f_long_poll_ok = false;
		next_short_poll_time = 0L;
		next_long_poll_time = 0L;

		return;
	}

}
