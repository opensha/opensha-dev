package scratch.aftershockStatistics.aafs;

import java.util.List;
import java.util.LinkedHashSet;

import java.io.IOException;
import java.io.PrintStream;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;
import scratch.aftershockStatistics.util.TimeSplitOutputStream;

import scratch.aftershockStatistics.ComcatException;
import scratch.aftershockStatistics.CompactEqkRupList;


/**
 * Task dispatcher for AAFS server.
 * Author: Michael Barall 03/18/2018.
 *
 * This class pulls tasks off the pending task queue and executes them.
 * It can run within its own thread.
 */
public class TaskDispatcher extends ServerComponent implements Runnable {

	//----- Dispatcher state variables -----

	// The time the task dispatcher was started.
	// Note: This is true time.

	private long start_time = 0L;

	// get_start_time - Get the time the task dispatcher was started.

	public long get_start_time () {
		return start_time;
	}




	// The last time that the task dispatcher was active.
	// Note: This is true time.

	private volatile long active_time = 0L;

	// get_active_time - Return the last time the dispatcher was active, or a special value.
	// Note: Can be called from any thread.

	public long get_active_time () {
		return active_time;
	}




	// The time the task dispatcher was last restarted, or 0L if never restarted.
	// Note: This is true time.

	private long restart_time = 0L;




	// Current dispatcher state, one of the special values below.

	private int dispatcher_state = 0;

	public static final int STATE_INITIAL         = 0;	// Initial state, dispatcher not started yet.
	public static final int STATE_FIRST_CONNECT   = 1;	// Connecting to MongoDB for the first time.
	public static final int STATE_RECONNECTING    = 2;	// Re-connecting to MongoDB.
	public static final int STATE_PROCESSING      = 3;	// Processing task.
	public static final int STATE_WAITING         = 4;	// Waiting for next task.
	public static final int STATE_POLLING         = 5;	// Polling queue.
	public static final int STATE_SHUTDOWN        = 6;	// Normal shutdown.

	// get_dispatcher_state - Get the task dispatcher state.
	// Note: This should only be called after the dispatcher exits.
	// The final state should be STATE_SHUTDOWN for normal exit.
	// Any other final state identifies when the failure occurred.

	public int get_dispatcher_state () {
		return dispatcher_state;
	}




	//----- Task context -----
	//
	// These variables are used by the task dispatcher to supply context for the currently-executing task.


	// Effective time at which the current task began to execute.

	private long dispatcher_time = 0L;

	// True time at which the current task began to execute.

	private long dispatcher_true_time = 0L;

	// Action configuration parameters for the current task.

	private ActionConfig dispatcher_action_config = null;




	// Get the effective time at which the current task began to execute.

	public long get_time () {
		return dispatcher_time;
	}




	// Get the true time at which the current task began to execute.

	public long get_true_time () {
		return dispatcher_true_time;
	}




	// Get the action configuration parameters for the current task.

	public ActionConfig get_action_config () {
		return dispatcher_action_config;
	}




	//----- Task results -----
	//
	// These variables are used by the currently-executing task to communicate results to the task dispatcher.


	// Time to insert in log entry for current task.
	// Defaults to dispatcher_time.

	private long taskres_log_time = 0L;

	// Remark to insert in log entry for current task.
	// Defaults to "".

	private String taskres_log_remark = "";

	// Execution time to use when staging current task.

	private long taskres_exec_time = 0L;

	// Stage to use when staging current task.

	private int taskres_stage = 0;

	// Event ID to use when staging current task, or null to leave event ID unchanged.

	private String taskres_event_id = null;




	// Set the remark to be used in a log entry, during task disposition.

	public void set_taskres_log (String log_remark) {
		taskres_log_remark = log_remark;
		return;
	}




	// Set and display the remark to be used in a log entry, during task disposition.

	public void set_display_taskres_log (String log_remark) {
		taskres_log_remark = log_remark;
		System.err.println (log_remark);
		return;
	}




	// Set the time and remark to be used in a log entry, during task disposition.

	public void set_taskres_log (long log_time, String log_remark) {
		taskres_log_time = log_time;
		taskres_log_remark = log_remark;
		return;
	}




	// Set the execution time and stage to be used when staging a task, during task disposition.

	public void set_taskres_stage (long exec_time, int stage) {
		taskres_exec_time = exec_time;
		taskres_stage = stage;
		taskres_event_id = null;
		return;
	}




	// Set the execution time and stage to be used when staging a task, during task disposition.

	public void set_taskres_stage (long exec_time, int stage, String event_id) {
		taskres_exec_time = exec_time;
		taskres_stage = stage;
		taskres_event_id = event_id;
		return;
	}




	// Display informational message, during task execution.

	public void display_taskinfo (String info) {
		System.out.println (info);
		return;
	}




	// Set dispatcher state to trigger shutdown, during task execution.

	public void set_state_shutdown () {
		dispatcher_state = STATE_SHUTDOWN;
		return;
	}




	//----- Dispatcher parameters -----

	// The polling delay, in milliseconds.

	private long polling_delay = 30000L;			// 30 seconds

	// The minimum polling delay, in milliseconds.

	private long polling_delay_min = 1000L;			// 1 second

	// The minimum delay before restarting after failure, in milliseconds.

	private long restart_delay_min = 20000L;		// 20 seconds

	// The maximum delay before restarting after failure, in milliseconds.

	private long restart_delay_max = 300000L;		// 5 minutes

	// True to enable messages at beginning and end of each task.

	private boolean dispatcher_verbose = true;




	//----- Idle time -----

	// The console log output stream, or null if none.

	private TimeSplitOutputStream console_log_tsop = null;

	// The summary log output stream, or null if none.

	private TimeSplitOutputStream summary_log_tsop = null;

	// List of all time split output streams.

	private LinkedHashSet<TimeSplitOutputStream> tsop_list;


	// Get the console log output stream, or null if none.

	public TimeSplitOutputStream get_console_log_tsop () {
		return console_log_tsop;
	}


	// Set the console log output stream, or null if none.

	public void set_console_log_tsop (TimeSplitOutputStream the_console_log_tsop) {
		remove_tsop (console_log_tsop);
		console_log_tsop = the_console_log_tsop;
		add_tsop (console_log_tsop);
		return;
	}


	// Get the summary log output stream, or null if none.

	public TimeSplitOutputStream get_summary_log_tsop () {
		return summary_log_tsop;
	}


	// Set the summary log output stream, or null if none.

	public void set_summary_log_tsop (TimeSplitOutputStream the_summary_log_tsop) {
		remove_tsop (summary_log_tsop);
		summary_log_tsop = the_summary_log_tsop;
		add_tsop (summary_log_tsop);

		// Set summary log destination

		sg.log_sup.set_summary_log_out (summary_log_tsop);

		return;
	}




	// Add a time split output stream.
	// If tsop is null, then perform no operation.
	// Note: The idle time code ignores any exceptions thrown by tsop.redirect().
	// If exception handling is required, you should subclass TimeSplitOutputStream
	// and override the redirect() method.

	public void add_tsop (TimeSplitOutputStream tsop) {
		if (tsop != null) {
			tsop_list.add (tsop);
		}
		return;
	}




	// Remove a time split output stream, if it is currently in the list.
	// If upstream is null, then perform no operation.

	public void remove_tsop (TimeSplitOutputStream tsop) {
		if (tsop != null) {
			tsop_list.remove (tsop);
		}
		return;
	}




	// Run idle time operations.
	// On entry, task context variables are set up:
	//  dispatcher_time, dispatcher_true_time, dispatcher_action_config

	private void exec_idle_time () {

		// Redirect time split output streams

		for (TimeSplitOutputStream tsop : tsop_list) {
			try {
				tsop.redirect (dispatcher_true_time);
			} catch (IOException e) {
			}
		}

		return;
	}




	//----- Construction -----

	public TaskDispatcher () {

		// Idle time initialization

		console_log_tsop = null;
		summary_log_tsop = null;
		tsop_list = new LinkedHashSet<TimeSplitOutputStream>();

		// Create and initialize the server group

		setup (new ServerGroup());
		sg.alloc_comp (this);
	}




	//----- Task posting functions -----




	/**
	 * post_task - Post a task.
	 * @param event_id = Event associated with this task, or "" if none. Cannot be null.
	 * @param sched_time = Time at which task should execute, in milliseconds
	 *                     since the epoch. Must be positive.
	 * @param submit_time = Time at which the task is submitted, in milliseconds
	 *                      since the epoch. Must be positive.
	 * @param submit_id = Person or entity submitting this task. Cannot be empty or null.
	 * @param opcode = Operation code used to dispatch the task.
	 * @param stage = Stage number, user-defined, effectively an extension of the opcode.
	 * @param details = Further details of this task. Can be null if there are none.
	 * @return
	 * Returns true if task is successfully posted, false if error.
	 * This function connects to MongoDB, puts the task on the queue, and then
	 * disconnects from MongoDB.
	 * If already connected to MongoDB, then use PendingTask.submit_task instead.
	 */
	public static boolean post_task (String event_id, long sched_time, long submit_time,
								String submit_id, int opcode, int stage, MarshalWriter details) {

		boolean result = true;

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Post the command

			PendingTask.submit_task (event_id, sched_time, submit_time,
									submit_id, opcode, stage, details);

			// Normal return

			result = true;

		// Abnormal return

        } catch (Exception e) {
			result = false;
            e.printStackTrace();
        } catch (Throwable e) {
			result = false;
            e.printStackTrace();
        }

		return result;
	}




	/**
	 * post_shutdown - Post a shutdown task to the queue.
	 * @param submit_id = Person or entity submitting this task. Cannot be empty or null.
	 * @return
	 * Returns true if shutdown is successfully posted, false if error.
	 */
	public static boolean post_shutdown (String submit_id) {

		boolean result = post_task (EVID_SHUTDOWN, EXEC_TIME_SHUTDOWN, ServerClock.get_true_time(),
								submit_id, OPCODE_SHUTDOWN, 0, null);

		return result;
	}




	//----- Task dispatching functions -----




	// Run the task dispatcher.

	@Override
	public void run() {

		// State = first connection

		dispatcher_state = STATE_FIRST_CONNECT;

		// Set up timers

		start_time = ServerClock.get_true_time();
		active_time = start_time;
		restart_time = start_time - restart_delay_max;

		// Restart loop, continue until shutdown or failure

		for (;;) {

			// Active task, null if none

			PendingTask task = null;

			// Connect to MongoDB

			try (
				MongoDBUtil mongo_instance = new MongoDBUtil();
			){

				// If first connection ...

				if (dispatcher_state == STATE_FIRST_CONNECT) {

					sg.log_sup.report_dispatcher_start ();
				
					// Remove any shutdown commands from the task queue

					List<PendingTask> shutdown_tasks = PendingTask.get_task_entry_range (0L, 0L, EVID_SHUTDOWN);

					for (PendingTask shutdown_task : shutdown_tasks) {
						if (shutdown_task.get_opcode() == OPCODE_SHUTDOWN) {
							PendingTask.delete_task (shutdown_task);
						}
					}

					// Initialize Comcat polling, begin disabled

					sg.poll_sup.init_polling_disabled();
				}

				else {
					sg.log_sup.report_dispatcher_restart ();
				}

				// Polling loop, continue until shutdown or exception

				while (dispatcher_state != STATE_SHUTDOWN) {

					// State = polling

					dispatcher_state = STATE_POLLING;

					// Get task time and configuration

					dispatcher_time = ServerClock.get_time();
					dispatcher_true_time = ServerClock.get_true_time();
					dispatcher_action_config = new ActionConfig();

					// Record the dispatcher active time

					active_time = dispatcher_true_time;

					// Get the next task on the pending queue, that's ready to execute, and activate it

					long cutoff_time = dispatcher_time;
					task = PendingTask.activate_first_ready_task (cutoff_time);

					// Check if there is a task

					if (task == null) {

						// State = waiting

						dispatcher_state = STATE_WAITING;

						// Execute idle time operations

						exec_idle_time();

						// Get polling delay, allowing for time consumed by idle time operations

						long eff_polling_delay = dispatcher_true_time + polling_delay - ServerClock.get_true_time();
						if (eff_polling_delay > polling_delay) {
							eff_polling_delay = polling_delay;
						}

						// Wait for the polling delay

						if (eff_polling_delay >= polling_delay_min) {
							try {
								Thread.sleep(eff_polling_delay);
							} catch (InterruptedException e) {
							}
						}

					} else {

						// State = processing

						dispatcher_state = STATE_PROCESSING;

						// Dispatch on opcode

						dispatch_task (task);

						// No active task

						task = null;
					}
				}

			// Operation failed with exception

			} catch (Exception e) {
				e.printStackTrace();
				if (task != null) {
					System.err.println ("Failing task: " + task.toString());
				}
				sg.log_sup.report_dispatcher_exception (task, e);
			} catch (Throwable e) {
				e.printStackTrace();
				if (task != null) {
					System.err.println ("Failing task: " + task.toString());
				}
				sg.log_sup.report_dispatcher_exception (task, e);
			}

			// If normal shutdown, exit the restart loop

			if (dispatcher_state == STATE_SHUTDOWN) {
				sg.log_sup.report_dispatcher_shutdown ();
				break;
			}

			// If first connection failed, exit the restart loop
			// (because this might be a configuration error)

			if (dispatcher_state == STATE_FIRST_CONNECT) {
				sg.log_sup.report_dispatcher_immediate_shutdown ();
				break;
			}

			// Calculate restart delay:
			// restart time is restart_delay_min after the current time,
			// or restart_delay_max after the last restart, whichever is later

			long restart_delay = Math.max (restart_delay_min,
						restart_delay_max + restart_time - ServerClock.get_true_time());

			// Wait until time for restart

			try {
				Thread.sleep(restart_delay);
			} catch (InterruptedException e) {
			}

			// New restart time

			restart_time = ServerClock.get_true_time();

			// State = reconnecting
		
			dispatcher_state = STATE_RECONNECTING;

		}

		return;
	}




	/**
	 * run_next_task - Run the next task.
	 * @param f_verbose = True to display the task to System.out.
	 * @param f_adjust_time = True to adjust time to the apparent exeuction time of the task.
	 * @return
	 * Returns true if task is executed, false if task queue is empty.
	 * This function connects to MongoDB, executes the next task on the queue, and then
	 * disconnects from MongoDB.
	 * This function also sets the clock to the task execution time.
	 */
	public boolean run_next_task (boolean f_verbose, boolean f_adjust_time) {

		boolean result = true;

		// Active task, null if none

		PendingTask task = null;

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

			// Get task time and configuration

			dispatcher_true_time = ServerClock.get_true_time();
			dispatcher_action_config = new ActionConfig();

			// Get the next task on the pending queue, and activate it

			long cutoff_time = EXEC_TIME_FAR_FUTURE;
			task = PendingTask.activate_first_ready_task (cutoff_time);

			// Check if there is a task

			if (task == null) {

				// No task

				result = false;

				// If verbose, write message

				if (f_verbose) {
					System.out.println ("No task available");
				}

			} else {

				// Got a task

				result = true;

				// Adjust the time

				if (f_adjust_time) {
					ServerClock.advance_frozen_time (task.get_apparent_time());
				}
				dispatcher_time = ServerClock.get_time();

				// If verbose, write message

				if (f_verbose) {
					System.out.println ("Executing task: " + task.toString());
				}

				// Dispatch on opcode

				dispatch_task (task);

				// No active task

				task = null;
			}

		// Abnormal return

        } catch (Exception e) {
            e.printStackTrace();
			if (task != null) {
				System.err.println ("Failing task: " + task.toString());
			}
        } catch (Throwable e) {
            e.printStackTrace();
			if (task != null) {
				System.err.println ("Failing task: " + task.toString());
			}
        }

		return result;
	}




	// Dispatch task based on its opcode.
	//
	// Each task must have an execution function that takes the task as its argument.
	// On entry to the execution function, the following are set:
	//   dispatcher_time = Effective time at which task is beginning to execute.
	//   dispatcher_true_time = True clock time at which task is beginning to execute.
	//   dispatcher_action_config = Configuration parameters.
	//   taskres_log_time = dispatcher_time.
	//   taskres_log_remark = "".
	// The execution function must return a result code (RESCODE_XXXXX).
	//   RESCODE_DELETE_XXXXX deletes the current task, with no other action.
	//   RESCODE_STAGE_XXXXX stages the current task, using taskres_exec_time,
	//    taskres_stage, and taskres_event_id.
	//   Any other return writes a log entry with the given result code, using
	//     taskres_log_time and taskres_log_remark; and then deletes the current task.
	// If the task is being restarted and a log entry has already been written,
	//   then the task is deleted and the execution function is not called.

	private void dispatch_task (PendingTask task) {

		// If restarting ...

		if (task.is_restarted()) {
		
			// If we wrote a log entry ...

			if (LogEntry.get_log_entry_for_key (task.get_record_key()) != null) {
			
				// Just remove the task from the queue

				PendingTask.delete_task (task);
				return;
			}
		}

		// Establish task context
		// (Note that dispatcher_time, dispatcher_true_time, and dispatcher_action_config
		// are established by our caller)

		taskres_log_time = dispatcher_time;
		taskres_log_remark = "";

		taskres_exec_time = 0L;
		taskres_stage = 0;
		taskres_event_id = null;

		// Say hello

		if (dispatcher_verbose) {

			if (task.is_restarted()) {

				System.out.println (LOG_SEPARATOR_LINE);
				display_taskinfo ("TASK-RESTART: " + SimpleUtils.time_to_string (dispatcher_time) + "\n"
					+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
					+ "event_id = " + task.get_event_id() + "\n"
					+ "stage = " + task.get_stage());

			} else {

				System.out.println (LOG_SEPARATOR_LINE);
				display_taskinfo ("TASK-BEGIN: " + SimpleUtils.time_to_string (dispatcher_time) + "\n"
					+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
					+ "event_id = " + task.get_event_id() + "\n"
					+ "stage = " + task.get_stage());

			}

		}

		if (task.is_restarted()) {
			sg.log_sup.report_task_restart (task);
		} else {
			sg.log_sup.report_task_begin (task);
		}

		// Invoke the execution function, depending on the opcode

		int opix = task.get_opcode();

		if (opix < OPCODE_MIN || opix > OPCODE_MAX) {
			opix = OPCODE_UNKNOWN;
		}

		int rescode = sg.dispatch_table[opix].exec_task (task);

		// Handle task disposition, depending on the result code

		switch (rescode) {

		default:

			// Display message

			if (dispatcher_verbose) {
				display_taskinfo ("TASK-END:\n"
					+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
					+ "rescode = " + get_rescode_as_string (rescode));
			}
			sg.log_sup.report_task_end (task, rescode);

			// Log the task

			LogEntry.submit_log_entry (task, taskres_log_time, rescode, taskres_log_remark);

			// Remove the task from the queue

			PendingTask.delete_task (task);

			break;

		case RESCODE_DELETE:
		case RESCODE_DELETE_TIMELINE_EXISTS:
		case RESCODE_DELETE_INTAKE_FILTERED:
		case RESCODE_DELETE_INTAKE_AGED:
		case RESCODE_DELETE_TIMELINE_NO_ALIAS:
		case RESCODE_DELETE_TIMELINE_BAD_STATE:
		case RESCODE_DELETE_NOT_IN_COMCAT:

			// Display message

			if (dispatcher_verbose) {
				display_taskinfo ("TASK-DELETE:\n"
					+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
					+ "rescode = " + get_rescode_as_string (rescode));
			}
			sg.log_sup.report_task_delete (task, rescode);

			// Remove the task from the queue

			PendingTask.delete_task (task);

			break;

		case RESCODE_STAGE:
		case RESCODE_STAGE_COMCAT_RETRY:
		case RESCODE_STAGE_PDL_RETRY:
		case RESCODE_STAGE_EVENT_ID:
		case RESCODE_STAGE_TIMELINE_ID:
		case RESCODE_STAGE_TOO_SOON:
		case RESCODE_STAGE_REPEATING_TASK:

			// Display message

			if (dispatcher_verbose) {
				display_taskinfo ("TASK-STAGE:\n"
					+ "opcode = " + get_opcode_as_string (task.get_opcode()) + "\n"
					+ "rescode = " + get_rescode_as_string (rescode) + "\n"
					+ "taskres_exec_time = " + SimpleUtils.time_to_string (taskres_exec_time) + "\n"
					+ "taskres_stage = " + taskres_stage + "\n"
					+ "taskres_event_id = " + ((taskres_event_id == null) ? "null" : taskres_event_id) );
			}
			sg.log_sup.report_task_stage (task, rescode, taskres_event_id, taskres_stage, taskres_exec_time);

			// Stage the task, so it will execute again

			PendingTask.stage_task (task, taskres_exec_time, taskres_stage, taskres_event_id);

			break;
		}

		return;
	}




	//----- Test functions -----




	// Set up the task context.
	// This sets up the task context as if a task were being dispatched.
	// It can be used when testing routines that expect the task context to be available.

	public void setup_task_context () {

		dispatcher_time = ServerClock.get_time();
		dispatcher_true_time = ServerClock.get_true_time();
		dispatcher_action_config = new ActionConfig();

		taskres_log_time = dispatcher_time;
		taskres_log_remark = "";

		taskres_exec_time = 0L;
		taskres_stage = 0;
		taskres_event_id = null;
	
		return;
	}




	// Get the server group.
	// This can be used for testing, to make calls directly into support routines.

	public ServerGroup get_server_group () {
		return sg;
	}




	// Execute the idle time operations.
	// This set up the task context and executes the idle time operations.
	// This is a test function.

	public void test_exec_idle_time () {

		dispatcher_time = ServerClock.get_time();
		dispatcher_true_time = ServerClock.get_true_time();
		dispatcher_action_config = new ActionConfig();
	
		exec_idle_time();

		return;
	}

}
