package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;

import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;

/**
 * Task dispatcher for AAFS server.
 * Author: Michael Barall 03/18/2018.
 *
 * This class pulls tasks off the pending task queue and executes them.
 * It can run within its own thread.
 */
public class TaskDispatcher implements Runnable {

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




	//----- Dispatcher parameters -----

	// The polling delay, in milliseconds.

	private long polling_delay = 30000L;			// 30 seconds

	// The minimum delay before restarting after failure, in milliseconds.

	private long restart_delay_min = 20000L;		// 20 seconds

	// The maximum delay before restarting after failure, in milliseconds.

	private long restart_delay_max = 300000L;		// 5 minutes




	// Opcodes.

	public static final int OPCODE_MIN = 0;					// Minimum allowed opcode
	public static final int OPCODE_NO_OP = 0;				// No operation
	public static final int OPCODE_SHUTDOWN = 1;			// Shut down the server
	public static final int OPCODE_CON_MESSAGE = 2;			// Write message to console
	public static final int OPCODE_DELAYED_SHUTDOWN = 3;	// Delayed shut down of server
	public static final int OPCODE_MAX = 3;					// Maximum allowed opcode




	// Special execution times.

	public static final long EXEC_TIME_ACTIVE = 0L;			// Task is active
	public static final long EXEC_TIME_SHUTDOWN = 1L;		// Execution time for shutdown task
	public static final long EXEC_TIME_FAR_FUTURE = 1000000000000000L;	// 10^15 ~ 30,000 years




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

		boolean result = post_task ("", EXEC_TIME_SHUTDOWN, ServerClock.get_true_time(),
								submit_id, OPCODE_SHUTDOWN, 0, null);

		return result;
	}




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

				// Polling loop, continue until shutdown or exception

				while (dispatcher_state != STATE_SHUTDOWN) {

					// State = polling

					dispatcher_state = STATE_POLLING;

					// Record the dispatcher active time

					active_time = ServerClock.get_true_time();

					// Get the next task on the pending queue, that's ready to execute, and activate it

					long cutoff_time = ServerClock.get_time();
					task = PendingTask.activate_first_ready_task (cutoff_time);

					// Check if there is a task

					if (task == null) {

						// State = waiting

						dispatcher_state = STATE_WAITING;

						// If none, wait for the polling delay

						try {
							Thread.sleep(polling_delay);
						} catch (InterruptedException e) {
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
			} catch (Throwable e) {
				e.printStackTrace();
				if (task != null) {
					System.err.println ("Failing task: " + task.toString());
				}
			}

			// If normal shutdown, exit the restart loop

			if (dispatcher_state == STATE_SHUTDOWN) {
				break;
			}

			// If first connection failed, exit the restart loop
			// (because this might be a configuration error)

			if (dispatcher_state == STATE_FIRST_CONNECT) {
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
	 * @return
	 * Returns true if task is executed, false if task queue is empty.
	 * This function connects to MongoDB, executes the next task on the queue, and then
	 * disconnects from MongoDB.
	 * This function also sets the clock to the task execution time.
	 */
	public boolean run_next_task (boolean f_verbose) {

		boolean result = true;

		// Active task, null if none

		PendingTask task = null;

		// Connect to MongoDB

		try (
			MongoDBUtil mongo_instance = new MongoDBUtil();
		){

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

				ServerClock.advance_frozen_time (task.get_apparent_time());

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

	private void dispatch_task (PendingTask task) {

		switch (task.get_opcode()) {

		default: exec_unknown (task); break;
				
		case OPCODE_NO_OP: exec_no_op (task); break;

		case OPCODE_SHUTDOWN: exec_shutdown (task); break;

		case OPCODE_CON_MESSAGE: exec_con_message (task); break;

		case OPCODE_DELAYED_SHUTDOWN: exec_delayed_shutdown (task); break;

		}

		return;
	}




	// Execute unknown opcode.

	private void exec_unknown (PendingTask task) {

		// Remove the task from the queue

		PendingTask.delete_task (task);

		// Throw exception

		throw new RuntimeException("TaskDispatcher: Invalid opcode\n" + task.toString());

		//return;	// would be unreachable
	}




	// Execute no operation.

	private void exec_no_op (PendingTask task) {

		// Remove the task from the queue

		PendingTask.delete_task (task);
		return;
	}




	// Execute shutdown.

	private void exec_shutdown (PendingTask task) {

		// Remove the task from the queue

		PendingTask.delete_task (task);

		// If the task was submitted after our start time, then signal polling loop to exit
		// (The purpose of this test is to discard any shutdown tasks that were already
		// on the queue prior to the dispatcher being started.)

		if (task.get_submit_time() > start_time) {
			dispatcher_state = STATE_SHUTDOWN;
		}

		return;
	}




	// Execute console message.
	// For testing, this message supports the following stages:
	//  0 = Write message normally.
	//  1 = If not restarting, throw exception before writing message.
	//  2 = If not restarting, throw exception after writing message.

	private void exec_con_message (PendingTask task) {

		// If restarting ...

		if (task.is_restarted()) {
		
			// If we wrote a log entry ...

			if (LogEntry.get_log_entry_for_key (task.get_record_key()) != null) {
			
				// Just remove the task from the queue

				PendingTask.delete_task (task);
				return;
			}
		}

		// If stage 1 and not restarting, throw exception

		if (task.get_stage() == 1 && !task.is_restarted()) {
			throw new RuntimeException("TaskDispatcher.exec_con_message: Pre-message exception");
		}

		// Write message

		System.out.println (task.get_details());

		// Log the task

		LogEntry.submit_log_entry (task, ServerClock.get_time(), 0, "");

		// If stage 2 and not restarting, throw exception

		if (task.get_stage() == 2 && !task.is_restarted()) {
			throw new RuntimeException("TaskDispatcher.exec_con_message: Post-message exception");
		}

		// Remove the task from the queue

		PendingTask.delete_task (task);

		return;
	}




	// Execute delayed shutdown.

	private void exec_delayed_shutdown (PendingTask task) {

		// If restarting ...

		if (task.is_restarted()) {
			
			// Just remove the task from the queue

			PendingTask.delete_task (task);
			return;
		}

		// Remove the task from the queue

		PendingTask.delete_task (task);

		// Signal polling loop to exit

		dispatcher_state = STATE_SHUTDOWN;

		return;
	}







}
