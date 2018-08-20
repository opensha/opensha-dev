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
 * Support functions for tasks.
 * Author: Michael Barall 06/23/2018.
 */
public class TaskSupport extends ServerComponent {




	//----- Task execution subroutines : General task management -----




	// Display an exception caused by an invalid task, and set the log remark.

	public void display_invalid_task (PendingTask task, Exception e) {

		// Display messages

		String trace = SimpleUtils.getStackTraceAsString(e);
		sg.task_disp.display_taskinfo (trace + "\nInvalid task: " + task.toString() + "\nDump:\n" + task.dump_details());

		sg.log_sup.report_invalid_task_exception (task, e);

		// Insert the stack trace into the log remark

		sg.task_disp.set_taskres_log (trace);
		return;
	}




	// Delete all waiting tasks with the given event id and any of the given opcodes.
	// Note: The currently active task is not deleted, even if it would match.

	public void delete_all_waiting_tasks (String event_id, int... opcodes) {

		// Get a list of all waiting tasks for the given event

		List<PendingTask> tasks = PendingTask.get_task_entry_range (EXEC_TIME_MIN_WAITING, 0L, event_id);

		// Delete tasks with the given opcode

		for (PendingTask task : tasks) {
			for (int opcode : opcodes) {
				if (task.get_opcode() == opcode) {
					PendingTask.delete_task (task);
					break;
				}
			}
		}
	
		return;
	}




	// Delete all tasks with the given event id and any of the given opcodes.
	// Note: The currently active task is deleted, if it matches.

	public void delete_all_tasks_for_event (String event_id, int... opcodes) {

		// Get a list of all waiting tasks for the given event

		List<PendingTask> tasks = PendingTask.get_task_entry_range (0L, 0L, event_id);

		// Delete tasks with the given opcode

		for (PendingTask task : tasks) {
			for (int opcode : opcodes) {
				if (task.get_opcode() == opcode) {
					PendingTask.delete_task (task);
					break;
				}
			}
		}
	
		return;
	}




	// Get the next prompt execution time to use.
	// Note: This is intended for tasks that execute reasonably quickly, without
	// contacting external services such as ComCat.  These tasks have priority over
	// shutdown, and execute last-in-first-out (although execution order should not
	// be considered guaranteed).

	public long get_prompt_exec_time () {

		// Get the next-up prompt task

		PendingTask prompt_task = PendingTask.get_first_task_entry (EXEC_TIME_MIN_PROMPT, EXEC_TIME_MAX_PROMPT, null);

		// If none, use the largest value

		if (prompt_task == null) {
			return EXEC_TIME_MAX_PROMPT;
		}

		// Otherwise, return one less that the next-up task

		return Math.max (EXEC_TIME_MIN_PROMPT, prompt_task.get_exec_time() - 1L);
	}




	//----- Construction -----


	// Default constructor.

	public TaskSupport () {}

}
