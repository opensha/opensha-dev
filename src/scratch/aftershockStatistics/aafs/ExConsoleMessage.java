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
 * Execute task: Console message.
 * Author: Michael Barall 06/25/2018.
 */
public class ExConsoleMessage extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_con_message (task);
	}




	// Execute console message.
	// For testing, this message supports the following stages:
	//  0 = Write message normally.
	//  1 = If not restarting, throw exception before writing message.
	//  2 = If not restarting, throw exception after writing message.

	private int exec_con_message (PendingTask task) {

		// Get payload and check for valid task

		OpConsoleMessage payload = new OpConsoleMessage();

		try {
			payload.unmarshal_task (task);

			if (task.get_stage() < 0 || task.get_stage() > 2) {
				throw new DBCorruptException("Invalid stage for console message task, stage = " + task.get_stage());
			}
		}

		// Invalid task

		catch (Exception e) {

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		// If stage 1 and not restarting, throw exception

		if (task.get_stage() == 1 && !task.is_restarted()) {
			throw new RuntimeException("ExConsoleMessage.exec_con_message: Pre-message exception");
		}

		// Write message

		System.out.println (payload.message);

		// Log the task
		// (We write this explicitly so we can throw an exception afterward)

		LogEntry.submit_log_entry (task, sg.task_disp.get_time(), RESCODE_SUCCESS, "");

		// If stage 2 and not restarting, throw exception

		if (task.get_stage() == 2 && !task.is_restarted()) {
			throw new RuntimeException("ExConsoleMessage.exec_con_message: Post-message exception");
		}

		// Remove the task from the queue

		return RESCODE_DELETE;
	}




	//----- Construction -----


	// Default constructor.

	public ExConsoleMessage () {}

}
