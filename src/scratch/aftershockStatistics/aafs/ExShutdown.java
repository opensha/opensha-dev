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
 * Execute task: Shutdown.
 * Author: Michael Barall 06/25/2018.
 */
public class ExShutdown extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_shutdown (task);
	}




	// Execute shutdown.

	private int exec_shutdown (PendingTask task) {

		// Signal polling loop to exit

		sg.task_disp.set_state_shutdown();

		// Remove the task from the queue

		return RESCODE_DELETE;
	}




	//----- Construction -----


	// Default constructor.

	public ExShutdown () {}

}
