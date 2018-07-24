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
 * Execute task: Unknown opcode.
 * Author: Michael Barall 06/25/2018.
 */
public class ExUnknown extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_unknown (task);
	}




	// Execute unknown opcode.

	private int exec_unknown (PendingTask task) {

		// Remove the task from the queue

		PendingTask.delete_task (task);

		// Throw exception

		throw new RuntimeException("Task has invalid opcode\n" + task.toString());

		//return RESCODE_DELETE;	// would be unreachable
	}




	//----- Construction -----


	// Default constructor.

	public ExUnknown () {}

}
