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

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.ComcatException;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

/**
 * Execute task: Stop Comcat polling.
 * Author: Michael Barall 08/05/2018.
 */
public class ExPollComcatStop extends ServerExecTask {


	//----- Task execution -----


	// Execute the task, called from the task dispatcher.
	// The parameter is the task to execute.
	// The return value is a result code.
	// Support functions, task context, and result functions are available through the server group.

	@Override
	public int exec_task (PendingTask task) {
		return exec_poll_comcat_stop (task);
	}




	// Stop Comcat polling.

	private int exec_poll_comcat_stop (PendingTask task) {

		//--- Get payload

		OpPollComcatStop payload = new OpPollComcatStop();

		try {
			payload.unmarshal_task (task);
		}

		// Invalid task

		catch (Exception e) {

			// Display the error and log the task

			sg.task_sup.display_invalid_task (task, e);
			return RESCODE_TASK_CORRUPT;
		}

		//--- Poll status

		// Disable polling

		sg.poll_sup.set_polling_disabled();

		//--- Final steps

		// Log the task

		return RESCODE_SUCCESS;
	}




	//----- Construction -----


	// Default constructor.

	public ExPollComcatStop () {}

}
