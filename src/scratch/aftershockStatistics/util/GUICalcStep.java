package scratch.aftershockStatistics.util;

/**
 * One calculation step in a series.
 * Author: Michael Barall 08/18/2018.
 *
 * Derived from code originally in AftershockStatsGUI.
 *
 * CAUTION: Any calculation step that does not run in the event dispatch thread
 * may not make calls to any Swing operations, except for SwingUtilities.invokeLater.
 * Be aware that calls to Swing can be non-obvious, for example, as a result of change
 * watchers or console redirection.
 */
public class GUICalcStep {
		
	private String title;
	private String progressMessage;
	private Runnable run;
	private boolean runInEDT;

	// Specify the progress monitor title, message, and operation to run.
		
	public GUICalcStep(String title, String progressMessage, Runnable run) {
		this(title, progressMessage, run, false);
	}

	// Specify the progress monitor title, message, operation to run,
	// and flag specifying if operation must run in the event dispatch thread.
		
	public GUICalcStep(String title, String progressMessage, Runnable run, boolean runInEDT) {
		this.title = title;
		this.progressMessage = progressMessage;
		this.run = run;
		this.runInEDT = runInEDT;
	}

	// Get the title.

	public String get_title () {
		return title;
	}

	// Get the progress message.

	public String get_progressMessage () {
		return progressMessage;
	}

	// Get the operation to run.

	public Runnable get_run () {
		return run;
	}

	// Get the flag to run in the event dispatch thread.

	public boolean get_runInEDT () {
		return runInEDT;
	}

}
