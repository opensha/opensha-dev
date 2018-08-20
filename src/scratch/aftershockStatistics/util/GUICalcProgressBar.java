package scratch.aftershockStatistics.util;

import java.util.ArrayDeque;

import java.awt.Component;
import java.awt.Dialog.ModalityType;

import javax.swing.SwingUtilities;

import org.opensha.sha.gui.infoTools.CalcProgressBar;


/**
 * Manages a simple progress bar that can be controlled from a worker thread.
 * Author: Michael Barall 08/18/2018.
 *
 * There are three operations:  req_init to initialize and create the progress bar,
 * req_update to update the title and progress message, req_dispose to destroy
 * the progress bar, other functions to pass calls to the underlying CalcProgressBar.
 *
 * Alternatively, you can use try-with-resources.  Construct the progress bar with
 * init_now = true to create the progress bar.  Or to use an existing progress bar,
 * invoke req_init and save the return value.  Closing the progress bar destroys
 * the progress bar.  Then, use req_update to update the title and progress message.
 *
 * Calls to req_init must be made from a worker thread.  The update and dispose
 * functions can be called from any thread.
 *
 * This is a wrapper around org.opensha.sha.gui.infoTools.CalcProgressBar,
 * although not with exactly the same set of functions.
 */
public class GUICalcProgressBar implements AutoCloseable {

	// The owner of the progress bar.
	// This is set in the constructor.

	private Component owner;

	// The initial title and info.
	// This is set in the constructor.

	String initial_title;
	String initial_info;

	// The progress bar.
	// This field is accessed only from the event dispatch thread.

	private CalcProgressBar progress;

	// True if we are active (init written to queue, dispose not written to queue).
	// Access to this field must be synchronized.

	private boolean is_active;

	// Class to perform an action on the progress bar.

	private abstract class BarAction {

		// True if waiting for the action to be performed.

		private boolean waiting = true;

		// Wait for the action to be performed.
		// Note: If this is used, then run_action must call notify_performed or else it will deadlock.

		public void wait_performed () {
			synchronized (this) {
				while (waiting) {
					try {
						wait();
					} catch (InterruptedException e) {}
				}
			}
		}

		// Notify that the action has been performed.

		protected void notify_performed () {
			synchronized (this) {
				waiting = false;
				notifyAll();
			}
			return;
		}

		// Perform the action, return true to exit the pump.

		public abstract boolean run_action ();
	}

	// Queue of actions.
	// Access to this field must be synchronized.

	private ArrayDeque<BarAction> deque;

	// Read from the queue and perform the actions.
	// This function must be called from the event dispatch thread.

	private void read_queue () {
		for (;;) {
			BarAction action;
			synchronized (this) {
				action = deque.pollFirst();
			}
			if (action == null) {
				break;
			}
			boolean exit_pump = action.run_action();
			if (exit_pump) {
				break;
			}
		}
		return;
	}

	// Write an action to the queue.

	private void write_queue (BarAction action) {
		write_queue (action, false, false);
	}

	private void write_queue (BarAction action, boolean is_init, boolean is_dispose) {

		synchronized (this) {

			// Must be inactive if this is init, active otherwise

			if (is_init) {
				if (is_active) {
					throw new IllegalStateException("GUICalcProgressBar.write_queue: Progress bar is already initialized");
				}
				is_active = true;
			}
			else {
				if (!( is_active )) {
					throw new IllegalStateException("GUICalcProgressBar.write_queue: Progress bar is already initialized");
				}
			}

			// If it's dispose, mark inactive

			if (is_dispose) {
				is_active = false;
			}

			// Place on queue

			deque.addLast(action);
		}

		// If called from the event dispatch queue, execute immediately

		if ( SwingUtilities.isEventDispatchThread() ) {
			read_queue();
			return;
		}

		// Otherwise, execute later

		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				read_queue();
			}
		});
		return;
	}

	// Construct object.
	// The visibility flag is ignored.

	public GUICalcProgressBar (Component owner, String title, String info, boolean visible) {
		this (owner, title, info, visible, false);
	}

	// Construct object, imediately initialize if flag is set.
	// The visibility flag is ignored.
	// If init_now is true, then must be called from an application thread.
	// If init_now is false, then can be called from any thread.

	public GUICalcProgressBar (Component owner, String title, String info, boolean visible, boolean init_now) {
		this.owner = owner;
		initial_title = title;
		initial_info = info;

		progress = null;
		is_active = false;
		deque = new ArrayDeque<BarAction>();

		if (init_now) {
			req_init();
		}
	}

	// Request initialization.
	// This must be called from an application thread.
	// It waits until initialization has occurred.

	public GUICalcProgressBar req_init () {

		// Error if on the event dispatch thread

		if ( SwingUtilities.isEventDispatchThread() ) {
			throw new IllegalStateException("GUICalcProgressBar.req_init: Called while on the event dispatch thread");
		}

		// The init action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress = new CalcProgressBar(owner, initial_title, initial_info, false);
				progress.setIndeterminate(true);
				progress.setModalityType(ModalityType.APPLICATION_MODAL);
				notify_performed();
				progress.setVisible(true);		// does not return until disposed or made not-visible
				return true;					// exit from action pump
			}
		};

		// Write init action to queue

		write_queue (action, true, false);

		// Wait for it to be performed

		action.wait_performed();
		return this;
	}

	// Request disposition.

	public void req_dispose () {

		// The dispose action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress.setVisible(false);
				progress.dispose();
				progress = null;
				notify_performed();
				return false;
			}
		};

		// Write dispose action to queue

		write_queue (action, false, true);

		// Wait for it to be performed, if not on the event dispatch thread

		if (!( SwingUtilities.isEventDispatchThread()) ) {
			action.wait_performed();
		}

		return;
	}

	// Request update of title and message.

	public void req_update (final String the_title, final String the_message) {

		// The action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress.setTitle(the_title);
				progress.setProgressMessage(the_message);
				progress.pack();
				notify_performed();
				return false;
			}
		};

		// Write action to queue

		write_queue (action);
		return;
	}

	// Update the progress bar to show portion completed.
	// If total is zero, then no proportion is shown.
	// If total is non-zero, the progress bar is switched to determinate mode.
	// If message is nul or omitted, then the message is not updated.

	public void updateProgress (final long count, final long total) {
		updateProgress(count, total, null);
	}
	
	public void updateProgress (final long count, final long total, final String message) {

		// The action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress.updateProgress(count, total, message);
				notify_performed();
				return false;
			}
		};

		// Write action to queue

		write_queue (action);
		return;
	}

	// Update the info message.

	public void setProgressMessage (final String s) {

		// The action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress.setProgressMessage(s);
				notify_performed();
				return false;
			}
		};

		// Write action to queue

		write_queue (action);
		return;
	}

	// Set the indeterminate flag.
	// If a non-null message is supplied, then the info message is updated too.

	public void setIndeterminate (final boolean indep) {
		setIndeterminate (indep, null);
	}

	public void setIndeterminate (final boolean indep, final String message) {

		// The action

		BarAction action = new BarAction() {
			@Override
			public boolean run_action() {
				progress.setIndeterminate(indep);
				if (message != null) {
					progress.setProgressMessage(message);
				}
				notify_performed();
				return false;
			}
		};

		// Write action to queue

		write_queue (action);
		return;
	}


	// Closing requests disposition.

	@Override
	public void close() {
		req_dispose ();
		return;
	}

}
