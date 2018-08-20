package scratch.aftershockStatistics.util;

import java.awt.Component;
import java.awt.Dialog.ModalityType;

import javax.swing.SwingUtilities;

import org.opensha.sha.gui.infoTools.CalcProgressBar;


/**
 * Manages a simple progress monitor that can be controlled from a worker thread.
 * Author: Michael Barall 08/18/2018.
 *
 * There are three operations:  req_init to initialize and create the progress bar,
 * req_update to update the title and progress message, and req_dispose to destroy
 * the progress bar.
 *
 * Alternatively, you can use try-with-resources.  Construct the progress manager with
 * init_now = true to create the progress bar.  Closing the progress manager destroys
 * the progress bar.  Then, use req_update to update the title and progress message.
 *
 * All calls into this class should be made from a worker thread.
 */
public class GUIProgressManager implements AutoCloseable {

	// The owner of the progress bar.

	private Component owner;

	// The progress bar.

	private CalcProgressBar progress;

	// The state: 0 = idle, 1 = init, 2 = update, 3 = dispose.

	private int state;

	public static final int STATE_IDLE = 0;
	public static final int STATE_INIT = 1;
	public static final int STATE_UPDATE = 2;
	public static final int STATE_DISPOSE = 3;

	// The title and message.

	private String title;
	private String message;

	// Construct object in idle state.

	public GUIProgressManager (Component the_owner) {
		this (the_owner, false);
	}

	// Construct object in idle state, imediately initialize if flag is set.

	public GUIProgressManager (Component the_owner, boolean init_now) {
		owner = the_owner;
		progress = null;
		state = STATE_IDLE;
		title = null;
		message = null;

		if (init_now) {
			req_init();
		}
	}

	// Request initialization.

	public void req_init () {
		synchronized (this) {
			while (state != STATE_IDLE) {
				try {
					wait();
				} catch (InterruptedException e) {}
			}
			state = STATE_INIT;
		}
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				run_next_op();
			}
		});
		return;
	}

	// Request update.

	public void req_update (String the_title, String the_message) {
		synchronized (this) {
			while (state != STATE_IDLE) {
				try {
					wait();
				} catch (InterruptedException e) {}
			}
			state = STATE_UPDATE;
			title = the_title;
			message = the_message;
		}
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				run_next_op();
			}
		});
		return;
	}

	// Request disposition.

	public void req_dispose () {
		synchronized (this) {
			while (state != STATE_IDLE) {
				try {
					wait();
				} catch (InterruptedException e) {}
			}
			state = STATE_DISPOSE;
		}
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				run_next_op();
			}
		});
		synchronized (this) {
			while (state != STATE_IDLE) {
				try {
					wait();
				} catch (InterruptedException e) {}
			}
		}
		return;
	}

	// Notify that the state is idle.

	public void notify_idle () {
		synchronized (this) {
			state = STATE_IDLE;
			notifyAll();
		}
		return;
	}

	// Run the next operation.

	public void run_next_op () {

		int my_state;
		String my_title;
		String my_message;
		synchronized (this) {
			my_state = state;
			my_title = title;
			my_message = message;
		}

		switch (my_state) {

		case STATE_INIT:
			if (progress == null) {
				progress = new CalcProgressBar(owner, "", "", false);
				progress.setIndeterminate(true);
			} else {
				if (progress.isVisible()) {
					progress.setVisible(false);
				}
			}
			progress.setModalityType(ModalityType.APPLICATION_MODAL);
			notify_idle();
			progress.setVisible(true);		// does not return until disposed or made not-visible
			break;

		case STATE_UPDATE:
			progress.setTitle(my_title);
			progress.setProgressMessage(my_message);
			progress.pack();
			notify_idle();
			break;

		case STATE_DISPOSE:
			progress.setVisible(false);
			progress.dispose();
			notify_idle();
			break;
		}

		return;
	}

	// Closing the manager requests disposition.

	@Override
	public void close() {
		req_dispose ();
		return;
	}

}
