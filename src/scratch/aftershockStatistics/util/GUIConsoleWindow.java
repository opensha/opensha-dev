package scratch.aftershockStatistics.util;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import java.util.ArrayDeque;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

/**
 * Redirects the console (System.out and System.err) to a JTextArea.
 * Author: Michael Barall 08/18/2018.
 *
 * This class is a replacement for org.opensha.commons.gui.ConsoleWindow.
 * The problem with ConsoleWindow is that it writes to the JTextArea any
 * time System.out or System.err is called, regardless of whether or not
 * it is on the event dispatch thread.  This class stores output in a queue,
 * and writes to the JTextArea only on the event dispatch thread.
 */
public class GUIConsoleWindow {
	
	private JDialog frame;
	private JTextArea text = new JTextArea();
	private JScrollPane scroll = new JScrollPane(text);

	private ArrayDeque<String> deque = new ArrayDeque<String>();
	
	public GUIConsoleWindow() {
		this(false);
	}
	
	public GUIConsoleWindow(boolean noFrame) {
		System.setErr(new PrintStream (new GUIOutputStream (System.err), true));
		System.setOut(new PrintStream (new GUIOutputStream (System.out), true));
		initGUI(noFrame);
	}
	
	public void initGUI(boolean noFrame) {
		if (!noFrame) {
			frame = new JDialog(new JFrame(), "Console Window");
			frame.setLocationRelativeTo(null);
			frame.setSize(800,500);
			frame.add(scroll);
		}
		text.setEditable(false);
		return;
	}
	
	public JScrollPane getScrollPane() {
		return scroll;
	}
	
	public JTextArea getTextArea() {
		return text;
	}
	
	public void setVisible(boolean show) {
		if (frame != null) {
			frame.setLocationRelativeTo(null);
			text.setCaretPosition(0);
			text.setCaretPosition(text.getText().length());
			frame.setVisible(show);
		}
		return;
	}

	// Write the entire contents of the deque to the window.
	// This method must be called on the event dispatcher thread!

	private void write_to_window () {
		for (;;) {
			String s;
			synchronized (this) {
				s = deque.pollFirst();
			}
			if (s == null) {
				break;
			}
			text.append(s);
			text.setCaretPosition(text.getText().length());
		}
		return;
	}




	// Inner class that serves as the argument to the PrintStream constructor.
	// It forwards all methods to both text and tee, except close() which it eats.
	// Data is sent to text through the deque.
	// Methods are synchronized because the console can be used by all threads.

	private class GUIOutputStream extends OutputStream {

		// Destination stream.

		private PrintStream tee;

		// Constructor.

		public GUIOutputStream (PrintStream tee) {
			this.tee = tee;
		}

		// Implementation of OutputStream.

		@Override
		public void write(int b) throws IOException {
			try {
				tee.write(b);
			} catch (Exception e) {
			}
			synchronized (GUIConsoleWindow.this) {
				deque.addLast(new String (new byte[]{(byte)b}));
			}
			return;
		}

		@Override
		public void write(byte[] b) throws IOException {
			try {
				tee.write(b);
			} catch (Exception e) {
			}
			synchronized (GUIConsoleWindow.this) {
				deque.addLast(new String (b));
			}
			return;
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			try {
				tee.write(b, off, len);
			} catch (Exception e) {
			}
			synchronized (GUIConsoleWindow.this) {
				deque.addLast(new String (b, off, len));
			}
			return;
		}

		@Override
		public void flush() throws IOException {
			try {
				tee.flush();
			} catch (Exception e) {
			}
			boolean empty;
			synchronized (GUIConsoleWindow.this) {
				empty = deque.isEmpty();
			}
			if (!( empty )) {
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {
						write_to_window();
					}
				});
			}
			return;
		}

		@Override
		public void close() throws IOException {
			return;
		}
	}

}
