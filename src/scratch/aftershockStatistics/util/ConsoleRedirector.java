package scratch.aftershockStatistics.util;

import java.io.Closeable;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.IOException;


/**
 * Redirects the console (System.out and System.err) to a file (or other destination).
 * Author: Michael Barall 08/09/2018.
 *
 * Closing this object restores the original System.out and System.err.
 */
public class ConsoleRedirector implements Closeable {

	// The destination streams, or null if not redirected, can be the same object.

	private OutputStream dest_out;
	private OutputStream dest_err;

	// True if the destination streams should be closed when this object is closed.

	private boolean f_close;

	// The original console streams, or null if not redirected.

	private PrintStream old_out;
	private PrintStream old_err;

	// The new console streams that we created, or null if not redirected.

	private PrintStream new_out;
	private PrintStream new_err;




	// Inner class that serves as the argument to the PrintStream constructor.
	// It forwards all methods to dest, except close() which it eats.
	// Methods are synchronized because the console can be used by all threads.

	private class ConRedOutputStream extends OutputStream {

		// Destination stream.

		private OutputStream dest;

		// Constructor saves the destination stream.

		public ConRedOutputStream (OutputStream dest) {
			this.dest = dest;
		}

		// Implementation of OutputStream.

		@Override
		public void write(int b) throws IOException {
			synchronized (ConsoleRedirector.this) {
				dest.write(b);
			}
			return;
		}

		@Override
		public void write(byte[] b) throws IOException {
			synchronized (ConsoleRedirector.this) {
				dest.write(b);
			}
			return;
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			synchronized (ConsoleRedirector.this) {
				dest.write(b, off, len);
			}
			return;
		}

		@Override
		public void flush() throws IOException {
			synchronized (ConsoleRedirector.this) {
				dest.flush();
			}
			return;
		}

		@Override
		public void close() throws IOException {
			return;
		}
	}




	// Inner class that serves as the argument to the PrintStream constructor.
	// It forwards all methods to both dest and tee, except close() which it eats.
	// Methods are synchronized because the console can be used by all threads.

	private class ConTeeOutputStream extends OutputStream {

		// Destination streams.

		private OutputStream dest;
		private PrintStream tee;

		// Constructor saves the destination stream.

		public ConTeeOutputStream (OutputStream dest, PrintStream tee) {
			this.dest = dest;
			this.tee = tee;
		}

		// Implementation of OutputStream.

		@Override
		public void write(int b) throws IOException {
			try {
				tee.write(b);
			} catch (Exception e) {
			}
			synchronized (ConsoleRedirector.this) {
				dest.write(b);
			}
			return;
		}

		@Override
		public void write(byte[] b) throws IOException {
			try {
				tee.write(b);
			} catch (Exception e) {
			}
			synchronized (ConsoleRedirector.this) {
				dest.write(b);
			}
			return;
		}

		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			try {
				tee.write(b, off, len);
			} catch (Exception e) {
			}
			synchronized (ConsoleRedirector.this) {
				dest.write(b, off, len);
			}
			return;
		}

		@Override
		public void flush() throws IOException {
			try {
				tee.flush();
			} catch (Exception e) {
			}
			synchronized (ConsoleRedirector.this) {
				dest.flush();
			}
			return;
		}

		@Override
		public void close() throws IOException {
			return;
		}
	}




	// Make the inner class.

	private OutputStream make_inner (OutputStream dest, PrintStream tee, boolean f_tee) {
		OutputStream result;
		if (f_tee) {
			result = new ConTeeOutputStream (dest, tee);
		} else {
			result = new ConRedOutputStream (dest);
		}
		return result;
	}




	// The close method restores the original streams, then flushes the new streams.

	@Override
	public void close() throws IOException {

		// Close exception

		Exception e_close = null;

		// Save the destination streams

		OutputStream my_dest_out = dest_out;
		OutputStream my_dest_err = dest_err;

		// If stdout is redirected ...

		if (dest_out != null) {

			// Restore the original console stream

			System.setOut (old_out);

			// Flush and close the new console stream (the flush is probably redundant)
			// (Ignore exceptions because PrintStream is not supposed to throw exceptions)

			try {
				new_out.flush();
			} catch (Exception e) {
			}

			try {
				new_out.close();
			} catch (Exception e) {
			}

			// Mark closed

			dest_out = null;
		}

		// If stderr is redirected ...

		if (dest_err != null) {

			// Restore the original console stream

			System.setErr (old_err);

			// Flush and close the new console stream (the flush is probably redundant)
			// (Ignore exceptions because PrintStream is not supposed to throw exceptions)

			try {
				new_err.flush();
			} catch (Exception e) {
			}

			try {
				new_err.close();
			} catch (Exception e) {
			}

			// Mark closed

			dest_err = null;
		}

		// If we want to close the streams ...
		// (Close must be delayed until after both are flushed, in case closing one closes the other)

		if (f_close) {

			// Close stdout destination, if it is redirected and not the same as stderr

			if (my_dest_out != null && my_dest_out != my_dest_err) {
				try {
					my_dest_out.close();
				} catch (Exception e) {
					if (e_close == null) {
						e_close = e;
					}
				}
			}

			// Close stderr destination, if it is redirected

			if (my_dest_err != null) {
				try {
					my_dest_err.close();
				} catch (Exception e) {
					if (e_close == null) {
						e_close = e;
					}
				}
			}
		}

		// If there was an exception while closing a stream, throw

		if (e_close != null) {
			throw new IOException ("ConsoleRedirector: Exception while closing stream", e_close);
		}
	
		return;
	}




	// Constructor redirects the console to the given output streams.
	// Parameters:
	//  dest_out = Destination stream for stdout, or null if no redirection is desired.
	//  dest_err = Destination stream for stderr, or null if no redirection is desired.
	//  f_close = True to close the destination streams when this object is closed.
	//  f_tee = True to write output to both destination and original streams.

	public ConsoleRedirector (OutputStream dest_out, OutputStream dest_err, boolean f_close, boolean f_tee) {

		// Save parameters and initialize fields

		this.dest_out = dest_out;
		this.dest_err = dest_err;
		this.f_close = f_close;

		old_out = null;
		old_err = null;
		new_out = null;
		new_err = null;

		boolean autoFlush = true;

		// If stdout redirection requested ...

		if (this.dest_out != null) {
	
			// Save existing console streams

			old_out = System.out;

			// Create new console streams, with auto-flush enabled

			new_out = new PrintStream (make_inner (this.dest_out, old_out, f_tee), autoFlush);

			// Install new console streams

			System.setOut (new_out);
		}

		// If stderr redirection requested ...

		if (this.dest_err != null) {
	
			// Save existing console streams

			old_err = System.err;

			// Create new console streams, with auto-flush enabled

			new_err = new PrintStream (make_inner (this.dest_err, old_err, f_tee), autoFlush);

			// Install new console streams

			System.setErr (new_err);
		}
	}




	// Constructor redirects the console to the given output stream.
	// Parameters:
	//  dest = Destination stream for stdout and stderr, or null if no redirection is desired.
	//  f_close = True to close the destination streams when this object is closed.
	//  f_tee = True to write output to both destination and original streams.

	public ConsoleRedirector (OutputStream dest, boolean f_close, boolean f_tee) {
		this (dest, dest, f_close, f_tee);
	}




	// Make a console redirector object, thereby redirecting the console.
	// Parameters:
	//  dest_out = Destination stream for stdout, or null if no redirection is desired.
	//  dest_err = Destination stream for stderr, or null if no redirection is desired.
	//  f_close = True to close the destination streams when this object is closed.
	//  f_tee = True to write output to both destination and original streams.
	// Note: You can also call the ConsoleRedirector constructor directly.
	// The difference is that this function returns null if dest_out and dest_err are both null.

	public static ConsoleRedirector make_redirector (OutputStream dest_out, OutputStream dest_err, boolean f_close, boolean f_tee) {
		if (dest_out == null && dest_err == null) {
			return null;
		}
		return new ConsoleRedirector (dest_out, dest_err, f_close, f_tee);
	}




	// Make a console redirector object, thereby redirecting the console.
	// Parameters:
	//  dest = Destination stream for stdout and stderr, or null if no redirection is desired.
	//  f_close = True to close the destination streams when this object is closed.
	//  f_tee = True to write output to both destination and original streams.
	// Note: You can also call the ConsoleRedirector constructor directly.
	// The difference is that this function returns null if dest is null.

	public static ConsoleRedirector make_redirector (OutputStream dest, boolean f_close, boolean f_tee) {
		if (dest == null) {
			return null;
		}
		return new ConsoleRedirector (dest, f_close, f_tee);
	}




	// Get the new stdout stream, if it is redirected.
	// If con_red is null, then return null.
	// If con_red is not null, and stdout is not redirected, then return null.
	// If con_red is not null, and stdout is redirected, then return the new stdout.

	public static PrintStream get_new_out (ConsoleRedirector con_red) {
		if (con_red == null) {
			return null;
		}
		return con_red.new_out;
	}




	// Get the new stderr stream, if it is redirected.
	// If con_red is null, then return null.
	// If con_red is not null, and stderr is not redirected, then return null.
	// If con_red is not null, and stderr is redirected, then return the new stderr.

	public static PrintStream get_new_err (ConsoleRedirector con_red) {
		if (con_red == null) {
			return null;
		}
		return con_red.new_err;
	}

}
