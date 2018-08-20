package scratch.aftershockStatistics.util;

import java.io.Closeable;
import java.io.Flushable;
import java.io.OutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.BufferedOutputStream;

import java.nio.file.Paths;
import java.nio.file.Path;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

import java.text.SimpleDateFormat;

import java.util.Date;
import java.util.TimeZone;
import java.util.ArrayList;
import java.util.LinkedHashSet;


/**
 * An output stream that is split into multiple files based on the current time.
 * Author: Michael Barall 08/10/2018.
 *
 */
public class TimeSplitOutputStream extends FilterOutputStream {

	// The current output stream, or null if none.

	//protected OutputStream out;			// inherited from FilterOutputStream

	// The pattern used to generate a filename, in the format of SimpleDateFormat.

	protected String filename_pattern;

	// The current filename to use for the output stream.

	protected String current_filename;

	// True if this object has been closed.

	protected boolean is_closed;

	// True if an attempt has been made to open the file.

	protected boolean f_open_attempted;

	// List of upstream objects to flush when underlying file is closed.

	protected LinkedHashSet<Flushable> upstream_list;

	// Formatter used to convert times to filenames

	protected SimpleDateFormat fmt;




	// Lazy-open the file.

	protected void lazy_open () throws IOException {

		// If this object is already closed, fail

		if (is_closed) {
			throw new IOException ("TimeSplitOutputStream: Attempt to write after stream is closed");
		}

		// If open already attempted, fail

		if (f_open_attempted) {
			throw new IOException ("TimeSplitOutputStream: Failed to open file: " + current_filename);
		}

		// Indicate open has been attempted

		f_open_attempted = true;

		// Get a path for the current filename, and its containing directory

		Path file_path;
		Path dir_path;

		try {
			file_path = Paths.get (current_filename);
			dir_path = file_path.getParent();
		} catch (Exception e) {
			throw new IOException ("TimeSplitOutputStream: Filename does not represent a valid filesystem path: " + current_filename, e);
		}

		// If there is a containing directory, create it if it doesn't already exist

		if (dir_path != null) {
			try {
				Files.createDirectories (dir_path);
			} catch (Exception e) {
				throw new IOException ("TimeSplitOutputStream: Unable to create parent directory for file: " + current_filename, e);
			}
		}

		// Create the stream, appending to an existing file, creating file if it does not exist

		OutputStream stream;

		try {
			stream = Files.newOutputStream (file_path, StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		} catch (Exception e) {
			throw new IOException ("TimeSplitOutputStream: Unable to open file: " + current_filename, e);
		}

		// Establish the stream

		out = stream;
		return;
	}




	// Flush all upstream objects.

	protected void flush_all_upstream () throws IOException {
		Exception e_flush = null;
		for (Flushable upstream : upstream_list) {
			try {
				upstream.flush();
			} catch (Exception e) {
				if (e_flush == null) {
					e_flush = e;
				}
			}
		}
		if (e_flush != null) {
			throw new IOException ("TimeSplitOutputStream: Exception while flushing upstream objects", e_flush);
		}
		return;
	}




	// Add an upstream object.
	// If upstream is null, then perform no operation.

	public void add_upstream (Flushable upstream) {
		if (upstream != null) {
			upstream_list.add (upstream);
		}
		return;
	}




	// Remove an upstream object, if it is currently in the list.
	// If upstream is null, then perform no operation.

	public void remove_upstream (Flushable upstream) {
		if (upstream != null) {
			upstream_list.remove (upstream);
		}
		return;
	}



	// This closeable inner class removes an upstream object when it is closed.

	protected class AutoRemoveUpstream implements Closeable {

		// The upstream object.

		private Flushable upstream;

		// Constructor.

		public AutoRemoveUpstream (Flushable upstream) {
			this.upstream = upstream;
		}

		// Remove the upstream object from the outer class.

		@Override
		public void close() throws IOException {
			if (upstream != null) {
				remove_upstream (upstream);
			}
			upstream = null;
			return;
		}
	}




	// Add an upstream object.
	// Returns a closeable object, which removes the upstream object when closed.
	// If upstream is null, then perform no operation and return null.
	// Note: This can be used in a try-with-resources statement to add an upstream
	// object immediately after it is created.

	public Closeable add_auto_upstream (Flushable upstream) {
		if (upstream == null) {
			return null;
		}
		upstream_list.add (upstream);
		return new AutoRemoveUpstream (upstream);
	}




	// Make the filename, given the currnet time (in milliseconds since the epoch).

	protected String make_filename (long the_time) {
		//SimpleDateFormat fmt = new SimpleDateFormat (filename_pattern);
		//fmt.setTimeZone (TimeZone.getTimeZone ("UTC"));
		return fmt.format (new Date (the_time));
	}




	// Redirect destination file, if necessary.
	// Parameters:
	//  the_time = Current time, in milliseconds since the epoch.

	public void redirect (long the_time) throws IOException {

		// If this object is already closed, fail

		if (is_closed) {
			throw new IOException ("TimeSplitOutputStream: Attempt to redirect after stream is closed");
		}

		// Get the new filename

		String new_filename = make_filename (the_time);

		// If filename is changed ...

		if (!( new_filename.equals (current_filename) )) {

			// An exception obtained while flushing, null if none

			Exception e_flush = null;

			// Flush all upstream objects
			
			for (Flushable upstream : upstream_list) {
				try {
					upstream.flush();
				} catch (Exception e) {
					if (e_flush == null) {
						e_flush = e;
					}
				}
			}

			// Flush the downstream object

			try {
				flush();
			} catch (Exception e) {
				if (e_flush == null) {
					e_flush = e;
				}
			}

			// Save the stream, then null it

			OutputStream stream = out;
			out = null;

			// Insert the new filename

			current_filename = new_filename;
			f_open_attempted = false;

			// Close the downstream object

			if (stream != null) {
				stream.close();
			}

			// If there was an exception during flushing, throw it now

			if (e_flush != null) {
				throw new IOException ("TimeSplitOutputStream: Exception while performing flush operation", e_flush);
			}
		}
	
		return;
	}




	// Write data to the stream.

	@Override
	public void write(int b) throws IOException {
		if (out == null) {
			lazy_open ();
		}
		out.write(b);
		return;
	}

	@Override
	public void write(byte[] b) throws IOException {
		if (out == null) {
			lazy_open ();
		}
		out.write(b);
		return;
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		if (out == null) {
			lazy_open ();
		}
		out.write(b, off, len);
		return;
	}




	// Flush the stream.

	@Override
	public void flush() throws IOException {

		// If this object is already closed, fail

		if (is_closed) {
			throw new IOException ("TimeSplitOutputStream: Attempt to flush after stream is closed");
		}

		// If a file is open, flush it
		// Note: This function must not flush upstream objects, because the flush functions
		// of the upstream objects are likely to call this function, creating a circularity.

		if (out != null) {
			out.flush();
		}
		return;
	}




	// Close the stream.

	@Override
	public void close() throws IOException {

		// If not already closed ...

		if (!( is_closed )) {

			// An exception obtained while flushing, null if none

			Exception e_flush = null;

			// Note: This function must not flush upstream objects, because the upstream
			// objects are likely to be already closed.

			try {
				flush();
			} catch (Exception e) {
				if (e_flush == null) {
					e_flush = e;
				}
			}

			// Save the stream, then null it

			OutputStream stream = out;
			out = null;
		
			// Mark it closed

			is_closed = true;

			// Close the downstream object

			if (stream != null) {
				stream.close();
			}

			// If there was an exception during flushing, throw it now

			if (e_flush != null) {
				throw new IOException ("TimeSplitOutputStream: Exception while performing flush operation", e_flush);
			}
		}

		return;
	}




	// Constructor.
	// Parameters:
	//  the_filename_pattern = The pattern used to generate a filename, in the format of SimpleDateFormat.
	//  the_time = Current time, in milliseconds since the epoch.

	public TimeSplitOutputStream (String the_filename_pattern, long the_time) {

		super (null);		// initializes out = null

		// Check parameters

		if (the_filename_pattern == null || the_filename_pattern.isEmpty()) {
			throw new IllegalArgumentException ("TimeSplitOutputStream: No filename pattern");
		}

		// Set up the formatter

		try {
			fmt = new SimpleDateFormat (the_filename_pattern);
			fmt.setTimeZone (TimeZone.getTimeZone ("UTC"));
		} catch (Exception e) {
			throw new IllegalArgumentException ("TimeSplitOutputStream: Bad filename pattern: " + filename_pattern, e);
		}

		// Save parameters and initialize fields

		out = null;
		filename_pattern = the_filename_pattern;
		current_filename = make_filename (the_time);
		is_closed = false;
		f_open_attempted = false;
		upstream_list = new LinkedHashSet<Flushable>();
	}




	// Make a time split output stream.
	// Parameters:
	//  the_filename_pattern = The pattern used to generate a filename, in the format of SimpleDateFormat.
	//  the_time = Current time, in milliseconds since the epoch.
	// Note: You can also call the TimeSplitOutputStream constructor directly.
	// The difference is that this function returns null if the filename pattern is null or empty.

	public static TimeSplitOutputStream make_tsop (String the_filename_pattern, long the_time) {
		if (the_filename_pattern == null || the_filename_pattern.isEmpty()) {
			return null;
		}
		return new TimeSplitOutputStream (the_filename_pattern, the_time);
	}




	// Add an upstream object.
	// Returns a closeable object, which removes the upstream object when closed.
	// If tsop or upstream is null, then perform no operation and return null.
	// Note: This can be used in a try-with-resources statement to add an upstream
	// object immediately after it is created.

	public static Closeable add_auto_upstream (TimeSplitOutputStream tsop, Flushable upstream) {
		if (tsop == null || upstream == null) {
			return null;
		}
		return tsop.add_auto_upstream (upstream);
	}



	// Test if a filename pattern is valid.
	// Note: This is a basic test examining only the text of the pattern,
	// to see if it acceptable to SimpleDateFormat.
	// It cannot determine if files can actually be created.

	public static boolean is_valid_pattern (String the_filename_pattern) {
		if (the_filename_pattern == null || the_filename_pattern.isEmpty()) {
			return false;
		}
		try {
			SimpleDateFormat the_fmt = new SimpleDateFormat (the_filename_pattern);
			the_fmt.setTimeZone (TimeZone.getTimeZone ("UTC"));
		} catch (Exception e) {
			return false;
		}
		return true;
	}

}
