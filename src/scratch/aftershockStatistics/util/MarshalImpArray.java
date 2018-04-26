package scratch.aftershockStatistics.util;

import java.util.Arrays;
import java.util.Set;
import java.util.HashSet;

import org.apache.commons.math3.distribution.UniformRealDistribution;

/**
 * Class for marshaling/unmarshaling parameters/data to/from array storage.
 * Author: Michael Barall 03/31/2018.
 */
public class MarshalImpArray implements MarshalReader, MarshalWriter {

	//----- Array storage -----

	// Storage of long.

	private long[] long_store;

	// Read and write indexes for long.

	private int long_read_index;
	private int long_write_index;

	// Storage of double.

	private double[] double_store;

	// Read and write indexes for double.

	private int double_read_index;
	private int double_write_index;

	// Storage of string.

	private String[] string_store;

	// Read and write indexes for string.

	private int string_read_index;
	private int string_write_index;

	//----- Context management -----

	// Class to hold current context.

	private static abstract class Context {

		// The previous context, null if this is the root context.

		protected Context previous;

		// The next context, null if this is the current context.

		protected Context next;

		// check_name - Check a name.

		public abstract void check_name (String name);

		// notify_child_begin - Notification that a child is beginning.

		public abstract void notify_child_begin (String name, Context child);

		// notify_child_end - Notification that a child is ending.

		public abstract void notify_child_end ();

		// close_map - Close a map context, return the previous context.

		public abstract Context close_map ();

		// close_array - Close an array context, return the previous context.

		public abstract Context close_array ();

		// Constructor.

		public Context (Context previous) {
			this.previous = previous;
			this.next = null;
		}
	}

	// Class to hold map context.

	private static class ContextMap extends Context {

		// Names currently in use.

		private Set<String> names;

		// Name of current child.

		private String child_name;

		// check_name - Check a name.

		@Override
		public void check_name (String name) {

			// Add the name, and throw exception if already in use

			if (name == null) {
				throw new MarshalException ("No name specified for element in map context");
			}
			if (!( names.add (name) )) {
				throw new MarshalException ("Duplicate element name in map context: name = " + name);
			}

			return;
		}

		// notify_child_begin - Notification that a child is beginning.

		@Override
		public void notify_child_begin (String name, Context child) {
			check_name (name);
			child_name = name;
			next = child;
			return;
		}

		// notify_child_end - Notification that a child is ending.

		@Override
		public void notify_child_end () {
			child_name = null;
			next = null;
			return;
		}

		// close_map - Close a map context, return the previous context.

		@Override
		public Context close_map () {
			previous.notify_child_end();
			return previous;
		}

		// close_array - Close an array context, return the previous context.

		@Override
		public Context close_array () {
			throw new MarshalException ("Attempt to end array context when in map context");
		}

		// Constructor.

		public ContextMap (String name, Context previous) {
			super (previous);
			this.names = new HashSet<String>();
			this.child_name = null;
			this.previous.notify_child_begin (name, this);
		}
	}

	// Class to hold array context.

	private static class ContextArray extends Context {

		// The array size.

		private int array_size;

		// The current index.

		private int array_index;

		// check_name - Check a name.

		@Override
		public void check_name (String name) {

			// Increment the index and check for overrun

			if (name != null) {
				throw new MarshalException ("Name specified for element in array context: name = " + name);
			}
			if (array_index == array_size) {
				throw new MarshalException ("Exceeded declared array size in array context: declared size = " + array_size);
			}
			++array_index;

			return;
		}

		// notify_child_begin - Notification that a child is beginning.

		@Override
		public void notify_child_begin (String name, Context child) {
			check_name (name);
			next = child;
			return;
		}

		// notify_child_end - Notification that a child is ending.

		@Override
		public void notify_child_end () {
			next = null;
			return;
		}

		// close_map - Close a map context, return the previous context.

		@Override
		public Context close_map () {
			throw new MarshalException ("Attempt to end map context when in array context");
		}

		// close_array - Close an array context, return the previous context.

		@Override
		public Context close_array () {
			if (array_index != array_size) {
				throw new MarshalException ("Array size mismatch in array context: declared size = " + array_size + ", actual size = " + array_index);
			}
			previous.notify_child_end();
			return previous;
		}

		// Constructor, specifies array size.

		public ContextArray (String name, Context previous, int array_size) {
			super (previous);
			if (array_size < 0) {
				throw new MarshalException ("Negative array size in array context: size = " + array_size);
			}
			this.array_size = array_size;
			this.array_index = 0;
			this.previous.notify_child_begin (name, this);
		}
	}

	// Class to hold root context.

	private static class ContextRoot extends Context {

		// True if a child of the root has been created.

		private boolean f_root_done;

		// Return true if a complete child has been processed, false if nothing processed, exception if in progress.

		public boolean get_root_status () {
			if (next != null) {
				throw new MarshalException ("Marshal/unmarshal is incomplete");
			}
			return f_root_done;
		}

		// check_name - Check a name.

		@Override
		public void check_name (String name) {

			// Throw exception

			if (name == null) {
				throw new MarshalException ("Attempt to add element in root context: name = null");
			}
			throw new MarshalException ("Attempt to add element in root context: name = " + name);
		}

		// notify_child_begin - Notification that a child is beginning.

		@Override
		public void notify_child_begin (String name, Context child) {
			if (next != null) {
				if (name == null) {
					throw new MarshalException ("Attempt to begin second child context when in root context: name = null");
				}
				throw new MarshalException ("Attempt to begin second child context when in root context: name = " + name);
			}
			if (f_root_done) {
				if (name == null) {
					throw new MarshalException ("Attempt to begin child context when in already-used root context: name = null");
				}
				throw new MarshalException ("Attempt to begin child context when in already-used root context: name = " + name);
			}
			if (name != null) {
				throw new MarshalException ("Attempt to add named child context when in root context: name = " + name);
			}
			next = child;
			return;
		}

		// notify_child_end - Notification that a child is ending.

		@Override
		public void notify_child_end () {
			if (next == null) {
				throw new MarshalException ("Attempt to end non-existent child context in root context");
			}
			next = null;
			f_root_done = true;
			return;
		}

		// close_map - Close a map context, return the previous context.

		@Override
		public Context close_map () {
			throw new MarshalException ("Attempt to end map context when in root context");
		}

		// close_array - Close an array context, return the previous context.

		@Override
		public Context close_array () {
			throw new MarshalException ("Attempt to end array context when in root context");
		}

		// Constructor.

		public ContextRoot () {
			super (null);
			f_root_done = false;
		}
	}

	// Root and current context for reading.

	private ContextRoot root_context_read;
	private Context current_context_read;

	// Root and current context for writing.

	private ContextRoot root_context_write;
	private Context current_context_write;

	//----- Implementation of MarshalReader -----

	/**
	 * Begin a map context.
	 */
	@Override
	public void unmarshalMapBegin (String name) {
		current_context_read = new ContextMap (name, current_context_read);
		return;
	}

	/**
	 * End a map context.
	 */
	@Override
	public void unmarshalMapEnd () {
		current_context_read = current_context_read.close_map();
		return;
	}

	/**
	 * Begin an array context, return the array size.
	 */
	@Override
	public int unmarshalArrayBegin (String name) {
		if (long_read_index == long_write_index) {
			throw new MarshalException ("Unmarshal long end-of-data: size = " + long_write_index);
		}
		long array_size_l = long_store[long_read_index++];
		if (array_size_l < 0L || array_size_l > (long)Integer.MAX_VALUE) {
			throw new MarshalException ("Unmarshal array size out-of-range: size = " + array_size_l);
		}
		int array_size = (int)array_size_l;
		current_context_read = new ContextArray (name, current_context_read, array_size);
		return array_size;
	}

	/**
	 * End an array context.
	 */
	@Override
	public void unmarshalArrayEnd () {
		current_context_read = current_context_read.close_array();
		return;
	}

	/**
	 * Unmarshal a long.
	 */
	@Override
	public long unmarshalLong (String name) {
		current_context_read.check_name (name);
		if (long_read_index == long_write_index) {
			throw new MarshalException ("Unmarshal long end-of-data: size = " + long_write_index);
		}
		return long_store[long_read_index++];
	}

	/**
	 * Unmarshal a double.
	 */
	@Override
	public double unmarshalDouble (String name) {
		current_context_read.check_name (name);
		if (double_read_index == double_write_index) {
			throw new MarshalException ("Unmarshal double end-of-data: size = " + double_write_index);
		}
		return double_store[double_read_index++];
	}

	/**
	 * Unmarshal a string.  (Null strings are not allowed.)
	 */
	@Override
	public String unmarshalString (String name) {
		current_context_read.check_name (name);
		if (string_read_index == string_write_index) {
			throw new MarshalException ("Unmarshal string end-of-data: size = " + string_write_index);
		}
		return string_store[string_read_index++];
	}

	//----- Implementation of MarshalWriter -----

	/**
	 * Begin a map context.
	 */
	@Override
	public void marshalMapBegin (String name) {
		current_context_write = new ContextMap (name, current_context_write);
		return;
	}

	/**
	 * End a map context.
	 */
	@Override
	public void marshalMapEnd () {
		current_context_write = current_context_write.close_map();
		return;
	}

	/**
	 * Begin an array context, specify the array size.
	 */
	@Override
	public void marshalArrayBegin (String name, int array_size) {
		if (long_write_index == long_store.length) {
			int new_capacity = Math.max (100, long_store.length * 2);
			long_store = Arrays.copyOf (long_store, new_capacity);
		}
		long_store[long_write_index++] = (long)array_size;
		current_context_write = new ContextArray (name, current_context_write, array_size);
		return;
	}

	/**
	 * End an array context.
	 */
	@Override
	public void marshalArrayEnd () {
		current_context_write = current_context_write.close_array();
		return;
	}

	/**
	 * Marshal a long.
	 */
	@Override
	public void marshalLong (String name, long x) {
		current_context_write.check_name (name);
		if (long_write_index == long_store.length) {
			int new_capacity = Math.max (100, long_store.length * 2);
			long_store = Arrays.copyOf (long_store, new_capacity);
		}
		long_store[long_write_index++] = x;
		return;
	}

	/**
	 * Marshal a double.
	 */
	@Override
	public void marshalDouble (String name, double x) {
		current_context_write.check_name (name);
		if (double_write_index == double_store.length) {
			int new_capacity = Math.max (100, double_store.length * 2);
			double_store = Arrays.copyOf (double_store, new_capacity);
		}
		double_store[double_write_index++] = x;
		return;
	}

	/**
	 * Marshal a string.  (Null strings are not allowed.)
	 */
	@Override
	public void marshalString (String name, String x) {
		current_context_write.check_name (name);
		if (string_write_index == string_store.length) {
			int new_capacity = Math.max (100, string_store.length * 2);
			string_store = Arrays.copyOf (string_store, new_capacity);
		}
		string_store[string_write_index++] = x;
		return;
	}

	//----- Construction -----

	/**
	 * Create an empty object, suitable for writing.
	 */
	public MarshalImpArray () {
		long_store = new long[100];
		long_read_index = 0;
		long_write_index = 0;

		double_store = new double[100];
		double_read_index = 0;
		double_write_index = 0;

		string_store = new String[100];
		string_read_index = 0;
		string_write_index = 0;

		root_context_read = new ContextRoot();
		current_context_read = root_context_read;

		root_context_write = new ContextRoot();
		current_context_write = root_context_write;
	}

	/**
	 * Create an object initialized with the given arrays, suitable for reading.
	 */
	public MarshalImpArray (long[] long_store, double[] double_store, String[] string_store) {
		this.long_store = long_store;
		long_read_index = 0;
		long_write_index = this.long_store.length;

		this.double_store = double_store;
		double_read_index = 0;
		double_write_index = this.double_store.length;

		this.string_store = string_store;
		string_read_index = 0;
		string_write_index = this.string_store.length;

		root_context_read = new ContextRoot();
		current_context_read = root_context_read;

		root_context_write = new ContextRoot();
		current_context_write = root_context_write;
	}

	//----- Control -----

	/**
	 * Set indexes to read from beginning of array.
	 */
	public void begin_read () {
		long_read_index = 0;
		double_read_index = 0;
		string_read_index = 0;

		root_context_read = new ContextRoot();
		current_context_read = root_context_read;
		return;
	}

	/**
	 * Set indexes to write to beginning of array.
	 */
	public void begin_write () {
		long_read_index = 0;
		long_write_index = 0;
		double_read_index = 0;
		double_write_index = 0;
		string_read_index = 0;
		string_write_index = 0;

		root_context_read = new ContextRoot();
		current_context_read = root_context_read;

		root_context_write = new ContextRoot();
		current_context_write = root_context_write;
		return;
	}

	/**
	 * Check read status, return true if read complete, false if nothing read, exception if in progress.
	 */
	public boolean check_read_complete () {
		return root_context_read.get_root_status();
	}

	/**
	 * Check write status, return true if write complete, false if nothing written, exception if in progress.
	 */
	public boolean check_write_complete () {
		return root_context_write.get_root_status();
	}

	/**
	 * Get the long store.
	 * Returns a copy, with length equal to the amount of data, but at least one element.
	 */
	public long[] get_long_store () {
		if (long_write_index == 0) {
			long[] result = new long[1];
			result[0] = 0L;
			return result;
		}
		return Arrays.copyOf (long_store, long_write_index);
	}

	/**
	 * Get the double store.
	 * Returns a copy, with length equal to the amount of data, but at least one element.
	 */
	public double[] get_double_store () {
		if (double_write_index == 0) {
			double[] result = new double[1];
			result[0] = 0.0;
			return result;
		}
		return Arrays.copyOf (double_store, double_write_index);
	}

	/**
	 * Get the string store.
	 * Returns a copy, with length equal to the amount of data, but at least one element.
	 */
	public String[] get_string_store () {
		if (string_write_index == 0) {
			String[] result = new String[1];
			result[0] = "";
			return result;
		}
		return Arrays.copyOf (string_store, string_write_index);
	}




	//----- Testing -----




	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("MarshalImpArray : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  num_long  num_double  num_string
		// Test marshaling and unmarshaling with the given numbers of long and double and string.
		// Data values are randomly generated.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Two additional arguments

			if (args.length != 4) {
				System.err.println ("MarshalImpArray : Invalid 'test1' subcommand");
				return;
			}
			int num_long = Integer.parseInt(args[1]);
			int num_double = Integer.parseInt(args[2]);
			int num_string = Integer.parseInt(args[3]);

			System.out.println (
				"num_long = " + num_long + "\n" +
				"num_double = " + num_double + "\n" +
				"num_string = " + num_string
			);

			// Random number generator

			UniformRealDistribution rangen = new UniformRealDistribution();

			// Generate random values

			System.out.println ("Generating random data ...");

			long[] long_data = new long[num_long];
			for (int i = 0; i < num_long; ++i) {
				long_data[i] = Math.round (rangen.sample() * 1.0e12);
				if (i < 10) {
					System.out.println ("long_data[" + i + "] = " + long_data[i]);
				}
			}

			double[] double_data = new double[num_double];
			for (int i = 0; i < num_double; ++i) {
				double_data[i] = rangen.sample() * 1.0e12;
				if (i < 10) {
					System.out.println ("double_data[" + i + "] = " + double_data[i]);
				}
			}

			String[] string_data = new String[num_string];
			for (int i = 0; i < num_string; ++i) {
				string_data[i] = "String" + Math.round (rangen.sample() * 1.0e12);
				if (i < 10) {
					System.out.println ("string_data[" + i + "] = " + string_data[i]);
				}
			}

			// Marshal the data

			System.out.println ("Marshaling data ...");

			MarshalImpArray writer = new MarshalImpArray();

			writer.marshalArrayBegin (null, num_long + num_double + num_string);

			for (int i = 0; i < num_long; ++i) {
				writer.marshalLong (null, long_data[i]);
			}

			for (int i = 0; i < num_double; ++i) {
				writer.marshalDouble (null, double_data[i]);
			}

			for (int i = 0; i < num_string; ++i) {
				writer.marshalString (null, string_data[i]);
			}

			writer.marshalArrayEnd ();

			if (!( writer.check_write_complete() )) {
				System.out.println ("Writer reports writing not complete");
				return;
			}

			long[] m_long_store = writer.get_long_store();
			double[] m_double_store = writer.get_double_store();
			String[] m_string_store = writer.get_string_store();

			writer = null;

			for (int i = 0; i < 11 && i < num_long + 1; ++i) {
				System.out.println ("m_long_store[" + i + "] = " + m_long_store[i]);
			}

			for (int i = 0; i < 10 && i < num_double; ++i) {
				System.out.println ("m_double_store[" + i + "] = " + m_double_store[i]);
			}

			for (int i = 0; i < 10 && i < num_string; ++i) {
				System.out.println ("m_string_store[" + i + "] = " + m_string_store[i]);
			}

			// Unmarshal and check the data

			System.out.println ("Unmarshaling data ...");

			MarshalImpArray reader = new MarshalImpArray (m_long_store, m_double_store, m_string_store);

			int array_size = reader.unmarshalArrayBegin (null);
			if (array_size != num_long + num_double + num_string) {
				System.out.println ("Reader reports incorrect array size: " + array_size);
				return;
			}

			int errors = 0;

			for (int i = 0; i < num_long; ++i) {
				long x = reader.unmarshalLong (null);
				if (x != long_data[i]) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched long: i = " + i + ", d = " + long_data[i] + ", x = " + x);
					}
				}
			}

			for (int i = 0; i < num_double; ++i) {
				double x = reader.unmarshalDouble (null);
				if (x != double_data[i]) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched double: i = " + i + ", d = " + double_data[i] + ", x = " + x);
					}
				}
			}

			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalString (null);
				if (!( x.equals(string_data[i]) )) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched string: i = " + i + ", d = " + string_data[i] + ", x = " + x);
					}
				}
			}

			reader.unmarshalArrayEnd ();

			if (!( reader.check_read_complete() )) {
				System.out.println ("Reader reports reading not complete");
				return;
			}

			System.out.println ("Error count: " + errors);

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("MarshalImpArray : Unrecognized subcommand : " + args[0]);
		return;

	}




}
