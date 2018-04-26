package scratch.aftershockStatistics.util;

import java.util.Arrays;
import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import java.io.Writer;
import java.io.Reader;
import java.io.IOException;

import org.json.simple.JSONAware;
import org.json.simple.JSONStreamAware;
import org.json.simple.ItemList;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.JSONValue;

import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.ContainerFactory;
import org.json.simple.parser.ContentHandler;
import org.json.simple.parser.Yytoken;

import org.apache.commons.math3.distribution.UniformRealDistribution;

/**
 * Class for marshaling parameters/data to JSON data structures.
 * Author: Michael Barall 04/06/2018.
 */
public class MarshalImpJsonWriter implements MarshalWriter {

	//----- Context management -----

	// Class to hold current context.

	private static abstract class Context {

		// The previous context, null if this is the root context.

		protected Context previous;

		// The next context, null if this is the current context.

		protected Context next;

		// check_name - Check a name, and supply the named object.

		public abstract void check_name (String name, Object o);

		// notify_child_begin - Notification that a child is beginning.

		public abstract void notify_child_begin (String name, Context child);

		// notify_child_end - Notification that a child is ending, and supply the child object.

		public abstract void notify_child_end (Object o);

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

		// The map.

		private JSONObject json_map;

		// Names currently in use.

		private Set<String> names;

		// Name of current child.

		private String child_name;

		// check_name - Check a name, and supply the named object.

		@Override
		public void check_name (String name, Object o) {

			// Add the name, and throw exception if already in use

			if (name == null) {
				throw new MarshalException ("No name specified for element in map context");
			}
			if (!( names.add (name) )) {
				throw new MarshalException ("Duplicate element name in map context: name = " + name);
			}
			json_map.put (name, o);

			return;
		}

		// notify_child_begin - Notification that a child is beginning.

		@Override
		public void notify_child_begin (String name, Context child) {
			child_name = name;
			next = child;
			return;
		}

		// notify_child_end - Notification that a child is ending, and supply the child object.

		@Override
		public void notify_child_end (Object o) {
			check_name (child_name, o);
			child_name = null;
			next = null;
			return;
		}

		// close_map - Close a map context, return the previous context.

		@Override
		public Context close_map () {
			previous.notify_child_end (json_map);
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
			this.json_map = new JSONObject();
			this.names = new HashSet<String>();
			this.child_name = null;
			this.previous.notify_child_begin (name, this);
		}
	}

	// Class to hold array context.

	private static class ContextArray extends Context {

		// The array.

		private JSONArray json_array;

		// The array size.

		private int array_size;

		// The current index.

		private int array_index;

		// check_name - Check a name, and supply the named object.

		@Override
		public void check_name (String name, Object o) {

			// Increment the index and check for overrun

			if (name != null) {
				throw new MarshalException ("Name specified for element in array context: name = " + name);
			}
			if (array_index == array_size) {
				throw new MarshalException ("Exceeded declared array size in array context: declared size = " + array_size);
			}
			json_array.add (o);
			++array_index;

			return;
		}

		// notify_child_begin - Notification that a child is beginning.

		@Override
		public void notify_child_begin (String name, Context child) {
			if (name != null) {
				throw new MarshalException ("Name specified for element in array context: name = " + name);
			}
			next = child;
			return;
		}

		// notify_child_end - Notification that a child is ending, and supply the child object.

		@Override
		public void notify_child_end (Object o) {
			check_name (null, o);
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
			previous.notify_child_end (json_array);
			return previous;
		}

		// Constructor, specifies array size.

		public ContextArray (String name, Context previous, int array_size) {
			super (previous);
			if (array_size < 0) {
				throw new MarshalException ("Negative array size in array context: size = " + array_size);
			}
			this.json_array = new JSONArray();
			this.array_size = array_size;
			this.array_index = 0;
			this.previous.notify_child_begin (name, this);
		}
	}

	// Class to hold root context.

	private static class ContextRoot extends Context {

		// The JSON container, can be either JSONObject or JSONArray, or null.

		private Object json_container;

		// True if a child of the root has been created.

		private boolean f_root_done;

		// Return true if a complete child has been processed, false if nothing processed, exception if in progress.

		public boolean get_root_status () {
			if (next != null) {
				throw new MarshalException ("Marshal is incomplete");
			}
			return f_root_done;
		}

		// Get the JSON container.

		public Object get_json_container () {
			return json_container;
		}

		// check_name - Check a name, and supply the named object.

		@Override
		public void check_name (String name, Object o) {

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

		// notify_child_end - Notification that a child is ending, and supply the child object.

		@Override
		public void notify_child_end (Object o) {
			if (next == null) {
				throw new MarshalException ("Attempt to end non-existent child context in root context");
			}
			next = null;
			json_container = o;
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
			json_container = null;
			f_root_done = false;
		}
	}

	// Root and current context for writing.

	private ContextRoot root_context_write;
	private Context current_context_write;

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
		current_context_write.check_name (name, new Long(x));
		return;
	}

	/**
	 * Marshal a double.
	 */
	@Override
	public void marshalDouble (String name, double x) {
		current_context_write.check_name (name, new Double(x));
		return;
	}

	/**
	 * Marshal a string.  (Null strings are not allowed.)
	 */
	@Override
	public void marshalString (String name, String x) {
		current_context_write.check_name (name, x);
		return;
	}

	/**
	 * Marshal a boolean.
	 */
	@Override
	public void marshalBoolean (String name, boolean x) {
		current_context_write.check_name (name, new Boolean(x));
		return;
	}

	/**
	 * Marshal a JSON string.  (Null strings are not allowed.)
	 * The string must contain a JSON object or array, or be an empty string.
	 * For JSON storage, the string is merged into the JSON instead of being
	 * embedded as string-valued data.  (An empty string becomes a JSON null.)
	 * The unmarshaled string may differ from the marshaled string due to JSON parsing.
	 * (Named element ordering, numeric formats, and spacing may be changed).
	 */
	@Override
	public void marshalJsonString (String name, String x) {
		Object json_container;

		if (x.equals("")) {
			json_container = null;
		}
		else {
			try {
				json_container = JSONValue.parseWithException (x);
			}
			catch (ParseException e) {
				throw new MarshalException ("Parsing error while parsing JSON string: name = " + ((name == null) ? "null" : name), e);
			}
			catch (Exception e) {
				throw new MarshalException ("Exception while parsing JSON string: name = " + ((name == null) ? "null" : name), e);
			}
			if (!( json_container instanceof JSONArray || json_container instanceof JSONObject )) {
				throw new MarshalException ("JSON String does not contain a JSON object or JSON array: name = " + ((name == null) ? "null" : name));
			}
		}

		current_context_write.check_name (name, json_container);
		return;
	}

	//----- Construction -----

	/**
	 * Create an empty object, suitable for writing.
	 */
	public MarshalImpJsonWriter () {
		root_context_write = new ContextRoot();
		current_context_write = root_context_write;
	}

	//----- Control -----

	/**
	 * Check write status, return true if write complete, false if nothing written, exception if in progress.
	 */
	public boolean check_write_complete () {
		return root_context_write.get_root_status();
	}

	/**
	 * Get the JSON container.
	 * It can be either JSONObject or JSONArray, or null.
	 */
	public Object get_json_container () {
		return root_context_write.get_json_container();
	}

	/**
	 * Get the JSON container, converted to a String.
	 */
	public String get_json_string () {
		Object json_container = root_context_write.get_json_container();
		String result;
		try {
			result = JSONValue.toJSONString(json_container);
		}
		catch (Exception e) {
			throw new MarshalException ("Exception while writing JSON string", e);
		}
		return result;
	}

	/**
	 * Write the JSON container to a file.
	 */
	public void write_json_file (Writer out) {
		Object json_container = root_context_write.get_json_container();
		try {
			JSONValue.writeJSONString(json_container, out);
		}
		catch (IOException e) {
			throw new MarshalException ("I/O error while writing JSON file", e);
		}
		catch (Exception e) {
			throw new MarshalException ("Exception while writing JSON string", e);
		}
		return;
	}




	//----- Testing -----




	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("MarshalImpJsonWriter : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  num_long  num_double  num_string
		// Test marshaling and unmarshaling with the given numbers of long and double and string.
		// Data values are randomly generated.

		if (args[0].equalsIgnoreCase ("test1")) {

			// Three additional arguments

			if (args.length != 4) {
				System.err.println ("MarshalImpJsonWriter : Invalid 'test1' subcommand");
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
				switch (i % 5) {
				case 0: double_data[i] = rangen.sample() * 1.0e12; break;
				case 1: double_data[i] = rangen.sample() * 2.0e7; break;
				case 2: double_data[i] = rangen.sample() * 3.0e0; break;
				case 3: double_data[i] = rangen.sample() * 2.0e-3; break;
				case 4: double_data[i] = rangen.sample() * 1.0e-12; break;
				}
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

			MarshalImpJsonWriter writer = new MarshalImpJsonWriter();

			writer.marshalArrayBegin (null, num_long + num_double + num_string + 2);

			writer.marshalBoolean (null, true);
			writer.marshalBoolean (null, false);

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

			Object json_container = writer.get_json_container();

			if (num_long + num_double + num_string <= 100) {
				System.out.println (writer.get_json_string());
			}

			writer = null;

			// Unmarshal and check the data

			System.out.println ("Unmarshaling data ...");

			MarshalImpJsonReader reader = new MarshalImpJsonReader (json_container);

			int array_size = reader.unmarshalArrayBegin (null);
			if (array_size != num_long + num_double + num_string + 2) {
				System.out.println ("Reader reports incorrect array size: " + array_size);
				return;
			}

			int errors = 0;

			if (!( reader.unmarshalBoolean (null) )) {
				++errors;
				System.out.println ("Mismatched boolean: expecting true, got false");
			}

			if ( reader.unmarshalBoolean (null) ) {
				++errors;
				System.out.println ("Mismatched boolean: expecting false, got true");
			}

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




		// Subcommand : Test #2
		// Command format:
		//  test2  num_long  num_double  num_string
		// Test marshaling and unmarshaling with the given numbers of long and double and string.
		// Data values are randomly generated.
		// Same as test #1 except data is passed from writer to reader thru a JSON string.

		if (args[0].equalsIgnoreCase ("test2")) {

			// Three additional arguments

			if (args.length != 4) {
				System.err.println ("MarshalImpJsonWriter : Invalid 'test2' subcommand");
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
				switch (i % 5) {
				case 0: double_data[i] = rangen.sample() * 1.0e12; break;
				case 1: double_data[i] = rangen.sample() * 2.0e7; break;
				case 2: double_data[i] = rangen.sample() * 3.0e0; break;
				case 3: double_data[i] = rangen.sample() * 2.0e-3; break;
				case 4: double_data[i] = rangen.sample() * 1.0e-12; break;
				}
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

			MarshalImpJsonWriter writer = new MarshalImpJsonWriter();

			writer.marshalArrayBegin (null, num_long + num_double + num_string + 2);

			writer.marshalBoolean (null, true);
			writer.marshalBoolean (null, false);

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

			String json_string = writer.get_json_string();

			if (num_long + num_double + num_string <= 100) {
				System.out.println (json_string);
			}

			writer = null;

			// Unmarshal and check the data

			System.out.println ("Unmarshaling data ...");

			MarshalImpJsonReader reader = new MarshalImpJsonReader (json_string);

			int array_size = reader.unmarshalArrayBegin (null);
			if (array_size != num_long + num_double + num_string + 2) {
				System.out.println ("Reader reports incorrect array size: " + array_size);
				return;
			}

			int errors = 0;

			if (!( reader.unmarshalBoolean (null) )) {
				++errors;
				System.out.println ("Mismatched boolean: expecting true, got false");
			}

			if ( reader.unmarshalBoolean (null) ) {
				++errors;
				System.out.println ("Mismatched boolean: expecting false, got true");
			}

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




		// Subcommand : Test #3
		// Command format:
		//  test3
		// Test marshaling and unmarshaling JSON strings.

		if (args[0].equalsIgnoreCase ("test3")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("MarshalImpJsonWriter : Invalid 'test3' subcommand");
				return;
			}

			// Set up data

			int num_string = 4;
			String[] string_data = new String[num_string];

			string_data[0] = "{\"one\":1, \"two\":2, \"three\":3, \"four\":4}";
			string_data[1] = "[5, 6, 7, 8, 9]";
			string_data[2] = "";
			string_data[3] = "{\"ten\":10, \"array\":[11, 12, 13], \"fourteen\":14}";

			System.out.println ("Original strings:");
			for (int i = 0; i < num_string; ++i) {
				System.out.println ("string_data[" + i + "] = " + string_data[i]);
			}

			// Marshal the data

			MarshalImpJsonWriter writer = new MarshalImpJsonWriter();

			writer.marshalMapBegin (null);

			for (int i = 0; i < num_string; ++i) {
				writer.marshalString ("embed" + i, string_data[i]);
			}

			for (int i = 0; i < num_string; ++i) {
				writer.marshalJsonString ("merge" + i, string_data[i]);
			}

			writer.marshalMapEnd ();

			if (!( writer.check_write_complete() )) {
				System.out.println ("Writer reports writing not complete");
				return;
			}

			String json_string = writer.get_json_string();

			writer = null;

			// Unmarshal the data

			MarshalImpJsonReader reader = new MarshalImpJsonReader (json_string);

			reader.unmarshalMapBegin (null);

			System.out.println ("Embedded unmarshaled strings:");
			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalString ("embed" + i);
				System.out.println ("embed[" + i + "] = " + x);
			}

			System.out.println ("Merged unmarshaled strings:");
			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalJsonString ("merge" + i);
				System.out.println ("merge[" + i + "] = " + x);
			}

			reader.unmarshalMapEnd ();

			if (!( reader.check_read_complete() )) {
				System.out.println ("Reader reports reading not complete");
				return;
			}

			System.out.println ("JSON string:");
			System.out.println (json_string);

			return;
		}




		// Subcommand : Test #4
		// Command format:
		//  test4
		// Test marshaling and unmarshaling JSON strings.
		// Same as test #3 except the outer container is an array.

		if (args[0].equalsIgnoreCase ("test4")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("MarshalImpJsonWriter : Invalid 'test4' subcommand");
				return;
			}

			// Set up data

			int num_string = 4;
			String[] string_data = new String[num_string];

			string_data[0] = "{\"one\":1, \"two\":2, \"three\":3, \"four\":4}";
			string_data[1] = "[5, 6, 7, 8, 9]";
			string_data[2] = "";
			string_data[3] = "{\"ten\":10, \"array\":[11, 12, 13], \"fourteen\":14}";

			System.out.println ("Original strings:");
			for (int i = 0; i < num_string; ++i) {
				System.out.println ("string_data[" + i + "] = " + string_data[i]);
			}

			// Marshal the data

			MarshalImpJsonWriter writer = new MarshalImpJsonWriter();

			writer.marshalArrayBegin (null, num_string * 2);

			for (int i = 0; i < num_string; ++i) {
				writer.marshalString (null, string_data[i]);
			}

			for (int i = 0; i < num_string; ++i) {
				writer.marshalJsonString (null, string_data[i]);
			}

			writer.marshalArrayEnd ();

			if (!( writer.check_write_complete() )) {
				System.out.println ("Writer reports writing not complete");
				return;
			}

			String json_string = writer.get_json_string();

			writer = null;

			// Unmarshal the data

			MarshalImpJsonReader reader = new MarshalImpJsonReader (json_string);

			int array_size = reader.unmarshalArrayBegin (null);
			if (array_size != num_string * 2) {
				System.out.println ("Reader reports incorrect array size: " + array_size);
				return;
			}

			System.out.println ("Embedded unmarshaled strings:");
			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalString (null);
				System.out.println ("embed[" + i + "] = " + x);
			}

			System.out.println ("Merged unmarshaled strings:");
			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalJsonString (null);
				System.out.println ("merge[" + i + "] = " + x);
			}

			reader.unmarshalArrayEnd ();

			if (!( reader.check_read_complete() )) {
				System.out.println ("Reader reports reading not complete");
				return;
			}

			System.out.println ("JSON string:");
			System.out.println (json_string);

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("MarshalImpJsonWriter : Unrecognized subcommand : " + args[0]);
		return;

	}




}
