package scratch.aftershockStatistics;

import java.util.Arrays;

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

	//----- Implementation of MarshalReader -----

	/**
	 * Unmarshal a long.
	 */
	@Override
	public long unmarshalLong () {
		if (long_read_index == long_write_index) {
			throw new MarshalException ("Unmarshal long end-of-data: size = " + long_write_index);
		}
		return long_store[long_read_index++];
	}

	/**
	 * Unmarshal a double.
	 */
	@Override
	public double unmarshalDouble () {
		if (double_read_index == double_write_index) {
			throw new MarshalException ("Unmarshal double end-of-data: size = " + double_write_index);
		}
		return double_store[double_read_index++];
	}

	/**
	 * Unmarshal a string.  (Null strings are not allowed.)
	 */
	@Override
	public String unmarshalString () {
		if (string_read_index == string_write_index) {
			throw new MarshalException ("Unmarshal string end-of-data: size = " + string_write_index);
		}
		return string_store[string_read_index++];
	}

	//----- Implementation of MarshalWriter -----

	/**
	 * Marshal a long.
	 */
	@Override
	public void marshalLong (long x) {
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
	public void marshalDouble (double x) {
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
	public void marshalString (String x) {
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
	}

	//----- Control -----

	/**
	 * Set indexes to read from beginning of array.
	 */
	public void begin_read () {
		long_read_index = 0;
		double_read_index = 0;
		string_read_index = 0;
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
		return;
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

			for (int i = 0; i < num_long; ++i) {
				writer.marshalLong (long_data[i]);
			}

			for (int i = 0; i < num_double; ++i) {
				writer.marshalDouble (double_data[i]);
			}

			for (int i = 0; i < num_string; ++i) {
				writer.marshalString (string_data[i]);
			}

			long[] m_long_store = writer.get_long_store();
			double[] m_double_store = writer.get_double_store();
			String[] m_string_store = writer.get_string_store();

			writer = null;

			// Unmarshal and check the data

			System.out.println ("Unmarshaling data ...");

			MarshalImpArray reader = new MarshalImpArray (m_long_store, m_double_store, m_string_store);

			int errors = 0;

			for (int i = 0; i < num_long; ++i) {
				long x = reader.unmarshalLong();
				if (x != long_data[i]) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched long: i = " + i + ", d = " + long_data[i] + ", x = " + x);
					}
				}
			}

			for (int i = 0; i < num_double; ++i) {
				double x = reader.unmarshalDouble();
				if (x != double_data[i]) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched double: i = " + i + ", d = " + double_data[i] + ", x = " + x);
					}
				}
			}

			for (int i = 0; i < num_string; ++i) {
				String x = reader.unmarshalString();
				if (!( x.equals(string_data[i]) )) {
					++errors;
					if (errors <= 10) {
						System.out.println ("Mismatched string: i = " + i + ", d = " + string_data[i] + ", x = " + x);
					}
				}
			}

			System.out.println ("Error count: " + errors);

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("MarshalImpArray : Unrecognized subcommand : " + args[0]);
		return;

	}




}
