package scratch.aftershockStatistics.util;

/**
 * Interface for unmarshaling parameters/data from the OAF database.
 * Author: Michael Barall 03/31/2018.
 */
public interface MarshalReader {

	/**
	 * Begin a map context.
	 */
	public void unmarshalMapBegin (String name);

	/**
	 * End a map context.
	 */
	public void unmarshalMapEnd ();

	/**
	 * Begin an array context, return the array size.
	 */
	public int unmarshalArrayBegin (String name);

	/**
	 * End an array context.
	 */
	public void unmarshalArrayEnd ();

	/**
	 * Unmarshal a long.
	 */
	public long unmarshalLong (String name);

	/**
	 * Unmarshal a double.
	 */
	public double unmarshalDouble (String name);

	/**
	 * Unmarshal a string.  (Null strings are not allowed.)
	 */
	public String unmarshalString (String name);

	/**
	 * Unmarshal a long, with required minimum value.
	 */
	public default long unmarshalLong (String name, long minValue) {
		long x = unmarshalLong (name);
		if (x < minValue) {
			throw new MarshalException ("Unmarshaled long out-of-range: value = " + x + ", min = " + minValue);
		}
		return x;
	}

	/**
	 * Unmarshal a long, with required minimum and maximum values.
	 */
	public default long unmarshalLong (String name, long minValue, long maxValue) {
		long x = unmarshalLong (name);
		if (x < minValue || x > maxValue) {
			throw new MarshalException ("Unmarshaled long out-of-range: value = " + x + ", min = " + minValue + ", max = " + maxValue);
		}
		return x;
	}

	/**
	 * Unmarshal an int.
	 */
	public default int unmarshalInt (String name) {
		long x = unmarshalLong (name);
		if (x < (long)Integer.MIN_VALUE || x > (long)Integer.MAX_VALUE) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + Integer.MIN_VALUE + ", max = " + Integer.MAX_VALUE);
		}
		return (int)x;
	}

	/**
	 * Unmarshal an int, with required minimum value.
	 */
	public default int unmarshalInt (String name, int minValue) {
		long x = unmarshalLong (name);
		if (x < (long)minValue || x > (long)Integer.MAX_VALUE) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + minValue + ", max = " + Integer.MAX_VALUE);
		}
		return (int)x;
	}

	/**
	 * Unmarshal an int, with required minimum and maximum values.
	 */
	public default int unmarshalInt (String name, int minValue, int maxValue) {
		long x = unmarshalLong (name);
		if (x < (long)minValue || x > (long)maxValue) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + minValue + ", max = " + maxValue);
		}
		return (int)x;
	}

	/**
	 * Unmarshal a boolean.
	 */
	public default boolean unmarshalBoolean (String name) {
		long x = unmarshalLong (name);
		if (x == 0L) {
			return false;
		}
		if (x == 1L) {
			return true;
		}
		throw new MarshalException ("Unmarshaled boolean out-of-range: value = " + x + ", min = 0, max = 1");
	}

	/**
	 * Unmarshal a JSON string.  (Null strings are not allowed.)
	 * The string must contain a JSON object or array, or be an empty string.
	 * For JSON storage, the string is merged into the JSON instead of being
	 * embedded as string-valued data.  (An empty string becomes a JSON null.)
	 * The unmarshaled string may differ from the marshaled string due to JSON parsing.
	 * (Named element ordering, numeric formats, and spacing may be changed).
	 */
	public default String unmarshalJsonString (String name) {
		return unmarshalString (name);
	}

	/**
	 * Unmarshal a long array.
	 */
	public default long[] unmarshalLongArray (String name) {
		int n = unmarshalArrayBegin (name);
		long[] x = new long[n];
		for (int i = 0; i < n; ++i) {
			x[i] = unmarshalLong (null);
		}
		unmarshalArrayEnd ();
		return x;
	}

	/**
	 * Unmarshal a double array.
	 */
	public default double[] unmarshalDoubleArray (String name) {
		int n = unmarshalArrayBegin (name);
		double[] x = new double[n];
		for (int i = 0; i < n; ++i) {
			x[i] = unmarshalDouble (null);
		}
		unmarshalArrayEnd ();
		return x;
	}

	/**
	 * Unmarshal a string array.  (Null strings are not allowed.)
	 */
	public default String[] unmarshalStringArray (String name) {
		int n = unmarshalArrayBegin (name);
		String[] x = new String[n];
		for (int i = 0; i < n; ++i) {
			x[i] = unmarshalString (null);
		}
		unmarshalArrayEnd ();
		return x;
	}

	/**
	 * Unmarshal an int array.
	 */
	public default int[] unmarshalIntArray (String name) {
		int n = unmarshalArrayBegin (name);
		int[] x = new int[n];
		for (int i = 0; i < n; ++i) {
			x[i] = unmarshalInt (null);
		}
		unmarshalArrayEnd ();
		return x;
	}

}
