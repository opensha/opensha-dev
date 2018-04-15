package scratch.aftershockStatistics.util;

/**
 * Interface for marshaling parameters/data to the OAF database.
 * Author: Michael Barall 03/31/2018.
 */
public interface MarshalWriter {

	/**
	 * Begin a map context.
	 */
	public void marshalMapBegin (String name);

	/**
	 * End a map context.
	 */
	public void marshalMapEnd ();

	/**
	 * Begin an array context, specify the array size.
	 */
	public void marshalArrayBegin (String name, int array_size);

	/**
	 * End an array context.
	 */
	public void marshalArrayEnd ();

	/**
	 * Marshal a long.
	 */
	public void marshalLong (String name, long x);

	/**
	 * Marshal a double.
	 */
	public void marshalDouble (String name, double x);

	/**
	 * Marshal a string.  (Null strings are not allowed.)
	 */
	public void marshalString (String name, String x);

	/**
	 * Marshal an int.
	 */
	public default void marshalInt (String name, int x) {
		marshalLong (name, (long)x);
		return;
	}

	/**
	 * Marshal a boolean.
	 */
	public default void marshalBoolean (String name, boolean x) {
		marshalLong (name, x ? 1L : 0L);
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
	public default void marshalJsonString (String name, String x) {
		marshalString (name, x);
		return;
	}

	/**
	 * Marshal a long array.
	 */
	public default void marshalLongArray (String name, long[] x) {
		int n = x.length;
		marshalArrayBegin (name, n);
		for (int i = 0; i < n; ++i) {
			 marshalLong (null, x[i]);
		}
		marshalArrayEnd ();
		return;
	}

	/**
	 * Marshal a double array.
	 */
	public default void marshalDoubleArray (String name, double[] x) {
		int n = x.length;
		marshalArrayBegin (name, n);
		for (int i = 0; i < n; ++i) {
			 marshalDouble (null, x[i]);
		}
		marshalArrayEnd ();
		return;
	}

	/**
	 * Marshal a string array.  (Null strings are not allowed.)
	 */
	public default void marshalStringArray (String name, String[] x) {
		int n = x.length;
		marshalArrayBegin (name, n);
		for (int i = 0; i < n; ++i) {
			 marshalString (null, x[i]);
		}
		marshalArrayEnd ();
		return;
	}

	/**
	 * Marshal an int array.
	 */
	public default void marshalIntArray (String name, int[] x) {
		int n = x.length;
		marshalArrayBegin (name, n);
		for (int i = 0; i < n; ++i) {
			 marshalInt (null, x[i]);
		}
		marshalArrayEnd ();
		return;
	}

}
