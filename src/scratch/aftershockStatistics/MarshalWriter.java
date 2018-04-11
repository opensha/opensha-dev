package scratch.aftershockStatistics;

/**
 * Interface for marshaling parameters/data to the OAF database.
 * Author: Michael Barall 03/31/2018.
 */
public interface MarshalWriter {

	/**
	 * Marshal a long.
	 */
	public void marshalLong (long x);

	/**
	 * Marshal a double.
	 */
	public void marshalDouble (double x);

	/**
	 * Marshal a string.  (Null strings are not allowed.)
	 */
	public void marshalString (String x);

	/**
	 * Marshal an int.
	 */
	public default void marshalInt (int x) {
		marshalLong ((long)x);
		return;
	}

	/**
	 * Marshal a boolean.
	 */
	public default void marshalBoolean (boolean x) {
		marshalLong (x ? 1L : 0L);
		return;
	}

}
