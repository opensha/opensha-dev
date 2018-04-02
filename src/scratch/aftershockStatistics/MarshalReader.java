package scratch.aftershockStatistics;

/**
 * Interface for unmarshaling parameters/data from the OAF database.
 * Author: Michael Barall 03/31/2018.
 */
public interface MarshalReader {

	/**
	 * Unmarshal a long.
	 */
	public long unmarshalLong ();

	/**
	 * Unmarshal a double.
	 */
	public double unmarshalDouble ();

	/**
	 * Unmarshal a long, with required minimum value.
	 */
	public default long unmarshalLong (long minValue) {
		long x = unmarshalLong ();
		if (x < minValue) {
			throw new MarshalException ("Unmarshaled long out-of-range: value = " + x + ", min = " + minValue);
		}
		return x;
	}

	/**
	 * Unmarshal a long, with required minimum and maximum values.
	 */
	public default long unmarshalLong (long minValue, long maxValue) {
		long x = unmarshalLong ();
		if (x < minValue || x > maxValue) {
			throw new MarshalException ("Unmarshaled long out-of-range: value = " + x + ", min = " + minValue + ", max = " + maxValue);
		}
		return x;
	}

	/**
	 * Unmarshal an int.
	 */
	public default int unmarshalInt () {
		long x = unmarshalLong ();
		if (x < (long)Integer.MIN_VALUE || x > (long)Integer.MAX_VALUE) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + Integer.MIN_VALUE + ", max = " + Integer.MAX_VALUE);
		}
		return (int)x;
	}

	/**
	 * Unmarshal an int, with required minimum value.
	 */
	public default int unmarshalInt (int minValue) {
		long x = unmarshalLong ();
		if (x < (long)minValue || x > (long)Integer.MAX_VALUE) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + minValue + ", max = " + Integer.MAX_VALUE);
		}
		return (int)x;
	}

	/**
	 * Unmarshal an int, with required minimum and maximum values.
	 */
	public default int unmarshalInt (int minValue, int maxValue) {
		long x = unmarshalLong ();
		if (x < (long)minValue || x > (long)maxValue) {
			throw new MarshalException ("Unmarshaled int out-of-range: value = " + x + ", min = " + minValue + ", max = " + maxValue);
		}
		return (int)x;
	}

	/**
	 * Unmarshal a boolean.
	 */
	public default boolean unmarshalBoolean () {
		long x = unmarshalLong ();
		if (x == 0L) {
			return false;
		}
		if (x == 1L) {
			return true;
		}
		throw new MarshalException ("Unmarshaled boolean out-of-range: value = " + x + ", min = 0, max = 1");
	}

}
