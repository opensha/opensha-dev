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

}
