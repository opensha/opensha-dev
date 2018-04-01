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

}
