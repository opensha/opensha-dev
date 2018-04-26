package scratch.aftershockStatistics.util;

/**
 * Exception class for parameter/data marshaling/unmarshaling.
 * Author: Michael Barall 03/31/2018.
 */
public class MarshalException extends RuntimeException {

	// Constructors.

	public MarshalException () {
		super ();
	}

	public MarshalException (String s) {
		super (s);
	}

	public MarshalException (String message, Throwable cause) {
		super (message, cause);
	}

	public MarshalException (Throwable cause) {
		super (cause);
	}

}
