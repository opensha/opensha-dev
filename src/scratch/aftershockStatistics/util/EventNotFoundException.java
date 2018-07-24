package scratch.aftershockStatistics.util;

/**
 * Exception class that can be used to signal that an event was not found.
 * Author: Michael Barall 06/09/2018.
 */
public class EventNotFoundException extends RuntimeException {

	// Constructors.

	public EventNotFoundException () {
		super ();
	}

	public EventNotFoundException (String s) {
		super (s);
	}

	public EventNotFoundException (String message, Throwable cause) {
		super (message, cause);
	}

	public EventNotFoundException (Throwable cause) {
		super (cause);
	}

}
