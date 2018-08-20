package scratch.aftershockStatistics;

/**
 * Exception class for Comcat access that fails to find expected data.
 * Author: Michael Barall 08/12/2018.
 *
 * This exception indicates that data which was previously found in Comcat
 * (or PDL) cannot be found now.  This may be a temporary condition, and so
 * the appropriate action is to wait a while and retry.
 *
 * This does not indicate a problem with accessing the Comcat service.  It also
 * does not indicate a deleted event (Comcat HTTP status 409), although event
 * deletion may lead to a failure to find expected data.  Note that deleted
 * events sometimes (usually?) re-appear within a few hours.
 */
public class ComcatRemovedException extends ComcatException {

	// Constructors.

	public ComcatRemovedException () {
		super ();
	}

	public ComcatRemovedException (String s) {
		super (s);
	}

	public ComcatRemovedException (String message, Throwable cause) {
		super (message, cause);
	}

	public ComcatRemovedException (Throwable cause) {
		super (cause);
	}

}
