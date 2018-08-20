package scratch.aftershockStatistics;

/**
 * Exception class for Comcat access soft errors.
 * Author: Michael Barall 08/12/2018.
 *
 * This exception indicates that inconsistent data was received during an
 * operation that involved multiple Comcat queries.  For example, two queries
 * that disagree about whether a given event exists, or that disagree about
 * which IDs refer to the same event.  This may indicate that the Comcat
 * database was changing during the operation, and so this may be a temporary
 * condition.  The appropriate action is to wait a while and retry.
 *
 * This does not indicate a problem with accessing the Comcat service.
 */
public class ComcatConflictException extends ComcatException {

	// Constructors.

	public ComcatConflictException () {
		super ();
	}

	public ComcatConflictException (String s) {
		super (s);
	}

	public ComcatConflictException (String message, Throwable cause) {
		super (message, cause);
	}

	public ComcatConflictException (Throwable cause) {
		super (cause);
	}

}
