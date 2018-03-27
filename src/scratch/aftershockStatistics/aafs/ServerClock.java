package scratch.aftershockStatistics.aafs;

/**
 * Clock for AAFS server.
 * Author: Michael Barall 03/18/2018.
 *
 * This class provides functions for getting and setting the time.
 * Normally, time is obtained from System.currentTimeMillis().
 * But, with this class it is possible to manipulate the clock for accelerated testing.
 */
public class ServerClock {

	// The frozen clock time, or 0 for normal clock operation, in milliseconds since the epoch.

	private static volatile long frozen_time = 0L;


	/**
	 * get_time - Get the current time, in milliseconds since the epoch..
	 */
	public static long get_time () {
		long the_time = frozen_time;
		if (the_time == 0L) {
			the_time = System.currentTimeMillis();
		}
		return the_time;
	}


	/**
	 * freeze_time - Freeze the clock at the given time.
	 * @param at_time = Time at which to freeze the clock, in milliseconds since the epoch.
	 * A parameter value of 0L will un-freeze the clock.
	 */
	public static void freeze_time (long at_time) {
		frozen_time = at_time;
		return;
	}


	/**
	 * advance_frozen_time - Freeze the clock at the given time,
	 *                       if the given time is later than the current frozen time.
	 * @param at_time = Time at which to freeze the clock, in milliseconds since the epoch.
	 * A parameter value of 0L will un-freeze the clock.
	 */
	public static void advance_frozen_time (long at_time) {
		if (at_time == 0L || frozen_time == 0L || at_time > frozen_time) {
			frozen_time = at_time;
		}
		return;
	}


	/**
	 * get_true_time - Get the true current time, in milliseconds since the epoch..
	 */
	public static long get_true_time () {
		return System.currentTimeMillis();
	}
}
