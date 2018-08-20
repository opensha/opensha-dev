package scratch.aftershockStatistics.util;

import java.util.Comparator;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

/**
 * Comparator to sort ruptures in increasing order by time (earliest rupture first).
 * Author: Michael Barall 07/29/2018.
 *
 * Events of equal time are sorted in decreasing order by magnitude.
 *
 * Events of equal time and magnitude are sorted in order by event ID.
 * Null event IDs are sorted before non-null event IDs.
 */
public class ObsEqkRupMinTimeComparator implements Comparator<ObsEqkRupture> {
	
	// Compares its two arguments for order. Returns a negative integer, zero, or a positive
	// integer as the first argument is less than, equal to, or greater than the second.

	@Override
    public int compare (ObsEqkRupture rupEvent1, ObsEqkRupture rupEvent2) {

		// Order by time, earliest first

		int result = Long.compare (rupEvent1.getOriginTime(), rupEvent2.getOriginTime());

		if (result == 0) {

			// Order by magnitude, largest first

			result = Double.compare (rupEvent2.getMag(), rupEvent1.getMag());

			if (result == 0) {

				// Order by event ID, lexicographically

				String eid1 = rupEvent1.getEventId();
				String eid2 = rupEvent2.getEventId();
				result = ( (eid1 == null)
							? ((eid2 == null) ? 0 : -1)
							: ((eid2 == null) ? 1 : (eid1.compareTo(eid2))) );
			}
		}

		return result;
    }

}
