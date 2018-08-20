package scratch.aftershockStatistics.util;

import java.util.Comparator;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

/**
 * Comparator to sort ruptures in increasing order by event ID.
 * Author: Michael Barall 07/29/2018.
 *
 * Null event IDs are sorted before non-null event IDs.
 *
 * Authors of prior versions: Nitin Gupta, Ned Field
 */
public class ObsEqkRupEventIdComparator implements Comparator<ObsEqkRupture> {
	
	// Compares its two arguments for order. Returns a negative integer, zero, or a positive
	// integer as the first argument is less than, equal to, or greater than the second.

	@Override
    public int compare (ObsEqkRupture rupEvent1, ObsEqkRupture rupEvent2) {

		// Order by event ID, lexicographically

		String eid1 = rupEvent1.getEventId();
		String eid2 = rupEvent2.getEventId();
		return (eid1 == null)
				? ((eid2 == null) ? 0 : -1)
				: ((eid2 == null) ? 1 : (eid1.compareTo(eid2)));
    }

}
