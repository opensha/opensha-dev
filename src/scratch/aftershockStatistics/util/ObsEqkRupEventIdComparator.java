package scratch.aftershockStatistics.util;

import java.util.Comparator;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

/**
 * <p>Title: ObsEqkRupEventIdComparator</p>
 *
 * <p>Description: This compares 2 observed ObsEqkRupture objects
 * based on their eventIds.
 * </p>
 *
 *
 * @author Nitin Gupta, rewritten by Ned Field
 * @version 1.0
 */
public class ObsEqkRupEventIdComparator
        implements Comparator<ObsEqkRupture>, java.io.Serializable {

    /**
     * Compares the event ids of the two arguments. Returns a negative integer, zero, or
     * a positive integer depending on whether the first event id is less than,
     * equal to, or greater than the second, respectively.
     */
    public int compare(ObsEqkRupture rupEvent1, ObsEqkRupture rupEvent2) {
        return rupEvent1.getEventId().compareTo(rupEvent2.getEventId());
    }

}