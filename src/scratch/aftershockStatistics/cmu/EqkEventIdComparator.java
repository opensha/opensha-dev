package scratch.aftershockStatistics.cmu;

import java.util.Comparator;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

/**
 * <p>Title: ObsEqkRupEventOriginTimeComparator</p>
 *
 * <p>Description: This compares 2 observed ObsEqkRupture objects
 * based on their eventIds.
 * </p>
 *
 *
 * @author Nitin Gupta, rewritten by Ned Field
 * @version 1.0
 */
public class EqkEventIdComparator
        implements Comparator<ObsEqkRupture>, java.io.Serializable {

    /**
     * Compares the origin times of the two arguments. Returns negative one, zero, or
     * positive one depending on whether the first origin time is less than,
     * equal to, or greater than the second, respectively.
     */
    public int compare(ObsEqkRupture rupEvent1, ObsEqkRupture rupEvent2) {
        return new Integer(rupEvent1.getEventId().compareTo(rupEvent2.getEventId()));
    }

}