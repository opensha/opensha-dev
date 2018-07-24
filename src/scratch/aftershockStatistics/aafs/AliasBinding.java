package scratch.aftershockStatistics.aafs;


/**
 * Interface that represents a relation between a Comcat ID and a timeline ID.
 * Author: Michael Barall 06/28/2018.
 */
public interface AliasBinding {

	// Get the Comcat ID.

	public String get_comcat_id ();

	// Get the timeline ID.

	public String get_timeline_id ();

	// Get the relation, returning one of the enumarated values.

	public static final int BINDREL_PRIMARY = 1;		// Comcat ID is currently the primary ID for the timeline.
	public static final int BINDREL_SECONDARY = 2;		// Comcat ID is currently a secondary ID for the timeline.
	public static final int BINDREL_REMOVED = 3;		// Comcat ID was previously associated with the timeline but has been removed.

	public int get_relation ();

}
