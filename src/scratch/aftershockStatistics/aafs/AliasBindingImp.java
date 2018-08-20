package scratch.aftershockStatistics.aafs;


/**
 * Simple, immutable implementation of the AliasBinding interface.
 * Author: Michael Barall 06/28/2018.
 */
public class AliasBindingImp implements AliasBinding {


	//----- Values -----

	// The Comcat ID.

	private String comcat_id;

	// The timeline ID.

	private String timeline_id;

	// The relation.

	private int relation;


	//----- Construction -----

	public AliasBindingImp (String comcat_id, String timeline_id, int relation) {
		this.comcat_id = comcat_id;
		this.timeline_id = timeline_id;
		this.relation = relation;
	}


	//----- Implementation of AliasBinding -----

	// Get the Comcat ID.

	@Override
	public String get_comcat_id () {
		return comcat_id;
	}

	// Get the timeline ID.

	@Override
	public String get_timeline_id () {
		return timeline_id;
	}

	// Get the relation, returning one of the enumarated values.

	@Override
	public int get_relation () {
		return relation;
	}

}
