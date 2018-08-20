package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Operation payload for intake of an event as a sync event.
 * Author: Michael Barall 05/20/2018.
 */
public class OpIntakeSync extends DBPayload {

	//----- Constants and variables -----

	// Parameters supplied by the analyst, or null if none.

	public AnalystOptions analyst_options;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpIntakeSync () {}


	// Set up the contents, for no analyst data.

	public void setup () {
		analyst_options = null;
		return;
	}


	// Set up the contents, with analyst data

	public void setup (AnalystOptions the_analyst_options) {
		analyst_options = the_analyst_options;
		return;
	}


	// Return the effective analyst parameters, or null if none.

	public ForecastParameters get_eff_analyst_params () {
		if (analyst_options != null) {
			return analyst_options.analyst_params;
		}
		return null;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 31001;

	private static final String M_VERSION_NAME = "OpIntakeSync";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		AnalystOptions.marshal_poly     (writer, "analyst_options"    , analyst_options    );

		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Superclass

		super.do_umarshal (reader);

		// Contents

		analyst_options     = AnalystOptions.unmarshal_poly     (reader, "analyst_options"    );

		return;
	}

	// Marshal object.

	@Override
	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	@Override
	public OpIntakeSync unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a pending task.

	@Override
	public OpIntakeSync unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

}
