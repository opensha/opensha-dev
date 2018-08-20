package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Operation payload for analyst intervention in an event.
 * Author: Michael Barall 05/22/2018.
 */
public class OpAnalystIntervene extends DBPayload {

	//----- Constants and variables -----

	// State change request values.

	public static final int ASREQ_MIN = 1;
	public static final int ASREQ_NONE = 1;			// Do not change state
	public static final int ASREQ_START = 2;		// Start or continue generating forecasts
	public static final int ASREQ_STOP = 3;			// Stop generating forecasts
	public static final int ASREQ_WITHDRAW = 4;		// Withdraw the timeline
	public static final int ASREQ_MAX = 4;

	// Requested state change.

	public int state_change;

	// Flag, true to create timeline if it doesn't exist.

	public boolean f_create_timeline;

	// Parameters supplied by the analyst, or null if none.

	public AnalystOptions analyst_options;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpAnalystIntervene () {}


	// Set up the contents, for no analyst data.

	public void setup (int the_state_change, boolean the_f_create_timeline) {
		state_change = the_state_change;
		f_create_timeline = the_f_create_timeline;
		analyst_options = null;
		return;
	}


	// Set up the contents, with analyst data

	public void setup (int the_state_change, boolean the_f_create_timeline,
						AnalystOptions the_analyst_options) {
		state_change = the_state_change;
		f_create_timeline = the_f_create_timeline;
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

	private static final int MARSHAL_VER_1 = 32001;

	private static final String M_VERSION_NAME = "OpAnalystIntervene";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalInt                       ("state_change"       , state_change       );
		writer.marshalBoolean                   ("f_create_timeline"  , f_create_timeline  );

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

		state_change        = reader.unmarshalInt                       ("state_change"       , ASREQ_MIN, ASREQ_MAX);
		f_create_timeline   = reader.unmarshalBoolean                   ("f_create_timeline"  );

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
	public OpAnalystIntervene unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a pending task.

	@Override
	public OpAnalystIntervene unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

}
