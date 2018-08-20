package scratch.aftershockStatistics.aafs;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Date;

import java.text.SimpleDateFormat;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;


/**
 * Operation payload for intake of an event discovered by polling.
 * Author: Michael Barall 05/20/2018.
 */
public class OpIntakePoll extends DBPayload {

	//----- Constants and variables -----

	// Event ID for delayed form of command, which can be a Comcat ID or a timeline ID.
	// Note: The delayed form of the command has the task's event_id equal to EVID_POLL.
	// It is executed by staging with event_id changed to delayed_event_id.  This makes
	// it possible to efficiently find and delete all delayed commands.  In non-delayed
	// commands, delayed_event_id should equal the task's event_id.

	public String delayed_event_id;

	// Parameters supplied by the analyst, or null if none.

	public AnalystOptions analyst_options;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpIntakePoll () {}


	// Set up the contents, for no analyst data.

	public void setup (String the_delayed_event_id) {
		delayed_event_id = the_delayed_event_id;
		analyst_options = null;
		return;
	}


	// Set up the contents, with analyst data

	public void setup (String the_delayed_event_id, AnalystOptions the_analyst_options) {
		delayed_event_id = the_delayed_event_id;
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

	private static final int MARSHAL_VER_1 = 41001;

	private static final String M_VERSION_NAME = "OpIntakePoll";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalString                    ("delayed_event_id"   , delayed_event_id   );
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

		delayed_event_id    = reader.unmarshalString                    ("delayed_event_id"   );
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
	public OpIntakePoll unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a pending task.

	@Override
	public OpIntakePoll unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

}
