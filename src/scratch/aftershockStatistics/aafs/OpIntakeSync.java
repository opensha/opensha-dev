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

	// Flag, true if this contains analyst data

	public boolean f_has_analyst;

	// Analyst that most recently reviewed this event, or "" if none.

	public String analyst_id;

	// Analyst remark for this event, or "" if none.

	public String analyst_remark;

	// Time at which analyst reviewed this event, in milliseconds since the epoch, or 0L if none.

	public long analyst_time;

	// Parameters supplied by the analyst.
	// This can be null, if the analyst has not supplied any parameters,
	// or if the analyst has intervened a second time to "unsupply" parameters.

	public ForecastParameters analyst_params;

	// Time lag at which an extra forecast is requested, in milliseconds since the mainshock.
	// The value is -1L if there has been no extra forecast requested.

	public long extra_forecast_lag;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpIntakeSync () {}


	// Set up the contents, for no analyst data.

	public void setup () {
		f_has_analyst = false;
		analyst_id = "";
		analyst_remark = "";
		analyst_time = 0L;
		analyst_params = null;
		extra_forecast_lag = -1L;
		return;
	}


	// Set up the contents, with analyst data

	public void setup (String the_analyst_id, String the_analyst_remark, long the_analyst_time,
						ForecastParameters the_analyst_params, long the_extra_forecast_lag) {
		f_has_analyst = true;
		analyst_id = the_analyst_id;
		analyst_remark = the_analyst_remark;
		analyst_time = the_analyst_time;
		analyst_params = the_analyst_params;
		extra_forecast_lag = the_extra_forecast_lag;
		return;
	}


	// Return the effective analyst parameters, or null if none.

	public ForecastParameters get_eff_analyst_params () {
		if (f_has_analyst) {
			return analyst_params;
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

		writer.marshalBoolean                   ("f_has_analyst"      , f_has_analyst      );
		writer.marshalString                    ("analyst_id"         , analyst_id         );
		writer.marshalString                    ("analyst_remark"     , analyst_remark     );
		writer.marshalLong                      ("analyst_time"       , analyst_time       );
		ForecastParameters.marshal_poly (writer, "analyst_params"     , analyst_params     );
		writer.marshalLong                      ("extra_forecast_lag" , extra_forecast_lag );

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

		f_has_analyst       = reader.unmarshalBoolean                   ("f_has_analyst"      );
		analyst_id          = reader.unmarshalString                    ("analyst_id"         );
		analyst_remark      = reader.unmarshalString                    ("analyst_remark"     );
		analyst_time        = reader.unmarshalLong                      ("analyst_time"       );
		analyst_params      = ForecastParameters.unmarshal_poly (reader, "analyst_params"     );
		extra_forecast_lag  = reader.unmarshalLong                      ("extra_forecast_lag" );

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
