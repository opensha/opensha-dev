package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Operation payload for generating a forecast.
 * Author: Michael Barall 05/07/2018.
 */
public class OpGenerateForecast extends DBPayload {

	//----- Constants and variables -----

	// Time stamp for the timeline entry when this was issued, in milliseconds since the epoch.

	public long action_time;

	// Time lag at which the last forecast occured, in milliseconds since the mainshock.
	// The value is -1L if there have been no prior forecasts.

	public long last_forecast_lag;

	// Time lag at which the next forecast will occur, in milliseconds since the mainshock.

	public long next_forecast_lag;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public OpGenerateForecast () {}


	// Set up the contents.

	public void setup (long the_action_time, long the_last_forecast_lag, long the_next_forecast_lag) {
		action_time = the_action_time;
		last_forecast_lag = the_last_forecast_lag;
		next_forecast_lag = the_next_forecast_lag;
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 28001;

	private static final String M_VERSION_NAME = "OpGenerateForecast";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalLong ("action_time"      , action_time      );
		writer.marshalLong ("last_forecast_lag", last_forecast_lag);
		writer.marshalLong ("next_forecast_lag", next_forecast_lag);

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

		action_time       = reader.unmarshalLong ("action_time"      );
		last_forecast_lag = reader.unmarshalLong ("last_forecast_lag");
		next_forecast_lag = reader.unmarshalLong ("next_forecast_lag");

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
	public OpGenerateForecast unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a pending task.

	@Override
	public OpGenerateForecast unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

}
