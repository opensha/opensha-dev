package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;


/**
 * Common base class for all timeline classes.
 * Author: Michael Barall 04/10/2018.
 *
 * The AAFS server maintains a timeline for each event it is tracking.
 */
public abstract class TimelineBase {

	//----- Constants and variables -----

	// Action code for this timeline entry.

	protected int actcode;

	public static final int ACTCODE_MIN        = 1;
	public static final int ACTCODE_TRACK      = 1;	// Initiate tracking the event
	public static final int ACTCODE_UNTRACK    = 2;	// Stop tracking the event
	public static final int ACTCODE_FORECAST   = 3;	// Generate a forecast
	public static final int ACTCODE_ANALYST    = 4;	// Supply analyst parameters
	public static final int ACTCODE_MAX        = 4;

	// Origin of this timeline.

	protected int origin;

	public static final int ORIGIN_MIN         = 1;
	public static final int ORIGIN_PDL         = 1;	// Injected from PDL notification
	public static final int ORIGIN_ANALYST     = 2;	// Analyst intervention
	public static final int ORIGIN_SYNC        = 3;	// Synchronizing to Comcat or backup
	public static final int ORIGIN_MAX         = 3;

	// Status of this timeline.

	protected int status;

	public static final int STATUS_MIN         = 1;
	public static final int STATUS_CREATED     = 1;	// Initial state, timeline created but no forecasts yet
	public static final int STATUS_RUNNING     = 2;	// Forecasts are being automatically generated
	public static final int STATUS_STOPPED     = 3;	// Timeline stopped, no further forecasts 
	public static final int STATUS_MAX         = 3;

	// Reason that the timeline has stopped.

	protected int stop_reason;

	public static final int SREASON_MIN        = 1;
	public static final int SREASON_NONE       = 1;	// Timeline not stopped yet
	public static final int SREASON_EXPIRED    = 2;	// Last forecast has been generated
	public static final int SREASON_ANALYST    = 3;	// Analyst has directed that event be untracked
	public static final int SREASON_FORESHOCK  = 4;	// The event is a foreshock
	public static final int SREASON_ERROR      = 5;	// Stopped because of unrecoverable errors
	public static final int SREASON_MAX        = 5;

	// Analyst that most recently reviewed this event, or "" if none.

	protected String analyst_id;

	// Analyst remark for this event, or "" if none.

	protected String analyst_remark;

	// Time at which analyst reviewed this event, in milliseconds since the epoch, or 0L if none.

	protected long analyst_time;

	// Parameters supplied by the analyst.
	// This can be null, if the analyst has not supplied any parameters.

	protected ForecastParameters analyst_params;

	// Parameters for the current forecast.
	// This cannot be null.
	// Contains the parameters used for the most recent forecast.
	// If no forecast has occurred yet, contains just mainshock parameters, with forecast_lag == 0L.

	protected ForecastParameters forecast_params;

	// Results for the current forecast.
	// Contains the results for the most recent forecast.
	// Can be null, if no forecast has occurred yet.

	protected ForecastResults forecast_results;

	// Time lag at which the next forecast will occur, in milliseconds since the mainshock.
	// The value is -1L if the timeline is stopped, indicating that no more forecasts will occur.

	protected long next_forecast_lag;

	// If this event is a foreshock, this is the event_id of the mainshock.
	// Otherwise, this is an empty string "".

	protected String foreshock_event_id;




	//----- Getters and setters -----

	public int get_actcode () {
		return actcode;
	}

	public int get_origin () {
		return origin;
	}

	public int get_status () {
		return status;
	}

	public int get_stop_reason () {
		return stop_reason;
	}

	public String get_analyst_id () {
		return analyst_id;
	}

	public String get_analyst_remark () {
		return analyst_remark;
	}

	public long get_analyst_time () {
		return analyst_time;
	}

	public ForecastParameters get_analyst_params () {
		return analyst_params;
	}

	public ForecastParameters get_forecast_params () {
		return forecast_params;
	}

	public ForecastResults get_forecast_results () {
		return forecast_results;
	}

	public long get_next_forecast_lag () {
		return next_forecast_lag;
	}

	public String get_foreshock_event_id () {
		return foreshock_event_id;
	}




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public TimelineBase () {}

	// Copy from another entry.

	public void copy_from (TimelineBase other) {
		actcode            = other.actcode;
		origin             = other.origin;
		status             = other.status;
		stop_reason        = other.stop_reason;
		analyst_id         = other.analyst_id;
		analyst_remark     = other.analyst_remark;
		analyst_time       = other.analyst_time;
		analyst_params     = other.analyst_params;
		forecast_params    = other.forecast_params;
		forecast_results   = other.forecast_results;
		next_forecast_lag  = other.next_forecast_lag;
		foreshock_event_id = other.foreshock_event_id;
		return;
	}


	///**
	// * Construct from given values.
	// */
	//public TimelineBase (int origin, String analyst_id, String analyst_remark) {
	//	if (!( origin >= ORIGIN_MIN && origin <= ORIGIN_MAX )) {
	//		throw new IllegalArgumentException ("TimelineBase: Origin out-of-range: origin = " + origin);
	//	}
	//	if (!( analyst_id != null )) {
	//		throw new NullPointerException ("TimelineBase: Parameter 'analyst_id' is null");
	//	}
	//	if (!( analyst_remark != null )) {
	//		throw new NullPointerException ("TimelineBase: Parameter 'analyst_remark' is null");
	//	}
	//
	//	this.origin = origin;
	//	this.analyst_id = analyst_id;
	//	this.analyst_remark = analyst_remark;
	//}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 15001;

	private static final String M_VERSION_NAME = "TimelineBase";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 15000;
	protected static final int MARSHAL_STATUS = 16001;
	protected static final int MARSHAL_UNTRACK = 17001;
	protected static final int MARSHAL_FORECAST = 18001;
	protected static final int MARSHAL_ANALYST = 19001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected abstract int get_marshal_type ();

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalInt                       ("actcode"           , actcode           );
		writer.marshalInt                       ("origin"            , origin            );
		writer.marshalInt                       ("status"            , status            );
		writer.marshalInt                       ("stop_reason"       , stop_reason       );
		writer.marshalString                    ("analyst_id"        , analyst_id        );
		writer.marshalString                    ("analyst_remark"    , analyst_remark    );
		writer.marshalLong                      ("analyst_time"      , analyst_time      );
		ForecastParameters.marshal_poly (writer, "analyst_params"    , analyst_params    );
		ForecastParameters.marshal_poly (writer, "forecast_params"   , forecast_params   );
		ForecastResults.marshal_poly    (writer, "forecast_results"  , forecast_results  );
		writer.marshalLong                      ("next_forecast_lag" , next_forecast_lag );
		writer.marshalString                    ("foreshock_event_id", foreshock_event_id);
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		actcode            = reader.unmarshalInt                       ("actcode"           );
		origin             = reader.unmarshalInt                       ("origin"            );
		status             = reader.unmarshalInt                       ("status"            );
		stop_reason        = reader.unmarshalInt                       ("stop_reason"       );
		analyst_id         = reader.unmarshalString                    ("analyst_id"        );
		analyst_remark     = reader.unmarshalString                    ("analyst_remark"    );
		analyst_time       = reader.unmarshalLong                      ("analyst_time"      );
		analyst_params     = ForecastParameters.unmarshal_poly (reader, "analyst_params"    );
		forecast_params    = ForecastParameters.unmarshal_poly (reader, "forecast_params"   );
		forecast_results   = ForecastResults.unmarshal_poly    (reader, "forecast_results"  );
		next_forecast_lag  = reader.unmarshalLong                      ("next_forecast_lag" );
		foreshock_event_id = reader.unmarshalString                    ("foreshock_event_id");

		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public TimelineBase unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, TimelineBase obj) {

		writer.marshalMapBegin (name);

		if (obj == null) {
			writer.marshalInt (M_TYPE_NAME, MARSHAL_NULL);
		} else {
			writer.marshalInt (M_TYPE_NAME, obj.get_marshal_type());
			obj.do_marshal (writer);
		}

		writer.marshalMapEnd ();

		return;
	}

	// Unmarshal object, polymorphic.

	public static TimelineBase unmarshal_poly (MarshalReader reader, String name) {
		TimelineBase result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("TimelineBase.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

//		case MARSHAL_STATUS:
//			result = new TimelineBaseTrack();
//			result.do_umarshal (reader);
//			break;
//
//		case MARSHAL_UNTRACK:
//			result = new TimelineBaseUntrack();
//			result.do_umarshal (reader);
//			break;
//
//		case MARSHAL_FORECAST:
//			result = new TimelineBaseForecast();
//			result.do_umarshal (reader);
//			break;
//
//		case MARSHAL_ANALYST:
//			result = new TimelineBaseAnalyst();
//			result.do_umarshal (reader);
//			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

}
