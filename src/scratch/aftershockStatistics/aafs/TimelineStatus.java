package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;

import scratch.aftershockStatistics.CompactEqkRupList;


/**
 * Status entry in the event timeline.
 * Author: Michael Barall 04/10/2018.
 *
 * The AAFS server maintains a timeline for each event it is tracking.
 *
 * The actcode defines the most recent action in the timeline.  It is copied
 * into TimelineEntry.  In principle, actcode could select amont different
 * classes, but currently TimelineStatus is used for all timeline entries.
 *
 * The action_time that appears in TimelineEntry is in all cases equal to the
 * mainshock_time plus forecast_lag (in forecast_parameters).  It is necessary
 * that the action_time be strictly increasing so that timeline entries are
 * properly ordered in the database.  The first timeline entry has forecast_lag
 * equal to zero.  The action_time should always be less than the current time.
 *
 * The event_id that appears in TimelineEntry is copied from forecast_parameters.
 *
 * ACTCODE_TRACK:  The system has just begun to track the event.  This may be
 * because the event has just been introduced (through PDL, analyst, or sync),
 * or because an analyst has re-started tracking for an event previously untracked.
 * The forecast_parameters must contain event_id, forecast_lag, and mainshock
 * parameters; they can optionally contain analyst-supplied parameters.
 * There are no forecast_results.
 * 
 */
public class TimelineStatus extends DBPayload {

	//----- Envelope information -----
	//
	// These variables are copied to/from TimelineEntry rather than being marshaled,
	// so they are exposed to database code.

	// Event ID for this timeline.

	public String event_id;

	// Action code for this timeline entry.

	public int actcode;

	public static final int ACTCODE_MIN                 = 1;
	public static final int ACTCODE_TRACK               = 1;
		// This is the first entry in the timeline, when event tracking begins.
		// The cause is recorded in fc_origin below.  It may also contain analsyt
		// parameter values.  Only the first timeline entry can be ACTCODE_TRACK.
	public static final int ACTCODE_FORECAST            = 2;
		// The entry contains a forecast.
	public static final int ACTCODE_ANALYST             = 3;
		// The entry contains an analyst intervention.  The analyst may assign
		// parameter values and turn tracking on or off.
	public static final int ACTCODE_PDL_UPDATE          = 4;
		// The entry contains an update to PDL status, usually the final result of a
		// delayed PDL report.  This is recorded in pdl_status below.
	public static final int ACTCODE_ERROR               = 5;
		// The entry is the last timeline entry due to an unrecoverable error.
		// In this case, the payload may not parse correctly.
	public static final int ACTCODE_WITHDRAWN           = 6;
		// The entry is the last timeline entry due to the automatic system
		// determining that the event is not eligible to be tracked.
	public static final int ACTCODE_COMCAT_FAIL         = 7;
		// The entry is the last timeline entry due to inability to retrieve
		// necessary data from ComCat.
	public static final int ACTCODE_FORESHOCK           = 8;
		// The entry is the last timeline entry because the automatic system
		// has determined it is a foreshock.
	public static final int ACTCODE_EXPIRED             = 9;
		// The entry is the last timeline entry because the automatic system
		// has determined it is expired. Note that an expired action may not
		// appear if expired status was set on the final forecast action.
	public static final int ACTCODE_STATUS_UPDATE       = 10;
		// The entry contains an update to forecast status.  This is recorded
		// in fc_status below.  There may also be a change to pdl_status.
	public static final int ACTCODE_MAX                 = 10;

	// Return a string describing the actcode.

	public String get_actcode_as_string () {
		switch (actcode) {
		case ACTCODE_TRACK: return "ACTCODE_TRACK";
		case ACTCODE_FORECAST: return "ACTCODE_FORECAST";
		case ACTCODE_ANALYST: return "ACTCODE_ANALYST";
		case ACTCODE_PDL_UPDATE: return "ACTCODE_PDL_UPDATE";
		case ACTCODE_ERROR: return "ACTCODE_ERROR";
		case ACTCODE_WITHDRAWN: return "ACTCODE_WITHDRAWN";
		case ACTCODE_COMCAT_FAIL: return "ACTCODE_COMCAT_FAIL";
		case ACTCODE_FORESHOCK: return "ACTCODE_FORESHOCK";
		case ACTCODE_EXPIRED: return "ACTCODE_EXPIRED";
		case ACTCODE_STATUS_UPDATE: return "ACTCODE_STATUS_UPDATE";
		}
		return "ACTCODE_INVALID(" + actcode + ")";
	}

	// Time stamp for this timeline entry, in milliseconds since the epoch.
	// This must be strictly monotone increasing within the timeline.
	// For forecasts, this is generally the nominal forecast time (mainshock_time + forecast_lag).
	// For other entries, this may be the creation time, or may be manipulated to preserve ordering.

	public long action_time;


	//----- Status variables -----
	//
	// These hold the current status of the timeline.
	// They are generally copied from one timeline entry to the next,
	// modified as necessary to accomodate state changes.

	// Time at which this entry was created, in milliseconds since the epoch.

	public long entry_time;

	// Origin of this timeline.

	public int fc_origin;

	public static final int FCORIG_MIN               = 1;
	public static final int FCORIG_PDL               = 1;
		// Timeline created in response to PDL notification.
	public static final int FCORIG_ANALYST           = 2;
		// Timeline created in response to analyst intervention.
	public static final int FCORIG_SYNC              = 3;
		// Timeline created while synchronizing with ComCat, backup, or another server.
	public static final int FCORIG_UNKNOWN           = 4;
		// Timeline created for unknown reason.  This should only be the case when
		// the timeline is in an error state.
	public static final int FCORIG_MAX               = 4;

	// Return a string describing the fc_origin.

	public String get_fc_origin_as_string () {
		switch (fc_origin) {
		case FCORIG_PDL: return "FCORIG_PDL";
		case FCORIG_ANALYST: return "FCORIG_ANALYST";
		case FCORIG_SYNC: return "FCORIG_SYNC";
		case FCORIG_UNKNOWN: return "FCORIG_UNKNOWN";
		}
		return "FCORIG_INVALID(" + fc_origin + ")";
	}

	// Status of this timeline.

	public int fc_status;

	public static final int FCSTAT_MIN               = 1;
	public static final int FCSTAT_ACTIVE_NORMAL     = 1;
		// Normal active state.  Forecasts are being generated automatically,
		// according to the configured schedule.
	public static final int FCSTAT_ACTIVE_INTAKE     = 2;
		// Active state, with tentative invocation.  At the time of the next forecast,
		// a decision must be made on whether to continue the active state or stop.
	public static final int FCSTAT_STOP_EXPIRED      = 3;
		// Stopped state.  All forecasts have been completed.
	public static final int FCSTAT_STOP_ANALYST      = 4;
		// Stopped state.  An analyst intervention requested that no more forecasts
		// be generated.
	public static final int FCSTAT_STOP_FORESHOCK    = 5;
		// Stopped state.  An aftershock was found with greater magnitude than the
		// mainshock, and so this timeline is now considered to be a foreshock.
	public static final int FCSTAT_STOP_ERROR        = 6;
		// Stopped state.  An error was encountered that makes further use of this
		// timeline impossible.
	public static final int FCSTAT_STOP_WITHDRAWN    = 7;
		// Stopped state.  The automatic system decided to stop tracking this event,
		// due to it being of insufficient magnitude or outside the intake region.
	public static final int FCSTAT_STOP_COMCAT_FAIL  = 8;
		// Stopped state.  The system is unable to retrieve information from ComCat.
	public static final int FCSTAT_MAX               = 8;

	// Return a string describing the fc_status.

	public String get_fc_status_as_string () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_NORMAL: return "FCSTAT_ACTIVE_NORMAL";
		case FCSTAT_ACTIVE_INTAKE: return "FCSTAT_ACTIVE_INTAKE";
		case FCSTAT_STOP_EXPIRED: return "FCSTAT_STOP_EXPIRED";
		case FCSTAT_STOP_ANALYST: return "FCSTAT_STOP_ANALYST";
		case FCSTAT_STOP_FORESHOCK: return "FCSTAT_STOP_FORESHOCK";
		case FCSTAT_STOP_ERROR: return "FCSTAT_STOP_ERROR";
		case FCSTAT_STOP_WITHDRAWN: return "FCSTAT_STOP_WITHDRAWN";
		case FCSTAT_STOP_COMCAT_FAIL: return "FCSTAT_STOP_COMCAT_FAIL";
		}
		return "FCSTAT_INVALID(" + fc_status + ")";
	}

	// PDL reporting status, for the current or most recent forecast.

	public int pdl_status;

	public static final int PDLSTAT_MIN              = 1;
	public static final int PDLSTAT_NONE             = 1;
		// No PDL report was requested for the forecast.
	public static final int PDLSTAT_SUCCESS          = 2;
		// Forecast was successfully reported to PDL.
	public static final int PDLSTAT_FAILURE          = 3;
		// Forecast was permanently not reported to PDL due to excessive retries.
	public static final int PDLSTAT_PENDING          = 4;
		// Report is pending.  One or more attempts to report the forecast failed,
		// but the attempt will be repeated.
	public static final int PDLSTAT_SECONDARY        = 5;
		// Forecast needs to be reported to PDL, but it was not done because this
		// is a secondary server.
	public static final int PDLSTAT_CONFIRMED        = 6;
		// Forecast needs to be reported to PDL, but it was not done because this
		// is a secondary server, and it has been confirmed that the primary server
		// sent the report.
	public static final int PDLSTAT_BYPASSED         = 7;
		// Forecast needs to be reported to PDL, but it was not done because it
		// was too late or some other blocking condition.
	public static final int PDLSTAT_MAX              = 7;

	// Return a string describing the pdl_status.

	public String get_pdl_status_as_string () {
		switch (pdl_status) {
		case PDLSTAT_NONE: return "PDLSTAT_NONE";
		case PDLSTAT_SUCCESS: return "PDLSTAT_SUCCESS";
		case PDLSTAT_FAILURE: return "PDLSTAT_FAILURE";
		case PDLSTAT_PENDING: return "PDLSTAT_PENDING";
		case PDLSTAT_SECONDARY: return "PDLSTAT_SECONDARY";
		case PDLSTAT_CONFIRMED: return "PDLSTAT_CONFIRMED";
		case PDLSTAT_BYPASSED: return "PDLSTAT_BYPASSED";
		}
		return "PDLSTAT_INVALID(" + pdl_status + ")";
	}

	// If this event is a foreshock, this is the event_id of the mainshock.
	// Otherwise, this is an empty string "".
	// This field should contain "" unless fc_status == FCSTAT_STOP_FORESHOCK.

	public String foreshock_event_id;


	//----- Analyst data -----
	//
	// If an analyst has intervened in the timeline, these variables hold the analyst info.
	// These are copied from one timeline entry to the next, to preserve the
	// analyst selections, unless the analyst intervenes again.

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


	//----- Forecast data -----

	// Parameters for the current forecast.
	// This must be non-null if the timeline entry contains a forecast, or null otherwise.

	public ForecastParameters forecast_params;

	// Results for the current forecast.
	// This must be non-null if the timeline entry contains a forecast, or null otherwise.

	public ForecastResults forecast_results;


	//----- Scheduling data -----

	// Last-obtained mainshock time, in milliseconds since the epoch.
	// The value is 0L if the mainshock time is unknown.

	public long last_mainshock_time;

	// Time lag at which the last forecast occured, in milliseconds since the mainshock.
	// The value is -1L if there have been no prior forecasts.

	public long last_forecast_lag;

	// Time lag at which an extra forecast is requested, in milliseconds since the mainshock.
	// The value is -1L if there has been no extra forecast requested.
	// This is set to -1L when any forecast is issued, meaning that the request is considered
	// to be granted when any forecast is issued.

	public long extra_forecast_lag;




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public TimelineStatus () {}

	// Copy from another entry.

	public void copy_from (TimelineStatus other) {

		event_id            = other.event_id;
		actcode             = other.actcode;
		action_time         = other.action_time;

		entry_time          = other.entry_time;
		fc_origin           = other.fc_origin;
		fc_status           = other.fc_status;
		pdl_status          = other.pdl_status;
		foreshock_event_id  = other.foreshock_event_id;

		analyst_id          = other.analyst_id;
		analyst_remark      = other.analyst_remark;
		analyst_time        = other.analyst_time;
		analyst_params      = other.analyst_params;

		forecast_params     = other.forecast_params;
		forecast_results    = other.forecast_results;

		last_mainshock_time = other.last_mainshock_time;
		last_forecast_lag   = other.last_forecast_lag;
		extra_forecast_lag  = other.extra_forecast_lag;
		return;
	}




	//----- Service functions: Query -----

	// Return true if this timeline is in a state that requires a future forecast.

	public boolean is_forecast_state () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_NORMAL:
		case FCSTAT_ACTIVE_INTAKE:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that requires intake to be
	// reconsidered when a forecast is generated.

	public boolean is_intake_state () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_INTAKE:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that requires a PDL report retry.

	public boolean is_pdl_retry_state () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_NORMAL:
		case FCSTAT_ACTIVE_INTAKE:
		case FCSTAT_STOP_EXPIRED:
			if (pdl_status == PDLSTAT_PENDING) {
				return true;
			}
			break;
		}
		return false;
	}

	// Return true if there should be an associated catalog snapshot.

	public boolean has_catalog_snapshot () {
		switch (actcode) {
		case ACTCODE_FORECAST:
			if (forecast_results.catalog_result_avail) {
				return true;
			}
			break;
		}
		return false;
	}

	// Get the associated catalog snapshot, or null if none.
	// Note: The catalog is not marshaled, so this will onl obtain the catalog
	// when the forecast is first generated.

	public CompactEqkRupList get_catalog_snapshot () {
		switch (actcode) {
		case ACTCODE_FORECAST:
			if (forecast_results.catalog_result_avail) {
				return forecast_results.catalog_aftershocks;
			}
			break;
		}
		return null;
	}

	// Return true if this is the first timeline entry.

	public boolean is_first_entry () {
		switch (actcode) {
		case ACTCODE_TRACK:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that can be converted to the
	// normal active state by changing fc_status.

	public boolean is_convertible_to_normal () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_INTAKE:
		case FCSTAT_STOP_WITHDRAWN:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that an analyst can
	// start by changing fc_status (with or without updating analyst data).

	public boolean can_analyst_start () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_INTAKE:
		case FCSTAT_STOP_ANALYST:
		case FCSTAT_STOP_WITHDRAWN:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that an analyst can
	// stop by changing fc_status (with or without updating analyst data).

	public boolean can_analyst_stop () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_NORMAL:
		case FCSTAT_ACTIVE_INTAKE:
		case FCSTAT_STOP_WITHDRAWN:
			return true;
		}
		return false;
	}

	// Return true if this timeline is in a state that an analyst can
	// update (withough changing fc_status).

	public boolean can_analyst_update () {
		switch (fc_status) {
		case FCSTAT_ACTIVE_NORMAL:
		case FCSTAT_ACTIVE_INTAKE:
		case FCSTAT_STOP_ANALYST:
		case FCSTAT_STOP_WITHDRAWN:
			return true;
		}
		return false;
	}




	//----- Service functions: Insertion -----

	// Set the fc_status variable.

	public void set_fc_status (int the_fc_status) {
		fc_status = the_fc_status;
		return;
	}

	// Set the pdl_status variable.

	public void set_pdl_status (int the_pdl_status) {
		pdl_status = the_pdl_status;
		return;
	}

	// Set the analyst data variables.
	// If the analyst data contains a mainshock time, save it too.

	public void set_analyst_data (String the_analyst_id, String the_analyst_remark,
		long the_analyst_time, ForecastParameters the_analyst_params, long the_extra_forecast_lag) {
		
		analyst_id = the_analyst_id;
		analyst_remark = the_analyst_remark;
		analyst_time = the_analyst_time;
		analyst_params = the_analyst_params;
		extra_forecast_lag = the_extra_forecast_lag;

		if (the_analyst_params != null) {
			if (the_analyst_params.mainshock_fetch_meth == ForecastParameters.FETCH_METH_ANALYST
					&& the_analyst_params.mainshock_avail) {
				last_mainshock_time = the_analyst_params.mainshock_time;
			}
		}

		return;
	}




	//----- Service functions: State setting -----

	// Set the state to error.

	public void set_state_error (String the_event_id, long the_action_time, long the_entry_time) {

		event_id            = the_event_id;
		actcode             = ACTCODE_ERROR;
		action_time         = the_action_time;

		entry_time          = the_entry_time;
		fc_origin           = FCORIG_UNKNOWN;
		fc_status           = FCSTAT_STOP_ERROR;
		pdl_status          = PDLSTAT_NONE;
		foreshock_event_id  = "";

		analyst_id          = "";
		analyst_remark      = "";
		analyst_time        = 0L;
		analyst_params      = null;

		forecast_params     = null;
		forecast_results    = null;

		last_mainshock_time = 0L;
		last_forecast_lag   = -1L;
		extra_forecast_lag  = -1L;
		return;
	}

	// Set the state to ComCat failure.

	public void set_state_comcat_fail (long the_entry_time) {

		//event_id            = kept;
		actcode             = ACTCODE_COMCAT_FAIL;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		fc_status           = FCSTAT_STOP_COMCAT_FAIL;
		pdl_status          = PDLSTAT_NONE;
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to foreshock.

	public void set_state_foreshock (long the_entry_time, String the_foreshock_event_id) {

		//event_id            = kept;
		actcode             = ACTCODE_FORESHOCK;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		fc_status           = FCSTAT_STOP_FORESHOCK;
		pdl_status          = PDLSTAT_NONE;
		foreshock_event_id  = the_foreshock_event_id;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to withdrawn.

	public void set_state_withdrawn (long the_entry_time, long the_last_mainshock_time) {

		//event_id            = kept;
		actcode             = ACTCODE_WITHDRAWN;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		fc_status           = FCSTAT_STOP_WITHDRAWN;
		pdl_status          = PDLSTAT_NONE;
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		last_mainshock_time = the_last_mainshock_time;
		//last_forecast_lag   = kept;
		extra_forecast_lag  = -1L;
		return;
	}

	// Set the state to forecast.
	//
	// Note: fc_status is set to FCSTAT_ACTIVE_NORMAL.  The caller can change fc_status
	// to FCSTAT_STOP_EXPIRED if there are no further forecasts scheduled.
	//
	// Note: pdl_status is set to either PDLSTAT_NONE or PDLSTAT_PENDING depending on the
	// return value of the_forecast_results.get_pdl_model().  The caller should call
	// the_forecast_results.pick_pdl_model() (or set the xxxx_pdl flags in some other way)
	// before calling this function.  If PDLSTAT_PENDING results, the caller will likely
	// want to take further action (such as an immediate attempt to send a report to PDL)
	// and make a futher update to pdl_status.

	public void set_state_forecast (long the_entry_time, ActionConfig action_config,
			ForecastParameters the_forecast_params, ForecastResults the_forecast_results) {

		//long min_lag = the_forecast_params.forecast_lag + action_config.get_forecast_min_gap();
		//long next_forecast_lag = action_config.get_next_forecast_lag (min_lag);
		long forecast_time = the_forecast_params.forecast_lag + the_forecast_params.mainshock_time;
		String pdl_json = the_forecast_results.get_pdl_model();

		//event_id            = kept;
		actcode             = ACTCODE_FORECAST;
		action_time         = action_config.floor_unit_lag (forecast_time, action_time);

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		//fc_status           = ((next_forecast_lag >= 0L) ? FCSTAT_ACTIVE_NORMAL : FCSTAT_STOP_EXPIRED);
		fc_status           = FCSTAT_ACTIVE_NORMAL;
		pdl_status          = ((pdl_json == null) ? PDLSTAT_NONE : PDLSTAT_PENDING);
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = the_forecast_params;
		forecast_results    = the_forecast_results;

		last_mainshock_time = the_forecast_params.mainshock_time;
		last_forecast_lag   = the_forecast_params.forecast_lag;
		extra_forecast_lag  = -1L;
		return;
	}

	// Set the state to PDL update.
	//
	// Note: This keeps the forecast parameters and results, so that it is possible
	// to send a PDL report for a forecast which originally didn't send a report.

	public void set_state_pdl_update (long the_entry_time, int the_pdl_status) {

		//event_id            = kept;
		actcode             = ACTCODE_PDL_UPDATE;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		//fc_status           = kpet;
		pdl_status          = the_pdl_status;
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		//forecast_params     = kept;
		//forecast_results    = kept;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to expired.

	public void set_state_expired (long the_entry_time) {

		//event_id            = kept;
		actcode             = ACTCODE_EXPIRED;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		if (fc_status == FCSTAT_ACTIVE_NORMAL || fc_status == FCSTAT_ACTIVE_INTAKE) {
			fc_status           = FCSTAT_STOP_EXPIRED;
		}
		if (pdl_status == PDLSTAT_PENDING) {
			pdl_status          = PDLSTAT_NONE;
		}
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to status update.
	//
	// Note: Sets pdl_status to PDLSTAT_NONE.

	public void set_state_status_update (long the_entry_time, int the_fc_status) {

		//event_id            = kept;
		actcode             = ACTCODE_STATUS_UPDATE;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		fc_status           = the_fc_status;
		pdl_status          = PDLSTAT_NONE;
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to analyst intervention.
	//
	// Note: Sets pdl_status to PDLSTAT_NONE.
	//
	// Note: The caller may use set_fc_status(), set_pdl_status(), and set_analyst_data()
	// to make changes to forecast status, PDL status, and analyst data.

	public void set_state_analyst_intervention (long the_entry_time) {

		//event_id            = kept;
		actcode             = ACTCODE_ANALYST;
		action_time         = action_time + 1L;

		entry_time          = the_entry_time;
		//fc_origin           = kept;
		//fc_status           = kept;
		pdl_status          = PDLSTAT_NONE;
		//foreshock_event_id  = kept;

		//analyst_id          = kept;
		//analyst_remark      = kept;
		//analyst_time        = kept;
		//analyst_params      = kept;

		forecast_params     = null;
		forecast_results    = null;

		//last_mainshock_time = kept;
		//last_forecast_lag   = kept;
		//extra_forecast_lag  = kept;
		return;
	}

	// Set the state to track.
	//
	// Note: Sets pdl_status to PDLSTAT_NONE.
	//
	// Note: The caller may use set_fc_status(), set_pdl_status(), and set_analyst_data()
	// to make changes to forecast status, PDL status, and analyst data.

	public void set_state_track (long the_entry_time, ActionConfig action_config,
			String the_event_id, long the_last_mainshock_time, int the_fc_origin, int the_fc_status) {

		event_id            = the_event_id;
		actcode             = ACTCODE_TRACK;
		action_time         = action_config.floor_unit_lag (the_last_mainshock_time);

		entry_time          = the_entry_time;
		fc_origin           = the_fc_origin;
		fc_status           = the_fc_status;
		pdl_status          = PDLSTAT_NONE;
		foreshock_event_id  = "";

		analyst_id          = "";
		analyst_remark      = "";
		analyst_time        = 0L;
		analyst_params      = null;

		forecast_params     = null;
		forecast_results    = null;

		last_mainshock_time = the_last_mainshock_time;
		last_forecast_lag   = -1L;
		extra_forecast_lag  = -1L;
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 27001;

	private static final String M_VERSION_NAME = "TimelineStatus";

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		writer.marshalLong                      ("entry_time"         , entry_time         );
		writer.marshalInt                       ("fc_origin"          , fc_origin          );
		writer.marshalInt                       ("fc_status"          , fc_status          );
		writer.marshalInt                       ("pdl_status"         , pdl_status         );
		writer.marshalString                    ("foreshock_event_id" , foreshock_event_id );

		writer.marshalString                    ("analyst_id"         , analyst_id         );
		writer.marshalString                    ("analyst_remark"     , analyst_remark     );
		writer.marshalLong                      ("analyst_time"       , analyst_time       );
		ForecastParameters.marshal_poly (writer, "analyst_params"     , analyst_params     );

		ForecastParameters.marshal_poly (writer, "forecast_params"    , forecast_params    );
		ForecastResults.marshal_poly    (writer, "forecast_results"   , forecast_results   );

		writer.marshalLong                      ("last_mainshock_time", last_mainshock_time);
		writer.marshalLong                      ("last_forecast_lag"  , last_forecast_lag  );
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

		entry_time          = reader.unmarshalLong                      ("entry_time"         );
		fc_origin           = reader.unmarshalInt                       ("fc_origin"          , FCORIG_MIN, FCORIG_MAX);
		fc_status           = reader.unmarshalInt                       ("fc_status"          , FCSTAT_MIN, FCSTAT_MAX);
		pdl_status          = reader.unmarshalInt                       ("pdl_status"         , PDLSTAT_MIN, PDLSTAT_MAX);
		foreshock_event_id  = reader.unmarshalString                    ("foreshock_event_id" );

		analyst_id          = reader.unmarshalString                    ("analyst_id"         );
		analyst_remark      = reader.unmarshalString                    ("analyst_remark"     );
		analyst_time        = reader.unmarshalLong                      ("analyst_time"       );
		analyst_params      = ForecastParameters.unmarshal_poly (reader, "analyst_params"     );

		forecast_params     = ForecastParameters.unmarshal_poly (reader, "forecast_params"    );
		forecast_results    = ForecastResults.unmarshal_poly    (reader, "forecast_results"   );

		last_mainshock_time = reader.unmarshalLong                      ("last_mainshock_time");
		last_forecast_lag   = reader.unmarshalLong                      ("last_forecast_lag"  );
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
	public TimelineStatus unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Unmarshal object, for a timeline entry.

	@Override
	public TimelineStatus unmarshal_timeline (TimelineEntry tentry) {
		try {
			if (tentry.get_actcode() == ACTCODE_ERROR) {
				throw new MarshalException("Timeline is in an unrecoverable error state, event_id = " + tentry.get_event_id());
			}

			unmarshal (tentry.get_details(), null);

			event_id = tentry.get_event_id();
			actcode = tentry.get_actcode();
			action_time = tentry.get_action_time();

			if (actcode < ACTCODE_MIN || actcode > ACTCODE_MAX) {
				throw new MarshalException("Timeline entry contains invalid action code, event_id = " + tentry.get_event_id() + ", actcode = " + actcode);
			}

		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling timeline entry payload\n" + tentry.toString() + "\nDump:\n" + tentry.dump_details(), e);
		}
		return this;
	}

}
