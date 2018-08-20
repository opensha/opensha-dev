package scratch.aftershockStatistics.aafs;

import java.util.GregorianCalendar;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;
import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.SphLatLon;
import scratch.aftershockStatistics.util.SphRegion;

import scratch.aftershockStatistics.AftershockStatsCalc;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.GenericRJ_Parameters;
import scratch.aftershockStatistics.MagCompPage_Parameters;
import scratch.aftershockStatistics.RJ_AftershockModel;
import scratch.aftershockStatistics.RJ_AftershockModel_Bayesian;
import scratch.aftershockStatistics.RJ_AftershockModel_Generic;
import scratch.aftershockStatistics.RJ_AftershockModel_SequenceSpecific;
import scratch.aftershockStatistics.RJ_Summary;
import scratch.aftershockStatistics.RJ_Summary_Bayesian;
import scratch.aftershockStatistics.RJ_Summary_Generic;
import scratch.aftershockStatistics.RJ_Summary_SequenceSpecific;
import scratch.aftershockStatistics.SeqSpecRJ_Parameters;
import scratch.aftershockStatistics.USGS_AftershockForecast;

import static scratch.aftershockStatistics.aafs.ForecastParameters.CALC_METH_AUTO_PDL;
import static scratch.aftershockStatistics.aafs.ForecastParameters.CALC_METH_AUTO_NO_PDL;
import static scratch.aftershockStatistics.aafs.ForecastParameters.CALC_METH_SUPPRESS;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;

import org.json.simple.JSONObject;

/**
 * Results of a forecast.
 * Author: Michael Barall 04/26/2018.
 *
 * All fields are public, since there is little benefit to having lots of getters and setters.
 */
public class ForecastResults {

	//----- Constants -----

	// Standard values of the advisory duration.

	public static final long ADVISORY_LAG_DAY   = 86400000L;	// 1 day
	public static final long ADVISORY_LAG_WEEK  = 604800000L;	// 1 week = 7 days
	public static final long ADVISORY_LAG_MONTH = 2592000000L;	// 1 month = 30 days
	public static final long ADVISORY_LAG_YEAR  = 31536000000L;	// 1 year = 365 days


	//----- Root parameters -----

	// Time results were prepared, in milliseconds since the epoch.
	// This is used as the start time for the forecasts, and typically equals the mainshock time plus the forecast lag.

	public long result_time = 0L;

	// Advisory duration, in milliseconds.

	public long advisory_lag = ADVISORY_LAG_WEEK;

	// Injectable text for PDL JSON files, or "" for none.

	public String injectable_text = "";


	//----- Catalog results -----

	// Catalog result available flag.

	public boolean catalog_result_avail = false;

	// Start time of aftershock sequence, in milliseconds since the epoch.
	// This is not the time of the first aftershock, because there may be an interval of no aftershock.

	public long catalog_start_time = 0L;

	// End time of aftershock sequence, in milliseconds since the epoch.
	// This is not the time of the last aftershock, because there may be an interval of no aftershock.

	public long catalog_end_time = 0L;

	// Number of aftershocks.

	public int catalog_eqk_count = 0;

	// Maximum magnitude of any aftershock (0.0 if there are no aftershocks).

	public double catalog_max_mag = 0.0;

	// Event id of the aftershock with maximum magnitude ("" if there are no aftershocks).

	public String catalog_max_event_id = "";

	// catalog_aftershocks - List of aftershocks.
	// Note: This field is not marshaled, because aftershock lists are stored in separate database records.

	public CompactEqkRupList catalog_aftershocks = null;

	// catalog_comcat_aftershocks - List of aftershocks, as returned by ComCat.
	// Note: This field is not marshaled, and is not rebuilt, because it is intended to be used
	// only immediately after the results are calculated;  mainly for finding foreshocks.

	public ObsEqkRupList catalog_comcat_aftershocks = null;

	// set_default_catalog_results - Set catalog results to default values.

	public void set_default_catalog_results () {
		catalog_start_time = 0L;
		catalog_end_time = 0L;
		catalog_eqk_count = 0;
		catalog_max_mag = 0.0;
		catalog_max_event_id = "";
		catalog_aftershocks = null;
		catalog_comcat_aftershocks = null;
		return;
	}

	// calc_catalog_results - Calculate catalog results.

	public void calc_catalog_results (ForecastMainshock fcmain, ForecastParameters params) {

		// Parameters must have mainshock, aftershock search region

		if (!( fcmain.mainshock_avail && params.aftershock_search_avail )) {
			set_default_catalog_results();
			catalog_result_avail = false;
			return;
		}

		// Retrieve list of aftershocks in the search region

		ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
		//ObsEqkRupList catalog_comcat_aftershocks;		// if this isn't an object field

		try {
			ComcatAccessor accessor = new ComcatAccessor();
			catalog_comcat_aftershocks = accessor.fetchAftershocks(mainshock, params.min_days, params.max_days,
				params.min_depth, params.max_depth, params.aftershock_search_region, false, params.min_mag);
		} catch (Exception e) {
			throw new RuntimeException("ForecastResults.calc_catalog_results: Comcat exception", e);
		}

		// Save catalog and info

		long eventTime = mainshock.getOriginTime();
		catalog_start_time = eventTime + (long)(params.min_days * ComcatAccessor.day_millis);
		catalog_end_time = eventTime + (long)(params.max_days * ComcatAccessor.day_millis);

		catalog_eqk_count = catalog_comcat_aftershocks.size();
		catalog_aftershocks = new CompactEqkRupList (catalog_comcat_aftershocks);

		if (catalog_eqk_count == 0) {
			catalog_max_mag = 0.0;
			catalog_max_event_id = "";
		} else {
			catalog_max_mag = -Double.MAX_VALUE;
			for (ObsEqkRupture rup : catalog_comcat_aftershocks) {
				double mag = rup.getMag();
				if (mag > catalog_max_mag) {
					catalog_max_mag = mag;
					catalog_max_event_id = rup.getEventId();
				}
			}
		}

		catalog_result_avail = true;
		return;
	}

	// rebuild_catalog_results - Rebuild transient catalog results.

	public void rebuild_catalog_results (ForecastMainshock fcmain, ForecastParameters params, CompactEqkRupList the_catalog_aftershocks) {

		// If there are results to rebuild ...

		if (catalog_result_avail) {

			// Parameters must have mainshock, aftershock search region

			if (!( fcmain.mainshock_avail && params.aftershock_search_avail )) {
				throw new RuntimeException("ForecastResults.rebuild_catalog_results: Invalid preconditions");
			}

			// Check for supplied catalog

			if (!( the_catalog_aftershocks != null )) {
				throw new RuntimeException("ForecastResults.rebuild_catalog_results: No aftershock catalog supplied");
			}

			if (!( the_catalog_aftershocks.size() == catalog_eqk_count )) {
				throw new RuntimeException("ForecastResults.rebuild_catalog_results: Aftershock catalog size mismatch, expecting " + catalog_eqk_count + ", got " + the_catalog_aftershocks.size());
			}

			// Save catalog

			catalog_aftershocks = the_catalog_aftershocks;
		}

		return;
	}


	//----- Generic results -----

	// Generic result available flag.

	public boolean generic_result_avail = false;

	// Generic results summary.

	public RJ_Summary_Generic generic_summary = null;

	// Generic results JSON.

	public String generic_json = "";

	// True if results sent to PDL.

	public boolean generic_pdl = false;

	// Generic aftershock model.
	// This field is not marshaled.

	RJ_AftershockModel_Generic generic_model = null;

	// set_default_generic_results - Set generic results to default values.

	public void set_default_generic_results () {
		generic_summary = null;
		generic_json = "";
		generic_pdl = false;
		generic_model = null;
		return;
	}

	// calc_generic_results - Calculate generic results.

	public void calc_generic_results (ForecastMainshock fcmain, ForecastParameters params) {

		// We need to have catalog results, mainshock parameters, and generic parameters

		if (!( (params.generic_calc_meth != CALC_METH_SUPPRESS)
				&& catalog_result_avail
				&& fcmain.mainshock_avail 
				&& params.generic_avail )) {
			set_default_generic_results();
			generic_result_avail = false;
			return;
		}

		try {

			// Build the generic model

			ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
			generic_model = new RJ_AftershockModel_Generic (mainshock.getMag(), params.generic_params);

			// Save the summary

			generic_summary = new RJ_Summary_Generic (generic_model);

			// Build the forecast

			GregorianCalendar eventDate = new GregorianCalendar();
			eventDate.setTimeInMillis(mainshock.getOriginTime());
			GregorianCalendar startDate = new GregorianCalendar();
			startDate.setTimeInMillis(Math.max (result_time, mainshock.getOriginTime() + 1000L));
			USGS_AftershockForecast forecast = new USGS_AftershockForecast (generic_model, catalog_aftershocks, eventDate, startDate);

			if (advisory_lag >= ADVISORY_LAG_YEAR) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_YEAR);
			} else if (advisory_lag >= ADVISORY_LAG_MONTH) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_MONTH);
			} else if (advisory_lag >= ADVISORY_LAG_WEEK) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_WEEK);
			} else {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_DAY);
			}

			String the_injectable_text = injectable_text;
			if (the_injectable_text.length() == 0) {
				the_injectable_text = null;		// convention for USGS_AftershockForecast
			}
			forecast.setInjectableText (the_injectable_text);

			// Get the JSON String

			JSONObject json = forecast.buildJSON(result_time);
			if (json == null) {
				throw new RuntimeException("ForecastResults.calc_generic_results: Unable to generate JSON");
			}
			generic_json = json.toJSONString();

		} catch (Exception e) {
			throw new RuntimeException("ForecastResults.calc_generic_results: Exception building generic forecast", e);
		}

		// Done

		generic_pdl = false;
		generic_result_avail = true;
		return;
	}

	// rebuild_generic_results - Rebuild transient generic results.

	public void rebuild_generic_results (ForecastMainshock fcmain, ForecastParameters params) {

		// If there are results to rebuild ...

		if (generic_result_avail) {

			// We need to have catalog results, mainshock parameters, and generic parameters

			if (!( (params.generic_calc_meth != CALC_METH_SUPPRESS)
					&& catalog_result_avail
					&& fcmain.mainshock_avail 
					&& params.generic_avail )) {
				throw new RuntimeException("ForecastResults.rebuild_generic_results: Invalid preconditions");
			}

			try {

				// Build the generic model

				ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
				generic_model = new RJ_AftershockModel_Generic (mainshock.getMag(), params.generic_params);

			} catch (Exception e) {
				throw new RuntimeException("ForecastResults.rebuild_generic_results: Exception building generic forecast", e);
			}
		}

		return;
	}


	//----- Sequence specific results -----

	// Sequence specific result available flag.

	public boolean seq_spec_result_avail = false;

	// Sequence specific results summary.

	public RJ_Summary_SequenceSpecific seq_spec_summary = null;

	// Sequence specific results JSON.

	public String seq_spec_json = "";

	// True if results sent to PDL.

	public boolean seq_spec_pdl = false;

	// Sequence specific aftershock model.
	// This field is not marshaled.

	RJ_AftershockModel_SequenceSpecific seq_spec_model = null;

	// set_default_seq_spec_results - Set sequence specific results to default values.

	public void set_default_seq_spec_results () {
		seq_spec_summary = null;
		seq_spec_json = "";
		seq_spec_pdl = false;
		seq_spec_model = null;
		return;
	}

	// calc_seq_spec_results - Calculate sequence specific results.

	public void calc_seq_spec_results (ForecastMainshock fcmain, ForecastParameters params, boolean f_seq_spec) {

		// We need to have catalog results, mainshock parameters, magnitude of completeness parameters, and sequence specific parameters

		if (!( f_seq_spec
				&& (params.seq_spec_calc_meth != CALC_METH_SUPPRESS)
				&& catalog_result_avail
				&& fcmain.mainshock_avail 
				&& params.mag_comp_avail
				&& params.seq_spec_avail )) {
			set_default_seq_spec_results();
			seq_spec_result_avail = false;
			return;
		}

		try {

			// Build the sequence specific model

			ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
			seq_spec_model = new RJ_AftershockModel_SequenceSpecific (mainshock, catalog_aftershocks,
				params.min_days, params.max_days, params.mag_comp_params, params.seq_spec_params);

			// Save the summary

			seq_spec_summary = new RJ_Summary_SequenceSpecific (seq_spec_model);

			// Build the forecast

			GregorianCalendar eventDate = new GregorianCalendar();
			eventDate.setTimeInMillis(mainshock.getOriginTime());
			GregorianCalendar startDate = new GregorianCalendar();
			startDate.setTimeInMillis(Math.max (result_time, mainshock.getOriginTime() + 1000L));
			USGS_AftershockForecast forecast = new USGS_AftershockForecast (seq_spec_model, catalog_aftershocks, eventDate, startDate);

			if (advisory_lag >= ADVISORY_LAG_YEAR) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_YEAR);
			} else if (advisory_lag >= ADVISORY_LAG_MONTH) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_MONTH);
			} else if (advisory_lag >= ADVISORY_LAG_WEEK) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_WEEK);
			} else {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_DAY);
			}

			String the_injectable_text = injectable_text;
			if (the_injectable_text.length() == 0) {
				the_injectable_text = null;		// convention for USGS_AftershockForecast
			}
			forecast.setInjectableText (the_injectable_text);

			// Get the JSON String

			JSONObject json = forecast.buildJSON(result_time);
			if (json == null) {
				throw new RuntimeException("ForecastResults.calc_seq_spec_results: Unable to generate JSON");
			}
			seq_spec_json = json.toJSONString();

		} catch (Exception e) {
			throw new RuntimeException("ForecastResults.calc_seq_spec_results: Exception building sequence specific forecast", e);
		}

		// Done

		seq_spec_pdl = false;
		seq_spec_result_avail = true;
		return;
	}

	// rebuild_seq_spec_results - Rebuild transient sequence specific results.

	public void rebuild_seq_spec_results (ForecastMainshock fcmain, ForecastParameters params) {

		// If there are results to rebuild ...

		if (seq_spec_result_avail) {

			// We need to have catalog results, mainshock parameters, magnitude of completeness parameters, and sequence specific parameters

			if (!( (params.seq_spec_calc_meth != CALC_METH_SUPPRESS)
					&& catalog_result_avail
					&& fcmain.mainshock_avail 
					&& params.mag_comp_avail
					&& params.seq_spec_avail )) {
				throw new RuntimeException("ForecastResults.rebuild_seq_spec_results: Invalid preconditions");
			}

			try {

				// Build the sequence specific model

				ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
				seq_spec_model = new RJ_AftershockModel_SequenceSpecific (mainshock, catalog_aftershocks,
					params.min_days, params.max_days, params.mag_comp_params, params.seq_spec_params);

			} catch (Exception e) {
				throw new RuntimeException("ForecastResults.rebuild_seq_spec_results: Exception building sequence specific forecast", e);
			}
		}

		return;
	}


	//----- Bayesian results -----

	// Bayesian result available flag.

	public boolean bayesian_result_avail = false;

	// Bayesian results summary.

	public RJ_Summary_Bayesian bayesian_summary = null;

	// Bayesian results JSON.

	public String bayesian_json = "";

	// True if results sent to PDL.

	public boolean bayesian_pdl = false;

	// Bayesian aftershock model.
	// This field is not marshaled.

	RJ_AftershockModel_Bayesian bayesian_model = null;

	// set_default_bayesian_results - Set bayesian results to default values.

	public void set_default_bayesian_results () {
		bayesian_summary = null;
		bayesian_json = "";
		bayesian_pdl = false;
		bayesian_model = null;
		return;
	}

	// calc_bayesian_results - Calculate bayesian results.

	public void calc_bayesian_results (ForecastMainshock fcmain, ForecastParameters params) {

		// We need to have catalog results, mainshock parameters, compatible generic and sequence specific models

		if (!( (params.bayesian_calc_meth != CALC_METH_SUPPRESS)
				&& catalog_result_avail
				&& fcmain.mainshock_avail 
				&& generic_result_avail
				&& seq_spec_result_avail
				&& RJ_AftershockModel_Bayesian.areModelsEquivalent(generic_model, seq_spec_model) )) {
			set_default_bayesian_results();
			bayesian_result_avail = false;
			return;
		}

		try {

			// Build the bayesian model

			ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
			bayesian_model = new RJ_AftershockModel_Bayesian (generic_model, seq_spec_model);

			// Save the summary

			bayesian_summary = new RJ_Summary_Bayesian (bayesian_model);

			// Build the forecast

			GregorianCalendar eventDate = new GregorianCalendar();
			eventDate.setTimeInMillis(mainshock.getOriginTime());
			GregorianCalendar startDate = new GregorianCalendar();
			startDate.setTimeInMillis(Math.max (result_time, mainshock.getOriginTime() + 1000L));
			USGS_AftershockForecast forecast = new USGS_AftershockForecast (bayesian_model, catalog_aftershocks, eventDate, startDate);

			if (advisory_lag >= ADVISORY_LAG_YEAR) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_YEAR);
			} else if (advisory_lag >= ADVISORY_LAG_MONTH) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_MONTH);
			} else if (advisory_lag >= ADVISORY_LAG_WEEK) {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_WEEK);
			} else {
				forecast.setAdvisoryDuration (USGS_AftershockForecast.Duration.ONE_DAY);
			}

			String the_injectable_text = injectable_text;
			if (the_injectable_text.length() == 0) {
				the_injectable_text = null;		// convention for USGS_AftershockForecast
			}
			forecast.setInjectableText (the_injectable_text);

			// Get the JSON String

			JSONObject json = forecast.buildJSON(result_time);
			if (json == null) {
				throw new RuntimeException("ForecastResults.calc_bayesian_results: Unable to generate JSON");
			}
			bayesian_json = json.toJSONString();

		} catch (Exception e) {
			//throw new RuntimeException("ForecastResults.calc_bayesian_results: Exception building bayesian forecast", e);

			// In case of any error, just don't to the Bayesian

			set_default_bayesian_results();
			bayesian_result_avail = false;
			return;
		}

		// Done

		bayesian_pdl = false;
		bayesian_result_avail = true;
		return;
	}

	// rebuild_bayesian_results - Rebuild transient bayesian results.

	public void rebuild_bayesian_results (ForecastMainshock fcmain, ForecastParameters params) {

		// If there are results to rebuild ...

		if (bayesian_result_avail) {

			// We need to have catalog results, mainshock parameters, compatible generic and sequence specific models

			if (!( (params.bayesian_calc_meth != CALC_METH_SUPPRESS)
					&& catalog_result_avail
					&& fcmain.mainshock_avail 
					&& generic_result_avail
					&& seq_spec_result_avail
					&& RJ_AftershockModel_Bayesian.areModelsEquivalent(generic_model, seq_spec_model) )) {
				throw new RuntimeException("ForecastResults.rebuild_bayesian_results: Invalid preconditions");
			}

			try {

				// Build the bayesian model

				ObsEqkRupture mainshock = fcmain.get_eqk_rupture();
				bayesian_model = new RJ_AftershockModel_Bayesian (generic_model, seq_spec_model);

			} catch (Exception e) {
				throw new RuntimeException("ForecastResults.rebuild_bayesian_results: Exception building bayesian forecast", e);
			}
		}

		return;
	}


	//----- Construction -----

	// Default constructor.

	public ForecastResults () {}

	// Calculate all results.
	// If f_seq_spec is false, then sequence specific results are not calculated.

	public void calc_all (long the_result_time, long the_advisory_lag, String the_injectable_text, ForecastMainshock fcmain, ForecastParameters params, boolean f_seq_spec) {
		result_time = the_result_time;
		advisory_lag = the_advisory_lag;
		injectable_text = ((the_injectable_text == null) ? "" : the_injectable_text);
		calc_catalog_results (fcmain, params);
		calc_generic_results (fcmain, params);
		calc_seq_spec_results (fcmain, params, f_seq_spec);
		calc_bayesian_results (fcmain, params);
		return;
	}

	// Rebuild all transient results.

	public void rebuild_all (ForecastMainshock fcmain, ForecastParameters params, CompactEqkRupList the_catalog_aftershocks) {
		rebuild_catalog_results (fcmain, params, the_catalog_aftershocks);
		rebuild_generic_results (fcmain, params);
		rebuild_seq_spec_results (fcmain, params);
		rebuild_bayesian_results (fcmain, params);
		return;
	}

	// Pick one of the models to be sent to PDL, and set the corresponding xxxx_pdl flag.
	// Models are picked in priority order: Bayesian, sequence specific, and generic.
	// Returns the JSON to be sent to PDL, or null if none.

	public String pick_pdl_model () {

		if (bayesian_result_avail) {
			if (bayesian_json.length() > 0) {
				bayesian_pdl = true;
				return bayesian_json;
			}
		}

		if (seq_spec_result_avail) {
			if (seq_spec_json.length() > 0) {
				seq_spec_pdl = true;
				return seq_spec_json;
			}
		}

		if (generic_result_avail) {
			if (generic_json.length() > 0) {
				generic_pdl = true;
				return generic_json;
			}
		}
	
		return null;
	}

	// Get the model prevously picked to be sent to PDL.
	// This function looks at the xxxx_pdl flags to find the model to return.
	// Returns the JSON to be sent to PDL, or null if none.

	public String get_pdl_model () {

		if (bayesian_pdl) {
			return bayesian_json;
		}

		if (seq_spec_pdl) {
			return seq_spec_json;
		}

		if (generic_pdl) {
			return generic_json;
		}
	
		return null;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();

		result.append ("ForecastResults:" + "\n");

		result.append ("result_time = " + result_time + "\n");
		result.append ("advisory_lag = " + advisory_lag + "\n");
		result.append ("injectable_text = " + injectable_text + "\n");

		result.append ("catalog_result_avail = " + catalog_result_avail + "\n");
		if (catalog_result_avail) {
			result.append ("catalog_start_time = " + catalog_start_time + "\n");
			result.append ("catalog_end_time = " + catalog_end_time + "\n");
			result.append ("catalog_eqk_count = " + catalog_eqk_count + "\n");
			result.append ("catalog_max_mag = " + catalog_max_mag + "\n");
			result.append ("catalog_max_event_id = " + catalog_max_event_id + "\n");
			result.append ("catalog_aftershocks = " + ((catalog_aftershocks == null) ? "null" : "available") + "\n");
			result.append ("catalog_comcat_aftershocks = " + ((catalog_comcat_aftershocks == null) ? "null" : "available") + "\n");
		}

		result.append ("generic_result_avail = " + generic_result_avail + "\n");
		if (generic_result_avail) {
			result.append ("generic_summary:\n" + generic_summary.toString() + "\n");
			result.append ("generic_json = " + generic_json + "\n");
			result.append ("generic_pdl = " + generic_pdl + "\n");
			result.append ("generic_model = " + ((generic_model == null) ? "null" : "available") + "\n");
		}

		result.append ("seq_spec_result_avail = " + seq_spec_result_avail + "\n");
		if (seq_spec_result_avail) {
			result.append ("seq_spec_summary:\n" + seq_spec_summary.toString() + "\n");
			result.append ("seq_spec_json = " + seq_spec_json + "\n");
			result.append ("seq_spec_pdl = " + seq_spec_pdl + "\n");
			result.append ("seq_spec_model = " + ((seq_spec_model == null) ? "null" : "available") + "\n");
		}

		result.append ("bayesian_result_avail = " + bayesian_result_avail + "\n");
		if (bayesian_result_avail) {
			result.append ("bayesian_summary:\n" + bayesian_summary.toString() + "\n");
			result.append ("bayesian_json = " + bayesian_json + "\n");
			result.append ("bayesian_pdl = " + bayesian_pdl + "\n");
			result.append ("bayesian_model = " + ((bayesian_model == null) ? "null" : "available") + "\n");
		}

		return result.toString();
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 23001;

	private static final String M_VERSION_NAME = "ForecastResults";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 23000;
	protected static final int MARSHAL_FCAST_RESULT = 23001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_FCAST_RESULT;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalLong   ("result_time"    , result_time    );
		writer.marshalLong   ("advisory_lag"   , advisory_lag   );
		writer.marshalString ("injectable_text", injectable_text);

		writer.marshalBoolean ("catalog_result_avail", catalog_result_avail);
		if (catalog_result_avail) {
			writer.marshalLong   ("catalog_start_time"  , catalog_start_time  );
			writer.marshalLong   ("catalog_end_time"    , catalog_end_time    );
			writer.marshalInt    ("catalog_eqk_count"   , catalog_eqk_count   );
			writer.marshalDouble ("catalog_max_mag"     , catalog_max_mag     );
			writer.marshalString ("catalog_max_event_id", catalog_max_event_id);
		}

		writer.marshalBoolean ("generic_result_avail", generic_result_avail);
		if (generic_result_avail) {
			generic_summary.marshal (writer, "generic_summary");
			writer.marshalJsonString ("generic_json", generic_json);
			writer.marshalBoolean    ("generic_pdl" , generic_pdl );
		}

		writer.marshalBoolean ("seq_spec_result_avail", seq_spec_result_avail);
		if (seq_spec_result_avail) {
			seq_spec_summary.marshal (writer, "seq_spec_summary");
			writer.marshalJsonString ("seq_spec_json", seq_spec_json);
			writer.marshalBoolean    ("seq_spec_pdl" , seq_spec_pdl );
		}

		writer.marshalBoolean ("bayesian_result_avail", bayesian_result_avail);
		if (bayesian_result_avail) {
			bayesian_summary.marshal (writer, "bayesian_summary");
			writer.marshalJsonString ("bayesian_json", bayesian_json);
			writer.marshalBoolean    ("bayesian_pdl" , bayesian_pdl );
		}
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		result_time     = reader.unmarshalLong   ("result_time"    );
		advisory_lag    = reader.unmarshalLong   ("advisory_lag"   );
		injectable_text = reader.unmarshalString ("injectable_text");

		catalog_result_avail = reader.unmarshalBoolean ("catalog_result_avail");
		if (catalog_result_avail) {
			catalog_start_time   = reader.unmarshalLong   ("catalog_start_time"  );
			catalog_end_time     = reader.unmarshalLong   ("catalog_end_time"    );
			catalog_eqk_count    = reader.unmarshalInt    ("catalog_eqk_count"   );
			catalog_max_mag      = reader.unmarshalDouble ("catalog_max_mag"     );
			catalog_max_event_id = reader.unmarshalString ("catalog_max_event_id");
			catalog_aftershocks = null;
			catalog_comcat_aftershocks = null;
		} else {
			set_default_catalog_results();
		}

		generic_result_avail = reader.unmarshalBoolean ("generic_result_avail");
		if (generic_result_avail) {
			generic_summary = (new RJ_Summary_Generic()).unmarshal (reader, "generic_summary");
			generic_json    = reader.unmarshalJsonString ("generic_json");
			generic_pdl     = reader.unmarshalBoolean    ("generic_pdl" );
			generic_model   = null;
		} else {
			set_default_generic_results();
		}

		seq_spec_result_avail = reader.unmarshalBoolean ("seq_spec_result_avail");
		if (seq_spec_result_avail) {
			seq_spec_summary = (new RJ_Summary_SequenceSpecific()).unmarshal (reader, "seq_spec_summary");
			seq_spec_json    = reader.unmarshalJsonString ("seq_spec_json");
			seq_spec_pdl     = reader.unmarshalBoolean    ("seq_spec_pdl" );
			seq_spec_model   = null;
		} else {
			set_default_seq_spec_results();
		}

		bayesian_result_avail = reader.unmarshalBoolean ("bayesian_result_avail");
		if (bayesian_result_avail) {
			bayesian_summary = (new RJ_Summary_Bayesian()).unmarshal (reader, "bayesian_summary");
			bayesian_json    = reader.unmarshalJsonString ("bayesian_json");
			bayesian_pdl     = reader.unmarshalBoolean    ("bayesian_pdl" );
			bayesian_model   = null;
		} else {
			set_default_bayesian_results();
		}

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

	public ForecastResults unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, ForecastResults obj) {

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

	public static ForecastResults unmarshal_poly (MarshalReader reader, String name) {
		ForecastResults result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("ForecastResults.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_FCAST_RESULT:
			result = new ForecastResults();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}




	//----- Testing -----

	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("ForecastResults : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1  event_id
		// Get parameters for the event, and display them.
		// Then get results for the event, and display them.

		if (args[0].equalsIgnoreCase ("test1")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastResults : Invalid 'test1' subcommand");
				return;
			}

			String the_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_event_id);

			System.out.println ("");
			System.out.println (fcmain.toString());

			// Set the forecast time to be 7 days after the mainshock

			long the_forecast_lag = Math.round(ComcatAccessor.day_millis * 7.0);

			// Get parameters

			ForecastParameters params = new ForecastParameters();
			params.fetch_all_params (the_forecast_lag, fcmain, null);

			// Display them

			System.out.println ("");
			System.out.println (params.toString());

			// Get results

			ForecastResults results = new ForecastResults();
			results.calc_all (fcmain.mainshock_time + the_forecast_lag, ADVISORY_LAG_WEEK, "test1 injectable.", fcmain, params, true);

			// Display them

			System.out.println ("");
			System.out.println (results.toString());

			return;
		}

		// Subcommand : Test #2
		// Command format:
		//  test2  event_id
		// Get parameters for the event, and display them.
		// Then get results for the event, and display them.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled results.
		// Then rebuild transient data, and display the results.

		if (args[0].equalsIgnoreCase ("test2")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastResults : Invalid 'test2' subcommand");
				return;
			}

			String the_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_event_id);

			System.out.println ("");
			System.out.println (fcmain.toString());

			// Set the forecast time to be 7 days after the mainshock

			long the_forecast_lag = Math.round(ComcatAccessor.day_millis * 7.0);

			// Get parameters

			ForecastParameters params = new ForecastParameters();
			params.fetch_all_params (the_forecast_lag, fcmain, null);

			// Display them

			System.out.println ("");
			System.out.println (params.toString());

			// Get results

			ForecastResults results = new ForecastResults();
			results.calc_all (fcmain.mainshock_time + the_forecast_lag, ADVISORY_LAG_WEEK, "", fcmain, params, true);

			// Display them

			System.out.println ("");
			System.out.println (results.toString());

			// Save catalog

			CompactEqkRupList saved_catalog_aftershocks = results.catalog_aftershocks;

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ForecastResults.marshal_poly (store, null, results);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			results = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			results = ForecastResults.unmarshal_poly (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (results.toString());

			// Rebuild transient data

			results.rebuild_all (fcmain, params, saved_catalog_aftershocks);

			System.out.println ("");
			System.out.println (results.toString());

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("ForecastResults : Unrecognized subcommand : " + args[0]);
		return;

	}

}
