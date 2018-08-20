package scratch.aftershockStatistics.aafs;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import scratch.aftershockStatistics.util.EventNotFoundException;
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
import scratch.aftershockStatistics.GenericRJ_Parameters;
import scratch.aftershockStatistics.GenericRJ_ParametersFetch;
import scratch.aftershockStatistics.MagCompPage_Parameters;
import scratch.aftershockStatistics.MagCompPage_ParametersFetch;
import scratch.aftershockStatistics.OAFTectonicRegime;
import scratch.aftershockStatistics.SeqSpecRJ_Parameters;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;

/**
 * Parameters for constructing a forecast.
 * Author: Michael Barall 04/21/2018.
 *
 * All fields are public, since there is little benefit to having lots of getters and setters.
 */
public class ForecastParameters {

	//----- Constants -----

	// Parameter fetch methods.

	public static final int FETCH_METH_MIN = 0;			// Minimum value
	public static final int FETCH_METH_AUTO = 0;		// Parameters determined by automatic system (default)
	public static final int FETCH_METH_ANALYST = 1;		// Parameters selected by analyst
	public static final int FETCH_METH_SUPPRESS = 2;	// Do not fetch parameters
	public static final int FETCH_METH_MAX = 2;			// Maximum value

	// Result calculation methods.

	public static final int CALC_METH_MIN = 0;			// Minimum value
	public static final int CALC_METH_AUTO_PDL = 0;		// Calculate by automatic system, eligible for PDL (default)
	public static final int CALC_METH_AUTO_NO_PDL = 1;	// Calculate by automatic system, not eligible for PDL
	public static final int CALC_METH_SUPPRESS = 2;		// Do not calculate result
	public static final int CALC_METH_MAX = 2;			// Maximum value

	// Special values of injectable text.

	public static final String INJ_TXT_USE_DEFAULT = "use-default";	// Use the configured default


	//----- Root parameters -----

	// Lag time of forecast, in milliseconds since the mainshock.
	// Restriction: Must be greater than the value in the previous forecast.

	public long forecast_lag = 0L;


	//----- Control parameters -----

	// Calculation method for generic result.

	public int generic_calc_meth = CALC_METH_AUTO_PDL;

	// Calculation method for sequence specific result.

	public int seq_spec_calc_meth = CALC_METH_AUTO_PDL;

	// Calculation method for bayesian result.

	public int bayesian_calc_meth = CALC_METH_AUTO_PDL;

	// Injectable text for PDL JSON files, or "" for none, or INJ_TXT_USE_DEFAULT for configured default.

	public String injectable_text = INJ_TXT_USE_DEFAULT;

	// Set control parameters to default.

	public void set_default_control_params () {
		generic_calc_meth = CALC_METH_AUTO_PDL;
		seq_spec_calc_meth = CALC_METH_AUTO_PDL;
		bayesian_calc_meth = CALC_METH_AUTO_PDL;
		injectable_text = INJ_TXT_USE_DEFAULT;
		return;
	}

	// Fetch control parameters.

	public void fetch_control_params (ForecastMainshock fcmain, ForecastParameters prior_params) {

		// If there are prior parameters, copy them

		if (prior_params != null) {
			generic_calc_meth  = prior_params.generic_calc_meth;
			seq_spec_calc_meth = prior_params.seq_spec_calc_meth;
			bayesian_calc_meth = prior_params.bayesian_calc_meth;
			injectable_text    = prior_params.injectable_text;
			return;
		}

		// Use defaults
	
		set_default_control_params();
		return;
	}

	// Get the effective injectable text.
	// Note: The return value is always non-null, and is "" if no injectable text is desired.

	public String get_eff_injectable_text (String def_injectable_text) {
		String result = injectable_text;
		if (result.equals (INJ_TXT_USE_DEFAULT)) {
			result = def_injectable_text;
			if (result == null) {
				result = "";
			}
		}
		return result;
	}


	//----- R&J generic parameters -----

	// Generic parameter fetch method.

	public int generic_fetch_meth = FETCH_METH_AUTO;

	// Generic parameter available flag.

	public boolean generic_avail = false;

	// Tectonic regime (null iff omitted).
	// (For analyst-supplied values, cannot be null, but can be an empty string,
	// and need not be the name of a known tectonic regime.)

	public String generic_regime = null;

	// Generic parameters (null iff omitted).

	public GenericRJ_Parameters generic_params = null;

	// Set generic parameters to default.

	public void set_default_generic_params () {
		generic_regime = null;
		generic_params = null;
		return;
	}

	// Fetch generic parameters.
	// Note: Mainshock parameters must be fetched first.

	public void fetch_generic_params (ForecastMainshock fcmain, ForecastParameters prior_params) {

		// Inherit fetch method from prior parameters, or use default

		if (prior_params != null) {
			generic_fetch_meth = prior_params.generic_fetch_meth;
		} else {
			generic_fetch_meth = FETCH_METH_AUTO;
		}

		// Handle non-auto fetch methods

		switch (generic_fetch_meth) {

		// Analyst, copy from prior parameters

		case FETCH_METH_ANALYST:
			generic_avail = prior_params.generic_avail;
			generic_regime = prior_params.generic_regime;
			generic_params = prior_params.generic_params;
			return;

		// Suppress, make not available

		case FETCH_METH_SUPPRESS:
			generic_avail = false;
			set_default_generic_params();
			return;
		}

		// If we don't have mainshock parameters, then we can't fetch generic parameters

		if (!( fcmain.mainshock_avail )) {
			generic_avail = false;
			set_default_generic_params();
			return;
		}

		// Fetch parameters based on mainshock location
		
		GenericRJ_ParametersFetch fetch = new GenericRJ_ParametersFetch();
		OAFTectonicRegime regime = fetch.getRegion (fcmain.get_eqk_location());

		generic_avail = true;
		generic_regime = regime.toString();
		generic_params = fetch.get(regime);

		return;
	}


	//----- R&J magnitude of completeness parameters -----

	// Magnitude of completeness parameter fetch method.

	public int mag_comp_fetch_meth = FETCH_METH_AUTO;

	// Magnitude of completeness parameter available flag.

	public boolean mag_comp_avail = false;

	// Tectonic regime (null iff omitted).
	// (For analyst-supplied values, cannot be null, but can be an empty string,
	// and need not be the name of a known tectonic regime.)

	public String mag_comp_regime = null;

	// Magnitude of completeness parameters (null iff omitted).

	public MagCompPage_Parameters mag_comp_params = null;

	// Set magnitude of completeness parameters to default.

	public void set_default_mag_comp_params () {
		mag_comp_regime = null;
		mag_comp_params = null;
		return;
	}

	// Fetch magnitude of completeness parameters.
	// Note: Mainshock parameters must be fetched first.

	public void fetch_mag_comp_params (ForecastMainshock fcmain, ForecastParameters prior_params) {

		// Inherit fetch method from prior parameters, or use default

		if (prior_params != null) {
			mag_comp_fetch_meth = prior_params.mag_comp_fetch_meth;
		} else {
			mag_comp_fetch_meth = FETCH_METH_AUTO;
		}

		// Handle non-auto fetch methods

		switch (mag_comp_fetch_meth) {

		// Analyst, copy from prior parameters

		case FETCH_METH_ANALYST:
			mag_comp_avail = prior_params.mag_comp_avail;
			mag_comp_regime = prior_params.mag_comp_regime;
			mag_comp_params = prior_params.mag_comp_params;
			return;

		// Suppress, make not available

		case FETCH_METH_SUPPRESS:
			mag_comp_avail = false;
			set_default_mag_comp_params();
			return;
		}

		// If we don't have mainshock parameters, then we can't fetch magnitude of completeness parameters

		if (!( fcmain.mainshock_avail )) {
			mag_comp_avail = false;
			set_default_mag_comp_params();
			return;
		}

		// Fetch parameters based on mainshock location
		
		MagCompPage_ParametersFetch fetch = new MagCompPage_ParametersFetch();
		OAFTectonicRegime regime = fetch.getRegion (fcmain.get_eqk_location());

		mag_comp_avail = true;
		mag_comp_regime = regime.toString();
		mag_comp_params = fetch.get(regime);

		return;
	}


	//----- R&J sequence specific parameters -----

	// Sequence specific parameter fetch method.

	public int seq_spec_fetch_meth = FETCH_METH_AUTO;

	// Sequence specific parameter available flag.

	public boolean seq_spec_avail = false;

	// Sequence specific parameters (null iff omitted).

	public SeqSpecRJ_Parameters seq_spec_params = null;

	// Set sequence specific parameters to default.

	public void set_default_seq_spec_params () {
		seq_spec_params = null;
		return;
	}

	// Fetch sequence specific parameters.
	// Note: Generic parameters must be fetched first.

	public void fetch_seq_spec_params (ForecastMainshock fcmain, ForecastParameters prior_params) {

		// Inherit fetch method from prior parameters, or use default

		if (prior_params != null) {
			seq_spec_fetch_meth = prior_params.seq_spec_fetch_meth;
		} else {
			seq_spec_fetch_meth = FETCH_METH_AUTO;
		}

		// Handle non-auto fetch methods

		switch (seq_spec_fetch_meth) {

		// Analyst, copy from prior parameters

		case FETCH_METH_ANALYST:
			seq_spec_avail = prior_params.seq_spec_avail;
			seq_spec_params = prior_params.seq_spec_params;
			return;

		// Suppress, make not available

		case FETCH_METH_SUPPRESS:
			seq_spec_avail = false;
			set_default_seq_spec_params();
			return;
		}

		// If we don't have generic parameters, then we can't fetch sequence specific parameters

		if (!( generic_avail )) {
			seq_spec_avail = false;
			set_default_seq_spec_params();
			return;
		}

		// Fetch parameters based on generic parameters

		seq_spec_avail = true;
		seq_spec_params = new SeqSpecRJ_Parameters(generic_params);

		return;
	}


	//----- Aftershock search parameters -----

	// Aftershock search parameter fetch method.

	public int aftershock_search_fetch_meth = FETCH_METH_AUTO;

	// Aftershock search parameter available flag.

	public boolean aftershock_search_avail = false;

	// Aftershock search region (null iff omitted).

	public SphRegion aftershock_search_region = null;

	// Minimum search time, in days after the mainshock.

	public double min_days = 0.0;

	// Maximum search time, in days after the mainshock.
	// Note: This get set to the current time, even if there are analyst-supplied values.

	public double max_days = 0.0;

	// Minimum search depth, in kilometers.

	public double min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;

	// Maximum search depth, in kilometers.
	// (Comcat has 1000 km maximum, OpenSHA has 700 km maximum.)

	public double max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;

	// Minimum magnitude to consider in search, or -10.0 if none.

	public double min_mag = -10.0;

	// Set aftershock search parameters to default.

	public void set_default_aftershock_search_params () {
		aftershock_search_region = null;
		min_days = 0.0;
		max_days = 0.0;
		min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;
		max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;
		min_mag = -10.0;
		return;
	}

	// Fetch aftershock search region.
	// Note: Mainshock parameters must be fetched first.

	public void fetch_aftershock_search_region (ForecastMainshock fcmain, ForecastParameters prior_params) {

		// Inherit fetch method from prior parameters, or use default

		if (prior_params != null) {
			aftershock_search_fetch_meth = prior_params.aftershock_search_fetch_meth;
		} else {
			aftershock_search_fetch_meth = FETCH_METH_AUTO;
		}

		// Handle non-auto fetch methods

		switch (aftershock_search_fetch_meth) {

		// Analyst, copy from prior parameters

		case FETCH_METH_ANALYST:
			aftershock_search_avail = prior_params.aftershock_search_avail;
			aftershock_search_region = prior_params.aftershock_search_region;
			min_days = prior_params.min_days;
			max_days = prior_params.max_days;
			min_depth = prior_params.min_depth;
			max_depth = prior_params.max_depth;
			min_mag = prior_params.min_mag;

			// Special handling for max_days, try to make it match forecast_lag if available

			max_days = ((double)forecast_lag)/ComcatAccessor.day_millis;

			return;

		// Suppress, make not available

		case FETCH_METH_SUPPRESS:
			aftershock_search_avail = false;
			set_default_aftershock_search_params();
			return;
		}

		// If we don't have mainshock and magnitude of completeness parameters, then we can't fetch aftershock search parameters

		if (!( fcmain.mainshock_avail && mag_comp_avail )) {
			aftershock_search_avail = false;
			set_default_aftershock_search_params();
			return;
		}

		// Get minimum magnitude and radius multiplier parameters

		double sample_min_mag = mag_comp_params.get_magSample();
		double sample_radius_mult = mag_comp_params.get_radiusSample();

		double centroid_min_mag = mag_comp_params.get_magCentroid();
		double centroid_radius_mult = mag_comp_params.get_radiusCentroid();

		// Get initial search radius from Wells and Coppersmith

		WC1994_MagLengthRelationship wcMagLen = new WC1994_MagLengthRelationship();
		double radius = wcMagLen.getMedianLength(fcmain.mainshock_mag);

		// The initial region is a circle centered at the epicenter

		SphRegion initial_region = SphRegion.makeCircle (fcmain.get_sph_eqk_location(), radius * centroid_radius_mult);

		// Time range used for sampling aftershocks, in days since the mainshock

		min_days = 0.0;

		//max_days = ((double)System.currentTimeMillis())/ComcatAccessor.day_millis;
		max_days = ((double)forecast_lag)/ComcatAccessor.day_millis;

		// Depth range used for sampling aftershocks, in kilometers

		min_depth = ComcatAccessor.DEFAULT_MIN_DEPTH;
		max_depth = ComcatAccessor.DEFAULT_MAX_DEPTH;

		// Minimum magnitude used for sampling aftershocks

		min_mag = sample_min_mag;

		// Retrieve list of aftershocks in the initial region

		ObsEqkRupList aftershocks;

		try {
			ComcatAccessor accessor = new ComcatAccessor();
			aftershocks = accessor.fetchAftershocks(fcmain.get_eqk_rupture(), min_days, max_days, min_depth, max_depth, initial_region, initial_region.getPlotWrap(), centroid_min_mag);
		} catch (Exception e) {
			throw new RuntimeException("ForecastParameters.fetch_aftershock_search_region: Comcat exception", e);
		}

		// Center of search region

		Location centroid;

		// If no aftershocks, use the hypocenter location

		if (aftershocks.isEmpty()) {
			centroid = fcmain.get_eqk_location();
		}

		// Otherwise, use the centroid of the aftershocks

		else {
			centroid = AftershockStatsCalc.getSphCentroid(fcmain.get_eqk_rupture(), aftershocks);
		}

		// Search region is a circle centered at the centroid (or hypocenter if no aftershocks)
			
		aftershock_search_region = SphRegion.makeCircle (new SphLatLon(centroid), radius * sample_radius_mult);

		aftershock_search_avail = true;
		return;
	}


	//----- Transient parameters -----

//	// The configured default injectable text, or "" if none, or null if not set.
//	// Note: This parameter is not marshaled/unmarshaled.
//
//	public String def_injectable_text = null;

	// Set transient parameters to default.

	public void set_default_transient_params () {
//		def_injectable_text = null;
		return;
	}

//	// Get the effective injectable text.
//	// Note: The return value is always non-null, and is "" if no injectable text is desired.
//
//	public String get_eff_injectable_text () {
//		String result = injectable_text;
//		if (result.equals (INJ_TXT_USE_DEFAULT)) {
//			result = def_injectable_text;
//			if (result == null) {
//				result = "";
//			}
//		}
//		return result;
//	}


	//----- Construction -----

	// Default constructor.

	public ForecastParameters () {}

	// Fetch all parameters.

	public void fetch_all_params (long the_forecast_lag, ForecastMainshock fcmain, ForecastParameters prior_params) {
		forecast_lag = the_forecast_lag;
		fetch_control_params (fcmain, prior_params);
		fetch_generic_params (fcmain, prior_params);
		fetch_mag_comp_params (fcmain, prior_params);
		fetch_seq_spec_params (fcmain, prior_params);
		fetch_aftershock_search_region (fcmain, prior_params);
		return;
	}

	// Set everything to default.
	// This is a useful starting point for setting up analyst parameters.

	public void setup_all_default () {
		forecast_lag = 0L;

		set_default_control_params();

		generic_fetch_meth = FETCH_METH_AUTO;
		generic_avail = false;
		set_default_generic_params();

		mag_comp_fetch_meth = FETCH_METH_AUTO;
		mag_comp_avail = false;
		set_default_mag_comp_params();

		seq_spec_fetch_meth = FETCH_METH_AUTO;
		seq_spec_avail = false;
		set_default_seq_spec_params();

		aftershock_search_fetch_meth = FETCH_METH_AUTO;
		aftershock_search_avail = false;
		set_default_aftershock_search_params();

		set_default_transient_params();
	
		return;
	}

	// Display our contents

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();

		result.append ("ForecastParameters:" + "\n");

		result.append ("forecast_lag = " + forecast_lag + "\n");

		result.append ("generic_calc_meth = " + generic_calc_meth + "\n");
		result.append ("seq_spec_calc_meth = " + seq_spec_calc_meth + "\n");
		result.append ("bayesian_calc_meth = " + bayesian_calc_meth + "\n");
		result.append ("injectable_text = " + injectable_text + "\n");

		result.append ("generic_fetch_meth = " + generic_fetch_meth + "\n");
		result.append ("generic_avail = " + generic_avail + "\n");
		if (generic_avail) {
			result.append ("generic_regime = " + generic_regime + "\n");
			result.append ("generic_params = " + generic_params.toString() + "\n");
		}

		result.append ("mag_comp_fetch_meth = " + mag_comp_fetch_meth + "\n");
		result.append ("mag_comp_avail = " + mag_comp_avail + "\n");
		if (mag_comp_avail) {
			result.append ("mag_comp_regime = " + mag_comp_regime + "\n");
			result.append ("mag_comp_params = " + mag_comp_params.toString() + "\n");
		}

		result.append ("seq_spec_fetch_meth = " + seq_spec_fetch_meth + "\n");
		result.append ("seq_spec_avail = " + seq_spec_avail + "\n");
		if (seq_spec_avail) {
			result.append ("seq_spec_params = " + seq_spec_params.toString() + "\n");
		}

		result.append ("aftershock_search_fetch_meth = " + aftershock_search_fetch_meth + "\n");
		result.append ("aftershock_search_avail = " + aftershock_search_avail + "\n");
		if (aftershock_search_avail) {
			result.append ("aftershock_search_region = " + aftershock_search_region.toString() + "\n");
			result.append ("min_days = " + min_days + "\n");
			result.append ("max_days = " + max_days + "\n");
			result.append ("min_depth = " + min_depth + "\n");
			result.append ("max_depth = " + max_depth + "\n");
			result.append ("min_mag = " + min_mag + "\n");
		}

//		if (def_injectable_text != null) {
//			result.append ("def_injectable_text = " + def_injectable_text + "\n");
//		}

		return result.toString();
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 22001;

	private static final String M_VERSION_NAME = "ForecastParameters";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 22000;
	protected static final int MARSHAL_FCAST_PARAM = 22001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected int get_marshal_type () {
		return MARSHAL_FCAST_PARAM;
	}

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalLong   ("forecast_lag"   , forecast_lag   );

		writer.marshalInt    ("generic_calc_meth" , generic_calc_meth );
		writer.marshalInt    ("seq_spec_calc_meth", seq_spec_calc_meth);
		writer.marshalInt    ("bayesian_calc_meth", bayesian_calc_meth);
		writer.marshalString ("injectable_text"   , injectable_text   );

		writer.marshalInt     ("generic_fetch_meth", generic_fetch_meth);
		writer.marshalBoolean ("generic_avail"     , generic_avail     );
		if (generic_avail) {
			writer.marshalString ("generic_regime", generic_regime);
			generic_params.marshal (writer, "generic_params");
		}

		writer.marshalInt     ("mag_comp_fetch_meth", mag_comp_fetch_meth);
		writer.marshalBoolean ("mag_comp_avail"     , mag_comp_avail     );
		if (mag_comp_avail) {
			writer.marshalString ("mag_comp_regime", mag_comp_regime);
			mag_comp_params.marshal (writer, "mag_comp_params");
		}

		writer.marshalInt     ("seq_spec_fetch_meth", seq_spec_fetch_meth);
		writer.marshalBoolean ("seq_spec_avail"     , seq_spec_avail     );
		if (seq_spec_avail) {
			seq_spec_params.marshal (writer, "seq_spec_params");
		}

		writer.marshalInt     ("aftershock_search_fetch_meth", aftershock_search_fetch_meth);
		writer.marshalBoolean ("aftershock_search_avail"     , aftershock_search_avail     );
		if (aftershock_search_avail) {
			SphRegion.marshal_poly (writer, "aftershock_search_region", aftershock_search_region);
			writer.marshalDouble ("min_days" , min_days );
			writer.marshalDouble ("max_days" , max_days );
			writer.marshalDouble ("min_depth", min_depth);
			writer.marshalDouble ("max_depth", max_depth);
			writer.marshalDouble ("min_mag"  , min_mag  );
		}
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		forecast_lag    = reader.unmarshalLong   ("forecast_lag"   );

		generic_calc_meth  = reader.unmarshalInt    ("generic_calc_meth" , CALC_METH_MIN, CALC_METH_MAX);
		seq_spec_calc_meth = reader.unmarshalInt    ("seq_spec_calc_meth", CALC_METH_MIN, CALC_METH_MAX);
		bayesian_calc_meth = reader.unmarshalInt    ("bayesian_calc_meth", CALC_METH_MIN, CALC_METH_MAX);
		injectable_text    = reader.unmarshalString ("injectable_text");

		generic_fetch_meth = reader.unmarshalInt     ("generic_fetch_meth", FETCH_METH_MIN, FETCH_METH_MAX);
		generic_avail      = reader.unmarshalBoolean ("generic_avail");
		if (generic_avail) {
			generic_regime = reader.unmarshalString ("generic_regime");
			generic_params = (new GenericRJ_Parameters()).unmarshal (reader, "generic_params");
		} else {
			set_default_generic_params();
		}

		mag_comp_fetch_meth = reader.unmarshalInt     ("mag_comp_fetch_meth", FETCH_METH_MIN, FETCH_METH_MAX);
		mag_comp_avail      = reader.unmarshalBoolean ("mag_comp_avail");
		if (mag_comp_avail) {
			mag_comp_regime = reader.unmarshalString ("mag_comp_regime");
			mag_comp_params = (new MagCompPage_Parameters()).unmarshal (reader, "mag_comp_params");
		} else {
			set_default_mag_comp_params();
		}

		seq_spec_fetch_meth = reader.unmarshalInt     ("seq_spec_fetch_meth", FETCH_METH_MIN, FETCH_METH_MAX);
		seq_spec_avail      = reader.unmarshalBoolean ("seq_spec_avail");
		if (seq_spec_avail) {
			seq_spec_params = (new SeqSpecRJ_Parameters()).unmarshal (reader, "seq_spec_params");
		} else {
			set_default_seq_spec_params();
		}

		aftershock_search_fetch_meth = reader.unmarshalInt     ("aftershock_search_fetch_meth", FETCH_METH_MIN, FETCH_METH_MAX);
		aftershock_search_avail      = reader.unmarshalBoolean ("aftershock_search_avail");
		if (aftershock_search_avail) {
			aftershock_search_region = SphRegion.unmarshal_poly (reader, "aftershock_search_region");
			if (aftershock_search_region == null) {
				throw new MarshalException ("Aftershock search region is null");
			}
			min_days  = reader.unmarshalDouble ("min_days" );
			max_days  = reader.unmarshalDouble ("max_days" );
			min_depth = reader.unmarshalDouble ("min_depth");
			max_depth = reader.unmarshalDouble ("max_depth");
			min_mag   = reader.unmarshalDouble ("min_mag"  );
		} else {
			set_default_aftershock_search_params();
		}

		set_default_transient_params();

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

	public ForecastParameters unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, ForecastParameters obj) {

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

	public static ForecastParameters unmarshal_poly (MarshalReader reader, String name) {
		ForecastParameters result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("ForecastParameters.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_FCAST_PARAM:
			result = new ForecastParameters();
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
			System.err.println ("ForecastParameters : Missing subcommand");
			return;
		}




		// Subcommand : Test #1
		// Command format:
		//  test1  query_event_id
		// Get parameters for the event, and display them.

		if (args[0].equalsIgnoreCase ("test1")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastParameters : Invalid 'test1' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_query_event_id);

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

			return;
		}




		// Subcommand : Test #2
		// Command format:
		//  test1  query_event_id
		// Get parameters for the event, and display them.
		// Then marshal to JSON, and display the JSON.
		// Then unmarshal, and display the unmarshaled parameters.

		if (args[0].equalsIgnoreCase ("test2")) {

			// One additional argument

			if (args.length != 2) {
				System.err.println ("ForecastParameters : Invalid 'test2' subcommand");
				return;
			}

			String the_query_event_id = args[1];

			// Fetch just the mainshock info

			ForecastMainshock fcmain = new ForecastMainshock();
			fcmain.setup_mainshock_only (the_query_event_id);

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

			// Marshal to JSON

			MarshalImpJsonWriter store = new MarshalImpJsonWriter();
			ForecastParameters.marshal_poly (store, null, params);
			store.check_write_complete ();
			String json_string = store.get_json_string();

			System.out.println ("");
			System.out.println (json_string);

			// Unmarshal from JSON
			
			params = null;

			MarshalImpJsonReader retrieve = new MarshalImpJsonReader (json_string);
			params = ForecastParameters.unmarshal_poly (retrieve, null);
			retrieve.check_read_complete ();

			System.out.println ("");
			System.out.println (params.toString());

			return;
		}




		// Unrecognized subcommand.

		System.err.println ("ForecastParameters : Unrecognized subcommand : " + args[0]);
		return;

	}

}
