package scratch.aftershockStatistics;

import java.io.IOException;
import java.net.URL;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;

import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;

import scratch.aftershockStatistics.OAFTectonicRegime;
import scratch.aftershockStatistics.OAFRegion;

public class GenericRJ_ParametersFetch {

	// dataMap - Maps tectonic domains to sets of R-J parameters.
	
	private static Map<OAFTectonicRegime, GenericRJ_Parameters> dataMap = null;

	// region_list - List of regions that override the Garcia regions.

	private static List<OAFRegion> region_list = null;

	// f_loaded - True if data has been successfully loaded.

	private static boolean f_loaded = false;

	// List of Garcia regions.

	private static final String[] garcia_regions = {
		"ANSR-DEEPCON",
		"ANSR-HOTSPOT",
		"ANSR-OCEANBD",
		"ANSR-SHALCON",
		"ANSR-ABSLDEC",
		"ANSR-ABSLOCB",
		"ANSR-ABSLSHC",
		"SCR-ABVSLAB",
		"SCR-GENERIC",
		"SOR-ABVSLAB",
		"SOR-GENERIC",
		"SZ-GENERIC",
		"SZ-INLBACK",
		"SZ-ONSHORE",
		"SZ-OUTERTR"
	};

	// load_table_int - Load an integer value for the tables.
	// An exception is thrown if there is no integer in the stream,
	// or if it lies outside the given minimum and maximum values.

	private static int load_table_int (Scanner sc, int minval, int maxval) {
	
		// Get the integer value

		int value;

		try {
			value = sc.nextInt();
		} catch (InputMismatchException e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Badly formatted integer value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Integer value out-of-range");
		}

		return value;
	}

	private static int load_table_int (Scanner sc, int minval) {
		return load_table_int (sc, minval, Integer.MAX_VALUE);
	}

	private static int load_table_int (Scanner sc) {
		return load_table_int (sc, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}

	// load_table_double - Load a double value for the tables.
	// An exception is thrown if there is no double in the stream,
	// or if it lies outside the given minimum and maximum values.

	private static double load_table_double (Scanner sc, double minval, double maxval) {
	
		// Get the double value

		double value;

		try {
			value = sc.nextDouble();
		} catch (InputMismatchException e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Badly formatted double value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Double value out-of-range");
		}

		return value;
	}

	private static double load_table_double (Scanner sc, double minval) {
		return load_table_double (sc, minval, Double.MAX_VALUE);
	}

	private static double load_table_double (Scanner sc) {
		return load_table_double (sc, -Double.MAX_VALUE, Double.MAX_VALUE);
	}

	// load_table_regime - Load a tectonic regime for the tables.

	private static OAFTectonicRegime load_table_regime (Scanner sc) {
	
		// Get a string with the name of the regime

		String regime_name;

		try {
			regime_name = sc.next();
		} catch (NoSuchElementException e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Unexpected end-of-file", e);
		}

		return OAFTectonicRegime.forName (regime_name);
	}

	// load_table_location - Load a location for the tables.

	private static Location load_table_location (Scanner sc) {
		double lat = load_table_double (sc, -90.0, 90.0);
		double lon = load_table_double (sc, -180.0, 180.0);
		return new Location(lat, lon);
	}

	// load_table_border_type - Load a border type for the tables.

	private static BorderType load_table_border_type (Scanner sc) {
		int value = load_table_int (sc, 1, 2);
		switch (value) {
		case 1:
			return BorderType.MERCATOR_LINEAR;
		case 2:
			return BorderType.GREAT_CIRCLE;
		}
		return null;
	}

	// load_data - Load parameters from the data file.
	// The data file format is:
	//	[int]		Number of tectonic regimes
	//	[repeated]	Repeated once for each tectonic regime:
	//		[string]	Name of tectonic regime
	//		[double]	p-value
	//		[double]	a-value mean
	//		[double]	a-value sigma
	//		[double]	a-value sigma1
	//		[double]	a-value sigma0
	//		[double]	b-value
	//		[double]	c-value
	//	[int]		Number of special regions
	//	[repeated]	Repeated once for each special region:
	//		[string]	Name of tectonic regime to apply in this region
	//		[double]	Minimum depth in km, positive down; use -1.0e10 if no bound
	//		[double]	Maximum depth in km, positive down; use 1.0e10 if no bound
	//		[int]		Polygon border type: 1 = Mercator Linear, 2 = Great Circle
	//		[int]		Number of polygon vertices
	//		[repeated]	Repeated once for each polygon vertex:
	//			[double]	Latitude in degrees, between -90 and 90
	//			[double]	Longitude in degrees, between -180 and 180
	// Notes:
	// The list of tectonic regimes must include at least the 15 Garcia regions.
	// A tectonic regime may appear in more than one special region.
	// In each special region, the first and last polygon vertices should be the same
	// (although, according to Region, this is not strictly necessary).
	// A special region cannot cross the date line.  If this is needed, use two regions.
	// If special regions overlap, the region listed first "wins".

	private static synchronized void load_data () {

		// If data is already loaded, do nothing

		if (f_loaded) {
			return;
		}
	
		// Working data

		Map<OAFTectonicRegime, GenericRJ_Parameters> wk_dataMap;
		List<OAFRegion> wk_region_list;

		// Any exception means load has failed

		try {

			// Open the data file

			InputStream stream = GenericRJ_ParametersFetch.class.getResourceAsStream ("GenericRJ_ParametersFetch.txt");
			if (stream == null) {
				throw new RuntimeException("GenericRJ_ParametersFetch: Cannot find data file GenericRJ_ParametersFetch.txt");
			}

			Scanner sc = new Scanner (stream);

			// Make working data

			wk_dataMap = new HashMap<OAFTectonicRegime, GenericRJ_Parameters>();
			wk_region_list = new ArrayList<OAFRegion>();

			// Number of tectonic regimes, must be at least 15 for the Garcia regions

			int regime_count = load_table_int (sc, 15);

			// For each tectonic regime ...

			for (int regime_i = 0; regime_i < regime_count; ++regime_i) {

				// Get the tectonic regime

				OAFTectonicRegime regime = load_table_regime (sc);

				// Get the R&J parameters
				
				double pValue = load_table_double (sc);
				double aValue_mean = load_table_double (sc);
				double aValue_sigma = load_table_double (sc);
				double aValue_sigma1 = load_table_double (sc);
				double aValue_sigma0 = load_table_double (sc);
				double bValue = load_table_double (sc);
				double cValue = load_table_double (sc);

				// Add to our table
				
				wk_dataMap.put(regime, new GenericRJ_Parameters(
						aValue_mean, aValue_sigma, aValue_sigma0, aValue_sigma1, bValue, pValue, cValue));
			}

			// Check that we have all the Garcia regions

			for (String garcia_region : garcia_regions) {
				if (!( wk_dataMap.containsKey (OAFTectonicRegime.forName (garcia_region)) ))
				{
					throw new RuntimeException("GenericRJ_ParametersFetch: No parameters defined for Garcia region : " + garcia_region);
				}
			}

			// Number of special regions

			int region_count = load_table_int (sc, 0);

			// For each special region ...

			for (int region_i = 0; region_i < region_count; ++region_i) {

				// Get the tectonic regime

				OAFTectonicRegime regime = load_table_regime (sc);

				if (!( wk_dataMap.containsKey (regime) ))
				{
					throw new RuntimeException("GenericRJ_ParametersFetch: No parameters defined for special region : " + regime);
				}

				// Get the depth range
				
				double min_depth = load_table_double (sc);
				double max_depth = load_table_double (sc);

				if (min_depth >= max_depth) {
					throw new RuntimeException("GenericRJ_ParametersFetch: Minimum and maximum depths are reversed");
				}

				// Get the polygon border type

				BorderType border_type = load_table_border_type (sc);

				// Get the number of polygon vertices

				int vertex_count = load_table_int (sc, 3);

				// Add every polygon vertex to a list

				LocationList locs = new LocationList();

				for (int vertex_i = 0; vertex_i < vertex_count; ++vertex_i) {
					locs.add (load_table_location (sc));
				}

				// Form the region

				OAFRegion region = new OAFRegion (regime, new Region(locs, border_type), min_depth, max_depth);
				
				// Add the region to the list

				wk_region_list.add (region);
			}

			// Close the file

			sc.close();

		} catch (Exception e) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Unable to load data file GenericRJ_ParametersFetch.txt", e);
		}

		// Save our working data into the static variables

		dataMap = wk_dataMap;
		region_list = wk_region_list;

		// Set flag indicating data is loaded

		f_loaded = true;
		return;
	}

	// Constructor loads the data if needed.
	
	public GenericRJ_ParametersFetch() {
		load_data();
	}

//	// Build the region that describes California.
//	
//	private static Region buildANSS_CA_Region() {
//		LocationList locs = new LocationList();
//		// region.scsn locations
//		locs.add(new Location(36.6847,   -117.793));
//		locs.add(new Location(35.8000,   -116.400));
//		locs.add(new Location(34.0815,   -114.472));
//		locs.add(new Location(32.0000,   -114.333));
//		locs.add(new Location(32.0000,   -120.500));
//		locs.add(new Location(34.6945,   -121.380));
//		// locs.add(new Location(36.6847,   -117.793)); // this was to close it, instead skip and continue on to NCSN
//		// region.ncsn locations, in reverse to match the ordering and continue on in SCSN
//		//locs.add(new Location(34.6945,   -121.380)); // ignore, duplicate with  SCSN
//		locs.add(new Location(40.0000,   -125.500));
//		locs.add(new Location(43.0200,   -125.000));
//		locs.add(new Location(42.0000,   -122.700));
//		locs.add(new Location(42.0000,   -121.417));
//		locs.add(new Location(39.5000,   -120.750));
//		locs.add(new Location(37.7500,   -119.500));
//		locs.add(new Location(37.7500,   -118.250));
//		locs.add(new Location(36.6847,   -117.793)); // this is also the first one from SCSN, and closes the region
//		//locs.add(new Location(34.6945,   -121.380)); // this was to close it, skip it (and we're already closed with SCSN)
//		
//		return new Region(locs, null);
//	}

// Comments for California parameters, retained here for historical reasons.
// via e-mail from Jeanne 2/8/17 and 2/9/17, subject "OAF To-Do List; Jan 31, 2017":
/*
	* I dug up the relevant email thread on the California regional parameters.  The summary is to use the modified
	* R&J parameters a=-1.85, b = 0.91, p = 1.08, and c = 0.05.  R&J don't define a-sigma, but I think we were okay
	* with using the global values of sigma0=0.49 and sigma1=750.  Use the same equation for the completeness, except with
	* Mcat= 2.5.  The California spatial region is given in the attached files that Andy provided (the union of
	* region.ncsn and region.scsn).
	*/
// aValue_sigma value of 1.76 is from the second e-mail, "Integrating equation 8 over M5 to M8, weighted by a MFD
// with b=1, the sigma value is 1.76. I know this is kind of large, but this is the extrapolation to lower
// magnitude of the observed relation between magnitude and sigma."
	
	/**
	 * Fetches a Garcia region for the given location and returns default omori parameters
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericRJ_Parameters get(Location loc) {
		OAFTectonicRegime region = getRegion(loc);
		return get(region);
	}
	
	/**
	 * Return parameters for a tectonic regime, throw exception if regime is unknown.
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericRJ_Parameters get(OAFTectonicRegime region) {
		GenericRJ_Parameters params = dataMap.get(region);
		if (params == null) {
			throw new RuntimeException("GenericRJ_ParametersFetch: Unknown tectonic regime : " + region);
		}
		return params;
	}
	
	/**
	 * Return parameters for a tectonic regime, or null if tectonic regime is unknown.
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericRJ_Parameters getOrNull(OAFTectonicRegime region) {
		GenericRJ_Parameters params = dataMap.get(region);
		return params;
	}
	
	public OAFTectonicRegime getRegion(Location loc) {

		// Location allows longitude -180 to 360, so bring it in range

		double lat = loc.getLatitude();
		double lon = loc.getLongitude();
		boolean f_changed = false;

		while (lon > 180.0) {
			lon -= 360.0;
			f_changed = true;
		}

		// For locations very close to the date line or pole, nudge them so that
		// they compare as expected to regions that end right at the date line or pole

		if (lon < -179.995) {
			lon = -179.995;
			f_changed = true;
		}

		if (lon > 179.995) {
			lon = 179.995;
			f_changed = true;
		}

		if (lat < -89.995) {
			lat = -89.995;
			f_changed = true;
		}

		if (lat > 89.995) {
			lat = 89.995;
			f_changed = true;
		}

		if (f_changed) {
			loc = new Location(lat, lon, loc.getDepth());
		}

		// I don't know if Region.contains is thread-safe, so synchronize this

		synchronized (GenericRJ_ParametersFetch.class) {

			// If the point is in a special region, use the tectonic regime for that region

			for (OAFRegion region : region_list) {
				if (region.contains (loc)) {
					return region.get_regime();
				}
			}

		}

		// Otherwise, use the tectonic regime table to get the Garcia region

		TectonicRegimeTable regime_table = new TectonicRegimeTable();
		return OAFTectonicRegime.forName (regime_table.get_strec_name (loc.getLatitude(), loc.getLongitude()));
	}
	
	public static void main(String[] args) {
		LocationList locs = new LocationList();

		locs.add(new Location(56.04640, -149.0728, 25.00000));
		locs.add(new Location(56.04640, 210.92720, 25.00000));
		locs.add(new Location(56.04640, 210.92720, 20.00000));
		locs.add(new Location(56.04640, 210.92720, 15.00000));
		locs.add(new Location(56.04640, 210.92720, 10.00000));
		locs.add(new Location(56.04640, 210.92720, 5.00000));
		locs.add(new Location(56.04640, 210.92720, 0.00000));

		locs.add(new Location(28.2305, 84.7314, 8.22));
		locs.add(new Location(35, -118, 7d));			// CALIFORNIA
		locs.add(new Location(35, -50, 7d));

//		for (int i=0; i<5; i++) {
//			double lat = 180d*Math.random()-90d;
//			double lon = 360d*Math.random()-180d;
//			double depth = 20d*Math.random();
//			locs.add(new Location(lat, lon, depth));
//		}

		locs.add(new Location(38.0, 80.0, 15.00000));	// ANSR_ABSLDEC
		locs.add(new Location(-15.0, 175.0, 15.00000));	// ANSR_ABSLOCB
		locs.add(new Location(60.0, -175.0, 15.00000));	// ANSR_ABSLSHC
		locs.add(new Location(32.0, 100.0, 15.00000));	// ANSR_DEEPCON
		locs.add(new Location(-45.0, 70.0, 15.00000));	// ANSR_HOTSPOT
		locs.add(new Location(-65.0, 175.0, 15.00000));	// ANSR_OCEANBD
		locs.add(new Location(68.0, -175.0, 15.00000));	// ANSR_SHALCON
		locs.add(new Location(0.0, -70.0, 15.00000));	// SCR_ABVSLAB
		locs.add(new Location(40.0, -90.0, 15.00000));	// SCR_GENERIC
		locs.add(new Location(-2.0, 158.0, 15.00000));	// SOR_ABVSLAB
		locs.add(new Location(35.0, -140.0, 15.00000));	// SOR_GENERIC
		locs.add(new Location(-6.0, 156.00, 15.00000));	// SZ_GENERIC
		locs.add(new Location(18.0, 145.0, 15.00000));	// SZ_INLBACK
		locs.add(new Location(45.0, -123.0, 15.00000));	// SZ_ONSHORE
		locs.add(new Location(50.0, -150.0, 15.00000));	// SZ_OUTERTR
		
		GenericRJ_ParametersFetch fetch = new GenericRJ_ParametersFetch();
		for (Location loc : locs) {
			OAFTectonicRegime regime = fetch.getRegion(loc);
			System.out.println(loc+", "+regime+": "+fetch.get(regime));
		}
	}

}
