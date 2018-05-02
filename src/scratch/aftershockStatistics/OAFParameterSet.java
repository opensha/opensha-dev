package scratch.aftershockStatistics;

import java.io.IOException;
import java.net.URL;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Properties;
import java.util.Set;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;

import scratch.aftershockStatistics.OAFTectonicRegime;
import scratch.aftershockStatistics.OAFRegion;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;

// Class for region-dependent parameters.
// Author: Michael Barall
//
// This class holds a set of parameters, which is determined by region.
// Type T is the class that holds parameter values.
//
// The primary function of this class is:  Given a location, return the corresponding
// parameter values, in an object of type T.  The Earth is paritioned into regions
// (which may also be tectonic regimes).  A query consists of finding the region that
// contains the given location, and then returning the parameter values assigned to
// that region.
//
// The parameter values come from a file, which is read in and stored in an object
// of this class.  This class contains the logic for reading the file.  The file
// contains the parameter values to be used within each region.  A region can be one
// of the 15 Garcia regions, or a region bounded by a user-supplied polygon.  Each
// region is identified by an OAFTectonicRegime, and the user-supplied polygons are
// represented by OAFRegion objects.
//
// This is an abstract class.  There must be a subclass for each distinct type T.
// In practice, there is also an outer class, which holds a singleton of the subclass.

public abstract class OAFParameterSet<T> {

	// dataMap - Maps tectonic domains to sets of parameters.
	
	private Map<OAFTectonicRegime, T> dataMap = null;

	// region_list - List of regions that override the Garcia regions.

	private List<OAFRegion> region_list = null;

	// f_world - True if the world region is included.

	private boolean f_world = false;

	// List of Garcia regions.

	protected static final String[] garcia_regions = {
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

	// The world region.

	protected static final String world_region = "WORLD";

	// load_table_int - Load an integer value for the tables.
	// An exception is thrown if there is no integer in the stream,
	// or if it lies outside the given minimum and maximum values.

	public static int load_table_int (Scanner sc, int minval, int maxval) {
	
		// Get the integer value

		int value;

		try {
			value = sc.nextInt();
		} catch (InputMismatchException e) {
			throw new RuntimeException("OAFParameterSet: Badly formatted integer value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("OAFParameterSet: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("OAFParameterSet: Integer value out-of-range");
		}

		return value;
	}

	public static int load_table_int (Scanner sc, int minval) {
		return load_table_int (sc, minval, Integer.MAX_VALUE);
	}

	public static int load_table_int (Scanner sc) {
		return load_table_int (sc, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}

	// load_table_long - Load a long value for the tables.
	// An exception is thrown if there is no long in the stream,
	// or if it lies outside the given minimum and maximum values.

	public static long load_table_long (Scanner sc, long minval, long maxval) {
	
		// Get the long value

		long value;

		try {
			value = sc.nextLong();
		} catch (InputMismatchException e) {
			throw new RuntimeException("OAFParameterSet: Badly formatted long integer value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("OAFParameterSet: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("OAFParameterSet: Long integer value out-of-range");
		}

		return value;
	}

	public static long load_table_long (Scanner sc, long minval) {
		return load_table_long (sc, minval, Long.MAX_VALUE);
	}

	public static long load_table_long (Scanner sc) {
		return load_table_long (sc, Long.MIN_VALUE, Long.MAX_VALUE);
	}

	// load_table_double - Load a double value for the tables.
	// An exception is thrown if there is no double in the stream,
	// or if it lies outside the given minimum and maximum values.

	public static double load_table_double (Scanner sc, double minval, double maxval) {
	
		// Get the double value

		double value;

		try {
			value = sc.nextDouble();
		} catch (InputMismatchException e) {
			throw new RuntimeException("OAFParameterSet: Badly formatted double value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("OAFParameterSet: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("OAFParameterSet: Double value out-of-range");
		}

		return value;
	}

	public static double load_table_double (Scanner sc, double minval) {
		return load_table_double (sc, minval, Double.MAX_VALUE);
	}

	public static double load_table_double (Scanner sc) {
		return load_table_double (sc, -Double.MAX_VALUE, Double.MAX_VALUE);
	}

	// load_table_boolean - Load a boolean value for the tables.
	// An exception is thrown if there is no boolean in the stream.

	public static boolean load_table_boolean (Scanner sc) {
	
		// Get the boolean value

		boolean value;

		try {
			value = sc.nextBoolean();
		} catch (InputMismatchException e) {
			throw new RuntimeException("OAFParameterSet: Badly formatted boolean value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("OAFParameterSet: Unexpected end-of-file", e);
		}

		return value;
	}

	// load_table_regime - Load a tectonic regime for the tables.

	public static OAFTectonicRegime load_table_regime (Scanner sc) {
	
		// Get a string with the name of the regime

		String regime_name;

		try {
			regime_name = sc.next();
		} catch (NoSuchElementException e) {
			throw new RuntimeException("OAFParameterSet: Unexpected end-of-file", e);
		}

		return OAFTectonicRegime.forName (regime_name);
	}

	// load_table_location - Load a location for the tables.

	public static Location load_table_location (Scanner sc) {
		double lat = load_table_double (sc, -90.0, 90.0);
		double lon = load_table_double (sc, -180.0, 180.0);
		return new Location(lat, lon);
	}

	// load_table_border_type - Load a border type for the tables.

	public static BorderType load_table_border_type (Scanner sc) {
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
	//		[any]		The parameter values, see load_param()
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
	// The list of tectonic regimes must include, at a minimum, either the 15 Garcia
	// regions or a world region (but not both).
	// A tectonic regime may appear in more than one special region.
	// In each special region, the first and last polygon vertices should be the same
	// (although, according to Region, this is not strictly necessary).
	// A special region cannot cross the date line.  If this is needed, use two regions.
	// If special regions overlap, the region listed first "wins".
	//
	// The caller must supply the Scanner from which parameters can be read.
	// In case of error, this function throws RuntimeException (or another unchecked exception).

	public void load_data (Scanner sc) {

		// Nothing loaded

		dataMap = null;
		region_list = null;

		// Make working data

		Map<OAFTectonicRegime, T> wk_dataMap = new HashMap<OAFTectonicRegime, T>();
		List<OAFRegion> wk_region_list = new ArrayList<OAFRegion>();

		// Number of tectonic regimes, must be at least 1

		int regime_count = load_table_int (sc, 1);

		// For each tectonic regime ...

		for (int regime_i = 0; regime_i < regime_count; ++regime_i) {

			// Get the tectonic regime

			OAFTectonicRegime regime = load_table_regime (sc);

			// Get the parameter values and add to our table
				
			wk_dataMap.put(regime, load_parameter_values (sc));
		}

		// Check if we have a world region

		f_world = wk_dataMap.containsKey (OAFTectonicRegime.forName (world_region));

		if (f_world) {

			// If we have a world region, then there should be no Garcia regions

			for (String garcia_region : garcia_regions) {
				if ( wk_dataMap.containsKey (OAFTectonicRegime.forName (garcia_region)) )
				{
					throw new RuntimeException("OAFParameterSet: Parameters defined for both World and Garcia region : " + garcia_region);
				}
			}

		} else {

			// Check that we have all the Garcia regions

			for (String garcia_region : garcia_regions) {
				if (!( wk_dataMap.containsKey (OAFTectonicRegime.forName (garcia_region)) ))
				{
					throw new RuntimeException("OAFParameterSet: No parameters defined for Garcia region : " + garcia_region);
				}
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
				throw new RuntimeException("OAFParameterSet: No parameters defined for special region : " + regime);
			}

			// Get the depth range
				
			double min_depth = load_table_double (sc);
			double max_depth = load_table_double (sc);

			if (min_depth >= max_depth) {
				throw new RuntimeException("OAFParameterSet: Minimum and maximum depths are reversed");
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

		// Save our working data into the variables

		dataMap = wk_dataMap;
		region_list = wk_region_list;

		return;
	}

	// load_parameter_values - Load parameter values for the tables.
	// This function should create a new object of type T, read the
	// parameter values from the Scanner, and return the object.
	// In case of error, this function should throw RuntimeException.

	protected abstract T load_parameter_values (Scanner sc);

	// Open a parameter data file.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.
	// Returns a stream that is open to the requested file.
	// This function first checks if a system property named "oafcfg" is defined.
	// If oafcfg is defined, its value is the name of a disk directory,
	//  and the file is opened within that directory.
	// If oafcfg is not defined, then the file is opened in the same directory
	//  as the class file of the requesting class (in practice, in the jar file).
	// Note: If oafcfg is defined and the file is not found on disk, this routine
	//  does NOT then fall back to reading from the directory of the class file.
	//  This is to avoid inadvertently using compiled-in parameters when a
	//  parameter file is missing.
	// If the file cannot be opened, then RuntimeException is thrown.

	public static InputStream open_param_file (String filename, Class<?> requester) {

		InputStream stream = null;

		// Get the oafcfg system property

		String oafcfg = System.getProperty("oafcfg");

		if (oafcfg == null) {

			// If not found, open in the class directory

			stream = requester.getResourceAsStream(filename);
			if (stream == null) {
				throw new RuntimeException("OAFParameterSet: Cannot find data file: " + filename);
			}

		} else {
		
			// If found, open on disk

			File pathname = new File(oafcfg, filename);

			try {
				stream = new FileInputStream(pathname);
			} catch (FileNotFoundException e) {
				throw new RuntimeException("OAFParameterSet: File not found: " + pathname.toString(), e);
			} catch (Exception e) {
				throw new RuntimeException("OAFParameterSet: Failed to open file: " + pathname.toString(), e);
			}

		}

		return stream;
	}

	// load_data - Load parameters from the data file.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.
	// This function first calls open_param_file to open the file,
	// then calls load_data above to read the date.
	// In case of error, this function throws RuntimeException (or another unchecked exception).

	public void load_data (String filename, Class<?> requester) {

		// Any exception means load has failed

		try (

			// Open the data file and make the scanner, it will be auto-closed on exit from the try block

			Scanner sc = new Scanner (open_param_file (filename, requester));
		){

			// Load the data

			load_data (sc);

		} catch (Exception e) {
			throw new RuntimeException("OAFParameterSet: Unable to load data file: " + filename, e);
		}

		return;
	}

	// load_properties - Load a set of properties from the data file.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.
	//  defaults = A set of properties that supply defaults for any properties that are
	//      not in the file.  It can be null if no defaults are needed.
	//  required = A list of strings that name properties that must be in the file (or
	//      in the defaults).  An exception is thrown if any listed string is not the
	//      name of a property.  It can be null if there are no required properties.
	// This function first calls open_param_file to open the file,
	// then uses java.util.Properties to parse the file.
	// The file format is as described in the documentation for java.util.Properties.
	// In case of error, this function throws RuntimeException (or another unchecked exception).

	public static Properties load_properties (String filename, Class<?> requester, Properties defaults, String[] required) {

		Properties prop;

		// Any exception means load has failed

		try (

			// Open the data file, it will be auto-closed on exit from the try block

			InputStream stream = open_param_file (filename, requester);
		){

			// Create an empty property list

			if (defaults == null) {
				prop = new Properties();
			} else {
				prop = new Properties(defaults);
			}

			// Load the data

			prop.load(stream);

		} catch (IOException e) {
			throw new RuntimeException("OAFParameterSet: Error reading from data file: " + filename, e);
		} catch (IllegalArgumentException e) {
			throw new RuntimeException("OAFParameterSet: Malformed text in data file: " + filename, e);
		} catch (Exception e) {
			throw new RuntimeException("OAFParameterSet: Unable to load data file: " + filename, e);
		}

		// Check that required properties are supplied

		if (required != null) {
			for (String key : required) {
				if (prop.getProperty(key) == null) {
					throw new RuntimeException("OAFParameterSet: Missing property '" + key + "' in data file: " + filename);
				}
			}
		}

		return prop;
	}

	// load_file_as_string - Load a data file into a string.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.
	// This function first calls open_param_file to open the file,
	// then reads the entire contents of the file into a string.
	// In case of error, this function throws RuntimeException (or another unchecked exception).

	public static String load_file_as_string (String filename, Class<?> requester) {

		String result;

		// Any exception means load has failed

		try (

			// Open the data file, it will be auto-closed on exit from the try block

			InputStream stream = open_param_file (filename, requester);
		){

			// Convert the InputStream to a Reader
			// (Could specify a charset in the InputStreamReader constructor)
			// (Don't need to close these because they are just wrappers around the InputStream)

			Reader reader = new BufferedReader (new InputStreamReader (stream));

			// Buffer to build string

			StringBuilder builder = new StringBuilder();

			// Buffer to read file

			char[] buffer = new char[16384];

			// Loop to read file

			for (;;) {
				int amount = reader.read (buffer, 0, buffer.length);
				if (amount <= 0) {
					break;
				}
				builder.append (buffer, 0, amount);
			}

			// Get the string

			result = builder.toString();

		} catch (IOException e) {
			throw new RuntimeException("OAFParameterSet: Error reading from data file: " + filename, e);
		} catch (Exception e) {
			throw new RuntimeException("OAFParameterSet: Unable to load data file: " + filename, e);
		}

		// Done

		return result;
	}

	// load_file_as_json - Load a data file as JSON.
	// Parameters:
	//  filename = Name of file (not including a path).
	//  requester = Class that is requesting the file.
	// This function first calls open_param_file to open the file,
	// then reads the entire contents of the file into a MarshalReader.
	// In case of error, this function throws RuntimeException (or another unchecked exception).

	public static MarshalReader load_file_as_json (String filename, Class<?> requester) {

		MarshalImpJsonReader result;

		// Any exception means load has failed

		try (

			// Open the data file, it will be auto-closed on exit from the try block

			InputStream stream = open_param_file (filename, requester);
		){

			// Convert the InputStream to a Reader
			// (Could specify a charset in the InputStreamReader constructor)
			// (Don't need to close these because they are just wrappers around the InputStream)

			Reader reader = new BufferedReader (new InputStreamReader (stream));

			// Read the file as JSON

			result = new MarshalImpJsonReader (reader);

		} catch (IOException e) {
			throw new RuntimeException("OAFParameterSet: Error reading from data file: " + filename, e);
		} catch (Exception e) {
			throw new RuntimeException("OAFParameterSet: Unable to load data file: " + filename, e);
		}

		// Done

		return result;
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

	/**
	 * Find the tectonic regime for the given location, and return its parameters.
	 * @param loc = Location.
	 * @return Object of type T containing parameters.
	 */
	public T get(Location loc) {
		OAFTectonicRegime region = getRegion(loc);
		return get(region);
	}
	
	/**
	 * Return parameters for a tectonic regime, throw exception if regime is unknown.
	 * @param region = Tectonic regime.
	 * @return Object of type T containing parameters.
	 * The function throws an exception if no parameters are defined for the region.
	 */
	public T get(OAFTectonicRegime region) {
		T params = dataMap.get(region);
		if (params == null) {
			throw new RuntimeException("OAFParameterSet: Unknown tectonic regime : " + region);
		}
		return params;
	}
	
	/**
	 * Return parameters for a tectonic regime, or null if tectonic regime is unknown.
	 * @param region = Tectonic regime.
	 * @return Object of type T containing parameters.
	 * The function returns null if no parameters are defined for the region.
	 */
	public T getOrNull(OAFTectonicRegime region) {
		T params = dataMap.get(region);
		return params;
	}
	
	/**
	 * Find the tectonic regime for the given location.
	 * @param loc = Location.
	 * @return Tectonic regime for the location.
	 */
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

		synchronized (this) {

			// If the point is in a special region, use the tectonic regime for that region

			for (OAFRegion region : region_list) {
				if (region.contains (loc)) {
					return region.get_regime();
				}
			}

		}

		// If we have a world region, return it

		if (f_world) {
			return OAFTectonicRegime.forName (world_region);
		}

		// Otherwise, use the tectonic regime table to get the Garcia region

		TectonicRegimeTable regime_table = new TectonicRegimeTable();
		return OAFTectonicRegime.forName (regime_table.get_strec_name (loc.getLatitude(), loc.getLongitude()));
	}
	
	/**
	 * Return a set containing the tectonic regimes.
	 */
	public Set<OAFTectonicRegime> getRegimeSet() {
		return dataMap.keySet();
	}
	
	/**
	 * Return a list of locations that can be used for testing.
	 * The list includes a point in each of the 15 Garcia regions,
	 * plus a point in California, plus a few others.
	 */
	public static LocationList getTestLocations() {
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

		return locs;
	}

}
