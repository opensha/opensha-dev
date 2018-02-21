package scratch.aftershockStatistics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;

import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;


// TectonicRegimeTable is used to obtain the tectonic domain (also known as
// the Garcia region) for any point on the earth.
//
// TectonicRegimeTable contains a built-in table so that it does not need to
// communicate with any outside servers.  The table has a resolution of 0.05
// degrees in latitude and longitude, which is approximately 5 km.
//
// There is a single public function:
//
//		public String get_strec_name (double lat, double lon)
//
// Given latitude and longitude, expressed in degrees, get_strec_name returns a
// String that identifies the tectonic domain.  It always returns one of the strings
// listed in array strec_name below.  These are the same names used in STREC.
//
// The method get_strec_name is thread-safe.
//
// To use this class, first create an object of type TectonicRegimeTable.
// The data tables are loaded during the constructor, if needed.  Tables are cached
// in memory, so that only the first TectonicRegimeTable object created performs
// the load.  Then, call get_strec_name on the object.
//
// Created by Michael Barall on 02/10/2018.

public class TectonicRegimeTable {

	// lat_min - Minimum latitude covered by the table.

	private static final double lat_min = -90.0;

	// lat_max - Maximum latitude covered by the table.

	private static final double lat_max = 90.0;

	// lat_count - Number of different latitude values in the table.

	private static final int lat_count = 3601;	// 0.05 degree spacing, 20*180+1

	// lon_min - Minimum longitude covered by the table.
	// This must be either -180.0 or else -180.0 + point_spacing.

	private static final double lon_min = -180.0;

	// lon_max - Maximum longitude covered by the table.
	// This must be either 180.0 or else 180.0 - point_spacing.

	private static final double lon_max = 180.0;

	// lon_count - Number of different longitude values in the table.

	private static final int lon_count = 7201;	// 0.05 degree spacing, 20*360+1

	// domain_count - Number of tectonic domains.

	private static final int domain_count = 15;

	// strec_name - Tectonic region names, as known to STREC.
	// The number of elements in this array equals domain_count.

	private static final String[] strec_name = {
		"ACR (deep)",
		"ACR (hot spot)",
		"ACR (oceanic boundary)",
		"ACR (shallow)",
		"ACR deep (above slab)",
		"ACR oceanic boundary (above slab)",
		"ACR shallow (above slab)",
		"SCR (above slab)",
		"SCR (generic)",
		"SOR (above slab)",
		"SOR (generic)",
		"SZ (generic)",
		"SZ (inland/back-arc)",
		"SZ (on-shore)",
		"SZ (outer-trench)"
	};

	// point_count - The number of grid points in each tectonic domain.
	// The number of elements in this array equals domain_count.

	private static int[] point_count = null;

	// group_col - Group columns for the table.
	// The number of rows in this array equals lat_count.
	// Each row in this array corresponds to one latitude value.
	// Rows appear in order of increasing latitude from lat_min to lat_max.
	// Within a row, each element is the column number where a group begins.
	// Groups are listed in order of increasing column number.
	// Each group consists of successive grid points that lie in the same tectonic domain.

	private static int[][] group_col = null;

	// domain_index - Domain index numbers for the table.
	// The table structure is the same as group_col.
	// Each element gives the index number of the tectonic domain for the corresponding group.
	// The index number can be applied to strec_names to get the name of the domain.

	private static int[][] domain_index = null;

	// f_loaded - Flag, indicating if the tables have been loaded from the data file.

	private static boolean f_loaded = false;

	// get_row_for_lat - Get the table row number for the given latitude, in degrees.

	private static int get_row_for_lat (double lat) {

		// Get latitude into range -90.0 to 90.0

		double my_lat = lat;
		if (my_lat > 90.0) {
			my_lat = 90.0;
		}
		if (my_lat < -90.0) {
			my_lat = -90.0;
		}

		// Round to integer

		int row = (int)Math.round ( ((my_lat - lat_min) / (lat_max - lat_min)) * ((double)(lat_count - 1)) );

		// Clip to table size

		if (row < 0) {
			row = 0;
		}
		else if (row >= lat_count) {
			row = lat_count - 1;
		}
	
		return row;
	}

	// get_lat_for_row - Get the latitude in degrees for a given table row.

	private static double get_lat_for_row (int row) {
		return ( (((double)row) * (lat_max - lat_min)) / ((double)(lat_count - 1)) ) + lat_min;
	}

	// get_col_for_lon - Get the table column number for the given longitude, in degrees.

	private static int get_col_for_lon (double lon) {

		// Get longitude into range -180.0 to 180.0

		double my_lon = lon;
		while (my_lon > 180.0) {
			my_lon -= 360.0;
		}
		while (my_lon < -180.0) {
			my_lon += 360.0;
		}

		// Round to integer

		int col = (int)Math.round ( ((my_lon - lon_min) / (lon_max - lon_min)) * ((double)(lon_count - 1)) );

		// Clip to table size, with wraparound

		if (col < 0) {
			col = lon_count - 1;
		}
		else if (col >= lon_count) {
			col = 0;
		}
	
		return col;
	}

	// get_lon_for_col - Get the longitude in degrees for a given table column.

	private static double get_lon_for_col (int col) {
		return ( (((double)col) * (lon_max - lon_min)) / ((double)(lon_count - 1)) ) + lon_min;
	}

	// get_index_for_name - Get the domain index, given a STREC name.
	// Returns -1 if the name is not found.

	private static int get_index_for_name (String name) {
	
		// Search the list of names
		// (We use simple linear search since the list is short and this is only used while building the tables)

		for (int index = 0; index < domain_count; ++index) {
			if (name.equalsIgnoreCase (strec_name[index])) {
				return index;
			}
		}

		// Didn't find the name

		return -1;
	}

	// get_domain_index - Get the domain index number for the given latitude and longitude, in degrees.

	private static int get_domain_index (double lat, double lon) {

		// Get the row and column numbers

		int row = get_row_for_lat (lat);
		int col = get_col_for_lon (lon);

		// Binary search to find the group that contains our column

		int lo = 0;
		int hi = group_col[row].length;

		while (hi - lo > 1) {
			int mid = (hi + lo)/2;
			if (group_col[row][mid] <= col) {
				lo = mid;
			} else {
				hi = mid;
			}
		}

		// Return the index number
	
		return domain_index[row][lo];
	}

	// The constructor loads the tables if needed.

	public TectonicRegimeTable () {
		load_tables();
	}

	// get_strec_name - Get the STREC name for the given latitude and longitude, in degrees.

	public String get_strec_name (double lat, double lon) {

		// Get the domain index number
		
		int index = get_domain_index (lat, lon);

		// return the name
		
		return strec_name[index];
	}

	// build_tables - Build the tables.
	// This function creates the required data file.
	// The given input file must contain on each line:
	//  latitude  longitude  name
	// The three fields on the line are separated by spaces or tabs.
	// The name is not enclosed in quotes, even though it contains spaces.
	// The data file is written to standard output.
	// The data file format is:
	//	[domain_count repetitions] For each tectonic domain:
	//		[int]		The number of grid points in the tectonic domain.
	//	[lat_count repetitions] For each latitude, in order from lat_min to lat_max:
	//		[int]		The number of groups.
	//		[N int]		The column where each group begins, in increasing order.
	//		[N int]		The tectonic domain index number for each group.

	private static boolean build_tables (String filename) {
	
		// Open the file

		File file = new File(filename);
		Scanner sc;
		try {
			sc = new Scanner(file);
		} catch (FileNotFoundException e) {
			System.err.println ("ERROR - Cannot open file : " + filename);
			return false;
		}

		// Create work array for point_count, initialize each element to 0

		int row;
		int col;
		int index;

		int[] wk_point_count = new int[domain_count];
		for (index = 0; index < domain_count; ++index) {
			wk_point_count[index] = 0;
		}

		// Create a full-sized work array for domain_index, initialize each element to -1

		int[][] wk_domain_index = new int[lat_count][];
		for (row = 0; row < lat_count; ++row) {
			wk_domain_index[row] = new int[lon_count];
			for (col = 0; col < lon_count; ++col) {
				wk_domain_index[row][col] = -1;
			}
		}

		// Line counter

		int line_number = 0;

		// Loop over lines in the file

		while (sc.hasNextLine()) {

			// Count and retrieve the line, trimming leading and trailing whitespace

			++line_number;
			String line = sc.nextLine().trim();

			// If the line is non-empty ...

			if (!( line.isEmpty() )) {

				// Split the line into latitude, longitude, name

				String[] words = line.split("\\s+", 3);

				if (words.length != 3) {
					System.err.println ("ERROR:" + line_number + " - Badly formatted line : " + line);
					return false;
				}
			
				// Convert the latitude

				double lat;
				try {
					lat = Double.valueOf(words[0]).doubleValue();
				} catch (NumberFormatException e) {
					System.err.println ("ERROR:" + line_number + " - Invalid latitude : " + line);
					return false;
				}

				if (lat < lat_min - 0.0001 || lat > lat_max + 0.0001) {
					System.err.println ("ERROR:" + line_number + " - Out-of-range latitude : " + line);
					return false;
				}
			
				// Convert the longitude

				double lon;
				try {
					lon = Double.valueOf(words[1]).doubleValue();
				} catch (NumberFormatException e) {
					System.err.println ("ERROR:" + line_number + " - Invalid longitude : " + line);
					return false;
				}

				if (lon < lon_min - 0.0001 || lon > lon_max + 0.0001) {
					System.err.println ("ERROR:" + line_number + " - Out-of-range longitude : " + line);
					return false;
				}

				// Get the domain index for this name

				index = get_index_for_name (words[2]);
				if (index < 0) {
					System.err.println ("ERROR:" + line_number + " - Unknown region : " + line);
					return false;
				}

				// Get row and column

				row = get_row_for_lat (lat);
				col = get_col_for_lon (lon);

				// Error if this point is already in the table

				if (wk_domain_index[row][col] != -1) {
					System.err.println ("ERROR:" + line_number + " - Duplicate point : " + line);
					return false;
				}

				// Save the index for this point

				wk_domain_index[row][col] = index;

				// Count use of this index

				wk_point_count[index]++;
			}
		}

		// Close the file

		sc.close();

		// Check that all the points in the table are filled

		for (row = 0; row < lat_count; ++row) {
			for (col = 0; col < lon_count; ++col) {
				if (wk_domain_index[row][col] == -1) {
					double missing_lat = get_lat_for_row (row);
					double missing_lon = get_lon_for_col (col);
					System.err.println ("ERROR - Missing point : lat=" + missing_lat + " lon=" + missing_lon);
					return false;
				}
			}
		}

		// For each domain, write the grid point counts

		for (index = 0; index < domain_count; ++index) {
			StringBuffer buf = new StringBuffer();
			buf.append (wk_point_count[index]);
			System.out.println (buf);
		}

		// For each latitude, write the group list

		for (row = 0; row < lat_count; ++row) {
			StringBuffer buf_group = new StringBuffer();
			StringBuffer buf_col = new StringBuffer();
			StringBuffer buf_index = new StringBuffer();

			// Scan columns, and break them into groups with the same domain index

			col = 0;
			buf_col.append (col);
			buf_index.append (wk_domain_index[row][col]);

			int group_count = 1;

			for (col = 1; col < lon_count; ++col) {
				if (wk_domain_index[row][col] != wk_domain_index[row][col - 1]) {
					buf_col.append (" ");
					buf_index.append (" ");
					buf_col.append (col);
					buf_index.append (wk_domain_index[row][col]);
					++group_count;
				}
			}

			buf_group.append (group_count);
			System.out.println (buf_group);
			System.out.println (buf_col);
			System.out.println (buf_index);
		}

		// Successful table construction

		return true;
	}

	// load_table_int - Load an integer value for the tables.
	// An exception is thrown if there is no integer in the stream,
	// or if it lies outside the given minimum and maximum values.

	private static int load_table_int (Scanner sc, int minval, int maxval) {
	
		// Get the integer value

//		if (!( sc.hasNextInt() )) {
//			throw new RuntimeException("TectonicRegimeTable: Unexpected end-of-file or badly formatted integer value");
//		}

		int value;

		try {
			value = sc.nextInt();
		} catch (InputMismatchException e) {
			throw new RuntimeException("TectonicRegimeTable: Badly formatted integer value", e);
		} catch (NoSuchElementException e) {
			throw new RuntimeException("TectonicRegimeTable: Unexpected end-of-file", e);
		}

		// Range checking

		if (value < minval || value > maxval) {
			throw new RuntimeException("TectonicRegimeTable: Integer value out-of-range");
		}

		return value;
	}

	// load_tables - Load the tables into memory.
	// Tables are assumed to be stored as a resource named TectonicRegimeTable.dat.
	// An exception is thrown if the tables cannot be loaded.
	// Performs no operation if the tables are already loaded.

	private static synchronized void load_tables () {

		// If tables are already loaded, do nothing

		if (f_loaded) {
			return;
		}

		// Working tables
	
		int[] wk_point_count;
		int[][] wk_group_col;
		int[][] wk_domain_index;

		// Any exception means load has failed

		try {

			// Open the data file

			InputStream stream = TectonicRegimeTable.class.getResourceAsStream ("TectonicRegimeTable.dat");
			if (stream == null) {
				throw new RuntimeException("TectonicRegimeTable: Cannot find data file TectonicRegimeTable.dat");
			}

			Scanner sc = new Scanner (stream);

			// Make working copies of our tables
	
			wk_point_count = new int[domain_count];
			wk_group_col = new int[lat_count][];
			wk_domain_index = new int[lat_count][];

			// Load the point counts

			int remaining = lat_count * lon_count;
			for (int index = 0; index < domain_count; ++index) {
				wk_point_count[index] = load_table_int (sc, (index == domain_count - 1) ? remaining : 0, remaining);
				remaining -= wk_point_count[index];
			}

			// Load the group lists

			for (int row = 0; row < lat_count; ++row) {

				// Get number of groups and allocate the arrays

				int group_count = load_table_int (sc, 1, lon_count);
				wk_group_col[row] = new int[group_count];
				wk_domain_index[row] = new int[group_count];

				// Load the group columns

				int group;
				for (group = 0; group < group_count; ++group) {
					wk_group_col[row][group] = load_table_int (
						sc,
						(group > 0) ? (wk_group_col[row][group - 1] + 1) : 0,
						(group > 0) ? (lon_count - 1) : 0 );
				}

				// Load the domain indexes

				for (group = 0; group < group_count; ++group) {
					wk_domain_index[row][group] = load_table_int (sc, 0, domain_count - 1);
				}
			}

			// Close the file

			sc.close();

		} catch (Exception e) {
			throw new RuntimeException("TectonicRegimeTable: Unable to load data file TectonicRegimeTable.dat", e);
		}

		// Save our working arrays into the static variables

		point_count = wk_point_count;
		group_col = wk_group_col;
		domain_index = wk_domain_index;

		// Set flag indicating tables are loaded

		f_loaded = true;
		return;
	}

	// test_tables - Test the tables.
	// The given input file must contain on each line:
	//  latitude  longitude  name
	// The three fields on the line are separated by spaced or tabs.
	// The name is not enclosed in quotes, even though it contains spaces.
	// Test results are written to standard output.
	// The file can be the same file that was used to build the tables.
	// You call this function to test the definitions for rep_count and domain_index.

	private static boolean test_tables (String filename) {
	
		// Open the file

		File file = new File(filename);
		Scanner sc;
		try {
			sc = new Scanner(file);
		} catch (FileNotFoundException e) {
			System.err.println ("ERROR - Cannot open file : " + filename);
			return false;
		}

		// Object to use for queries

		TectonicRegimeTable regime_table = new TectonicRegimeTable();

		// Line counter

		int line_number = 0;
		int test_count = 0;
		int error_count = 0;

		// Loop over lines in the file

		while (sc.hasNextLine()) {

			// Count and retrieve the line, trimming leading and trailing whitespace

			++line_number;
			String line = sc.nextLine().trim();

			// If the line is non-empty ...

			if (!( line.isEmpty() )) {

				// Split the line into latitude, longitude, name

				String[] words = line.split("\\s+", 3);

				if (words.length != 3) {
					System.err.println ("ERROR:" + line_number + " - Badly formatted line : " + line);
					return false;
				}
			
				// Convert the latitude

				double lat;
				try {
					lat = Double.valueOf(words[0]).doubleValue();
				} catch (NumberFormatException e) {
					System.err.println ("ERROR:" + line_number + " - Invalid latitude : " + line);
					return false;
				}
			
				// Convert the longitude

				double lon;
				try {
					lon = Double.valueOf(words[1]).doubleValue();
				} catch (NumberFormatException e) {
					System.err.println ("ERROR:" + line_number + " - Invalid longitude : " + line);
					return false;
				}

				// Get the name for this latitude and longitude

				String name = regime_table.get_strec_name (lat, lon);
				++test_count;
				if (!( name.equalsIgnoreCase (words[2]) )) {
					++error_count;
					if (error_count < 11) {
						System.out.println ("ERROR:" + line_number + " - Wrong result : " + line);
						System.out.println ("Got name : " + name);
					}
					if (error_count == 11) {
						System.out.println ("Additional error messages suppressed");
					}
				}
			}
		}

		// Close the file

		sc.close();

		// Show results

		System.out.println ("Number of tests : " + test_count);
		System.out.println ("Number of errors : " + error_count);

		if (error_count == 0 && test_count > 0) {
			System.out.println ("SUCCESS - All tests completed successfully");
		}

		return true;
	}

	// Entry point.
	
	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("TectonicRegimeTable : Missing subcommand");
			return;
		}

		// Subcommand : Build the data tables
		// Command format:
		//  build  filename

		if (args[0].equalsIgnoreCase ("build")) {

			// Second argument is filename

			if (args.length != 2) {
				System.err.println ("TectonicRegimeTable : Invalid 'build' subcommand");
				return;
			}

			// Invoke the operation

			build_tables (args[1]);

			return;
		}
		
		// Subcommand : Test the data tables
		// Command format:
		//  test  filename

		if (args[0].equalsIgnoreCase ("test")) {

			// Second argument is filename

			if (args.length != 2) {
				System.err.println ("TectonicRegimeTable : Invalid 'test' subcommand");
				return;
			}

			// Invoke the operation

			test_tables (args[1]);

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("TectonicRegimeTable : Unrecognized subcommand : " + args[0]);
		return;
	}

}
