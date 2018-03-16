package scratch.aftershockStatistics;

import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Set;

import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;

import scratch.aftershockStatistics.OAFTectonicRegime;
import scratch.aftershockStatistics.OAFRegion;
import scratch.aftershockStatistics.OAFParameterSet;

public class MagCompPage_ParametersFetch {

	// parameter_set - The parameter set.

	private static OAFParameterSet<MagCompPage_Parameters> parameter_set = null;

	// load_data - Load parameters from the data file.
	// The data file format is:
	//	[int]		Number of tectonic regimes
	//	[repeated]	Repeated once for each tectonic regime:
	//		[string]	Name of tectonic regime
	//		[double]	magCat
	//		[double]	capG
	//		[double]	capH
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

	private static synchronized void load_data () {

		// If data is already loaded, do nothing

		if (parameter_set != null) {
			return;
		}

		// Working data

		OAFParameterSet<MagCompPage_Parameters> wk_parameter_set = new OAFParameterSet<MagCompPage_Parameters>(){
			
			// load_parameter_values - Load parameter values for the tables.
			// This function should create a new object of type T, read the
			// parameter values from the Scanner, and return the object.
			// In case of error, this function should throw RuntimeException.

			@Override
			protected MagCompPage_Parameters load_parameter_values (Scanner sc) {

				// Get the Page parameters
				
				double magCat = load_table_double (sc);
				double capG = load_table_double (sc);
				double capH = load_table_double (sc);

				// Make the parameter object

				return new MagCompPage_Parameters(magCat, capG, capH);
			}
		};

		// Load the data

		wk_parameter_set.load_data ("MagCompPage_ParametersFetch.txt", MagCompPage_ParametersFetch.class);

		// Save our working data into the static variable

		parameter_set = wk_parameter_set;
		return;
	}

	// Constructor loads the data if needed.
	
	public MagCompPage_ParametersFetch() {
		load_data();
	}
	
	/**
	 * Find the tectonic regime for the given location, and return its parameters.
	 * @param loc = Location.
	 * @return Object of type T containing parameters.
	 */
	public MagCompPage_Parameters get(Location loc) {
		return parameter_set.get(loc);
	}
	
	/**
	 * Return parameters for a tectonic regime, throw exception if regime is unknown.
	 * @param region = Tectonic regime.
	 * @return Object of type T containing parameters.
	 * The function throws an exception if no parameters are defined for the region.
	 */
	public MagCompPage_Parameters get(OAFTectonicRegime region) {
		return parameter_set.get(region);
	}
	
	/**
	 * Return parameters for a tectonic regime, or null if tectonic regime is unknown.
	 * @param region = Tectonic regime.
	 * @return Object of type T containing parameters.
	 * The function returns null if no parameters are defined for the region.
	 */
	public MagCompPage_Parameters getOrNull(OAFTectonicRegime region) {
		return parameter_set.getOrNull(region);
	}
	
	/**
	 * Find the tectonic regime for the given location.
	 * @param loc = Location.
	 * @return Tectonic regime for the location.
	 */
	public OAFTectonicRegime getRegion(Location loc) {
		return parameter_set.getRegion(loc);
	}
	
	/**
	 * Return a set containing the tectonic regimes.
	 */
	public Set<OAFTectonicRegime> getRegimeSet() {
		return parameter_set.getRegimeSet();
	}



	
	public static void main(String[] args) {

		// There needs to be at least one argument, which is the subcommand

		if (args.length < 1) {
			System.err.println ("MagCompPage_ParametersFetch : Missing subcommand");
			return;
		}

		// Subcommand : Test #1
		// Command format:
		//  test1
		// Display tectonic regime and parameters for a selected list of points.

		if (args[0].equalsIgnoreCase ("test1")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("MagCompPage_ParametersFetch : Invalid 'test1' subcommand");
				return;
			}

			// Display info for each point in a list of locations

			LocationList locs = OAFParameterSet.getTestLocations();;
		
			MagCompPage_ParametersFetch fetch = new MagCompPage_ParametersFetch();
			for (Location loc : locs) {
				OAFTectonicRegime regime = fetch.getRegion(loc);
				System.out.println(loc+", "+regime+": "+fetch.get(regime));
			}

			return;
		}

		// Subcommand : Test #2
		// Command format:
		//  test1
		// List all regimes, and the properties assigned to each.

		if (args[0].equalsIgnoreCase ("test2")) {

			// No additional arguments

			if (args.length != 1) {
				System.err.println ("MagCompPage_ParametersFetch : Invalid 'test2' subcommand");
				return;
			}

			// Display info for each regime
		
			MagCompPage_ParametersFetch fetch = new MagCompPage_ParametersFetch();
			Set<OAFTectonicRegime> regimes = fetch.getRegimeSet();
			for (OAFTectonicRegime regime : regimes) {
				System.out.println(regime+": "+fetch.get(regime));
			}

			return;
		}

		// Unrecognized subcommand.

		System.err.println ("MagCompPage_ParametersFetch : Unrecognized subcommand : " + args[0]);
		return;

	}

}
