package scratch.aftershockStatistics;

import java.util.Map;

import org.opensha.commons.data.siteData.impl.TectonicRegime;

import java.util.HashMap;

// OAFTectonicRegime is intended to be a drop-in replacement for TectonicRegime
// (in package org.opensha.commons.data.siteData.impl) for use in OAF code.
//
// The important feature of OAFTectonicRegime is that the set of tectonic
// regimes can be enlarged dynamically at run-time.  The 15 Garcia regions
// are added automatically, each with three names:  the Garcia name (with an
// underscore in place of a dash), the Garcia name (with the dash), and
// the STREC name.  The first name is the one that appears when converting
// the regime to a string (in order to match TectonicRegime).  Dynamically-
// added regimes have only one name.
//
// OAFTectonicRegime uses the flyweight design pattern.  That means there is
// only one object in the system with a given value.  So, objects x and y refer
// to the same tectonic regime if and only if x==y.  In many respects,
// OAFTectonicRegime objects can be used as if they were enums.

public final class OAFTectonicRegime {

	// mappings - Maps names to tectonic regimes.

	private static Map<String, OAFTectonicRegime> mappings = null;

	// regime_name - The name of this tectonic regime.

	private String regime_name;

	// Conversion to string returns our regime name.

	@Override
	public String toString() {
		return regime_name;
	}

	// Private constructor.

	private OAFTectonicRegime (String name) {
		regime_name = name;
	}
	
	private static void add_garcia_region (TectonicRegime regime) {
		add_garcia_region(regime.getNames());
	}

	// add_garcia_region - Add a Garcia region to the table.
	// It is assumed that a caller has synchronized.

	private static void add_garcia_region (String... names) {
		OAFTectonicRegime regime = new OAFTectonicRegime (names[0]);
		for (String name : names)
			mappings.put (name, regime);
		return;
	}

	// setup_mappings - Create and set up the mappings table.
	// It is assumed that a caller has synchronized.

	private static void setup_mappings () {

		// Create the table

		mappings = new HashMap<String, OAFTectonicRegime>();

		// Add the Garcia regions

		add_garcia_region (TectonicRegime.ANSR_DEEPCON);
		add_garcia_region (TectonicRegime.ANSR_HOTSPOT);
		add_garcia_region (TectonicRegime.ANSR_OCEANBD);
		add_garcia_region (TectonicRegime.ANSR_SHALCON);
		add_garcia_region (TectonicRegime.ANSR_ABSLDEC);
		add_garcia_region (TectonicRegime.ANSR_ABSLOCB);
		add_garcia_region (TectonicRegime.ANSR_ABSLSHC);
		add_garcia_region (TectonicRegime.SCR_ABVSLAB);
		add_garcia_region (TectonicRegime.SCR_GENERIC);
		add_garcia_region (TectonicRegime.SOR_ABVSLAB);
		add_garcia_region (TectonicRegime.SOR_GENERIC);
		add_garcia_region (TectonicRegime.SZ_GENERIC);
		add_garcia_region (TectonicRegime.SZ_INLBACK);
		add_garcia_region (TectonicRegime.SZ_ONSHORE);
		add_garcia_region (TectonicRegime.SZ_OUTERTR);

		return;
	}

	// forName - Get the tectonic regime with the given name.
	// If there is no regime for the given name, create one.
	
	public static synchronized OAFTectonicRegime forName (String name) {

		// Set up the table if needed

		if (mappings == null) {
			setup_mappings();
		}

		// Get the regime from the table

		OAFTectonicRegime regime = mappings.get(name);

		// If we didn't find a regime, create one

		if (regime == null) {
			regime = new OAFTectonicRegime (name);
			mappings.put (name, regime);
		}

		return regime;
	}

}
