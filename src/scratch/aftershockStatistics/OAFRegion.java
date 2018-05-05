package scratch.aftershockStatistics;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;

import scratch.aftershockStatistics.OAFTectonicRegime;

// OAFRegion is a region for defining a tectonic regime.
// It includes a 3D region on and below the earth's surface,
// which may be limited in depth as well as in horizontal extent.
// It also includes the tectonic regime.

public abstract class OAFRegion {

	// MIN_DEPTH_UNBOUNDED - Value to use for minimum depth when there is no bound.

	public static final double MIN_DEPTH_UNBOUNDED = -1.0e10;

	// MAX_DEPTH_UNBOUNDED - Value to use for maximum depth when there is no bound.

	public static final double MAX_DEPTH_UNBOUNDED = 1.0e10;

	// regime - The tectonic regime.

	private OAFTectonicRegime regime;

	// get_regime - Gets the tectonic regime.

	public OAFTectonicRegime get_regime () {
		return regime;
	}

	// Constructor saves the regime.

	protected OAFRegion (OAFTectonicRegime the_regime) {
		regime = the_regime;
	}

	// contains - Determine whether the given location is inside the region.

	public abstract boolean contains (Location loc);

}
