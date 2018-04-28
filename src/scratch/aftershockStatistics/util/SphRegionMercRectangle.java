package scratch.aftershockStatistics.util;

import java.util.List;
import java.util.ArrayList;

import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.BorderType;

import static java.lang.Math.PI;
import static org.opensha.commons.geo.GeoTools.TWOPI;
import static org.opensha.commons.geo.GeoTools.TO_DEG;
import static org.opensha.commons.geo.GeoTools.TO_RAD;
import static org.opensha.commons.geo.GeoTools.EARTH_RADIUS_MEAN;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;


/**
 * Region for locating aftershocks -- Mercator Rectangle.
 * Author: Michael Barall 04/10/2018.
 *
 * The AAFS server locates aftershocks by drawing a region surrounding the
 * mainshock.  Earthquakes within that region are considered to be aftershocks.
 * This class represents a rectangle on a Mercator map.
 * The rectangle is specified by giving two diagonally opposite corners.
 * A rectangle is limited to 180 degrees of longitude.
 * Attempting to make a rectangle spanning more than 180 degrees of longitude
 * will make a rectangle that goes around the Earth the "other way".
 */
public class SphRegionMercRectangle extends SphRegion {

	//----- Querying -----

	/**
	 * contains - Test if the region contains the given location.
	 * @param loc = Location to check.
	 * @return
	 * Returns true if loc is inside the region, false if loc is outside the region.
	 * Note: Due to rounding errors, it may be indeterminate whether points exactly on,
	 * or very close to, the boundary of the region are considered inside or outside.
	 */
	@Override
	public boolean contains (SphLatLon loc) {
		return contains (loc.get_lat(), loc.get_lon());
	}

	/**
	 * contains - Test if the region contains the given location.
	 * @param loc = Location to check.
	 * @return
	 * Returns true if loc is inside the region, false if loc is outside the region.
	 * Note: Due to rounding errors, it may be indeterminate whether points exactly on,
	 * or very close to, the boundary of the region are considered inside or outside.
	 */
	@Override
	public boolean contains (Location loc) {
		return contains (loc.getLatitude(), loc.getLongitude());
	}

	/**
	 * contains - Test if the region contains the given location.
	 * @param lat = Latitude to check.
	 * @param lon = Longitude to check, can be -180 to +360.
	 * @return
	 * Returns true if loc is inside the region, false if loc is outside the region.
	 * Note: Due to rounding errors, it may be indeterminate whether points exactly on,
	 * or very close to, the boundary of the region are considered inside or outside.
	 * Implementation note: The function uses plane geometry, in the domain selected
	 * by plot_wrap.
	 */
	private boolean contains (double lat, double lon) {

		// Coerce longitude according to our wrapping domain.

		if (plot_wrap) {
			if (lon < 0.0) {
				lon += 360.0;
			}
		} else {
			if (lon > 180.0) {
				lon -= 360.0;
			}
		}

		// Now just compare to box limits

		return lon >= min_lon && lon <= max_lon && lat >= min_lat && lat <= max_lat;
	}




	//----- Plotting -----

	// This function is called when it is necessary to build plot_border.
	// The subclass should supply a function to build the border.
	// A subclass might choose to build the border in its constructor and
	// do_marshal method, in which case this function can just throw an exception.

	@Override
	protected void make_plot_border () {
		plot_border = new LocationList();
	
		plot_border.add(new Location(min_lat, min_lon));
		plot_border.add(new Location(min_lat, max_lon));
		plot_border.add(new Location(max_lat, max_lon));
		plot_border.add(new Location(max_lat, min_lon));

		return;
	}




	//----- Special regions -----

	/**
	 * Return true if this is a rectangular region in a Mercator projection.
	 * If this function returns true, then the region is exactly the box
	 * given by min_lat, etc.
	 */
	@Override
	public boolean isRectangular() {
		return true;
	}




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public SphRegionMercRectangle () {}


	/**
	 * Construct from given corners.
	 * @param corner_1 = Rectangle corner #1.
	 * @param corner_2 = Rectangle corner #2, must be diagonally opposite to corner #1.
	 */
	public SphRegionMercRectangle (SphLatLon corner_1, SphLatLon corner_2) {
		setup (corner_1, corner_2);
	}


	/**
	 * Set up the region.
	 */
	private void setup (SphLatLon corner_1, SphLatLon corner_2) {

		// Get the latitudes and longitudes

		double lat_1 = corner_1.get_lat();
		double lon_1 = corner_1.get_lon();

		double lat_2 = corner_2.get_lat();
		double lon_2 = corner_2.get_lon();

		// Set up a box in the -180 to +180 domain

		plot_wrap = false;

		min_lat = Math.min (lat_1, lat_2);
		max_lat = Math.max (lat_1, lat_2);

		min_lon = Math.min (lon_1, lon_2);
		max_lon = Math.max (lon_1, lon_2);

		// If it spans more than 180 degrees, use the 0 to 360 domain to go around the other way

		if (max_lon - min_lon > 180.0) {
			plot_wrap = true;
		
			min_lon = Math.max (lon_1, lon_2);
			max_lon = Math.min (lon_1, lon_2) + 360.0;
		}

		plot_border = null;

		return;
	}

	// Display our contents

	@Override
	public String toString() {
		return "SphRegionMercRectangle:" + "\n"
		+ "plot_wrap = " + plot_wrap + "\n"
		+ "min_lat = " + min_lat + "\n"
		+ "max_lat = " + max_lat + "\n"
		+ "min_lon = " + min_lon + "\n"
		+ "max_lon = " + max_lon;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_HWV_1 = 1;		// human-writeable version
	private static final int MARSHAL_VER_1 = 20001;

	private static final String M_VERSION_NAME = "SphRegionMercRectangle";

	// Get the type code.

	@Override
	protected int get_marshal_type () {
		return MARSHAL_CIRCLE;
	}

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents


		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME);

		switch (ver) {

		default:
			throw new MarshalException ("SphRegionMercRectangle.do_umarshal: Unknown version code: version = " + ver);
		
		// Human-writeable version

		case MARSHAL_HWV_1: {

			// Get corners

			SphLatLon corner_1 = (new SphLatLon()).unmarshal (reader, "corner_1");;
			SphLatLon corner_2 = (new SphLatLon()).unmarshal (reader, "corner_2");;

			// Set up region

			try {
				setup (corner_1, corner_2);
			}
			catch (Exception e) {
				throw new MarshalException ("SphRegionMercRectangle.do_umarshal: Failed to set up region", e);
			}
		}
		break;

		// Machine-written version

		case MARSHAL_VER_1: {

			// Superclass

			super.do_umarshal (reader);

			// Contents

		}
		break;

		}

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
	public SphRegionMercRectangle unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, SphRegionMercRectangle obj) {

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

	public static SphRegionMercRectangle unmarshal_poly (MarshalReader reader, String name) {
		SphRegionMercRectangle result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("SphRegionMercRectangle.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_MERC_RECTANGLE:
			result = new SphRegionMercRectangle();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

}
