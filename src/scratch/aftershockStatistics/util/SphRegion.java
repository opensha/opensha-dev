package scratch.aftershockStatistics.util;

import java.util.List;
import java.util.ArrayList;

import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GeoTools;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;


/**
 * Region on a sphere.
 * Author: Michael Barall 04/12/2018.
 *
 * The OpenSHA class Region performs plane geometry on an equirectangular
 * projection of the Earth.  Attempting to use Region in AAFS has proven to be
 * error-prone.  This is an abstract class that represents a region on a sphere,
 * with calculations done using spherical geometry.  Subclasses define specific
 * types of region (circle, rectangle, polygon, etc.).
 */
public abstract class SphRegion {

	//----- Querying -----

	/**
	 * contains - Test if the region contains the given location.
	 * @param loc = Location to check.
	 * @return
	 * Returns true if loc is inside the region, false if loc is outside the region.
	 * Note: Due to rounding errors, it may be indeterminate whether points exactly on,
	 * or very close to, the boundary of the region are considered inside or outside.
	 */
	public abstract boolean contains (SphLatLon loc);

	/**
	 * contains - Test if the region contains the given location.
	 * @param loc = Location to check.
	 * @return
	 * Returns true if loc is inside the region, false if loc is outside the region.
	 * Note: Due to rounding errors, it may be indeterminate whether points exactly on,
	 * or very close to, the boundary of the region are considered inside or outside.
	 */
	public abstract boolean contains (Location loc);




	//----- Plotting -----

	// OpenSHA allows longitude to range from -180 to +360, which means that
	// each point in the western hemisphere has two possible longitudes.  Since
	// SphRegion does spherical geometry, SphRegion and related classes only permit
	// longitude to range from -180 to +180.  This can be a problem when plotting
	// regions that cross the date line.  We address this by defining two plotting
	// domains: one that extends from -180 to +180, and one that extends from 0 to +360.
	// Regions that cross the date line are plotted in the second domain.  Note that
	// both domains lie within the OpenSHA permitted longitude range.
	//
	// If plot_wrap is false, then plots should be drawn by coercing longitudes
	// to lie between -180 and +180.  If plot_wrap is true, then plots should
	// be drawn by coercing longitudes to lie between 0 and +360.
	// This field should be set up by the subclass.

	protected boolean plot_wrap;

	/**
	 * Returns the plot domain.
	 */
	public boolean getPlotWrap() {
		return plot_wrap;
	}

	// These fields define a latitude/longitude bounding box for the region.
	// When plotting, this box can be used to set the limits of the plot.
	// The latitude values must satisfy:
	//  -90 <= min_lat <= max_lat <= +90
	// If plot_wrap is false, then the longitude values must satisfy:
	//  -180 <= min_lon <= max_lon <= +180
	// If plot_wrap is true, then the longitude values must satisfy:
	//  -0 <= min_lon <= max_lon <= +360
	// These fields should be set up by the subclass.  A subclass may also
	// choose to use these fields for a "quick reject" of points that are
	// well outside the region.

	protected double min_lat;
	protected double max_lat;
	protected double min_lon;
	protected double max_lon;

	/**
	 * Returns the minimum latitude.
	 */
	public double getMinLat() {
		return min_lat;
	}

	/**
	 * Returns the maximum latitude.
	 */
	public double getMaxLat() {
		return max_lat;
	}

	/**
	 * Returns the minimum longitude.
	 */
	public double getMinLon() {
		return min_lon;
	}

	/**
	 * Returns the maximum longitude.
	 */
	public double getMaxLon() {
		return max_lon;
	}

	// A list of points that is used for drawing the region's border on a plot.
	// Each point must be within the bounding box defined by min_lat, etc.
	// This is initialized to null, and built when needed.
	// The last point should not coincide with the first point.
	//
	// Note: This field is typically not marshaled because it can be re-built
	// when needed. The subclass constructor and do_unmarshal method must either
	// set plot_border to null, or else build the list.  If they build the list,
	// than make_plot_border can just throw an exception.  This class's
	// constructor and do_unmarshal method set plot_border to null.

	protected LocationList plot_border = null;

	// This function is called when it is necessary to build plot_border.
	// The subclass should supply a function to build the border.
	// A subclass might choose to build the border in its constructor and
	// do_marshal method, in which case this function can just throw an exception.

	protected abstract void make_plot_border ();

	// The subclass can call this function to add a point to the border.
	// The longitude is adjusted to the plotting domain, then the latitude
	// and longitude are clipped to the bounding box.

	protected void add_border_point (SphLatLon loc) {
		double lat = loc.get_lat();
		double lon = loc.get_lon(plot_wrap);

		if (lat < min_lat) {
			lat = min_lat;
		} else if (lat > max_lat) {
			lat = max_lat;
		}

		if (lon < min_lon) {
			lon = min_lon;
		} else if (lon > max_lon) {
			lon = max_lon;
		}

		plot_border.add (new Location (lat, lon));
		return;
	}

	/**
	 * Returns an unmodifiable view of the border.
	 */
	public LocationList getBorder() {
		if (plot_border == null) {
			make_plot_border();
		}
		return plot_border.unmodifiableList();
	}




	//----- Special regions -----

	/**
	 * Return true if this is a rectangular region in a Mercator projection.
	 * If this function returns true, then the region is exactly the box
	 * given by min_lat, etc.
	 */
	public boolean isRectangular() {
		return false;
	}

	/**
	 * Return true if this is a circular region on the sphere.
	 * If this function returns true, then the center and radius can be
	 * obtained from getCircleCenter and getCircleRadiusDeg.
	 */
	public boolean isCircular() {
		return false;
	}

	/**
	 * If this is a circular region on the sphere, then get the center.
	 * Otherwise, throw an exception.
	 */
	public SphLatLon getCircleCenter() {
		throw new UnsupportedOperationException ("The region is not a circle");
	}

	/**
	 * If this is a circular region on the sphere, then get the radius in degrees.
	 * The returned value ranges from 0 to +180.
	 */
	public double getCircleRadiusDeg() {
		throw new UnsupportedOperationException ("The region is not a circle");
	}

	/**
	 * Return true if this region is the entire world.
	 * If this function returns true, then the region is exactly the box
	 * given by min_lat, etc., and that box has latitude -90 to +90
	 * and longitude -180 to +180.
	 */
	public boolean isWorld() {
		return false;
	}




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public SphRegion () {}

	// Display our contents

	@Override
	public String toString() {
		return "SphRegion:" + "\n"
		+ "plot_wrap = " + plot_wrap + "\n"
		+ "min_lat = " + min_lat + "\n"
		+ "max_lat = " + max_lat + "\n"
		+ "min_lon = " + min_lon + "\n"
		+ "max_lon = " + max_lon;
	}




	//----- Factory methods -----


	/**
	 * Construct a circle from given center location and radius.
	 * @param center = Center of the circle.
	 * @param radius = Radius of the circle, in kilometers, must be positive and span no more than 180 degrees.
	 */
	public static SphRegion makeCircle (SphLatLon center, double radius) {
		return new SphRegionCircle (center, radius);
	}


	/**
	 * Construct a Mercator rectangle from given corners.
	 * @param corner_1 = Rectangle corner #1.
	 * @param corner_2 = Rectangle corner #2, must be diagonally opposite to corner #1.
	 */
	public static SphRegion makeMercRectangle (SphLatLon corner_1, SphLatLon corner_2) {
		return new SphRegionMercRectangle (corner_1, corner_2);
	}


	/**
	 * Construct a Mercator polygon from given list of vertices.
	 * @param vertex_list = List of vertices.
	 */
	public static SphRegion makeMercPolygon (List<SphLatLon> vertex_list) {
		return new SphRegionMercPolygon (vertex_list);
	}


	/**
	 * Construct a region that is the entire world.
	 */
	public static SphRegion makeWorld () {
		return new SphRegionWorld ();
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 12001;

	private static final String M_VERSION_NAME = "SphRegion";

	// Marshal type code.

	protected static final int MARSHAL_NULL = 12000;
	protected static final int MARSHAL_CIRCLE = 13001;
	protected static final int MARSHAL_MERC_POLYGON = 14001;
	protected static final int MARSHAL_MERC_RECTANGLE = 20001;
	protected static final int MARSHAL_GC_POLYGON = 21001;
	protected static final int MARSHAL_WORLD = 33001;

	protected static final String M_TYPE_NAME = "ClassType";

	// Get the type code.

	protected abstract int get_marshal_type ();

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		writer.marshalBoolean ("plot_wrap", plot_wrap);
		writer.marshalDouble  ("min_lat"  , min_lat  );
		writer.marshalDouble  ("max_lat"  , max_lat  );
		writer.marshalDouble  ("min_lon"  , min_lon  );
		writer.marshalDouble  ("max_lon"  , max_lon  );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		plot_wrap = reader.unmarshalBoolean ("plot_wrap");
		min_lat   = reader.unmarshalDouble  ("min_lat"  );
		max_lat   = reader.unmarshalDouble  ("max_lat"  );
		min_lon   = reader.unmarshalDouble  ("min_lon"  );
		max_lon   = reader.unmarshalDouble  ("max_lon"  );

		plot_border = null;

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

	public SphRegion unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, SphRegion obj) {

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

	public static SphRegion unmarshal_poly (MarshalReader reader, String name) {
		SphRegion result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("SphRegion.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_CIRCLE:
			result = new SphRegionCircle();
			result.do_umarshal (reader);
			break;

		case MARSHAL_MERC_POLYGON:
			result = new SphRegionMercPolygon();
			result.do_umarshal (reader);
			break;

		case MARSHAL_MERC_RECTANGLE:
			result = new SphRegionMercRectangle();
			result.do_umarshal (reader);
			break;

		case MARSHAL_WORLD:
			result = new SphRegionWorld();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

}
