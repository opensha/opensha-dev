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
public class SphRegionMercPolygon extends SphRegion {

	//----- Region definition -----

	// Latitude and longitude of the vertices, in degrees.
	// There must be at least 3 vertices.

	private ArrayList<SphLatLon> vertex_list;




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
	 * by plot_wrap.  It walks the vertex list, keeping track of how many edges intersect
	 * a ray extending due north from the query point.  The point is considered to be
	 * inside if there are an odd number of intersections.
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

		// If outside the box limits, we can return false.
		// This check is necessary so that we don't wrap around the end of the domain
		// while finding intersections. The algorithm needs to know if each vertext is
		// to the left or to the right of the query point.

		if (!( lon >= min_lon && lon <= max_lon && lat >= min_lat && lat <= max_lat )) {
			return false;
		}

		// Now introduce an (x,y) coordinate system where x increases to east from the query
		// point, and y increases to north from the query point.  Also, s is the sign of x.

		double x1;
		double y1;
		int s1;
		double x2;
		double y2;
		int s2;
		double x_last;
		double y_last;
		int s_last;

		int n = vertex_list.size();

		// Let (x1,y1) and (x2,y2) be the coordinates of two successive vertices.
		// It counts as a crossing if the edge connecting the vertices intersects the
		// positive y-axis.  For that to happen, x1 and x2 must have opposite signs,
		// and the y-coordinate where it intersects the y-axis must be positive.
		// The standard linear interpolation formula says that the y-coordinate where
		// it intersects the y-axis (x == 0) is:
		//  y = (x2*y1 - x1*y2) / (x2 - x1)
		// So, if x2 > x1, then a crossing from left to right occurs if x2*y1 - x1*y2 > 0.
		// If x2 < x1, then a crossing from right to left occurs if x2*y1 - x1*y2 < 0.

		// Get (x,y) from the last point in the list

		x2 = vertex_list.get(n-1).get_lat() - lat;
		y2 = vertex_list.get(n-1).get_lon(plot_wrap) - lon;
		s2 = ((x2 >= 0.0) ? 1 : -1);

		x_last = x2;
		y_last = y2;
		s_last = s2;

		// Winding number, counts number of left-to-right crossings

		int winding = 0;

		// Loop over vertices, each edge runs from prior vertex to this vertex

		for (int i = 0; i < n; ++i) {

			// Endpoints of the current edge

			x1 = x2;
			y1 = y2;
			s1 = s2;

			if (i == n-1) {
				x2 = x_last;
				y2 = y_last;
				s2 = s_last;
			} else {
				x2 = vertex_list.get(i).get_lat() - lat;
				y2 = vertex_list.get(i).get_lon(plot_wrap) - lon;
				s2 = ((x2 >= 0.0) ? 1 : -1);
			}

			// If possible left-to-right crossing ...

			if (s2 > s1) {

				// Count if the crossing point is above the query latitude

				if (x2 * y1 > x1 * y2) {
					++winding;
				}
			}

			// If possible right-to-left crossing ...

			if (s2 < s1) {

				// Count if the crossing point is above the query latitude

				if (x2 * y1 < x1 * y2) {
					--winding;
				}
			}
		}

		return (winding % 2) != 0;
	}




	//----- Plotting -----

	// This function is called when it is necessary to build plot_border.
	// The subclass should supply a function to build the border.
	// A subclass might choose to build the border in its constructor and
	// do_marshal method, in which case this function can just throw an exception.

	@Override
	protected void make_plot_border () {
		plot_border = new LocationList();

		for (SphLatLon vertex : vertex_list) {
			add_border_point (vertex);
		}

		return;
	}




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public SphRegionMercPolygon () {}


	/**
	 * Construct from given list of vertices.
	 * @param vertex_list = List of vertices.
	 */
	public SphRegionMercPolygon (List<SphLatLon> vertex_list) {
		setup (vertex_list);
	}


	/**
	 * Set up the region.
	 */
	private void setup (List<SphLatLon> the_vertex_list) {

		// Must be at least 3 points

		int n = the_vertex_list.size();
		if (n < 3) {
			throw new IllegalArgumentException ("SphRegionMercPolygon: Must have at least 3 vertices: n = " + n);
		}

		// If the first and last vertices are the same, drop the last vertex

		if (the_vertex_list.get(0).get_lat() == the_vertex_list.get(n-1).get_lat()
			&& the_vertex_list.get(0).get_lon() == the_vertex_list.get(n-1).get_lon()) {
			--n;
			if (n < 3) {
				throw new IllegalArgumentException ("SphRegionMercPolygon: Must have at least 3 vertices: n = " + n);
			}
		}

		// Copy the list

		vertex_list = new ArrayList<SphLatLon>();

		for (int i = 0; i < n; ++i) {
			vertex_list.add (the_vertex_list.get(i));
		}

		// Scan list to find min and max latitude, longitude, and wrapped longitude

		plot_wrap = false;

		min_lat = Double.MAX_VALUE;
		max_lat = -Double.MAX_VALUE;

		min_lon = Double.MAX_VALUE;
		max_lon = -Double.MAX_VALUE;

		double min_wrap_lon = Double.MAX_VALUE;
		double max_wrap_lon = -Double.MAX_VALUE;

		for (SphLatLon loc : vertex_list) {
			double lat = loc.get_lat();
			double lon = loc.get_lon();
			double wrap_lon = loc.get_lon(true);

			min_lat = Math.min(min_lat, lat);
			max_lat = Math.max(max_lat, lat);

			min_lon = Math.min(min_lon, lon);
			max_lon = Math.max(max_lon, lon);

			min_wrap_lon = Math.min(min_wrap_lon, wrap_lon);
			max_wrap_lon = Math.max(max_wrap_lon, wrap_lon);
		}

		// If longitude spans more than 180 degress, use the wrapped longitude instead

		if (max_lon - min_lon > 180.0) {
			plot_wrap = true;
			min_lon = min_wrap_lon;
			max_lon = max_wrap_lon;

			// If longitude still spans more than 180 degrees, it's an error

			if (max_lon - min_lon > 180.0) {
				throw new IllegalArgumentException ("SphRegionMercPolygon: Region spans more than 180 degrees in longitude");
			}
		}

		plot_border = null;

		return;
	}

	// Display our contents

	@Override
	public String toString() {
		return "SphRegionMercPolygon:" + "\n"
		+ "plot_wrap = " + plot_wrap + "\n"
		+ "min_lat = " + min_lat + "\n"
		+ "max_lat = " + max_lat + "\n"
		+ "min_lon = " + min_lon + "\n"
		+ "max_lon = " + max_lon + "\n"
		+ "vertex_list size = " + vertex_list.size();
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_HWV_1 = 1;		// human-writeable version
	private static final int MARSHAL_VER_1 = 14001;

	private static final String M_VERSION_NAME = "SphRegionMercPolygon";

	// Get the type code.

	@Override
	protected int get_marshal_type () {
		return MARSHAL_MERC_POLYGON;
	}

	// Marshal object, internal.

	@Override
	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Superclass

		super.do_marshal (writer);

		// Contents

		SphLatLon.marshal_list (writer, "vertex_list", vertex_list);

		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME);

		switch (ver) {

		default:
			throw new MarshalException ("SphRegionMercPolygon.do_umarshal: Unknown version code: version = " + ver);
		
		// Human-writeable version

		case MARSHAL_HWV_1: {

			// Get vertex list

			ArrayList<SphLatLon> the_vertex_list = SphLatLon.unmarshal_list (reader, "vertex_list");

			// Set up region

			try {
				setup (the_vertex_list);
			}
			catch (Exception e) {
				throw new MarshalException ("SphRegionMercPolygon.do_umarshal: Failed to set up region", e);
			}
		}
		break;

		// Machine-written version

		case MARSHAL_VER_1: {

			// Superclass

			super.do_umarshal (reader);

			// Contents

			vertex_list = SphLatLon.unmarshal_list (reader, "vertex_list");
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
	public SphRegionMercPolygon unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, SphRegionMercPolygon obj) {

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

	public static SphRegionMercPolygon unmarshal_poly (MarshalReader reader, String name) {
		SphRegionMercPolygon result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("SphRegionMercPolygon.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_MERC_POLYGON:
			result = new SphRegionMercPolygon();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

}
