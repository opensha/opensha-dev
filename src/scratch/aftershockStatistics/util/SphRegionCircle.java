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
 * Region for locating aftershocks -- Circle.
 * Author: Michael Barall 04/10/2018.
 *
 * The AAFS server locates aftershocks by drawing a region surrounding the
 * mainshock.  Earthquakes within that region are considered to be aftershocks.
 * This class represents a circle.
 */
public class SphRegionCircle extends SphRegion {

	//----- Region definition -----

	// Latitude and longitude of the center, in degrees.

	private SphLatLon center;

	// Radius of the circle, in kilometers, must be positive, and no more than 180.0001 degrees.

	private double radius;




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
		return SphLatLon.horzDistance(center, loc) <= radius;
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
		return SphLatLon.horzDistance(center, loc) <= radius;
	}




	//----- Plotting -----

	// This function is called when it is necessary to build plot_border.
	// The subclass should supply a function to build the border.
	// A subclass might choose to build the border in its constructor and
	// do_marshal method, in which case this function can just throw an exception.

	@Override
	protected void make_plot_border () {
		plot_border = new LocationList();

		// Latitude, longitude, and distance in radians

		double lat_rad = center.get_lat_rad();
		double lon_rad = center.get_lon_rad();
		double dist_rad = SphLatLon.distance_to_rad (radius);

		// Determine number of subdivisions
		// (Constant 230.0 chosen so radius = 1000 km corresponds to 36 divisions)

		int divs = (int)Math.round (230.0 * Math.sin(dist_rad));

		if (divs < 36) {
			divs = 36;
		} else if (divs > 150) {
			divs = 150;
		}

		double wedge = TWOPI / ((double)divs);

		// Form the circle by traveling along the great circle in each direction

		for (int i = 0; i < divs; ++i) {
			double az_rad = wedge * ((double)i);
			SphLatLon vertex = SphLatLon.gc_travel_rad (lat_rad, lon_rad, az_rad, dist_rad);
			add_border_point (vertex);
		}

		return;
	}




	//----- Special regions -----

	/**
	 * Return true if this is a circular region on the sphere.
	 * If this function returns true, then the center and radius can be
	 * obtained from getCircleCenter and getCircleRadiusDeg.
	 */
	@Override
	public boolean isCircular() {
		return true;
	}

	/**
	 * If this is a circular region on the sphere, then get the center.
	 * Otherwise, throw an exception.
	 */
	@Override
	public SphLatLon getCircleCenter() {
		return center;
	}

	/**
	 * If this is a circular region on the sphere, then get the radius in degrees.
	 * The returned value ranges from 0 to +180.
	 */
	@Override
	public double getCircleRadiusDeg() {
		return Math.min (180.0, SphLatLon.distance_to_deg (radius));
	}




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public SphRegionCircle () {}


	/**
	 * Construct from given center location and radius.
	 * @param center = Center of the circle.
	 * @param radius = Radius of the circle, in kilometers, must be positive and span no more than 180 degrees.
	 */
	public SphRegionCircle (SphLatLon center, double radius) {
		setup (center, radius);
	}


	/**
	 * Set up the region.
	 */
	private void setup (SphLatLon the_center, double the_radius) {

		// Save parameters

		center = the_center;
		radius = the_radius;

		// Radius in degrees

		double radius_deg = SphLatLon.distance_to_deg (radius);
		if (!( radius_deg > 0.0 && radius_deg <= 180.0001 )) {
			throw new IllegalArgumentException ("SphRegionCircle: Radius out-of-range: radius = " + radius);
		}

		// Use it to compute minimum and maximum latitude

		min_lat = center.get_lat() - radius_deg;
		max_lat = center.get_lat() + radius_deg;

		// If either exceeds 89 degrees, construct a polar region

		if (min_lat < -89.0 || max_lat > 89.0) {
			min_lat = Math.max (-90.0, min_lat);
			max_lat = Math.min ( 90.0, max_lat);

			// If the center is closer to the date line than the prime meridian, use the 0 to +360 domain

			if (center.get_lat() < -90.0 || center.get_lat() > 90.0) {
				min_lon = 0.0;
				max_lon = 360.0;
				plot_wrap = true;
			}

			// Otherwise use the -180 to +180 domain

			else {
				min_lon = -180.0;
				max_lon =  180.0;
				plot_wrap = false;
			}
		}

		// Otherwise, draw a box around the Circle

		else {

			// This is the half-angle subtended by the circle at the north or south pole.
			// It is derived from the law of sines for spherical triangles.

			double theta = Math.asin(Math.sin(radius_deg * TO_RAD) / Math.sin((90.0 - center.get_lat()) * TO_RAD)) * TO_DEG;

			// Use it to compute the minimum and maximum longitude

			min_lon = center.get_lon() - theta;
			max_lon = center.get_lon() + theta;

			// If it extends below -180, shift it

			if (min_lon < -180.0) {
				min_lon += 360.0;
				max_lon += 360.0;
				plot_wrap = true;
			}

			// If it extends above +180, no need to shift but the plot needs to be wrapped

			else if (max_lon > 180.0) {
				plot_wrap = true;
			}

			// Otherwise it's in the -180 to +180 domain

			else {
				plot_wrap = false;
			}
		
		}

		plot_border = null;

		return;
	}

	// Display our contents

	@Override
	public String toString() {
		return "SphRegionCircle:" + "\n"
		+ "plot_wrap = " + plot_wrap + "\n"
		+ "min_lat = " + min_lat + "\n"
		+ "max_lat = " + max_lat + "\n"
		+ "min_lon = " + min_lon + "\n"
		+ "max_lon = " + max_lon + "\n"
		+ "center_lat = " + center.get_lat() + "\n"
		+ "center_lon = " + center.get_lon() + "\n"
		+ "radius = " + radius;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_HWV_1 = 1;		// human-writeable version
	private static final int MARSHAL_VER_1 = 13001;

	private static final String M_VERSION_NAME = "SphRegionCircle";

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

		center.marshal (writer, "center");
		writer.marshalDouble ("radius", radius);

		return;
	}

	// Unmarshal object, internal.

	@Override
	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME);

		switch (ver) {

		default:
			throw new MarshalException ("SphRegionCircle.do_umarshal: Unknown version code: version = " + ver);
		
		// Human-writeable version

		case MARSHAL_HWV_1: {

			// Get center and radius

			SphLatLon the_center = (new SphLatLon()).unmarshal (reader, "center");;
			double the_radius = reader.unmarshalDouble ("radius");

			// Set up region

			try {
				setup (the_center, the_radius);
			}
			catch (Exception e) {
				throw new MarshalException ("SphRegionCircle.do_umarshal: Failed to set up region", e);
			}
		}
		break;

		// Machine-written version

		case MARSHAL_VER_1: {

			// Superclass

			super.do_umarshal (reader);

			// Contents

			center = (new SphLatLon()).unmarshal (reader, "center");
			radius = reader.unmarshalDouble ("radius");
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
	public SphRegionCircle unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, polymorphic.

	public static void marshal_poly (MarshalWriter writer, String name, SphRegionCircle obj) {

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

	public static SphRegionCircle unmarshal_poly (MarshalReader reader, String name) {
		SphRegionCircle result;

		reader.unmarshalMapBegin (name);
	
		// Switch according to type

		int type = reader.unmarshalInt (M_TYPE_NAME);

		switch (type) {

		default:
			throw new MarshalException ("SphRegionCircle.unmarshal_poly: Unknown class type code: type = " + type);

		case MARSHAL_NULL:
			result = null;
			break;

		case MARSHAL_CIRCLE:
			result = new SphRegionCircle();
			result.do_umarshal (reader);
			break;
		}

		reader.unmarshalMapEnd ();

		return result;
	}

}
