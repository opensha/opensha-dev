package scratch.aftershockStatistics.util;

import java.util.List;
import java.util.ArrayList;

import org.opensha.commons.geo.Location;

import static java.lang.Math.PI;
import static org.opensha.commons.geo.GeoTools.TWOPI;
import static org.opensha.commons.geo.GeoTools.TO_DEG;
import static org.opensha.commons.geo.GeoTools.TO_RAD;
import static org.opensha.commons.geo.GeoTools.EARTH_RADIUS_MEAN;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;


/**
 * Class to hold a (lat, lon) pair in degrees, denoting a point on the sphere.
 * Author: Michael Barall 04/13/2018.
 *
 * This is an immutable object.
 */
public class SphLatLon {

	//----- Contents -----

	// Latitude, in degress, in range -180 to +180.

	private double lat;

	// Longitude, in degress, in range -90 to +90.

	private double lon;




	//----- Construction -----

	// Set the latitude and longitude, return true if success, false if out-of-range.
	// Values slightly outside the range are coerced to be in range.

	private boolean set (double the_lat, double the_lon) {
		lat = the_lat;
		lon = the_lon;

		if (lat < -90.0) {
			if (lat < -90.0001) {
				return false;
			}
			lat = -90.0;
		} else if (lat > 90.0) {
			if (lat > 90.0001) {
				return false;
			}
			lat = 90.0;
		}

		if (lon < -180.0) {
			if (lon < -180.0001) {
				return false;
			}
			lon = -180.0;
		} else if (lon > 180.0) {
			if (lon > 180.0001) {
				return false;
			}
			lon = 180.0;
		}

		return true;
	}

	// Default constructor.

	public SphLatLon () {}

	// Construct from given latitude and longitude.

	public SphLatLon (double lat, double lon) {
		if (!( set (lat, lon) )) {
			throw new IllegalArgumentException ("SphLatLon: Coordinates out-of-range: lat = " + lat + ", lon = " + lon);
		}
	}

	// Construct from a Location.

	public SphLatLon (Location loc) {
		double the_lat = loc.getLatitude();
		double the_lon = loc.getLongitude();

		// Location permits longitudes -180 to +360, so do the wrapping
		if (the_lon > 180.0) {
			the_lon -= 360.0;
		}

		if (!( set (the_lat, the_lon) )) {
			throw new IllegalArgumentException ("SphLatLon: Coordinates out-of-range: lat = " + the_lat + ", lon = " + the_lon);
		}
	}




	//----- Getters -----

	// Get latitude.

	public double get_lat () {
		return lat;
	}

	// Get longitude.

	public double get_lon () {
		return lon;
	}

	// Get longitude.
	// If f_wrap_lon is true, the result lies between 0 and +360.
	// If f_wrap_lon is false, the result lies between -180 and +180.

	public double get_lon (boolean f_wrap_lon) {
		return (f_wrap_lon && lon < 0.0) ? (lon + 360.0) : lon;
	}

	// Get latitude, in radians.

	public double get_lat_rad () {
		return lat * TO_RAD;
	}

	// Get longitude, in radians.

	public double get_lon_rad () {
		return lon * TO_RAD;
	}

	// Get longitude, in radians.
	// If f_wrap_lon is true, the result lies between 0 and +360.
	// If f_wrap_lon is false, the result lies between -180 and +180.

	public double get_lon_rad (boolean f_wrap_lon) {
		return ((f_wrap_lon && lon < 0.0) ? (lon + 360.0) : lon) * TO_RAD;
	}




	//----- Services for Location -----

	// Get location.

	public Location get_location () {
		return new Location (get_lat(), get_lon());
	}

	// Get location.
	// If f_wrap_lon is true, the longitude lies between 0 and +360.
	// If f_wrap_lon is false, the longitude lies between -180 and +180.

	public Location get_location (boolean f_wrap_lon) {
		return new Location (get_lat(), get_lon(f_wrap_lon));
	}

	// Get longitude from a Location.
	// If f_wrap_lon is true, the result lies between 0 and +360.
	// If f_wrap_lon is false, the result lies between -180 and +180.

	public static double get_lon (Location loc, boolean f_wrap_lon) {
		double lon = loc.getLongitude();
		if (f_wrap_lon) {
			if (lon < 0.0) {
				lon += 360.0;
			}
		}
		else {
			if (lon > 180.0) {
				lon -= 360.0;
			}
		}
		return lon;
	}




	//----- Spherical geometry -----

	/**
	 * Calculates the angle between two points using the haversine formula.
	 * http://en.wikipedia.org/wiki/Haversine_formula
	 * This method properly handles values spanning 180 degrees longitude.
	 * See http://www.edwilliams.org/avform.htm Aviation Formulary for source.
	 * Result is returned in radians.
	 * Note: Same code as in LocationUtils.java.
	 * 
	 * @param p1 the first point
	 * @param p2 the second point
	 * @return the angle between the points (in radians)
	 */
	public static double angle_rad (SphLatLon p1, SphLatLon p2) {
		double lat1 = p1.get_lat_rad();
		double lat2 = p2.get_lat_rad();
		double sinDlatBy2 = Math.sin((lat2 - lat1) / 2.0);
		double sinDlonBy2 = Math.sin((p2.get_lon_rad() - p1.get_lon_rad()) / 2.0);
		// half length of chord connecting points
		double c = (sinDlatBy2 * sinDlatBy2) +
			(Math.cos(lat1) * Math.cos(lat2) * sinDlonBy2 * sinDlonBy2);
		return 2.0 * Math.atan2(Math.sqrt(c), Math.sqrt(1 - c));
	}

	/**
	 * Version that accepts one SphLatLon and one Location.
	 */
	public static double angle_rad (SphLatLon p1, Location p2) {
		double lat1 = p1.get_lat_rad();
		double lat2 = p2.getLatRad();
		double sinDlatBy2 = Math.sin((lat2 - lat1) / 2.0);
		double sinDlonBy2 = Math.sin((p2.getLonRad() - p1.get_lon_rad()) / 2.0);
		// half length of chord connecting points
		double c = (sinDlatBy2 * sinDlatBy2) +
			(Math.cos(lat1) * Math.cos(lat2) * sinDlonBy2 * sinDlonBy2);
		return 2.0 * Math.atan2(Math.sqrt(c), Math.sqrt(1 - c));
	}

	/**
	 * Calculates the great circle surface distance between two
	 * points using the haversine formula for computing the
	 * angle between two points.
	 * Note: Same code as in LocationUtils.java.
	 * 
	 * @param p1 the first point
	 * @param p2 the second point
	 * @return the distance between the points in km
	 */
	public static double horzDistance (SphLatLon p1, SphLatLon p2) {
		return EARTH_RADIUS_MEAN * angle_rad (p1, p2);
	}

	/**
	 * Version that accepts one SphLatLon and one Location.
	 */
	public static double horzDistance (SphLatLon p1, Location p2) {
		return EARTH_RADIUS_MEAN * angle_rad (p1, p2);
	}

	/**
	 * Convert horizontal distance to radians.
	 */
	public static double distance_to_rad (double distance) {
		return distance / EARTH_RADIUS_MEAN;
	}

	/**
	 * Convert horizontal distance to degrees.
	 */
	public static double distance_to_deg (double distance) {
		return (distance / EARTH_RADIUS_MEAN) * TO_DEG;
	}

	/**
	 * Calculates the endpoint when traveling along a great circle,
	 * with a given starting point, azimuth, and angular distance traveled.
	 * Note: Derived from code in LocationUtils.java.
	 * 
	 * @param lat_rad = Starting point latitude, in radians.
	 * @param lon_rad = Starting point longitude, in radians.
	 * @param az_rad = Initial azimuth, in radians; 0 = north, pi/2 = east.
	 * @param dist_rad = Angular distance traveled, in radians.
	 * @return
	 * Returns the ending point.
	 */
	public static SphLatLon gc_travel_rad (double lat_rad, double lon_rad, double az_rad, double dist_rad) {

		double sinLat1 = Math.sin(lat_rad);
		double cosLat1 = Math.cos(lat_rad);
		double sinD = Math.sin(dist_rad);
		double cosD = Math.cos(dist_rad);

		double lat2 = Math.asin(sinLat1 * cosD + cosLat1 * sinD * Math.cos(az_rad));

		double lon2 = lon_rad +
			Math.atan2(Math.sin(az_rad) * sinD * cosLat1,
				cosD - sinLat1 * Math.sin(lat2));

		lat2 = lat2 * TO_DEG;
		lon2 = lon2 * TO_DEG;

		while (lon2 < -180.0) {
			lon2 += 360.0;
		}
		while (lon2 > 180.0) {
			lon2 -= 360.0;
		}

		return new SphLatLon (lat2, lon2);
	}

	/**
	 * Calculates the endpoint when traveling along a great circle,
	 * with a given starting point, azimuth, and kilometer distance traveled.
	 * Note: Derived from code in LocationUtils.java.
	 * 
	 * @param lat_rad = Starting point latitude, in radians.
	 * @param lon_rad = Starting point longitude, in radians.
	 * @param az_rad = Initial azimuth, in radians; 0 = north, pi/2 = east.
	 * @param dist_km = Distance traveled, in kilometers.
	 * @return
	 * Returns the ending point.
	 */
	public static SphLatLon gc_travel_km (double lat_rad, double lon_rad, double az_rad, double dist_km) {
		return gc_travel_rad (lat_rad, lon_rad, az_rad, dist_km / EARTH_RADIUS_MEAN);
	}




	//----- Marshaling -----

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		writer.marshalDouble ("lat", lat);
		writer.marshalDouble ("lon", lon);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public SphLatLon unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		double the_lat = reader.unmarshalDouble ("lat");
		double the_lon = reader.unmarshalDouble ("lon");
		if (!( set (the_lat, the_lon) )) {
			throw new MarshalException ("SphLatLon: Coordinates out-of-range: lat = " + the_lat + ", lon = " + the_lon);
		}
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal a list of objects.

	public static void marshal_list (MarshalWriter writer, String name, List<SphLatLon> obj_list) {
		int n = obj_list.size();
		writer.marshalArrayBegin (name, n);
		for (SphLatLon obj : obj_list) {
			obj.marshal (writer, null);
		}
		writer.marshalArrayEnd ();
		return;
	}

	// Unmarshal a list of objects.

	public static ArrayList<SphLatLon> unmarshal_list (MarshalReader reader, String name) {
		ArrayList<SphLatLon> obj_list = new ArrayList<SphLatLon>();
		int n = reader.unmarshalArrayBegin (name);
		for (int i = 0; i < n; ++i) {
			obj_list.add ((new SphLatLon()).unmarshal (reader, null));
		}
		reader.unmarshalArrayEnd ();
		return obj_list;
	}
}
