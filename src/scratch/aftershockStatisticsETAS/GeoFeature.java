package scratch.aftershockStatisticsETAS;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;

import com.google.common.io.Files;

public class GeoFeature {
	// define a class to hold city information for the map
	String type;
	String name;
	Location loc;
	int mapLevel;

	public GeoFeature(String type, String name, double lat, double lon, int mapLevel){
		this(type, name, new Location(lat, lon), mapLevel);
	}

	public GeoFeature(String type, String name, Location loc, int mapLevel){
		this.type = type;
		this.name = name;
		this.loc =loc;
		this.mapLevel = mapLevel;
	}

	public double distanceTo(GeoFeature otherFeature){
		double R = 6371; // kilometres
		double lat0 = loc.getLatRad();
		double lat = otherFeature.loc.getLatRad();
		double lon0 = loc.getLonRad();
		double lon = otherFeature.loc.getLonRad();

		double deltaLat = (lat-lat0);
		double deltaLon = (lon-lon0);

		double a = Math.sin(deltaLat/2) * Math.sin(deltaLat/2) +
				Math.cos(lat0) * Math.cos(lat) *
				Math.sin(deltaLon/2) * Math.sin(deltaLon/2);
		double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));

		return R * c;
	}

	public boolean isInside(Region region){
		return region.contains(this.loc);
	}
	
	@Override
	public String toString(){
		return name + " " + loc.getLatitude() + " " + loc.getLongitude() + " " + mapLevel;
	}

	
		
}
