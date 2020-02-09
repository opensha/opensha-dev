package scratch.kevin.simulators.ruptures.azimuthal;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;

public class SimpleGC2SiteLocationCalc {
	
	private static final double half_pi = 0.5*Math.PI;
	
	public static Location gc2ToLoc(Location p1, Location p2, double rx, double ry) {
		double strike = LocationUtils.azimuthRad(p1, p2);
		
		// first move along the x axis
		Location xPt = LocationUtils.location(p1, strike+half_pi, rx);
		// now move along the y axis
		return LocationUtils.location(xPt, strike, ry);
	}
	
	public static double[] locToGC2(Location p1, Location p2, Location siteLoc) {
		double rx = LocationUtils.distanceToLineFast(p1, p2, siteLoc);
		double strike = LocationUtils.azimuthRad(p1, p2);
		// create a point that is left in the reference frame of the rupture
		// thus for a line from p1 to this point, values on the right will have a positive
		// distance to line (rY)
		Location perpLeftLoc = LocationUtils.location(p1, strike-half_pi, 10d);
		double ry = LocationUtils.distanceToLineFast(p1, perpLeftLoc, siteLoc);
		
		return new double[] { rx, ry };
	}

	public static void main(String[] args) {
		Location p1 = new Location(34d, -118d);
		Location p2 = new Location(35d, -118d);
		
		Location siteLoc = new Location(35d, -117d);
		double[] gc2 = locToGC2(p1, p2, siteLoc);
		System.out.println("rX="+gc2[0]+"\trY="+gc2[1]);
		Location recoveredSite = gc2ToLoc(p1, p2, gc2[0], gc2[1]);
		System.out.println("original site: "+siteLoc);
		System.out.println("recovered site: "+recoveredSite);
		double[] gc2_2 = locToGC2(p1, p2, recoveredSite);
		System.out.println("recovered rX="+gc2_2[0]+"\trY="+gc2_2[1]);
		
	}

}
