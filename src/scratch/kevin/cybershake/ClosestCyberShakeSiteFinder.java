package scratch.kevin.cybershake;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;

public class ClosestCyberShakeSiteFinder {

	public static void main(String[] args) {
		Location loc = new Location(34.04288, -118.36893);
		double maxDist = 2d; // km
		
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(Cybershake_OpenSHA_DBApplication.getDB());
		
		for (CybershakeSite site : sites2db.getAllSitesFromDB()) {
			Location testLoc = site.createLocation();
			double dist = LocationUtils.horzDistanceFast(loc, testLoc);
			if (dist < maxDist)
				System.out.println((float)dist+" km: "+site);
		}
	}

}
