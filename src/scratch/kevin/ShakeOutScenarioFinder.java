package scratch.kevin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.ERF2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.AbstractEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

public class ShakeOutScenarioFinder {
	
	public static void findMatches(ProbEqkRupture rup, int sourceID, int rupID,
			Location loc, double targetMag, ERF2DB erf2db, String sourceName) {
		if (rup.getMag() < (targetMag - 0.7) || rup.getMag() > (targetMag + 0.7))
			return;
		boolean closeEnough = false;
		Iterator<Location> it = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface().iterator();
		while (it.hasNext()) {
			Location surfLoc = it.next();
			
			double dist = LocationUtils.linearDistanceFast(loc, surfLoc);
			if (dist < 10d) {
				closeEnough = true;
				break;
			}
		}
		if (!closeEnough)
			return;
		
		HashMap<Integer, Location> hypos = erf2db.getHypocenters(35, sourceID, rupID, 3);
		
		if (hypos.size() == 0)
			return;
		
		double minDist = Double.MAX_VALUE;
		ArrayList<Integer> ids = null;
		Location closest = null;
		
		for (int id : hypos.keySet()) {
			Location hypo = hypos.get(id);
			double dist = LocationUtils.linearDistance(loc, hypo);
			if (dist < minDist) {
				closest = hypo;
				ids = new ArrayList<Integer>();
				minDist = dist;
			}
			if (dist == minDist) {
				ids.add(id);
			}
		}
		
		if (minDist > 15)
			return;
		
		String idString = null;
		for (int id : ids) {
			if (idString == null)
				idString = "";
			else
				idString += ",";
			idString += id;
		}
		
		System.out.println(sourceName + "\tmin hypo dist: " + minDist + "\t src:\t" + sourceID + "\trup: " + rupID + "\tmag: " + rup.getMag() + "\thypo: " + closest + "\tids: " + idString);
	}
	
	public static void main(String[] args) {
		Location loc = new Location(34.02, -117.95, 12.5);
		double mag = 7.1;
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		
		ERF2DB erf2db = new ERF2DB(db);
		
		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = source.getRupture(rupID);
				
				findMatches(rup, sourceID, rupID, loc, mag, erf2db, source.getName());
			}
		}
		System.exit(0);
	}

}
