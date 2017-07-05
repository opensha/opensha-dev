package scratch.kevin.cybershake;

import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;

public class CloseSitesFinder {

	public static void main(String[] args) {
		HazardCurveFetcher fetch = new HazardCurveFetcher(
				Cybershake_OpenSHA_DBApplication.getDB(), 35, 21);
		
		List<CybershakeSite> sites = fetch.getCurveSites();
		
		// remove test site(s)
		for (int i=sites.size(); --i>=0;)
			if (sites.get(i).type_id == CybershakeSite.TYPE_TEST_SITE)
				sites.remove(i);
		
		double[] cutoffs = { 0.5, 1d, 2d, 3d };
		
		double[][] dists = new double[sites.size()][sites.size()];
		
		System.out.println("Calculating dists for "+sites.size()+" sites");
		for (int i=0; i<sites.size(); i++) {
			Location loc1 = sites.get(i).createLocation();
			dists[i][i] = Double.POSITIVE_INFINITY;
			for (int j=i+1; j<sites.size(); j++) {
				Location loc2 = sites.get(j).createLocation();
				double dist = LocationUtils.linearDistanceFast(loc1, loc2);
				dists[i][j] = dist;
				dists[j][i] = dist;
			}
		}
		
		int[] counts = new int[cutoffs.length];
		
		for (int n=0; n<counts.length; n++) {
			for (int i=0; i<sites.size(); i++) {
				for (int j=i+1; j<sites.size(); j++) {
					if (dists[i][j] <= cutoffs[n]) {
						counts[n]++;
						System.out.println(sites.get(i).short_name+" and "+sites.get(j).short_name
								+" are within "+cutoffs[n]+" km ("+dists[i][j]+")");
					}
				}
			}
		}
		
		for (int i=0; i<counts.length; i++)
			System.out.println(cutoffs[i]+" km: "+counts[i]+" pairs");
		
		Cybershake_OpenSHA_DBApplication.getDB().destroy();
	}

}
