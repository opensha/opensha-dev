package scratch.kevin.cybershake;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.InterpolatedEvenlyGriddedSurface;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

public class BatchSiteDistUpdate {

	public static void main(String[] args) throws IOException {
		DBAccess db = Cybershake_OpenSHA_DBApplication.getAuthenticatedDBAccess(
				true, true, Cybershake_OpenSHA_DBApplication.PRODUCTION_HOST_NAME);
		
//		int erfID = 36;
//		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2_200mERF(true);
		int erfID = 41;
		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2_200mERF(false);
		erf.updateForecast();
		
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
		List<CybershakeSite> sites = sites2db.getAllSitesFromDB();
		
		siteLoop:
		for (CybershakeSite site : sites) {
			String select = "SELECT Source_ID,Rupture_ID FROM CyberShake_Site_Ruptures WHERE CS_Site_ID="+site.id
					+" AND ERF_ID="+erfID+" AND Site_Rupture_Dist IS NULL";
			
			System.out.println("Handling "+site.id+". "+site.short_name);
			Location loc = site.createLocation();
			
			List<DistCalcTask> tasks = Lists.newArrayList();
			try {
				ResultSet rs = db.selectData(select);
				
				if (!rs.first()) {
					System.out.println("\talready done.");
					continue;
				}
				
				while (!rs.isAfterLast()) {
					int sourceID = rs.getInt("Source_ID");	
					int rupID = rs.getInt("Rupture_ID");
					
					ProbEqkSource source = erf.getSource(sourceID);
					int numRups = source.getNumRuptures();
					Preconditions.checkState(numRups > rupID, "Bad rupture ID for source %s, have %s but ID %s encountered. Source: %s",
							sourceID, numRups, rupID, source.getName());
					
					tasks.add(new DistCalcTask(erf, sourceID, rupID, loc));
					
					rs.next();
				}
				rs.close();
			} catch (Exception e) {
				e.printStackTrace();
				break;
			}
			
			if (tasks.isEmpty()) {
				System.out.println("\talready done.");
				continue;
			}
			
			System.out.println("\tcalculating "+tasks.size()+" missing distances");
			
			Stopwatch watch = Stopwatch.createStarted();
			try {
				new ThreadedTaskComputer(tasks, true).computeThreaded();
			} catch (InterruptedException e) {
				e.printStackTrace();
				break;
			}
			watch.stop();
			long secs = watch.elapsed(TimeUnit.SECONDS);
			System.out.println("\tcalc took "+secs+" seconds");
			
			MinMaxAveTracker distTrack = new MinMaxAveTracker();
			
			for (DistCalcTask task : tasks) {
				distTrack.addValue(task.minDist);
				
				String sql = "UPDATE CyberShake_Site_Ruptures SET Site_Rupture_Dist="+task.minDist+" WHERE CS_Site_ID="+site.id
					+" AND ERF_ID="+erfID+" AND Source_ID="+task.sourceID+" AND Rupture_ID="+task.rupID;
				
				try {
					db.insertUpdateOrDeleteData(sql);
				} catch (SQLException e) {
					System.out.println("Error on statement: "+sql);
					e.printStackTrace();
					break siteLoop;
				}
			}
			
			System.out.println("\tdone inserting. distance range: "+distTrack);
		}
		
		db.destroy();
		
		System.exit(0);
	}
	
	private static class DistCalcTask implements Task {
		
		private AbstractERF erf;
		private int sourceID, rupID;
		private Location loc;
		
		private double minDist = Double.NaN;

		public DistCalcTask(AbstractERF erf, int sourceID, int rupID, Location loc) {
			super();
			this.erf = erf;
			this.sourceID = sourceID;
			this.rupID = rupID;
			this.loc = loc;
		}

		@Override
		public void compute() {
			ProbEqkRupture rup = erf.getRupture(sourceID, rupID);
			
			EvenlyGriddedSurface surf = (EvenlyGriddedSurface) rup.getRuptureSurface();
			if (surf instanceof InterpolatedEvenlyGriddedSurface)
				surf = ((InterpolatedEvenlyGriddedSurface)surf).getLowResSurface();
			
//			int mod = 1;
//			if (surf.getAveGridSpacing() < 0.6)
//				mod = (int)Math.floor(0.6/surf.getAveGridSpacing());
//			mod = 1;
//			
//			System.out.println("Mod: "+mod);
			
			minDist = Double.POSITIVE_INFINITY;
			
			for (int row=0; row<surf.getNumRows(); row++) {
//				if (row % mod != 0)
//					continue;
				for (int col=0; col<surf.getNumCols(); col++) {
//					if (col % mod != 0)
//						continue;
					minDist = Math.min(minDist, LocationUtils.linearDistanceFast(loc, surf.get(row, col)));
				}
			}
//			System.out.println("Min dist: "+minDist);
		}
		
	}

}