package scratch.kevin.cybershake;

import java.io.File;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM_Vs30;
import org.opensha.commons.data.siteData.impl.CVM_Vs30.CVM;
import org.opensha.commons.geo.Location;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.SiteInfo2DB;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;

public class NicoRecordFetch {

	public static void main(String[] args) throws SQLException, IOException {
		int datasetID = 61; // study 15.12
		int erfID = 36;
		SiteData<Double> vs30Prov = new CVM_Vs30(CVM.CVMS4i26);
		CyberShakeComponent comp = CyberShakeComponent.GEOM_MEAN;
		
		int numPerSite = 3;
		int refIM = 21;
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB(Cybershake_OpenSHA_DBApplication.ARCHIVE_HOST_NAME);
		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		SiteInfo2DB siteDB = new SiteInfo2DB(db);
		
		HazardCurveFetcher fetcher = new HazardCurveFetcher(db, datasetID, refIM);
		
		List<Integer> runIDs = fetcher.getRunIDs();
		List<CybershakeSite> sites = fetcher.getCurveSites();
		
		Random r = new Random();
		
		HashSet<Integer> imTypeIDs = new HashSet<Integer>();
		
		List<Result> results = Lists.newArrayList();
		
		for (int i=0; i<sites.size(); i++) {
			CybershakeSite site = sites.get(i);
			System.out.println("Site: "+site);
			if (site.type_id == CybershakeSite.TYPE_TEST_SITE)
				continue;
			List<Integer> srcIdList = siteDB.getSrcIdsForSite(site.id, erfID);
			
			Location loc = site.createLocation();
			double vs30 = vs30Prov.getValue(loc);
			
			for (int n=0; n<numPerSite; n++) {
				int sourceID = srcIdList.get(r.nextInt(srcIdList.size()));
				int rupID = r.nextInt(erf.getNumRuptures(sourceID));

				ProbEqkRupture rup = erf.getRupture(sourceID, rupID);

				double dist = rup.getRuptureSurface().getDistanceRup(loc);
				
				String sql = "SELECT Rup_Var_ID,IM_Type_ID,IM_Value FROM PeakAmplitudes "
						+ "WHERE Source_ID="+sourceID+" AND Rupture_ID="+rupID+" AND Run_ID="+runIDs.get(i);
				ResultSet rs = db.selectData(sql);
				
				// RV, IMT, IM
				Table<Integer, Integer, Double> siteResults = HashBasedTable.create();
				
				int maxRV = 0;
				
				rs.first();
				while(!rs.isAfterLast()){
					try {
						int rvID = rs.getInt("Rup_Var_ID");
						if (rvID > maxRV)
							maxRV = rvID;
						int imt = rs.getInt("IM_Type_ID");
						imTypeIDs.add(imt);
						double im = rs.getDouble("IM_Value");
						siteResults.put(rvID, imt, im);
						rs.next();
					} catch (Exception e) {
						System.out.println("No values!");
						break;
					}
				}
				if (siteResults.isEmpty())
					continue;
				
				int rvID = r.nextInt(maxRV+1);
				Map<Integer, Double> rvResults = siteResults.row(rvID);
				results.add(new Result(site.id, runIDs.get(i), sourceID, rupID, rvID, dist, rup.getMag(), vs30, rvResults));
				
				rs.close();
			}
		}
		
		HazardCurve2DB curves2db = new HazardCurve2DB(db);
		List<CybershakeIM> ims = Lists.newArrayList();
		for (int imTypeID : imTypeIDs) {
			CybershakeIM im;
			try {
				im = curves2db.getIMFromID(imTypeID);
			} catch (Exception e) {
				continue;
			}
			if (im.getComponent() == comp && im.getMeasure() == IMType.SA)
				ims.add(im);
		}
		Collections.sort(ims, new Comparator<CybershakeIM>() {

			@Override
			public int compare(CybershakeIM o1, CybershakeIM o2) {
				return Double.compare(o1.getVal(), o2.getVal());
			}
		});
		
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Site ID", "Run ID", "Source ID", "Rup ID", "RV ID", "RRup", "Mag", "Vs30");
		for (CybershakeIM im : ims)
			header.add((float)im.getVal()+"s");
		csv.addLine(header);
		
		for (Result rs : results) {
			List<String> line = Lists.newArrayList(rs.siteID+"", rs.runID+"", rs.sourceID+"", rs.rupID+"", rs.rvID+"", rs.dist+"",
					rs.mag+"", rs.vs30+"");
			for (CybershakeIM im : ims) {
				double val = rs.imMap.get(im.getID());
				line.add(val/HazardCurveComputation.CONVERSION_TO_G+"");
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/cs_amps_for_nico.csv"));
		
		db.destroy();
	}
	
	private static class Result {
		private double dist, mag, vs30;
		private Map<Integer, Double> imMap;
		private int siteID;
		private int runID;
		private int sourceID;
		private int rupID;
		private int rvID;
		
		public Result(int siteID, int runID, int sourceID, int rupID, int rvID, double dist, double mag, double vs30,
				Map<Integer, Double> imMap) {
			super();
			this.siteID = siteID;
			this.runID = runID;
			this.sourceID = sourceID;
			this.rupID = rupID;
			this.rvID = rvID;
			this.dist = dist;
			this.mag = mag;
			this.vs30 = vs30;
			this.imMap = imMap;
		}
	}

}
