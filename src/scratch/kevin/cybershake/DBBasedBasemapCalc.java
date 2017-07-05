package scratch.kevin.cybershake;

import java.sql.SQLException;
import java.util.HashMap;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.geo.Location;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.cybershake.db.AttenRelCurves2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;

public class DBBasedBasemapCalc {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		
		int datasetID = 1;
		int imTypeID = 21;
		
//		boolean isProbAt_IML = false;
//		double level = 0.0004;
		boolean isProbAt_IML = true;
		double level = 0.5;
		
		AttenRelCurves2DB curves2db = new AttenRelCurves2DB(db);
		try {
//			HashMap<Location, Integer> ids = curves2db.getCurveIDs(datasetID, imTypeID);
//			System.out.println("Got "+ids.size()+" IDs!");
//			HashMap<Location, ArbitrarilyDiscretizedFunc> curves =
//				curves2db.fetchCurves(datasetID, imTypeID);
			
			GeoDataSet xyz = curves2db.fetchMap(datasetID, imTypeID, isProbAt_IML, level, true, null);
			System.out.println("Got "+xyz.size()+" values!");
			
			ArbDiscrGeoDataSet.writeXYZFile(xyz, "/tmp/map.xyz");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		db.destroy();
		System.exit(0);
	}

}
