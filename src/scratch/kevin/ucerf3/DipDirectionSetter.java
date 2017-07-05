package scratch.kevin.ucerf3;

import java.sql.SQLException;
import java.util.ArrayList;

import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionData;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

public class DipDirectionSetter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadWriteConn();
		DB_ConnectionPool.authenticateDBConnection(true, true);
		
		try {
			PrefFaultSectionDataDB_DAO pref2db = new PrefFaultSectionDataDB_DAO(db);
			ArrayList<FaultSectionPrefData> prefData = pref2db.getAllFaultSectionPrefData();
//			for (FaultSectionPrefData data : pref2db.getAllFaultSectionPrefData()) {
//				double dir = data.getDipDirection();
//				System.out.println("DIP DIRECTION: "+dir);
////				if (dir != 0 && !Double.isNaN(dir)) {
////					System.out.println("CUSTOM DIP DIRECTION: "+dir);
////				}
//			}
			
			FaultSectionVer2_DB_DAO fs2db = new FaultSectionVer2_DB_DAO(db);
//			System.out.println("getting sections...");
//			ArrayList<FaultSectionData> sects = fs2db.getAllFaultSections();
			System.out.println("setting dip dirs!");
			for (FaultSectionPrefData data : prefData) {
				double newDir = data.getFaultTrace().getDipDirection();
				System.out.println("setting dip dir. "+data.getSectionName()+": "+(float)newDir);
				fs2db.updateDipDirection(data.getSectionId(), (float)newDir);
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				db.destroy();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}
		
		System.out.println("DONE.");
		
		System.exit(0);
	}

}
