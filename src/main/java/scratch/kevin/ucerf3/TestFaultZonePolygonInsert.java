package scratch.kevin.ucerf3;

import java.sql.SQLException;

import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionData;

public class TestFaultZonePolygonInsert {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			DB_AccessAPI db = DB_ConnectionPool.getDB3ReadWriteConn();
			DB_ConnectionPool.authenticateDBConnection(false, false);
			
			FaultSectionVer2_DB_DAO fs2db = new FaultSectionVer2_DB_DAO(db);
			
			FaultSectionData fault = fs2db.getFaultSection(861);
			System.out.println(fault.getZonePolygon());
			
//			PrefFaultSectionDataDB_DAO pref2db = new PrefFaultSectionDataDB_DAO(db);
//			
//			pref2db.getAllFaultSectionPrefData();
			
//			FaultSectionData orig = fs2db.getFaultSection(114);
//			FaultSectionData newData = orig.clone();
//			newData.setSectionId(-1);
//			newData.setSectionName("Kevin Test Section!");
//			newData.setZonePolygon(new Region(newData.getFaultTrace(), 10d));
//			fs2db.removeFaultSection(942);
			
//			fs2db.addFaultSection(newData);
			
			db.destroy();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
