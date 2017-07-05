package scratch.kevin;

import java.sql.SQLException;

import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionData;

public class FaultDBDipDirTest {
	
	private static void removeAllDipDirs(DB_AccessAPI db) throws SQLException {
		// warning this will delete ALL of them!
		String sql = "UPDATE "+FaultSectionVer2_DB_DAO.TABLE_NAME
						+" SET "+FaultSectionVer2_DB_DAO.DIP_DIRECTION+"=NULL";
		db.insertUpdateOrDeleteData(sql);
	}

	/**
	 * @param args
	 * @throws SQLException 
	 */
	public static void main(String[] args) throws SQLException {
		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
		
		removeAllDipDirs(db);
		System.exit(0);
		
		FaultSectionVer2_DB_DAO secs2db = new FaultSectionVer2_DB_DAO(db);
		for (FaultSectionData fault : secs2db.getAllFaultSections()) {
			float dipDir = fault.getDipDirection();
			if (!Float.isNaN(dipDir)) {
				double calculatedDir = fault.getFaultTrace().getDipDirection();
				if (dipDir != (float)calculatedDir) {
					System.out.println(fault.getSectionId() + ". " + fault.getSectionName()
							+ " differs!\tdb: " + dipDir + "\tcalc: " + calculatedDir
							+ "\tdiff: " + Math.abs(dipDir - calculatedDir));
				}
			}
		}
		try {
			db.destroy();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.exit(0);
	}

}
