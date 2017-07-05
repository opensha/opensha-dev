package scratch.kevin.ucerf3;

import java.sql.SQLException;
import java.util.ArrayList;

import org.opensha.commons.data.estimate.Estimate;
import org.opensha.commons.data.estimate.MinMaxPrefEstimate;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.gui.addEdit.faultSection.EditFaultSection;
import org.opensha.refFaultParamDb.vo.EstimateInstances;
import org.opensha.refFaultParamDb.vo.FaultSectionData;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

public class RakeRangeFixer {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		DB_AccessAPI db = DB_ConnectionPool.getLatestReadWriteConn();
		DB_ConnectionPool.authenticateDBConnection(true, false);
		
		PrefFaultSectionDataDB_DAO prefDB = new PrefFaultSectionDataDB_DAO(db);
		FaultSectionVer2_DB_DAO fsDB = new FaultSectionVer2_DB_DAO(db);
		
		ArrayList<FaultSectionPrefData> datas = prefDB.getAllFaultSectionPrefData();
		
		for (FaultSectionPrefData data : datas) {
			double rake = data.getAveRake();
			if (Double.isNaN(rake))
				continue;
			double origRake = rake;
			while (rake > 180)
				rake -=360;
			while (rake < -180)
				rake += 360;
			if (rake != origRake) {
				System.out.println("Updating "+data.getName()+" ("+origRake+" => "+rake+")");
				FaultSectionData sect = fsDB.getFaultSection(data.getSectionId());
//				MinMaxPrefEstimate rakeEst = (MinMaxPrefEstimate)sect.getAveRakeEst().getEstimate();
				MinMaxPrefEstimate rakeEst = new MinMaxPrefEstimate(Double.NaN, Double.NaN, rake, Double.NaN, Double.NaN, Double.NaN);
				sect.setAveRakeEst(new EstimateInstances(rakeEst, EditFaultSection.RAKE_UNITS));
				fsDB.update(sect);
			}
		}
		
		System.out.println("DONE!");
		try {
			db.destroy();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.exit(0);
	}

}
