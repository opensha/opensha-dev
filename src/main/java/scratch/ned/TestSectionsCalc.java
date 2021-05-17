package scratch.ned;

import java.util.List;

import org.opensha.commons.geo.LocationUtils;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionData;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

public class TestSectionsCalc {
	
	public TestSectionsCalc() {
		
		System.out.println("loading data...");
		// this get's the DB accessor (version 3)
	    DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
	    PrefFaultSectionDataDB_DAO faultSectionDB_DAO = new PrefFaultSectionDataDB_DAO(db);
	    List<FaultSectionPrefData> sections = faultSectionDB_DAO.getAllFaultSectionPrefData();
		
//		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
//		FaultSectionVer2_DB_DAO fs2db = new FaultSectionVer2_DB_DAO(db);
//		List<FaultSectionData> sections = fs2db.getAllFaultSections();
		
		System.out.println("starting search...");
		for(FaultSectionPrefData secData: sections) {
			FaultTrace trace = secData.getFaultTrace();
			for(int i=0; i<trace.size()-1;i++) {
				double dist = LocationUtils.horzDistance(trace.get(i), trace.get(i+1));
				if(dist>7 && (i==0 || i==trace.size()-2)) {
					System.out.println(Math.round(dist)+"\t"+i+"\t"+(trace.size()-2)+"\t"+secData.getName());
				}
			}
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TestSectionsCalc test = new TestSectionsCalc();
		
	}

}
