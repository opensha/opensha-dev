package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultModelDB_DAO;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

public class AseismicityFileWriter {
	
	private static void write(DB_AccessAPI db, File file, int fmID) throws IOException {
		PrefFaultSectionDataDB_DAO fs2db = new PrefFaultSectionDataDB_DAO(db);
		FaultModelDB_DAO fm2db = new FaultModelDB_DAO(db);
		
		ArrayList<Integer> fm = fm2db.getFaultSectionIdList(fmID);
		
		ArrayList<FaultSectionPrefData> sects = fs2db.getAllFaultSectionPrefData();
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine(Lists.newArrayList("ID", "Name", "Aseismicity Factor"));
		
		for (FaultSectionPrefData sect : sects) {
			if (!fm.contains(sect.getSectionId()))
				continue;
			
			csv.addLine(Lists.newArrayList(sect.getSectionId()+"", sect.getSectionName(), sect.getAseismicSlipFactor()+""));
		}
		
		csv.writeToFile(file);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
		File outputFile = new File("/tmp/ucerf3_aseis.csv");
		int fmID = 101;
		
//		DB_AccessAPI db = DB_ConnectionPool.getDB2ReadOnlyConn();
//		File outputFile = new File("/tmp/ucerf2_aseis.csv");
//		int fmID = 41;
		
		try {
			write(db, outputFile, fmID);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			try {
				db.destroy();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		System.exit(0);
	}

}
