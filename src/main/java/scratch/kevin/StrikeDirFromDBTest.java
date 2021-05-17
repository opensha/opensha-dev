package scratch.kevin;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.PrefFaultSectionDataDB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

public class StrikeDirFromDBTest {

	private static final DecimalFormat df = new DecimalFormat("000.0000");
	
	public static ArrayList<FaultSectionPrefData> getDBData() {
		
		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
		PrefFaultSectionDataDB_DAO prefDB = new PrefFaultSectionDataDB_DAO(db);
		ArrayList<FaultSectionPrefData> datas = prefDB.getAllFaultSectionPrefData();
		
		try {
			db.destroy();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return datas;
	}
	
	public static double getDiff(double strikeDir, double aveStrike) {
		double diff = Math.abs(strikeDir - aveStrike);
		if (diff < 20)
			return diff;
		if (strikeDir > 180)
			strikeDir -= 360;
		if (aveStrike > 180)
			aveStrike -= 360;
		diff = Math.abs(strikeDir - aveStrike);
		return diff;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int calcTries = 100000;
		
		System.out.println("Loading data...");
		ArrayList<FaultSectionPrefData> datas = getDBData();
		
		Object[][] table = new Object[datas.size()][];
		
		System.out.println("Doing calculations...");
		for (int i=0; i<datas.size(); i++) {
			FaultSectionPrefData data = datas.get(i);
			FaultTrace trace = data.getFaultTrace();
			
			double strikeDir = trace.getStrikeDirection();
			long start = System.currentTimeMillis();
			for (int j=0; j<calcTries; j++) {
				trace.getStrikeDirection();
			}
			long strikeDirMilis = System.currentTimeMillis()-start;
			
			double aveStrike = trace.getAveStrike();
			start = System.currentTimeMillis();
			for (int j=0; j<calcTries; j++) {
				trace.getAveStrike();
			}
			long aveStrikeMilis = System.currentTimeMillis()-start;
			double diff = getDiff(strikeDir, aveStrike);
			
			String name = data.getName();
			Integer pts = trace.getNumLocations();
			Long dirTime = strikeDirMilis;
			Long aveTime = aveStrikeMilis;
			Double dirVal = strikeDir;
			Double aveVal = aveStrike;
			Double deltaVal = diff;
			
			Object[] row = { name, pts, dirTime, aveTime, dirVal, aveVal, deltaVal }; 
			table[i] = row;
		}
		String[] header = { "= '''name''' =", "= '''pts''' =", "= '''dirTime''' =", "= '''aveTime''' =",
							"= '''dirVal''' =", "= '''aveVal''' =", "= '''deltaVal''' =" }; 
		String[] footer = { "'''AVERAGE'''", "pts", "dirTime", "aveTime", "dirVal", "aveVal", "deltaVal" };
		
		double[] aves = new double[footer.length-1];
		
		ArrayList<Object[]> newTable = new ArrayList<Object[]>();
		for (int i=0; i<table.length; i++) {
			if (i % 20 == 0)
				newTable.add(header);
			Object[] tableRow = table[i];
			String[] row = new String[tableRow.length];
			for (int j=0; j<tableRow.length; j++) {
				Object obj = tableRow[j];
				String str;
				if (obj instanceof Double) {
					str = df.format((Double)obj);
				} else {
					str = obj.toString();
				}
				if (j>0) {
					if (obj instanceof Double)
						aves[j-1] += (Double)obj;
					else if (obj instanceof Integer)
						aves[j-1] += (Integer)obj;
					else if (obj instanceof Long)
						aves[j-1] += (Long)obj;
				}
				row[j] = str;
			}
			newTable.add(row);
		}
		if (newTable.get(newTable.size()-1) != header)
			newTable.add(header);
		for (int i=0; i<aves.length; i++) {
			footer[i+1] = df.format(aves[i] / (double)table.length);
		}
		newTable.add(footer);
		
		System.out.println(TracTableCreator.getTracTableString(newTable));
		
		System.exit(0);
	}

}
