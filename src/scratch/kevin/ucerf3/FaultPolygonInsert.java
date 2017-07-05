package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.RegionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.vo.FaultSectionData;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.refFaultParamDb.vo.FaultSectionSummary;

import scratch.UCERF3.enumTreeBranches.FaultModels;

import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Maps;

public class FaultPolygonInsert {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File file = new File("D:\\Documents\\temp\\Fault_zones_pts_v3.csv");
		CSVFile<String> csv = CSVFile.readFile(file, true);
		
		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadWriteConn();
		DB_ConnectionPool.authenticateDBConnection(true, true);
		
		try {
			FaultSectionVer2_DB_DAO fs2db = new FaultSectionVer2_DB_DAO(db);
			ArrayList<FaultSectionSummary> summaries = fs2db.getAllFaultSectionsSummary();
			
//			File txtFile = new File(file.getParentFile(), "fault_names.txt");
//			File csvFile = new File(file.getParentFile(), "fault_names.csv");
//			CSVFile<String> namesCSV = new CSVFile<String>(true);
//			HashSet<Integer> dones = new HashSet<Integer>();
//			ArrayList<FaultSectionPrefData> datas = new ArrayList<FaultSectionPrefData>();
//			datas.addAll(FaultModels.FM3_1.fetchFaultSections());
//			datas.addAll(FaultModels.FM3_2.fetchFaultSections());
//			FileWriter fw = new FileWriter(txtFile);
//			for (FaultSectionPrefData data : datas) {
//				int id = data.getSectionId();
//				String name = data.getSectionName();
//				if (dones.contains(id))
//					continue;
//				dones.add(id);
//				fw.write(id+"\t\""+name+"\"\n");
//				namesCSV.addLine(id+"", name);
//			}
//			fw.close();
//			namesCSV.writeToFile(csvFile);
			
			HashMap<Integer, String> nameIDMap = Maps.newHashMap();
			
			for (FaultSectionSummary summ : summaries)
				nameIDMap.put(summ.getSectionId(), summ.getSectionName());
			
			HashMap<Integer, LocationList> locLists = Maps.newHashMap();
			
			for (int row=1; row<csv.getNumRows(); row++) {
				String name = csv.get(row, 0);
				int id;
				try {
					id = Integer.parseInt(csv.get(row, 1));
				} catch (Exception e) {
					id = (int)Double.parseDouble(csv.get(row, 1));
				}
				double lon = Double.parseDouble(csv.get(row, 2));
				double lat = Double.parseDouble(csv.get(row, 3));
				Location loc = new Location(lat, lon);
				
				if (!locLists.containsKey(id)) {
					locLists.put(id, new LocationList());
					// now make sure that the name is accurate
					String realName = nameIDMap.get(id);
					if (realName == null)
						System.out.println("UH OH: name is null!");
					else {
						int diff = StringUtils.getLevenshteinDistance(name, realName);
//						if (diff > realName.length()/10) {
						if (diff > 1) {
							System.out.println("WARNING on name match: "+name+" => "
									+realName+" ("+diff+" steps)");
						}
					}
				}
				LocationList locs = locLists.get(id);
				locs.add(loc);
			}
			int nameMismatches = 0;
			for (int id : locLists.keySet()) {
				if (!nameIDMap.containsKey(id)) {
					System.out.println("Unkown ID: "+id);
					nameMismatches++;
				}
//				Preconditions.checkState(nameIDMap.containsKey(name), "Unkown name: "+name);
				LocationList list = locLists.get(id);
				String name = nameIDMap.get(id);
				Preconditions.checkState(list.size()>2, "too few locs for "+name+": "+list.size());
				
				try {
					Region reg = new Region(list, BorderType.MERCATOR_LINEAR);
					System.out.println("Updating polygon for: "+id+". "+name);
					fs2db.updateZonePolygon(id, reg);
				} catch (Exception e) {
					System.out.println("ERROR for: "+name);
					e.printStackTrace();
					RegionUtils.locListToKML(list, "error_"+id+".kml", Color.RED);
				}
			}
			System.out.println("Name mismatches: "+nameMismatches+"/"+locLists.keySet().size());
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				db.destroy();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}
		System.exit(0);
	}

}
