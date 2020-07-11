package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class SubSectProbWithGeomGen {

	public static void main(String[] args) throws IOException {
//		String probColName = "Supra-Seismogenic 30y Participation Prob";
//		File probCSVFile = new File("/tmp/sub_section_probabilities_supra.csv");
//		String probColName = "M>=6.7 30y Participation Prob";
//		File probCSVFile = new File("/tmp/sub_section_probabilities_m6.7.csv");
		String probColName = "M>=7.7 30y Participation Prob";
		File probCSVFile = new File("/tmp/sub_section_probabilities_m7.7.csv");
		File outputFile = new File(probCSVFile.getAbsolutePath().replaceAll(".csv", "")+"_with_geom.csv");
		
		CSVFile<String> probCSV = CSVFile.readFile(probCSVFile, true);
		int probCol = 3; // TD prob
		HashMap<String, Double> probsMap = new HashMap<>();
		
		for (int row=1; row<probCSV.getNumRows(); row++)
			probsMap.put(probCSV.get(row, 0), Double.parseDouble(probCSV.get(row, probCol)));
		
		List<? extends FaultSection> fm31 = RSQSimUtils.getUCERF3SubSectsForComparison(FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		List<? extends FaultSection> fm32 = RSQSimUtils.getUCERF3SubSectsForComparison(FaultModels.FM3_2, DeformationModels.GEOLOGIC);
		
		HashMap<String, FaultSection> sectNameMap = new HashMap<>();
		for (FaultSection sect : fm31)
			sectNameMap.put(sect.getName(), sect);
		for (FaultSection sect : fm32)
			sectNameMap.put(sect.getName(), sect);
		
		List<FaultSection> allSects = new ArrayList<>(sectNameMap.values());
		allSects.sort(new SubSectNameComparator());
		
		CSVFile<String> outCSV = new CSVFile<>(true);
		outCSV.addLine("Parent Section Name", "Parent Section ID", "Sub Section Number",
				"Start Latitude", "Start Longitude", "End Latitude", "End Longitude", probColName);
		
		for (FaultSection sect : allSects) {
			double prob = probsMap.get(sect.getName());
			FaultTrace trace = sect.getFaultTrace();
			String parentName = sect.getParentSectionName();
			int parentID = sect.getParentSectionId();
			int sectNum = getSubsectionNumber(sect.getSectionName());
			for (int i=0; i<trace.size()-1; i++) {
				Location p1 = trace.get(i);
				Location p2 = trace.get(i+1);
				outCSV.addLine(parentName, parentID+"", sectNum+"", (float)p1.getLatitude()+"", (float)p1.getLongitude()+"",
						(float)p2.getLatitude()+"", (float)p2.getLongitude()+"", (float)prob+"");
			}
		}
		
		outCSV.writeToFile(outputFile);
	}
	
	private static final String sect_key = "Subsection ";
	private static final int sect_key_len = sect_key.length();
	
	private static int getSubsectionNumber(String name) {
		int ss1_index = name.indexOf(sect_key)+sect_key_len;
		return Integer.parseInt(name.substring(ss1_index));
	}
	
	public static class SubSectNameComparator implements Comparator<FaultSection> {
		@Override
		public int compare(FaultSection s1, FaultSection s2) {
			String o1 = s1.getName();
			String o2 = s2.getName();
			int ss1_index = o1.indexOf(sect_key)+sect_key_len;
			int ss2_index = o2.indexOf(sect_key)+sect_key_len;
			int ret = o1.substring(0, ss1_index).compareTo(o2.substring(0, ss2_index));
			if (ret != 0)
				return ret;
			// this means same parent section
			int ss1 = getSubsectionNumber(o1);
			int ss2 = getSubsectionNumber(o2);
			Preconditions.checkState(ss1 >= 0 && ss2 >= 0);
			return new Integer(ss1).compareTo(ss2);
		}
		
	}

}
