package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class SubSectMidpointFinder {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		int subSectIndex = 2166;
		FaultModels fm = FaultModels.FM3_2;
		File file = new File("D:\\Documents\\temp\\fm3_2_subSect_midpts.csv");
		
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Sub Section Index", "Sub Section Name", "Parent Section Name", "Parent Section ID",
				"Midpoint Lat", "Midpoint Lon");
		
		ArrayList<FaultSectionPrefData> subSects = new DeformationModelFetcher(
				fm, DeformationModels.GEOLOGIC,
				UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1).getSubSectionList();
		
		for (FaultSectionPrefData subSect : subSects) {
			Location midPt = FaultUtils.resampleTrace(subSect.getFaultTrace(), 11).get(5);
			csv.addLine(subSect.getSectionId()+"", subSect.getName(), subSect.getParentSectionName(),
					subSect.getParentSectionId()+"",midPt.getLatitude()+"", midPt.getLongitude()+"");
		}
		
		csv.writeToFile(file);
//		FaultSectionPrefData subSect = subSects.get(subSectIndex);
//		
//		System.out.println("Name: "+subSect.getName());
//		System.out.println("Parent Name: "+subSect.getParentSectionName());
//		System.out.println("Parent ID: "+subSect.getParentSectionId());
//		Location midPt = FaultUtils.resampleTrace(subSect.getFaultTrace(), 11).get(5);
//		
//		System.out.println("Midpoint: "+midPt.getLatitude()+", "+midPt.getLongitude());
	}

}
