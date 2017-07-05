package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class SubSectCenterCalc {
	
	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_2;
		
		ArrayList<FaultSectionPrefData> subSects = new DeformationModelFetcher(
				fm, DeformationModels.GEOLOGIC, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1)
				.getSubSectionList();
		
		CSVFile<String> output = new CSVFile<String>(true);
		
		output.addLine("ID", "Name", "Center Lat", "Center Lon");
		
		for (FaultSectionPrefData subSect : subSects) {
			FaultTrace trace = subSect.getFaultTrace();
			trace = FaultUtils.resampleTrace(trace, 101);
			Location center = trace.get(50);
			
			output.addLine(subSect.getSectionId()+"", subSect.getName(),
					center.getLatitude()+"", center.getLongitude()+"");
		}
		
		output.writeToFile(new File("/tmp/fm3_2_subsect_centers.csv"));
	}

}
