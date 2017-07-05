package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class TomParsonsSubSectGen {
	
	public static void main(String[] args) throws IOException {
		FaultModels[] fms = { FaultModels.FM2_1, FaultModels.FM3_1, FaultModels.FM3_2 };
		
		DeformationModels[] fm2dms = { DeformationModels.UCERF2_ALL };
		DeformationModels[] fm3dms = { DeformationModels.GEOLOGIC,
				DeformationModels.NEOKINEMA, DeformationModels.ABM, DeformationModels.ZENG };
		
		for (FaultModels fm : fms) {
			DeformationModels[] dms;
			if (fm == FaultModels.FM2_1)
				dms = fm2dms;
			else
				dms = fm3dms;
			
			for (DeformationModels dm : dms) {
				CSVFile<String> csv = new CSVFile<String>(true);
				
				csv.addLine("ID", "Start Lat", "Start Lon", "End Lat", "End Lon", "Slip Rate (mm/yr)", "Rake");
				
				DeformationModelFetcher fetch = new DeformationModelFetcher(
						fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1d);
				
				ArrayList<FaultSectionPrefData> subSects = fetch.getSubSectionList();
				
				for (FaultSectionPrefData subSect : subSects) {
					FaultTrace trace = subSect.getFaultTrace();
					Location loc1 = trace.get(0);
					Location loc2 = trace.get(trace.size()-1);
					
					List<String> line = Lists.newArrayList();
					
					line.add(subSect.getSectionId()+"");
					line.add(loc1.getLatitude()+"");
					line.add(loc1.getLongitude()+"");
					line.add(loc2.getLatitude()+"");
					line.add(loc2.getLongitude()+"");
					line.add(subSect.getOrigAveSlipRate()+"");
					line.add(subSect.getAveRake()+"");
					
					csv.addLine(line);
				}
				
				String prefix = fm.name()+"_"+dm.name()+"_sub_sects";
				csv.writeToFile(new File("/tmp/"+prefix+".csv"));
			}
		}
	}

}
