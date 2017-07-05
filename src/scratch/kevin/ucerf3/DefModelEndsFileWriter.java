package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.collect.Lists;

import scratch.UCERF3.analysis.DeformationModelsCalc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class DefModelEndsFileWriter {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultModels[] fms = {FaultModels.FM3_1, FaultModels.FM3_2};
		DeformationModels dm = DeformationModels.GEOLOGIC;
		
		int[] maxAways = { 0, 1 };
		
		File dir = new File("/tmp");
		
		for (FaultModels fm : fms) {
			DeformationModelFetcher fetch = new DeformationModelFetcher(fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR,
					InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE);
			
			Map<IDPairing, Double> distances = fetch.getSubSectionDistanceMap(5d);
			
				for (int maxAwayFromEnd : maxAways) {
				
				List<FaultSectionPrefData> isolated =
						DeformationModelsCalc.getIsolatedEndpoints(fetch.getSubSectionList(), distances, maxAwayFromEnd);
				
				File file = new File(dir, fm.getShortName()+"_isolated_"+maxAwayFromEnd+"_from_end.csv");
				
				CSVFile<String> csv = new CSVFile<String>(true);
				
				csv.addLine(Lists.newArrayList(
						"Sub Section ID", "Parent Section ID", "Sub Section Name", "Endpoint Lat", "Endpoint Lon"));
				for (FaultSectionPrefData sect : isolated) {
					String[] split = sect.getName().split(" ");
					int subSectID = Integer.parseInt(split[split.length-1]);
					FaultTrace trace = sect.getFaultTrace();
					Location loc;
					if (subSectID == 0)
						loc = trace.get(0);
					else
						loc = trace.get(trace.size()-1);
					
					csv.addLine(Lists.newArrayList(sect.getSectionId()+"", sect.getParentSectionId()+"",
							sect.getName(), (float)loc.getLatitude()+"", (float)loc.getLongitude()+""));
				}
				
				csv.writeToFile(file);
				csv.writeToTabSeparatedFile(new File(dir, file.getName().replaceAll(".csv", ".txt")), 1);
			}
		}
	}

}
