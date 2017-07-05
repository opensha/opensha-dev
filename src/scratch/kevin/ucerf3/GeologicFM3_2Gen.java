package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class GeologicFM3_2Gen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		Map<Integer, DeformationSection> fm3_1 =
				DeformationModelFileParser.load(UCERF3_DataUtils.locateResource(
						"DeformationModels", "geologic_slip_rake_2012_02_21.csv"));
		
		File outputFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/DeformationModels/geologic_slip_rake_fm3pt2_2012_02_27.csv");
		
		ArrayList<FaultSectionPrefData> fm3_2_sects = FaultModels.FM3_1.fetchFaultSections();
		
		ArrayList<DeformationSection> fm3_2 = new ArrayList<DeformationModelFileParser.DeformationSection>();
		
		for (FaultSectionPrefData data : fm3_2_sects) {
			int id = data.getSectionId();
			DeformationSection def = fm3_1.get(id);
			if (def == null) {
				// this is a FM 3.2 only section!
				System.out.println("Handling FM 3.2 only section: "+data.getName());
				def = new DeformationSection(id);
				FaultTrace trace = data.getFaultTrace();
				
				for (int i=1; i<trace.size(); i++) {
					Location loc1 = trace.get(i-1);
					Location loc2 = trace.get(i);
					double slip = data.getOrigAveSlipRate();
					double rake = data.getAveRake();
					
					def.add(loc1, loc2, slip, rake);
				}
			}
			
			fm3_2.add(def);
		}
		
		DeformationModelFileParser.write(fm3_2, outputFile);
		
		System.exit(0);
	}

}
