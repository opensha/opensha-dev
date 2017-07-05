package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class ParkfieldAseisFileGen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		
		FaultModels fm = FaultModels.FM3_1;
		
		List<DeformationModels> dms = DeformationModels.forFaultModel(fm);
		dms.remove(dms.indexOf(DeformationModels.GEOLOGIC_UPPER));
		dms.remove(dms.indexOf(DeformationModels.GEOLOGIC_LOWER));
		
		ArrayList<String> header = Lists.newArrayList("Section Index", "Section Name", "Length (KM)");
		for (DeformationModels dm : dms) {
			header.add(dm.getShortName()+" Orig Slip");
			header.add(dm.getShortName()+" Mo Rate Reduction");
		}
		
		csv.addLine(header);
		
		ArrayList<FaultSectionPrefData> sects = fm.fetchFaultSections();
		
		ArrayList<Integer> parentSects = Lists.newArrayList(657, 658, 32, 285, 300);
		
		List<Map<Integer, DeformationSection>> dmDatas = Lists.newArrayList();
		
		for (DeformationModels dm : dms) {
			Map<Integer, DeformationSection> def = DeformationModelFileParser.load(dm.getDataFileURL(fm));
			DeformationModelFileParser.applyMomentReductions(def, Double.POSITIVE_INFINITY);
			dmDatas.add(def);
		}
		
		for (int sectID : parentSects) {
			String name = null;
			FaultTrace trace = null;
			for (FaultSectionPrefData data : sects) {
				if (data.getSectionId() == sectID) {
					name = data.getSectionName();
					trace = data.getFaultTrace();
				}
			}
			int numSubs = dmDatas.get(0).get(sectID).getSlips().size();
			
			for (int i=0; i<numSubs; i++) {
				ArrayList<String> line = Lists.newArrayList();
				int[] mini = {sectID, i+1};
				line.add(DeformationModelFileParser.getMinisectionString(mini));
				line.add(name+" (mini sect "+(i+1)+")");
				double length = LocationUtils.horzDistance(trace.get(i), trace.get(i+1));
				line.add(length+"");
				for (int j=0; j<dms.size(); j++) {
					Map<Integer, DeformationSection> dmData = dmDatas.get(j);
					DeformationSection def = dmData.get(sectID);
					double slip = def.getSlips().get(i);
					double momRed = def.getMomentReductions().get(i);
					
					line.add(slip+"");
					line.add(momRed+"");
				}
				
				csv.addLine(line);
			}
		}
		
		csv.writeToFile(new File("/tmp/parkfield_mom_reds.csv"));
		csv.writeToTabSeparatedFile(new File("/tmp/parkfield_mom_reds.txt"), 1);
	}

}
