package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class SpecialFaultSlipRates {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_1;
		
		int[] ids = {184, 196, 240};
		Map<DeformationModels, double[]> slips = Maps.newHashMap();
		
		List<FaultSectionPrefData> sects = fm.fetchFaultSections();
		Map<Integer, FaultSectionPrefData> sectsMap = Maps.newHashMap();
		for (FaultSectionPrefData sect : sects)
			sectsMap.put(sect.getSectionId(), sect);
		
		for (DeformationModels dm : DeformationModels.forFaultModel(fm)) {
			if (dm.getRelativeWeight(null) == 0 && !(dm == DeformationModels.GEOLOGIC_LOWER || dm == DeformationModels.GEOLOGIC_UPPER))
				continue;
			Map<Integer, DeformationSection> dmSects = DeformationModelFileParser.load(dm.getDataFileURL(fm));
			
			double[] dmSlips = new double[ids.length];
			for (int i=0; i<ids.length; i++) {
				int id = ids[i];
				DeformationSection dmSect = dmSects.get(id);
				
				double avg = 0d;
				for (double slip : dmSect.getSlips())
					avg += slip;
				avg /= dmSect.getSlips().size();
				dmSlips[i] = avg;
			}
			
			slips.put(dm, dmSlips);
		}
		
		List<DeformationModels> dms = Lists.newArrayList(slips.keySet());
		
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> line = Lists.newArrayList("Section Name", "Section ID");
		for (DeformationModels dm : dms)
			line.add(dm.getName());
		csv.addLine(line);
		
		for (int i=0; i<ids.length; i++) {
			int id = ids[i];
			FaultSectionPrefData sect = sectsMap.get(id);
			line = Lists.newArrayList(sect.getSectionName(), sect.getSectionId()+"");
			for (DeformationModels dm : dms)
				line.add((float)slips.get(dm)[i]+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/dm_slips.csv"));
	}

}
