package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.geo.Location;

import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class DoformationModelCombiner {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR.getParentFile(), "DeformationModels");
		Map<Integer, DeformationSection> geologicModel = DeformationModelFileParser.load(
				new File(dir, "geologic_slip_rake_fm_3_1_2012_05_29.csv"));
		Map<Integer, DeformationSection> abmModel = DeformationModelFileParser.load(
				new File(dir, "ABM_slip_rake_fm_3_1_2012_06_08.csv"));
		
		ArrayList<DeformationSection> combined = new ArrayList<DeformationModelFileParser.DeformationSection>();
		
		for (int id : geologicModel.keySet()) {
			DeformationSection geologic = geologicModel.get(id);
			DeformationSection abm = abmModel.get(id);
			
			DeformationSection comb = new DeformationSection(id);
			
			for (int i=0; i<geologic.getLocs1().size(); i++) {
				Location loc1 = geologic.getLocs1().get(i);
				Location loc2 = geologic.getLocs2().get(i);
				double slip1 = geologic.getSlips().get(i);
				if (Double.isNaN(slip1))
					slip1 = 0;
				double slip2 = abm.getSlips().get(i);
				if (Double.isNaN(slip2))
					slip2 = 0;
				double rake = geologic.getRakes().get(i); // keep rake from geologic
				
				double slip = 0.5*(slip1+slip2);
				comb.add(loc1, loc2, slip, rake);
			}
			combined.add(comb);
		}
		
		DeformationModelFileParser.write(combined,
				new File(dir, "geologic_plus_ABM_slip_rake_fm_3_1_2012_06_08.csv"));
	}

}
