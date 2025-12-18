package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.ComparablePairing;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.faultSurface.FaultSection;

public class LargestMminSearch {
	
	public static void main(String[] args) throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		DecimalFormat df = new DecimalFormat("0.0");
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		Map<String, Double> parentMinMags = new HashMap<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			double minMag = rupSet.getMinMagForSection(s);
			String name = sect.getParentSectionName();
			if (!parentMinMags.containsKey(name) || minMag < parentMinMags.get(name))
				parentMinMags.put(name, minMag);
		}
		
		List<String> sorted = ComparablePairing.getSortedData(parentMinMags);
		Collections.reverse(sorted);
		for (int i=0; i<20; i++) {
			String name = sorted.get(i);
			System.out.println(name+":\t"+df.format(parentMinMags.get(name)));
		}
	}

}
