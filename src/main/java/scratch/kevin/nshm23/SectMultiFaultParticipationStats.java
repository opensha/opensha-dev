package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;

public class SectMultiFaultParticipationStats {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		List<? extends FaultSection> sects = rupSet.getFaultSectionDataList();
		
		int[] parentIDs = {
				// Art's mendocino to not actually elsinore example
//				FaultSectionUtils.findParentSectionID(sects, "Mendocino"),
////				FaultSectionUtils.findParentSectionID(sects, "San Felipe"),
//				339, // san felipe
//				FaultSectionUtils.findParentSectionID(sects, "Big", "Bend"),
				
				// Art's "all SAF" which actually goes to Rodgers creek
//				FaultSectionUtils.findParentSectionID(sects, "Andreas", "Bernardino", "north"),
//				FaultSectionUtils.findParentSectionID(sects, "Rodgers", "Healdsburg"),
				
				FaultSectionUtils.findParentSectionID(sects, "Cucamonga"),
		};
		
		HashSet<Integer> matches = null;
		
		for (int parentID : parentIDs) {
			if (matches == null)
				matches = new HashSet<>(rupSet.getRupturesForParentSection(parentID));
			else
				matches.retainAll(rupSet.getRupturesForParentSection(parentID));
		}
		
		System.out.println("Found "+matches.size()+" ruptures with all parents");
		
		double rate = 0d;
		for (int rupIndex : matches)
			rate += sol.getRateForRup(rupIndex);
		
		System.out.println("Total rate: "+rate);
		
		if (parentIDs.length == 1) {
			// see how often it ruptures with anything else
			double coruptureRate = 0d;
			for (int rupIndex : matches) {
				for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex)) {
					if (sect.getParentSectionId() != parentIDs[0]) {
						coruptureRate += sol.getRateForRup(rupIndex);
						break;
					}
				}
			}
			DecimalFormat pDF = new DecimalFormat("0.00%");
			System.out.println("Corupture rate: "+coruptureRate+" ("+pDF.format(coruptureRate/rate)+")");
		}
	}

}
