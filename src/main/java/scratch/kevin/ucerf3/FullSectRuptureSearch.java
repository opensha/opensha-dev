package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class FullSectRuptureSearch {
	
	public static void main(String[] args) throws IOException {
		int[] parentIDs = {
				285, // cholame
				282, // SB N
		};
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip"));
		
		Map<Integer, List<FaultSection>> parentSectsMap = rupSet.getFaultSectionDataList().stream().collect(
				Collectors.groupingBy(S -> S.getParentSectionId()));
		
		BitSet potentialRups = null;
		for (int i=0; i<parentIDs.length; i++) {
			BitSet myBitSet = new BitSet(rupSet.getNumRuptures());
			for (int r : rupSet.getRupturesForParentSection(parentIDs[i]))
				myBitSet.set(r);
			if (i == 0) {
				potentialRups = myBitSet;
			} else {
				potentialRups.and(myBitSet);
			}
		}
		System.out.println("Found "+potentialRups.cardinality()+" ruptures that break all parents");
		
		// filter down to fully breaking
		for (int parentID : parentIDs) {
			List<FaultSection> sects = parentSectsMap.get(parentID);
			for (FaultSection sect : sects) {
				BitSet myBitSet = new BitSet(rupSet.getNumRuptures());
				for (int r : rupSet.getRupturesForSection(sect.getSectionId()))
					myBitSet.set(r);
				potentialRups.and(myBitSet);
			}
		}
		System.out.println("Found "+potentialRups.cardinality()+" ruptures that fully break all parents");
		
		// find the smallest
		double minMag = Double.POSITIVE_INFINITY;
		int smallest = -1;
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (potentialRups.get(r)) {
				double mag = rupSet.getMagForRup(r);
				if (mag < minMag) {
					smallest = r;
					minMag = mag;
				}
			}
		}
		
		System.out.println("Smallest is a M"+(float)minMag+": "+smallest);
	}

}
