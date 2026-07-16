package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class BranchMFDExample {

	public static void main(String[] args) throws IOException {
		File solfile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solfile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		BranchAveragingOrder branches = sol.requireModule(BranchAveragingOrder.class);
		BranchSectParticMFDs sectParticMFDs = sol.requireModule(BranchSectParticMFDs.class);
		BranchParentSectParticMFDs parentParticMFDs = sol.requireModule(BranchParentSectParticMFDs.class);
		int numBranches = branches.getNumBranches();
		int numSects = rupSet.getNumSections();
		System.out.println(numBranches+" branches");
		double sumWeight = 0d;
		for (int b=0; b<numBranches; b++)
			sumWeight += branches.getBranchWeight(b);
		for (int b=0; b<numBranches; b++) {
			double weight = branches.getBranchWeight(b) / sumWeight;
			System.out.println("Branch "+b+" (wt="+(float)weight+"):\t"+branches.getBranchFileName(b));
			for (int s=0; s<numSects; s++) {
				IncrementalMagFreqDist sectParticMFD = sectParticMFDs.getSectionMFD(b, s);
				Preconditions.checkNotNull(sectParticMFD);
				int parentID = rupSet.getFaultSectionData(s).getParentSectionId();
				IncrementalMagFreqDist parentParticMFD = parentParticMFDs.getSectionMFD(b, parentID);
				Preconditions.checkNotNull(parentParticMFD);
			}
		}
	}

}