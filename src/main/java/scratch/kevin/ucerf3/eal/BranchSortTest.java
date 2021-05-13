package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.ZipException;

import com.google.common.base.Preconditions;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMM_Epistemic;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMMs;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_ProbModels;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_Vs30Model;

public class BranchSortTest {
	
	public static void main(String[] args) throws ZipException, IOException {
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(new File(
				"/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		int numTests = 10;
		
		List<List<U3_EAL_LogicTreeBranch>> branchesList = new ArrayList<>();
		
		int numEach = -1;
		
		for (int i=0; i<numTests; i++) {
			List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>();
			for (U3_EAL_ProbModels probModel : U3_EAL_ProbModels.values())
				for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values())
					for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values())
						for (U3_EAL_Vs30Model vs30 : U3_EAL_Vs30Model.values())
							for (LogicTreeBranch tiBranch : cfss.getBranches())
								branches.add(new U3_EAL_LogicTreeBranch(tiBranch, probModel, gmm, gmmEpi, vs30));
			Collections.shuffle(branches);
			Collections.sort(branches);
			if (branchesList.isEmpty())
				numEach = branches.size();
			else
				Preconditions.checkState(branches.size() == numEach);
			branchesList.add(branches);
		}
		
		for (int i=0; i<numEach; i++) {
			U3_EAL_LogicTreeBranch branch0 = branchesList.get(0).get(i);
			for (int j=0; j<numTests; j++) {
				U3_EAL_LogicTreeBranch testBranch = branchesList.get(j).get(i);
				Preconditions.checkState(testBranch.equals(branch0), "Branch mismatch for test %s, index %s\n\tBranch 0: %s\n\tBranch %s: %s",
						j, i, branch0, i, testBranch);
				Preconditions.checkState(testBranch.buildFileName().equals(branch0.buildFileName()));
			}
		}
	}

}
