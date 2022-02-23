package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.U3InversionConfigFactory;

public class OrigU3BACreator {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2021_11_30-u3_branches-orig_calcs-5h");
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		
		slt.setProcessor(new U3InversionConfigFactory.UCERF3_SolutionProcessor());
		
		LogicTree<?> tree = slt.getLogicTree();
		
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(tree.getWeightProvider());
		
		int count = 0;
		for (LogicTreeBranch<?> branch : tree) {
			if (branch.hasValue(FaultModels.FM3_1)) {
				System.out.println("Branch "+(count++)+": "+branch);
				baCreator.addSolution(slt.forBranch(branch), branch);
			}
		}
		
		FaultSystemSolution baSol = baCreator.build();
		
		baSol.write(new File(dir, "results_FM3_1_branch_averaged.zip"));
	}

}
