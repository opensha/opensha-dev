package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.FileBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class FakeBALogicTreeGen {

	public static void main(String[] args) throws IOException {
		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_U3_HybridLogicTreeBranch.U3_FM;
		
		List<LogicTreeNode> nodes = new ArrayList<>();
		List<FaultSystemSolution> nodeSols = new ArrayList<>();
		
		boolean stripRupMFDs = true;
		
		nodes.add(FaultModels.FM3_1);
		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip")));
		
		nodes.add(FaultModels.FM3_2);
		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_2_branch_averaged.zip")));
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2022_09_27-u3-ba_only");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		FileBuilder builder = new SolutionLogicTree.FileBuilder(new File(outputDir, "results.zip"));
		
		for (int i=0; i<nodes.size(); i++) {
			LogicTreeNode node = nodes.get(i);
			FaultSystemSolution sol = nodeSols.get(i);
			
			if (stripRupMFDs)
				sol.removeModuleInstances(RupMFDsModule.class);
			
			LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(List.of(level));
			branch.setValue(node);
			
			builder.solution(sol, branch);
		}
		
		builder.build();
	}

}
