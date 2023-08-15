package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_EpistemicLogicTreeNode;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_LogicTreeNode;

import com.google.common.base.Preconditions;

public class GMM_BranchLogicTreeAdd {

	public static void main(String[] args) throws IOException {
//		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//		File inputTreeFile = new File(dir, "logic_tree_full_gridded.json");
		
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only");
		File inputTreeFile = new File(dir, "logic_tree.json");
		
		boolean includeIdriss = false;
		
		boolean allNGAs = true;
		boolean addEpi = false;
		File outputFile = new File(dir, "logic_tree_full_gridded-nga_w2s.json");
		
//		boolean allNGAs = true;
//		boolean addEpi = true;
//		File outputFile = new File(dir, "logic_tree_full_gridded-nga_w2s-gmm_add_epi.json");
		
		LogicTree<?> inputLogicTree = LogicTree.read(inputTreeFile);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> origLevels = new ArrayList<>();
		origLevels.addAll(inputLogicTree.getLevels());
		
		List<LogicTreeLevel<? extends LogicTreeNode>> addLevels = new ArrayList<>();
		if (allNGAs) {
			LogicTreeLevel<NGAW2_LogicTreeNode> level = LogicTreeLevel.forEnum(
					NGAW2_LogicTreeNode.class, "NGA-W2 GMM", "NGA-W2");
			addLevels.add(level);
		}
		if (addEpi) {
			LogicTreeLevel<NGAW2_EpistemicLogicTreeNode> level = LogicTreeLevel.forEnum(
					NGAW2_EpistemicLogicTreeNode.class, "NGA-W2 Additional Epistemic Uncertainty", "AddEpi");
			addLevels.add(level);
		}
		
		LogicTree<?> addTree = LogicTree.buildExhaustive(addLevels, true);
		if (!includeIdriss && allNGAs)
			addTree = addTree.matchingNone(NGAW2_LogicTreeNode.IDRISS_2014);
		
		int numBranches = addTree.size()*inputLogicTree.size();
		System.out.println("Have "+addTree.size()+" new sub-branches, total branches: "+numBranches);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> allLevels = new ArrayList<>();
		allLevels.addAll(origLevels);
		allLevels.addAll(addLevels);
		
		List<LogicTreeBranch<LogicTreeNode>> allBranches = new ArrayList<>(numBranches);
		for (LogicTreeBranch<?> origBranch : inputLogicTree) {
			for (LogicTreeBranch<?> addBranch : addTree) {
				LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(allLevels);
				for (LogicTreeNode node : origBranch)
					branch.setValue(node);
				for (LogicTreeNode node : addBranch)
					branch.setValue(node);
				allBranches.add(branch);
			}
		}
		Preconditions.checkState(allBranches.size() == numBranches);
		
		LogicTree<?> fullTree = LogicTree.fromExisting(allBranches.get(0).getLevels(), allBranches);
		fullTree.write(outputFile);
	}

}
