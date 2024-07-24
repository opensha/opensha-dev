package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.AffectsNone;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.EnumBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.imr.AttenRelSupplier;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_AttenRelSupplier;
import org.opensha.sha.imr.logicTree.ScalarIMR_LogicTreeNode;
import org.opensha.sha.imr.logicTree.ScalarIMR_ParamsLogicTreeNode;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.gmm.Gmm;

public class GMMLogicTreeWriter {

	public static void main(String[] args) throws IOException {
		File sourceDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_01-prvi25_crustal_branches-dmSample5x");
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_01-prvi25_crustal_branches-dmSample5x-gmTreeCalcs");
		
		LogicTree<?> tree = LogicTree.read(new File(sourceDir, "logic_tree.json"));
		System.out.println("Read "+tree.size()+" branches");
		
		EnumBackedLevel<NGAW2_Node> gmmLevel = new EnumBackedLevel<>("GMM", "GMM", NGAW2_Node.class);
		
		System.out.println("GMM level check worked? "+isGMMLevel(gmmLevel));
		
		List<LogicTreeLevel<? extends LogicTreeNode>> combLevels = new ArrayList<>();
		combLevels.addAll(tree.getLevels());
		combLevels.add(gmmLevel);
		
		List<LogicTreeBranch<LogicTreeNode>> combBranches = new ArrayList<>();
		
		for (LogicTreeBranch<?> branch : tree) {
			for (NGAW2_Node gmmNode : NGAW2_Node.values()) {
				LogicTreeBranch<LogicTreeNode> combBranch = new LogicTreeBranch<>(combLevels);
				for (LogicTreeNode node : branch)
					combBranch.setValue(node);
				combBranch.setValue(gmmNode);
				combBranches.add(combBranch);
			}
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		LogicTree<?> combTree = LogicTree.fromExisting(combLevels, combBranches);
		combTree.write(new File(outputDir, "logic_tree.json"));
		System.out.println("Write "+combTree.size()+" branches");
		
		combTree = LogicTree.read(new File(outputDir, "logic_tree.json"));
		
		LogicTreeBranch<?> prevBranch = null;
		for (int i=0; i<combTree.size(); i++) {
			LogicTreeBranch<?> branch = combTree.getBranch(i);
			if (prevBranch != null) {
				boolean onlyGmmDifferent = true;
				for (int l=0; l<branch.size(); l++) {
					if (!branch.getValue(l).equals(prevBranch.getValue(l))) {
						if (!isGMMLevel(branch.getLevel(l))) {
							onlyGmmDifferent = false;
							break;
						}
					}
				}
				if (onlyGmmDifferent)
					System.out.println("Only GMM branches differ for "+i+" (relative to "+(i-1)+")");
			}
			prevBranch = branch;
		}
	}
	
	private static boolean isGMMLevel(LogicTreeLevel<?> level) {
		return ScalarIMR_LogicTreeNode.class.isAssignableFrom(level.getType()) ||
				ScalarIMR_ParamsLogicTreeNode.class.isAssignableFrom(level.getType());
	}
	
	@AffectsNone
	public static enum NGAW2_Node implements ScalarIMR_LogicTreeNode {
		ASK(Gmm.ASK_14_BASE),
		BSSA(Gmm.BSSA_14_BASE),
		CB(Gmm.CB_14_BASE),
		CY(Gmm.CY_14_BASE);
		
		private Gmm gmm;

		private NGAW2_Node(Gmm gmm) {
			this.gmm = gmm;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return 1d;
		}

		@Override
		public String getFilePrefix() {
			return gmm.name();
		}

		@Override
		public AttenRelSupplier getSupplier() {
			return new NSHMP_AttenRelSupplier(gmm, false);
		}
		
	}

}
