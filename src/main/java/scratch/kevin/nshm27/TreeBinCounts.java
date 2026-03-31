package scratch.kevin.nshm27;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeLevel.BinnableLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.BinnedLevel;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

import gov.usgs.earthquake.nshmp.erf.logicTree.TectonicRegionBranchTreeNode;

public class TreeBinCounts {

	public static void main(String[] args) throws IOException {
		LogicTree<?> tree = LogicTree.read(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2026_03_23-nshm26-AMSAM-1000samples-gridded/logic_tree.json"));
//				+ "2026_03_25-nshm26-AMSAM-2000samples-gridded/logic_tree.json"));
//				+ "2026_03_27-nshm26-AMSAM-2000samples-gridded/logic_tree.json"));
				+ "2026_03_27-nshm26-GNMI-2000samples-gridded/logic_tree.json"));
		
		System.out.println("Stats for "+tree.size()+" branches");
		
		for (int t=0; t<tree.getLevels().size(); t++) {
			TectonicRegionBranchTreeNode.Level trtLevel = (TectonicRegionBranchTreeNode.Level)tree.getLevels().get(t);
			System.out.println(trtLevel.getName());
			
			int levels = trtLevel.getNodes().get(0).getValue().size();
			int branches = trtLevel.getNodes().size();
			for (int l=0; l<levels; l++) {
				LogicTreeLevel<?> level = null;
				boolean binned = false;
				Map<LogicTreeNode, Double> nodeWeights = new LinkedHashMap<>();
				double sumWeight = 0d;
				for (int i=0; i<branches; i++) {
					LogicTreeBranch<?> branch = trtLevel.getNodes().get(i).getValue();
					double weight = branch.getOrigBranchWeight();
					sumWeight += weight;
					if (i == 0) {
						level = branch.getLevel(l);
						if (level instanceof BinnableLevel<?,?,?>) {
							level = ((BinnableLevel<?,?,?>)level).toBinnedLevel();
							binned = true;
						}
						for (LogicTreeNode node : level.getNodes())
							nodeWeights.put(node, 0d);
					}
					LogicTreeNode node = branch.getValue(l);
					if (binned)
						node = ((BinnedLevel<?, ?>)level).getBinUnchecked(node);
					if (nodeWeights.containsKey(node))
						nodeWeights.put(node, nodeWeights.get(node) + weight);
					else
						nodeWeights.put(node, weight);
				}
				System.out.println(level.getName());
				for (LogicTreeNode node : nodeWeights.keySet()) {
					double weight = nodeWeights.get(node);
					double nodeWeight = node.getNodeWeight(null);
					if (weight > 0d || nodeWeight > 0d)
						System.out.println("\t"+node.getName()+":\t"+(float)weight
								+"\t(nodeWeight="+(float)nodeWeight+")");
				}
			}
		}
	}

}
