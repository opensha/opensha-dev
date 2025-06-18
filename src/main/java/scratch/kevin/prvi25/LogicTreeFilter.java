package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;

import com.google.common.base.Preconditions;

public class LogicTreeFilter {

	public static void main(String[] args) throws IOException {
//		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2025_05_09-prvi25_crustal_subduction_combined_branches/");
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2025_05_09-prvi25_crustal_branches-dmSample10x/");
		
		File inputFile = new File(dir, "logic_tree_full_gridded.json");
		LogicTree<?> tree = LogicTree.read(inputFile);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		for (LogicTreeLevel<?> level : tree.getLevels())
			levels.add(level);
		List<LogicTreeBranch<LogicTreeNode>> branches = new ArrayList<>();
		
		boolean evenWeight = true;
		File outputFile = new File(dir, "logic_tree_full_gridded_3bvalues_even.json");
//		boolean evenWeight = false;
//		File outputFile = new File(dir, "logic_tree_full_gridded_3bvalues_25_50_25.json");
		for (LogicTreeBranch<?> branch : tree) {
			boolean skip = false;
			List<LogicTreeNode> bValNodes = new ArrayList<>();
			for (LogicTreeNode node : branch) {
				if (node.getShortName().startsWith("b=")) {
					bValNodes.add(node);
					if (node.getShortName().contains("b=0.25") || node.getShortName().contains("b=0.75"))
						skip = true;
				}
			}
			if (!skip) {
				Preconditions.checkState(!bValNodes.isEmpty());
				LogicTreeBranch<LogicTreeNode> modBranch = cleanBranch(branch, levels);
				if (!evenWeight) {
					double branchWeight = branch.getOrigBranchWeight();
					for (LogicTreeNode node : bValNodes) {
						// clear out old weight
						branchWeight /= node.getNodeWeight(branch);
						if (node.getShortName().contains("b=0.5"))
							branchWeight *= 0.5;
						else
							branchWeight *= 0.25;
					}
					modBranch.setOrigBranchWeight(branchWeight);
				}
				branches.add(modBranch);
			}
		}
		
		
//		File outputFile = new File(dir, "logic_tree_full_gridded_mue_with_car.json");
//		int mueIndex = -1;
//		int carIndex = -1;
//		for (int l=0; l<levels.size(); l++) {
//			if (levels.get(l).getName().contains("Muertos"))
//				mueIndex = l;
//			else if (levels.get(l).getName().contains("Caribbean"))
//				carIndex = l;
//		}
//		double origLowWeight = 0d;
//		double origPrefWeight = 0d;
//		double origHighWeight = 0d;
//		double modLowWeight = 0d;
//		double modPrefWeight = 0d;
//		double modHighWeight = 0d;
//		for (LogicTreeBranch<?> branch : tree) {
//			String carVal = branch.getValue(carIndex).getShortName();
//			String mueVal = branch.getValue(mueIndex).getShortName();
//			if (carVal.equals("Low") || mueVal.equals("Low"))
//				origLowWeight += tree.getBranchWeight(branch);
//			if (carVal.equals("Preferred") || mueVal.equals("Preferred"))
//				origPrefWeight += tree.getBranchWeight(branch);
//			if (carVal.equals("High") || mueVal.equals("High"))
//				origHighWeight += tree.getBranchWeight(branch);
//			if (carVal.equals(mueVal)) {
//				if (carVal.equals("Low"))
//					modLowWeight += tree.getBranchWeight(branch);
//				if (carVal.equals("Preferred"))
//					modPrefWeight += tree.getBranchWeight(branch);
//				if (carVal.equals("High"))
//					modHighWeight += tree.getBranchWeight(branch);
//				branches.add(cleanBranch(branch, levels));
//			}
//		}
//		// now re-weight
//		double origFractLow = origLowWeight/(origLowWeight+origPrefWeight+origHighWeight);
//		double origFractPref = origPrefWeight/(origLowWeight+origPrefWeight+origHighWeight);
//		double origFractHigh = origHighWeight/(origLowWeight+origPrefWeight+origHighWeight);
//		double modFractLow = modLowWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		double modFractPref = modPrefWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		double modFractHigh = modHighWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		double lowScalar = origFractLow/modFractLow;
//		double prefScalar = origFractPref/modFractPref;
//		double highScalar = origFractHigh/modFractHigh;
//		System.out.println("Low had "+(float)origFractLow+" weight, reduced to "+(float)modFractLow+", will scale by "+(float)lowScalar);
//		System.out.println("Pref had "+(float)origFractPref+" weight, reduced to "+(float)modFractPref+", will scale by "+(float)prefScalar);
//		System.out.println("High had "+(float)origFractHigh+" weight, reduced to "+(float)modFractHigh+", will scale by "+(float)highScalar);
//		modLowWeight = 0d;
//		modPrefWeight = 0d;
//		modHighWeight = 0d;
//		for (int b=0; b<branches.size(); b++) {
//			LogicTreeBranch<LogicTreeNode> branch = branches.get(b);
//			String carVal = branch.getValue(carIndex).getShortName();
//			if (carVal.equals("Low")) {
//				double weight = branch.getOrigBranchWeight()*lowScalar;
//				modLowWeight += weight;
//				branch.setOrigBranchWeight(weight);
//			} else if (carVal.equals("High")) {
//				double weight = branch.getOrigBranchWeight()*highScalar;
//				modHighWeight += weight;
//				branch.setOrigBranchWeight(weight);
//			} else {
//				double weight = branch.getOrigBranchWeight()*prefScalar;
//				modPrefWeight += weight;
//				branch.setOrigBranchWeight(weight);
//			}
//		}
//		modFractLow = modLowWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		modFractPref = modPrefWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		modFractHigh = modHighWeight/(modLowWeight+modPrefWeight+modHighWeight);
//		System.out.println("Final low weight: "+(float)modFractLow);
//		System.out.println("Final pref weight: "+(float)modFractPref);
//		System.out.println("Final high weight: "+(float)modFractHigh);
		
		System.out.println("Retained "+branches.size()+"/"+tree.size()+" branches");
		LogicTree<?> modTree = LogicTree.fromExisting(levels, branches);
		modTree.setWeightProvider(new BranchWeightProvider.OriginalWeights());
		
		modTree.write(outputFile);
	}
	
	private static LogicTreeBranch<LogicTreeNode> cleanBranch(LogicTreeBranch<?> branch, List<LogicTreeLevel<? extends LogicTreeNode>> levels) {
		LogicTreeBranch<LogicTreeNode> ret = new LogicTreeBranch<>(levels);
		for (int i=0; i<branch.size(); i++)
			ret.setValue(i, branch.getValue(i));
		ret.setOrigBranchWeight(branch.getOrigBranchWeight());
		return ret;
	}

}
