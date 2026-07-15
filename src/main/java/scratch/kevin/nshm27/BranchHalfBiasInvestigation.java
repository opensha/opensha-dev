package scratch.kevin.nshm27;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;

import gov.usgs.earthquake.nshmp.erf.logicTree.TectonicRegionBranchTreeNode;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceObsSeisDMAdjustment;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceCouplingDepthModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceHingedBValue.CombinedSampledType;

public class BranchHalfBiasInvestigation {

	public static void main(String[] args) throws IOException {
//		File treeFile = new File("/home/kevin/OpenSHA/fss_inversions/2026_06_27-nshm27-AMSAM-3000samples-lhs_pairwise/logic_tree.json");
		File treeFile = new File("/home/kevin/OpenSHA/fss_inversions/2026_07_13-nshm27-AMSAM-5000samples-lhs_pairwise/logic_tree.json");
//		File treeFile = new File("/home/kevin/OpenSHA/fss_inversions/2026_07_11-nshm27-AMSAM-2000samples-lhs_pairwise/logic_tree.json");
		
		LogicTree<?> tree = LogicTree.read(treeFile);
		LogicTree<LogicTreeNode> interfaceTree = null;
		for (int l=0; l<tree.getLevels().size(); l++) {
			LogicTreeLevel<?> level = tree.getLevels().get(l);
			if (level instanceof TectonicRegionBranchTreeNode.Level) {
				interfaceTree = ((TectonicRegionBranchTreeNode.Level)level).getTree();
				break;
			}
		}
		
		int total = tree.size();
		int halfIndex = tree.size()/2;
		int numHinged1 = 0;
		int numHinged2 = 0;
		double sumB1 = 0d;
		int numB1 = 0;
		double sumB2 = 0d;
		int numB2 = 0;
		
		int numExtraploate1 = 0;
		int numExtraploate2 = 0;
		
		int numTaperDouble1 = 0;
		int numTaperNone1 = 0;
		int numTaperDouble2 = 0;
		int numTaperNone2 = 0;
		
		for (int i=0; i<total; i++) {
			LogicTreeBranch<LogicTreeNode> branch = interfaceTree.getBranch(i);
			CombinedSampledType node = branch.requireValue(CombinedSampledType.class);
			boolean first = i < halfIndex;
			if (node.isHinged()) {
				if (first)
					numHinged1++;
				else
					numHinged2++;
			} else {
				double b = node.getB(null, null);
				if (first) {
					sumB1 += b;
					numB1++;
				} else {
					sumB2 += b;
					numB2++;
				}
			}
			if (branch.hasValue(NSHM27_InterfaceObsSeisDMAdjustment.EXTRAPOLATE)) {
				if (first)
					numExtraploate1++;
				else
					numExtraploate2++;
			}
			if (branch.hasValue(NSHM27_InterfaceCouplingDepthModels.DOUBLE_TAPER)) {
				if (first)
					numTaperDouble1++;
				else
					numTaperDouble2++;
			}
			if (branch.hasValue(NSHM27_InterfaceCouplingDepthModels.NONE)) {
				if (first)
					numTaperNone1++;
				else
					numTaperNone2++;
			}
		}
		
		double avgB1 = sumB1 / (double)numB1;
		double avgB2 = sumB2 / (double)numB2;
		System.out.println("First half:");
		System.out.println("\t"+numHinged1+" hinged");
		System.out.println("\t"+avgB1+" average b");
		System.out.println("\t"+numExtraploate1+" extrapolate");
		System.out.println("\t"+numTaperDouble1+" taper-double");
		System.out.println("\t"+numTaperNone1+" taper-none");
		System.out.println("Second half:");
		System.out.println("\t"+numHinged2+" hinged");
		System.out.println("\t"+avgB2+" average b");
		System.out.println("\t"+numExtraploate2+" extrapolate");
		System.out.println("\t"+numTaperDouble2+" taper-double");
		System.out.println("\t"+numTaperNone2+" taper-none");
	}

}
