package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.HazardMapPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;

import com.google.common.base.Preconditions;

public class LogicTreeBranchAverageWriter {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_16-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-ZENGBB-Shaw09Mod-DsrUni-TotNuclRate-SubB1");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
		File mainDir = new File(invDir, "2022_01_16-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-U3RupSet-SubB1-5000ip");
		File resultsFile = new File(mainDir, "results.zip");
		File fullBAFile = new File(mainDir, "results_FM3_1_U3RupSet_branch_averaged.zip");
		
		HazardMapPlot.SPACING_DEFAULT = 0.2;
		
		File outputDir = new File(mainDir, "node_branch_averaged");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		FaultSystemSolution fullBA = fullBAFile == null ? null : FaultSystemSolution.load(fullBAFile);
		PlotLevel plt = PlotLevel.DEFAULT;
		boolean compWithLoaded = false;
		
		LogicTree<?> tree = slt.getLogicTree();
		
//		tree = tree.matchingAll(SupraSeisBValues.B_0p0, DeformationModels.GEOLOGIC,
//				SubSectConstraintModels.TOT_NUCL_RATE, SegmentationModels.SHAW_R0_3);
//		System.out.println("Test sub-tree size: "+tree.size());
//		Preconditions.checkState(tree.size() > 0);
//		outputDir = new File(outputDir, "test_subset");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		plt = PlotLevel.LIGHT;
//		compWithLoaded = true;
		
		Map<LogicTreeLevel<?>, HashSet<LogicTreeNode>> levelNodes = new HashMap<>();
		List<? extends LogicTreeLevel<?>> levels = tree.getLevels();
		for (LogicTreeLevel<?> level : levels)
			levelNodes.put(level, new HashSet<>());
		
		Map<LogicTreeNode, LogicTreeLevel<?>> nodeLevels = new HashMap<>();
		
		for (LogicTreeBranch<?> branch : tree) {
			for (int i=0; i<levels.size(); i++) {
				LogicTreeNode node = branch.getValue(i);
				levelNodes.get(levels.get(i)).add(node);
				if (!nodeLevels.containsKey(node))
					nodeLevels.put(node, levels.get(i));
			}
		}
		
		Map<LogicTreeNode, BranchAverageSolutionCreator> nodeBACreators = new HashMap<>();
		
		for (LogicTreeLevel<?> level : levels) {
			HashSet<LogicTreeNode> nodes = levelNodes.get(level);
			if (nodes.size() > 1) {
				System.out.println("Building "+nodes.size()+" BAs for "+level.getName());
				for (LogicTreeNode node : nodes) {
					BranchAverageSolutionCreator creator = new BranchAverageSolutionCreator(tree.getWeightProvider());
//					creator.skipModule(SupraSeisBValInversionTargetMFDs.class);
					creator.skipModule(InversionTargetMFDs.class);
					nodeBACreators.put(node, creator);
				}
			}
		}
		
		for (LogicTreeBranch<?> branch : tree) {
			FaultSystemSolution sol = slt.forBranch(branch);
			
			for (LogicTreeNode node : branch)
				if (nodeBACreators.containsKey(node))
					nodeBACreators.get(node).addSolution(sol, branch);
		}
		
		RupSetMetadata comparison = fullBA == null ? null : new RupSetMetadata("Full Tree", fullBA);
		
		for (LogicTreeNode node : nodeBACreators.keySet()) {
			System.out.println("Building BA for "+node);
			BranchAverageSolutionCreator creator = nodeBACreators.get(node);
			FaultSystemSolution ba = creator.build();
			
			LogicTreeLevel<?> level = nodeLevels.get(node);
			
			String prefix = level.getShortName().replaceAll("\\W+", "_")+"_"+node.getFilePrefix();
			
			File solFile = new File(outputDir, prefix+".zip");
			ba.write(solFile);
			
			File reportDir = new File(outputDir, prefix+"_report");
			Preconditions.checkState(reportDir.exists() || reportDir.mkdir());
			System.out.println("Writing report to "+reportDir.getAbsolutePath());
			RupSetMetadata primary = new RupSetMetadata(level.getShortName()+": "+node.getShortName(), ba);
			
			if (compWithLoaded)
				// test to compare serialized and regular
				comparison = new RupSetMetadata("Loaded", FaultSystemSolution.load(solFile));
			
			ReportMetadata meta = new ReportMetadata(primary, comparison);
			ReportPageGen pageGen = new ReportPageGen(meta, reportDir, ReportPageGen.getDefaultSolutionPlots(plt));
			pageGen.setReplot(true);
			
			pageGen.generatePage();
		}
	}

}
