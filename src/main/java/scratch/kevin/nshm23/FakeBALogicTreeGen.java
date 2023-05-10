package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.FileBuilder;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.SolutionProcessor;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_LogicTreeBranch;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class FakeBALogicTreeGen {

	public static void main(String[] args) throws IOException {
		
		List<LogicTreeNode> nodes = new ArrayList<>();
		List<FaultSystemSolution> nodeSols = new ArrayList<>();
		
		boolean stripRupMFDs = true;

		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_U3_HybridLogicTreeBranch.U3_FM;
		nodes.add(FaultModels.FM3_1);
		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip")));
		
		nodes.add(FaultModels.FM3_2);
		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_2_branch_averaged.zip")));
		
//		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2023_05_08-u3-both_fms-ba_only");
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2023_05_08-u3-both_fms-ba_only-nshm23_gridded");
//		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2023_01_04-u3-FM3_1-ba_only-nshm23_gridded");
		GridSourceProvider gridProv = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip")).getGridSourceProvider();
		for (FaultSystemSolution sol : nodeSols)
			sol.setGridSourceProvider(gridProv);
		
		SolutionProcessor processor = null;
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		FileBuilder builder = new SolutionLogicTree.FileBuilder(processor, new File(outputDir, "results.zip"));
		
//		File inputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
//				+ "2022_03_24-u3_branches-FM3_1-2000ip");
//		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_U3_HybridLogicTreeBranch.U3_FM;
//		nodes.add(FaultModels.FM3_1);
//		nodeSols.add(FaultSystemSolution.load(new File(inputDir, "results_FM3_1_branch_averaged.zip")));
//		SolutionProcessor processor = null;
//		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only");
		
//		File inputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
////				+ "2022_12_20-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
////				+ "2023_01_25-nshm18_branches-new_scale-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED-TotNuclRate-NoRed-ThreshAvgIterRelGR");
////				+ "2023_04_13-nshm18_branches-new_scale-u3_paleo-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED-TotNuclRate-NoRed-ThreshAvgIterRelGR");
////				+ "2023_04_14-nshm18_branches-wc_94-u3_paleo-NSHM18_WUS_PlusU3_FM_3p1-CoulombRupSet-BRANCH_AVERAGED-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//				+ "2023_04_14-nshm23_u3_hybrid_branches-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//////		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_U3_HybridLogicTreeBranch.U3_FM;
//////		nodes.add(FaultModels.FM3_1);
////		LogicTreeLevel<? extends LogicTreeNode> level = NSHM18_LogicTreeBranch.FM;
////		nodes.add(NSHM18_FaultModels.NSHM18_WUS_PlusU3_FM_3p1);
//////		FaultSystemSolution tmpSol = FaultSystemSolution.load(new File(inputDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip"));
//////		FaultSystemSolution tmpSol = FaultSystemSolution.load(new File(inputDir, "results_NSHM18_WUS_PlusU3_FM_3p1_CoulombRupSet_branch_averaged.zip"));
////		FaultSystemSolution tmpSol = FaultSystemSolution.load(new File(inputDir, "node_branch_averaged/SegModel_Classic.zip"));
//////		FaultSystemSolution tmpSol = FaultSystemSolution.load(new File(inputDir, "results_FM3_1_CoulombRupSet_NoClassic_branch_averaged.zip"));
//////		tmpSol.setGridSourceProvider(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip")).getGridSourceProvider());
////		tmpSol.setGridSourceProvider(FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
////				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
////				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip")).getGridSourceProvider());
////		nodeSols.add(tmpSol);
//////
//////		SolutionProcessor processor = null;
////		SolutionProcessor processor = new NSHM23_InvConfigFactory.NSHM23SolProcessor();
//////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only");
//////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-no_classic");
//////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-nshm23_gridded");
////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-nshm23_gridded-classic_only");
//		
//		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_U3_HybridLogicTreeBranch.U3_FM;
//		nodes.add(FaultModels.FM3_1);
//		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2022_12_20-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_FM3_1_CoulombRupSet_branch_averaged.zip")));
//		nodes.add(FaultModels.FM3_2);
//		nodeSols.add(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2023_02_07-nshm23_u3_hybrid_branches-FM3_2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_FM3_2_CoulombRupSet_branch_averaged.zip")));
//		SolutionProcessor processor = null;
//		
////		GridSourceProvider gridProv = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
////				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
////				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip")).getGridSourceProvider();
////		for (FaultSystemSolution sol : nodeSols)
////			sol.setGridSourceProvider(gridProv);
////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-nshm23_gridded");
//		
//		nodeSols.get(0).setGridSourceProvider(FaultSystemSolution.load(
//				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip")).getGridSourceProvider());
//		nodeSols.get(1).setGridSourceProvider(FaultSystemSolution.load(
//				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_2_branch_averaged.zip")).getGridSourceProvider());
//		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-u3_gridded");
//
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		FileBuilder builder = new SolutionLogicTree.FileBuilder(processor, new File(outputDir, "results.zip"));
		
//		File inputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
//				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//		LogicTreeLevel<? extends LogicTreeNode> level = NSHM23_LogicTreeBranch.FM;
//		nodes.add(NSHM23_FaultModels.NSHM23_v2);
//		nodeSols.add(FaultSystemSolution.load(new File(inputDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip")));
//		
//		SolutionProcessor processor = new NSHM23_InvConfigFactory.NSHM23SolProcessor();
//		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only");
////		nodes.add(NSHM23_FaultModels.NSHM23_v2);
////		nodeSols.add(FaultSystemSolution.load(new File(inputDir, "node_branch_averaged/SegModel_Classic.zip")));
////		
////		SolutionProcessor processor = new NSHM23_InvConfigFactory.NSHM23SolProcessor();
////		File outputDir = new File(inputDir.getParentFile(), inputDir.getName()+"-ba_only-classic_only");
//		
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		FileBuilder builder = new SolutionLogicTree.FileBuilder(processor, new File(outputDir, "results.zip"));
		
		for (int i=0; i<nodes.size(); i++) {
			LogicTreeNode node = nodes.get(i);
			FaultSystemSolution sol = nodeSols.get(i);
			
			System.out.println("Writing "+node+" (has gridded? "+(sol.getGridSourceProvider() != null)+")");
			
			if (stripRupMFDs)
				sol.removeModuleInstances(RupMFDsModule.class);
			
			LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(List.of(level));
			branch.setValue(node);
			
			builder.solution(sol, branch);
		}
		
		builder.build();
	}

}
