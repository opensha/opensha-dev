package scratch.kevin.nshm23.segModelTests;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.BranchWeightProvider.CurrentWeights;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;

public class MaxDistReweight {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File mainDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
		
		File resultsFile = new File(mainDir, "results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		MaxJumpDistModels.WEIGHT_TARGET_R0 = 3;
		
		System.out.println("Weights for "+MaxJumpDistModels.WEIGHT_TARGET_R0);
		for (MaxJumpDistModels m : MaxJumpDistModels.values())
			System.out.println(m+": "+m.getNodeWeight(null));
		
		BranchAverageSolutionCreator baCreate =
				new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		
		baCreate.skipModule(InversionTargetMFDs.class);
		
		LogicTree<?> tree = slt.getLogicTree();
		
		for (int i=0; i<tree.size(); i++) {
			LogicTreeBranch<?> branch = tree.getBranch(i);
			System.out.println("Processing branch "+i+"/"+tree.size()+": "+branch);
			baCreate.addSolution(slt.forBranch(branch), branch);
		}
		
		FaultSystemSolution sol = baCreate.build();
		
		sol.write(new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged_reweight_r0_"
				+(float)MaxJumpDistModels.WEIGHT_TARGET_R0+".zip"));
		File outputDir = new File(mainDir, "branch_averaged_reweight_r0_"+(float)MaxJumpDistModels.WEIGHT_TARGET_R0);
		
		ReportPageGen pageGen = new ReportPageGen(
				new ReportMetadata(new RupSetMetadata("Reweight R0="+(float)MaxJumpDistModels.WEIGHT_TARGET_R0, sol)),
				outputDir, ReportPageGen.getDefaultSolutionPlots(PlotLevel.FULL));
		pageGen.skipSectBySect();
		
		pageGen.generatePage();
		
		System.exit(0);
	}

}
