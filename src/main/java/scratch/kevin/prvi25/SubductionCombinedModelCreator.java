package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.util.TrueMeanSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;

import com.google.common.base.Preconditions;

public class SubductionCombinedModelCreator {

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("USAGE: <subduction-dir>");
			System.exit(2);
		}
		File subductionDir = new File(args[0]);
		String suffix = "_branch_averaged_gridded.zip";
		File outputFile = new File(subductionDir, "results_PRVI_SUB_FMs_combined"+suffix);
		
		Map<PRVI25_SubductionFaultModels, FaultSystemSolution> subductionBASols = new HashMap<>();
		
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			File subductionBA = new File(subductionDir, "results_"+fm.getFilePrefix()+suffix);
			if (subductionBA.exists()) {
				FaultSystemSolution sol = FaultSystemSolution.load(subductionBA);
				GridSourceProvider gridProv = sol.getGridSourceProvider();
				if (gridProv == null) {
					// could be pre-abstarct/precomputed refactor
					gridProv = sol.getArchive().loadUnlistedModule(GridSourceList.Precomputed.class, FaultSystemSolution.NESTING_PREFIX);
					Preconditions.checkNotNull(gridProv);
					sol.setGridSourceProvider(gridProv);
				}
				subductionBASols.put(fm, sol);
			}
		}
		Preconditions.checkState(!subductionBASols.isEmpty());
		
		FaultSystemSolution trueMean = combine(subductionBASols);
		trueMean.write(outputFile);
	}

	public static FaultSystemSolution combine(
			Map<PRVI25_SubductionFaultModels, FaultSystemSolution> subductionBASols) {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		List<LogicTreeBranch<LogicTreeNode>> branches = new ArrayList<>();
		
		levels.add(PRVI25_LogicTreeBranch.SUB_FM);
		
		for (PRVI25_SubductionFaultModels subductionFM : subductionBASols.keySet())
			branches.add(new LogicTreeBranch<>(levels, List.of(subductionFM)));
		
		LogicTree<?> tree = LogicTree.fromExisting(levels, branches);
		
		TrueMeanSolutionCreator creator = new TrueMeanSolutionCreator(tree);
		creator.setDoGridProv(true);
		List<BranchRegionalMFDs> regionalMFDs = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		for (LogicTreeBranch<?> branch : tree) {
			FaultSystemSolution subductionSol = subductionBASols.get(branch.requireValue(PRVI25_SubductionFaultModels.class));
			weights.add(tree.getBranchWeight(branch));
			
			creator.addSolution(subductionSol, branch);
			
			if (regionalMFDs != null) {
				BranchRegionalMFDs regMFDs = subductionSol.getModule(BranchRegionalMFDs.class);
				if (regMFDs == null)
					regionalMFDs = null;
				else
					regionalMFDs.add(regMFDs);
			}
		}
		
		FaultSystemSolution trueMean = creator.build();
		RegionsOfInterest roi = subductionBASols.values().iterator().next().getRupSet().getModule(RegionsOfInterest.class);
		if (roi != null)
			trueMean.getRupSet().addModule(roi);
		if (regionalMFDs != null)
			trueMean.addModule(BranchRegionalMFDs.combine(regionalMFDs, weights));
		return trueMean;
	}

}
