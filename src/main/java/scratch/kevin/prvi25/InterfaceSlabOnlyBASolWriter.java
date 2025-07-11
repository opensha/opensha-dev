package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.SolModuleStripper;
import org.opensha.sha.earthquake.faultSysSolution.util.TrueMeanSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.util.TectonicRegionType;

public class InterfaceSlabOnlyBASolWriter {

	public static void main(String[] args) throws IOException {
		File dir;
		if (args.length == 0)
			dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_09_04-prvi25_subduction_branches");
		else
			dir = new File(args[0]);
		File largeInputFile = new File(dir, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip");
		File smallInputFile = new File(dir, "results_PRVI_SUB_FM_SMALL_branch_averaged_gridded.zip");
		File slabOutputFile = new File(dir, "results_PRVI_SLAB_ONLY_branch_averaged_gridded.zip");
		File interfaceOutputFile = new File(dir, "results_PRVI_INTERFACE_ONLY_branch_averaged_gridded.zip");

		FaultSystemSolution largeSol = FaultSystemSolution.load(largeInputFile);
		FaultSystemSolution smallSol = FaultSystemSolution.load(smallInputFile);
		
		GridSourceList orig = largeSol.requireModule(GridSourceList.class);
		double[] zeroRates = new double[largeSol.getRupSet().getNumRuptures()];
		FaultSystemSolution slabSol = new FaultSystemSolution(largeSol.getRupSet(), zeroRates);
		slabSol = SolModuleStripper.stripModules(slabSol, 0);
		
		List<List<GriddedRupture>> rupLists = new ArrayList<>(orig.getNumLocations());
		for (int i=0; i<orig.getNumLocations(); i++)
			rupLists.add(orig.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, i));
		GridSourceList slabOnly = new GridSourceList.Precomputed(
				orig.getGriddedRegion(), TectonicRegionType.SUBDUCTION_SLAB, rupLists);
		slabSol.setGridSourceProvider(slabOnly);
		slabSol.write(slabOutputFile);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		List<LogicTreeBranch<LogicTreeNode>> branches = new ArrayList<>();
		List<FaultSystemSolution> sols = new ArrayList<>();
		
		levels.add(PRVI25_LogicTree.SUB_FM);
		
		sols.add(largeSol);
		branches.add(new LogicTreeBranch<>(levels, List.of(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE)));
		sols.add(smallSol);
		branches.add(new LogicTreeBranch<>(levels, List.of(PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL)));
		
		LogicTree<?> tree = LogicTree.fromExisting(levels, branches);
		
		TrueMeanSolutionCreator creator = new TrueMeanSolutionCreator(tree);
		creator.setDoGridProv(true);
		for (int s=0; s<sols.size(); s++) {
			FaultSystemSolution subductionSol = sols.get(s);
			
			orig = subductionSol.requireModule(GridSourceList.class);
			rupLists = new ArrayList<>(orig.getNumLocations());
			for (int i=0; i<orig.getNumLocations(); i++)
				rupLists.add(orig.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, i));
			subductionSol.setGridSourceProvider(new GridSourceList.Precomputed(
					orig.getGriddedRegion(), TectonicRegionType.SUBDUCTION_INTERFACE, rupLists));
			
			creator.addSolution(subductionSol, branches.get(s));
		}
		
		FaultSystemSolution trueMean = creator.build();
		trueMean.write(interfaceOutputFile);
	}

}
