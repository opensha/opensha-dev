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
import org.opensha.sha.earthquake.faultSysSolution.hazard.AbstractLogicTreeHazardCombiner;
import org.opensha.sha.earthquake.faultSysSolution.hazard.AbstractLogicTreeHazardCombiner.CombinedRupSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ProxyFaultSectionInstances;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.util.SolModuleStripper;
import org.opensha.sha.earthquake.faultSysSolution.util.TrueMeanSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;

import com.google.common.base.Preconditions;

public class CrustalSubductionTrueMeanCreator {

	public static void main(String[] args) throws IOException {
		if (args.length < 3 || args.length > 4) {
			System.err.println("USAGE: <crustal-dir> <subduction-dir> <output-file> [<gridded?>]");
			System.exit(2);
		}
		File crustalDir = new File(args[0]);
		File subductionDir = new File(args[1]);
		File outputFile = new File(args[2]);
		boolean gridded = args.length < 4 ? true : Boolean.parseBoolean(args[3]);
		
		Map<PRVI25_CrustalFaultModels, FaultSystemSolution> crustalBASols = new HashMap<>();
		Map<PRVI25_SubductionFaultModels, FaultSystemSolution> subductionBASols = new HashMap<>();
		
		String simplifiedSuffix = gridded ? "_branch_averaged_gridded_simplified.zip" : "_branch_averaged_simplified.zip";
		String suffix = gridded ? "_branch_averaged_gridded.zip" : "_branch_averaged.zip";
		
		for (PRVI25_CrustalFaultModels fm : PRVI25_CrustalFaultModels.values()) {
			File crustalBA = new File(crustalDir, "results_"+fm.getFilePrefix()+simplifiedSuffix);
			boolean simplify = false;
			if (!crustalBA.exists()) {
				// try non-simplified
				crustalBA = new File(crustalDir, "results_"+fm.getFilePrefix()+suffix);
				simplify = true;
			}
			if (crustalBA.exists()) {
				FaultSystemSolution sol = FaultSystemSolution.load(crustalBA);
				if (gridded) {
					GridSourceProvider gridProv = sol.getGridSourceProvider();
					if (gridProv == null) {
						// could be pre-abstarct/precomputed refactor
						gridProv = sol.getArchive().loadUnlistedModule(GridSourceList.Precomputed.class, FaultSystemSolution.NESTING_PREFIX);
						Preconditions.checkNotNull(gridProv);
						sol.setGridSourceProvider(gridProv);
					}
					if (gridProv instanceof MFDGridSourceProvider) {
						gridProv = GridSourceList.convert(
								(MFDGridSourceProvider)sol.getGridSourceProvider(),
								sol.getRupSet().requireModule(FaultGridAssociations.class),
								new NSHM23_SingleRegionGridSourceProvider.NSHM23_WUS_FiniteRuptureConverter());
						sol.setGridSourceProvider(gridProv);
					}
				}
				if (simplify)
					// simplify it (removes proxies)
					sol = SolModuleStripper.stripModules(sol, 5d, true, false);
				crustalBASols.put(fm, sol);
			}
		}
		Preconditions.checkState(!crustalBASols.isEmpty());
		
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			File subductionBA = new File(subductionDir, "results_"+fm.getFilePrefix()+simplifiedSuffix);
			if (!subductionBA.exists())
				// try non-simplified
				subductionBA = new File(subductionDir, "results_"+fm.getFilePrefix()+suffix);
			if (subductionBA.exists()) {
				FaultSystemSolution sol = FaultSystemSolution.load(subductionBA);
				if (gridded) {
					GridSourceProvider gridProv = sol.getGridSourceProvider();
					if (gridProv == null) {
						// could be pre-abstarct/precomputed refactor
						gridProv = sol.getArchive().loadUnlistedModule(GridSourceList.Precomputed.class, FaultSystemSolution.NESTING_PREFIX);
						Preconditions.checkNotNull(gridProv);
						sol.setGridSourceProvider(gridProv);
					}
				}
				subductionBASols.put(fm, sol);
			}
		}
		Preconditions.checkState(subductionBASols.size() == 2);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		List<LogicTreeBranch<LogicTreeNode>> branches = new ArrayList<>();
		
		levels.add(PRVI25_LogicTreeBranch.CRUSTAL_FM);
		levels.add(PRVI25_LogicTreeBranch.SUB_FM);
		
		for (PRVI25_CrustalFaultModels crustalFM : crustalBASols.keySet())
			for (PRVI25_SubductionFaultModels subductionFM : subductionBASols.keySet())
				branches.add(new LogicTreeBranch<>(levels, List.of(crustalFM, subductionFM)));
		
		LogicTree<?> tree = LogicTree.fromExisting(levels, branches);
		
		TrueMeanSolutionCreator creator = new TrueMeanSolutionCreator(tree);
		creator.setDoGridProv(gridded);
		for (LogicTreeBranch<?> branch : tree) {
			FaultSystemSolution crustalSol = crustalBASols.get(branch.requireValue(PRVI25_CrustalFaultModels.class));
			Preconditions.checkState(crustalSol.getRupSet().hasModule(RupSetTectonicRegimes.class), "Crustal solution doesn't have TRTs");
			FaultSystemSolution subductionSol = subductionBASols.get(branch.requireValue(PRVI25_SubductionFaultModels.class));
			Preconditions.checkState(subductionSol.getRupSet().hasModule(RupSetTectonicRegimes.class), "Subduction solution doesn't have TRTs");
			
			FaultSystemSolution combined = AbstractLogicTreeHazardCombiner.combineSols(crustalSol, subductionSol, true);
			Preconditions.checkState(combined.getRupSet().hasModule(RupSetTectonicRegimes.class), "Combined solution doesn't have TRTs");
			if (gridded) {
				GridSourceList crustalGridded = crustalSol.requireModule(GridSourceList.class);
				GridSourceList subductionGridded = subductionSol.requireModule(GridSourceList.class);
				CombinedRupSetMappings mappings = combined.getRupSet().requireModule(CombinedRupSetMappings.class);
				crustalGridded = GridSourceList.remapAssociations(crustalGridded, mappings.getInnerSectMappings());
				subductionGridded = GridSourceList.remapAssociations(subductionGridded, mappings.getOuterSectMappings());
				combined.setGridSourceProvider(GridSourceList.combine(subductionGridded, crustalGridded));
			}
			
			creator.addSolution(combined, branch);
		}
		
		FaultSystemSolution trueMean = creator.build();
		Preconditions.checkState(trueMean.getRupSet().hasModule(RupSetTectonicRegimes.class), "True mean solution doesn't have TRTs");
		trueMean.write(outputFile);
	}

}
