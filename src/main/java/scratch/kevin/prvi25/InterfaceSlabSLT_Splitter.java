package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class InterfaceSlabSLT_Splitter {

	public static void main(String[] args) throws IOException {
		File dir;
		boolean griddedOnly;
		
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			System.out.println("HARDCODED");
			File baseDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions");
			dir = new File(baseDir, "2024_08_16-prvi25_subduction_branches");
			griddedOnly = true;
		} else if (args.length == 2) {
			dir = new File(args[0]);
			griddedOnly = Boolean.parseBoolean(args[1]);
		} else {
			throw new IllegalArgumentException("USAGE: <dir> <gridded-only?>");
		}
		System.out.println("Input directory: "+dir.getAbsolutePath());
		System.out.println("Gridded only? "+griddedOnly);
		Preconditions.checkState(dir.exists());
		
		if (griddedOnly) {
			File intputFile = new File(dir, "results_gridded_branches.zip");
			File slabOutputFile = new File(dir, "results_gridded_branches_slab_only.zip");
			File interfaceOutputFile = new File(dir, "results_gridded_branches_interface_only.zip");
			
			SolutionLogicTree slt = SolutionLogicTree.load(intputFile);
			LogicTree<?> tree = slt.getLogicTree();
			
			SolutionLogicTree.FileBuilder slabWriter = new SolutionLogicTree.FileBuilder(slabOutputFile);
			SolutionLogicTree.FileBuilder interfaceWriter = new SolutionLogicTree.FileBuilder(interfaceOutputFile);
			// we'll manually write them
			slabWriter.setSerializeGridded(false);
			interfaceWriter.setSerializeGridded(false);
			for (int b=0; b<tree.size(); b++) {
				LogicTreeBranch<?> branch = tree.getBranch(b);
				System.out.println("Processing "+b+"/"+tree.size()+": "+branch);
				FaultSystemSolution sol = slt.forBranch(branch);
				GridSourceList gridSources = sol.requireModule(GridSourceList.class);
				
				List<List<GriddedRupture>> slabRupturesList = new ArrayList<>(gridSources.getNumLocations());
				List<List<GriddedRupture>> interfaceRupturesList = new ArrayList<>(gridSources.getNumLocations());
				
				for (int l=0; l<gridSources.getNumLocations(); l++) {
					slabRupturesList.add(gridSources.getRuptures(TectonicRegionType.SUBDUCTION_SLAB, l));
					interfaceRupturesList.add(gridSources.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l));
				}

				GridSourceList slabSources = new GridSourceList.Precomputed(gridSources.getGriddedRegion(),
						TectonicRegionType.SUBDUCTION_SLAB, slabRupturesList);
				GridSourceList interfaceSources = new GridSourceList.Precomputed(gridSources.getGriddedRegion(),
						TectonicRegionType.SUBDUCTION_INTERFACE, interfaceRupturesList);
				
				slabWriter.solution(sol, branch);
				slabWriter.writeGridProvToArchive(slabSources, branch);
				interfaceWriter.solution(sol, branch);
				interfaceWriter.writeGridProvToArchive(interfaceSources, branch);
			}
			slabWriter.close();
			interfaceWriter.close();
			
			tree.write(new File(dir, "logic_tree_gridded_only.json"));
		} else {
			File intputDir = new File(dir, "results");
			File outputFile = new File(dir, "results_full_gridded_interface_only.zip");
			LogicTree<?> tree = LogicTree.read(new File(dir, "logic_tree_full_gridded.json"));
			
			SolutionLogicTree slt = new SolutionLogicTree.ResultsDirReader(intputDir, tree);
			
			SolutionLogicTree.FileBuilder interfaceWriter = new SolutionLogicTree.FileBuilder(outputFile);
			// we'll manually write them
			interfaceWriter.setSerializeGridded(false);
			for (int b=0; b<tree.size(); b++) {
				LogicTreeBranch<?> branch = tree.getBranch(b);
				System.out.println("Processing "+b+"/"+tree.size()+": "+branch);
				FaultSystemSolution sol = slt.forBranch(branch);
				GridSourceList gridSources = (GridSourceList)sol.getGridSourceProvider();
				
				List<List<GriddedRupture>> interfaceRupturesList = new ArrayList<>(gridSources.getNumLocations());
				
				for (int l=0; l<gridSources.getNumLocations(); l++) {
					interfaceRupturesList.add(gridSources.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l));
				}
				
				GridSourceList interfaceSources = new GridSourceList.Precomputed(gridSources.getGriddedRegion(),
						TectonicRegionType.SUBDUCTION_INTERFACE, interfaceRupturesList);
				
				interfaceWriter.solution(sol, branch);
				interfaceWriter.writeGridProvToArchive(interfaceSources, branch);
			}
			interfaceWriter.close();
		}
		
		System.out.println("DONE");
	}

}
