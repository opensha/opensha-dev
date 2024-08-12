package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.geo.CubedGriddedRegion;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModelRegion;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_CombinedRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_MaxMagOffFault;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;

public class GriddedAddCEUS {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File baSolFile = new File(invDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File outputFile = new File(invDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded_with_ceus.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(baSolFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		MFDGridSourceProvider wusProv = sol.requireModule(MFDGridSourceProvider.class);
		FaultCubeAssociations wusCubes = rupSet.requireModule(FaultCubeAssociations.class);
		
		// average across Mmax
		LogicTreeBranch<LogicTreeNode> defaultBranch = NSHM23_LogicTreeBranch.DEFAULT_OFF_FAULT;
		
		FaultCubeAssociations ceusCubes = null;
		AveragingAccumulator<GridSourceProvider> ceusProvAccumulator = null;
		for (NSHM23_MaxMagOffFault mMax : NSHM23_MaxMagOffFault.values()) {
			if (mMax.getNodeWeight(null) == 0d)
				continue;
			System.out.println("Building CEUS for "+mMax.getShortName());
			LogicTreeBranch<LogicTreeNode> branch = defaultBranch.copy();
			branch.setValue(mMax);
			if (ceusCubes == null) {
				System.out.println("Initializing CEUS cube associations");
				ceusCubes = NSHM23_InvConfigFactory.buildFaultCubeAssociations(
						rupSet, branch, SeismicityRegions.CONUS_EAST.load());
			}
			GridSourceProvider prov = NSHM23_InvConfigFactory.buildGridSourceProv(
					sol, branch, List.of(SeismicityRegions.CONUS_EAST), ceusCubes);
			if (ceusProvAccumulator == null)
				ceusProvAccumulator = prov.averagingAccumulator();
			ceusProvAccumulator.process(prov, branch.getBranchWeight());
		}
		
		System.out.println("Building averaged CEUS grid prov");
		MFDGridSourceProvider ceusProv = (MFDGridSourceProvider)ceusProvAccumulator.getAverage();
		
		System.out.println("Building stitched cubes");
		GriddedRegion stitchedReg = NSHM23_InvConfigFactory.getGriddedSeisRegion(
				List.of(SeismicityRegions.CONUS_WEST, SeismicityRegions.CONUS_EAST));
		CubedGriddedRegion cgr = new CubedGriddedRegion(stitchedReg);
		FaultCubeAssociations stitchedCubes = FaultCubeAssociations.stitch(cgr, List.of(wusCubes, ceusCubes));
		
		NSHM23_CombinedRegionGridSourceProvider stitchedProv = new NSHM23_CombinedRegionGridSourceProvider(
				sol, stitchedCubes, List.of(wusProv, ceusProv));
		
		rupSet.addModule(stitchedCubes);
		rupSet.addModule(new ModelRegion(NSHM23_RegionLoader.loadFullConterminousUS()));
		sol.setGridSourceProvider(stitchedProv);
		
		sol.write(outputFile);
	}

}
