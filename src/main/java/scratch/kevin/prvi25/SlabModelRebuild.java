package scratch.kevin.prvi25;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class SlabModelRebuild {

	public static void main(String[] args) throws IOException {
//		File inputFile = COMBINED_SOL;
//		File inputFile = SUBDUCTION_SOL_LARGE;
		File inputFile = SUBDUCTION_SOL_SMALL;
		
		File outputFile = new File(inputFile.getParentFile(), inputFile.getName().replace(".zip", "")+"_slab_rebuild.zip");
		FaultSystemSolution origSol = FaultSystemSolution.load(inputFile);
		
		GridSourceList origGridList = origSol.requireModule(GridSourceList.class);
		GriddedRegion gridReg = origGridList.getGriddedRegion();
		
		// add in everything but slab
		EnumMap<TectonicRegionType, List<List<GriddedRupture>>> modRupsTRT = new EnumMap<>(TectonicRegionType.class);
		for (TectonicRegionType trt : origGridList.getTectonicRegionTypes()) {
			if (trt == TectonicRegionType.SUBDUCTION_SLAB)
				continue;
			List<List<GriddedRupture>> rupLists = new ArrayList<>(gridReg.getNodeCount());
			for (int l=0; l<gridReg.getNodeCount(); l++)
				rupLists.add(origGridList.getRuptures(trt, l));
			modRupsTRT.put(trt, rupLists);
		}
		
		GridSourceList gridListWithoutSlab = new GridSourceList.Precomputed(gridReg, modRupsTRT);
		
		LogicTree<?> gridTree = LogicTree.buildExhaustive(PRVI25_LogicTree.levelsSubductionGridded, true);
		System.out.println("Built "+gridTree.size()+" branches");
		
		List<CompletableFuture<GridSourceList>> futures = new ArrayList<>();
		for (LogicTreeBranch<?> branch : gridTree) {
			futures.add(CompletableFuture.supplyAsync(new Supplier<GridSourceList>() {

				@Override
				public GridSourceList get() {
					try {
						return PRVI25_GridSourceBuilder.buildSlabGridSourceList(branch);
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
			}));
		}
		
		AveragingAccumulator<GridSourceProvider> gridAvg = null;
		for (int i=0; i<futures.size(); i++) {
			double weight = gridTree.getBranchWeight(i);
			GridSourceList slab = futures.get(i).join();
			System.out.println("Done branch "+i+"/"+gridTree.size());
			if (gridAvg == null)
				gridAvg = slab.averagingAccumulator();
			gridAvg.process(slab, weight);
		}
		
		GridSourceList avgSlab = (GridSourceList)gridAvg.getAverage();
		GridSourceList modGridList = GridSourceList.combine(gridListWithoutSlab, avgSlab);
		
		origSol.setGridSourceProvider(modGridList);
		origSol.write(outputFile);
		
//		ModuleArchive<OpenSHA_Module> slabOnlyArchive = new ModuleArchive<>();
//		slabOnlyArchive.addModule(avgSlab);
//		slabOnlyArchive.write(new File(COMBINED_DIR, "slab_gridded_rebuild.zip"));
	}

}
