package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.ImmutableList;

public class ConvertToTruePointSource {

	public static void main(String[] args) throws IOException {
		File inputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_31-prvi25_subduction_branches");
		File inputFile = new File(inputDir, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip");
		
//		EnumSet<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.SUBDUCTION_SLAB);
//		File outputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_31-prvi25_subduction_branches-ba_only-LARGE-slab_pt_src");
		EnumSet<TectonicRegionType> trts = EnumSet.allOf(TectonicRegionType.class);
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_31-prvi25_subduction_branches-ba_only-LARGE-true_pt_src");
		File outputFile = new File(outputDir, inputFile.getName());
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		GridSourceList gridSources = sol.requireModule(GridSourceList.class);
		
		EnumMap<TectonicRegionType, List<List<GriddedRupture>>> trtRuptureLists = new EnumMap<>(TectonicRegionType.class);
		for (TectonicRegionType trt : gridSources.getTectonicRegionTypes()) {
			List<List<GriddedRupture>> ruptureLists = new ArrayList<>(gridSources.getNumLocations());
			for (int gridIndex=0; gridIndex<gridSources.getNumLocations(); gridIndex++) {
				ImmutableList<GriddedRupture> ruptures = gridSources.getRuptures(trt, gridIndex);
				if (ruptures.isEmpty()) {
					ruptureLists.add(null);
					continue;
				} else {
					List<GriddedRupture> ptRuptures = new ArrayList<>(ruptures.size());
					for (GriddedRupture rup : ruptures) {
						if (trts.contains(rup.properties.tectonicRegionType)) {
							GriddedRuptureProperties props = new GriddedRuptureProperties(gridIndex, rup.properties.location, rup.properties.magnitude,
									rup.properties.rake, rup.properties.dip, Double.NaN, null, rup.properties.upperDepth, rup.properties.upperDepth, 0d,
									Double.NaN, Double.NaN, trt);
							ptRuptures.add(new GriddedRupture(props, rup.rate));
						} else {
							ptRuptures.add(rup);
						}
					}
					ruptureLists.add(ptRuptures);
				}
			}
			trtRuptureLists.put(trt, ruptureLists);
		}
		GridSourceList pointSources = new GridSourceList.Precomputed(gridSources.getGriddedRegion(), trtRuptureLists);
		sol.setGridSourceProvider(pointSources);
		
		sol.write(outputFile);
	}

}
