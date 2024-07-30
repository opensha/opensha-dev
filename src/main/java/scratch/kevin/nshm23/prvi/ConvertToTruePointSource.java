package scratch.kevin.nshm23.prvi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.ImmutableList;

public class ConvertToTruePointSource {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_07_26-prvi25_subduction_branches");
		File inputFile = new File(dir, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip");
		File outputFile = new File(dir, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded_true_pt_src.zip");
		
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
					for (GriddedRupture rup : ruptures)
						ptRuptures.add(new GriddedRupture(gridIndex, rup.location, rup.magnitude, rup.rate,
								rup.rake, rup.dip, Double.NaN, null, rup.upperDepth, rup.upperDepth, 0d,
								Double.NaN, Double.NaN, trt));
					ruptureLists.add(ptRuptures);
				}
			}
			trtRuptureLists.put(trt, ruptureLists);
		}
		GridSourceList pointSources = new GridSourceList(gridSources.getGriddedRegion(), trtRuptureLists);
		sol.setGridSourceProvider(pointSources);
		
		sol.write(outputFile);
	}

}
