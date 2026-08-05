package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.util.TectonicRegionType;

public class RevisedGridProvUpdate2026 {

	public static void main(String[] args) throws IOException {
		// updates corresponding to the 2026 revision paper
		
		// from Peter via e-mail 3/9/2026:
		// * Set all previously stable grid nodes to 1/3 stable and 2/3 active
		// * Mark all fault-system ruptures as ACTIVE
		
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3");
		
//		File inFile = new File(dir, "results_WUS_FM_v3_branch_averaged.zip");
//		File inFile = new File(dir, "results_WUS_FM_v3_branch_averaged_gridded.zip");
		File inFile = new File(dir, "results_WUS_FM_v3_branch_averaged_gridded_simplified.zip");
		
		File outFile = new File(dir, inFile.getName().substring(0, inFile.getName().indexOf(".zip"))+"_revised2026.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		if (rupSet.hasModule(RupSetTectonicRegimes.class))
			rupSet.addModule(RupSetTectonicRegimes.constant(rupSet, TectonicRegionType.ACTIVE_SHALLOW));
		
		if (sol.hasModule(GridSourceProvider.class)) {
			GridSourceList origGridList = sol.requireModule(GridSourceList.class);
			sol.setGridSourceProvider(updateGridList(origGridList));
		}
		
		sol.write(outFile);
	}
	
	private static GridSourceList updateGridList(GridSourceList origGridList) {
		double ACTIVE_FRACT = 2d/3d;
		double STABLE_FRACT = 1d/3d;
		
		
		EnumMap<TectonicRegionType, List<? extends List<GriddedRupture>>> trtRuptureLists = new EnumMap<>(TectonicRegionType.class);
		List<List<GriddedRupture>> activeList = new ArrayList<>();
		trtRuptureLists.put(TectonicRegionType.ACTIVE_SHALLOW, activeList);
		List<List<GriddedRupture>> stableList = new ArrayList<>();
		trtRuptureLists.put(TectonicRegionType.STABLE_SHALLOW, stableList);
		
		Map<GriddedRuptureProperties, GriddedRuptureProperties> activePropsCache = new HashMap<>();
		
		for (int l=0; l<origGridList.getNumLocations(); l++) {
			List<GriddedRupture> activeRups = origGridList.getRuptures(TectonicRegionType.ACTIVE_SHALLOW, l);
			List<GriddedRupture> stableRups = origGridList.getRuptures(TectonicRegionType.STABLE_SHALLOW, l);
			
			if (stableRups != null) {
				activeRups = new ArrayList<>(activeRups);
				List<GriddedRupture> modStableRups = new ArrayList<>(stableRups.size());
				
				for (GriddedRupture rup : stableRups) {
					GriddedRuptureProperties props = rup.properties;
					GriddedRuptureProperties activeProps = new GriddedRuptureProperties(
							props.magnitude, props.rake, props.dip, props.strike, props.strikeRange,
							props.upperDepth, props.lowerDepth, props.length, props.hypocentralDepth,
							props.hypocentralDAS, TectonicRegionType.ACTIVE_SHALLOW);
					if (activePropsCache.containsKey(activeProps))
						activeProps = activePropsCache.get(activeProps);
					else
						activePropsCache.put(activeProps, activeProps);
					
					activeRups.add(new GriddedRupture(l, rup.location, activeProps,
							rup.rate*ACTIVE_FRACT, rup.associatedSections, rup.associatedSectionFracts));
					modStableRups.add(new GriddedRupture(l, rup.location, props,
							rup.rate*STABLE_FRACT, rup.associatedSections, rup.associatedSectionFracts));
				}
				
				stableRups = modStableRups;
			}
			
			activeList.add(activeRups);
			stableList.add(stableRups);
		}
		
		return new GridSourceList.Precomputed(origGridList.getGriddedRegion(), trtRuptureLists);
	}

}
