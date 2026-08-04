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

import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider.NSHM23_WUS_FiniteRuptureConverter;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupturePropertiesBuilder;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupturePropertiesCache;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

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

		double stableToActiveFract = 2d/3d;
		double stableStaysStableFract = 1d/3d;
		boolean updateRakes = true;
		boolean applyStablePropsToOverlapActive = true;
		
		File outFile = new File(dir, inFile.getName().substring(0, inFile.getName().indexOf(".zip"))+"_revised2026.zip");
		
//		updateRakes = false;
//		File outFile = new File(dir, inFile.getName().substring(0, inFile.getName().indexOf(".zip"))+"_revised2026_origRakes.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		if (rupSet.hasModule(RupSetTectonicRegimes.class))
			rupSet.addModule(RupSetTectonicRegimes.constant(rupSet, TectonicRegionType.ACTIVE_SHALLOW));
		
		if (updateRakes) {
			// update rakes to match Peter's rounding, which unfortunately matters because it pushes ruptures across GMM
			// fault style bins, e.g., -135 is normal but -135.4 is SS. (although Peter's bins are also wrong, it should be
			// 150 for most NGA-W2's).
			// also apply minor rounding to mags via the build(true) call below 
			double[] rakes = new double[rupSet.getNumRuptures()];
			for (int i=0; i<rakes.length; i++)
				rakes[i] = DataUtils.roundFixed(rupSet.getAveRakeForRup(i), 0);
			rupSet = FaultSystemRupSet.buildFromExisting(rupSet, true).rupRakes(rakes).build(true);
			List<OpenSHA_Module> prevSolModules = sol.getModules(true);
			sol = new FaultSystemSolution(rupSet, sol.getRateForAllRups());
			for (OpenSHA_Module module : prevSolModules)
				sol.addModule(module);
		}
		
		// remove rup MFDs (not used by Peter)
		sol.removeModuleInstances(RupMFDsModule.class);
		
		// now update the grid list
		if (sol.getGridSourceProvider() != null) {
			GridSourceList origGridList = sol.requireModule(GridSourceList.class);
			System.out.println("Orig grid list Mmin="+(float)minMag(origGridList));
			sol.setGridSourceProvider(updateGridList(origGridList, stableToActiveFract, stableStaysStableFract, applyStablePropsToOverlapActive));
		}
		
		sol.write(outFile);
	}
	
	private static double minMag(GridSourceList gridList) {
		double minMag = Double.POSITIVE_INFINITY;
		for (int l=0; l<gridList.getNumLocations(); l++)
			for (GriddedRupture rup : gridList.getRuptures(null, l))
				minMag = Math.min(minMag, rup.properties.magnitude);
		return minMag;
	}
	
	public static GridSourceList updateGridList(GridSourceList origGridList, double stabletoActiveFract,
			double stableStaysStableFract, boolean applyStablePropsToOverlapActive) {
		EnumMap<TectonicRegionType, List<? extends List<GriddedRupture>>> trtRuptureLists = new EnumMap<>(TectonicRegionType.class);
		List<List<GriddedRupture>> activeList = new ArrayList<>();
		trtRuptureLists.put(TectonicRegionType.ACTIVE_SHALLOW, activeList);
		List<List<GriddedRupture>> stableList = new ArrayList<>();
		trtRuptureLists.put(TectonicRegionType.STABLE_SHALLOW, stableList);
		
		NSHM23_WUS_FiniteRuptureConverter converter = new NSHM23_WUS_FiniteRuptureConverter();
		GriddedRupturePropertiesCache cache = new GriddedRupturePropertiesCache();
		
		FocalMech mech = FocalMech.STRIKE_SLIP;
		
		for (int l=0; l<origGridList.getNumLocations(); l++) {
			double origSumRate = origGridList.getRuptures(null, l).stream().mapToDouble(R->R.rate).sum();
			List<GriddedRupture> activeRups = origGridList.getRuptures(TectonicRegionType.ACTIVE_SHALLOW, l);
			List<GriddedRupture> stableRups = origGridList.getRuptures(TectonicRegionType.STABLE_SHALLOW, l);
			
			if (!stableRups.isEmpty()) {
				double sumStable = stableRups.stream().mapToDouble(R->R.rate).sum();
				double sumActive = activeRups.stream().mapToDouble(R->R.rate).sum();
				Preconditions.checkState(activeRups.isEmpty(), "Already has both active and stable?"
						+ "\n\tStable:\t%s rups\trate=%s\n\tActive:\t%s rups\trate=%s",
						stableRups.size(), (float)sumStable, activeRups.size(), (float)sumActive);
				activeRups = new ArrayList<>(stableRups.size());
				List<GriddedRupture> modStableRups = new ArrayList<>(stableRups.size());
				
				for (GriddedRupture rup : stableRups) {
					Preconditions.checkState((float)mech.rake() == (float)rup.properties.rake);
					Preconditions.checkState((float)mech.dip() == (float)rup.properties.dip);
					if (stabletoActiveFract > 0d) {
						if (applyStablePropsToOverlapActive) {
							// NSHMP-haz treats the "active" sources in the overlap zone as a GMM-override only, i.e.,
							// it continues to use the stable grid properties with Ztor=5 for all magnitudes
							
							// build it as stable first
							GriddedRupture tempRup = converter.buildFiniteRupture(l, rup.location, rup.properties.magnitude, rup.rate*stabletoActiveFract,
									mech, TectonicRegionType.STABLE_SHALLOW, rup.associatedSections, rup.associatedSectionFracts, cache);
							
							// now convert to active
							GriddedRuptureProperties props = new GriddedRupturePropertiesBuilder(tempRup.properties)
									.tectonicRegionType(TectonicRegionType.ACTIVE_SHALLOW).build();
							props = cache.getCached(props);
							activeRups.add(new GriddedRupture(l, tempRup.location, props, tempRup.rate, tempRup.associatedSections, tempRup.associatedSectionFracts));
						} else {
							// this is probably the right way, but not how Peter does it
							activeRups.add(converter.buildFiniteRupture(l, rup.location, rup.properties.magnitude, rup.rate*stabletoActiveFract,
									mech, TectonicRegionType.ACTIVE_SHALLOW, rup.associatedSections, rup.associatedSectionFracts, cache));
						}
					}
					// this will also correct the zTOR for nshmp-haz that puts everything at 5km (even M>6.5)
					if (stableStaysStableFract > 0d)
						modStableRups.add(converter.buildFiniteRupture(l, rup.location, rup.properties.magnitude, rup.rate*stableStaysStableFract,
								mech, TectonicRegionType.STABLE_SHALLOW, rup.associatedSections, rup.associatedSectionFracts, cache));
				}
				
				stableRups = modStableRups;
			}
			
			activeList.add(activeRups);
			stableList.add(stableRups);
			double newSumRate = activeRups.stream().mapToDouble(R->R.rate).sum() + stableRups.stream().mapToDouble(R->R.rate).sum();
			Preconditions.checkState((float)newSumRate == (float)origSumRate);
		}
		
		return new GridSourceList.Precomputed(origGridList.getGriddedRegion(), trtRuptureLists);
	}

}
