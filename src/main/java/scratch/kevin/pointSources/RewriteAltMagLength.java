package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRuptureProperties;
import org.opensha.sha.util.TectonicRegionType;

public class RewriteAltMagLength {
	
	/**
	 * Leonard (2010)
	 * M = a log L + b
	 * log L = (M - b)/a
	 * L = 10^((M - b)/a)
	 */
	public static final Function<GriddedRuptureProperties, Double> LEONARD =
			P -> Math.pow(10.0, (P.magnitude - (P.dip < 90 ? 4.24 : 4.17))/1.67);
	
	/**
	 * WC 94 RLD (subsurface rupture length) formula
	 */
	public static final Function<GriddedRuptureProperties, Double> WC94_RLD =
			P -> Math.pow(10.0,-2.44+0.59*P.magnitude);

	public static void main(String[] args) throws IOException {
		File inputSol = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
				+ "results_WUS_FM_v3_branch_averaged_gridded.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSol);
		
		sol.removeModuleInstances(BranchRegionalMFDs.class);
		sol.removeModuleInstances(BranchSectBVals.class);
		sol.removeModuleInstances(BranchSectNuclMFDs.class);
		sol.removeModuleInstances(BranchSectParticMFDs.class);
		sol.removeModuleInstances(BranchParentSectParticMFDs.class);
		
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
//		Function<GriddedRuptureProperties, Double> magLen = WC94_RLD;
//		File outputSol = new File("/tmp/results_WUS_FM_v3_branch_averaged_gridded_mod_wc_lengths.zip");
		
		Function<GriddedRuptureProperties, Double> magLen = LEONARD;
		File outputSol = new File("/tmp/results_WUS_FM_v3_branch_averaged_gridded_mod_leonard_lengths.zip");
		
		GridSourceList modGridList = new GridSourceList.DynamicallyBuilt(gridList.getTectonicRegionTypes(),
				gridList.getGriddedRegion(), gridList.getRefMFD()) {
			
			@Override
			public int getNumSources() {
				return gridList.getNumSources();
			}
			
			@Override
			public int getLocationIndexForSource(int sourceIndex) {
				return gridList.getLocationIndexForSource(sourceIndex);
			}
			
			@Override
			public TectonicRegionType tectonicRegionTypeForSourceIndex(int sourceIndex) {
				return gridList.tectonicRegionTypeForSourceIndex(sourceIndex);
			}
			
			@Override
			public Set<Integer> getAssociatedGridIndexes(int sectionIndex) {
				return gridList.getAssociatedGridIndexes(sectionIndex);
			}
			
			@Override
			protected List<GriddedRupture> buildRuptures(TectonicRegionType tectonicRegionType, int gridIndex) {
				List<GriddedRupture> orig = gridList.getRuptures(tectonicRegionType, gridIndex);
				if (orig == null || orig.isEmpty())
					return orig;
				List<GriddedRupture> ret = new ArrayList<>(orig.size());
				for (GriddedRupture rup : orig) {
					GriddedRuptureProperties props = rup.properties;
					double mag = props.magnitude;
					double dipRad = Math.toRadians(props.dip);
					double depth = (float)mag <= 6.5f ? 5d : 1d;
					double length = magLen.apply(props);
					double aspectWidth = length / 1.5;
					double ddWidth = (14.0 - depth) / Math.sin(dipRad);
					ddWidth = Math.min(aspectWidth, ddWidth);
					double lower = depth + ddWidth * Math.sin(dipRad);
					props = new GriddedRuptureProperties(mag, props.rake, props.dip,
							props.strike, props.strikeRange, depth, lower, length,
							Double.NaN, Double.NaN, tectonicRegionType);
					rup = new GriddedRupture(rup.gridIndex, rup.location, props, rup.rate, rup.associatedSections, rup.associatedSectionFracts);
					ret.add(rup);
				}
				return ret;
			}
		};
		
		sol.setGridSourceProvider(modGridList);
		
		// don't copy unkonwn files
		sol.write(ArchiveOutput.getDefaultOutput(outputSol, sol.getArchive().getInput()), false);
	}

}
