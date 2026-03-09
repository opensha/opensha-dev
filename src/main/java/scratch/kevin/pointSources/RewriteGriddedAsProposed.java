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

public class RewriteGriddedAsProposed {

	public static void main(String[] args) throws IOException {
		File inputSol = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
				+ "results_WUS_FM_v3_branch_averaged_gridded.zip");
		File outputSol = new File("/tmp/results_WUS_FM_v3_branch_averaged_gridded_mod_proposed.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSol);
		
		sol.removeModuleInstances(BranchRegionalMFDs.class);
		sol.removeModuleInstances(BranchSectBVals.class);
		sol.removeModuleInstances(BranchSectNuclMFDs.class);
		sol.removeModuleInstances(BranchSectParticMFDs.class);
		sol.removeModuleInstances(BranchParentSectParticMFDs.class);
		
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
		Function<GriddedRuptureProperties, Double> magLen = RewriteAltMagLength.LEONARD;
		
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
					double length = magLen.apply(props);
					double[] depths = RewriteMagDepthProps.calcDepths(mag, length, dipRad);
					double upper = depths[0];
					double lower = depths[1];
					props = new GriddedRuptureProperties(mag, props.rake, props.dip,
							props.strike, props.strikeRange, upper, lower, length,
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
