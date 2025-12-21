package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.util.Interpolate;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.util.TectonicRegionType;

public class RewriteMagDepthProps {
	
	private static WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();

	public static void main(String[] args) throws IOException {
		File inputSol = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip");
		File outputSol = new File("/tmp/results_WUS_FM_v3_branch_averaged_gridded_mod_grid_depths.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSol);
		
		sol.removeModuleInstances(BranchRegionalMFDs.class);
		sol.removeModuleInstances(BranchSectBVals.class);
		sol.removeModuleInstances(BranchSectNuclMFDs.class);
		sol.removeModuleInstances(BranchSectParticMFDs.class);
		sol.removeModuleInstances(BranchParentSectParticMFDs.class);
		
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
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
					double length = WC94.getMedianLength(mag);
					double[] depths = calcDepths(mag, length, dipRad);
					double upper = depths[0];
					double lower = depths[1];
					props = new GriddedRuptureProperties(mag, props.rake, props.dip,
							props.strike, props.strikeRange, upper, lower, props.length,
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
	private static final double zCenter = 5d;	// km
	private static final double zTopMin = 1d;	// km
	private static final double zBotMax = 14d;	// km
	private static final double maxVertThickness = zBotMax - zTopMin;	// 13 km
	
	public static double[] calcDepths(double mag, double length, double dipRad) {
		// Same aspect control as before
		double aspectWidthDD = length / 1.5;	// down-dip width (km)

		// Cap by maximum possible down-dip width given seismogenic thickness
		double sinDip = Math.sin(dipRad);
		double maxWidthDD = maxVertThickness / sinDip;
		double ddWidth = Math.min(aspectWidthDD, maxWidthDD);

		// Convert to vertical thickness (km)
		double vertThickness = ddWidth * sinDip;

		// Center vertically at 5 km
		double upper = zCenter - 0.5 * vertThickness;
		double lower = zCenter + 0.5 * vertThickness;

		// Enforce saturation bounds by shifting the interval if needed.
		// (Keep thickness fixed unless it exceeds the max thickness.)
		if (vertThickness >= maxVertThickness) {
			upper = zTopMin;
			lower = zBotMax;
		} else if (upper < zTopMin) {
			upper = zTopMin;
			lower = upper + vertThickness;
		} else if (lower > zBotMax) {
			lower = zBotMax;
			upper = lower - vertThickness;
		}
		return new double[] {upper, lower};
	}

}
