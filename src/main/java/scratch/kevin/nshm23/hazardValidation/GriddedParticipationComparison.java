package scratch.kevin.nshm23.hazardValidation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.BitSet;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.SourceTree;

public class GriddedParticipationComparison {

	public static void main(String[] args) throws IOException {
		File inputSolFile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
		double[] minMags = {0d, 6d, 7d};
		
		File outputDir = new File("/data/kevin/nshm23/nshmp-haz-models/gridded_partic_debug");
		
		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		// we want direct mappings
		sol.getRupSet().removeModuleInstances(FaultGridAssociations.class);
		
		HazardModel model = HazardModel.load(modelDir.toPath());
		
		GriddedGeoDataSet[] fssXYZs = new GriddedGeoDataSet[minMags.length];
		GriddedGeoDataSet[] nhXYZs = new GriddedGeoDataSet[minMags.length];
		for (int m=0; m<minMags.length; m++) {
			fssXYZs[m] = new GriddedGeoDataSet(gridReg);
			nhXYZs[m] = new GriddedGeoDataSet(gridReg);
		}
		
		// do FSS
		if (bgOp == IncludeBackgroundOption.INCLUDE || bgOp == IncludeBackgroundOption.EXCLUDE) {
			FaultSystemRupSet rupSet = sol.getRupSet();
			BitSet[] sectBitSets = new BitSet[rupSet.getNumSections()];
			for (int s=0; s<rupSet.getNumSections(); s++) {
				FaultSection sect = rupSet.getFaultSectionData(s);
				sectBitSets[s] = new BitSet(gridReg.getNodeCount());
				for (Location loc : sect.getFaultSurface(1d).getEvenlyDiscritizedListOfLocsOnSurface()) {
					int node = gridReg.indexForLocation(loc);
					if (node >= 0)
						sectBitSets[s].set(node);
				}
			}
			
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				BitSet rupBits = null;
				for (int s : rupSet.getSectionsIndicesForRup(r)) {
					if (rupBits == null)
						rupBits = (BitSet)sectBitSets[s].clone();
					else
						rupBits.or(sectBitSets[s]);
				}
				double rate = sol.getRateForRup(r);
				double mag = rupSet.getMagForRup(r);
				for (int m=0; m<minMags.length; m++)
					if ((float)mag >= (float)minMags[m])
						for (int i = rupBits.nextSetBit(0); i >= 0; i = rupBits.nextSetBit(i + 1))
							fssXYZs[m].add(i, rate);
			}
		}
		if (bgOp == IncludeBackgroundOption.INCLUDE || bgOp == IncludeBackgroundOption.ONLY) {
			GridSourceList gridList = sol.requireModule(GridSourceList.class);
			
			for (int l=0; l<gridList.getNumLocations(); l++) {
				Location loc = gridList.getLocation(l);
				int locIndex = gridReg.indexForLocation(loc);
				if (locIndex >= 0) {
					for (GriddedRupture rup : gridList.getRuptures(trt, l)) {
						for (int m=0; m<minMags.length; m++) {
							if ((float)rup.properties.magnitude >= (float)minMags[m])
								fssXYZs[m].add(locIndex, rup.rate);
						}
					}
				}
			}
		}
		
		// do NSHMP-haz
		for (SourceTree tree : model) {
			// TODO
		}
	}
	
	

}
