package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider.AbstractPrecomputed;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_AbstractGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_GridFocalMechs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

public class GridUpdateFocalRegions {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File solFileGridded = new File(dir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		File solFile = new File(dir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		
		FaultSystemSolution solGridded = FaultSystemSolution.load(solFileGridded);
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		// first do tectonic regimes
		Region stableReg = NSHM23_RegionLoader.GridSystemRegions.CEUS_STABLE.load();
		Map<Region, TectonicRegionType> regRegimes = Map.of(stableReg, TectonicRegionType.STABLE_SHALLOW);
		RupSetTectonicRegimes regimes = RupSetTectonicRegimes.forRegions(
				sol.getRupSet(), regRegimes, TectonicRegionType.ACTIVE_SHALLOW, 0.5);
		
		solGridded.getRupSet().addModule(regimes);
		sol.getRupSet().addModule(regimes);
		
		// now do focal mechanisms
		AbstractPrecomputed gridProv = (AbstractPrecomputed)solGridded.getGridSourceProvider();
		System.out.println("Grid prov is of type: "+gridProv.getClass().getName());
		
		SeismicityRegions region = SeismicityRegions.CONUS_WEST;
		GriddedRegion gridReg = gridProv.getGriddedRegion();
		
		double[] fractStrikeSlip = NSHM23_GridFocalMechs.getFractStrikeSlip(region, gridReg);
		double[] fractReverse = NSHM23_GridFocalMechs.getFractReverse(region, gridReg);
		double[] fractNormal = NSHM23_GridFocalMechs.getFractNormal(region, gridReg);
		
		NSHM23_AbstractGridSourceProvider.Precomputed updated = new NSHM23_AbstractGridSourceProvider.Precomputed(
				gridReg, gridProv.getNodeSubSeisMFDs(), gridProv.getNodeUnassociatedMFDs(),
				fractStrikeSlip, fractNormal, fractReverse);
		
		solGridded.setGridSourceProvider(updated);
		solGridded.write(new File(dir, "updated_"+solFileGridded.getName()));
		sol.write(new File(dir, "updated_"+solFile.getName()));
	}

}
