package scratch.kevin.prvi25.figures;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.util.TectonicRegionType;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;

public class GriddedStatsCalc {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions");
		
		List<File> solFiles = new ArrayList<>();
		List<String> names = new ArrayList<>();
		
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		solFiles.add(new File(dir, "2024_11_19-prvi25_crustal_branches-dmSample5x/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		names.add("Do Nothing");
		
		solFiles.add(new File(dir, "2024_11_20-prvi25_crustal_branches-dmSample5x-grided_rate_balancing/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		names.add("Rate Balancing");
		
		solFiles.add(new File(dir, "2024_11_20-prvi25_crustal_branches-dmSample5x-limit_below_obs_constraint/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		names.add("MFD Constraint");
		
		solFiles.add(new File(dir, "2024_11_20-prvi25_crustal_branches-dmSample5x-limit_below_obs_constraint-grided_rate_balancing/results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		names.add("MFD Constraint + Rate Balancing");
		
		List<FaultSystemSolution> sols = new ArrayList<>();
		for (File solFile : solFiles) {
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			sols.add(sol);
		}
		
		DecimalFormat twoDigits = new DecimalFormat("0.00");
		double[] mags = {6d, 6.5, 7d};
		for (double mag : mags) {
			System.out.println("M>"+(float)mag);
			for (int i=0; i<sols.size(); i++) {
				FaultSystemSolution sol = sols.get(i);
				FaultSystemRupSet rupSet = sol.getRupSet();
				GridSourceList gridProv = sol.requireModule(GridSourceList.class);
				RupMFDsModule mfds = sol.requireModule(RupMFDsModule.class);
				
				double faultRate = 0d;
				for (int r=0; r<rupSet.getNumRuptures(); r++) {
					DiscretizedFunc mfd = mfds.getRuptureMFD(r);
					if (mfd == null) {
						if (rupSet.getMagForRup(r) >= mag)
							faultRate += sol.getRateForRup(r);
					} else {
						for (Point2D pt : mfd)
							if (pt.getX() >= mag)
								faultRate += pt.getY();
					}
				}
				
				double griddedRate = 0d;
				for (int g=0; g<gridProv.getNumLocations(); g++) {
					for (GriddedRupture rup : gridProv.getRuptures(TectonicRegionType.ACTIVE_SHALLOW, g))
						if (rup.properties.magnitude >= mag)
							griddedRate += rup.rate;
				}
				System.out.println("\t"+names.get(i));
				System.out.println("\t\tFault rate:"+(float)faultRate);
				System.out.println("\t\tGridded rate:"+(float)griddedRate);
				System.out.println("\t\t"+twoDigits.format(faultRate/griddedRate)+"x more likely on fault");
			}
		}
		
		FaultSystemSolution subModel = FaultSystemSolution.load(new File(dir,
				"2024_11_19-prvi25_subduction_branches/results_PRVI_SUB_FMs_combined_branch_averaged_gridded.zip"));
		
		GridSourceList gridSources = subModel.requireModule(GridSourceList.class);
		PRVI25_SeismicityRegions[] interfaceRegions = {
				PRVI25_SeismicityRegions.CAR_INTERFACE,
				PRVI25_SeismicityRegions.MUE_INTERFACE
		};
		for (PRVI25_SeismicityRegions seisReg : interfaceRegions) {
			Region reg = seisReg.load();
			
			double faultTotalRate = 0d;
			for (FaultSection sect : subModel.getRupSet().getFaultSectionDataList()) {
				boolean muertos = sect.getSectionName().contains("Muertos");
				if (seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE) {
					if (muertos)
						faultTotalRate += subModel.calcNucleationRateForSect(sect.getSectionId(), 0d, 10d);
				} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE) {
					if (!muertos)
						faultTotalRate += subModel.calcNucleationRateForSect(sect.getSectionId(), 0d, 10d);
				}
			}
			
			double griddedRateM5 = 0d;
			double griddedRateM6 = 0d;
			for (int i=0; i<gridSources.getNumLocations(); i++) {
				Location loc = gridSources.getLocation(i);
				for (GriddedRupture rup : gridSources.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, i)) {
					if (reg.contains(loc)) {
						if (rup.properties.magnitude >= 5d)
							griddedRateM5 += rup.rate;
						if (rup.properties.magnitude >= 6d)
							griddedRateM6 += rup.rate;
					}
				}
			}
			
			DecimalFormat pDF = new DecimalFormat("0.##%");
			System.out.println(seisReg);
			System.out.println("\tFault rate: "+(float)faultTotalRate);
			System.out.println("\tGridded rate >M5: "+(float)griddedRateM5+";\tfaultFract="+pDF.format(faultTotalRate/griddedRateM5));
			System.out.println("\tGridded rate >M6: "+(float)griddedRateM6+";\tfaultFract="+pDF.format(faultTotalRate/griddedRateM6));
		}
	}

}
