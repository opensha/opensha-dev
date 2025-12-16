package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

public class SolGriddedToleranceCompare {

	public static void main(String[] args) throws IOException {
		if (args.length != 4) {
			System.err.println("USAGE: <sol1.zip> <sol2.zip> <abolute-tolerance> <relative-tolerance>");
			System.exit(2);
		}
		File sol1 = new File(args[0]);
		File sol2 = new File(args[1]);
		double absTol = Double.parseDouble(args[2]);
		double relTol = Double.parseDouble(args[3]);
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		System.out.println("Loading GridSourceList1 from "+sol1.getAbsolutePath());
		GridSourceList gridList1 = FaultSystemSolution.load(sol1).requireModule(GridSourceList.class);
		System.out.println("Loading GridSourceList2 from "+sol2.getAbsolutePath());
		GridSourceList gridList2 = FaultSystemSolution.load(sol2).requireModule(GridSourceList.class);
		
		compare(gridList1, gridList2, absTol, relTol);
	}
	
	public static void compare(GridSourceList gridList1, GridSourceList gridList2, double absTol, double relTol) {
		Preconditions.checkState(gridList1.getNumLocations() == gridList2.getNumLocations(), "Location count mismatch");
		Preconditions.checkState(gridList1.getTectonicRegionTypes().equals(gridList2.getTectonicRegionTypes()), "TRTs mismatch");
		
		Set<TectonicRegionType> trts = gridList1.getTectonicRegionTypes();
		
		System.out.println("Validating "+gridList1.getNumLocations()+" location and "+trts.size()+" TRTs");
		System.out.println("\tUsing absTol="+absTol+" and relTol="+relTol);
		for (int l=0; l<gridList1.getNumLocations(); l++) {
			for (TectonicRegionType trt : trts) {
				ImmutableList<GriddedRupture> rups1 = gridList1.getRuptures(trt, l);
				ImmutableList<GriddedRupture> rups2 = gridList2.getRuptures(trt, l);
				if (rups1 != null && !rups1.isEmpty()) {
					Preconditions.checkState(rups2 != null && !rups2.isEmpty(),
							"Have trt=%s rups for %s in 1 but not 2", trt, l);
					
					Preconditions.checkState(rups1.size() == rups2.size(),
							"Rupture count mismatch for trt=%s l=%s: %s != %s", trt, l, rups1.size(), rups2.size());
					for (int r=0; r<rups1.size(); r++) {
						GriddedRupture rup1 = rups1.get(r);
						GriddedRupture rup2 = rups2.get(r);
						if (!rup1.equals(rup2)) {
							String rupsStr = "Location "+l+", "+trt+", rupture "+r
									+"\n\t1: "+gridList1.buildRuptureCSVLine(rup1)+"\n\t2: "+gridList2.buildRuptureCSVLine(rup2);
							// chech within tolerance
							if (!rup1.properties.equals(rup2.properties)) {
								check(rup1.properties.magnitude, rup2.properties.magnitude, "Magnitude", rupsStr, absTol, relTol);
								check(rup1.properties.rake, rup2.properties.rake, "Rake", rupsStr, absTol, relTol);
								check(rup1.properties.dip, rup2.properties.dip, "Dip", rupsStr, absTol, relTol);
								check(rup1.properties.upperDepth, rup2.properties.upperDepth, "Upper Depth", rupsStr, absTol, relTol);
								check(rup1.properties.length, rup2.properties.length, "Length", rupsStr, absTol, relTol);
							}
							check(rup1.rate, rup2.rate, "Rate", rupsStr, absTol, relTol);
						}
						
					}
				} else {
					Preconditions.checkState(rups2 == null || rups2.isEmpty(),
							"Have trt=%s rups for %s in 2 but not 1", trt, l);
				}
			}
		}
		System.out.println("DONE; all valid");
	}
	
	private static void check(double val1, double val2, String property, String rupsStr, double absTol, double relTol) {
		checkAbs(val1, val2, absTol, property, rupsStr);
		checkRel(val1, val2, relTol, property, rupsStr);
	}
	
	private static void checkAbs(double val1, double val2, double absTol, String property, String rupsStr) {
		Preconditions.checkState(Precision.equals(val1, val2, absTol),
				"%s mismatch: %s != %s to absTol=%s\n%s", property, val1, val2, absTol, rupsStr);
	}
	
	private static void checkRel(double val1, double val2, double relTol, String property, String rupsStr) {
		if (Math.abs(val1) < 1e-10) {
			Preconditions.checkState(Math.abs(val1) <= 2e-10);
			return;
		}
		double tol = Math.abs(val1)*relTol;
		Preconditions.checkState(Precision.equals(val1, val2, tol),
				"%s mismatch: %s != %s to relTol=%s->%s\n%s", property, val1, val2, relTol, tol, rupsStr);
	}

}
