package scratch.kevin.nshm23.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.PlausibilityFilterPlot;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.PlausibilityFilterPlot.RupSetPlausibilityResult;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityResult;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.DirectPathPlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.MinSectsPerParentFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.path.NucleationClusterEvaluator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.path.PathPlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.SectCountAdaptiveRuptureGrowingStrategy.ConnPointCleanupFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupCartoonGenerator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class ObsRupFiniteFaultPlausibilityComparison {

	public static void main(String[] args) throws IOException {
		
		File outputDir = new File("/tmp/obs_rup_compare");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<String> rupNames = new ArrayList<>();
		List<int[]> rupSectIDs = new ArrayList<>();
		
//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
//		
//		rupNames.add("Landers");
////		rupSectIDs.add(new int[] {243, 242, 557, 556, 555, 880, 879, 878, 1025, 1024, 991, 990, 989, 560, 559, 175});
//		rupSectIDs.add(new int[] {243, 242, 557, 556, 555, 880, 879, 878, 1025, 1024, 991, 990, 989, 560, 559});
//		
//		rupNames.add("Hector Mine");
//		rupSectIDs.add(new int[] {1583, 1584, 1585, 846, 845, 844, 843});
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
				+ "NSHM23_v1p4_plausibleMulti15km_adaptive6km_direct_cmlRake360_jumpP0.001_slipP0.05incrCapDist_"
				+ "cff0.75IntsPos_comb2Paths_cffFavP0.01_cffFavRatioN2P0.5_sectFractGrow0.1.zip"));
		rupNames.add("Ridgecrest M6.4");
//		rupSectIDs.add(new int[] {3794, 3793, 3792, 3791, 3283, 3284});
		rupSectIDs.add(new int[] {3794, 3793, 3792, 3283, 3284});
		
		PlausibilityConfiguration config = rupSet.requireModule(PlausibilityConfiguration.class);
		
		RuptureConnectionSearch connSearch = new RuptureConnectionSearch(rupSet, config.getDistAzCalc(), 15d, false);
		
//		System.out.println(config.getDistAzCalc().getDistance(989, 879));
//		System.exit(0);
		
		ExecutorService exec = Executors.newSingleThreadExecutor();
		
		List<PlausibilityFilter> inputFilters = new ArrayList<>();
		for (PlausibilityFilter filter : config.getFilters()) {
			if (filter instanceof MinSectsPerParentFilter)
				continue;
			if (filter instanceof ConnPointCleanupFilter)
				continue;
			if (filter instanceof DirectPathPlausibilityFilter)
				continue;
			inputFilters.add(filter);
		}
		
		for (int r=0; r<rupNames.size(); r++) {
			String rupName = rupNames.get(r);
			int[] sectIDs = rupSectIDs.get(r);
			
			System.out.println("Working on rupture: "+rupName);
			
			List<FaultSection> sects = new ArrayList<>();
			for (int sectID : sectIDs)
				sects.add(rupSet.getFaultSectionData(sectID));
			List<FaultSubsectionCluster> clusters = connSearch.calcClusters(sects, false);
			List<Jump> jumps = connSearch.calcRuptureJumps(clusters, false);
//			System.exit(0);
			ClusterRupture rup = connSearch.buildClusterRupture(clusters, jumps, false);
			
			double maxJumpDist = 0d;
			for (Jump jump : rup.getJumpsIterable())
				maxJumpDist = Math.max(maxJumpDist, jump.distance);
			System.out.println("Maximum jump distance: "+(float)maxJumpDist);
			
			RupCartoonGenerator.plotRupture(outputDir, rupName, rup, rupName, false, true);
			
			RupSetPlausibilityResult results = PlausibilityFilterPlot.testRupSetPlausibility(
					List.of(rup), inputFilters, config, connSearch, exec);
			
			List<PlausibilityFilter> filters = results.filters;
			
			PlausibilityResult netResult = PlausibilityResult.PASS;
			
			for (int f=0; f<filters.size(); f++) {
				PlausibilityFilter filter  = filters.get(f);
				PlausibilityResult result = results.filterResults.get(f).get(0);
				
				netResult = netResult.logicalAnd(result);
				
				if (results.scalarVals.get(f) != null) {
					Double scalar = results.scalarVals.get(f).get(0);
					System.out.println("\t"+filter.getName()+": "+result+"  ("+scalar+")");
				} else {
					System.out.println("\t"+filter.getName()+": "+result);
				}
//				if (!result.isPass())
//					filter.apply(rup, true);
				if (filter instanceof PathPlausibilityFilter && ((PathPlausibilityFilter)filter).getEvaluators().length > 1) {
					NucleationClusterEvaluator[] pathEvals = ((PathPlausibilityFilter)filter).getEvaluators();
					for (NucleationClusterEvaluator eval : pathEvals) {
						if (eval instanceof NucleationClusterEvaluator.Scalar<?>) {
							PathPlausibilityFilter.Scalar<?> indvFilter = new PathPlausibilityFilter.Scalar<>(
									(NucleationClusterEvaluator.Scalar<?>)eval);
							PlausibilityResult indvResult = indvFilter.apply(rup, false);
							System.out.println("\t\t"+eval.getName()+": "+indvResult+" ("+indvFilter.getValue(rup)+")");
						} else {
							PathPlausibilityFilter indvFilter = new PathPlausibilityFilter(eval);
							PlausibilityResult indvResult = indvFilter.apply(rup, false);
							System.out.println("\t\t"+eval.getName()+": "+indvResult);
						}
					}
				}
			}
			System.out.println("Net result: "+netResult);
		}
		
		exec.shutdown();
	}

}
