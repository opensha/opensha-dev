package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.InputJumpsOrDistClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;

public class RSQSimFSSConvert {

	public static void main(String[] args) throws IOException {
		File inputFile = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_old/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5.zip");
		File outputFile = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		SectionDistanceAzimuthCalculator distCalc = new SectionDistanceAzimuthCalculator(rupSet.getFaultSectionDataList());
		RuptureConnectionSearch search = new RuptureConnectionSearch(rupSet, distCalc, 100d, false);
		rupSet.addModule(ClusterRuptures.instance(rupSet, search));
		
		HashSet<Jump> jumps = new HashSet<>();
		int maxSplays = 0;
		for (ClusterRupture rup : rupSet.getModule(ClusterRuptures.class)) {
			for (Jump jump : rup.getJumpsIterable()) {
				jumps.add(jump);
				jumps.add(jump.reverse());
			}
			maxSplays = Integer.max(maxSplays, rup.getTotalNumSplays());
		}
		
		ClusterConnectionStrategy connStrat = new InputJumpsOrDistClusterConnectionStrategy(
				distCalc.getSubSections(), distCalc, 15d, jumps);
		PlausibilityConfiguration config = new PlausibilityConfiguration(new ArrayList<>(), maxSplays, connStrat, distCalc);
		rupSet.addModule(config);
		
		sol.write(outputFile);
		
		
	}

}
