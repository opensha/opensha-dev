package scratch.kevin.nshm23.segModelTests;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.util.IDPairing;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;

public class RSQSimComparableConnStratAttach {

	public static void main(String[] args) throws IOException {
		File inputRSQSimSol = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(inputRSQSimSol);
		File outputSolFile = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5_plus_coulomb_conns.zip");
		
		FaultSystemRupSet refRupSet = FaultSystemRupSet.load(
				new File("/home/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
		ClusterConnectionStrategy refConnStrat = refRupSet.requireModule(PlausibilityConfiguration.class).getConnectionStrategy();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		ClusterRuptures inputRups = rupSet.getModule(ClusterRuptures.class);
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(rupSet.getFaultSectionDataList());
		if (inputRups == null) {
			RuptureConnectionSearch connSearch = new RuptureConnectionSearch(rupSet, distAzCalc);
			inputRups = ClusterRuptures.instance(rupSet, connSearch);
			
		}
		HashSet<Jump> inputJumps = new HashSet<>();
		for (ClusterRupture rup : inputRups)
			for (Jump jump : rup.getJumpsIterable())
				inputJumps.add(jump);
		InputJumpsOrExternalClusterConnectionStrategy newConnStrat = new InputJumpsOrExternalClusterConnectionStrategy(
				rupSet.getFaultSectionDataList(), distAzCalc, inputJumps, refConnStrat);
		
		PlausibilityConfiguration config = new PlausibilityConfiguration(null, 0, newConnStrat, distAzCalc);
		
		rupSet.addModule(config);
		sol.write(outputSolFile);
	}
	
	public static class InputJumpsOrExternalClusterConnectionStrategy extends ClusterConnectionStrategy {

		private SectionDistanceAzimuthCalculator distCalc;
		private double maxJumpDist;
		private HashSet<IDPairing> allowedSectConnections;

		public InputJumpsOrExternalClusterConnectionStrategy(List<? extends FaultSection> subSects,
				SectionDistanceAzimuthCalculator distCalc, Collection<Jump> inputJumps,
				ClusterConnectionStrategy externalConnStrat) {
			super(subSects, distCalc);
			this.maxJumpDist = externalConnStrat.getMaxJumpDist();
			this.distCalc = distCalc;
			initAllowedSectConnections(inputJumps, externalConnStrat);
		}
		
		private void initAllowedSectConnections(Collection<Jump> jumps, ClusterConnectionStrategy externalConnStrat) {
			allowedSectConnections = new HashSet<>();
			HashSet<IDPairing> allowedParentConnections = new HashSet<>();
			for (Jump jump : jumps) {
				IDPairing pair = new IDPairing(jump.fromSection.getSectionId(), jump.toSection.getSectionId());
				allowedSectConnections.add(pair);
				allowedSectConnections.add(pair.getReversed());
				IDPairing parentPair = new IDPairing(jump.fromSection.getParentSectionId(), jump.toSection.getParentSectionId());
				allowedParentConnections.add(parentPair);
				allowedParentConnections.add(parentPair.getReversed());
			}
			// now add any external jumps from non-connected parents
			for (Jump jump : externalConnStrat.getAllPossibleJumps()) {
				IDPairing parentPair = new IDPairing(jump.fromSection.getParentSectionId(), jump.toSection.getParentSectionId());
				if (!allowedParentConnections.contains(parentPair)) {
					// connections between these faults are not in the input catalog, add external
					IDPairing pair = new IDPairing(jump.fromSection.getSectionId(), jump.toSection.getSectionId());
					allowedSectConnections.add(pair);
					allowedSectConnections.add(pair.getReversed());
				}
			}
		}

		@Override
		protected List<Jump> buildPossibleConnections(FaultSubsectionCluster from, FaultSubsectionCluster to) {
			List<Jump> ret = new ArrayList<>();
			for (FaultSection s1 : from.subSects) {
				for (FaultSection s2 : to.subSects) {
					double dist = distCalc.getDistance(s1, s2);
					if (allowedSectConnections.contains(new IDPairing(s1.getSectionId(), s2.getSectionId())))
						// it's an input jump
						ret.add(new Jump(s1, from, s2, to, dist));
				}
			}
			if (ret.isEmpty())
				return null;
			return ret;
		}

		@Override
		public String getName() {
			return "InputPlusDist: maxDist="+(float)maxJumpDist+" km";
		}

		@Override
		public double getMaxJumpDist() {
			return maxJumpDist;
		}

	}

}
