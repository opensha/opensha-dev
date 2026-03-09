package scratch.kevin.segmentation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Supplier;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchAlignment;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchLocation;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessDistribution;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;

import com.google.common.base.Preconditions;

public class RelativeCFFCalculator {
	
	private SubSectStiffnessCalculator stiffnessCalc;
	
	private Map<Jump, Double> cffRefRatios;

	public RelativeCFFCalculator(FaultSystemRupSet rupSet, SubSectStiffnessCalculator stiffnessCalc, int threads) {
		this.stiffnessCalc = stiffnessCalc;
		
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		
		// for each jump, find the set of destination clusters on each side; keep those that are the largest or that
		// terminate at subsequent jumps.
		Map<Jump, DestinationClusters> jumpDestinations = new HashMap<>();
		System.out.println("Detecting jump destinations for "+cRups.size()+" ruptures");
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			ClusterRupture rup = cRups.get(rupIndex);
			if (rup.getTotalNumJumps() == 0)
				continue;
			RuptureTreeNavigator nav = rup.getTreeNavigator();
//			for (FaultSubsectionCluster cluster : rup.getClustersIterable()) {
//				Jump jumpTo = nav.getJumpTo(cluster);
//				Jump jumpFrom = nav.get
//			}
			for (Jump origJump : rup.getJumpsIterable()) {
				for (Jump jump : new Jump[] {origJump, origJump.reverse()}) {
//					jump.toCluster;
					boolean intermediate = !nav.getDescendants(jump.toCluster).isEmpty();
					DestinationClusters destinations = jumpDestinations.get(jump);
					if (destinations == null) {
						destinations = new DestinationClusters();
						jumpDestinations.put(jump, destinations);
					}
					destinations.add(jump.toCluster, jump.toSection, intermediate);
				}
			}
		}
		
		if (threads == 1) {
			cffRefRatios = new HashMap<>();
			for (Jump jump : jumpDestinations.keySet()) {
				DestinationClusters sources = jumpDestinations.get(jump.reverse());
				DestinationClusters destinations = jumpDestinations.get(jump);
				double val;
				try {
					val = new JumpRelativeSupplier(jump, sources, destinations).get();
				} catch (Exception e) {
					val = new JumpRelativeSupplier(jump, sources, destinations, true).get();
				}
				cffRefRatios.put(jump, val);
			}
		} else {
			ExecutorService exec = Executors.newFixedThreadPool(threads);
			
			System.out.println("Calculationg CFF reference ratios for "+jumpDestinations.size()+" jumps");
			Map<Jump, CompletableFuture<Double>> futures = new HashMap<>(jumpDestinations.size());
			for (Jump jump : jumpDestinations.keySet()) {
				DestinationClusters sources = jumpDestinations.get(jump.reverse());
				DestinationClusters destinations = jumpDestinations.get(jump);
				futures.put(jump, CompletableFuture.supplyAsync(new JumpRelativeSupplier(jump, sources, destinations), exec));
			}

			cffRefRatios = new HashMap<>();
			for (Jump jump : futures.keySet())
				cffRefRatios.put(jump, futures.get(jump).join());
			
			exec.shutdown();
		}
	}
	
	public double getRelativeCFF(Jump jump) {
		Preconditions.checkState(cffRefRatios.containsKey(jump), "Unexpected jump: "+jump);
		return cffRefRatios.get(jump);
	}
	
	private static class DestinationClusters {
		private FaultSubsectionCluster largestSingle;
		private FaultSubsectionCluster largestForward;
		private FaultSubsectionCluster largestBackward;
		private Set<FaultSubsectionCluster> intermediateClusters;
		
		public DestinationClusters() {
			this.intermediateClusters = new HashSet<>();
		}
		
		public void add(FaultSubsectionCluster toCluster, FaultSection toSect, boolean isIntermediate) {
			int toAlong = toSect.getSubSectionIndexAlong();
			Preconditions.checkState(toAlong >= 0, "Bad indexAlong=%s for %s", toAlong, toSect);
			int minAlong = toAlong;
			int maxAlong = toAlong;
			List<FaultSection> forwardSects = new ArrayList<>();
			List<FaultSection> backwardSects = new ArrayList<>();
			for (FaultSection sect : toCluster.subSects) {
				int along = sect.getSubSectionIndexAlong();
				if (along >= toAlong) {
					forwardSects.add(sect);
					maxAlong = Integer.max(maxAlong, along);
				}
				if (along <= toAlong) {
					backwardSects.add(sect);
					minAlong = Integer.min(minAlong, along);
				}
			}
			if (minAlong == maxAlong) {
				if (largestSingle == null || toCluster.subSects.size() > largestSingle.subSects.size())
					largestSingle = toCluster;
				
				if (isIntermediate)
					intermediateClusters.add(toCluster);
			} else {
				if (!forwardSects.isEmpty()) {
					FaultSubsectionCluster forwardCluster = forwardSects.size() == toCluster.subSects.size() ?
							toCluster : new FaultSubsectionCluster(forwardSects);
					if (largestForward == null || forwardSects.size() > largestForward.subSects.size())
						largestForward = forwardCluster;
					
					if (isIntermediate)
						intermediateClusters.add(forwardCluster);
				}
				
				if (!backwardSects.isEmpty()) {
					FaultSubsectionCluster backwardCluster = backwardSects.size() == toCluster.subSects.size() ?
							toCluster : new FaultSubsectionCluster(backwardSects);
					if (largestBackward == null || backwardSects.size() > largestBackward.subSects.size())
						largestBackward = backwardCluster;
					
					if (isIntermediate)
						intermediateClusters.add(backwardCluster);
				}
			}
		}
		
		public Set<Set<FaultSection>> getUnique() {
			Set<Set<FaultSection>> uniques = new HashSet<>();
			if (intermediateClusters != null)
				for (FaultSubsectionCluster cluster : intermediateClusters)
					uniques.add(Set.copyOf(cluster.subSects));
			
			if (largestForward == null && largestBackward == null && largestSingle != null) {
				// only include single if we don't have a path in either direction
				uniques.add(Set.copyOf(largestSingle.subSects));
			} else {
				if (largestForward != null)
					uniques.add(Set.copyOf(largestForward.subSects));
				if (largestBackward != null)
					uniques.add(Set.copyOf(largestBackward.subSects));
			}
			
			return uniques;
		}
	}
	
	private class JumpRelativeSupplier implements Supplier<Double> {
		
		private Jump jump;
		private DestinationClusters sources;
		private DestinationClusters destinations;
		
		private final boolean D;

		public JumpRelativeSupplier(Jump jump, DestinationClusters sources, DestinationClusters destinations) {
			this(jump, sources, destinations, false);
		}

		public JumpRelativeSupplier(Jump jump, DestinationClusters sources, DestinationClusters destinations, boolean D) {
			this.jump = jump;
			this.sources = sources;
			this.destinations = destinations;
			this.D = D;
		}

		@Override
		public Double get() {
			Set<Set<FaultSection>> sourcesSet = sources.getUnique();
			Set<Set<FaultSection>> destinationsSet = destinations.getUnique();
			
			if (D) System.out.println("Calculating for jump "+jump+" with "+sourcesSet.size()
						+" sourceSets and "+destinationsSet.size()+" destinationSets");
			
			double best = 0d;
			for (Set<FaultSection> sources : sourcesSet) {
				int minAlong = Integer.MAX_VALUE;
				int maxAlong = Integer.MIN_VALUE;
				for (FaultSection sect : sources) {
					int along = sect.getSubSectionIndexAlong();
					minAlong = Integer.min(minAlong, along);
					maxAlong = Integer.max(maxAlong, along);
				}
				int jumpAlong = jump.fromSection.getSubSectionIndexAlong();
				Preconditions.checkState(jumpAlong == minAlong || jumpAlong == maxAlong);
				
				List<FaultSection> jumpEdgeSects = new ArrayList<>(1);
				for (FaultSection sect : sources)
					if (sect.getSubSectionIndexAlong() == jumpAlong)
						jumpEdgeSects.add(sect);
				
				for (Set<FaultSection> destinations : destinationsSet) {
					// this propery accounts for down-dip sections
					double len = FaultSystemRupSet.calculateLength(destinations) * 1e-3; // m -> km
					
					// figure out the best reference CFF value that is of the same length as this destination
					if (D) System.out.println("\tFinding reference for destination of lengh "+(float)len+": "+destinations);
					double bestRef = 0d;
					for (boolean forward : new boolean[] {true,false}) {
						if (minAlong != maxAlong) {
							// see what direction we're going
							boolean actualForward = jumpAlong > minAlong;
							if (actualForward != forward)
								continue;
						}
						if (D) System.out.println("\t\tBuilding source reference of lengh "+(float)len
								+", forward="+forward+" for: "+sources+" with edges: "+jumpEdgeSects);
						// extend the edge in a planar manner in the average strike direction of the subesction
						// to get a reference CFF
						List<FaultSection> jumpToSects = new ArrayList<>(jumpEdgeSects.size());
						for (FaultSection edge : jumpEdgeSects) {
							FaultTrace origTrace = edge.getFaultTrace();
							double direction = forward ? origTrace.getStrikeDirection() : 180d + origTrace.getStrikeDirection();
							Location edgeLoc = forward ? origTrace.last() : origTrace.first();
							GeoJSONFaultSection extension = new GeoJSONFaultSection.Builder(
									100000+edge.getSectionId(), edge.getName()+" extension",
										FaultTrace.of(edgeLoc, LocationUtils.location(edgeLoc, Math.toRadians(direction), len)))
									.upperDepth(edge.getOrigAveUpperDepth())
									.lowerDepth(edge.getAveLowerDepth())
									.dip(edge.getAveDip()).rake(edge.getAveRake())
									.dipDir(edge.getDipDirection())
									.aseismicity(edge.getAseismicSlipFactor())
									.build();
//							if (D) System.out.println("\t\t\tExtension trace: "+extension.getFaultTrace());
							jumpToSects.add(extension);
						}
						double refSum = 0d;
						for (FaultSection source : sources) {
							List<PatchLocation> sourcePatches = stiffnessCalc.getPatches(source);
							for (FaultSection receiver : jumpToSects) {
								List<PatchLocation> receiverPatches = stiffnessCalc.buildPatches(receiver);
								double[] selfStiffness = null;
								if (stiffnessCalc.getSelfStiffnessCap() > 0d)
									selfStiffness = stiffnessCalc.calcSelfStiffness(receiverPatches);
								StiffnessDistribution dist = stiffnessCalc.calcStiffnessDistribution(sourcePatches, receiverPatches, selfStiffness);
								for (double[] vals : dist.get(StiffnessType.CFF))
									for (double val : vals)
										refSum += val;
							}
						}
						if (D) System.out.println("\t\t\tRef CFF="+(float)refSum);
						if (refSum > bestRef)
							bestRef = refSum;
					}
					
					if (D) System.out.println("\t\tBestRef CFF="+(float)bestRef);
					
					double mySum = 0d;
					for (FaultSection source : sources) {
						for (FaultSection receiver : destinations) {
							StiffnessDistribution dist = stiffnessCalc.calcStiffnessDistribution(source, receiver);
							for (double[] vals : dist.get(StiffnessType.CFF))
								for (double val : vals)
									mySum += val;
						}
					}
					
					if (D) System.out.println("\t\tMyCFF="+(float)mySum+";\tBestRefCFF="+(float)bestRef+";\tRatio="+(float)(mySum/bestRef));
					
					Preconditions.checkState(bestRef > 0, "Couldn't find any positive reference CFF value: %s", bestRef);
					best = Math.max(best, mySum/bestRef);
				}
			}
			return best;
		}
		
	}
	
//	private FaultSection buildExtension
	
	public static void main(String[] args) throws IOException {
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
		FaultSystemRupSet rs = factory.buildRuptureSet(NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT, FaultSysTools.defaultNumThreads());
		ClusterConnectionStrategy connStrat = rs.requireModule(PlausibilityConfiguration.class).getConnectionStrategy();
		
		SubSectStiffnessCalculator stiffnessCalc = new SubSectStiffnessCalculator(
				rs.getFaultSectionDataList(), 2d, 3e4, 3e4, 0.5, PatchAlignment.FILL_OVERLAP, 1d);
		
		int threads = 1;
//		int threads = FaultSysTools.defaultNumThreads();
		RelativeCFFCalculator calc = null;
		try {
			calc = new RelativeCFFCalculator(rs, stiffnessCalc, threads);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		
		int fromID = 1424;
		for (Jump jump : connStrat.getJumpsFrom(rs.getFaultSectionData(fromID))) {
			System.out.println(jump);
			System.out.println("\tFROM:\t"+jump.fromCluster);
			System.out.println("\tTO:\t"+jump.toCluster);
			System.out.println("\tForward relCFF:\t"+(float)calc.getRelativeCFF(jump));
			System.out.println("\tReverse relCFF:\t"+(float)calc.getRelativeCFF(jump.reverse()));
		}
	}

}
