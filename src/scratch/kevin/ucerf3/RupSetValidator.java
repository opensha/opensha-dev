package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.SectionCluster;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester.TestType;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class RupSetValidator {
	
	private static final boolean D = true;
	
	private static boolean checkSingleParentRups(
			int parent, int numSects, Set<Integer> rups, FaultSystemRupSet rupSet) {
		int expected = (numSects + 1) * numSects / 2 - numSects;
		int actual = rups.size();
		if (actual != expected) {
			if (D) {
				System.out.println("Single Parent Rup Test Failed! "
						+parent+". expected: "+expected+"\tactual: "+actual);
				for (FaultSectionPrefData data : rupSet.getFaultSectionDataList())
					if (data.getParentSectionId() == parent) {
						System.out.println(data.getSectionName());
						break;
					}
			}
			return false;
		}
		return true;
	}
	
	private static boolean checkAllSingleParentRups(FaultSystemRupSet rupSet) {
		int prevParent = -1;
		int sectsForParent = 0;
		HashSet<Integer> rupsForParent = null;
//		System.out.println("Starting the test!");
		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			FaultSectionPrefData data = rupSet.getFaultSectionData(sectIndex);
			if (rupsForParent == null || prevParent != data.getParentSectionId()) {
				if (rupsForParent != null) {
					// do the test
					if (!checkSingleParentRups(prevParent, sectsForParent, rupsForParent, rupSet))
						return false;
				}
				sectsForParent = 0;
				rupsForParent = new HashSet<Integer>();
				prevParent = data.getParentSectionId();
			}
			
			sectsForParent++;
			for (int r : rupSet.getRupturesForSection(sectIndex))
				if (!isMultiFault(r, rupSet) && !rupsForParent.contains(r))
					rupsForParent.add(r);
		}
		if (!checkSingleParentRups(prevParent, sectsForParent, rupsForParent, rupSet))
			return false;
		return true;
	}
	
	private static boolean isMultiFault(int rup, FaultSystemRupSet rupSet) {
		int parentID = -1;
		for (int sectIndex : rupSet.getSectionsIndicesForRup(rup)) {
			FaultSectionPrefData data = rupSet.getFaultSectionData(sectIndex);
			if (parentID < 0)
				parentID = data.getParentSectionId();
			else if (parentID != data.getParentSectionId())
				return true;
		}
		return false;
	}
	
	private static boolean doMultiFaultParentChecks(InversionFaultSystemRupSet rupSet) {
		FaultModels fm = rupSet.getFaultModel();
		
		// now check to make sure some faults rupture together
		ArrayList<ArrayList<Integer>> parentChecks = new ArrayList<ArrayList<Integer>>();
		// this checks SAF wall to wall the hard way through San Gorgonio Pass
		parentChecks.add(Lists.newArrayList(653, 283, 284, 295));
		// this checks Landers
		parentChecks.add(Lists.newArrayList(91, 92, 723, 90, 724));
		// SAF to San Jacinto
		parentChecks.add(Lists.newArrayList(286, 293));
		// All of San Jacinto
		parentChecks.add(Lists.newArrayList(119, 28));
		// Whittier to end of Elsinore
		if (fm == FaultModels.FM3_1)
			parentChecks.add(Lists.newArrayList(236, 102));
		else
			parentChecks.add(Lists.newArrayList(237, 102));
		// all of Death Valley
		parentChecks.add(Lists.newArrayList(45, 246));
		// all of Garlock
		parentChecks.add(Lists.newArrayList(49, 341, 48));
		
		ArrayList<ArrayList<Integer>> rupParentsList = new ArrayList<ArrayList<Integer>>();
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			ArrayList<Integer> parents = new ArrayList<Integer>();
			for (int s : rupSet.getSectionsIndicesForRup(r)) {
				int parent = rupSet.getFaultSectionData(s).getParentSectionId();
				if (!parents.contains(parent))
					parents.add(parent);
			}
			Preconditions.checkState(parents.size() >= 1);
			rupParentsList.add(parents);
		}
		
		for (ArrayList<Integer> parentCheck : parentChecks) {
			boolean found = false;
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				ArrayList<Integer> rupParents = rupParentsList.get(r);
				boolean match = true;
				for (int expectedParent : parentCheck) {
					if (!rupParents.contains(expectedParent)) {
						match = false;
						break;
					}
				}
				if (match) {
					found = true;
					break;
				}
			}
			if (!found) {
				if (D) {
					System.out.println("Failed parent check for: "+Joiner.on(", ").join(parentCheck));
				}
				return false;
			}
		}
		
		return true;
	}
	
	private static boolean doMultiFaultSubSectChecks(InversionFaultSystemRupSet rupSet) {
		ArrayList<ArrayList<Integer>> subSectChecks = new ArrayList<ArrayList<Integer>>();
		
		switch (rupSet.getFaultModel()) {
		case FM3_1:
			// full length SAF the hard way
			subSectChecks.add(Lists.newArrayList(1893, 1943, 1791));
//			subSectChecks.add(Lists.newArrayList(1900, 1943, 165));
			
			// full length SAF the easy way
			subSectChecks.add(Lists.newArrayList(1900, 1840, 1791));
//			subSectChecks.add(Lists.newArrayList(1900, 1840, 165));
			
			// full length Garlock
			subSectChecks.add(Lists.newArrayList(621, 620));
			
			// full length Elsinore
			subSectChecks.add(Lists.newArrayList(2586, 510));
			break;
		case FM3_2:
			// full length SAF the hard way
			subSectChecks.add(Lists.newArrayList(1949, 1999, 1847));
			
			// full length SAF the easy way
			subSectChecks.add(Lists.newArrayList(1956, 1896, 1847));
			
			// full length Garlock
			subSectChecks.add(Lists.newArrayList(643, 642));
			
			// full length Elsinore
			subSectChecks.add(Lists.newArrayList(2642, 535));
			break;

		default:
			throw new IllegalStateException("Unkown FM: "+rupSet.getFaultModel());
		}
		
		for (ArrayList<Integer> subSectCheck : subSectChecks) {
			boolean found = false;
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				boolean match = true;
				List<Integer> rupSects = rupSet.getSectionsIndicesForRup(r);
				for (int expectedSubSect : subSectCheck) {
					if (!rupSects.contains(expectedSubSect)) {
						match = false;
						break;
					}
				}
				if (match) {
					found = true;
					break;
				}
			}
			if (!found) {
				if (D) {
					System.out.println("Failed subSect check for: "+Joiner.on(", ").join(subSectCheck));
				}
				return false;
			}
		}
		
		return true;
	}
	
	private static boolean validateRupSet(InversionFaultSystemRupSet rupSet) {
		// first do the single parent rup purmutations check
		if (!checkAllSingleParentRups(rupSet))
			return false;
		
		if (!doMultiFaultParentChecks(rupSet))
			return false;
		
		if (!doMultiFaultSubSectChecks(rupSet))
			return false;
		
		return true;
	}
	
	private static ArrayList<LaughTestFilter> generateFilters() {
		ArrayList<LaughTestFilter> filters = new ArrayList<LaughTestFilter>();
		
//		double maxJumpDist = 5;
//		double[] maxAzChanges = { 45, 75, 90 };
//		double[] maxTotAzChanges = { 45, 75, 90 };
//		double maxRakeDiff = Double.POSITIVE_INFINITY;
////		double[] maxCmlJumpDists = { 5, 7.5, 10 };
//		double[] maxCmlJumpDists = { 10 };
//		int minNumSectInRup = 2;
//		double[] maxCmlRakeChanges = { 270, 360, 450 };
//		double[] maxCmlAzChanges = { 360, 450, 540 };
//		
//		double[] minAverageProbs = { 0.05, 0.1, 0.2 };
//		double[] minIndividualProbs = { 0.02, 0.05, 0.075, 0.1 };
////		double minimumStressExclusionCeiling = 1d;
//		double[] minimumStressExclusionCeilings = { Double.POSITIVE_INFINITY, 1d, 2.5d };
		
		
		double maxJumpDist = 5;
//		double[] maxAzChanges = { 45, 60, 75, 90 };
//		double[] maxAzChanges = { 60 };
		double[] maxAzChanges = { 60 };
//		double[] maxAzChanges = { 90 };
//		double[] maxTotAzChanges = { 45, 60, 75, 90 };
//		double[] maxTotAzChanges = { 90 };
		double[] maxTotAzChanges = { 60 };
//		double[] maxCmlJumpDists = { 5, 7.5, 10 };
//		double[] maxCmlJumpDists = { 10 };
		double[] maxCmlJumpDists = { 5 };
		int minNumSectInRup = 2;
//		double[] maxCmlRakeChanges = { 180, 225, 270, 315, 360, 405, 450 };
//		double[] maxCmlRakeChanges = { 450 };
		double[] maxCmlRakeChanges = { 180 };
//		double[] maxCmlAzChanges = { 360, 450, 540 };
//		double[] maxCmlAzChanges = { 630 };
//		double[] maxCmlAzChanges = { 540, 550, 560 };
		double[] maxCmlAzChanges = { 560 };
		
//		double[] minAverageProbs = { 0.05, 0.1, 0.2 };
//		double[] minAverageProbs = { 0.02, 0.05, 0.075, 0.1, 0.15, 0.2 };
//		double[] minAverageProbs = { 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275 };
		double[] minAverageProbs = { 0.1 };
//		double[] minAverageProbs = { 0.4 };
//		double[] minIndividualProbs = { 0.02, 0.035, 0.05, 0.075, 0.1 };
		double[] minIndividualProbs = { 0.1 };
//		double[] minIndividualProbs = { 0.05, 0.075, 0.1, 0.125, 0.15, 0.175 };
//		double minimumStressExclusionCeiling = 1d;
//		double[] minimumStressExclusionCeilings = { 1, 1.25, 1.5, 1.75, 2, 2.15  };
		double[] minimumStressExclusionCeilings = { 1.5 };
		boolean applyBranchesOnly = true;
		
		for (double maxAzimuthChange : maxAzChanges) {
			for (double maxTotAzimuthChange : maxTotAzChanges) {
				for (double maxCumJumpDist : maxCmlJumpDists) {
					for (double maxCmlRakeChange : maxCmlRakeChanges) {
						for (double maxCmlAzimuthChange : maxCmlAzChanges) {
							for (double minAverageProb : minAverageProbs) {
								for (double minIndividualProb : minIndividualProbs) {
									for (double minimumStressExclusionCeiling : minimumStressExclusionCeilings) {
										if (minIndividualProb > minAverageProb)
											continue;
										CoulombRatesTester coulombFilter =
											new CoulombRatesTester(TestType.COULOMB_STRESS, minAverageProb,
													minIndividualProb, minimumStressExclusionCeiling, applyBranchesOnly, true);
										filters.add(new LaughTestFilter(maxJumpDist, maxAzimuthChange,
												maxTotAzimuthChange, maxCumJumpDist,
												maxCmlRakeChange, maxCmlAzimuthChange, minNumSectInRup,
												false, coulombFilter));
									}
								}
							}
						}
					}
				}
			}
		}
		
		return filters;
	}
	
	private static class ValidationTask implements Task {
		
		private FaultModels faultModel;
		private DeformationModels defModel;
		private LaughTestFilter filter;
		private List<FaultSectionPrefData> faultSectionData;
		private Map<IDPairing, Double> subSectionDistances;
		private Map<IDPairing, Double> subSectionAzimuths;
		private int numRups;
		private Boolean passes = null;
		
		public ValidationTask(FaultModels faultModel, DeformationModels defModel, LaughTestFilter filter,
				List<FaultSectionPrefData> faultSectionData, Map<IDPairing, Double> subSectionDistances, Map<IDPairing, Double> subSectionAzimuths) {
			this.faultModel = faultModel;
			this.defModel = defModel;
			this.filter = filter;
			this.faultSectionData = faultSectionData;
			this.subSectionDistances = subSectionDistances;
			this.subSectionAzimuths = subSectionAzimuths;
		}

		@Override
		public void compute() {
			System.out.println("Building a rup set.");
			SectionClusterList clusters = new SectionClusterList(faultModel, defModel, filter,
					faultSectionData, subSectionDistances, subSectionAzimuths);
			InversionFaultSystemRupSet rupSet = new FakeRupSet(clusters);
			numRups = rupSet.getNumRuptures();
			System.out.println("Done building a rup set.");
			if (numRups < currentMin) {
				passes = validateRupSet(rupSet);
				if (passes)
					updateCurrentMin(numRups, filter);
				System.out.println("Checked for "+numRups+" rups. passed? "+passes);
			} else
				System.out.println("Didn't bother checking for "+numRups+" rups.");
		}
		
	}
	
	private static int currentMin = Integer.MAX_VALUE;
	private static synchronized void updateCurrentMin(int rups, LaughTestFilter filter) {
		if (rups < currentMin) {
			currentMin = rups;
			System.out.println("We have a new min ("+currentMin+"): "+filter);
		}
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		try {
			ArrayList<LaughTestFilter> filters = generateFilters();
			Collections.shuffle(filters);
			System.out.println("Num filters: "+filters.size());
			
			LaughTestFilter best = null;
			FaultModels faultModel = FaultModels.FM3_2;
			DeformationModels defModel = DeformationModels.GEOLOGIC;
			DeformationModelFetcher fetch = new DeformationModelFetcher(faultModel,
					defModel, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE);
			Map<IDPairing, Double> subSectionDistances =
				fetch.getSubSectionDistanceMap(filters.get(0).getMaxJumpDist());
			Map<IDPairing, Double> subSectionAzimuths =
				fetch.getSubSectionAzimuthMap(subSectionDistances.keySet());
			ArrayList<FaultSectionPrefData> faultSectionData = fetch.getSubSectionList();
			
			ArrayList<ValidationTask> tasks = new ArrayList<RupSetValidator.ValidationTask>();
			for (LaughTestFilter filter : filters) {
				tasks.add(new ValidationTask(faultModel, defModel, filter,
						faultSectionData, subSectionDistances, subSectionAzimuths));
			}
			
			ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks);
			comp.computeThreaded(6);
			
			int minRups = Integer.MAX_VALUE;
			int cnt = 0;
			for (ValidationTask task : tasks) {
				System.out.println((cnt++)+". "+task.numRups+" ? "
						+task.passes+": "+task.filter);
				// == true here because passes is a Boolean that can be null (if it wasn't checked due to size)
				if (task.passes != null && task.passes && task.numRups < minRups) {
					minRups = task.numRups;
					best = task.filter;
				}
			}
			
			System.out.println("Best found is "+minRups+" rups: "+best);
			
//			int minRups = Integer.MAX_VALUE;
//			int numChecked = 0;
//			for (LaughTestFilter filter : filters) {
////			FaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(FaultModels.FM3_1,
////					DeformationModels.GEOLOGIC, MagAreaRelationships.ELL_B,
////					AveSlipForRupModels.ELLSWORTH_B, SlipAlongRuptureModels.TAPERED, filter);
//				SectionClusterList clusters = new SectionClusterList(faultModel, defModel, filter,
//						faultSectionData, subSectionDistances, subSectionAzimuths);
//				FaultSystemRupSet rupSet = new FakeRupSet(clusters);
//				numChecked++;
//				System.out.println(numChecked+"/"+filters.size()+": size="+rupSet.getNumRuptures());
//				if (rupSet.getNumRuptures() < minRups) {
//					// then it's worth checking
//					if (validateRupSet(rupSet)) {
//						minRups = rupSet.getNumRuptures();
//						best = filter;
//						System.out.println(numChecked+"/"+filters.size()+": Got a new best ("+minRups+"): "+best);
//					}
//				}
//			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}
	
	private static class FakeRupSet extends InversionFaultSystemRupSet {
		
		private List<FaultSectionPrefData> data;
		private List<List<Integer>> rups;
		private int numRups;
		private FaultModels fm;
		
		public FakeRupSet(SectionClusterList clusters) {
			super(LogicTreeBranch.fromValues(false, clusters.getFaultModel()), clusters, clusters.getFaultSectionData());
			this.data = clusters.getFaultSectionData();
			rups = new ArrayList<List<Integer>>();
			for (SectionCluster c : clusters) {
				rups.addAll(c.getSectionIndicesForRuptures());
			}
			numRups = rups.size();
			fm = clusters.getFaultModel();
		}

		@Override
		public int getNumRuptures() {
			return numRups;
		}
		
		@Override
		public List<Integer> getSectionsIndicesForRup(int rupIndex) {
			return rups.get(rupIndex);
		}
		
		@Override
		public FaultSectionPrefData getFaultSectionData(int sectIndex) {
			return data.get(sectIndex);
		}
		
		@Override
		public int getNumSections() {
			return data.size();
		}
		
		@Override
		public FaultModels getFaultModel() {
			return fm;
		}

		@Override
		public List<List<Integer>> getSectionIndicesForAllRups() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double[] getMagForAllRups() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getMagForRup(int rupIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double[] getAveSlipForAllRups() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getAveSlipForRup(int rupIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public SlipAlongRuptureModels getSlipAlongRuptureModel() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double[] getAveRakeForAllRups() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getAveRakeForRup(int rupIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double[] getAreaForAllRups() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getAreaForRup(int rupIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double[] getAreaForAllSections() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getAreaForSection(int sectIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public List<FaultSectionPrefData> getFaultSectionDataList() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public List<FaultSectionPrefData> getFaultSectionDataForRupture(
				int rupIndex) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getSlipRateForSection(int sectIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double[] getSlipRateForAllSections() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double getSlipRateStdDevForSection(int sectIndex) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public double[] getSlipRateStdDevForAllSections() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public String getInfoString() {
			// TODO Auto-generated method stub
			return null;
		}
		
		public void setInfoString(String info) {
			
		}

		@Override
		public List<Integer> getCloseSectionsList(int sectIndex) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public List<List<Integer>> getCloseSectionsListList() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public int getNumClusters() {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public int getNumRupturesForCluster(int index) {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public List<Integer> getSectionsForCluster(int index) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public List<Integer> getRupturesForCluster(int index)
				throws IndexOutOfBoundsException {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public DeformationModels getDeformationModel() {
			// TODO Auto-generated method stub
			return null;
		}
		
		// net added the following to parent, but doesn't know what to do here
		@Override
		public double getLengthForRup(int rupIndex) {
			throw new RuntimeException("not yet implemented");
		}

		
	}
	
	

}
