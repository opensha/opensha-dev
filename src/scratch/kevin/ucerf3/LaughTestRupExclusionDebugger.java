package scratch.kevin.ucerf3;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.ClassUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.laughTest.AbstractLaughTest;
import scratch.UCERF3.inversion.laughTest.AzimuthChangeFilter;
import scratch.UCERF3.inversion.laughTest.BuggyCoulombFilter;
import scratch.UCERF3.inversion.laughTest.CoulombFilter;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class LaughTestRupExclusionDebugger {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		int[] sects = {612, 613, 1495, 1496};
//		int[] sects = {2343, 2342, 1962, 1925, 1926};
//		int[] sects = {1137, 1136, 503, 502};
//		List<Integer> sectsList = Lists.newArrayList();
//		for (int i=1502; i>=1498; i--)
//			sectsList.add(i);
//		for (int i=613; i>=594; i--)
//			sectsList.add(i);
//		for (int i=635; i>=622; i--)
//			sectsList.add(i);
//		for (int i=1832; i<=1836; i++)
//			sectsList.add(i);
//		int[] sects = Ints.toArray(sectsList);
//		int[] sects = {1399, 1398, 1949, 1950};
		int[] sects = {2593, 2592, 334, 335};
		
		boolean applyGarlockPintoMtnFix = true;
		
		LaughTestFilter filter = LaughTestFilter.getDefault();
		filter.setAllowSingleSectDuringJumps(true);
		
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		DeformationModelFetcher fetch = new DeformationModelFetcher(fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1);
		List<FaultSectionPrefData> datas = fetch.getSubSectionList();
		
		Map<IDPairing, Double> subSectionDistances = fetch.getSubSectionDistanceMap(filter.getMaxJumpDist());
		Map<IDPairing, Double> subSectionAzimuths = fetch.getSubSectionAzimuthMap(subSectionDistances.keySet());
		
		List<FaultSectionPrefData> rupture = Lists.newArrayList();
		for (int sect : sects)
			rupture.add(datas.get(sect));
		
		CoulombRates coulombRates = CoulombRates.loadUCERF3CoulombRates(fm);
		
		List<List<Integer>> sectionConnectionsListList = SectionClusterList.computeCloseSubSectionsListList(
				datas, subSectionDistances, filter.getMaxJumpDist(), coulombRates);
		
		List<AbstractLaughTest> laughTests = filter.buildLaughTests(subSectionAzimuths, subSectionDistances, null, coulombRates,
				applyGarlockPintoMtnFix, sectionConnectionsListList, datas);
		
		for (int i=1; i<rupture.size(); i++) {
			if (rupture.get(i).getParentSectionId() != rupture.get(i-1).getParentSectionId()) {
				IDPairing pairing = new IDPairing(rupture.get(i-1).getSectionId(), rupture.get(i).getSectionId());
				if (!sectionConnectionsListList.get(pairing.getID1()).contains(pairing.getID2())) {
					System.out.println("Pairing doesn't exist in connections list: "+pairing);
					System.out.println("\tPossibilities for "+pairing.getID1()+": "
							+Joiner.on(",").join(sectionConnectionsListList.get(pairing.getID1())));
					System.out.println("\tPossibilities for "+pairing.getID2()+": "
							+Joiner.on(",").join(sectionConnectionsListList.get(pairing.getID2())));
				}
			}
		}
		
		boolean failedCoulomb = false;
		for (AbstractLaughTest test : laughTests) {
			if (!test.doesRupturePass(rupture)) {
				System.out.println("FAILED: "+ClassUtils.getClassNameWithoutPackage(test.getClass()));
				if (test instanceof CoulombFilter || test instanceof BuggyCoulombFilter)
					failedCoulomb = true;
			}
		}
		
		if (failedCoulomb) {
			System.out.println("Debugging coulomb");
			// print out coulomb at junctions
			int prevID = -1;
			int prevParent = -1;
			for (FaultSectionPrefData sect : rupture) {
				int id = sect.getSectionId();
				int parent = sect.getParentSectionId();
				if (prevParent != -1 && prevParent != parent) {
					// junction
					IDPairing pair = new IDPairing(id, prevID);
					System.out.println("Junction: "+pair.getReversed());
					System.out.println("\t"+coulombRates.get(pair.getReversed()));
					System.out.println("\t"+coulombRates.get(pair));
				}
				prevParent = parent;
				prevID = id;
			}
		}
	}

}
