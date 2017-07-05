package scratch.kevin.ucerf3.inversion;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;

public class NamedFaultsTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		for (FaultModels fm : FaultModels.values()) {
			if (fm != FaultModels.FM2_1)
				continue;
			Map<Integer, List<Integer>> named = fm.getNamedFaultsMap();
			
			ArrayList<FaultSectionPrefData> sects = fm.fetchFaultSections();
			
			for (FaultSectionPrefData fault : sects) {
				int id = fault.getSectionId();
				if (id == 402)
					System.out.println("402: "+fault.getSectionName());
				if (!named.containsKey(id))
					System.out.println(fm+" missing named fault ID: "+id);
			}
			FaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.cachedForBranch(
					true, fm, DeformationModels.forFaultModel(fm).get(0));
			for (FaultSectionPrefData subSect : rupSet.getFaultSectionDataList()) {
				int parentID = subSect.getParentSectionId();
				if (!named.containsKey(parentID)) {
					System.out.println(fm+" missing named faults for parent: "+parentID);
					System.out.println("sub sect: "+subSect.getName()+" ("+subSect.getSectionId()+")");
				}
			}
		}
	}

}
