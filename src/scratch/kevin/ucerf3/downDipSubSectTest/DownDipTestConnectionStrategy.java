package scratch.kevin.ucerf3.downDipSubSectTest;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.IDPairing;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.UCERF3.inversion.SectionConnectionStrategy;

public class DownDipTestConnectionStrategy implements SectionConnectionStrategy {
	
	private DownDipSubSectBuilder downDipBuilder;
	private int downDipParent;
	private double maxJumpDist;

	public DownDipTestConnectionStrategy(DownDipSubSectBuilder downDipBuilder, double maxJumpDist) {
		this.downDipBuilder = downDipBuilder;
		this.maxJumpDist = maxJumpDist;
		this.downDipParent = downDipBuilder.getParentID();
	}

	@Override
	public List<List<Integer>> computeCloseSubSectionsListList(List<? extends FaultSection> faultSectionData,
			Map<IDPairing, Double> subSectionDistances) {
		List<List<Integer>> connections = new ArrayList<>();
		
		// first add neighbors within each parent section
		for (int s=0; s<faultSectionData.size(); s++) {
			FaultSection sect = faultSectionData.get(s);
			Preconditions.checkState(sect.getSectionId() == s);
			int parent = sect.getParentSectionId();
			List<Integer> myConnections = new ArrayList<>();
			connections.add(myConnections);
			
			if (parent == downDipParent) {
				// special down-dip rules
				myConnections.addAll(downDipBuilder.getNeighbors(sect));
			} else {
				// sequential (within section) rules
				int sBefore = s-1;
				if (s > 0 && faultSectionData.get(sBefore).getParentSectionId() == parent)
					myConnections.add(sBefore);
				int sAfter = s+1;
				if (s < faultSectionData.size()-1 &&
						faultSectionData.get(sAfter).getParentSectionId() == parent)
					myConnections.add(sAfter);
			}
		}
		
		// now add closest connection to each parent (if within max jump dist)
		
		// first bundle by parent section IDs
		Map<Integer, List<FaultSection>> parentSectsMap = new HashMap<>();
		for (FaultSection sect : faultSectionData) {
			List<FaultSection> parentSects = parentSectsMap.get(sect.getParentSectionId());
			if (parentSects == null) {
				parentSects = new ArrayList<>();
				parentSectsMap.put(sect.getParentSectionId(), parentSects);
			}
			parentSects.add(sect);
		}
		// now find and add closest pairings
		for (int p1 : parentSectsMap.keySet()) {
			List<FaultSection> sects1 = parentSectsMap.get(p1);
			for (int p2 : parentSectsMap.keySet()) {
				if (p1 == p2)
					continue;
				List<FaultSection> sects2 = parentSectsMap.get(p2);
				double minDist = Double.POSITIVE_INFINITY;
				IDPairing minPairing = null;
				for (int s1=0; s1<sects1.size(); s1++) {
					for (int s2=0; s2<sects2.size(); s2++) {
						IDPairing pair = new IDPairing(sects1.get(s1).getSectionId(),
								sects2.get(s2).getSectionId());
						Double dist = subSectionDistances.get(pair);
						if (dist != null && dist < minDist) {
							minDist = dist;
							minPairing = pair;
						}
					}
				}
				Preconditions.checkNotNull(minPairing);
				if (minDist <= maxJumpDist)
					connections.get(minPairing.getID1()).add(minPairing.getID2());
			}
		}
		
		return connections;
	}

}
