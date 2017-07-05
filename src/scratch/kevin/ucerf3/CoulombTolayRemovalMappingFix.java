package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesRecord;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class CoulombTolayRemovalMappingFix {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		int origNum = 9;
		int newNum = 12;
		int newCnt = newNum-origNum;
		
		double maxDistance = 5d;
		
		int imperialParent = 97;
		List<FaultSectionPrefData> newImperialSubSects = null;
		
		FaultModels[] fms = { FaultModels.FM3_1, FaultModels.FM3_2 };
		
		for (FaultModels fm : fms) {
			if (newImperialSubSects == null) {
				FaultSectionPrefData section = fm.fetchFaultSectionsMap().get(imperialParent);
				double ddw = section.getOrigDownDipWidth();
				newImperialSubSects = section.getSubSectionsList(ddw*0.5, 0, 2);
			}
			DeformationModelFetcher fetch = new DeformationModelFetcher(
					fm, DeformationModels.GEOLOGIC, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1);
			List<FaultSectionPrefData> origSubSects = fetch.getSubSectionList();
			
			List<FaultSectionPrefData> newSubSects = Lists.newArrayList();
			// now remap sections
			int firstImperialIndex = -1;
			int firstPostImperialIndex = -1;
			int newSectIndex = 0;
			for (int i=0; i<origSubSects.size(); i++) {
				FaultSectionPrefData sect = origSubSects.get(i);
				if (sect.getParentSectionId() == 97) {
					if (origSubSects.get(i-1).getParentSectionId() != 97) {
						// this means first imperial
						firstImperialIndex = newSubSects.size();
						for (int j=0; j<newImperialSubSects.size(); j++) {
							FaultSectionPrefData newImperialSect = new FaultSectionPrefData();
							newImperialSect.setFaultSectionPrefData(newImperialSubSects.get(j));
							newImperialSect.setSectionId(newSectIndex++);
							newSubSects.add(newImperialSect);
						}
						firstPostImperialIndex = newSubSects.size();
					}
				} else {
					FaultSectionPrefData newSect = new FaultSectionPrefData();
					newSect.setFaultSectionPrefData(origSubSects.get(i));
					newSect.setSectionId(newSectIndex++);
					newSubSects.add(newSect);
				}
			}
			Preconditions.checkState(firstPostImperialIndex-newCnt-origNum == firstImperialIndex,
					"Somethings not right: "+firstPostImperialIndex+"-"+newCnt+"-"+origNum+" != "+firstImperialIndex);
			Preconditions.checkState(newSectIndex == newSubSects.size());
			Preconditions.checkState(newSubSects.size() == origSubSects.size()+newCnt);
			
			Map<IDPairing, Double> newDistsMap = DeformationModelFetcher.calculateDistances(maxDistance, newSubSects);
			for (IDPairing pairing : Lists.newArrayList(newDistsMap.keySet()))
				newDistsMap.put(pairing.getReversed(), newDistsMap.get(pairing));
			
			CoulombRates origRates = CoulombRates.loadUCERF3CoulombRates(fm);
			Map<IDPairing, CoulombRatesRecord> newCoulombRatesMap = Maps.newHashMap();
			
			for (IDPairing pairing : origRates.keySet()) {
				int id1 = getRemappedID(pairing.getID1(), firstImperialIndex, origNum, newNum);
				int id2 = getRemappedID(pairing.getID2(), firstImperialIndex, origNum, newNum);
				
				newCoulombRatesMap.put(new IDPairing(id1, id2), origRates.get(pairing));
			}
			// now we have to add the fake records for the new values
			for (int i=firstImperialIndex; i<firstImperialIndex+newCnt; i++) {
				IDPairing pair1 = new IDPairing(i, i+1);
				IDPairing pair2 = new IDPairing(i+1, i);
				
				CoulombRatesRecord rec1, rec2;
				if (i == firstImperialIndex)
					rec1 = new CoulombRatesRecord(pair1, 10d, 1d, 10d, 1d);
				else
					rec1 = new CoulombRatesRecord(pair1, 10d, 0.5, 10d, 0.5);
				
				if (i == firstImperialIndex+(newCnt-1)) {
					// last one, fix things up
					IDPairing nextPairing = new IDPairing(i+1, i+2);
					CoulombRatesRecord nextFWD = newCoulombRatesMap.get(nextPairing);
					double nextDS = nextFWD.getShearStressChange();
					double nextDCFF = nextFWD.getCoulombStressChange();
					
					double nextTotDS = nextDS + 10;
					double nextTotDCFF = nextDCFF + 10;
					
					rec2 = new CoulombRatesRecord(pair2, 10d, 10d/nextTotDS, 10d, 10d/nextTotDCFF);
					// now fix the old first section
					newCoulombRatesMap.remove(nextPairing);
					newCoulombRatesMap.put(nextPairing, new CoulombRatesRecord(
							nextPairing, nextDS, nextDS/nextTotDS, nextDCFF, nextDCFF/nextTotDCFF));
				} else {
					rec2 = new CoulombRatesRecord(pair2, 10d, 0.5, 10d, 0.5);
				}
				newCoulombRatesMap.put(pair1, rec1);
				newCoulombRatesMap.put(pair2, rec2);
			}
			
			// check for completeness
			Preconditions.checkState(newCoulombRatesMap.size() == origRates.size()+(newCnt*2));
//			for (IDPairing pair : newDistsMap.keySet())
//				Preconditions.checkNotNull(newCoulombRatesMap.get(pair), "New coulomb missing mapping: "+pair);
			for (IDPairing pair : newCoulombRatesMap.keySet())
				Preconditions.checkNotNull(newDistsMap.get(pair),
						"New coulomb has extra mapping: "+pair+" (orig has it? "+origRates.containsKey(pair)+")");
			
			CoulombRates.writeExcelFile(new CoulombRates(newCoulombRatesMap), newDistsMap,
					new File("/tmp/2013_04_08-Stress_Table-"+fm.getShortName().replaceAll("_", ".")+".xls"));
		}
	}
	
	private static int getRemappedID(int id, int minImperialID, int origNumImperial, int newNumImperial) {
		if (id < minImperialID)
			// before imperial, simple case
			return id;
		// just add the difference. this works for Imperial IDs as well as we want to add the nubsections at the south
		// end which is the first index.
		return id + (newNumImperial-origNumImperial);
	}

}
