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

public class CoulombImperialMappingFix {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		double maxDistance = 5d;
		
		int tolayParent = 761;
		
		FaultModels[] fms = { FaultModels.FM3_1, FaultModels.FM3_2 };
		
		for (FaultModels fm : fms) {
			DeformationModelFetcher fetch = new DeformationModelFetcher(
					fm, DeformationModels.GEOLOGIC, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1);
			List<FaultSectionPrefData> origSubSects = fetch.getSubSectionList();
			
			List<FaultSectionPrefData> newSubSects = Lists.newArrayList();
			// now remap sections
			int origTolayIndex = -1;
			int origTolayNum = 0;
			int newSectIndex = 0;
			for (int i=0; i<origSubSects.size(); i++) {
				FaultSectionPrefData sect = origSubSects.get(i);
				if (sect.getParentSectionId() == tolayParent) {
					if (origSubSects.get(i-1).getParentSectionId() != 97) {
						// this means first tolay
						origTolayIndex = newSubSects.size();
					}
					origTolayNum++;
				} else {
					FaultSectionPrefData newSect = new FaultSectionPrefData();
					newSect.setFaultSectionPrefData(origSubSects.get(i));
					newSect.setSectionId(newSectIndex++);
					newSubSects.add(newSect);
				}
			}
			Preconditions.checkState(newSectIndex == newSubSects.size());
			Preconditions.checkState(newSubSects.size() == origSubSects.size()-origTolayNum);
			
			Map<IDPairing, Double> newDistsMap = DeformationModelFetcher.calculateDistances(maxDistance, newSubSects);
			for (IDPairing pairing : Lists.newArrayList(newDistsMap.keySet()))
				newDistsMap.put(pairing.getReversed(), newDistsMap.get(pairing));
			
			CoulombRates origRates = CoulombRates.loadUCERF3CoulombRates(fm);
			Map<IDPairing, CoulombRatesRecord> newCoulombRatesMap = Maps.newHashMap();
			
			for (IDPairing pairing : origRates.keySet()) {
				// skip if either is Tolay
				if (isTolay(pairing.getID1(), origTolayIndex, origTolayNum)
						|| isTolay(pairing.getID2(), origTolayIndex, origTolayNum))
					continue;
				int id1 = getRemappedID(pairing.getID1(), origTolayIndex, origTolayNum);
				int id2 = getRemappedID(pairing.getID2(), origTolayIndex, origTolayNum);
				
				newCoulombRatesMap.put(new IDPairing(id1, id2), origRates.get(pairing));
			}
			
			// now fix records where Tolay used to be an option
			for (int id1=0; id1<newSubSects.size(); id1++) {
				int origID1;
				if (id1 < origTolayIndex)
					origID1 = id1;
				else
					origID1 = id1 + origTolayNum;
				List<IDPairing> pairings = Lists.newArrayList();
				boolean recalc = false;
				for (IDPairing pair : origRates.keySet()) {
					if (pair.getID1() != origID1)
						continue;
					// check if tolay
					if (isTolay(pair.getID2(), origTolayIndex, origTolayNum))
						recalc = true;
					else
						pairings.add(pair);
				}
				// do the calcs anyway just to check
				double totDS = 0;
				double totDCFF = 0;
				for (IDPairing pair : pairings) {
					totDS += origRates.get(pair).getShearStressChange();
					totDCFF += origRates.get(pair).getCoulombStressChange();
				}
				for (IDPairing pair : pairings) {
					CoulombRatesRecord  rec = origRates.get(pair);
					double pds = rec.getShearStressChange() / totDS;
					double pdcff = rec.getCoulombStressChange() / totDCFF;
					if (recalc) {
						IDPairing remapped = new IDPairing(getRemappedID(pair.getID1(), origTolayIndex, origTolayNum),
								getRemappedID(pair.getID2(), origTolayIndex, origTolayNum));
						System.out.println("Fixing coulomb rates for orig pair="+pair+", new pair="+remapped);
						newCoulombRatesMap.remove(remapped);
						rec = new CoulombRatesRecord(remapped, rec.getShearStressChange(), pds, rec.getCoulombStressChange(), pdcff);
						newCoulombRatesMap.put(remapped, rec);
					} else {
						// lets just double check them
						Preconditions.checkState(Math.abs(pds - rec.getShearStressProbability()) < 0.0001,
								pds+" != "+rec.getShearStressProbability()+" ("+pair+")");
						Preconditions.checkState(Math.abs(pdcff - rec.getCoulombStressProbability()) < 0.0001,
								pdcff+" != "+rec.getCoulombStressProbability()+" ("+pair+")");
					}
				}
			}
			
			System.out.println("Orig coulomb count: "+origRates.size());
			System.out.println("New coulomb count: "+newCoulombRatesMap.size());
			
			// check for completeness
//			Preconditions.checkState(newCoulombRatesMap.size() == origRates.size()+(newCnt*2));
//			for (IDPairing pair : newDistsMap.keySet())
//				Preconditions.checkNotNull(newCoulombRatesMap.get(pair), "New coulomb missing mapping: "+pair);
			for (IDPairing pair : newCoulombRatesMap.keySet())
				Preconditions.checkNotNull(newDistsMap.get(pair),
						"New coulomb has extra mapping: "+pair+" (orig has it? "+origRates.containsKey(pair)+")");
			
			CoulombRates.writeExcelFile(new CoulombRates(newCoulombRatesMap), newDistsMap,
					new File("/tmp/2013_04_08-Stress_Table-"+fm.getShortName().replaceAll("_", ".")+".xls"));
		}
	}
	private static boolean isTolay(int id, int origTolayIndex, int origTolayNum) {
		return id >= origTolayIndex && id < origTolayIndex+origTolayNum;
	}
	
	private static int getRemappedID(int id, int origTolayIndex, int origTolayNum) {
		if (id < origTolayIndex)
			// before imperial, simple case
			return id;
		return id - origTolayNum;
	}

}
