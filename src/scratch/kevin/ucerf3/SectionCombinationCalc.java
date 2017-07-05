package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.utils.FaultSystemIO;

public class SectionCombinationCalc {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		Table<Integer, Integer, String> sectMap = HashBasedTable.create();
		sectMap.put(651, 0, "RC"); // Rodgers Creek - Healdsburg 2011 CFM
		sectMap.put(639, 1, "HN"); // Hayward (No) 2011 CFM
		sectMap.put(638, 2, "HS"); // Hayward (So) 2011 CFM
		sectMap.put(637, 3, "HSE"); // Hayward (So) extension 2011 CFM
		sectMap.put(601, 4, "CN"); // Calaveras (No) 2011 CFM
		sectMap.put(602, 5, "CC"); // Calaveras (Central) 2011 CFM
		sectMap.put(603, 6, "CS"); // Calaveras (So) 2011 CFM
		sectMap.put(621, 7, "CSE"); // Calaveras (So) - Paicines extension 2011 CFM"
		
		HashSet<Integer> potentialRups = new HashSet<Integer>();
		for (int sect : sectMap.rowKeySet())
			potentialRups.addAll(rupSet.getRupturesForParentSection(sect));
		Map<String, Double> maxMagForCombos = Maps.newHashMap();
		for (int rup : potentialRups) {
			boolean allOnSects = true;
			HashSet<Integer> mySects = new HashSet<Integer>();
			for (FaultSectionPrefData sect : rupSet.getFaultSectionDataForRupture(rup)) {
				if (!sectMap.containsRow(sect.getParentSectionId())) {
					allOnSects = false;
					break;
				}
				mySects.add(sect.getParentSectionId());
			}
			if (allOnSects) {
				String key = getKey(mySects, sectMap);
				double myMag = rupSet.getMagForRup(rup);
				Double maxMag = maxMagForCombos.get(key);
				if (maxMag == null || myMag > maxMag)
					maxMagForCombos.put(key, myMag);
			}
		}
		List<String> keys = Lists.newArrayList(maxMagForCombos.keySet());
		Collections.sort(keys);
		
		for (String key : keys)
			System.out.println(key+": "+maxMagForCombos.get(key).floatValue());
	}
	
	private static Joiner j = Joiner.on("+");
	private static String getKey(HashSet<Integer> sects, Table<Integer, Integer, String> sectMap) {
		List<String> labels = Lists.newArrayList();
		List<Integer> orders = Lists.newArrayList();
		for (int sect : sects) {
			Map<Integer, String> row = sectMap.row(sect);
			int order = row.keySet().iterator().next();
			String label = row.values().iterator().next();
			labels.add(label);
			orders.add(order);
		}
		List<String> sortedLabeles = ComparablePairing.getSortedData(orders, labels);
		return j.join(sortedLabeles);
	}

}
