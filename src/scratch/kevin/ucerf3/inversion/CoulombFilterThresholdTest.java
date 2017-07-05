package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesRecord;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester.TestType;
import scratch.UCERF3.inversion.laughTest.CoulombFilter;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.IDPairing;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class CoulombFilterThresholdTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_1;
		CoulombRates rates = CoulombRates.loadUCERF3CoulombRates(fm);
		
		File file = new File("/tmp/F2F_section_connection_review.csv");
		List<BiasiCoulombRecord> biasi = loadBiasiTable(file);
		
		File outputFile = new File("/tmp/coulomb_thresholds.csv");
		
		double[] minStresses = { 0.5, 0.75, 1, 1.25, 1.5, Double.POSITIVE_INFINITY };
		double[] pdcffs = { 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 1 };
		
		LaughTestFilter filter = LaughTestFilter.getDefault();
		filter.setAllowSingleSectDuringJumps(true);
		filter.setCoulombFilter(new CoulombRatesTester(TestType.COULOMB_STRESS, 0d, 0d, 0d, true, true));
		
		FaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(
				filter, InversionFaultSystemRupSetFactory.DEFAULT_ASEIS_VALUE, fm);
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine("PDCFF", "Min Stress For Exclusion", "Total Rups", "Rups Excluded",
				"Biasi 'Likely' Included", "Biasi 'Likely' Excluded", "Biasi 'Possible' Included", "Biasi 'Possible' Excluded",
				"Biasi 'No' Included", "Biasi 'No' Excluded", "Likely's Excluded", "Possibles Excluded", "No's Included");
		
		for (double minStress : minStresses) {
			for (double pdcff : pdcffs) {
				CoulombRatesTester tester = new CoulombRatesTester(TestType.COULOMB_STRESS, pdcff, pdcff, minStress, true, true);
				
				CoulombFilter cf = new CoulombFilter(rates, tester);
				
				int filteredCount = 0;
				for (int r=0; r<rupSet.getNumRuptures(); r++) {
					List<FaultSectionPrefData> rupture = rupSet.getFaultSectionDataForRupture(r);
					if (cf.doesRupturePass(rupture))
						filteredCount++;
				}
				
				Map<Likelihood, List<BiasiCoulombRecord>> includedMap = Maps.newHashMap();
				Map<Likelihood, List<BiasiCoulombRecord>> excludedMap = Maps.newHashMap();
				
				for (BiasiCoulombRecord rec : biasi) {
					List<CoulombRatesRecord> forwardRates = Lists.newArrayList(rates.get(new IDPairing(rec.id1, rec.id2)));
					List<CoulombRatesRecord> backwardRates = Lists.newArrayList(rates.get(new IDPairing(rec.id2, rec.id1)));
					
					Likelihood likelihood = rec.likelihood;
					
					boolean included = tester.doesRupturePass(forwardRates, backwardRates);
					
					Map<Likelihood, List<BiasiCoulombRecord>> map;
					if (included)
						map = includedMap;
					else
						map = excludedMap;
					
					List<BiasiCoulombRecord> list = map.get(likelihood);
					if (list == null) {
						list = Lists.newArrayList();
						map.put(likelihood, list);
					}
					list.add(rec);
				}
				
				List<String> line = Lists.newArrayList();
				line.add((float)pdcff+"");
				line.add((float)minStress+"");
				line.add(filteredCount+"");
				line.add((rupSet.getNumRuptures() - filteredCount)+"");
				line.add(sizeNullZero(includedMap.get(Likelihood.LIKELY))+"");
				line.add(sizeNullZero(excludedMap.get(Likelihood.LIKELY))+"");
				line.add(sizeNullZero(includedMap.get(Likelihood.POSSIBLE))+"");
				line.add(sizeNullZero(excludedMap.get(Likelihood.POSSIBLE))+"");
				line.add(sizeNullZero(includedMap.get(Likelihood.NO))+"");
				line.add(sizeNullZero(excludedMap.get(Likelihood.NO))+"");
				List<FaultSectionPrefData> datas = rupSet.getFaultSectionDataList();
				line.add(getPairingsString(excludedMap.get(Likelihood.LIKELY), datas, rates));
				line.add(getPairingsString(excludedMap.get(Likelihood.POSSIBLE), datas, rates));
				line.add(getPairingsString(includedMap.get(Likelihood.NO), datas, rates));
				
				csv.addLine(line);
			}
		}
		
		csv.writeToFile(outputFile);
	}
	
	private static int sizeNullZero(List<?> list) {
		if (list == null)
			return 0;
		return list.size();
	}
	
	private static String getPairingsString(List<BiasiCoulombRecord> recs, List<FaultSectionPrefData> datas, CoulombRates rates) {
		if (recs == null || recs.isEmpty())
			return "";
		List<String> list = Lists.newArrayList();
		for (BiasiCoulombRecord rec : recs) {
			CoulombRatesRecord fwd = rates.get(new IDPairing(rec.id1, rec.id2));
			CoulombRatesRecord bkw = rates.get(new IDPairing(rec.id2, rec.id1));
			list.add(datas.get(rec.id1).getSectionName()+" ["+(float)fwd.getCoulombStressProbability()
					+"/"+fwd.getCoulombStressChange()+"] <==> ["+(float)bkw.getCoulombStressProbability()+"/"
					+bkw.getCoulombStressChange()+"] "+datas.get(rec.id2).getSectionName());
		}
		return Joiner.on(",").join(list);
	}
	
	private static List<BiasiCoulombRecord> loadBiasiTable(File file) throws IOException {
		List<BiasiCoulombRecord> recs = Lists.newArrayList();
		
		CSVFile<String> csv = CSVFile.readFile(file, false);
		
		for (List<String> line : csv) {
			if (line.isEmpty())
				continue;
			String classification = line.get(0).trim();
			if (classification.isEmpty())
				continue;
			classification = classification.toLowerCase();
			Likelihood likelihood;
			if (classification.startsWith("no"))
				likelihood = Likelihood.NO;
			else if (classification.startsWith("possible"))
				likelihood = Likelihood.POSSIBLE;
			else if (classification.startsWith("likely"))
				likelihood = Likelihood.LIKELY;
			else
				continue;
			
			int id1 = Integer.parseInt(line.get(1));
			int id2 = Integer.parseInt(line.get(2));
			
			recs.add(new BiasiCoulombRecord(likelihood, id1, id2));
		}
		
		
		return recs;
	}
	
	private static enum Likelihood {
		LIKELY,
		POSSIBLE,
		NO;
	}
	
	private static class BiasiCoulombRecord {
		private Likelihood likelihood;
		private int id1;
		private int id2;
		public BiasiCoulombRecord(Likelihood likelihood, int id1, int id2) {
			super();
			this.likelihood = likelihood;
			this.id1 = id1;
			this.id2 = id2;
		}
	}

}
