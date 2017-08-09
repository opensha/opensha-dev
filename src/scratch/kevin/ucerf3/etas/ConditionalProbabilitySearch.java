package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ConditionalProbabilitySearch {

	public static void main(String[] args) {
		File resultsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin");
		
		double[] durations = { 1d/365.25, 7d/365.25, 1d/12d, 1d, 10d };
		String[] labels = { "1 day", "1 week", "1 month", "1 year", "10 years" };
		
		int[] numInitials = new int[durations.length];
		int[] counts = new int[durations.length];
		int[] triggeredCounts = new int[durations.length];
		
		double initialMag = 7d;
		double targetMag = 7d;
		
		for (List<ETAS_EqkRupture> catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(resultsFile, 0d)) {
			Map<Integer, ETAS_EqkRupture> catalogMap = new HashMap<>();
			for (ETAS_EqkRupture rup : catalog)
				catalogMap.put(rup.getID(), rup);
			long maxTime = catalog.get(catalog.size()-1).getOriginTime();
			
			for (int i=0; i<catalog.size(); i++) {
				ETAS_EqkRupture initial = catalog.get(i);
				if (initial.getMag() >= initialMag) {
					long ot = initial.getOriginTime();
					
					for (int d=0; d<durations.length; d++) {
						long maxOT = (long)(ot + durations[d]*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
						if (maxOT > maxTime)
							// window cut off by end of catalog
							continue;
						numInitials[d]++;
						for (int j=i+1; j<catalog.size(); j++) {
							ETAS_EqkRupture target = catalog.get(j);
							if (target.getOriginTime() > maxOT)
								break;
							if (target.getMag() >= targetMag) {
								// it's a match
								counts[d]++;
								break; // so we don't double count
							}
						}
						// actually triggered
						for (int j=i+1; j<catalog.size(); j++) {
							ETAS_EqkRupture target = catalog.get(j);
							if (target.getOriginTime() > maxOT)
								break;
							if (target.getMag() >= targetMag) {
								// it's a match
								// is it triggered?
								if (descendsFrom(catalogMap, initial, target)) {
									triggeredCounts[d]++;
									break; // so we don't double count
								}
							}
						}
					}
				}
			}
		}
		
		System.out.println("Full Catalog Probabilities");
		for (int d=0; d<durations.length; d++) {
			double percent = 100d*(double)counts[d]/(double)numInitials[d];
			System.out.println(labels[d]+":\t"+(float)percent+" %");
		}
		System.out.println();
		System.out.println("Triggered Catalog Probabilities");
		for (int d=0; d<durations.length; d++) {
			double percent = 100d*(double)triggeredCounts[d]/(double)numInitials[d];
			System.out.println(labels[d]+":\t"+(float)percent+" %");
		}
	}
	
	private static boolean descendsFrom(Map<Integer, ETAS_EqkRupture> catalog, ETAS_EqkRupture initial, ETAS_EqkRupture target) {
		ETAS_EqkRupture curRup = target;
		while (true) {
			int parentID = curRup.getParentID();
			
			if (parentID == initial.getID())
				return true;
			if (parentID < 0 || !catalog.containsKey(parentID))
				return false;
			
			curRup = catalog.get(parentID);
		}
	}

}
