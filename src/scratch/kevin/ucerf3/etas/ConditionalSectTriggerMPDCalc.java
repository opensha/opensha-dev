package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;

public class ConditionalSectTriggerMPDCalc {

	public static void main(String[] args) throws IOException, DocumentException {
		File catalogsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined/results_descendents_m5_preserve.bin");
				+ "2016_10_27-2016_bombay_swarm-10yr-full_td-no_ert-combined/results.bin");
		int parentID = 295;
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogsFile);
		
		long maxOT = calcMaxOTforDuration(catalogs, 7/365.25);
		
		EvenlyDiscretizedFunc mpd = calcCumulativeMPD(catalogs, sol, parentID, 6.5, 8.5, 0.1, maxOT);
		System.out.println(mpd);
	}
	
	public static long calcMaxOTforDuration(Iterable<List<ETAS_EqkRupture>> catalogs, double durationYears) {
		long minOT = Long.MAX_VALUE;
		
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			if (catalog.isEmpty())
				continue;
			long myFirst = catalog.get(0).getOriginTime();
			if (myFirst < minOT)
				minOT = myFirst;
		}
		
		return (long)(minOT + durationYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
	}
	
	public static EvenlyDiscretizedFunc calcCumulativeMPD(Iterable<List<ETAS_EqkRupture>> catalogs, FaultSystemSolution sol,
			int parentID, double minMag, double maxMag, double deltaMag, long maxOT) {
		int numMag = (int)((maxMag - minMag)/deltaMag + 0.5) + 1;
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(minMag, numMag, deltaMag);
		double minMagForInclusion = minMag-0.5*deltaMag;
		
		HashSet<Integer> rupIndexes = new HashSet<Integer>(sol.getRupSet().getRupturesForParentSection(parentID));
//		HashSet<Integer> rupIndexes = new HashSet<Integer>();
//		for (int i=0; i<sol.getRupSet().getNumRuptures(); i++)
//			rupIndexes.add(i);
		
		int numCatalogs = 0;
		
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			numCatalogs++;
			
			double maxMatchMag = 0d;
			
			for (ETAS_EqkRupture rup : catalog) {
				if (maxOT > 0l && rup.getOriginTime() > maxOT)
					break;
				if (rupIndexes.contains(rup.getFSSIndex())) {
					// it's a match, on the fault!
					maxMatchMag = Math.max(maxMatchMag, rup.getMag());
				}
			}
			
			for (int i=0; i<func.size(); i++) {
				if (maxMatchMag < func.getX(i))
					break;
				func.add(i, 1d);
			}
			
//			for (ETAS_EqkRupture rup : catalog) {
//				if (maxOT > 0l && rup.getOriginTime() > maxOT)
//					break;
//				if (rupIndexes.contains(rup.getFSSIndex())) {
//					// it's a match, on the fault!
//					if (rup.getMag() < minMagForInclusion)
//						continue;
//					int ind = func.getClosestXIndex(rup.getMag());
//					for (int i=0; i<=ind; i++)
//						func.add(i, 1d);
//				}
//			}
		}
		
		func.scale(1d/(double)numCatalogs);
		
//		for (int i=0; i<func.size(); i++) {
//			double mag = func.getX(i)-0.5*deltaMag;
//			double prob = ETAS_MultiSimAnalysisTools.calcETASParticProb((List<List<ETAS_EqkRupture>>)catalogs,
//					Long.MIN_VALUE, maxOT, rupIndexes, null, mag);
//			func.set(i, prob);
//		}
		
		return func;
	}

}
