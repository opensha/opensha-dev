package scratch.kevin.cybershake.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class BiasiDownsampleMapper {
	
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/biasi_downsample_tests");
		FaultSystemSolution origSol = FaultSystemIO.loadSol(new File(dir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip"));
		
		FaultSystemRupSet downsampledRupSet = FaultSystemIO.loadRupSet(new File(dir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledBoth_mean_sol.zip"));
		File outputFile = new File(dir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledBothMapped_mean_sol.zip");
//		FaultSystemRupSet downsampledRupSet = FaultSystemIO.loadRupSet(new File(dir,
//				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledStarts_mean_sol.zip"));
//		File outputFile = new File(dir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledStartsMapped_mean_sol.zip");
//		FaultSystemRupSet downsampledRupSet = FaultSystemIO.loadRupSet(new File(dir,
//				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledEnds_mean_sol.zip"));
//		File outputFile = new File(dir, "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledEndsMapped_mean_sol.zip");
		
		FaultSystemRupSet origRupSet = origSol.getRupSet();
		
		List<HashSet<Integer>> downsampledRupSets = Lists.newArrayList();
		for (int r=0; r<downsampledRupSet.getNumRuptures(); r++)
			downsampledRupSets.add(new HashSet<Integer>(downsampledRupSet.getSectionsIndicesForRup(r)));
		
		HistogramFunction hist = new HistogramFunction(0d, 21, 1d);
		
		double[] downsampledRates = new double[downsampledRupSet.getNumRuptures()];
		
		int numNoParentsMatch = 0;
		int numNoMatch = 0;
		
		for (int rupIndex=0; rupIndex<origRupSet.getNumRuptures(); rupIndex++) {
			if (rupIndex % 5000 == 0)
				System.out.println("Rupture "+rupIndex);
			List<Integer> sects = origRupSet.getSectionsIndicesForRup(rupIndex);
			double mag = origRupSet.getMagForRup(rupIndex);
			HashSet<Integer> sectsSet = new HashSet<Integer>(sects);
			
			List<Integer> parents = origRupSet.getParentSectionsForRup(rupIndex);
			HashSet<Integer> possibles = new HashSet<Integer>(downsampledRupSet.getRupturesForParentSection(parents.get(0)));
			for (int i=1; i<parents.size(); i++)
				possibles.retainAll(downsampledRupSet.getRupturesForParentSection(parents.get(i)));
			if (possibles.isEmpty()) {
//				System.err.println("No ruptures exist with these parents:");
//				String prevParent = null;
//				for (FaultSectionPrefData sect : origRupSet.getFaultSectionDataForRupture(rupIndex)) {
//					String parentName = sect.getParentSectionName();
//					if (prevParent == null || !prevParent.equals(parentName))
//						System.err.println("\t"+parentName);
//					prevParent = parentName;
//				}
//				System.exit(0); 
				numNoParentsMatch++;
//				possibles = allDownsampledRups;
//				// consider all ruptures which touch any
//				possibles = new HashSet<Integer>(downsampledRupSet.getRupturesForParentSection(parents.get(0)));
//				for (int i=1; i<parents.size(); i++)
//					possibles.addAll(downsampledRupSet.getRupturesForParentSection(parents.get(i)));
				// consider all ruptures which touch the middle section
				possibles = new HashSet<Integer>(downsampledRupSet.getRupturesForSection(sects.get(sects.size()/2)));
			}
			Preconditions.checkState(!possibles.isEmpty(), "");
			double minMag = mag - 0.2;
			double maxMag = mag + 0.2;
			
			List<Candidate> candidates = Lists.newArrayList();
			
			double minDiff = Double.MAX_VALUE;
			Candidate best = null;
			
//			for (int r : downsampledRupSet.getRupturesForSection(middleSect)) {
			for (int r : possibles) {
				double oMag = downsampledRupSet.getMagForRup(r);
				if (oMag < minMag || oMag > maxMag)
					continue;
				Candidate candidate = new Candidate(r, oMag, mag);
				HashSet<Integer> oSects = downsampledRupSets.get(r);
				for (Integer s : sectsSet)
					if (!oSects.contains(s))
						candidate.sectsMissing++;
				for (Integer s : oSects)
					if (!sectsSet.contains(s))
						candidate.sectsExtra++;
				candidates.add(candidate);
				
				double diff = candidate.diff();
				if (diff < minDiff) {
					minDiff = diff;
					best = candidate;
					if (diff == 0)
						// exact match
						break;
				}
			}
			if (best == null) {
				numNoMatch++;
				System.out.println("WARNING: NO MATCH FOUND WITHIN MAG TOLERANCE!!!!!!");
				continue;
			}
			Preconditions.checkNotNull(best);
			if (best.diff() > 20) {
				System.out.println("Terrible match at rupIndex="+rupIndex);
				System.out.print("Orig:\t");
				printRup(sects);
				System.out.print("Match:\t");
				printRup(downsampledRupSet.getSectionsIndicesForRup(best.rupIdnex));
			}
//			System.out.println("Best match diff="+minDiff+", origMag="+mag+", mappedMag="+best.mag);
			downsampledRates[best.rupIdnex] += origSol.getRateForRup(rupIndex);
			hist.add(hist.getClosestXIndex(best.diff()), 1d);
		}
		System.out.println(hist);
		System.out.println(numNoMatch+"/"+origRupSet.getNumRuptures()+" have NO MATCH");
		System.out.println(numNoParentsMatch+"/"+origRupSet.getNumRuptures()+" have no parent section match");
		FaultSystemIO.writeSol(new FaultSystemSolution(downsampledRupSet, downsampledRates), outputFile);
	}
	
	private static void printRup(List<Integer> rup) {
//		rup = Lists.newArrayList(rup);
//		Collections.sort(rup);
		System.out.println(Joiner.on(",").join(rup));
	}
	
	private static class Candidate {
		int rupIdnex;
		double mag;
		double magDelta;
		int sectsMissing;
		int sectsExtra;
		
		public Candidate(int rupIndex, double mag, double origMag) {
			this.rupIdnex = rupIndex;
			this.mag = mag;
			magDelta = Math.abs(origMag - mag);
		}
		
		private double diff() {
			return sectsMissing + sectsExtra + magDelta;
		}
	}

}
