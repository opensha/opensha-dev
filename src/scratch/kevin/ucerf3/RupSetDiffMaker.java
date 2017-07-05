package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.laughTest.AzimuthChangeFilter;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.FaultSystemIO;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

public class RupSetDiffMaker {

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
//		File rupSet1File = new File("/tmp/GarlockPintoMtnFix_RupSet.zip");
////		File rupSet2File = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/" +
////				"scratch/InversionSolutions/FM3_1_ZENG_Shaw09Mod_DsrTap_CharConst_M5Rate8.7" +
////				"_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
//		boolean oldRups = true;
//		
//		File diffFile;
//		if (oldRups)
//			diffFile = new File("/tmp/garlockOldRups.zip");
//		else
//			diffFile = new File("/tmp/garlockNewRups.zip");
//		
//		FaultSystemRupSet rupSet1 = SimpleFaultSystemRupSet.fromZipFile(rupSet1File);
////		FaultSystemRupSet rupSet2 = SimpleFaultSystemRupSet.fromZipFile(rupSet2File);
//		FaultSystemRupSet rupSet2 = InversionFaultSystemRupSetFactory.forBranch(FaultModels.FM3_1);
//		if (oldRups) {
//			FaultSystemRupSet tmp = rupSet1;
//			rupSet1 = rupSet2;
//			rupSet2 = tmp;
//		}
//		
//		writeDiffs(diffFile, rupSet1, rupSet2);
		
		LaughTestFilter laughTest = LaughTestFilter.getDefault();
//		laughTest.setCoulombFilter(null);
//		laughTest.setMaxCmlRakeChange(Double.POSITIVE_INFINITY);
		
//		LaughTestFilter.USE_BUGGY_COULOMB = false;
//		CoulombRatesTester.BUGGY_MIN_STRESS = false;
//		CumulativeAzimuthChangeFilter.USE_BUGGY_AZ_CHANGE = false;
//		AzimuthChangeFilter.INCLUDE_UCERF3p3_NEW_LL = true;
//		laughTest.setAllowSingleSectDuringJumps(true);
//		laughTest.getCoulombFilter().setMinIndividualProb(0.1);
//		laughTest.getCoulombFilter().setMinAverageProb(0.1);
//		laughTest.getCoulombFilter().setMinimumStressExclusionCeiling(1.5);
		Stopwatch watch = Stopwatch.createStarted();
		InversionFaultSystemRupSet rupSet1 = InversionFaultSystemRupSetFactory.forBranch(laughTest, 0.1, FaultModels.FM3_1);
		watch.stop();
		FaultSystemIO.writeRupSet(rupSet1, new File("/tmp/rupSet1.zip"));
		double secsNew = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
		rupSet1.setInfoString("");
//		laughTest.clearLaughTests();
//		laughTest = LaughTestFilter.getUCERF3p2Filter();
//		laughTest.setMaxCmlAzimuthChange(Double.POSITIVE_INFINITY);
//		LaughTestFilter.USE_BUGGY_COULOMB = false;
//		laughTest.getCoulombFilter().setMinAverageProb(0.04d);
//		laughTest.getCoulombFilter().setMinIndividualProb(0.04d);
//		laughTest.setMaxAzimuthChange(90d);
//		CoulombRatesTester.BUGGY_MIN_STRESS = false;
//		CumulativeAzimuthChangeFilter.USE_BUGGY_AZ_CHANGE = false;
//		AzimuthChangeFilter.INCLUDE_UCERF3p3_NEW_LL = false;
//		laughTest.setAllowSingleSectDuringJumps(true);
//		laughTest.getLaughTest(AzimuthChangeFilter.class).setTotAzChangeAtJunctionsOnly(true);
//		SectionCluster.NEW_ADD_RUPS = false;
//		CoulombRatesTester.BUGGY_MIN_STRESS = true;
//		laughTest.setAllowSingleSectDuringJumps(false);
		watch.reset();
		watch.start();
		InversionFaultSystemRupSet rupSet2 = InversionFaultSystemRupSetFactory.forBranch(laughTest, 0.1, FaultModels.FM3_1);
		watch.stop();
		double secsOld = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
		
		System.out.println("New method: "+rupSet1.getNumRuptures()+" rups ("+(float)secsNew+" s)");
		System.out.println("Old method: "+rupSet2.getNumRuptures()+" rups ("+(float)secsOld+" s)");
		
		writeDiffs(new File("/tmp/new_rups_in.zip"), rupSet1, rupSet2);
		writeDiffs(new File("/tmp/new_rups_out.zip"), rupSet2, rupSet1);
	}

	public static void writeDiffs(File diffFile, InversionFaultSystemRupSet rupSet1,
			FaultSystemRupSet rupSet2) throws IOException {
		HashSet<Rup> rups2 = new HashSet<Rup>();
		for (int r=0; r<rupSet2.getNumRuptures(); r++) {
			rups2.add(new Rup(rupSet2.getSectionsIndicesForRup(r)));
		}
		
		List<Integer> newRups = Lists.newArrayList();
		
		for (int r=0; r<rupSet1.getNumRuptures(); r++) {
			Rup rup = new Rup(rupSet1.getSectionsIndicesForRup(r));
			if (!rups2.contains(rup))
				newRups.add(r);
		}
		
		System.out.println("Found "+newRups.size()+" new rups ("
				+rupSet1.getNumRuptures()+" => "+rupSet2.getNumRuptures()+")");
		
		if (newRups.isEmpty()) {
			System.out.println("No rups!");
			return;
		}
		
//		// verify
//		for (Integer r : newRups) {
//			HashSet<Integer> rup = new HashSet<Integer>(rupSet1.getSectionsIndicesForRup(r));
//			rup2loop:
//			for (int r2=0; r2<rupSet2.getNumRuptures(); r2++) {
//				List<Integer> rup2 = rupSet2.getSectionsIndicesForRup(r2);
//				if (rup2.size() == rup.size()) {
//					for (Integer s : rup2)
//						if (!rup.contains(s))
//							continue rup2loop;
//					System.out.println("Equals ? "+new Rup(rupSet1.getSectionsIndicesForRup(r))
//								.equals(new Rup(rupSet2.getSectionsIndicesForRup(r2))));
//					throw new IllegalStateException("Found a match, wtf??");
//				}
//			}
//		}
		
		double[] mags = new double[newRups.size()];
		double[] rupAveSlips = new double[newRups.size()];
		double[] rakes = new double[newRups.size()];
		double[] rupAreas = new double[newRups.size()];
		double[] rupLenghts = new double[newRups.size()];
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		
		for (int i=0; i<newRups.size(); i++) {
			int r = newRups.get(i);
			mags[i] = rupSet1.getMagForRup(r);
			rupAveSlips[i] = rupSet1.getAveSlipForRup(r);
			rakes[i] = rupSet1.getAveRakeForRup(r);
			rupAreas[i] = rupSet1.getAreaForRup(r);
			rupLenghts[i] = rupSet1.getLengthForRup(r);
			sectionForRups.add(rupSet1.getSectionsIndicesForRup(r));
		}
		
		FaultSystemRupSet diffSet = new FaultSystemRupSet(
				rupSet1.getFaultSectionDataList(), rupSet1.getSlipRateForAllSections(),
				rupSet1.getSlipRateStdDevForAllSections(), rupSet1.getAreaForAllSections(),
				sectionForRups, mags, rakes, rupAreas, rupLenghts,
				rupSet1.getInfoString());
		diffSet = new InversionFaultSystemRupSet(diffSet, rupSet1.getLogicTreeBranch(),
				null, rupAveSlips, rupSet1.getCloseSectionsListList(),
				rupSet1.getRupturesForClusters(), rupSet1.getSectionsForClusters());
		
		FaultSystemIO.writeRupSet(diffSet, diffFile);
	}
	
	private static class Rup {
		private List<Integer> sects;
		public Rup(List<Integer> sects) {
			this.sects = Lists.newArrayList(sects);
			Collections.sort(this.sects);
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((sects == null) ? 0 : sects.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Rup other = (Rup) obj;
			if (sects == null) {
				if (other.sects != null)
					return false;
			} else if (!sects.equals(other.sects))
				return false;
			return true;
		}
	}

}
