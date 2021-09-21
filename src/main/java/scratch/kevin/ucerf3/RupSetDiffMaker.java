package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;

import scratch.UCERF3.U3FaultSystemRupSet;
import scratch.UCERF3.SlipEnabledRupSet;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.laughTest.AzimuthChangeFilter;
import scratch.UCERF3.inversion.laughTest.UCERF3PlausibilityConfig;
import scratch.UCERF3.utils.U3FaultSystemIO;

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
////		File rupSet1File = new File("/tmp/GarlockPintoMtnFix_RupSet.zip");
//////		File rupSet2File = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/" +
//////				"scratch/InversionSolutions/FM3_1_ZENG_Shaw09Mod_DsrTap_CharConst_M5Rate8.7" +
//////				"_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
////		boolean oldRups = true;
////		
////		File diffFile;
////		if (oldRups)
////			diffFile = new File("/tmp/garlockOldRups.zip");
////		else
////			diffFile = new File("/tmp/garlockNewRups.zip");
////		
////		FaultSystemRupSet rupSet1 = SimpleFaultSystemRupSet.fromZipFile(rupSet1File);
//////		FaultSystemRupSet rupSet2 = SimpleFaultSystemRupSet.fromZipFile(rupSet2File);
////		FaultSystemRupSet rupSet2 = InversionFaultSystemRupSetFactory.forBranch(FaultModels.FM3_1);
////		if (oldRups) {
////			FaultSystemRupSet tmp = rupSet1;
////			rupSet1 = rupSet2;
////			rupSet2 = tmp;
////		}
////		
////		writeDiffs(diffFile, rupSet1, rupSet2);
//		
//		LaughTestFilter laughTest = LaughTestFilter.getDefault();
////		laughTest.setCoulombFilter(null);
////		laughTest.setMaxCmlRakeChange(Double.POSITIVE_INFINITY);
//		
////		LaughTestFilter.USE_BUGGY_COULOMB = false;
////		CoulombRatesTester.BUGGY_MIN_STRESS = false;
////		CumulativeAzimuthChangeFilter.USE_BUGGY_AZ_CHANGE = false;
////		AzimuthChangeFilter.INCLUDE_UCERF3p3_NEW_LL = true;
////		laughTest.setAllowSingleSectDuringJumps(true);
////		laughTest.getCoulombFilter().setMinIndividualProb(0.1);
////		laughTest.getCoulombFilter().setMinAverageProb(0.1);
////		laughTest.getCoulombFilter().setMinimumStressExclusionCeiling(1.5);
//		Stopwatch watch = Stopwatch.createStarted();
//		InversionFaultSystemRupSet rupSet1 = InversionFaultSystemRupSetFactory.forBranch(laughTest, 0.1, FaultModels.FM3_1);
//		watch.stop();
//		FaultSystemIO.writeRupSet(rupSet1, new File("/tmp/rupSet1.zip"));
//		double secsNew = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
//		rupSet1.setInfoString("");
////		laughTest.clearLaughTests();
////		laughTest = LaughTestFilter.getUCERF3p2Filter();
////		laughTest.setMaxCmlAzimuthChange(Double.POSITIVE_INFINITY);
////		LaughTestFilter.USE_BUGGY_COULOMB = false;
////		laughTest.getCoulombFilter().setMinAverageProb(0.04d);
////		laughTest.getCoulombFilter().setMinIndividualProb(0.04d);
////		laughTest.setMaxAzimuthChange(90d);
////		CoulombRatesTester.BUGGY_MIN_STRESS = false;
////		CumulativeAzimuthChangeFilter.USE_BUGGY_AZ_CHANGE = false;
////		AzimuthChangeFilter.INCLUDE_UCERF3p3_NEW_LL = false;
////		laughTest.setAllowSingleSectDuringJumps(true);
////		laughTest.getLaughTest(AzimuthChangeFilter.class).setTotAzChangeAtJunctionsOnly(true);
////		SectionCluster.NEW_ADD_RUPS = false;
////		CoulombRatesTester.BUGGY_MIN_STRESS = true;
////		laughTest.setAllowSingleSectDuringJumps(false);
//		watch.reset();
//		watch.start();
//		InversionFaultSystemRupSet rupSet2 = InversionFaultSystemRupSetFactory.forBranch(laughTest, 0.1, FaultModels.FM3_1);
//		watch.stop();
//		double secsOld = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
		
//		File rsDir = new File("/home/kevin/Simulators/catalogs/rundir4983_stitched/fss");
//		File rsRupsFile = new File(rsDir, "rsqsim_sol_m6.5_skip5000_sectArea0.2.zip");
//		File u3File = new File("/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/"
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
//		
//		FaultSystemRupSet rsRupSet = FaultSystemIO.loadRupSet(rsRupsFile);
//		FaultSystemRupSet u3RupSet = FaultSystemIO.loadRupSet(u3File);
//		
//		System.out.println("RSQSim has: "+rsRupSet.getNumRuptures()+" rups");
//		System.out.println("UCERF3 has: "+u3RupSet.getNumRuptures()+" rups");
//		
//		writeDiffs(new File(rsDir, "new_rups_in.zip"), rsRupSet, u3RupSet);
//		writeDiffs(new File(rsDir, "new_rups_out.zip"), u3RupSet, rsRupSet);
		
//		File testFile = new File("/tmp/test_rup_set.zip");
//		File u3File = new File("/home/kevin/workspace/opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/"
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		File rupSetsDir = new File("/home/kevin/OpenSHA/UCERF4/rup_sets");
		File refRupSetFile = new File(rupSetsDir, "fm3_1_plausibleMulti10km_direct_slipP0.05incr_cff0.75IntsPos_comb2Paths_cffFavP0.02_cffFavRatioN2P0.5_sectFractPerm0.05.zip");
		File testRupSetFile = new File(rupSetsDir, "fm3_1_plausibleMulti10km_direct_slipP0.05incr_cff0.75IntsPos_comb2Paths_cffFavP0.02_cffFavRatioN2P0.5_sectFractPerm0.05_comp/add_CumulativeAzimuth.zip");
		
		U3FaultSystemRupSet testRupSet = U3FaultSystemIO.loadRupSet(testRupSetFile);
		U3FaultSystemRupSet refRupSet = U3FaultSystemIO.loadRupSet(refRupSetFile);
		
		System.out.println("Test has: "+testRupSet.getNumRuptures()+" rups");
		System.out.println("Ref has: "+refRupSet.getNumRuptures()+" rups");
		
		writeDiffs(new File("/tmp/test_new_rups_in.zip"), testRupSet, refRupSet);
		writeDiffs(new File("/tmp/test_new_rups_out.zip"), refRupSet, testRupSet);
	}

	public static void writeDiffs(File diffFile, U3FaultSystemRupSet rupSet1,
			U3FaultSystemRupSet rupSet2) throws IOException {
		HashSet<UniqueRupture> rups2 = new HashSet<UniqueRupture>();
		for (int r=0; r<rupSet2.getNumRuptures(); r++) {
			rups2.add(UniqueRupture.forIDs(rupSet2.getSectionsIndicesForRup(r)));
		}

		HashSet<UniqueRupture> newRupsSet = new HashSet<UniqueRupture>();
		List<Integer> newRups = Lists.newArrayList();
		
		for (int r=0; r<rupSet1.getNumRuptures(); r++) {
			UniqueRupture rup = UniqueRupture.forIDs(rupSet1.getSectionsIndicesForRup(r));
			if (!rups2.contains(rup) && !newRupsSet.contains(rup)) {
				newRups.add(r);
				newRupsSet.add(rup);
			}
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
		double[] rupLenghts = rupSet1.getLengthForAllRups() == null || rupSet2.getLengthForAllRups() == null
				? null : new double[newRups.size()];
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		
		for (int i=0; i<newRups.size(); i++) {
			int r = newRups.get(i);
			mags[i] = rupSet1.getMagForRup(r);
			if (rupSet1 instanceof SlipEnabledRupSet)
				rupAveSlips[i] = ((SlipEnabledRupSet)rupSet1).getAveSlipForRup(r);
			rakes[i] = rupSet1.getAveRakeForRup(r);
			rupAreas[i] = rupSet1.getAreaForRup(r);
			if (rupLenghts != null)
				rupLenghts[i] = rupSet1.getLengthForRup(r);
			sectionForRups.add(rupSet1.getSectionsIndicesForRup(r));
		}
		
		U3FaultSystemRupSet diffSet = new U3FaultSystemRupSet(
				rupSet1.getFaultSectionDataList(), rupSet1.getSlipRateForAllSections(),
				rupSet1.getSlipRateStdDevForAllSections(), rupSet1.getAreaForAllSections(),
				sectionForRups, mags, rakes, rupAreas, rupLenghts,
				rupSet1.getInfoString());
		if (rupSet1 instanceof InversionFaultSystemRupSet) {
			InversionFaultSystemRupSet invSet = (InversionFaultSystemRupSet)rupSet1;
			diffSet = new InversionFaultSystemRupSet(diffSet, invSet.getLogicTreeBranch(),
				null, rupAveSlips, invSet.getCloseSectionsListList(),
				invSet.getRupturesForClusters(), invSet.getSectionsForClusters());
		}
		
		U3FaultSystemIO.writeRupSet(diffSet, diffFile);
	}

}
