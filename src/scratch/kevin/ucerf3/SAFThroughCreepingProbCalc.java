package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.util.DataUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class SAFThroughCreepingProbCalc {
	
	private static void calcForRay(FaultSystemSolution sol) throws IOException {
		final FaultSystemRupSet rupSet = sol.getRupSet();
		
		List<HashSet<Integer>> rupsToConsider = Lists.newArrayList();
		List<String> setNames = Lists.newArrayList();
		
		// paleo sites
		int[] parentIDs = new int[] {32, 285, 300, 287, 286, 301, 282, 283, 284, 295, 294}; // parkfield to coachella
		HashSet<Integer> parentIDSet = new HashSet<Integer>(Ints.asList(parentIDs));
		for (PaleoRateConstraint constr :
			UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList())) {
			int sectID = constr.getSectionIndex();
			int parentID = rupSet.getFaultSectionData(sectID).getParentSectionId();
			if (!parentIDSet.contains(parentID))
				continue;
			
			rupsToConsider.add(new HashSet<Integer>(rupSet.getRupturesForSection(sectID)));
			setNames.add(constr.getPaleoSiteName().replaceAll("\t", "").trim());
		}
		
//		HashSet<Integer> rups = new HashSet<Integer>();
////		for (int parentID : FaultModels.FM3_1.getNamedFaultsMapAlt().get("San Andreas"))
////			rups.addAll(rupSet.getRupturesForParentSection(parentID));
//		int[] parentIDs = new int[] {32, 285, 300, 287, 286, 301, 282, 283, 284, 295, 294}; // parkfield to coachella
////		int[] parentIDs = new int[] {32, 285, 300, 287, 286, 301}; // parkfield to mojave s
//		for (int parentID : parentIDs)
//			rups.addAll(rupSet.getRupturesForParentSection(parentID));
//		rupsToConsider.add(rups);
//		setNames.add("SSAF");
		
		for (int i=0; i<rupsToConsider.size(); i++) {
			double safRate = 0d;
			List<Double> magsAbove7 = Lists.newArrayList();
			
			HashSet<Integer> rupsSet = rupsToConsider.get(i);
			System.out.println("*** "+setNames.get(i)+" ***");
			List<Integer> rupsSorted = Lists.newArrayList(rupsSet);
			// sort by magnitude
			Collections.sort(rupsSorted, new Comparator<Integer>() {

				@Override
				public int compare(Integer o1, Integer o2) {
					return Double.compare(rupSet.getMagForRup(o1), rupSet.getMagForRup(o2));
				}
			});
			
			for (int r : rupsSorted) {
				double mag = rupSet.getMagForRup(r);
				if (mag >= 7d) {
					safRate += sol.getRateForRup(r);
					magsAbove7.add(mag);
				}
			}
			
			// now find weighted median
			double halfRate = safRate * 0.5;
			double cmlRate = 0d;
			for (int r : rupsSorted) {
				if (rupSet.getMagForRup(r) >= 7d)
					cmlRate += sol.getRateForRup(r);
				if (cmlRate >= halfRate) {
					System.out.println("Weighted Median Mag: "+(float)rupSet.getMagForRup(r));
					break;
				}
			}
			System.out.println("Unweighted median: "+(float)DataUtils.median(Doubles.toArray(magsAbove7)));
			System.out.println("M >= 7 rate: "+(float)safRate);
		}
	}

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution meanSol = FaultSystemIO.loadSol(
				new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = meanSol.getRupSet();
		
		// quick calc for Ray
//		calcForRay(meanSol);
		
		double[] minMags = { 6.5d, 7d, 7.5, 8d };
		int duration = 30;
		
//		int peninsulaParent = 655;
//		int carrizoParent = 300;

		int peninsulaParent = 653; // actually SAF offshore
//		int peninsulaParent = 13; // actually mendocino
		int carrizoParent = 295; // actually coachella
		
		HashSet<Integer> peninsulaRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(peninsulaParent));
		HashSet<Integer> carrizoRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(carrizoParent));
		
		HashSet<Integer> bothRups = new HashSet<Integer>();
		for (Integer peninsulaRup : peninsulaRups)
			if (carrizoRups.contains(peninsulaRup))
				bothRups.add(peninsulaRup);
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(meanSol);
		
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setDuration(duration);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.updateForecast();
		
		for (double minMag : minMags) {
			int peninsulaAbove = 0;
			int carrizoAbove = 0;
			int bothAbove = 0;
			
			List<Double> peninsulaProbs = Lists.newArrayList();
			List<Double> carrizoProbs = Lists.newArrayList();
			List<Double> bothProbs = Lists.newArrayList();
			
			for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
				int rupIndex = erf.getFltSysRupIndexForSource(sourceID);
				double prob = erf.getSource(sourceID).computeTotalProbAbove(minMag);
				if (prob == 0)
					continue;
				if (peninsulaRups.contains(rupIndex)) {
					peninsulaProbs.add(prob);
					peninsulaAbove++;
				}
				if (carrizoRups.contains(rupIndex)) {
					carrizoProbs.add(prob);
					carrizoAbove++;
				}
				if (bothRups.contains(rupIndex)) {
					bothProbs.add(prob);
					bothAbove++;
				}
			}
			
			double peninsulaProb = FaultSysSolutionERF_Calc.calcSummedProbs(peninsulaProbs);
			double carrizoProb = FaultSysSolutionERF_Calc.calcSummedProbs(carrizoProbs);
			double bothProb = FaultSysSolutionERF_Calc.calcSummedProbs(bothProbs);
			
			System.out.println("Peninsula "+duration+"yr M>="+(float)minMag+" prob: "+(float)peninsulaProb
					+" ("+peninsulaAbove+" rups)");
			System.out.println("Carrizo "+duration+"yr M>="+(float)minMag+" prob: "+(float)carrizoProb
					+" ("+carrizoAbove+" rups)");
			System.out.println("Together "+duration+"yr M>="+(float)minMag+" prob: "+(float)bothProb
					+" ("+bothAbove+" rups)");
			System.out.println();
		}
	}

}

