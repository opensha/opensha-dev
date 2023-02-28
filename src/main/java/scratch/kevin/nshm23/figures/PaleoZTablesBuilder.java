package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;

import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint.U3AveSlipProbModel;

class PaleoZTablesBuilder {

	public static void main(String[] args) throws IOException {
		File u3SolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged_full_modules.zip");
		File nshm23MethodologyDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_06-nshm23_u3_hybrid_branches-no_paleo_slip-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File nshm23ModelDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_12_07-nshm23_branches-no_paleo_slip-mod_dm_weights-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("Model", "Average z-score", "Average Absolute z-score");
		DecimalFormat zDF = new DecimalFormat("0.00");
		
		List<FaultSystemSolution> caSols = new ArrayList<>();
		List<String> caNames = new ArrayList<>();
		
		caSols.add(FaultSystemSolution.load(u3SolFile));
		caNames.add("UCERF3");
		
		caSols.add(FaultSystemSolution.load(new File(nshm23MethodologyDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip")));
		caNames.add("NSHM23 Methodology");
		
		for (NSHM23_PaleoUncertainties uncert : NSHM23_PaleoUncertainties.values()) {
			File solFile = new File(nshm23MethodologyDir, "node_branch_averaged/PaleoUncert_"+uncert.getFilePrefix()+".zip");
			caSols.add(FaultSystemSolution.load(solFile));
			caNames.add("NSHM23 Methodology, "+uncert.getShortName().replace("F", "-F"));
		}
		
		List<FaultSystemSolution> wasatchSols = new ArrayList<>();
		List<String> wasatchNames = new ArrayList<>();
		
		caSols.add(FaultSystemSolution.load(new File(nshm23ModelDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip")));
		caNames.add("NSHM23 CA");
		wasatchSols.add(caSols.get(caSols.size()-1));
		wasatchNames.add("NSHM23 Wasatch");
		
		for (NSHM23_PaleoUncertainties uncert : NSHM23_PaleoUncertainties.values()) {
			File solFile = new File(nshm23ModelDir, "node_branch_averaged/PaleoUncert_"+uncert.getFilePrefix()+".zip");
			caSols.add(FaultSystemSolution.load(solFile));
			caNames.add("NSHM23 CA, "+uncert.getShortName().replace("F", "-F"));
			wasatchSols.add(caSols.get(caSols.size()-1));
			wasatchNames.add("NSHM23 Wasatch, "+uncert.getShortName().replace("F", "-F"));
		}
		
		Region ca = new CaliforniaRegions.RELM_TESTING();
		
		for (int i=0; i<caSols.size(); i++) {
			FaultSystemSolution sol = caSols.get(i);
			double avg = avgZScore(sol, false, ca);
			double avgAbs = avgZScore(sol, true, ca);
			
			csv.addLine(caNames.get(i), zDF.format(avg), zDF.format(avgAbs));
		}
		
		Region ut = NSHM23_SingleStates.UT.loadRegion();
		for (int i=0; i<wasatchSols.size(); i++) {
			FaultSystemSolution sol = wasatchSols.get(i);
			double avg = avgZScore(sol, false, ut);
			double avgAbs = avgZScore(sol, true, ut);
			
			csv.addLine(wasatchNames.get(i), zDF.format(avg), zDF.format(avgAbs));
		}
		
		csv.writeToFile(new File(outputDir, "paleo_z_scores.csv"));
		
		// now paleo slips (only U3)
		csv = new CSVFile<>(true);
		
		csv.addLine("Model", "Average z-score", "Average Absolute z-score");
		
		FaultSystemSolution u3Sol = caSols.get(0);
		double avg = avgPaleoSlipZScore(u3Sol, false);
		double avgAbs = avgPaleoSlipZScore(u3Sol, true);
		csv.addLine("UCERF3", zDF.format(avg), zDF.format(avgAbs));
		csv.writeToFile(new File(outputDir, "paleo_slip_z_scores.csv"));
	}
	
	private static double avgZScore(FaultSystemSolution sol, boolean abs, Region region) {
		PaleoseismicConstraintData data = sol.getRupSet().requireModule(PaleoseismicConstraintData.class);
		
		List<? extends SectMappedUncertainDataConstraint> origConstraints = data.getPaleoRateConstraints();
		List<SectMappedUncertainDataConstraint> reginalConstraints = new ArrayList<>();
		for (SectMappedUncertainDataConstraint constraint : origConstraints)
			if (region.contains(constraint.dataLocation) ||
					region.contains(sol.getRupSet().getFaultSectionData(constraint.sectionIndex).getFaultTrace().first()))
				reginalConstraints.add(constraint);
		
		double avg = 0d;
		for (SectMappedUncertainDataConstraint constraint : reginalConstraints) {
			double solRate = sol.calcTotPaleoVisibleRateForSect(constraint.sectionIndex, data.getPaleoProbModel());
			
			double z = (solRate - constraint.bestEstimate)/constraint.getPreferredStdDev();
			if (abs)
				avg += Math.abs(z);
			else
				avg += z;
		}
		avg /= reginalConstraints.size();
		return avg;
	}
	
	private static double avgPaleoSlipZScore(FaultSystemSolution sol, boolean abs) throws IOException {
		FaultSystemRupSet rupSet = sol.getRupSet();
		AveSlipModule aveSlipModule = rupSet.requireModule(AveSlipModule.class);
		SlipAlongRuptureModel slipAlongModule = rupSet.requireModule(SlipAlongRuptureModel.class);
		// use original U3 data
		List<U3AveSlipConstraint> rawConstraints = U3AveSlipConstraint.load(sol.getRupSet().getFaultSectionDataList());
		List<SectMappedUncertainDataConstraint> constraints =
				PaleoseismicConstraintData.inferRatesFromSlipConstraints(rupSet, rawConstraints, false);
		U3AveSlipProbModel probModel = new U3AveSlipProbModel();
		
		double avg = 0d;
		for (SectMappedUncertainDataConstraint constraint : constraints) {
//			sol.calctot
			double solRate = 0d;
			for (int rupIndex : sol.getRupSet().getRupturesForSection(constraint.sectionIndex)) {
				int sectIndexInRup = rupSet.getSectionsIndicesForRup(rupIndex).indexOf(constraint.sectionIndex);
				double slip = slipAlongModule.calcSlipOnSectionsForRup(rupSet, aveSlipModule, rupIndex)[sectIndexInRup];
				solRate += sol.getRateForRup(rupIndex) * probModel.getProbabilityOfObservedSlip(slip);
			}
			
			double z = (solRate - constraint.bestEstimate)/constraint.getPreferredStdDev();
			if (abs)
				avg += Math.abs(z);
			else
				avg += z;
		}
		avg /= constraints.size();
		return avg;
	}

}
