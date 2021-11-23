package scratch.kevin.nshm23;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDLaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SubSectMFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.U3MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.MisfitStdDevCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.nshm23.targetMFDs.DraftModelConstraintBuilder;

public class BatchInversionScriptWriter {
	
	private static final DecimalFormat oDF = new DecimalFormat("0.###");

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteToalThreads = 20;
		int remoteTotalMemGB = 55;
		BatchScriptWriter scriptWrite = new USC_CARC_ScriptWriter();
		String queue = "scec";
		
		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
//		String dirName = "2021_11_09";
		
		List<InversionConfiguration> configs = new ArrayList<>();
		List<String> subDirNames = new ArrayList<>();
		
		File rsDir = new File("/home/kevin/markdown/inversions/");
		
		File rupSetFile;
		String rsPrefix;
		
		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_reproduce_ucerf3.zip");
		rsPrefix = "reproduce-ucerf3-ref_branch-uniform";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_coulomb.zip");
//		rsPrefix = "coulomb-ref_branch-uniform";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_geol_uniform_reproduce_ucerf3.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-geol-uniform";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_tapered_reproduce_ucerf3.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-tapered";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_reproduce_ucerf3_fractGrow0.1.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-uniform-grow0.1";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_tapered_reproduce_ucerf3.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-tapered";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_coulomb.zip");
//		rsPrefix = "coulomb-fm31-ref_branch-uniform";
		
//		File remoteMeanCompFile = new File(remoteMainDir,
//				"2021_10_18-reproduce-ucerf3-ref_branch-uniform-new_anneal-5x_avg-try_zero-var_perturb-noWL-5h/mean_solution.zip");
//		String remoteMeanCompareName = "U3-New-Anneal-NoWL";
//		File remoteMeanCompFile = new File(remoteMainDir,
//				"2021_10_25-reproduce-ucerf3-ref_branch-uniform-new_anneal-uncert_weighted-mfd_sd_0.1-minimize10000-smooth1000-5h/mean_solution.zip");
//		String remoteMeanCompareName = "U3-Uncert-Wtd";
		File remoteMeanCompFile = null;
		String remoteMeanCompareName = null;
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(rupSetFile);
		
		// this is applied to the first configuration supplied
//		PlotLevel primaryPlotLevel = PlotLevel.FULL;
		PlotLevel primaryPlotLevel = PlotLevel.DEFAULT;
//		PlotLevel primaryPlotLevel = PlotLevel.LIGHT;
		// this is applied to everything else
//		PlotLevel allPlotLevel = PlotLevel.FULL;
//		PlotLevel allPlotLevel = PlotLevel.DEFAULT;
//		PlotLevel allPlotLevel = PlotLevel.LIGHT;
		PlotLevel allPlotLevel = null;
		// this is applied to the average job (if avgJob==true)
		PlotLevel avgPlotLevel = PlotLevel.FULL;
		
		boolean avgJob = false;
		
		/*
		 * UCERF3 as was, 10 times
		 */
////		dirName += "-"+rsPrefix+"-u3Iters";
//		dirName += "-reproduce-ucerf3-ref_branch-uniform-40h";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints,
//				TimeCompletionCriteria.getInHours(40))
////				new IterationCompletionCriteria(22088044l))
//				.subCompletion(new IterationCompletionCriteria(1227))
//				.threads(5)
//				.nonNegativity(NonnegativityConstraintType.LIMIT_ZERO_RATES)
//				.perturbation(GenerationFunctionType.UNIFORM_0p001)
//				.waterLevel(u3Gen.getWaterLevelRates()).build();
//		for (int i=0; i<10; i++) {
//			configs.add(config);
//			subDirNames.add("u3_reproduce_run_"+i);
//		}
//		avgJob = true;
		
		/*
		 * UCERF3 convergence test, 200 times
		 */
//		dirName += "-"+rsPrefix+"-convergence-u3Iters";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints,
////				TimeCompletionCriteria.getInHours(5))
//				new IterationCompletionCriteria(22088044l))
//				.subCompletion(new IterationCompletionCriteria(1227))
//				.threads(5)
//				.nonNegativity(NonnegativityConstraintType.LIMIT_ZERO_RATES)
//				.perturbation(GenerationFunctionType.UNIFORM_0p001)
//				.waterLevel(u3Gen.getWaterLevelRates()).build();
//		for (int i=0; i<200; i++) {
//			configs.add(config);
//			subDirNames.add("u3_converge_run_"+i);
//		}
//		avgJob = true;
		
		/*
		 * UCERF3 as was with different thread counts
		 */
//		dirName += "-"+rsPrefix+"-thread_test";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints,
//				TimeCompletionCriteria.getInHours(5))
////				new IterationCompletionCriteria(22088044l))
////				.subCompletion(new IterationCompletionCriteria(1227))
//				.threads(5)
//				.nonNegativity(NonnegativityConstraintType.LIMIT_ZERO_RATES)
//				.perturbation(GenerationFunctionType.UNIFORM_0p001)
//				.waterLevel(u3Gen.getWaterLevelRates()).build();
//		for (int i=1; i<=remoteToalThreads; i++) {
//			configs.add(InversionConfiguration.builder(config).threads(i).build());
//			subDirNames.add("u3_reproduce_"+i+"_threads");
//		}
//		avgJob = false;
		
		/*
		 * new annealing defaults, 10 times
		 */
//		dirName += "-"+rsPrefix+"-u3_constr-new_anneal";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		
//		dirName += "-no_paleo";
//		removeConstraintsByType(u3Constraints, PaleoSlipInversionConstraint.class);
//		removeConstraintsByType(u3Constraints, PaleoRateInversionConstraint.class);
//		
//		dirName += "-no_parkfield";
//		removeConstraintsByType(u3Constraints, ParkfieldInversionConstraint.class);
//		
//		dirName += "-no_u2_ss_mfds";
//		removeConstraintsByType(u3Constraints, U3MFDSubSectNuclInversionConstraint.class);
//		
//		dirName += "-no_smooth";
//		removeConstraintsByType(u3Constraints, MFDLaplacianSmoothingInversionConstraint.class);
//		
//		dirName += "-single_mfd_region";
//		for (int c=0; c<u3Constraints.size(); c++) {
//			if (u3Constraints.get(c) instanceof MFDInversionConstraint) {
//				MFDInversionConstraint orig = (MFDInversionConstraint)u3Constraints.get(c);
//				List<? extends IncrementalMagFreqDist> mfds = orig.getMFDs();
//				Preconditions.checkState(mfds.size() == 2);
//				SummedMagFreqDist sumMFD = null;
//				for (IncrementalMagFreqDist mfd : mfds) {
//					if (sumMFD == null)
//						sumMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
//					sumMFD.addIncrementalMagFreqDist(mfd);
//				}
//				sumMFD.setRegion(new CaliforniaRegions.RELM_TESTING());
//				u3Constraints.set(c, new MFDInversionConstraint(rupSet, orig.getWeight(), orig.isInequality(),
//						orig.getWeightingType(), List.of(sumMFD), orig.getExcludeRupIndexes()));
//			}
//		}
//		
//		System.out.println("Constraint list:");
//		for (InversionConstraint constr : u3Constraints)
//			System.out.println("\t"+constr.getName()+": wt="+(float)constr.getWeight()+"\twtType="+constr.getWeightingType());
//		
////		boolean waterLevel = true; dirName += "-u3WL";
//		boolean waterLevel = false; dirName += "-noWL";
//		
////		CompletionCriteria completion = TimeCompletionCriteria.getInHours(10); dirName += "-10h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5); dirName += "-5h";
////		CompletionCriteria completion = TimeCompletionCriteria.getInHours(2); dirName += "-2h";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(5);
////		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(2);
////		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(1);
////		CompletionCriteria avgCompletion = null; dirName += "-noAvg";
//		
//		InversionConfiguration.Builder builder = InversionConfiguration.builder(
//				u3Constraints, completion).threads(remoteToalThreads);
//		
//		if (avgCompletion != null)
//			builder.avgThreads(remoteToalThreads/4, avgCompletion);
//		if (waterLevel)
//			builder.waterLevel(u3Gen.getWaterLevelRates());
////				.nonNegativity(NonnegativityConstraintType.LIMIT_ZERO_RATES)
//		
//		InversionConfiguration config = builder.build();
//		for (int i=0; i<10; i++) {
//			configs.add(config);
//			subDirNames.add("u3_new_anneal_run_"+i);
//		}
//		avgJob = true;
		
		/*
		 * new NSHM23 draft scheme
		 */
		dirName += "-"+rsPrefix+"-nshm23_draft";
		double bVal = 0.8;
		dirName += "-supra_b_"+(float)bVal;
		
		double slipWeight=0d, paleoWeight=0d, parkWeight=0d, mfdWeight=0d, nuclWeight=0d, nuclMFDWeight=0d, paleoSmoothWeight=0d;
		
		// weights are zero (disabled) unless uncommented here
		slipWeight = 1d;
		paleoWeight = 5;
		parkWeight = 100;
		mfdWeight = 10;
		nuclWeight = 0.5;
		nuclMFDWeight = 0.01;
		paleoSmoothWeight = paleoWeight > 0 ? 10000 : 0;
		
		boolean skipBelow = true;
		boolean applyDefModelUncertaintiesToNucl = true;
		boolean addSectCountUncertaintiesToMFD = false;
		boolean adjustForIncompatibleData = paleoWeight > 0 || parkWeight > 0;
		
//		List<InversionConstraint> u3Constraints = InversionsCLI.getU3Constraints(rupSet);
		
//		if (applyDefModelUncertaintiesToNucl)
//			dirName += "-dm_uncert_nucl";
		
		if (adjustForIncompatibleData && (mfdWeight > 0d || nuclWeight > 0d || nuclMFDWeight > 0d))
			dirName += "-adj_ucert_for_data";
		
		if (addSectCountUncertaintiesToMFD && mfdWeight > 0d)
			dirName += "-mfd_count_uncert";

//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_03-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-2h/mean_solution.zip");
//		remoteMeanCompareName = "All-New-Constr-b=0.8";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_08-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-only-slip-mfd-mfd_wt_10-skipBelow-2h/mean_solution.zip");
//		remoteMeanCompareName = "New-MFD-Constr-b=0.8-No-Nucl";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_18-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_"+(float)bVal+"-adj_ucert_for_data-paleo_wt_5-parkfield_wt_100-mfd_wt_10-sect_wt_0.5-smooth_paleo_wt_10000-skipBelow-2h/mean_solution.zip");
//		remoteMeanCompareName = "Tot-Sect-Nucl-Wt-Adjusted";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_18-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_"+(float)bVal+"-adj_ucert_for_data-paleo_wt_5-parkfield_wt_100-mfd_wt_10-no_sect_rate-sect_nucl_mfd_0.01-smooth_paleo_wt_10000-skipBelow-2h/mean_solution.zip");
//		remoteMeanCompareName = "Sect-Nucl-MFD-Wt-Adjusted";
		
		int num = 5;
		
		if (num == 1)
			primaryPlotLevel = avgPlotLevel;

		DraftModelConstraintBuilder constrBuilder = new DraftModelConstraintBuilder(rupSet, bVal,
				applyDefModelUncertaintiesToNucl, addSectCountUncertaintiesToMFD, adjustForIncompatibleData);
		
//		dirName += "-u3_supra_reduction";
//		constrBuilder.useExistingTargetSlipRates();
		
		if (slipWeight > 0d) {
			if (slipWeight != 1d)
				dirName += "-slip_wt_"+oDF.format(slipWeight);
			constrBuilder.slipRates().weight(slipWeight);
		} else {
			dirName += "-no_slip";
		}
		
		if (paleoWeight > 0d) {
			dirName += "-paleo_wt_"+oDF.format(paleoWeight);
			constrBuilder.paleoRates().weight(paleoWeight);
			constrBuilder.paleoSlips().weight(paleoWeight);
		} else {
			dirName += "-no_paleo";
		}
		
		if (parkWeight > 0d) {
			dirName += "-parkfield_wt_"+oDF.format(parkWeight);
			constrBuilder.parkfield().weight(parkWeight);
		} else {
			dirName += "-no_parkfield";
		}
		
		if (mfdWeight > 0d) {
			dirName += "-mfd_wt_"+oDF.format(mfdWeight);
			constrBuilder.supraBValMFDs().weight(mfdWeight);
		} else {
			dirName += "-no_mfd";
		}
		
		if (nuclWeight > 0d) {
			dirName += "-sect_wt_"+oDF.format(nuclWeight);
			constrBuilder.sectSupraRates().weight(nuclWeight);
		} else {
			dirName += "-no_sect_rate";
		}
		
		if (nuclMFDWeight > 0d) {
			dirName += "-sect_nucl_mfd_"+oDF.format(nuclMFDWeight);
			constrBuilder.sectSupraNuclMFDs().weight(nuclMFDWeight);
		} else {
//			dirName += "-no_sect_mfd";
		}
		
//		dirName += "-smooth_all";
//		constrBuilder.supraSmooth();
//		constrBuilder.weight(10000); dirName += "_wt_10000";
		
//		int parentWt = 10000;
//		dirName += "-parent_smooth_wt_"+parentWt;
//		constrBuilder.add(new ParentSectSmoothnessConstraint(rupSet, parentWt, true));
		
		if (paleoSmoothWeight > 0d) {
			dirName += "-smooth_paleo";
			constrBuilder.supraPaleoSmooth();
			constrBuilder.weight(paleoSmoothWeight); dirName += "_wt_"+oDF.format(paleoSmoothWeight);
		}
		
//		dirName += "-u3_target_mfds";
////		boolean consolidateRegion = true; dirName += "-single_mfd_region";
//		boolean consolidateRegion = false;
//		for (int c=0; c<u3Constraints.size(); c++) {
//			if (u3Constraints.get(c) instanceof MFDInversionConstraint) {
//				MFDInversionConstraint orig = (MFDInversionConstraint)u3Constraints.get(c);
//				List<? extends IncrementalMagFreqDist> mfds = orig.getMFDs();
//				Preconditions.checkState(mfds.size() == 2);
//				if (consolidateRegion) {
//					SummedMagFreqDist sumMFD = null;
//					for (IncrementalMagFreqDist mfd : mfds) {
//						if (sumMFD == null)
//							sumMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
//						sumMFD.addIncrementalMagFreqDist(mfd);
//					}
//					sumMFD.setRegion(new CaliforniaRegions.RELM_TESTING());
//					constrBuilder.add(new MFDInversionConstraint(rupSet, orig.getWeight(), orig.isInequality(),
//							orig.getWeightingType(), List.of(sumMFD), orig.getExcludeRupIndexes()));
//				} else {
//					constrBuilder.add(orig);
//				}
//			}
//		}
		
		IntegerPDF_FunctionSampler sampler = null;
		if (skipBelow) {
			sampler = constrBuilder.getSkipBelowMinSampler();
			dirName += "-skipBelow";
			constrBuilder.except(RupRateMinimizationConstraint.class);
		} else {
			dirName += "-minimizeBelow";
			constrBuilder.minimizeBelowSectMinMag();
		}
		
		List<InversionConstraint> constraints = constrBuilder.build();
		
		System.out.println("Constraint list:");
		for (InversionConstraint constr : constraints)
			System.out.println("\t"+constr.getName()+": wt="+(float)constr.getWeight()+"\twtType="+constr.getWeightingType());
		
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(40); dirName += "-40h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(20); dirName += "-20h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(10); dirName += "-10h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5); dirName += "-5h";
		CompletionCriteria completion = TimeCompletionCriteria.getInHours(2); dirName += "-2h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(1); dirName += "-1h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(30); dirName += "-30m";
//		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(20); dirName += "-20m";
//		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(10); dirName += "-10m";
		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(5);
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(2);
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(1);
//		CompletionCriteria avgCompletion = null; dirName += "-noAvg";
		
		if (num > 10)
			dirName += "-"+num+"x";
		
		InversionConfiguration.Builder builder = InversionConfiguration.builder(constraints, completion)
				.threads(remoteToalThreads).sampler(sampler);
		if (avgCompletion != null)
			builder.avgThreads(remoteToalThreads/4, avgCompletion);
		
//		dirName += "-sampler";
//		builder.sampler(Inversions.getDefaultVariablePerturbationBasis(rupSet));
		
//		dirName += "-sampler";
//		builder.sampler(Inversions.getDefaultVariablePerturbationBasis(rupSet));
		
//		dirName += "-initial";
//		builder.initialSolution(Inversions.getDefaultVariablePerturbationBasis(rupSet));
		
//		dirName += "-simple_exp_perturb";
//		builder.perturbation(GenerationFunctionType.EXPONENTIAL_SCALE);
		
		InversionConfiguration config = builder.build();
		for (int i=0; i<num; i++) {
			configs.add(config);
			subDirNames.add("run_"+i);
		}
		avgJob = true;
		allPlotLevel = null;
		
		/*
		 * U3 constraints, but swapping in new nucleation and smoothing scheme
		 */
//		dirName += "-"+rsPrefix+"-u3_constr-except";
//		//	dirName += "-"+rsPrefix+"-new_anneal-no_avg-try_zero-var_perturb-noWL-5h";
//		//	dirName += "-"+rsPrefix+"-new_anneal-no_avg-limit_zero-var_perturb-noWL-5h";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		
//		FaultSystemRupSet altRupSet = FaultSystemRupSet.buildFromExisting(rupSet, true).build();
//		
//		DraftModelConstraintBuilder newBuilder = new DraftModelConstraintBuilder(altRupSet);
//		
////		for (int c=u3Constraints.size(); --c>=0;)
////			if (u3Constraints.get(c) instanceof MFDSubSectNuclInversionConstraint)
////				u3Constraints.remove(c);
////		dirName += "-no-sect-mfd";
//		
//		for (int c=u3Constraints.size(); --c>=0;)
//			if (u3Constraints.get(c) instanceof MFDSubSectNuclInversionConstraint)
//				u3Constraints.remove(c);
//		double nuclWt = 0.001;
//		newBuilder.u2NuclBVals(true, false).weight(nuclWt);
//		dirName += "-sect_nucl_wt_"+oDF.format(nuclWt);
//		
//		for (int c=u3Constraints.size(); --c>=0;)
//			if (u3Constraints.get(c) instanceof MFDLaplacianSmoothingInversionConstraint)
//				u3Constraints.remove(c);
//		double smoothWt = 10000;
//		newBuilder.supraPaleoSmooth().weight(smoothWt);
//		dirName += "-supra_paleo_smooth_wt_"+oDF.format(smoothWt);
//		
//		u3Constraints.addAll(newBuilder.build());
//		
////		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5); dirName += "-5h";
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(2); dirName += "-2h";
////		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(20); dirName += "-20m";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(5);
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints, completion)
//				.threads(remoteToalThreads)
//				.avgThreads(remoteToalThreads/4, avgCompletion)
//				.build();
//		for (int i=0; i<10; i++) {
//			configs.add(config);
//			subDirNames.add("run_"+i);
//		}
//		avgJob = true;
		
		
		/*
		 * uncertainty weighted, 200 runs misfit-targeted
		 */
//		dirName += "-"+rsPrefix+"-misfit_std_dev_targeted";
//		dirName += "-slip_only";
//		DoubleUnaryOperator mfdStdDevFunc = null;
////		DoubleUnaryOperator mfdStdDevFunc = M->0.1; dirName += "-mfd_sd_0.1";
////		DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5)); dirName += "-mfd_sd_0.1xMmin5";
////		DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-mfd_sd_0.1pls10powMmin8";
//		double slipWeight = 1d;
//		double mfdWeight = 0d;
//		double paleoWeight = 0d;
//		double parkfieldWeight = 0d;
//		double minimizeWeight = 10000d;
//		double mfdSmoothWeight = 0d;
//		double supraSmoothWeight = 1000d;
//		double u2NuclWeight = 0d;
//		double parentSmoothWeight = 0;
//		
//		int num = 500;
//
//		if (mfdWeight > 0 && mfdWeight != 1d && mfdStdDevFunc != null)
//			dirName += "_wt"+oDF.format(mfdWeight);
//		if (slipWeight > 0 && slipWeight != 1d)
//			dirName += "-slip"+oDF.format(slipWeight);
//		if (paleoWeight > 0 && paleoWeight != 1d)
//			dirName += "-paleo"+oDF.format(paleoWeight);
//		if (parkfieldWeight > 0 && parkfieldWeight != 1d)
//			dirName += "-prakfield"+oDF.format(parkfieldWeight);
//		if (minimizeWeight > 0)
//			dirName += "-minimize"+oDF.format(minimizeWeight);
//		if (supraSmoothWeight > 0)
//			dirName += "-supra_smooth"+oDF.format(supraSmoothWeight);
//		if (mfdSmoothWeight > 0)
//			dirName += "-mfd_smooth"+oDF.format(mfdSmoothWeight);
//		if (parentSmoothWeight > 0)
//			dirName += "-parent_smooth"+oDF.format(parentSmoothWeight);
//		if (u2NuclWeight > 0)
//			dirName += "-u2Nucl"+oDF.format(u2NuclWeight);
//		List<InversionConstraint> u3Constraints = InversionsCLI.getStdDevWeightedU3Constraints(
//				rupSet, slipWeight, mfdWeight, mfdStdDevFunc, paleoWeight, parkfieldWeight,
//				minimizeWeight, u2NuclWeight, supraSmoothWeight, mfdSmoothWeight);
//		if (parentSmoothWeight > 0d)
//			u3Constraints.add(new ParentSectSmoothnessConstraint(rupSet, parentSmoothWeight, true));
//		
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_10_26-reproduce-ucerf3-ref_branch-uniform-new_anneal-uncert_weighted-only-slip_rates-minimize10000-smooth1000-5h/mean_solution.zip");
//		remoteMeanCompareName = "Only-Slip-Overfit";
//		
////		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5); dirName += "-5h";
////		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(20);
//
////		double[] startingModel = null;
//		
//		dirName += "-start_smooth";
//		double[] startingModel = Inversions.getDefaultVariablePerturbationBasis(rupSet);
//		if (minimizeWeight > 0d) {
//			// have it never sample reates that are being minimized
//			for (InversionConstraint constraint : u3Constraints) {
//				if (constraint instanceof RupRateMinimizationConstraint) {
//					for (int rupIndex : ((RupRateMinimizationConstraint)constraint).getRupIndexes())
//						startingModel[rupIndex] = 0d;
//				}
//			}
//		}
//		
//		double[] samplerRates = null;
//		
////		dirName += "-rup_sampler";
////		double[] samplerRates = Inversions.getDefaultVariablePerturbationBasis(rupSet);
////		if (minimizeWeight > 0d) {
////			// have it never sample reates that are being minimized
////			for (InversionConstraint constraint : u3Constraints) {
////				if (constraint instanceof RupRateMinimizationConstraint) {
////					for (int rupIndex : ((RupRateMinimizationConstraint)constraint).getRupIndexes())
////						samplerRates[rupIndex] = 0d;
////				}
////			}
////		}
//		
//		dirName += "-sd1";
//		CompletionCriteria completion = new MisfitStdDevCompletionCriteria(
//				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, 1d);
//		
//		InversionConfiguration config = InversionConfiguration.builder(u3Constraints, completion)
//				.threads(remoteToalThreads).sampler(samplerRates).initialSolution(startingModel).build();
//		for (int i=0; i<num; i++) {
//			configs.add(config);
//			subDirNames.add("misfit_sd1_run_"+i);
//		}
//		avgJob = true;
//		allPlotLevel = null;
		
		/*
		 * uncertainty weighted with MFD alternatives
		 */
//		double b = 1d;
////		double bWeight = 0;
//		double bWeight = 1000;
//		double rateWeight = 1d;
//		double[] rateMags = { 0d };
//		dirName += "-"+rsPrefix+"-uncert_weighted";
//		
//		double slipWeight = 1d;
//		double paleoWeight = 1d;
//		double parkfieldWeight = 1d;
//		double minimizeWeight = 10000d;
////		double minimizeWeight = 0d;
////		double mfdSmoothWeight = 1000d;
//		double mfdSmoothWeight = 0d;
//		double supraSmoothWeight = 1000d;
////		double supraSmoothWeight = 0d;
////		double u2NuclWeight = 0.01d;
//		double u2NuclWeight = 0d;
//		
//		int num = 5;
//
//		if (bWeight > 0)
//			dirName += "-bVal"+(float)b;
//		if (rateWeight > 0d) {
//			dirName += "-rate";
//			for (int m=0; m<rateMags.length; m++) {
//				if (m > 0)
//					dirName += "_";
//				if (rateMags[m] > 0)
//					dirName += "M"+oDF.format(rateMags[m]);
//				else
//					dirName += "Supra";
//			}
//		}
//		if (minimizeWeight > 0)
//			dirName += "-minimize"+oDF.format(minimizeWeight);
//		if (supraSmoothWeight > 0)
//			dirName += "-supra_smooth"+oDF.format(supraSmoothWeight);
//		if (mfdSmoothWeight > 0)
//			dirName += "-mfd_smooth"+oDF.format(mfdSmoothWeight);
//		if (u2NuclWeight > 0)
//			dirName += "-u2Nucl"+oDF.format(u2NuclWeight);
//		dirName += "-5h";
//		List<InversionConstraint> u3Constraints = InversionsCLI.getStdDevWeightedU3Constraints(
//				rupSet, slipWeight, 0d, null, paleoWeight, parkfieldWeight,
//				minimizeWeight, u2NuclWeight, supraSmoothWeight, mfdSmoothWeight);
//		if (bWeight > 0d)
//			u3Constraints.add(0, new RelativeBValueConstraint(rupSet, b, bWeight));
//		if (rateWeight > 0d) {
//			InversionTargetMFDs mfds = rupSet.requireModule(InversionTargetMFDs.class);
//			EvenlyDiscretizedFunc cmlMFD = mfds.getTotalOnFaultSupraSeisMFD().getCumRateDistWithOffset();
//			for (int m=rateMags.length; --m>=0;) {
//				double target = cmlMFD.getInterpolatedY_inLogYDomain(rateMags[m]);
//				u3Constraints.add(0, new TotalRateInversionConstraint(rateWeight, target, rupSet, rateMags[m],
//						ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, 0.1*target));
//			}
//		}
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints, TimeCompletionCriteria.getInHours(5))
//				.threads(remoteToalThreads)
//				.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
//				.build();
//		for (int i=0; i<num; i++) {
//			configs.add(config);
//			subDirNames.add("uncert_weight_run_"+i);
//		}
//		avgJob = true;
		
		/*
		 * UCERF3 with each combination of constraints
		 */
//		dirName += "-"+rsPrefix+"-u3_constraint_by_constraint";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		Set<Set<InversionConstraint>> powerSet = Sets.powerSet(new HashSet<>(u3Constraints));
//		// for consistent ordering
//		HashMap<String, Integer> nameOrders = new HashMap<>();
//		for (int i=0; i<u3Constraints.size(); i++)
//			nameOrders.put(u3Constraints.get(i).getName(), i);
//		Comparator<InversionConstraint> constrComp = new Comparator<InversionConstraint>() {
//
//			@Override
//			public int compare(InversionConstraint o1, InversionConstraint o2) {
//				int i0 = nameOrders.get(o1.getName());
//				int i1 = nameOrders.get(o2.getName());
//				return Integer.compare(i0, i1);
//			}
//		};
//		for (Set<InversionConstraint> subSet : powerSet) {
//			if (subSet.isEmpty())
//				continue;
//			List<InversionConstraint> myConstraints = new ArrayList<>(subSet);
//			Collections.sort(myConstraints, constrComp);
//			String name = null;
//			for (InversionConstraint constr : myConstraints) {
//				if (name == null)
//					name = "";
//				else
//					name += "-";
//				name += constr.getShortName().replaceAll("\\W+", "");
//			}
//			if (!name.contains("RateMinimize"))
//				// always keep the minimization constraint
//				continue;
//			if (name.contains("SlipRate") && !name.contains("NormSlipRate-SlipRate"))
//				// bundle both slip rate constraints together
//				continue;
//			if ((name.contains("MFDEquality") || name.contains("MFDInequality"))
//					&& !name.contains("MFDEquality-MFDInequality"))
//				// bundle both MFD constraints together
//				continue;
//			subDirNames.add(name);
//			InversionConfiguration config = InversionConfiguration.builder(
//					myConstraints, TimeCompletionCriteria.getInHours(5))
//					.threads(remoteToalThreads)
//					.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
//					.build();
//			configs.add(config);
//		}
//		avgJob = false;
//		primaryPlotLevel = PlotLevel.DEFAULT;
//		allPlotLevel = PlotLevel.DEFAULT;
		
		/*
		 * Constraint removal one at a time, 10x each
		 */
////		Class<? extends InversionConstraint> removeClass = SlipRateInversionConstraint.class;
////		Class<? extends InversionConstraint> removeClass = PaleoRateInversionConstraint.class;
////		Class<? extends InversionConstraint> removeClass = PaleoSlipInversionConstraint.class;
////		Class<? extends InversionConstraint> removeClass = RupRateMinimizationConstraint.class;
////		Class<? extends InversionConstraint> removeClass = MFDInversionConstraint.class;
////		Class<? extends InversionConstraint> removeClass = MFDSubSectNuclInversionConstraint.class;
////		Class<? extends InversionConstraint> removeClass = MFDLaplacianSmoothingInversionConstraint.class;
//		Class<? extends InversionConstraint> removeClass = ParkfieldInversionConstraint.class;
//		dirName += "-"+rsPrefix+"-u3_without";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		List<InversionConstraint> myConstraints = new ArrayList<>();
//		for (InversionConstraint constraint : u3Constraints) {
//			if (removeClass.isAssignableFrom(constraint.getClass()))
//				dirName += "_"+constraint.getShortName().replaceAll("\\W+", "");
//			else
//				myConstraints.add(constraint);
//		}
//		InversionConfiguration config = InversionConfiguration.builder(
//				myConstraints, TimeCompletionCriteria.getInHours(5))
//				.threads(remoteToalThreads)
//				.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
//				.build();
//		for (int i=0; i<10; i++) {
//			subDirNames.add("removal_run_"+i);
//			configs.add(config);
//		}
//		avgJob = true;
//		primaryPlotLevel = PlotLevel.FULL;
//		allPlotLevel = PlotLevel.LIGHT;
		
		Preconditions.checkState(!configs.isEmpty());
		Preconditions.checkState(configs.size() == subDirNames.size());
		
		SlipAlongRuptureModel dsr = rupSet.getModule(SlipAlongRuptureModel.class);
		if (dirName.contains("uniform"))
			Preconditions.checkState(dsr instanceof SlipAlongRuptureModel.Uniform,
					"Directory name (%s) indicates uniform, but Dsr module is %s", dirName, dsr.getName());
		else if (dirName.contains("tapered"))
			Preconditions.checkState(dsr instanceof SlipAlongRuptureModel.Tapered,
					"Directory name (%s) indicates tapered, but Dsr module is %s", dirName, dsr.getName());
		
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		System.out.println("Writing "+configs.size()+" configurations to "+localDir.getName());
		File remoteDir = new File(remoteMainDir, dirName);
		File localRupSet = new File(localDir, "rupture_set.zip");
		rupSet.write(localRupSet);
		File remoteRupSet = new File(remoteDir, localRupSet.getName());
		File remoteJar = new File(remoteDir, "opensha-all.jar");
		
		String java = "java -Djava.awt.headless=true -Xmx"+remoteTotalMemGB+"G -cp "+remoteJar.getAbsolutePath();
		
		List<File> outputSolutions = new ArrayList<>();
		
		boolean allConfigsEqual = true;
		InversionConfiguration config0 = configs.get(0);
		for (int i=1; i<configs.size(); i++) {
			if (config0 != configs.get(i)) {
				allConfigsEqual = false;
				break;
			}
		}
		
		File equalConfigsRemoteFile = null;
		if (allConfigsEqual) {
			System.out.println("All configurations are equal, only writing out once");
			File equalConfigsLocalFile = new File(localDir, "config.json");
			equalConfigsRemoteFile = new File(remoteDir, equalConfigsLocalFile.getName());
			System.out.println("Writing configuration to "+equalConfigsRemoteFile.getAbsolutePath());
			InversionConfiguration.writeJSON(config0, equalConfigsLocalFile);
		}
		
		for (int i=0; i<configs.size(); i++) {
			String name = subDirNames.get(i);
			InversionConfiguration subConfig = configs.get(i);
			
			File localSubDir = new File(localDir, name);
			Preconditions.checkState(localSubDir.exists() || localSubDir.mkdir());
			File remoteSubDir = new File(remoteDir, name);
			
			int mins;
			CompletionCriteria invCompletion = subConfig.getCompletionCriteria();
			if (invCompletion instanceof TimeCompletionCriteria) {
				long invertMillis = ((TimeCompletionCriteria)subConfig.getCompletionCriteria()).getMillis();
				mins = (int)((double)invertMillis/(1000d*60d));
			} else if (invCompletion instanceof IterationCompletionCriteria) {
				long iters = ((IterationCompletionCriteria)invCompletion).getMinIterations();
				long secs = iters/5000; // conservative: 5000 iterations per second
				mins = (int)(secs/60);
				System.out.println("Estimated "+mins+"m inversion runtime");
			} else if (invCompletion instanceof MisfitStdDevCompletionCriteria){
				mins = 20*60;
			} else {
				throw new IllegalStateException("Cannot estimate job runtime for completion criteria: "+invCompletion);
			}
			// add an extra hour buffer
			mins += 60;
			
			// write the configuration
			File remoteConfig;
			if (allConfigsEqual) {
				remoteConfig = equalConfigsRemoteFile;
			} else {
				File localConfig = new File(localSubDir, "config.json");
				System.out.println("Writing configuration to "+localConfig.getAbsolutePath());
				InversionConfiguration.writeJSON(subConfig, localConfig);
				remoteConfig = new File(remoteSubDir, localConfig.getName());
			}
			
			File remoteOutput = new File(remoteSubDir, "solution.zip");
			outputSolutions.add(remoteOutput);
			
			List<String> script = new ArrayList<>();
			
			PlotLevel plotLevel = i == 0 ? primaryPlotLevel : allPlotLevel;
			
			script.add("#!/bin/bash");
			script.add("");
			script.add("# run the inversion");
			script.add("echo \"Running the Inversion\"");
			script.add(java+" "+Inversions.class.getName()
				+" --rupture-set "+remoteRupSet.getAbsolutePath()
				+" --output-file "+remoteOutput.getAbsolutePath()
				+" --config-json "+remoteConfig.getAbsolutePath());
			if (plotLevel != null) {
				script.add("");
				script.add("if [[ -e "+remoteOutput.getAbsolutePath()+" ]];then");
				script.add("    # build a report");
				script.add("    echo \"Building Report\"");
				script.add("    "+java+" "+ReportPageGen.class.getName()
					+" --input-file "+remoteOutput.getAbsolutePath()
					+" --plot-level "+plotLevel.name()
					+" --name "+name
					+" --output-dir "+remoteSubDir.getAbsolutePath());
				script.add("else");
				script.add("    echo \"Inversion failed, see errors above\"");
				script.add("    exit 1");
				script.add("fi");
			}
			script.add("");
			script.add("# open up permissions");
			script.add("chmod -R go+rX "+remoteSubDir.getAbsolutePath());
			
			script = scriptWrite.buildScript(script, mins, 1, remoteToalThreads, queue);
			
			File localScript = new File(localSubDir, name+".slurm");
			System.out.println("Writing "+localScript.getAbsolutePath());
			FileWriter fw = new FileWriter(localScript);
			for (String line : script)
				fw.write(line+"\n");
			fw.close();
		}
		
		if (avgJob) {
			writeMeanJob(remoteToalThreads, scriptWrite, queue, remoteMeanCompFile, remoteMeanCompareName,
					avgPlotLevel, localDir, remoteDir, java, outputSolutions);
		}
	}

	public static void writeMeanJob(int remoteToalThreads, BatchScriptWriter scriptWrite, String queue,
			File remoteMeanCompFile, String remoteMeanCompareName, PlotLevel plotLevel, File localDir,
			File remoteDir, String java, List<File> outputSolutions) throws IOException {
		List<String> script = new ArrayList<>();

		File remoteOutput = new File(remoteDir, "mean_solution.zip");
		File plotDir = new File(remoteDir, "mean_solution");
		
		script.add("#!/bin/bash");
		script.add("");
		script.add("# averaging");
		script.add("echo \"Averaging "+outputSolutions.size()+" solutions\"");
		String argz = java+" "+AverageSolutionCreator.class.getName()+" "+remoteOutput.getAbsolutePath();
		for (File output : outputSolutions)
			argz += " "+output.getAbsolutePath();
		script.add(argz);
		script.add("");
		script.add("if [[ -e "+remoteOutput.getAbsolutePath()+" ]];then");
		script.add("    # build a report");
		script.add("    echo \"Building Report\"");
		script.add("    export FST_HAZARD_SPACING=0.2");
		String reportCommand = "    "+java+" "+ReportPageGen.class.getName()
			+ " --input-file "+remoteOutput.getAbsolutePath()
			+" --plot-level "+plotLevel.name()
			+" --output-dir "+plotDir.getAbsolutePath();
		if (remoteMeanCompFile != null) {
			reportCommand += " --compare-to "+remoteMeanCompFile.getAbsolutePath();
			if (remoteMeanCompareName != null)
				reportCommand += " --comp-name \""+remoteMeanCompareName+"\"";
		}
		script.add(reportCommand);
		script.add("else");
		script.add("    echo \"Averaging failed, see errors above\"");
		script.add("    exit 1");
		script.add("fi");
		script.add("");
		script.add("# open up permissions");
		script.add("chmod -R go+rX "+remoteDir.getAbsolutePath());
		
		script = scriptWrite.buildScript(script, 60, 1, remoteToalThreads, queue);
		
		File localScript = new File(localDir, "mean_solution.slurm");
		System.out.println("Writing "+localScript.getAbsolutePath());
		FileWriter fw = new FileWriter(localScript);
		for (String line : script)
			fw.write(line+"\n");
		fw.close();
	}
	
	private static void removeConstraintsByType(List<InversionConstraint> constraints, Class<? extends InversionConstraint> type) {
		int num = 0;
		for (int i=constraints.size(); --i>=0;) {
			if (type.isAssignableFrom(constraints.get(i).getClass())) {
				constraints.remove(i);
				num++;
			}
		}
		Preconditions.checkState(num > 0);
	}

}
