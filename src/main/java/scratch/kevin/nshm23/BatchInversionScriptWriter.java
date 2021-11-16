package scratch.kevin.nshm23;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.APrioriInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDLaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParentSectSmoothnessConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RelativeBValueConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.TotalRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.MisfitStdDevCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.A_PrioriRupRates;
import org.opensha.sha.magdist.gui.MagFreqDistAppWindow;

import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;

import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.nshm23.targetMFDs.DraftModelConstraintBuilder;
import scratch.nshm23.targetMFDs.InversionTargetMFDsFromBValAndDefModel;

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
//		dirName += "-"+rsPrefix+"-new_anneal-5x_avg-try_zero-var_perturb-noWL-5h";
////		dirName += "-"+rsPrefix+"-new_anneal-no_avg-try_zero-var_perturb-noWL-5h";
////		dirName += "-"+rsPrefix+"-new_anneal-no_avg-limit_zero-var_perturb-noWL-5h";
//		UCERF3InversionInputGenerator u3Gen = InversionsCLI.getU3Generator(rupSet);
//		List<InversionConstraint> u3Constraints = u3Gen.getConstraints();
//		InversionConfiguration config = InversionConfiguration.builder(
//				u3Constraints, TimeCompletionCriteria.getInHours(5))
//				.threads(remoteToalThreads)
//				.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
////				.nonNegativity(NonnegativityConstraintType.LIMIT_ZERO_RATES)
//				.build();
//		for (int i=0; i<10; i++) {
//			configs.add(config);
//			subDirNames.add("new_anneal_run_"+i);
//		}
//		avgJob = true;
		
		/*
		 * uncertainty weighted
		 */
//		dirName += "-"+rsPrefix+"-new_anneal-uncert_weighted";
////		DoubleUnaryOperator mfdStdDevFunc = null;
//		DoubleUnaryOperator mfdStdDevFunc = M->0.1; dirName += "-mfd_sd_0.1";
////		DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5)); dirName += "-mfd_sd_0.1xMmin5";
////		DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-mfd_sd_0.1pls10powMmin8";
//		double slipWeight = 1d;
//		double mfdWeight = 5d;
//		double paleoWeight = 5d;
//		double parkfieldWeight = 5d;
//		double minimizeWeight = 10000d;
////		double minimizeWeight = 0d;
////		double mfdSmoothWeight = 1000d;
//		double mfdSmoothWeight = 0d;
//		double supraSmoothWeight = 1000d;
////		double supraSmoothWeight = 0d;
////		double u2NuclWeight = 0.01d;
//		double u2NuclWeight = 0d;
//		double parentSmoothWeight = 0;
//		
//		int num = 10;
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
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5); dirName += "-5h";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(20);
//		
//		InversionConfiguration.Builder builder = InversionConfiguration.builder(u3Constraints, completion)
//				.threads(remoteToalThreads).avgThreads(remoteToalThreads/4, avgCompletion);
//		
//		InversionConfiguration config = builder.build();
//		for (int i=0; i<num; i++) {
//			configs.add(config);
//			subDirNames.add("uncert_weight_run_"+i);
//		}
//		avgJob = true;
//		allPlotLevel = null;
		
		/*
		 * new NSHM23 draft scheme
		 */
		dirName += "-"+rsPrefix+"-nshm23_draft";
		double bVal = 0.8;
		dirName += "-supra_b_"+(float)bVal;

//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_03-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-2h/mean_solution.zip");
//		remoteMeanCompareName = "All-New-Constr-b=0.8";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_08-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-only-slip-mfd-mfd_wt_10-skipBelow-2h/mean_solution.zip");
//		remoteMeanCompareName = "New-MFD-Constr-b=0.8-No-Nucl";
		remoteMeanCompFile = new File(remoteMainDir,
				"2021_11_15-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_100-sect_wt_0.5-smooth_paleo_wt_10000-skipBelow-2h/mean_solution.zip");
		remoteMeanCompareName = "Weighted-Draft-b=0.8";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_10-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-sect_wt_0.5-smooth_paleo_wt_10000-skipBelow-10h-20x/mean_solution.zip");
//		remoteMeanCompareName = "New-Draft-0.8-10hr-20x";
//		remoteMeanCompFile = new File(remoteMainDir,
//		"2021_11_10-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-no_sect_rate-smooth_paleo_wt_10000-skipBelow-10h-20x/mean_solution.zip");
//		remoteMeanCompareName = "New-Draft-NoSect-0.8-10hr-20x";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_12-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-no_sect_rate-smooth_all_wt_10000-skipBelow-10h/mean_solution.zip");
//		remoteMeanCompareName = "New-Draft-NoSect-0.8-SmoothAll-10hr";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_13-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-no_sect_rate-parent_smooth_wt_10000-smooth_paleo_wt_10000-skipBelow-10h/mean_solution.zip");
//		remoteMeanCompareName = "ParentSmooth-10hr";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_12-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-no_sect_rate-smooth_paleo_wt_10000-skipBelow-10h-simple_exp_perturb/mean_solution.zip");
//		remoteMeanCompareName = "10hr-NoSect-SimplePerturb";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_13-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_10-paleo_wt_5-parkfield_wt_10-no_sect_rate-smooth_paleo_wt_10000-parentFairSampler-10h/mean_solution.zip");
//		remoteMeanCompareName = "10hr-ParentFairSampler";
//		remoteMeanCompFile = new File(remoteMainDir,
//				"2021_11_14-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_0.8-mfd_wt_1-paleo_wt_5-parkfield_wt_10-no_sect_rate-smooth_paleo_wt_10000-skipBelow-10h/mean_solution.zip");
//		remoteMeanCompareName = "10hr-MFDWt1";
		
		int num = 5;

		DraftModelConstraintBuilder constrBuilder = new DraftModelConstraintBuilder(rupSet);
		
		constrBuilder.defaultConstraints(bVal);
		
//		dirName += "-no_mfd";
//		constrBuilder.except(MFDInversionConstraint.class);
		
//		dirName += "-only-slip-mfd-sect_rate";
//		constrBuilder.slipRates().supraBValMFDs(bVal).sectSupraRates(bVal).defaultMetaConstraints();
		
//		dirName += "-only-slip-mfd";
//		constrBuilder.slipRates().supraBValMFDs(bVal).defaultMetaConstraints();
		
//		dirName += "-only-slip-sect_rate";
//		constrBuilder.slipRates().sectSupraRates(bVal).defaultMetaConstraints();
		
		double mfdWeight = 10;
		dirName += "-mfd_wt_"+oDF.format(mfdWeight);
		constrBuilder.weight(MFDInversionConstraint.class, mfdWeight);
		
		double paleoWeight = 5;
		dirName += "-paleo_wt_"+oDF.format(paleoWeight);
		constrBuilder.weight(PaleoRateInversionConstraint.class, paleoWeight);
		constrBuilder.weight(PaleoSlipInversionConstraint.class, paleoWeight);
		
//		dirName += "-no_paleo";
//		constrBuilder.except(PaleoRateInversionConstraint.class).except(PaleoSlipInversionConstraint.class);
		
		double parkWeight = 100;
		dirName += "-parkfield_wt_"+oDF.format(parkWeight);
		constrBuilder.weight(ParkfieldInversionConstraint.class, parkWeight);
		
//		dirName += "-no_parkfield";
//		constrBuilder.except(ParkfieldInversionConstraint.class);
		
		double nuclWeight = 0.5;
		dirName += "-sect_wt_"+oDF.format(nuclWeight);
		constrBuilder.weight(SectionTotalRateConstraint.class, nuclWeight);
		
//		dirName += "-no_sect_rate";
//		constrBuilder.except(SectionTotalRateConstraint.class);
		
//		dirName += "_u2";
//		boolean aFaults = true;
////		boolean aFaults = false; dirName += "_all_as_bflts"
//		boolean reproduce = false;
//		constrBuilder.except(SectionTotalRateConstraint.class);
//		constrBuilder.u2NuclBVals(aFaults, reproduce).weight(nuclWeight);
		
//		FaultSystemSolution refSol = FaultSystemSolution.load(
//				new File(localMainDir, "2021_10_18-reproduce-ucerf3-ref_branch-uniform-new_anneal-5x_avg"
//						+ "-try_zero-var_perturb-noWL-5h/mean_solution.zip"));
//		constrBuilder.except(SectionTotalRateConstraint.class);
//		constrBuilder.testSameBVals(refSol); dirName += "_u3_same";
//		constrBuilder.weight(nuclWeight);
		
//		dirName += "-smooth_all";
//		constrBuilder.supraSmooth();
//		constrBuilder.weight(10000); dirName += "_wt_10000";
		
//		int parentWt = 10000;
//		dirName += "-parent_smooth_wt_"+parentWt;
//		constrBuilder.add(new ParentSectSmoothnessConstraint(rupSet, parentWt, true));
		
		dirName += "-smooth_paleo";
		constrBuilder.supraPaleoSmooth();
		constrBuilder.weight(10000); dirName += "_wt_10000";
		
//		FaultSystemSolution refSol = FaultSystemSolution.load(
//				new File(new File(localMainDir, remoteMeanCompFile.getParentFile().getName()),
//						remoteMeanCompFile.getName()));
//		constrBuilder.testFlipBVals(refSol, bVal); dirName += "-test_flip_b";
//		constrBuilder.testSameBVals(refSol); dirName += "-test_same_b";
//		constrBuilder.weight(10d); dirName += "_wt_10";
//		IntegerPDF_FunctionSampler sampler = constrBuilder.testGetSampleAllNew(refSol, true); dirName += "-test_sample_new";
//		constrBuilder.except(RupRateMinimizationConstraint.class);
		
		IntegerPDF_FunctionSampler sampler = constrBuilder.getSkipBelowMinSampler();
		dirName += "-skipBelow";
		constrBuilder.except(RupRateMinimizationConstraint.class);
		
//		IntegerPDF_FunctionSampler sampler = constrBuilder.testSampleFaultsEqually(true, 100d);
//		dirName += "-parentFairSampler";
//		constrBuilder.except(RupRateMinimizationConstraint.class);
		
//		dirName += "-sampleMohawkSurprise10x";
//		for (int r : rupSet.getRupturesForParentSection(686)) // mohawk valley
//			sampler.set(r, sampler.getY(r)*10);
//		for (int r : rupSet.getRupturesForParentSection(708)) // surprise valley
//			sampler.set(r, sampler.getY(r)*10);
		
//		dirName += "-redo";
		
//		IntegerPDF_FunctionSampler sampler = null;
		
		List<InversionConstraint> constraints = constrBuilder.build();
		
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

}
