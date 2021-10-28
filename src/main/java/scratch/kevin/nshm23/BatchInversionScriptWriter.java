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

import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.APrioriInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDLaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldUncertaintyWeightedInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint.WeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.A_PrioriRupRates;

import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;

import scratch.UCERF3.inversion.InversionTargetMFDs;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;

public class BatchInversionScriptWriter {
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteToalThreads = 20;
		int remoteTotalMemGB = 55;
		BatchScriptWriter scriptWrite = new USC_CARC_ScriptWriter();
		String queue = "scec";
		
		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
//		String dirName = "2021_10_12";
		
		List<InversionConfiguration> configs = new ArrayList<>();
		List<String> subDirNames = new ArrayList<>();
		
		File rsDir = new File("/home/kevin/markdown/inversions/");
		
		File rupSetFile;
		String rsPrefix;
		
		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_reproduce_ucerf3.zip");
		rsPrefix = "reproduce-ucerf3-ref_branch-uniform";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_reproduce_ucerf3_fractGrow0.1.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-uniform-grow0.1";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_tapered_reproduce_ucerf3.zip");
//		rsPrefix = "reproduce-ucerf3-ref_branch-tapered";
		
//		rupSetFile = new File(rsDir, "fm3_1_u3ref_uniform_coulomb.zip");
//		rsPrefix = "coulomb-fm31-ref_branch-uniform";
		
//		File remoteMeanCompFile = new File(remoteMainDir,
//				"2021_10_18-reproduce-ucerf3-ref_branch-uniform-new_anneal-5x_avg-try_zero-var_perturb-noWL-5h/mean_solution.zip");
//		String remoteMeanCompareName = "U3-New-Anneal";
		File remoteMeanCompFile = new File(remoteMainDir,
				"2021_10_25-reproduce-ucerf3-ref_branch-uniform-new_anneal-uncert_weighted-mfd_sd_0.1-minimize10000-smooth1000-5h/mean_solution.zip");
		String remoteMeanCompareName = "U3-Uncert-Wtd";
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(rupSetFile);
		
		// this is applied to the first configuration supplied
//		PlotLevel primaryPlotLevel = PlotLevel.FULL;
		PlotLevel primaryPlotLevel = PlotLevel.DEFAULT;
//		PlotLevel primaryPlotLevel = PlotLevel.LIGHT;
		// this is applied to everything else
//		PlotLevel allPlotLevel = PlotLevel.FULL;
//		PlotLevel allPlotLevel = PlotLevel.DEFAULT;
		PlotLevel allPlotLevel = PlotLevel.LIGHT;
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
		dirName += "-"+rsPrefix+"-new_anneal-uncert_weighted";
		DoubleUnaryOperator mfdStdDevFunc = M->0.1; dirName += "-mfd_sd_0.1";
//		DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5)); dirName += "-mfd_sd_0.1xMmin5";
//		DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-mfd_sd_0.1pls10powMmin8";
		double minimizeWeight = 10000d;
//		double minimizeWeight = 0d;
//		double mfdSmoothWeight = 1000d;
		double mfdSmoothWeight = 0d;
		double supraSmoothWeight = 1000d;
//		double supraSmoothWeight = 0d;
//		double u2NuclWeight = 0.01d;
		double u2NuclWeight = 0d;
		
		if (minimizeWeight > 0)
			dirName += "-minimize"+oDF.format(minimizeWeight);
		if (supraSmoothWeight > 0)
			dirName += "-supra_smooth"+oDF.format(supraSmoothWeight);
		if (mfdSmoothWeight > 0)
			dirName += "-mfd_smooth"+oDF.format(mfdSmoothWeight);
		if (u2NuclWeight > 0)
			dirName += "-u2Nucl"+oDF.format(u2NuclWeight);
		dirName += "-5h";
		List<InversionConstraint> u3Constraints = InversionsCLI.getStdDevWeightedU3Constraints(
				rupSet, mfdStdDevFunc, minimizeWeight, u2NuclWeight, supraSmoothWeight, mfdSmoothWeight);
		InversionConfiguration config = InversionConfiguration.builder(
				u3Constraints, TimeCompletionCriteria.getInHours(5))
				.threads(remoteToalThreads)
				.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
				.build();
		for (int i=0; i<10; i++) {
			configs.add(config);
			subDirNames.add("uncert_weight_run_"+i);
		}
		avgJob = true;
		
		/*
		 * uncertainty weighted, individual
		 */
//		dirName += "-"+rsPrefix+"-new_anneal-uncert_weighted-only";
//		
//		boolean slipRates = true; dirName += "-slip_rates";
////		boolean slipRates = false;
//		
////		boolean parkfield = true; dirName += "-parkfield";
//		boolean parkfield = false;
//		
////		boolean paleoRate = true; dirName += "-paleoRate";
//		boolean paleoRate = false;
//		
////		boolean paleoSlip = true; dirName += "-paleoSlip";
//		boolean paleoSlip = false;
//		
//		DoubleUnaryOperator mfdStdDevFunc = null;
////		DoubleUnaryOperator mfdStdDevFunc = M->0.1; dirName += "-mfd_sd_0.1";
////		DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5)); dirName += "-mfd_sd_0.1xMmin5";
////		DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-mfd_sd_0.1pls10powMmin8";'
//		
//		double minimizeWeight = 10000d;
////		double minimizeWeight = 0d;
////		double mfdSmoothWeight = 1000d;
//		double mfdSmoothWeight = 0d;
//		double supraSmoothWeight = 1000d;
////		double supraSmoothWeight = 0d;
////		double u2NuclWeight = 0.01d;
//		double u2NuclWeight = 0d;
//		
//		List<InversionConstraint> constraints = InversionsCLI.getStdDevWeightedU3Constraints(
//				rupSet, mfdStdDevFunc, minimizeWeight, u2NuclWeight, supraSmoothWeight, mfdSmoothWeight);
//		for (int i=constraints.size(); --i>=0;) {
//			InversionConstraint constr = constraints.get(i);
//			if (!slipRates && constr instanceof SlipRateInversionConstraint)
//				constraints.remove(i);
//			else if (!parkfield && constr instanceof ParkfieldUncertaintyWeightedInversionConstraint)
//				constraints.remove(i);
//			else if (!paleoRate && constr instanceof PaleoRateInversionConstraint)
//				constraints.remove(i);
//			else if (!paleoSlip && constr instanceof PaleoSlipInversionConstraint)
//				constraints.remove(i);
//			else
//				System.out.println("Keeping constraint: "+constr.getName());
//		}
//		
//		if (minimizeWeight > 0)
//			dirName += "-minimize"+oDF.format(minimizeWeight);
//		if (supraSmoothWeight > 0)
//			dirName += "-supra_smooth"+oDF.format(supraSmoothWeight);
//		if (mfdSmoothWeight > 0)
//			dirName += "-mfd_smooth"+oDF.format(mfdSmoothWeight);
//		if (u2NuclWeight > 0)
//			dirName += "-u2Nucl"+oDF.format(u2NuclWeight);
//		dirName += "-5h";
//		
//		InversionConfiguration config = InversionConfiguration.builder(
//				constraints, TimeCompletionCriteria.getInHours(5))
//				.threads(remoteToalThreads)
//				.avgThreads(remoteToalThreads/4, TimeCompletionCriteria.getInMinutes(20))
//				.build();
//		for (int i=0; i<10; i++) {
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
		
		for (int i=0; i<configs.size(); i++) {
			String name = subDirNames.get(i);
			InversionConfiguration subConfig = configs.get(i);
			
			File localSubDir = new File(localDir, name);
			Preconditions.checkState(localSubDir.exists() || localSubDir.mkdir());
			File remoteSubDir = new File(remoteDir, name);
			
			int mins;
			CompletionCriteria completion = subConfig.getCompletionCriteria();
			if (completion instanceof TimeCompletionCriteria) {
				long invertMillis = ((TimeCompletionCriteria)subConfig.getCompletionCriteria()).getMillis();
				mins = (int)((double)invertMillis/(1000d*60d));
			} else if (completion instanceof IterationCompletionCriteria) {
				long iters = ((IterationCompletionCriteria)completion).getMinIterations();
				long secs = iters/800; // conservative: 800 iterations per second
				mins = (int)(secs/60);
				System.out.println("Estimated "+mins+"m inversion runtime");
			} else {
				throw new IllegalStateException("Cannot estimate job runtime for completion criteria: "+completion);
			}
			// add an extra hour buffer
			mins += 60;
			
			// write the configuration
			File localConfig = new File(localSubDir, "config.json");
			System.out.println("Writing configuration to "+localConfig.getAbsolutePath());
			InversionConfiguration.writeJSON(subConfig, localConfig);
			File remoteConfig = new File(remoteSubDir, localConfig.getName());
			
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
