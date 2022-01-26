package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.function.DoubleUnaryOperator;

import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.JumpProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.JumpProbabilityConstraint.InitialModelParticipationRateEstimator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.LaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDLaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.U3MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSystemRupSetCalc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.U3InversionConfigFactory;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.U3SectionMFD_constraint;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class InversionsCLI {

	private static final DecimalFormat oDF = new DecimalFormat("0.##");

	public static void main(String[] args) throws InterruptedException, IOException {
		Date date = new Date();
		
//		System.out.println("Yawn...");
//		long minute = 1000l*60l;
//		long hour = minute*60l;
//		Thread.sleep(4l*hour + 20l*minute);
//		System.out.println("Im awake! "+new Date());
		
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		// run this if I need to attach UCERF3 modules after rebuilding default rupture sets
//		reprocessDefaultRupSets(parentDir);
//		System.exit(0);

		List<String> argz = new ArrayList<>();

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(date);

		dirName += "-coulomb-u3";
		File origRupSetFile = new File(parentDir, "fm3_1_u3ref_uniform_coulomb.zip");

//		dirName += "-u3rs";
//		File origRupSetFile = new File(parentDir, "fm3_1_u3ref_uniform_reproduce_ucerf3.zip");
		
		File rupSetFile = origRupSetFile;
		
//		dirName += "-coulomb-nshm23";
//		File origRupSetFile = new File(parentDir, "nshm23_geo_dm_coulomb.zip");
		
		File compSol = null;

		argz.add("--threads");
		argz.add("16");
		argz.add("--avg-threads");
		argz.add("4");

		dirName += "-slip_constr";
		argz.add("--slip-constraint");
////		argz.add("--slip-weight"); argz.add("1");
////		argz.add("--norm-slip-weight"); argz.add("0.01");
		argz.add("--uncertain-slip-weight");
		argz.add("1");
		dirName += "_uncertain";

//		double b = 0.8;
//		dirName += "-rel_gr_b"+oDF.format(b);
//		argz.add("--rel-gr-constraint");
//		argz.add("--b-value"); argz.add(oDF.format(b));
////		argz.add("--mfd-weight"); argz.add("100"); dirName += "_wt100";
////		argz.add("--mfd-weight"); argz.add("1000"); dirName += "_wt1000";
////		argz.add("--mfd-ineq"); dirName += "_ineq";

//		dirName += "-u3_mfd";
//		argz.add("--mfd-constraint");
//		argz.add("--mfd-transition-mag"); argz.add("7.8"); dirName += "_trans7.8";

		double b = 0.5;
		dirName += "-infer_gr_b"+oDF.format(b);
		argz.add("--mfd-constraint");
		argz.add("--infer-target-gr");
		argz.add("--b-value"); argz.add(oDF.format(b));
		argz.add("--mfd-ineq"); dirName += "_ineq";
//		argz.add("--mfd-transition-mag"); argz.add("7.8"); dirName += "_trans7.8";

//		dirName += "-smooth";
//		argz.add("--smooth");
//		argz.add("--smooth-weight");
//		argz.add("1000");

//		dirName += "-minimize_below";
//		argz.add("--minimize-below-sect-min");
//		argz.add("--minimize-weight");
//		argz.add("10000");
		
		boolean u3Constraints = false;
//		boolean u3Constraints = true;
		boolean u3StdDevConstraints = false;
//		boolean u3StdDevConstraints = true;
//		boolean nshmDraftConstraints = true;
		boolean nshmDraftConstraints = false;
		
		FaultSystemRupSet rupSet = null;

		List<InversionConstraint> extraConstraints = new ArrayList<>();
		
		if (u3Constraints) {
			rupSet = FaultSystemRupSet.load(origRupSetFile);
			extraConstraints.addAll(getU3Constraints(rupSet));
			dirName += "-u3_constraints";
		}
		
		if (u3StdDevConstraints) {
			rupSet = FaultSystemRupSet.load(origRupSetFile);
			DoubleUnaryOperator mfdStdDevFunc = M->0.1;
//			DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5));
//			DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-aggresivePowMSD";
			double slipWeight = 1d;
			double mfdWeight = 5d;
//			double paleoWeight = 1d;
			double paleoWeight = 5d;
			double parkfieldWeight = 10d;
			extraConstraints.addAll(getStdDevWeightedU3Constraints(rupSet, slipWeight, mfdWeight, mfdStdDevFunc,
					paleoWeight, parkfieldWeight, 0, 0, 0, 0d));
			dirName += "-u3_std_dev_tests";
			if (paleoWeight != 1d && paleoWeight > 0)
				dirName += "-paleo"+oDF.format(paleoWeight);
			if (parkfieldWeight != 1d && paleoWeight > 0)
				dirName += "-parkfield"+oDF.format(parkfieldWeight);
		}
		
		if (nshmDraftConstraints) {
			double supraBVal = 0.8;
			rupSet = FaultSystemRupSet.load(origRupSetFile);
			rupSet = new U3InversionConfigFactory().updateRuptureSetForBranch(rupSet,
					rupSet.requireModule(U3LogicTreeBranch.class));
			dirName += "-nshm23_draft-supra_b_"+oDF.format(supraBVal);
			
			boolean applyDefModelUncertaintiesToNucl = true;
			boolean addSectCountUncertaintiesToMFD = false;
			boolean adjustForIncompatibleData = true;
			boolean randWeight = true;

			NSHM23_ConstraintBuilder constrBuilder = new NSHM23_ConstraintBuilder(rupSet, supraBVal,
					applyDefModelUncertaintiesToNucl, addSectCountUncertaintiesToMFD, adjustForIncompatibleData);
			constrBuilder.defaultConstraints();
			
			if (randWeight) {
				dirName += "-randWeightChange";
				constrBuilder.weight(SlipRateInversionConstraint.class, 1*(Math.random()+0.5));
				constrBuilder.weight(MFDInversionConstraint.class, 10*(Math.random()+0.5));
				constrBuilder.weight(PaleoRateInversionConstraint.class, 5*(Math.random()+0.5));
				constrBuilder.weight(PaleoSlipInversionConstraint.class, 5*(Math.random()+0.5));
				constrBuilder.weight(ParkfieldInversionConstraint.class, 100*(Math.random()+0.5));
				constrBuilder.weight(ParkfieldInversionConstraint.class, 100*(Math.random()+0.5));
				constrBuilder.weight(SectionTotalRateConstraint.class, 0.5*(Math.random()+0.5));
			}
			
//			double mfdWeight = 10;
//			dirName += "-mfd_wt_"+oDF.format(mfdWeight);
//			constrBuilder.weight(MFDInversionConstraint.class, mfdWeight);
//			
//			double paleoWeight = 5;
//			dirName += "-paleo_wt_"+oDF.format(paleoWeight);
//			constrBuilder.weight(PaleoRateInversionConstraint.class, paleoWeight);
//			constrBuilder.weight(PaleoSlipInversionConstraint.class, paleoWeight);
//			
////			dirName += "-no_paleo";
////			constrBuilder.except(PaleoRateInversionConstraint.class).except(PaleoSlipInversionConstraint.class);
//			
//			double parkWeight = 10;
//			dirName += "-parkfield_wt_"+oDF.format(parkWeight);
//			constrBuilder.weight(ParkfieldInversionConstraint.class, parkWeight);
//			
////			dirName += "-no_parkfield";
////			constrBuilder.except(ParkfieldInversionConstraint.class);
//			
//			double nuclWeight = 0.5;
//			dirName += "-sect_wt_"+oDF.format(nuclWeight);
//			constrBuilder.weight(SectionTotalRateConstraint.class, nuclWeight);
			
			dirName += "-no_sect_rate";
			constrBuilder.except(SectionTotalRateConstraint.class);
			
			extraConstraints.addAll(constrBuilder.build());
			
			// write out new ruptures set with the new target MFDs
			rupSetFile = File.createTempFile("fst_cli_temp", "rup_set.zip");
			rupSetFile.deleteOnExit();
			rupSet.write(rupSetFile);
		}
		
//		if (rupSet == null)
//			rupSet = FaultSystemRupSet.load(rupSetFile);
//		double r0 = 6d;
//		extraConstraints.add(new RelativeRupJumpDistConstraint(rupSet, r0, 1d, 1000d, false));
//		dirName += "-rup_jump_r0_"+oDF.format(r0);
		
//		if (rupSet == null)
//			rupSet = FaultSystemRupSet.load(rupSetFile);
//		double r0 = 3d;
//		double weight = 1d;
//		JumpProbabilityCalc jumpProbCalc = new Shaw07JumpDistProb(1d, r0);
//		
//		boolean ineq = true;
//		
////		extraConstraints.add(new JumpProbabilityConstraint.ProxySlip(weight, ineq, rupSet, jumpProbCalc));
////		dirName += "-rel_slip_jump_r0_"+oDF.format(r0);
//		
//		extraConstraints.add(new JumpProbabilityConstraint.RelativeRate(weight, ineq, rupSet, jumpProbCalc,
//				new InitialModelParticipationRateEstimator(rupSet, Inversions.getDefaultVariablePerturbationBasis(rupSet))));
//		dirName += "-rel_rate_jump_r0_"+oDF.format(r0);
//		
//		if (weight != 1d)
//			dirName += "_wt"+oDF.format(weight);
//		if (ineq)
//			dirName += "_ineq";
//		
//		compSol = new File(new File(parentDir, "2021_12_10-coulomb-u3-nshm23_draft-supra_b_0.8-10m"), "solution.zip");
		
//		double b = 0.8;
//		dirName += "rel_gr_b"+oDF.format(b);
//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(origRupSetFile);
////		DoubleUnaryOperator relMagStdDev = (M)->0.1; dirName += "_sd0.1";
////		DoubleUnaryOperator relMagStdDev = (M)->0.1+0.5*Math.pow(Math.abs(M-6.5), 2); dirName += "_sd0.1_xCustom";
////		DoubleUnaryOperator relMagStdDev = (M)->0.1+0.25*Math.pow(Math.abs(M-6.5), 2); dirName += "_sd0.1_xCustom2";
//		DoubleUnaryOperator relMagStdDev = (M)->0.1+0.5*Math.pow(Math.abs(M-7), 2); dirName += "_sd0.1_xCustom3";
//		extraConstraints.add(new RelativeBValueConstraint(rupSet, b, 1, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, relMagStdDev));

//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(origRupSetFile);
//		boolean slipAdj = true;
//		double pWeight = 500;
//		extraConstraints.add(new ParentSectSmoothnessConstraint(rupSet, pWeight, slipAdj));
//		dirName += "_pSmooth";
//		if (slipAdj)
//			dirName += "SlipAdj";
//		dirName += oDF.format(pWeight);
		
////		double[] minMags = { 0d };
////		double[] minMags = { 6.5d };
//		double[] minMags = { 7d };
////		double[] minMags = { 0d, 6.5d };
////		double[] minMags = { 0d, 6.5d, 7.5d };
////		double[] minMags = { 0d, 6.5d, 7d, 7.5d };
//		ConstraintWeightingType rateWeightType = ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY;
//		double rateRelStdDev = 0.1;
//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(origRupSetFile);
//		EvenlyDiscretizedFunc supraCmlMFD = rupSet.requireModule(InversionTargetMFDs.class)
//				.getTotalOnFaultSupraSeisMFD().getCumRateDistWithOffset();
//		dirName += rateWeightType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY ? "-uncertTotRate" : "-totRate";
//		for (int i=0; i<minMags.length; i++) {
//			double minMag = minMags[i];
//			double rate = supraCmlMFD.getInterpolatedY_inLogYDomain(minMag);
//			System.out.println("Target rate of M"+(float)minMag+" events: "+rate);
//			extraConstraints.add(new TotalRateInversionConstraint(100d, rate, rupSet, minMag, rateWeightType, rate*rateRelStdDev));
//			if (i > 0)
//				dirName += "_";
//			if (minMag > 0)
//				dirName += "M"+oDF.format(minMag);
//			else
//				dirName += "Supra";
//		}

//		dirName += "-5h";
//		argz.add("--completion"); argz.add("5h");
//		argz.add("--avg-completion"); argz.add("5m");
//		dirName += "-2h";
//		argz.add("--completion"); argz.add("2h");
//		argz.add("--avg-completion"); argz.add("5m");
//		dirName += "-1h";
//		argz.add("--completion"); argz.add("1h");
//		argz.add("--avg-completion"); argz.add("5m");
//		dirName += "-30m";
//		argz.add("--completion"); argz.add("30m");
//		argz.add("--avg-completion"); argz.add("1m");
		dirName += "-10m";
		argz.add("--completion"); argz.add("10m");
		argz.add("--avg-completion"); argz.add("1m");
//		dirName += "-sd1";
//		argz.add("--completion-sd"); argz.add("1");
//		argz.add("--completion-sd-type"); argz.add(ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY.name());
//		argz.add("--avg-completion"); argz.add("1m");

		argz.add("--rupture-set");
		argz.add(rupSetFile.getAbsolutePath());

		File dir = new File(parentDir, dirName);
		System.out.println("Output directory: " + dir.getAbsolutePath());
		if (dir.exists())
			System.err.println("WARNING: directory already exists");
		else
			Preconditions.checkState(dir.mkdir());

		argz.add("--output-file");
		argz.add(new File(dir, "solution.zip").getAbsolutePath());
		
		argz.add("--write-config"); argz.add(new File(dir, "config.json").getAbsolutePath());

		PlotLevel plots = PlotLevel.FULL;

		try {
			FaultSystemSolution sol = Inversions.run(argz.toArray(new String[0]), extraConstraints);

			if (plots != null) {
				// plot it
				String name = dirName.substring(11); // skip date
				File reportDir = new File(dir, "sol_report");
				RupSetMetadata compMeta = null;
				if (compSol != null) {
					FaultSystemSolution cSol = FaultSystemSolution.load(compSol);
					compMeta = new RupSetMetadata(compSol.getParentFile().getName().substring(11), cSol);
				}
				ReportMetadata meta = new ReportMetadata(new RupSetMetadata(name, sol), compMeta);
				ReportPageGen report = new ReportPageGen(meta, reportDir, ReportPageGen.getDefaultSolutionPlots(plots));
				report.setReplot(true);
				report.generatePage();
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static void reprocessDefaultRupSets(File parentDir) throws IOException {
		File[] files = {
				new File(parentDir, "fm3_1_u3ref_uniform_coulomb.zip"),
				new File(parentDir, "fm3_1_u3ref_uniform_reproduce_ucerf3.zip"),
				new File(parentDir, "fm3_1_u3ref_geol_uniform_reproduce_ucerf3.zip"),
				new File(parentDir, "fm3_1_u3ref_uniform_reproduce_ucerf3_fractGrow0.1.zip"),
				new File(parentDir, "fm3_1_u3ref_tapered_reproduce_ucerf3.zip"),
				new File(parentDir, "nshm23_geo_dm_coulomb.zip") };

		U3LogicTreeBranch branch = U3LogicTreeBranch.fromValues(SlipAlongRuptureModels.UNIFORM);
		DeformationModels defaultDM = branch.getValue(DeformationModels.class);
		System.out.println("Default logic tree branch: "+defaultDM);

		for (File file : files) {
			System.out.println("Attaching default UCERF3 modules to " + file.getAbsolutePath());
			FaultSystemRupSet rupSet = FaultSystemRupSet.load(file);
			U3LogicTreeBranch myBranch = branch.copy();
			if (file.getName().startsWith("fm3_1")) {
				myBranch.setValue(FaultModels.FM3_1);
				if (file.getName().contains("_geol_"))
					myBranch.setValue(DeformationModels.GEOLOGIC);
				else
					myBranch.setValue(defaultDM);
			} else {
				myBranch.clearValue(FaultModels.class);
				myBranch.clearValue(DeformationModels.class);
			}
			if (file.getName().contains("tapered"))
				myBranch.setValue(SlipAlongRuptureModels.TAPERED);
			else if (file.getName().contains("uniform"))
				myBranch.setValue(SlipAlongRuptureModels.UNIFORM);
			rupSet = new U3InversionConfigFactory().updateRuptureSetForBranch(rupSet, branch);
			rupSet.write(file);
		}
	}
	
	public static List<InversionConstraint> getU3Constraints(FaultSystemRupSet rupSet) throws IOException {
		return getU3Generator(rupSet).getConstraints();
	}
	
	public static UCERF3InversionInputGenerator getU3Generator(FaultSystemRupSet rupSet) throws IOException {
		LogicTreeBranch<?> branch = rupSet.requireModule(LogicTreeBranch.class);
		
		InversionTargetMFDs targetMFDs = rupSet.requireModule(InversionTargetMFDs.class);
		FaultModels fm = branch.getValue(FaultModels.class);
		UCERF3InversionConfiguration config = UCERF3InversionConfiguration.forModel(
				branch.getValue(InversionModels.class), rupSet, fm, targetMFDs);
		
		// get the paleo rate constraints
		List<U3PaleoRateConstraint> paleoRateConstraints = CommandLineInversionRunner.getPaleoConstraints(
				fm, rupSet);

		// get the improbability constraints
		double[] improbabilityConstraint = null; // null in UCERF3

		// paleo probability model
		PaleoProbabilityModel paleoProbabilityModel = UCERF3InversionInputGenerator.loadDefaultPaleoProbabilityModel();

		List<U3AveSlipConstraint> aveSlipConstraints = U3AveSlipConstraint.load(rupSet.getFaultSectionDataList());

		return new UCERF3InversionInputGenerator(rupSet, config, paleoRateConstraints, aveSlipConstraints,
				improbabilityConstraint, paleoProbabilityModel);
	}
	
	public static List<InversionConstraint> getStdDevWeightedU3Constraints(
			FaultSystemRupSet rupSet, double slipWeight, double mfdWeight, DoubleUnaryOperator relMFDstdDevFunc,
			double paleoWeight, double parkfieldWeight, double minimizeWeight,
			double u2NuclWeight, double supraSmoothWeight, double mfdSmoothWeight) throws IOException {
		InversionTargetMFDs targetMFDs = rupSet.requireModule(InversionTargetMFDs.class);
		
		List<InversionConstraint> constraints = new ArrayList<>();
		
		// slip rate constraints
		if (slipWeight > 0d) 
			constraints.add(new SlipRateInversionConstraint(slipWeight, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, rupSet));
		
		if (paleoWeight > 0d) {
			// paleo event rate
			List<U3PaleoRateConstraint> paleoRateConstraints =
					UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList());
			constraints.add(new PaleoRateInversionConstraint(
					rupSet, paleoWeight, paleoRateConstraints, UCERF3_PaleoProbabilityModel.load()));
			
			// paleo slip
			List<U3AveSlipConstraint> aveSlipConstraints = U3AveSlipConstraint.load(rupSet.getFaultSectionDataList());
			constraints.add(new PaleoSlipInversionConstraint(
					rupSet, paleoWeight, aveSlipConstraints, U3AveSlipConstraint.slip_prob_model, true));
		}
		
		// parkfield
		if (parkfieldWeight > 0d) {
			double parkfieldMeanRate = 1.0/25.0; // Bakun et al. (2005)
			double parkfieldStdDev = 0.1d*parkfieldMeanRate; // TODO made up
			
			// Find Parkfield M~6 ruptures
			List<Integer> parkfieldRups = UCERF3InversionInputGenerator.findParkfieldRups(rupSet);
			constraints.add(new ParkfieldInversionConstraint(parkfieldWeight, parkfieldMeanRate, parkfieldRups,
					ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, parkfieldStdDev));
		}
		
		if (relMFDstdDevFunc != null && mfdWeight > 0d) {
			List<UncertainIncrMagFreqDist> mfds = new ArrayList<>();
			for (IncrementalMagFreqDist mfd : targetMFDs.getMFD_Constraints())
				mfds.add(UncertainIncrMagFreqDist.relStdDev(mfd, relMFDstdDevFunc));
			constraints.add(new MFDInversionConstraint(rupSet, mfdWeight, false,
					ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, mfds));
		}
		
		if (minimizeWeight > 0d) {
			ModSectMinMags modMinMags = rupSet.requireModule(ModSectMinMags.class);
			Preconditions.checkNotNull(modMinMags, "Rupture set must supply ModSectMinMags if minimization constraint is enabled");
			List<Integer> belowMinIndexes = new ArrayList<>();
			for (int r=0; r<rupSet.getNumRuptures(); r++)
//				if (rupSet.isRuptureBelowSectMinMag(r))
				if (FaultSystemRupSetCalc.isRuptureBelowSectMinMag(rupSet, r, modMinMags))
					belowMinIndexes.add(r);
			constraints.add(new RupRateMinimizationConstraint(minimizeWeight, belowMinIndexes));
		}
		
		ArrayList<U3SectionMFD_constraint> sectMFDConstraints = null;
		if (u2NuclWeight > 0d) {
			sectMFDConstraints = FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(rupSet);
			constraints.add(new U3MFDSubSectNuclInversionConstraint(rupSet, u2NuclWeight, sectMFDConstraints));
		}
		
		if (supraSmoothWeight > 0d)
			constraints.add(new LaplacianSmoothingInversionConstraint(rupSet, supraSmoothWeight));
		
		if (mfdSmoothWeight > 0d) {
			if (sectMFDConstraints == null)
				sectMFDConstraints = FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(rupSet);
			
			constraints.add(new MFDLaplacianSmoothingInversionConstraint(rupSet, mfdSmoothWeight, sectMFDConstraints));
		}

		return constraints;
	}

}
