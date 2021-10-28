package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.LaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDLaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDUncertaintyWeightedInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldUncertaintyWeightedInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint.WeightingType;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSystemRupSetCalc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionTargetMFDs;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.MFD_InversionConstraint;
import scratch.UCERF3.utils.MFD_WeightedInversionConstraint;
import scratch.UCERF3.utils.SectionMFD_constraint;
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
//		Thread.sleep(5l*hour + 30l*minute);
//		System.out.println("Im awake! "+new Date());
		
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		// run this if I need to attach UCERF3 modules after rebuilding default rupture sets
//		reprocessDefaultRupSets(parentDir);
//		System.exit(0);

		List<String> argz = new ArrayList<>();

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(date);

//		dirName += "-coulomb-u3";
//		File origRupSetFile = new File(parentDir, "fm3_1_coulomb.zip");

		dirName += "-u3rs";
		File origRupSetFile = new File(parentDir, "fm3_1_u3ref_uniform_reproduce_ucerf3.zip");
		
//		dirName += "-coulomb-nshm23";
//		File origRupSetFile = new File(parentDir, "nshm23_geo_dm_coulomb.zip");

		argz.add("--rupture-set");
		argz.add(origRupSetFile.getAbsolutePath());

		argz.add("--threads");
		argz.add("16");
		argz.add("--avg-threads");
		argz.add("4");

		dirName += "-slip_constr";
		argz.add("--slip-constraint");
//		argz.add("--slip-weight"); argz.add("1");
//		argz.add("--norm-slip-weight"); argz.add("0.01");
		argz.add("--uncertain-slip-weight");
		argz.add("1");
		dirName += "_uncertain";

//		double b = 1;
//		dirName += "-rel_gr_b"+oDF.format(b);
//		argz.add("--rel-gr-constraint");
//		argz.add("--b-value"); argz.add(oDF.format(b));
////		argz.add("--mfd-weight"); argz.add("1000"); dirName += "_wt1000";
////		argz.add("--mfd-ineq"); dirName += "_ineq";

//		dirName += "-u3_mfd";
//		argz.add("--mfd-constraint");
//		argz.add("--mfd-transition-mag"); argz.add("7.8"); dirName += "_trans7.8";

//		double b = 1;
//		dirName += "-infer_gr_b"+oDF.format(b);
//		argz.add("--mfd-constraint");
//		argz.add("--infer-target-gr");
//		argz.add("--b-value"); argz.add(oDF.format(b));
////		argz.add("--mfd-ineq"); dirName += "_ineq";
//		argz.add("--mfd-transition-mag"); argz.add("7.8"); dirName += "_trans7.8";

		dirName += "-smooth";
		argz.add("--smooth");
		argz.add("--smooth-weight");
		argz.add("1000");

//		dirName += "-minimize_below";
//		argz.add("--minimize-below-sect-min");
//		argz.add("--minimize-weight");
//		argz.add("10000");
		
		boolean u3Constraints = false;
//		boolean u3Constraints = true; dirName += "-u3_constraints";
		boolean u3StdDevConstraints = false;
//		boolean u3StdDevConstraints = true; dirName += "-u3_std_dev_tests";

//		dirName += "-5h";
//		argz.add("--completion"); argz.add("5h");
//		argz.add("-avg-completion"); argz.add("5m");
//		dirName += "-1h";
//		argz.add("--completion"); argz.add("1h");
//		argz.add("-avg-completion"); argz.add("5m");
//		dirName += "-30m";
//		argz.add("--completion"); argz.add("30m");
//		argz.add("-avg-completion"); argz.add("1m");
		dirName += "-10m";
		argz.add("--completion"); argz.add("10m");
		argz.add("-avg-completion"); argz.add("1m");

		List<InversionConstraint> extraConstraints = new ArrayList<>();
		
		if (u3Constraints) {
			FaultSystemRupSet rupSet = FaultSystemRupSet.load(origRupSetFile);
			extraConstraints.addAll(getU3Constraints(rupSet));
		}
		
		if (u3StdDevConstraints) {
			FaultSystemRupSet rupSet = FaultSystemRupSet.load(origRupSetFile);
			DoubleUnaryOperator mfdStdDevFunc = M->0.1;
//			DoubleUnaryOperator mfdStdDevFunc = M->Math.max(0.1, 0.1*(M-5));
//			DoubleUnaryOperator mfdStdDevFunc = M->0.1+Math.pow(10, M-8); dirName += "-aggresivePowMSD";
			extraConstraints.addAll(getStdDevWeightedU3Constraints(rupSet, mfdStdDevFunc, 0, 0, 0, 0d));
		}

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
				ReportPageGen report = new ReportPageGen(sol.getRupSet(), sol, name, reportDir,
						ReportPageGen.getDefaultSolutionPlots(plots));
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
				myBranch.setValue(defaultDM);
			} else {
				myBranch.clearValue(FaultModels.class);
				myBranch.clearValue(DeformationModels.class);
			}
			if (file.getName().contains("tapered"))
				myBranch.setValue(SlipAlongRuptureModels.TAPERED);
			else if (file.getName().contains("uniform"))
				myBranch.setValue(SlipAlongRuptureModels.UNIFORM);
			rupSet = FaultSystemRupSet.buildFromExisting(rupSet).forU3Branch(myBranch).build();
			rupSet.write(file);
		}
	}
	
	public static List<InversionConstraint> getU3Constraints(FaultSystemRupSet rupSet) throws IOException {
		return getU3Generator(rupSet).getConstraints();
	}
	
	public static UCERF3InversionInputGenerator getU3Generator(FaultSystemRupSet rupSet) throws IOException {
		U3LogicTreeBranch branch = rupSet.requireModule(U3LogicTreeBranch.class);
		
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
			FaultSystemRupSet rupSet, DoubleUnaryOperator relMFDstdDevFunc, double minimizeWeight,
			double u2NuclWeight, double supraSmoothWeight, double mfdSmoothWeight) throws IOException {
		InversionTargetMFDs targetMFDs = rupSet.requireModule(InversionTargetMFDs.class);
		
		List<InversionConstraint> constraints = new ArrayList<>();
		
		final double weightEach = 1d;
		
		// slip rate constraints
		constraints.add(new SlipRateInversionConstraint(weightEach, WeightingType.NORMALIZED_BY_UNCERTAINTY, rupSet));
		
		// paleo event rate
		List<U3PaleoRateConstraint> paleoRateConstraints =
				UCERF3_PaleoRateConstraintFetcher.getConstraints(rupSet.getFaultSectionDataList());
		constraints.add(new PaleoRateInversionConstraint(
				rupSet, weightEach, paleoRateConstraints, UCERF3_PaleoProbabilityModel.load()));
		
		// paleo slip
		List<U3AveSlipConstraint> aveSlipConstraints = U3AveSlipConstraint.load(rupSet.getFaultSectionDataList());
		constraints.add(new PaleoSlipInversionConstraint(
				rupSet, weightEach, aveSlipConstraints, U3AveSlipConstraint.slip_prob_model, true));
		
		// parkfield
		double parkfieldMeanRate = 1.0/25.0; // Bakun et al. (2005)
		double parkfieldStdDev = 0.1d*parkfieldMeanRate; // TODO made up
		
		// Find Parkfield M~6 ruptures
		List<Integer> parkfieldRups = UCERF3InversionInputGenerator.findParkfieldRups(rupSet);
		constraints.add(new ParkfieldUncertaintyWeightedInversionConstraint(
				weightEach, parkfieldMeanRate, parkfieldStdDev, parkfieldRups));
		
		if (relMFDstdDevFunc != null)
			constraints.add(new MFDInversionConstraint(rupSet, weightEach, false,
					MFDInversionConstraint.WeightingType.NORMALIZED_BY_UNCERTAINTY, targetMFDs.getMFD_Constraints(),
					MFDInversionConstraint.calcStdDevsFromRelativeFunc(targetMFDs.getMFD_Constraints(), relMFDstdDevFunc)));
		
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
		
		ArrayList<SectionMFD_constraint> sectMFDConstraints = null;
		if (u2NuclWeight > 0d) {
			sectMFDConstraints = FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(rupSet);
			constraints.add(new MFDSubSectNuclInversionConstraint(rupSet, u2NuclWeight, sectMFDConstraints));
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
