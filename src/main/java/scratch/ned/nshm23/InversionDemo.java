package scratch.ned.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.*;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ThreadedSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.ProgressTrackingCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.InitialSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SubSeismoOnFaultMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.WaterLevelRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.griddedSeismicity.UCERF3_GridSourceGenerator;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionTargetMFDs;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;

public class InversionDemo {

	public static void main(String[] args) throws IOException {
		// output directory, change this
		File outputDir = new File("/tmp/inversion_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		// number of annealing threads
		int threads = Runtime.getRuntime().availableProcessors();
		
		// if false, and rupture_set.zip already exists in the output directory, will just load it instead of rebuilding it
		boolean rebuildRupSet = true;
		
		// UCERF3 logic tree branch, starts with our reference branch
		U3LogicTreeBranch branch = U3LogicTreeBranch.DEFAULT;
		// I've been using the uniform slip model
		branch.setValue(SlipAlongRuptureModels.UNIFORM);
		
		FaultModels fm = branch.getValue(FaultModels.class);
		ScalingRelationships scale = branch.getValue(ScalingRelationships.class);
		
//		// this will generate a new Coulomb (draft U4) rupture set, and may take a little while
//		RupSetConfig rsConfig = new RuptureSets.CoulombRupSetConfig(fm , scale);
		// this will generate a UCERF3-style rupture set
		RupSetConfig rsConfig = new RuptureSets.U3RupSetConfig(fm, scale);
		
		File rsFile = new File(outputDir, "rupture_set.zip");
		FaultSystemRupSet rupSet;
		if (!rebuildRupSet && rsFile.exists()) {
			// load existing rupture set
			rupSet = FaultSystemRupSet.load(rsFile);
		} else {
			// build the rupture set
			rupSet = rsConfig.build(threads);
			
			// now add all of the things it needs to run a UCERF3-style inversion (e.g., polygons, target MFDs)
			rupSet = FaultSystemRupSet.buildFromExisting(rupSet).forU3Branch(branch).build();
			
			// write it out
			rupSet.write(new File(outputDir, "rupture_set.zip"));
		}
		
		// choose your SA parameters
		
		// Perturbation function:
		// what we used in UCREF3, not that great
//		GenerationFunctionType perturb = GenerationFunctionType.UNIFORM_NO_TEMP_DEPENDENCE;
		// new (better) exponential scale version that chooses values across a wide (log-spaced) range
//		GenerationFunctionType perturb = GenerationFunctionType.EXPONENTIAL_SCALE;
		// smarter version that uses the starting solution to choose random perturbations that are more relevant
		// to rupture, but needs an extra step below
		GenerationFunctionType perturb = GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE;
		
		// Water level:
		// In UCERF3, water level was applied as 1e-2*smoothStartingModel
//		double wlFract = 1e-2;
		// or we can disable it
		double wlFract = 0d;
		
		// Non-negativity constraint:
		
		// this is the one that we used in UCERF3. it does it's stated job of limiting zero rates, but that really just means
		// that a lot of ruptures get stuck at near-zero rates when the inversion wants them to be zero
//		NonnegativityConstraintType nonNeg = NonnegativityConstraintType.LIMIT_ZERO_RATES;
		// this one allows the inversion to set rupture rates back to zero
		NonnegativityConstraintType nonNeg = NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN;
		
		// completion criteria
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(5);
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(1);
		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(30);
//		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(5);
		
		// now configure inversion inputs
		
		// this is the UCERF3-style target MFDs
		InversionTargetMFDs targetMFDs = rupSet.requireModule(InversionTargetMFDs.class);
		// UCERF3 inversion parameters, will be replaced for NSHM23
		UCERF3InversionConfiguration config = UCERF3InversionConfiguration.forModel(
				branch.getValue(InversionModels.class), rupSet, fm, targetMFDs);
		config.setMinimumRuptureRateFraction(wlFract);
		// you can set other inversion options or change weights here. e.g., to disable paleoseismic rates:
//		config.setPaleoRateConstraintWt(0d);
//		config.setPaleoSlipWt(0d);
		
		// various data inputs that we needed in UCERF3. in the future, these could be attached to the rupture set as modules
		
		// get the paleo rate constraints
		List<U3PaleoRateConstraint> paleoRateConstraints = CommandLineInversionRunner.getPaleoConstraints(
				fm, rupSet);

		// get the improbability constraints
		double[] improbabilityConstraint = null; // this wasn't used in UCERF3

		// paleo probability model
		PaleoProbabilityModel paleoProbabilityModel = UCERF3InversionInputGenerator.loadDefaultPaleoProbabilityModel();

		List<U3AveSlipConstraint> aveSlipConstraints = U3AveSlipConstraint.load(rupSet.getFaultSectionDataList());

		// this actually builds the constraint matrices
		UCERF3InversionInputGenerator inputGen = new UCERF3InversionInputGenerator(
				rupSet, config, paleoRateConstraints, aveSlipConstraints, improbabilityConstraint, paleoProbabilityModel);
		
		System.out.println("Generating inversion inputs");
		inputGen.generateInputs();
		// this makes the inversion faster
		inputGen.columnCompress();
		
		// wrapping the completion criteria in this allows us to track inversion progress
		ProgressTrackingCompletionCriteria progress = new ProgressTrackingCompletionCriteria(completion);
		// time between across-thread synchronization, if threads > 1
		TimeCompletionCriteria subCompletion = TimeCompletionCriteria.getInSeconds(1);
		
		// build the SA calculator
		SimulatedAnnealing sa;
		if (threads > 1)
			sa = new ThreadedSimulatedAnnealing(inputGen.getA(), inputGen.getD(),
					inputGen.getInitialSolution(), 0d, inputGen.getA_ineq(), inputGen.getD_ineq(), threads, subCompletion);
		else
			sa = new SerialSimulatedAnnealing(inputGen.getA(), inputGen.getD(),
					inputGen.getInitialSolution(), 0d, inputGen.getA_ineq(), inputGen.getD_ineq());
		
		// if we enabled the variable perturbation function, we need to tell it about the smooth initial solution
		// it will choose perturbations that are +/- 2 orders of magnitude around these values 
		if (perturb == GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE)
			sa.setVariablePerturbationBasis(config.getMinimumRuptureRateBasis());
		
		// keep track of information on individual constraints
		progress.setConstraintRanges(inputGen.getConstraintRowRanges());
		sa.setConstraintRanges(inputGen.getConstraintRowRanges());
		
		// set SA params
		sa.setPerturbationFunc(perturb);
		sa.setNonnegativeityConstraintAlgorithm(nonNeg);
		
		// actually run the inversion
		sa.iterate(progress);
		
		double[] rawSol = sa.getBestSolution();
		// this does nothing if we don't have a water level
		double[] rates = inputGen.adjustSolutionForWaterLevel(rawSol);
		
		FaultSystemSolution sol = new FaultSystemSolution(rupSet, rates);
		
		// This solution has everything we need for standard fault-based hazard calculations, but we'll add extra modules
		// here, including gridded seismicity and some inversion metadata
		
		// add sub-seismo MFDs to the solution (used in some plots, or by ETAS)
		sol.addModule(targetMFDs.getOnFaultSubSeisMFDs());
		// add grid source provider
		sol.setGridSourceProvider(new UCERF3_GridSourceGenerator(sol, branch.getValue(SpatialSeisPDF.class),
				branch.getValue(MomentRateFixes.class), targetMFDs,
				sol.requireModule(SubSeismoOnFaultMFDs.class),
				branch.getValue(MaxMagOffFault.class).getMaxMagOffFault(),
				rupSet.requireModule(FaultGridAssociations.class)));
		// add inversion progress
		sol.addModule(progress.getProgress());
		// add water level rates
		if (inputGen.getWaterLevelRates() != null)
			sol.addModule(new WaterLevelRates(inputGen.getWaterLevelRates()));
		if (inputGen.hasInitialSolution())
			sol.addModule(new InitialSolution(inputGen.getInitialSolution()));
		
		// write solution
		sol.write(new File(outputDir, "solution.zip"));
		
		// now write a solution report
		PlotLevel level = PlotLevel.FULL; // full set of plots, can be slow (mostly due to fault-by-fault pages)
//		PlotLevel level = PlotLevel.DEFAULT; // quicker set of plots
		ReportPageGen solReport = new ReportPageGen(rupSet, sol, "Test Inversion", outputDir,
				ReportPageGen.getDefaultSolutionPlots(level));
	 
		solReport.generatePage();
		
		System.out.println("DONE");
		System.exit(0);
	}

}
