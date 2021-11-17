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
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SubSeismoOnFaultMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.WaterLevelRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

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
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;

public class InversionDemo {

	public static void main(String[] args) throws IOException {
		
		
		
		GutenbergRichterMagFreqDist gr1 = new GutenbergRichterMagFreqDist(6.05, 8.25, 23);
		
		gr1.setAllButBvalue(6.05, 6.75, 3.7838799937196456E14, 2.0749095531599396E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.75, 4.3127491137557506E14, 2.2976760901106445E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.55, 4.3842269539330035E15, 6.702614738621997E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 7.25, 5.2373404698224336E16, 7.135858032044776E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 7.25, 1.15140045094023472E17, 0.0016908489705185934); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.65, 1.1428453352452396E16, 0.0013149259891796574); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.65, 1.1428453352452396E16, 0.0013149259891796574); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.65, 1.0778898967764078E16, 0.0011653382563877287); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.75, 8.935645534028222E13, 6.950777655367173E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.75, 9.270216796769267E13, 7.412982416872361E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.634969197884943E15, 0.0010499640528587298); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.696994789822025E15, 0.001070784775408472); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.802500050112605E15, 0.0011053092240358366); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.735216226569275E15, 0.0010722499993214023); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.515896677451394E15, 0.0010050635975455081); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.55, 6.439272005933971E15, 9.834252838071306E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.130839919549189E14, 5.92432189446333E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.203043943411271E14, 7.340182629539844E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.2226659695346742E14, 7.107652717151501E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.2240791557929788E14, 7.307239492717735E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.228264243183299E14, 7.501830796028453E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.2049163438800512E14, 6.920870423981034E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.8500000000000005, 1.193210418516083E14, 6.4412167984721925E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.55, 3.4560126871563685E15, 4.3106629027510936E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.55, 3.713817753851802E15, 4.984073469126277E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.55, 3.713817753851802E15, 4.984073469126277E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.45, 8.618233991387049E15, 0.0019882238349078775); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.45, 7.84129178797964E15, 0.0017084300688041296); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.55, 5.052629024040283E15, 7.469748070921114E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.15, 7.65, 1.71956884071431744E17, 6.461451302818004E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.15, 7.65, 1.71974067776472928E17, 6.473642185443594E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.15, 7.65, 1.71742471688573152E17, 6.478963948196094E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.15, 7.65, 1.7167982182730688E17, 6.503235398276339E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.15, 7.65, 1.7168911871366512E17, 6.539782869859572E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.75, 7.3500000000000005, 6.0104921391823432E16, 5.62262706726246E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.75, 7.3500000000000005, 6.0213856361001632E16, 5.664962908936068E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.75, 7.3500000000000005, 4.8813401293454208E16, 4.0668833643636395E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.75, 7.3500000000000005, 4.4358829837110568E16, 3.547926891138277E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 6.75, 8.515776456550969E15, 7.219893648958316E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.45, 5.862244902414524E15, 0.001040401436677825); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.45, 5.862244902414524E15, 0.001040401436677825); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.3500000000000005, 6.45, 5.568783939499474E15, 9.485690620761399E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.95, 9.934526585968634E13, 4.218543473669787E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 6.65, 2.923885303861126E15, 3.299029879613995E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 6.65, 3.113614729526637E15, 3.657481180554488E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.15, 3.6827807922915144E14, 1.8913530819756328E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.15, 4.2002263229122644E14, 2.3153302187783949E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.15, 4.2002263229122644E14, 2.3153302187783949E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.05, 6.15, 3.6590557333680525E14, 1.8719135816657454E-4); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 7.3500000000000005, 2.0416729918605388E14, 1.9077663329550877E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 7.3500000000000005, 2.1995352513149644E14, 2.2922099719415666E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 7.3500000000000005, 2.1708147043380566E14, 1.9836983897712592E-6); System.out.println(gr1.get_bValue());
		gr1.setAllButBvalue(6.25, 6.45, 5.422876084333603E14, 1.1540123576378476E-4); System.out.println(gr1.get_bValue());
		
		System.exit(0);

		
		
		
		double[] bArray = {-1, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 0.896, 0.894, -0.894, -0.896};
//		double[] bArray = {0.5 };
		
		for(double b:bArray) {
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(b, 1.0 , 6.05, 8.25, 23);
//			double moRate = gr.getTotalMomentRate();
//			System.out.println(gr);
//			System.out.println((float)moRate);
//			gr.setAllButBvalue(6.05, 8.25, moRate, 1.0);
//			double b_computed = gr.get_bValue();
			
			double b_computed = gr.compute_bValueAlt(6.05, 8.25);
			
//			double highRate = gr.getCumRate(7.25);
//			b_computed = Math.log10((gr.getTotalIncrRate()-highRate)/highRate);
			
			System.out.println(b+"\t"+(float)b_computed +"\t"+(float)(b_computed/b)+"\t"+(double)(gr.getTotalIncrRate()/1.0));
		}
		
		System.exit(0);

		
		
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(6.35,8.55, 23);
		System.out.println(mfd);

		mfd.set(6.35, 5.9339914E-6);
		mfd.set(6.45, 0.0);
		mfd.set(6.55, 4.1930953E-6);
		mfd.set(6.65, 2.0332966E-6);
		mfd.set(6.75, 7.737724E-6);
		mfd.set(6.85, 1.6636932E-5);
		mfd.set(6.95, 2.1538113E-5);
		mfd.set(7.05, 2.5683565E-5);
		mfd.set(7.15, 1.8583287E-5);
		mfd.set(7.25, 1.6224169E-5);
		mfd.set(7.35, 1.0909239E-5);
		mfd.set(7.45, 2.4610978E-5);
		mfd.set(7.55, 4.4499055E-5);
		mfd.set(7.65, 4.9235724E-4);
		mfd.set(7.75, 4.7000538E-4);
		mfd.set(7.85, 9.2351076E-4);
		mfd.set(7.95, 9.881476E-4);
		mfd.set(8.05, 8.211161E-4);
		mfd.set(8.15, 6.1357504E-4);
		mfd.set(8.25, 3.7129354E-4);
		mfd.set(8.35, 3.7080667E-5);
		mfd.set(8.45, 3.1506759E-6);
		mfd.set(8.55, 1.6095449E-6);
		
		System.out.println(mfd.getGR_fit(6.35,8.55));
		
		System.out.println(mfd.getTotalMomentRate()+"\t"+mfd.getGR_fit(6.35,8.55).getTotalMomentRate());

		
		
		System.exit(0);
		
		
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
