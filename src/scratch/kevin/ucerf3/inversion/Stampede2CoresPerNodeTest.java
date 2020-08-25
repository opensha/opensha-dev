package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.inversion.laughTest.UCERF3PlausibilityConfig;
import scratch.UCERF3.simulatedAnnealing.ThreadedSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.completion.CompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.ProgressTrackingCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;
import scratch.UCERF3.utils.aveSlip.AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;

public class Stampede2CoresPerNodeTest extends MPJTaskCalculator {
	
	int[] thread_counts;
	
	CompletionCriteria criteria;
	CompletionCriteria subCompetionCriteria;
	
	DoubleMatrix2D A;
	double[] d;
	DoubleMatrix2D A_ineq;
	double[] d_ineq;
	double[] initial;
	double[] minimumRuptureRates;
	
	double relativeSmoothnessWt;
	
	File outputDir;

	public Stampede2CoresPerNodeTest(CommandLine cmd, File outputDir) {
		super(cmd);
		
		this.outputDir = outputDir;
		if (rank == 0)
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		if (rank == 0)
			thread_counts = new int[] { 56 };
		else if (rank == 1 || rank == 2)
			thread_counts = new int[] { 1, 1, 2, 2, 16, 32 };
		else if (rank == 3)
			thread_counts = new int[] { 4, 4, 4, 4, 8, 8, 8, 8 };
		
		UCERF3PlausibilityConfig filter = UCERF3PlausibilityConfig.getDefault();
		InversionModels inversionModel = InversionModels.CHAR_CONSTRAINED;
		
		double defaultAseis = 0.1;
		InversionFaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(
				filter, defaultAseis, inversionModel, FaultModels.FM3_1);
		
		UCERF3InversionConfiguration config = UCERF3InversionConfiguration.forModel(inversionModel, rupSet);
		
		// get the paleo rate constraints
		List<PaleoRateConstraint> paleoRateConstraints = null;
		try {
			paleoRateConstraints = CommandLineInversionRunner.getPaleoConstraints(
					rupSet.getFaultModel(), rupSet);
		} catch (IOException e1) {
			e1.printStackTrace();
			// exit
			System.exit(1);
		}
		
		// get the improbability constraints
		double[] improbabilityConstraint = null; // null for now
//		improbabilityConstraint = getCoulombWeights(faultSystemRupSet.getNumRuptures(), CoulombWeightType.MEAN_SIGMA, precomputedDataDir);
		
		// paleo probability model
		PaleoProbabilityModel paleoProbabilityModel = null;
		try {
			paleoProbabilityModel = UCERF3InversionInputGenerator.loadDefaultPaleoProbabilityModel();
		} catch (IOException e) {
			e.printStackTrace();
			// exit
			System.exit(1);
		}
		
		List<AveSlipConstraint> aveSlipConstraints = null;
		try {
			aveSlipConstraints = AveSlipConstraint.load(rupSet.getFaultSectionDataList());
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// create the input generator
		UCERF3InversionInputGenerator gen = new UCERF3InversionInputGenerator(rupSet, config, paleoRateConstraints, aveSlipConstraints,
				improbabilityConstraint, paleoProbabilityModel);
		
		// generate the inputs
		gen.generateInputs();
		// optionally we can specify the class we want to use for the A matrix:
//		gen.generateInputs(SparseDoubleMatrix2D.class);
		
		// column compress it for fast annealing!
		gen.columnCompress();
		
		// fetch matrices
		A = gen.getA();
		d = gen.getD();
		A_ineq = gen.getA_ineq();
		d_ineq = gen.getD_ineq();
		initial = gen.getInitialSolution();
		minimumRuptureRates = gen.getWaterLevelRates();
		
		// use one of these to run it for a set amount of time: 
//		criteria = TimeCompletionCriteria.getInMinutes(90);
		criteria = TimeCompletionCriteria.getInHours(8);
//		criteria = TimeCompletionCriteria.getInSeconds(30);
		// or use this to run until a set amount of iterations have been completed
//		criteria = new IterationCompletionCriteria(1);
		// this is the "sub completion criteria" - the amount of time (or iterations) between synchronization
		subCompetionCriteria = TimeCompletionCriteria.getInSeconds(1); // 1 second;
		
		relativeSmoothnessWt = config.getSmoothnessWt();
	}

	@Override
	protected int getNumTasks() {
		return size;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		if (thread_counts == null)
			return;
		
		ExecutorService exec = Executors.newFixedThreadPool(thread_counts.length+1);
		
		List<Future<?>> futures = new ArrayList<>();
		
		int myIndex = 0;
		for (int numThreads : thread_counts) {
			String fileName = "node_"+rank+"_calc_"+(myIndex++)+"_threads_"+numThreads+".csv";
			File outputFile = new File(outputDir, fileName);
			
			ProgressTrackingCompletionCriteria criteria =
					new ProgressTrackingCompletionCriteria(this.criteria, outputFile);
			
			final ThreadedSimulatedAnnealing tsa = new ThreadedSimulatedAnnealing(A, d, initial, relativeSmoothnessWt,
					A_ineq, d_ineq, minimumRuptureRates, numThreads, subCompetionCriteria);
			
			futures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					tsa.iterate(criteria);
				}
			}));
		}
		
		for (Future<?> f : futures) {
			f.get();
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		
	}

	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, Stampede2CoresPerNodeTest.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(Stampede2CoresPerNodeTest.class)
						+" <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			
			Stampede2CoresPerNodeTest driver = new Stampede2CoresPerNodeTest(cmd, outputDir);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
