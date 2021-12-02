package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Supplier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.modules.ArchivableModule;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.AbstractExternalFetcher;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.AsyncPostBatchHook;
import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_LogicTreeInversionRunner extends MPJTaskCalculator {

	private File outputDir;
	
	private int runsPerBranch = 1;
	
	private InversionConfigurationFactory factory;

	private LogicTree<LogicTreeNode> tree;

	private CommandLine cmd;
	
	private int annealingThreads;
	private int runsPerBundle = 1;

	public MPJ_LogicTreeInversionRunner(CommandLine cmd) throws IOException {
		super(cmd);
		this.cmd = cmd;
		this.annealingThreads = Integer.parseInt(cmd.getOptionValue("annealing-threads"));
		Preconditions.checkState(annealingThreads >= 1);
		
		if (cmd.hasOption("runs-per-bundle"))
			runsPerBundle = Integer.parseInt(cmd.getOptionValue("runs-per-bundle"));
		
		this.shuffle = false;
		
		tree = LogicTree.read(new File(cmd.getOptionValue("logic-tree")));
		if (rank == 0)
			debug("Loaded "+tree.size()+" tree nodes");
		
		outputDir = new File(cmd.getOptionValue("output-dir"));
		
		if (rank == 0)
			waitOnDir(outputDir, 5, 1000);
		
		if (cmd.hasOption("runs-per-branch"))
			runsPerBranch = Integer.parseInt(cmd.getOptionValue("runs-per-branch"));
		
		try {
			@SuppressWarnings("unchecked")
			Class<? extends InversionConfigurationFactory> factoryClass = (Class<? extends InversionConfigurationFactory>)
					Class.forName(cmd.getOptionValue("inversion-factory"));
			factory = factoryClass.getDeclaredConstructor().newInstance();
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		debug("Factory type: "+factory.getClass().getName());
		
		if (rank == 0) {
			SolutionLogicTree treeWriter = factory.initSolutionLogicTree(tree);
			if (treeWriter != null) {
				this.postBatchHook = new AsyncLogicTreeWriter(treeWriter, cmd.hasOption("branch-average"));
			}
		}
	}
	
	private class AsyncLogicTreeWriter extends AsyncPostBatchHook {
		
		private AsyncSolutionLogicTree asyncWriter;
		private CompletableFuture<?> archiveFuture;
		private CompletableFuture<?> branchAvgFuture;
		private boolean branchAverage;
		
		public AsyncLogicTreeWriter(SolutionLogicTree treeWriter, boolean branchAverage) {
			super(1);
			this.branchAverage = branchAverage;
			asyncWriter = new AsyncSolutionLogicTree(treeWriter);
		}

		@Override
		protected void batchProcessedAsync(int[] batch, int processIndex) {
			debug("AsyncLogicTree: beginning async call with batch size "
					+batch.length+" from "+processIndex+": "+getCountsString());
			for (int index : batch)
				asyncWriter.calcDone(index);
			
			if (archiveFuture == null) {
				synchronized (this) {
					if (archiveFuture == null) {
						File outputFile = new File(outputDir+".zip");
						debug("AsyncLogicTree: submitting async write future");
						archiveFuture = CompletableFuture.runAsync(new AsyncWriteRunnable(asyncWriter, outputFile));
						if (branchAverage) {
							File baFile = new File(outputDir.getParent(), "branch_averaged.zip");
							branchAvgFuture = CompletableFuture.runAsync(new AsyncBranchAvgRunnable(asyncWriter, baFile));
						}
						debug("AsyncLogicTree: submitted async write future");
					}
				}
			}
			debug("AsyncLogicTree: exiting async process, stats: "+getCountsString());
		}

		@Override
		public void shutdown() {
			super.shutdown();
			
			debug("AsyncLogicTree: waiting on async write future");
			try {
				archiveFuture.get();
				if (branchAvgFuture != null)
					branchAvgFuture.get();
			} catch (InterruptedException | ExecutionException e) {
				abortAndExit(e, 1);
			}
		}
		
	}
	
	private class AsyncWriteRunnable implements Runnable {

		private AsyncSolutionLogicTree asyncWriter;
		private File outputFile;

		public AsyncWriteRunnable(AsyncSolutionLogicTree asyncWriter, File outputFile) {
			this.asyncWriter = asyncWriter;
			this.outputFile = outputFile;
		}

		@Override
		public void run() {
			debug("AsyncLogicTree: beginning write to "+outputFile.getAbsolutePath());
			ModuleArchive<AsyncSolutionLogicTree> compoundArchive = new ModuleArchive<>();
			compoundArchive.addModule(asyncWriter);
			try {
				compoundArchive.write(outputFile);
			} catch (IOException e) {
				abortAndExit(e, 1);
			}
			debug("AsyncLogicTree: DONE writing");
		}
		
	}
	
	private class AsyncBranchAvgRunnable implements Runnable {

		private AsyncSolutionLogicTree asyncWriter;
		private File outputFile;

		public AsyncBranchAvgRunnable(AsyncSolutionLogicTree asyncWriter, File outputFile) {
			this.asyncWriter = asyncWriter;
			this.outputFile = outputFile;
		}

		@Override
		public void run() {
			debug("AsyncLogicTree: beginning write to "+outputFile.getAbsolutePath());
			FaultSystemSolution avgSol;
			try {
				avgSol = asyncWriter.calcBranchAveraged();
			} catch (Exception e) {
				debug("Failed to build branch averaged solution, see error below");
				System.out.flush();
				e.printStackTrace();
				System.err.flush();
				return;
			}
			try {
				avgSol.write(outputFile);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			debug("AsyncLogicTree: DONE writing branch averaged");
		}
		
	}
	
	private class AsyncSolutionLogicTree extends AbstractExternalFetcher {

		private SolutionLogicTree treeWriter;
		
		private boolean[] dones;
		private Map<LogicTreeBranch<?>, CompletableFuture<Callable<FaultSystemSolution>>> branchFutureMap;

		protected AsyncSolutionLogicTree(SolutionLogicTree treeWriter) {
			super(treeWriter.getLogicTree());
			this.treeWriter = treeWriter;
			
			dones = new boolean[getNumTasks()];
			branchFutureMap = new HashMap<>();
			
			for (LogicTreeBranch<?> branch : treeWriter.getLogicTree())
				branchFutureMap.put(branch, new CompletableFuture<>());
		}
		
		void calcDone(int index) {
			int branchIndex = branchForCalcIndex(index);
			LogicTreeBranch<?> branch = getLogicTree().getBranch(branchIndex);
			
			debug("AsyncLogicTree: calcDone "+index+" = branch "+branchIndex+": "+branch);
			
			synchronized (dones) {
				dones[index] = true;
				List<File> solFiles = new ArrayList<>();
				for (int run=0; run<runsPerBranch; run++) {
					int doneIndex = indexForBranchRun(branchIndex, run);
					if (dones[doneIndex]) {
						solFiles.add(getSolFile(branch, run));
					} else {
						// not all runs for this branch are done
						debug("AsyncLogicTree: not ready, waiting on run "+run+" for branch "+branchIndex
								+" (origIndex="+index+", checkIndex="+doneIndex+"): "+branch);
						return;
					}
				}
				
				Callable<FaultSystemSolution> call;
				if (runsPerBranch > 1) {
					call = new Callable<FaultSystemSolution>() {
						
						@Override
						public FaultSystemSolution call() throws Exception {
							FaultSystemSolution[] inputs = new FaultSystemSolution[solFiles.size()];
							for (int i=0; i<inputs.length; i++)
								inputs[i] = FaultSystemSolution.load(solFiles.get(i));
							return AverageSolutionCreator.buildAverage(inputs);
						}
					};
				} else {
					call = new Callable<FaultSystemSolution>() {

						@Override
						public FaultSystemSolution call() throws Exception {
							return FaultSystemSolution.load(solFiles.get(0));
						}
						
					};
				}
				
				CompletableFuture<Callable<FaultSystemSolution>> future = branchFutureMap.get(branch);
				Preconditions.checkNotNull(future);
				debug("AsyncLogicTree: calcDone "+index+" = branch "+branchIndex+" is READY: "+branch);
				future.complete(call);
				debug("AsyncLogicTree: calcDone "+index+" = branch "+branchIndex+" is COMPLETE: "+branch);
			}
		}

		@Override
		protected FaultSystemSolution loadExternalForBranch(LogicTreeBranch<?> branch) throws IOException {
			try {
				debug("AsyncLogicTree: waiting on "+branch);
				Callable<FaultSystemSolution> solCall = branchFutureMap.get(branch).get();
				debug("AsyncLogicTree: loading "+branch);
				return solCall.call();
			} catch (Exception e) {
				System.err.println("FAILED to load Async!!");
				e.printStackTrace();
				abortAndExit(e, 2);
				throw ExceptionUtils.asRuntimeException(e);
			}
		}

		@Override
		public List<? extends LogicTreeLevel<?>> getLevelsAffectingFile(String fileName) {
			return treeWriter.getLevelsAffectingFile(fileName);
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForFaultSections() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForRuptureSectionIndices() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForRuptureProperties() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForRuptureRates() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForGridRegion() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForGridMechs() {
			throw new IllegalStateException();
		}

		@Override
		protected List<? extends LogicTreeLevel<?>> getLevelsForGridMFDs() {
			throw new IllegalStateException();
		}

		@Override
		public Class<? extends ArchivableModule> getLoadingClass() {
			return treeWriter.getLoadingClass();
		}
		
	}
	
	private static void waitOnDir(File dir, int maxRetries, long sleepMillis) {
		int retry = 0;
		while (!(dir.exists() || dir.mkdir())) {
			try {
				Thread.sleep(sleepMillis);
			} catch (InterruptedException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (retry++ > maxRetries)
				throw new IllegalStateException("Directory doesn't exist and couldn't be created after "
						+maxRetries+" retries: "+dir.getAbsolutePath());
		}
	}

	@Override
	protected int getNumTasks() {
		return tree.size()*runsPerBranch;
	}
	
	protected File getSolFile(LogicTreeBranch<?> branch, int run) {
		String dirName = branch.buildFileName();
		if (runsPerBranch > 1)
			dirName += "_run "+run;
		File runDir = new File(outputDir, dirName);
		Preconditions.checkState(runDir.exists() || runDir.mkdir());
		
		return new File(runDir, "solution.zip");
	}
	
	private int branchForCalcIndex(int index) {
		return index / runsPerBranch;
	}
	
	private int runForCalcIndex(int index) {
		return index % runsPerBranch;
	}
	
	private int indexForBranchRun(int branchIndex, int runIndex) {
		return branchIndex * runsPerBranch + runIndex;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		if (runsPerBundle > 1) {
			ExecutorService exec = Executors.newFixedThreadPool(runsPerBundle);
			
			List<Future<?>> futures = new ArrayList<>();
			for (int index : batch)
				futures.add(exec.submit(new CalcRunnable(index)));
			
			for (Future<?> future : futures) {
				try {
					future.get();
				} catch (InterruptedException | ExecutionException e) {
					exec.shutdown();
					throw e;
				}
			}
			
			exec.shutdown();
		} else {
			for (int index : batch) {
				new CalcRunnable(index).run();
			}
		}
	}

	private class CalcRunnable implements Runnable {
		
		private int index;

		public CalcRunnable(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			int branchIndex = branchForCalcIndex(index);
			int run = runForCalcIndex(index);
			
			LogicTreeBranch<LogicTreeNode> branch = tree.getBranch(branchIndex);
			
			debug("index "+index+" is branch "+branchIndex+" run "+run+": "+branch);
			
			File solFile = getSolFile(branch, run);
			
			if (solFile.exists()) {
				debug(solFile.getAbsolutePath()+" exists, testing loading...");
				try {
					FaultSystemSolution.load(solFile);
					debug("skipping "+index+" (already done)");
					return;
				} catch (Exception e) {
					debug("Failed to load, re-inverting: "+e.getMessage());
				}
			}
			
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, annealingThreads);
			rupSet.addModule(branch);
			
			InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, cmd, annealingThreads);
			
			debug("Running inversion for task "+index+" with "+config.getConstraints().size()+" constraints");
			FaultSystemSolution sol = Inversions.run(rupSet, config);
			
			try {
				sol.write(solFile);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		if (rank == 0) {
			debug("waiting for any post batch hook operations to finish");
			((AsyncLogicTreeWriter)postBatchHook).shutdown();
			debug("post batch hook done");
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		ops.addRequiredOption("lt", "logic-tree", true, "Path to logic tree JSON file");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addRequiredOption("at", "annealing-threads", true, "Number of annealing threads per inversion");
		ops.addOption("rpb", "runs-per-branch", true, "Runs per branch (default is 1)");
		ops.addOption("rpb", "runs-per-bundle", true, "Simultaneous runs to executure (default is 1)");
		ops.addRequiredOption("ifc", "inversion-factory", true, "Inversion configuration factory classname");
		ops.addOption("ba", "branch-average", false, "Flag to also build a branch-averaged solution");
		
		for (Option op : InversionConfiguration.createSAOptions().getOptions())
			ops.addOption(op);
		
		return ops;
	}

	public static void main(String[] args) {
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_LogicTreeInversionRunner.class);
			
			MPJ_LogicTreeInversionRunner driver = new MPJ_LogicTreeInversionRunner(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
