package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.AbstractAsyncLogicTreeWriter;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeInversionRunner;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InfoModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.SolutionProcessor;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.ShawSegmentationModels;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_StrictSegLogicTreeTranslation extends MPJTaskCalculator {

	private File inputDir;
	private File outputDir;
	
	double minDistance;
	
	private InversionConfigurationFactory factory;

	private LogicTree<LogicTreeNode> inputTree;
	private LogicTree<LogicTreeNode> outputTree;

	public MPJ_StrictSegLogicTreeTranslation(CommandLine cmd) throws IOException {
		super(cmd);
		
		this.shuffle = false;
		
		inputTree = LogicTree.read(new File(cmd.getOptionValue("input-logic-tree")));
		outputTree = LogicTree.read(new File(cmd.getOptionValue("output-logic-tree")));
		if (rank == 0)
			debug("Loaded "+outputTree.size()+" tree nodes");

		inputDir = new File(cmd.getOptionValue("input-dir"));
		outputDir = new File(cmd.getOptionValue("output-dir"));
		
		if (cmd.hasOption("min-distance"))
			minDistance = Double.parseDouble(cmd.getOptionValue("min-distance"));
		
		if (rank == 0)
			waitOnDir(outputDir, 5, 1000);
		
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
			this.postBatchHook = new AsyncLogicTreeWriter(factory.getSolutionLogicTreeProcessor());
		}
	}
	
	private class AsyncLogicTreeWriter extends AbstractAsyncLogicTreeWriter {
		
		public AsyncLogicTreeWriter(SolutionProcessor processor) {
			super(outputDir, processor, outputTree.getWeightProvider());
		}

		@Override
		public int getNumTasks() {
			return MPJ_StrictSegLogicTreeTranslation.this.getNumTasks();
		}

		@Override
		public void debug(String message) {
			MPJ_StrictSegLogicTreeTranslation.this.debug(message);
		}

		@Override
		public LogicTreeBranch<?> getBranch(int calcIndex) {
			return outputTree.getBranch(calcIndex);
		}

		@Override
		public FaultSystemSolution getSolution(LogicTreeBranch<?> branch, int calcIndex) throws IOException {
			File solFile = getOutputSolFile(branch);
			return FaultSystemSolution.load(solFile);
		}

		@Override
		public void abortAndExit(Throwable t, int exitCode) {
			MPJ_LogicTreeInversionRunner.abortAndExit(t, exitCode);
		}
		
	}
	
	private void memoryDebug(String info) {
		if (info == null || info.isBlank())
			info = "";
		else
			info += "; ";
	    
	    System.gc();
		Runtime rt = Runtime.getRuntime();
		long totalMB = rt.totalMemory() / 1024 / 1024;
		long freeMB = rt.freeMemory() / 1024 / 1024;
		long usedMB = totalMB - freeMB;
		debug(info+"mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
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
		return outputTree.size();
	}
	
	protected File getInputSolFile(LogicTreeBranch<?> branch) {
		String dirName = branch.buildFileName();
		File runDir = new File(inputDir, dirName);
		Preconditions.checkState(runDir.exists() || runDir.mkdir());
		
		return new File(runDir, "solution.zip");
	}
	
	protected File getOutputSolFile(LogicTreeBranch<?> branch) {
		String dirName = branch.buildFileName();
		File runDir = new File(outputDir, dirName);
		Preconditions.checkState(runDir.exists() || runDir.mkdir());
		
		return new File(runDir, "solution.zip");
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		int threads = getNumThreads();
		if (threads > 1) {
			ExecutorService exec = Executors.newFixedThreadPool(threads);
			
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
			LogicTreeBranch<LogicTreeNode> branch = outputTree.getBranch(index);
			
			debug("index "+index+" is branch: "+branch);
			
			File destSolFile = getOutputSolFile(branch);
			
			if (destSolFile.exists()) {
				debug(destSolFile.getAbsolutePath()+" exists, testing loading...");
				try {
					FaultSystemSolution.load(destSolFile);
					debug("skipping "+index+" (already done)");
					return;
				} catch (Exception e) {
					debug("Failed to load, re-inverting: "+e.getMessage());
				}
			}
			
			memoryDebug("Loading strict cutoff solutions for "+index);
			
			try {
				List<FaultSystemSolution> sols = new ArrayList<>();
				List<LogicTreeBranch<?>> branches = new ArrayList<>();
				List<MaxJumpDistModels> maxDists = new ArrayList<>();
				List<CSVFile<String>> passthroughCSVs = new ArrayList<>();
				
				ClusterRuptures cRups = null;
				ClusterConnectionStrategy connStrat = null;
				SectionDistanceAzimuthCalculator distAzCalc = null;
				
				File tempDir = Files.createTempDir();
				
				for (MaxJumpDistModels maxDist : MaxJumpDistModels.values()) {
					if ((float)maxDist.getMaxDist() < (float)minDistance)
						continue;
					
					ImmutableList<LogicTreeLevel<? extends LogicTreeNode>> levels = inputTree.getLevels();
					LogicTreeBranch<LogicTreeNode> inputBranch = new LogicTreeBranch<>(levels);
					for (LogicTreeNode node : branch) {
						for (int i=0; i<levels.size(); i++) {
							if (levels.get(i).isMember(node)) {
								inputBranch.setValue(i, node);
								break;
							}
						}
					}
					inputBranch.setValue(maxDist);
					Preconditions.checkState(inputBranch.isFullySpecified());
					File inputSolFile = getInputSolFile(inputBranch);
					Preconditions.checkState(inputSolFile.exists());
					FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
					
					if (cRups == null) {
						FaultSystemRupSet rupSet = sol.getRupSet();
						cRups = rupSet.requireModule(ClusterRuptures.class);
						PlausibilityConfiguration plausibility = rupSet.requireModule(PlausibilityConfiguration.class);
						connStrat = plausibility.getConnectionStrategy();
						distAzCalc = plausibility.getDistAzCalc();
					}
					
					SegmentationCalculator calc = new SegmentationCalculator(
							sol, cRups.getAll(), connStrat, distAzCalc, new double[] {0d});
					String prefix = maxDist.name();
					calc.plotDistDependComparison(tempDir, maxDist.name(), false, RateCombiner.MIN);
					File csvFile = new File(tempDir, prefix+"_supra_seis.csv");
					
					maxDists.add(maxDist);
					sols.add(sol);
					branches.add(inputBranch);
					passthroughCSVs.add(CSVFile.readFile(csvFile, true));
				};
				
				FileUtils.deleteRecursive(tempDir);
				
				SegmentationModelBranchNode segModel = branch.requireValue(
						SegmentationModelBranchNode.class);
				double r0;
				if (segModel instanceof ShawSegmentationModels) {
					Shaw07JumpDistProb shawProb = (Shaw07JumpDistProb)segModel.getModel(sols.get(0).getRupSet(), null);
					r0 = shawProb.getR0();
				} else if (segModel instanceof NSHM23_SegmentationModels) {
					r0 = ((NSHM23_SegmentationModels)segModel).getShawR0();
				} else {
					// just try it
					JumpProbabilityCalc model = segModel.getModel(sols.get(0).getRupSet(), null);
					if (!(model instanceof Shaw07JumpDistProb))
						model = segModel.getModel(null, null);
					Preconditions.checkState(model instanceof Shaw07JumpDistProb,
							"Couldn't get a Shaw07 from seg model: %s", segModel);
					r0 = ((Shaw07JumpDistProb)model).getR0();
				}
				
				List<Double> weights = new ArrayList<>();
				Map<LogicTreeBranch<?>, Double> weightsMap = new HashMap<>();
				synchronized (MaxJumpDistModels.class) {
					// recalculate weights
					MaxJumpDistModels.invertForWeights(maxDists.toArray(new MaxJumpDistModels[0]), passthroughCSVs, r0);
					for (int i=0; i<maxDists.size(); i++) {
						MaxJumpDistModels maxDist = maxDists.get(i);
						double weight = maxDist.getNodeWeight(null);
						weights.add(weight);
						weightsMap.put(branches.get(i), weight);
					}
				}
				
//				FaultSystemSolution avgSol = AverageSolutionCreator.buildAverage(
//						sols.toArray(new FaultSystemSolution[0]), Doubles.toArray(weights));
				
				// need to do BA sol, some things like slip rates can vary across max dists (due to sub-seis reduction)
				BranchAverageSolutionCreator baCreate = new BranchAverageSolutionCreator(
						new BranchWeightProvider.HardcodedWeights(weightsMap));
				for (int i=0; i<sols.size(); i++) {
					if (weights.get(i) == 0d) {
						System.out.println("Warning: inversion supplied zero weight for maxDist="+maxDists.get(i)+" on branch "+branches.get(i));
						continue;
					}
					baCreate.addSolution(sols.get(i), branches.get(i));
				}
				
				FaultSystemSolution avgSol = baCreate.build();
				
				String info = "Post-processed segmentation branch from "+sols.size()+" strict cutoff runs. Weights:\n\n";
				for (int i=0; i<sols.size(); i++)
					info += "\t"+maxDists.get(i)+": "+weights.get(i);
				avgSol.addModule(new InfoModule(info));
				
				avgSol.write(destSolFile);
				avgSol = null;
				sols = null;
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			memoryDebug("DONE "+index);
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		if (rank == 0) {
			memoryDebug("waiting for any post batch hook operations to finish");
			((AsyncLogicTreeWriter)postBatchHook).shutdown();
			memoryDebug("post batch hook done");
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		ops.addRequiredOption("ilt", "input-logic-tree", true, "Path to input logic tree JSON file");
		ops.addRequiredOption("olt", "output-logic-tree", true, "Path to output logic tree JSON file");
		ops.addRequiredOption("id", "input-dir", true, "Path to input directory");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addRequiredOption("ifc", "inversion-factory", true, "Inversion configuration factory classname");
		ops.addOption("md", "min-distance", true, "Minimum jump distance to include (default uses all)");
		
		return ops;
	}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_LogicTreeInversionRunner.class);
			
			MPJ_StrictSegLogicTreeTranslation driver = new MPJ_StrictSegLogicTreeTranslation(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
