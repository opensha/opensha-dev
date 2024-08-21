package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.ExecutorUtils;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAverageableModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.reports.AbstractRupSetPlot;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.HazardMapPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupsThroughCreepingSect;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.ShawSegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSectConstraintModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class LogicTreeBranchAverageWriter {

	public static void main(String[] args) throws IOException {
		System.setProperty("java.awt.headless", "true");
		
		File fullBAFile = null;
		File outputDir;
		
		HashSet<Class<? extends LogicTreeNode>> restrictBAClasses = null;
		LogicTreeNode[] restrictNodes = null;
		
		List<Class<? extends BranchAverageableModule<?>>> skipModules = null;
		
		int totThreads = FaultSysTools.defaultNumThreads();
		int asyncThreads = -1;
		
		SolutionLogicTree slt;
		
		PlotLevel plt = ReportPageGen.PLOT_LEVEL_DEFAULT;
		boolean skipSectBySect = false;
		
		if (args.length == 0) {
			File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
			
//			File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2021_12_16-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-ZENGBB-Shaw09Mod-DsrUni-TotNuclRate-SubB1");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2022_01_16-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-U3RupSet-SubB1-5000ip");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_U3RupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2022_01_19-nshm23_branches-reweighted_even_fit-CoulombRupSet-DsrUni-SubB1-ShawR0_3-5000ip");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2022_01_25-nshm23_u3_hybrid_branches-max_dist-CoulombRupSet-U3_ZENG-Shaw09Mod-DsrUni-SubB1-2000ip");
//			File resultsFile = new File(mainDir, "results.zip");
//			File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
			
//			File mainDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//			File mainDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-seg_bin_dist_capped_distr-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//			File mainDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//			File mainDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//			File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-strict_cutoff_seg-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1");
//			File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//			File mainDir = new File(invDir, "2022_06_01-nshm23_u3_hybrid_branches-cluster_specific_inversion-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ShawR0_3-Shift2km-ThreshAvg");
//			File mainDir = new File(invDir, "2022_05_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg");
//			File mainDir = new File(invDir, "2022_07_25-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
//			File mainDir = new File(invDir, "2022_07_28-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
//			File mainDir = new File(invDir, "2022_07_29-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR");
//			File mainDir = new File(invDir, "2022_07_29-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR");
//			File mainDir = new File(invDir, "2022_08_08-nshm23_branches-wide_seg_branches-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-SubB1-ThreshAvgIterRelGR");
//			File mainDir = new File(invDir, "2022_08_05-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-SubB1-ThreshAvgIterRelGR");
			File mainDir = new File(invDir, "2022_08_22-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
			File resultsFile = new File(mainDir, "results_gridded_branches.zip");
//			File resultsFile = new File(mainDir, "results.zip");
			fullBAFile = new File(mainDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
			
//			HashSet<Class<? extends LogicTreeNode>> restrictBAClasses = null;
//			restrictBAClasses = new HashSet<>();
//			restrictBAClasses.add(NSHM23_DeformationModels.class);
//			restrictBAClasses.add(NSHM23_ScalingRelationships.class);
//			restrictBAClasses.add(NSHM23_SegmentationModels.class);
//			restrictBAClasses.add(MaxJumpDistModels.class);
//			restrictBAClasses.add(SegmentationModels.class);
//			restrictBAClasses.add(SubSectConstraintModels.class);
//			restrictBAClasses.add(SupraSeisBValues.class);
//			restrictBAClasses.add(RupsThroughCreepingSect.class);
			
			skipModules = new ArrayList<>();
			skipModules.add(InversionTargetMFDs.class);
			
//			restrictNodes = new LogicTreeNode[] {
//					FaultModels.FM3_1
//			};
			
//			plt = PlotLevel.FULL;
//			skipSectBySect = true;
			
			outputDir = new File(mainDir, "node_branch_averaged");
			
			slt = SolutionLogicTree.load(resultsFile);
			
			HazardMapPlot.SPACING_DEFAULT = 0.2;
		} else {
			CommandLine cmd = FaultSysTools.parseOptions(createOptions(), args, ReportPageGen.class);
			
			File inputFile = new File(cmd.getOptionValue("input-file"));
			Preconditions.checkArgument(inputFile.exists(), "Input file doesn't exist: %s", inputFile.getAbsolutePath());
			
			if (inputFile.isDirectory()) {
				Preconditions.checkArgument(cmd.hasOption("logic-tree"), "Must supply logic tree file if input-file is"
						+ " a results directory");
				File logicTreeFile = new File(cmd.getOptionValue("logic-tree"));
				Preconditions.checkArgument(logicTreeFile.exists(), "Logic tree file doesn't exist: %s",
						logicTreeFile.getAbsolutePath());
				LogicTree<?> tree = LogicTree.read(logicTreeFile);
				
				slt = new SolutionLogicTree.ResultsDirReader(inputFile, tree);
			} else {
				// it should be SolutionLogicTree zip file
				slt = SolutionLogicTree.load(inputFile);
			}
			
			if (cmd.hasOption("branch-averaged-file"))
				fullBAFile = new File(cmd.getOptionValue("branch-averaged-file"));
			
			outputDir = new File(cmd.getOptionValue("output-dir"));
			
			totThreads = FaultSysTools.getNumThreads(cmd);
			if (cmd.hasOption("async-threads"))
				asyncThreads = Integer.parseInt(cmd.getOptionValue("async-threads"));
			
			if (cmd.hasOption("plot-level"))
				plt = PlotLevel.valueOf(cmd.getOptionValue("plot-level").trim().toUpperCase());
			
			skipSectBySect = cmd.hasOption("skip-sect-by-sect");
			
			if (cmd.hasOption("level-class")) {
				try {
					Class<? extends LogicTreeNode> clazz = (Class<? extends LogicTreeNode>) Class.forName(cmd.getOptionValue("level-class"));
					restrictBAClasses = new HashSet<>();
					restrictBAClasses.add(clazz);
				} catch (ClassNotFoundException | ClassCastException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				
			}
		}
		
		if (asyncThreads < 1) {
			if (totThreads > 20)
				asyncThreads = 8;
			else
				asyncThreads = Integer.min(4, totThreads);
		}
		boolean replot = true;
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemSolution fullBA = fullBAFile == null ? null : FaultSystemSolution.load(fullBAFile);
		boolean compWithLoaded = false;
		
		LogicTree<?> tree = slt.getLogicTree();
		
		if (restrictNodes != null && restrictNodes.length > 0)
			tree = tree.matchingAll(restrictNodes);
		
//		tree = tree.matchingAll(SupraSeisBValues.B_0p0, DeformationModels.GEOLOGIC,
//				SubSectConstraintModels.TOT_NUCL_RATE, SegmentationModels.SHAW_R0_3);
//		System.out.println("Test sub-tree size: "+tree.size());
//		Preconditions.checkState(tree.size() > 0);
//		outputDir = new File(outputDir, "test_subset");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		plt = PlotLevel.LIGHT;
//		compWithLoaded = true;
		
		Map<LogicTreeLevel<?>, HashSet<LogicTreeNode>> levelNodes = new HashMap<>();
		List<? extends LogicTreeLevel<?>> levels = tree.getLevels();
		for (LogicTreeLevel<?> level : levels)
			levelNodes.put(level, new HashSet<>());
		
		Map<LogicTreeNode, LogicTreeLevel<?>> nodeLevels = new HashMap<>();
		
		for (LogicTreeBranch<?> branch : tree) {
			for (int i=0; i<levels.size(); i++) {
				LogicTreeNode node = branch.getValue(i);
				levelNodes.get(levels.get(i)).add(node);
				if (!nodeLevels.containsKey(node))
					nodeLevels.put(node, levels.get(i));
			}
		}
		
		Map<LogicTreeNode, BranchAverageSolutionCreator> nodeBACreators = new HashMap<>();
		
		for (LogicTreeLevel<?> level : levels) {
			HashSet<LogicTreeNode> nodes = levelNodes.get(level);
			if (nodes.size() > 1) {
				if (restrictBAClasses != null && !restrictBAClasses.contains(level.getType())) {
					System.out.println("Skipping "+nodes.size()+" BAs for "+level.getName());
					continue;
				}
				System.out.println("Building "+nodes.size()+" BAs for "+level.getName());
				for (LogicTreeNode node : nodes) {
					BranchAverageSolutionCreator creator = new BranchAverageSolutionCreator(tree.getWeightProvider());
					if (skipModules != null)
						for (Class<? extends BranchAverageableModule<?>> moduleClass : skipModules)
							creator.skipModule(moduleClass);
					nodeBACreators.put(node, creator);
				}
			}
		}
		
		// this dicates how far ahead we will read in solutions, needs to be > asyncThreads to make use of parallel
		// reading, but should be small enough soas to not use too much memory storing to-be-processed solutions
		int maxTasks = Integer.min(asyncThreads * 2, asyncThreads + 2);
		ExecutorService exec = ExecutorUtils.newBlockingThreadPool(asyncThreads, maxTasks);
		
		int count = 0;
		List<Future<?>> futures = new ArrayList<>();
		for (LogicTreeBranch<?> branch : tree) {
			System.out.println("Processing branch "+(count++)+"/"+tree.size()+": "+branch);
			FaultSystemSolution sol = slt.forBranch(branch);
			
			futures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					for (LogicTreeNode node : branch)
						if (nodeBACreators.containsKey(node))
							nodeBACreators.get(node).addSolution(sol, branch);
				}
			}));
		}
		
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		exec.shutdown();
		
		RupSetMetadata comparison = fullBA == null ? null : new RupSetMetadata("Full Tree", fullBA);
		
		int plotThreads =  Math.min(asyncThreads, Math.min(totThreads, nodeBACreators.size()));
		
		int threadsEach = Integer.max(1, totThreads / plotThreads);
		System.out.println("Building "+nodeBACreators.size()+" BAs with "+plotThreads
				+" plot threads ("+threadsEach+" threads per report)");
		
		exec = null;
		futures = null;
		
		if (plotThreads > 1) {
			exec = Executors.newFixedThreadPool(plotThreads);
			futures = new ArrayList<>();
		}
		
		List<AbstractRupSetPlot> plots = ReportPageGen.getDefaultSolutionPlots(plt);
		
		final boolean doSkipSectBySect = skipSectBySect;
		
		for (LogicTreeNode node : nodeBACreators.keySet()) {
			Runnable run = new Runnable() {
				
				@Override
				public void run() {
					try {
						System.out.println("Building BA for "+node);
						BranchAverageSolutionCreator creator = nodeBACreators.get(node);
						FaultSystemSolution ba = creator.build();
						
						LogicTreeLevel<?> level = nodeLevels.get(node);
						
						String prefix = level.getShortName().replaceAll("\\W+", "_")+"_"+node.getFilePrefix();
						
						File solFile = new File(outputDir, prefix+".zip");
						ba.write(solFile);
						
						File reportDir = new File(outputDir, prefix+"_report");
						Preconditions.checkState(reportDir.exists() || reportDir.mkdir());
						System.out.println("Writing report to "+reportDir.getAbsolutePath());
						RupSetMetadata primary = new RupSetMetadata(level.getShortName()+": "+node.getShortName(), ba);
						
						RupSetMetadata myComp = comparison;
						
						if (compWithLoaded)
							// test to compare serialized and regular
							myComp = new RupSetMetadata("Loaded", FaultSystemSolution.load(solFile));
						
						ReportMetadata meta = new ReportMetadata(primary, myComp);
						ReportPageGen pageGen = new ReportPageGen(meta, reportDir, plots);
						if (doSkipSectBySect)
							pageGen.skipSectBySect();
						pageGen.setReplot(replot);
						pageGen.setNumThreads(threadsEach);
						
						pageGen.generatePage();
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
			};
			if (plotThreads > 1)
				futures.add(exec.submit(run));
			else
				run.run();
		}
		
		if (futures != null) {
			System.out.println("Waiting on "+futures.size()+" futures");
			for (Future<?> future : futures) {
				try {
					future.get();
				} catch (Exception e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			
			exec.shutdown();
		}
	}
	
	public static Options createOptions() {
		Options ops = new Options();

		ops.addRequiredOption("if", "input-file", true, "Input SolutionLogicTree zip file or results directory. If a "
				+ "results directory is supplied, you must also specify --logic-tree <file>");
		ops.addOption("baf", "branch-averaged-file", true, "Optional path to a full branch-averaged solution, used "
				+ "for comparison in reports");
		ops.addOption("lt", "logic-tree", true, "Path to logic tree JSON file, required if a results directory is "
				+ "supplied with --input-file");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addOption(FaultSysTools.threadsOption());
		ops.addOption("at", "async-threads", true, "Maximum number of asynchronous load/process threads, lower to "
				+ "reduce memory usage or I/O load");
		ops.addOption("pl", "plot-level", true, "This determins which set of plots should be included. One of: "
						+FaultSysTools.enumOptions(PlotLevel.class)+". Default: "+ReportPageGen.PLOT_LEVEL_DEFAULT.name());
		ops.addOption("ssbs", "skip-sect-by-sect", false,
				"Flag to skip section-by-section plots, regardless of selected plot level");
		ops.addOption("lc", "level-class", true,
				"Flag to limit aggregation to only the given level, specified as a fully-qualified class name");
		
		return ops;
	}

}
