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

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
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
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2021_12_16-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-ZENGBB-Shaw09Mod-DsrUni-TotNuclRate-SubB1");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2022_01_16-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-U3RupSet-SubB1-5000ip");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_U3RupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2022_01_19-nshm23_branches-reweighted_even_fit-CoulombRupSet-DsrUni-SubB1-ShawR0_3-5000ip");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2022_01_25-nshm23_u3_hybrid_branches-max_dist-CoulombRupSet-U3_ZENG-Shaw09Mod-DsrUni-SubB1-2000ip");
//		File resultsFile = new File(mainDir, "results.zip");
//		File fullBAFile = new File(mainDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
//		File mainDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//		File mainDir = new File(invDir, "2022_02_08-nshm23_u3_hybrid_branches-seg_bin_dist_capped_distr-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//		File mainDir = new File(invDir, "2022_01_28-nshm23_u3_hybrid_branches-max_dist-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip");
//		File mainDir = new File(invDir, "2022_02_27-nshm23_u3_hybrid_branches-strict_cutoff_seg-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-2000ip");
//		File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-strict_cutoff_seg-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1");
//		File mainDir = new File(invDir, "2022_05_09-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
//		File mainDir = new File(invDir, "2022_06_01-nshm23_u3_hybrid_branches-cluster_specific_inversion-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ShawR0_3-Shift2km-ThreshAvg");
//		File mainDir = new File(invDir, "2022_05_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg");
//		File mainDir = new File(invDir, "2022_07_25-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
		File mainDir = new File(invDir, "2022_07_28-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
		File resultsFile = new File(mainDir, "results.zip");
		File fullBAFile = new File(mainDir, "results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
		
//		HashSet<Class<? extends LogicTreeNode>> restrictBAClasses = null;
		HashSet<Class<? extends LogicTreeNode>> restrictBAClasses = new HashSet<>();
		restrictBAClasses.add(NSHM23_DeformationModels.class);
//		restrictBAClasses.add(NSHM23_ScalingRelationships.class);
		restrictBAClasses.add(NSHM23_SegmentationModels.class);
//		restrictBAClasses.add(MaxJumpDistModels.class);
//		restrictBAClasses.add(SegmentationModels.class);
//		restrictBAClasses.add(SubSectConstraintModels.class);
//		restrictBAClasses.add(SupraSeisBValues.class);
//		restrictBAClasses.add(RupsThroughCreepingSect.class);
		
		LogicTreeNode[] restrictNodes = {
//				FaultModels.FM3_1
		};
		
		int totThreads = FaultSysTools.defaultNumThreads();
		
		int asyncThreads = Integer.min(8, totThreads);
		boolean replot = true;
		
		HazardMapPlot.SPACING_DEFAULT = 0.2;
		
		File outputDir = new File(mainDir, "node_branch_averaged");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		FaultSystemSolution fullBA = fullBAFile == null ? null : FaultSystemSolution.load(fullBAFile);
//		PlotLevel plt = PlotLevel.FULL;
		PlotLevel plt = PlotLevel.DEFAULT;
		boolean compWithLoaded = false;
		boolean skipSectBySect = true;
		
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
//					creator.skipModule(SupraSeisBValInversionTargetMFDs.class);
					creator.skipModule(InversionTargetMFDs.class);
					nodeBACreators.put(node, creator);
				}
			}
		}
		
		int maxTasks = asyncThreads * 2;
		ExecutorService exec = new ThreadPoolExecutor(asyncThreads, asyncThreads,
                0L, TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<Runnable>(maxTasks), new ThreadPoolExecutor.CallerRunsPolicy());
		
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
						ReportPageGen pageGen = new ReportPageGen(meta, reportDir, ReportPageGen.getDefaultSolutionPlots(plt));
						if (skipSectBySect)
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

}
