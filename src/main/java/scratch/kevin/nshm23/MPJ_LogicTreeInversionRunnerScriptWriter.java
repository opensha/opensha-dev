package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.inversion.ClusterSpecificInversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeInversionRunner;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.DistDependSegShift;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupsThroughCreepingSect;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSectConstraintModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSeisMoRateReductions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.U3_UncertAddDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_LogicTreeBranch;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.inversion.U3InversionConfigFactory;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;

public class MPJ_LogicTreeInversionRunnerScriptWriter {
	
	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		List<String> extraArgs = new ArrayList<>();
		
		boolean strictSeg = false;
		double segTransMaxDist = 3d;
		boolean hazardGridded = false;
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteInversionsPerBundle = 1;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
		double itersPerSec = 200000;
		int runsPerBranch = 1;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
//		File remoteMainDir = new File("/work/00950/kevinm/stampede2/nshm23/batch_inversions");
//		int remoteTotalThreads = 48;
//		int remoteInversionsPerBundle = 3;
//		int remoteTotalMemGB = 100;
//		String queue = "skx-normal";
//		int nodes = 100;
//		double itersPerSec = 300000;
//		int runsPerBranch = 1;
////		String queue = "skx-dev";
////		int nodes = 4;
//		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
//				StampedeScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, StampedeScriptWriter.FMPJ_HOME);
//		BatchScriptWriter pbsWrite = new StampedeScriptWriter();

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
//		String dirName = "2022_02_27";
		
		/*
		 * UCERF3 logic tree
		 */
//		List<LogicTreeLevel<? extends U3LogicTreeBranchNode<?>>> levels = U3LogicTreeBranch.getLogicTreeLevels();
//		LogicTree<U3LogicTreeBranchNode<?>> logicTree = LogicTree.buildExhaustive(levels, true);
//		dirName += "-u3_branches";
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.class;
////		int avgNumRups = 250000;
//		
//		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.ScaleLowerDepth1p3.class;
//		int avgNumRups = 250000;
//		dirName += "-scaleLowerDepth1.3";
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.OriginalCalcParamsNewAvgConverged.class;
////		dirName += "-orig_calc_params-new_avg-converged";
////		int avgNumRups = 250000;
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.OriginalCalcParamsNewAvgNoWaterLevel.class;
////		dirName += "-orig_calc_params-new_avg-converged-noWL";
////		int avgNumRups = 250000;
//		
////		dirName += "-new_perturb";
////		extraArgs.add("--perturb "+GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE.name());
//		
////		dirName += "-try_zero_often";
////		extraArgs.add("--non-negativity "+NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN.name());
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.ThinnedRupSet.class;
////		int avgNumRups = 150000;
////		dirName += "-thinned_0.1";
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.OriginalCalcParams.class;
////		dirName += "-orig_calc_params";
////		int avgNumRups = 250000;
////		remoteInversionsPerBundle = 4;
////		runsPerBranch = 10;
////		String completionArg = null;
////		int invMins = 30;
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.OriginalCalcParamsNewAvg.class;
////		dirName += "-orig_calc_params-new_avg";
////		int avgNumRups = 250000;
////		String completionArg = null;
////		int invMins = 30;
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.NoPaleoParkfieldSingleReg.class;
////		dirName += "-no_paleo-no_parkfield-single_mfd_reg";
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.CoulombRupSet.class;
////		dirName += "-coulomb";
////		int avgNumRups = 325000;
//		
////		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.CoulombBilateralRupSet.class;
////		dirName += "-coulomb-bilateral";
////		int avgNumRups = 500000;
//		
////		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1, DeformationModels.GEOLOGIC,
////				ScalingRelationships.SHAW_2009_MOD, TotalMag5Rate.RATE_7p9 };
////		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1, DeformationModels.ZENGBB,
////				ScalingRelationships.SHAW_2009_MOD };
//		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1 };
////		U3LogicTreeBranchNode<?>[] required = {  };
//		Class<? extends LogicTreeNode> sortBy = null;
		/*
		 * END UCERF3 logic tree
		 */

		/*
		 * NSHM23 logic tree
		 */
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_U3_HybridLogicTreeBranch.levels;
		dirName += "-nshm23_u3_hybrid_branches";
		double avgNumRups = 300000;
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levels;
//		dirName += "-nshm23_branches";
//		double avgNumRups = 600000;
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM18_LogicTreeBranch.levels;
//		dirName += "-nshm18_branches";
		
		levels = new ArrayList<>(levels);
		for (int i=levels.size(); --i>=0;)
			if (levels.get(i).getType().isAssignableFrom(SegmentationModels.class)
					|| levels.get(i).getType().isAssignableFrom(SegmentationMFD_Adjustment.class)
					|| levels.get(i).getType().isAssignableFrom(DistDependSegShift.class))
				levels.remove(i);
		dirName += "-no_seg";
//		levels.add(NSHM23_LogicTreeBranch.MAX_DIST);
//		dirName += "-strict_cutoff_seg"; strictSeg = true;
		
//		dirName += "-reweight_seg_2_3_4";
		
//		levels = new ArrayList<>(levels);
//		int origSize = levels.size();
//		for (int i=levels.size(); --i>=0;)
//			if (levels.get(i).getType().isAssignableFrom(ScalingRelationships.class))
//				levels.remove(i);
//		Preconditions.checkState(levels.size() < origSize);
//		levels.add(LogicTreeLevel.forEnum(HardcodedTestScaleRels.class, "Hardocded Test Scaling Rels", "TestScales"));
//		dirName += "-test_scale_rels";
			
		LogicTree<LogicTreeNode> logicTree = LogicTree.buildExhaustive(levels, true);
		
		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.class;
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.MFDUncert0p1.class;
//		dirName += "-mfd_uncert_0p1";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.FullSysInv.class;
//		dirName += "-full_sys_inv";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.ClusterSpecific.class;
//		dirName += "-cluster_specific_inversion";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.SegWeight100.class;
//		dirName += "-seg_weight_100";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.HardcodedPrevWeightAdjust.class;
//		dirName += "-no_reweight_use_prev";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.HardcodedOrigWeights.class;
//		dirName += "-no_reweight_use_orig";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.HardcodedPrevAvgWeights.class;
//		dirName += "-no_reweight_use_prev_avg";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.NoPaleoParkfield.class;
//		dirName += "-no_paleo_parkfield";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.NoMFDScaleAdjust.class;
//		dirName += "-no_scale_adj_mfds";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.NoIncompatibleDataAdjust.class;
//		dirName += "-no_mfd_sigma_data_adj";
		
//		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.ScaleLowerDepth1p3.class;
//		dirName += "-scaleLowerDepth1.3";
		
		LogicTreeNode[] required = {
				FaultModels.FM3_1,
//				NSHM18_FaultModels.NSHM18_WUS_NoCA,
//				NSHM23_FaultModels.NSHM23_v1p4,

//				RupturePlausibilityModels.COULOMB,
				RupturePlausibilityModels.COULOMB_5km,
//				RupturePlausibilityModels.AZIMUTHAL,
//				RupturePlausibilityModels.SEGMENTED,
//				RupturePlausibilityModels.UCERF3,
//				RupturePlausibilityModels.UCERF3_REDUCED,
				
//				U3_UncertAddDeformationModels.U3_ZENG,
				
//				ScalingRelationships.SHAW_2009_MOD,
				
				SlipAlongRuptureModels.UNIFORM,
//				SlipAlongRuptureModels.TAPERED,
				
				SubSectConstraintModels.TOT_NUCL_RATE,
//				SubSectConstraintModels.NUCL_MFD,
				
				SubSeisMoRateReductions.SUB_B_1,
//				SubSeisMoRateReductions.NONE,
//				SubSeisMoRateReductions.SYSTEM_AVG,
//				SubSeisMoRateReductions.SYSTEM_AVG_SUB_B_1,
				
//				SupraSeisBValues.B_0p8,
				
//				SegmentationModels.SHAW_R0_3,
				
//				DistDependSegShift.NONE,
//				DistDependSegShift.ONE_KM,
//				DistDependSegShift.TWO_KM,
//				DistDependSegShift.THREE_KM,
				
//				SegmentationMFD_Adjustment.NONE,
//				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG,
//				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG,
//				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG_ITERATIVE,
//				SegmentationMFD_Adjustment.CAPPED_REDIST,
//				SegmentationMFD_Adjustment.CAPPED_REDIST_SELF_CONTAINED,
//				SegmentationMFD_Adjustment.GREEDY,
//				SegmentationMFD_Adjustment.GREEDY_SELF_CONTAINED,
//				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG_MATCH_STRICT,
				
				RupsThroughCreepingSect.INCLUDE,
//				RupsThroughCreepingSect.EXCLUDE,
				};
//		LogicTreeNode[] required = { FaultModels.FM3_1, SubSeisMoRateReductionNode.SYSTEM_AVG };
//		LogicTreeNode[] required = { FaultModels.FM3_1, SubSeisMoRateReductionNode.FAULT_SPECIFIC };
		Class<? extends LogicTreeNode> sortBy = SubSectConstraintModels.class;
		/*
		 * END NSHM23 logic tree
		 */
		
		int rounds = 2000;
		String completionArg = null;
		
////		int rounds = 5000;
//		int rounds = 10000;
//		String completionArg = rounds+"ip";
		
		double numIters = avgNumRups*rounds;
		double invSecs = numIters/itersPerSec;
		int invMins = (int)(invSecs/60d + 0.5);
		if (ClusterSpecificInversionConfigurationFactory.class.isAssignableFrom(factoryClass))
			invMins *= 2;
		System.out.println("Estimate "+invMins+" minues per inversion");
		
//		String completionArg = "1m"; int invMins = 1;
//		String completionArg = "10m"; int invMins = 10;
//		String completionArg = "30m"; int invMins = 30;
//		String completionArg = "2h"; int invMins = 2*60;
//		String completionArg = "5h"; int invMins = 5*60;
//		String completionArg = null; int invMins = defaultInvMins;
		
//		int numSamples = nodes*5;
//		int numSamples = nodes*4;
//		long randSeed = System.currentTimeMillis();
		long randSeed = 12345678l;
		int numSamples = 0;
		
		if (required != null && required.length > 0) {
			for (LogicTreeNode node : required)
				dirName += "-"+node.getFilePrefix();
			logicTree = logicTree.matchingAll(required);
		}
		
		if (numSamples > 0) {
			if (numSamples < logicTree.size()) {
				System.out.println("Reducing tree of size "+logicTree.size()+" to "+numSamples+" samples");
				dirName += "-"+numSamples+"_samples";
				logicTree = logicTree.sample(numSamples, false, new Random(randSeed));
			} else {
				System.out.println("Won't sample logic tree, as tree has "+logicTree.size()+" values, which is fewer "
						+ "than the specified "+numSamples+" samples.");
			}
		}
		
		if (sortBy != null) {
			Comparator<LogicTreeBranch<?>> groupingComparator = new Comparator<LogicTreeBranch<?>>() {

				@SuppressWarnings("unchecked")
				@Override
				public int compare(LogicTreeBranch<?> o1, LogicTreeBranch<?> o2) {
					LogicTreeNode v1 = o1.getValue(sortBy);
					LogicTreeNode v2 = o2.getValue(sortBy);
					if (v1.equals(v2))
						return ((LogicTreeBranch<LogicTreeNode>)o1).compareTo((LogicTreeBranch<LogicTreeNode>)o2);
					
					return v1.getShortName().compareTo(v2.getShortName());
				}
			};
			logicTree = logicTree.sorted(groupingComparator);
		}
		
		System.out.println("Built "+logicTree.size()+" branches");
		Preconditions.checkState(logicTree.size() > 0, "No matching branches");
		
		int origNodes = nodes;
		nodes = Integer.min(nodes, logicTree.size());
		
		int numCalcs = logicTree.size()*runsPerBranch;
		int nodeRounds = (int)Math.ceil((double)numCalcs/(double)(nodes*remoteInversionsPerBundle));
		double calcNodes = Math.ceil((double)numCalcs/(double)(nodeRounds*remoteInversionsPerBundle));
		System.out.println("Implies "+(float)calcNodes+" nodes for "+nodeRounds+" rounds");
		nodes = Integer.min(nodes, (int)calcNodes);
		if (origNodes != nodes)
			System.out.println("Adjusted "+origNodes+" to "+nodes+" nodes to evenly divide "+numCalcs+" calcs");
		
		if (completionArg != null)
			dirName += "-"+completionArg;
		
		if (hazardGridded)
			dirName += "-gridded";
		
		System.out.println("Directory name: "+dirName);
		
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		File remoteDir = new File(remoteMainDir, dirName);
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		
		File localLogicTree = new File(localDir, "logic_tree.json");
		logicTree.write(localLogicTree);
		File remoteLogicTree = new File(remoteDir, localLogicTree.getName());
		
		mpjWrite.setClasspath(classpath);
		if (mpjWrite instanceof MPJExpressShellScriptWriter)
			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		else if (mpjWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		int annealingThreads = remoteTotalThreads/remoteInversionsPerBundle;
		
		File resultsDir = new File(remoteDir, "results");
		String argz = "--logic-tree "+remoteLogicTree.getAbsolutePath();
		argz += " --output-dir "+resultsDir.getAbsolutePath();
		argz += " --inversion-factory '"+factoryClass.getName()+"'"; // surround in single quotes to escape $'s
		argz += " --annealing-threads "+annealingThreads;
		argz += " --cache-dir "+new File(remoteMainDir, "cache").getAbsolutePath();
		if (remoteInversionsPerBundle > 0)
			argz += " --runs-per-bundle "+remoteInversionsPerBundle;
		if (completionArg != null)
			argz += " --completion "+completionArg;
		if (runsPerBranch > 1)
			argz += " --runs-per-branch "+runsPerBranch;
		for (String arg : extraArgs)
			argz += " "+arg;
		argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(remoteInversionsPerBundle).build();
		List<String> script = mpjWrite.buildScript(MPJ_LogicTreeInversionRunner.class.getName(), argz);
		
		int mins = (int)(nodeRounds*invMins);
		mins += Integer.max(60, invMins);
		mins += (nodeRounds-1)*10; // add a little extra for overhead associated with each round
		System.out.println("Total job time: "+mins+" mins = "+(float)((double)mins/60d)+" hours");
		
		pbsWrite.writeScript(new File(localDir, "batch_inversion.slurm"), script, mins, nodes, remoteTotalThreads, queue);
		
		// now write hazard script
		argz = "--input-file "+new File(resultsDir.getAbsolutePath()+".zip");
		argz += " --output-dir "+resultsDir.getAbsolutePath();
		if (hazardGridded)
			argz += " --gridded-seis INCLUDE";
		argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(1).threads(remoteTotalThreads).build();
		script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
		
		mins = (int)60*10;
		nodes = Integer.min(40, nodes);
		if (queue != null && queue.equals("scec"))
			// run hazard in the high priority queue
			queue = "scec_hiprio";
		pbsWrite.writeScript(new File(localDir, "batch_hazard.slurm"), script, mins, nodes, remoteTotalThreads, queue);
		
		if (strictSeg) {
			// write strict segmentation branch translation job
			String modDirName = localDir.getName()+"-branch-translated";
			if (segTransMaxDist > 0d)
				modDirName += "-min"+new DecimalFormat("0.#").format(segTransMaxDist)+"km";
			File modLocalDir = new File(localMainDir, modDirName);
			Preconditions.checkState(modLocalDir.exists() || modLocalDir.mkdir());
			File modRemoteDir = new File(remoteMainDir, modLocalDir.getName());
			
			classpath = new ArrayList<>();
			classpath.add(new File(modRemoteDir, "opensha-dev-all.jar"));
			mpjWrite.setClasspath(classpath);
			
			List<LogicTreeLevel<? extends LogicTreeNode>> modLevels = new ArrayList<>();
			for (LogicTreeLevel<? extends LogicTreeNode> level : levels)
				if (!level.matchesType(MaxJumpDistModels.class))
					modLevels.add(level);
			modLevels.add(NSHM23_LogicTreeBranch.SEG);
			
			List<LogicTreeBranch<LogicTreeNode>> modBranches = new ArrayList<>();
			
			HashSet<String> branchNameHash = new HashSet<>();
			for (LogicTreeBranch<?> branch : logicTree) {
				LogicTreeBranch<LogicTreeNode> modBranch = new LogicTreeBranch<LogicTreeNode>(modLevels);
				for (LogicTreeNode node : branch)
					if (!(node instanceof MaxJumpDistModels))
						modBranch.setValue(node);
				String baseName = modBranch.toString();
				if (branchNameHash.contains(baseName))
					continue;
				branchNameHash.add(baseName);
				for (SegmentationModels segModel : SegmentationModels.values()) {
					if (segModel.getNodeWeight(modBranch) > 0d) {
						LogicTreeBranch<LogicTreeNode> fullBranch = modBranch.copy();
						fullBranch.setValue(segModel);
						modBranches.add(fullBranch);
					}
				}
			}
			
			System.out.println("Translating "+logicTree.size()+" branches to "+modBranches.size());
			LogicTree<LogicTreeNode> modTree = LogicTree.fromExisting(modLevels, modBranches);
			File modLocalLogicTree = new File(modLocalDir, remoteLogicTree.getName());
			modTree.write(modLocalLogicTree);
			File modRemoteLogicTree = new File(modRemoteDir, modLocalLogicTree.getName());
			
			File modOutputDir = new File(modRemoteDir, "results");
			argz = "--input-dir "+resultsDir.getAbsolutePath();
			argz += " --input-logic-tree "+remoteLogicTree.getAbsolutePath();
			argz += " --output-logic-tree "+modRemoteLogicTree.getAbsolutePath();
			argz += " --inversion-factory '"+factoryClass.getName()+"'"; // surround in single quotes to escape $'s
			argz += " --output-dir "+modOutputDir.getAbsolutePath();
			if (segTransMaxDist > 0d)
				argz += " --min-distance "+segTransMaxDist;
			argz += " "+MPJTaskCalculator.argumentBuilder().threads(Integer.min(10, remoteTotalThreads)).build();
			script = mpjWrite.buildScript(MPJ_StrictSegLogicTreeTranslation.class.getName(), argz);
			
			int transNodes = Integer.min(16, nodes);
			mins = (int)60*10;
			pbsWrite.writeScript(new File(modLocalDir, "batch_strict_branch_translate.slurm"), script, mins, transNodes,
					remoteTotalThreads, queue);
			
			// now write hazard script
			argz = "--input-file "+new File(modOutputDir.getAbsolutePath()+".zip");
			argz += " --output-dir "+modOutputDir.getAbsolutePath();
			argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(1).threads(remoteTotalThreads).build();
			script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
			
			mins = (int)60*10;
			nodes = Integer.min(40, nodes);
			pbsWrite.writeScript(new File(modLocalDir, "batch_hazard.slurm"), script, mins, nodes, remoteTotalThreads, queue);
		}
	}

}
