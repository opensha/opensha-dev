package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeInversionRunner;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSectConstraintModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSeisMoRateReductions;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class InvFileDeadlockScriptWriter {

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		int numRuns = 100;
		int minsEach = 60;
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteInversionsPerBundle = 1;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 2;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		String testDirName = "2022_11_08-inv_file_write_deadlock_tests";
		
		File localBatchDir = new File(localMainDir, testDirName);
		Preconditions.checkState(localBatchDir.exists() || localBatchDir.mkdirs());
		File remoteBatchDir = new File(remoteMainDir, testDirName);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levelsOnFault;
		levels = new ArrayList<>(levels);
		levels.add(NSHM23_LogicTreeBranch.SINGLE_STATES);
		
		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.class;
		
		LogicTreeNode[] required = {
				// FAULT MODELS
				NSHM23_FaultModels.WUS_FM_v2,
				
//				// SINGLE STATE
				NSHM23_SingleStates.CA,

				// RUPTURE SETS
				RupturePlausibilityModels.AZIMUTHAL,
				
				// DEFORMATION MODELS
				NSHM23_DeformationModels.GEOLOGIC,
				
				// SCALING RELATIONSHIPS
				NSHM23_ScalingRelationships.AVERAGE,
				
				// SLIP ALONG RUPTURE
				
				// SUB-SECT CONSTRAINT
				SubSectConstraintModels.TOT_NUCL_RATE,
//				SubSectConstraintModels.NUCL_MFD,
				
				// SUB-SEIS MO REDUCTION
//				SubSeisMoRateReductions.SUB_B_1,
				SubSeisMoRateReductions.NONE,
//				SubSeisMoRateReductions.SYSTEM_AVG,
//				SubSeisMoRateReductions.SYSTEM_AVG_SUB_B_1,
				
				// SUPRA-SEIS-B
//				SupraSeisBValues.B_0p5,
				
				// PALEO UNCERT
				NSHM23_PaleoUncertainties.EVEN_FIT,
				
				// SEGMENTATION
//				SegmentationModels.SHAW_R0_3,
//				NSHM23_SegmentationModels.AVERAGE,
//				NSHM23_SegmentationModels.MID,
				NSHM23_SegmentationModels.CLASSIC,
				
				// SEG-SHIFT
//				DistDependSegShift.NONE,
//				DistDependSegShift.ONE_KM,
//				DistDependSegShift.TWO_KM,
//				DistDependSegShift.THREE_KM,
				
				// SEG ADJUSTMENT
//				SegmentationMFD_Adjustment.NONE,
//				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG,
//				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG_SINGLE_ITER,
				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG,
//				SegmentationMFD_Adjustment.CAPPED_REDIST,
//				SegmentationMFD_Adjustment.CAPPED_REDIST_SELF_CONTAINED,
//				SegmentationMFD_Adjustment.GREEDY,
//				SegmentationMFD_Adjustment.GREEDY_SELF_CONTAINED,
//				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG_MATCH_STRICT,
				
				// CREEPING SECTION
//				RupsThroughCreepingSect.INCLUDE,
//				RupsThroughCreepingSect.EXCLUDE,
				};
		/*
		 * END NSHM23 logic tree
		 */
		
		LogicTree<LogicTreeNode> logicTree = LogicTree.buildExhaustive(levels, true, new BranchWeightProvider.NodeWeightOverrides(required, 1d), required);

		System.out.println("Built "+logicTree.size()+" branches:");
		
		for (LogicTreeBranch<?> branch : logicTree)
			System.out.println("\t"+branch);
		
		File localLogicTree = new File(localBatchDir, "logic_tree.json");
		logicTree.write(localLogicTree);
		File remoteLogicTree = new File(remoteBatchDir, localLogicTree.getName());
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteBatchDir, "opensha-dev-all.jar"));
		mpjWrite.setClasspath(classpath);
		
		for (int run=0; run<numRuns; run++) {
			File localDir = new File(localBatchDir, "run_"+run);
			Preconditions.checkState(localDir.exists() || localDir.mkdir());
			File remoteDir = new File(remoteBatchDir, localDir.getName());
			
			File resultsDir = new File(remoteDir, "results");
			String argz = "--logic-tree "+remoteLogicTree.getAbsolutePath();
			argz += " --output-dir "+resultsDir.getAbsolutePath();
			argz += " --inversion-factory '"+factoryClass.getName()+"'"; // surround in single quotes to escape $'s
			argz += " --annealing-threads "+remoteTotalThreads;
			argz += " --cache-dir "+new File(remoteMainDir, "cache").getAbsolutePath();
			argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(remoteInversionsPerBundle).build();
			List<String> script = mpjWrite.buildScript(MPJ_LogicTreeInversionRunner.class.getName(), argz);
			
			pbsWrite.writeScript(new File(localDir, localDir.getName()+".slurm"), script, minsEach, nodes, remoteTotalThreads, queue);
		}
	}

}
