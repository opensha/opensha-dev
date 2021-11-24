package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;
import scratch.nshm23.logicTree.DraftNSHM23InvConfigFactory;
import scratch.nshm23.logicTree.DraftNSHM23LogicTreeBranch;

public class MPJ_LogicTreeInversionRunnerScriptWriter {
	
	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteToalThreads = 20;
		int remoteTotalMemGB = 55;
		BatchScriptWriter scriptWrite = new USC_CARC_ScriptWriter();
		String queue = "scec";

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		int nodes = 38;
		
//		LogicTree<U3LogicTreeBranchNode<?>> logicTree = LogicTree.buildExhaustive(
//				U3LogicTreeBranch.getLogicTreeLevels(), true);
//		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.class;
////		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1, DeformationModels.GEOLOGIC,
////				ScalingRelationships.SHAW_2009_MOD, TotalMag5Rate.RATE_7p9 };
//		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1 };
//		dirName += "-u3_branches";
		
		LogicTree<LogicTreeNode> logicTree = LogicTree.buildExhaustive(
				DraftNSHM23LogicTreeBranch.levels, true);
		Class<? extends InversionConfigurationFactory> factoryClass = DraftNSHM23InvConfigFactory.class;
//		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1, DeformationModels.GEOLOGIC,
//				ScalingRelationships.SHAW_2009_MOD, TotalMag5Rate.RATE_7p9 };
		U3LogicTreeBranchNode<?>[] required = { FaultModels.FM3_1 };
		dirName += "-nshm23_draft_branches";
		
		for (LogicTreeNode node : required)
			dirName += "-"+node.getFilePrefix();
		logicTree = logicTree.matchingAll(required);
		
		System.out.println("Built "+logicTree.size()+" branches");
		
		int origNodes = nodes;
		nodes = Integer.min(nodes, logicTree.size());
		
		int nodeRounds = (int)Math.ceil((double)logicTree.size()/(double)nodes);
		double calcNodes = Math.ceil((double)logicTree.size()/(double)nodeRounds);
		System.out.println("Implies "+(float)calcNodes+" nodes for "+nodeRounds+" rounds");
		nodes = Integer.min(nodes, (int)calcNodes);
		if (origNodes != nodes)
			System.out.println("Ajusted "+origNodes+"to "+nodes+" nodes to evenly divide "+logicTree.size()+" branches");
		
//		String completionArg = "10m"; int invMins = 10;
//		String completionArg = "5h"; int invMins = 5*60;
		String completionArg = null; int invMins = 5*60;
		
		if (completionArg != null)
			dirName += "-"+completionArg;
		
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		File remoteDir = new File(remoteMainDir, dirName);
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		
		File localLogicTree = new File(localDir, "logic_tree.json");
		logicTree.write(localLogicTree);
		File remoteLogicTree = new File(remoteDir, localLogicTree.getName());
		
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, classpath, USC_CARC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		String argz = "--logic-tree "+remoteLogicTree.getAbsolutePath();
		argz += " --output-dir "+new File(remoteDir, "results").getAbsolutePath();
		argz += " --inversion-factory "+factoryClass.getName();
		if (completionArg != null)
			argz += " --completion "+completionArg;
		argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(1).build();
		List<String> script = mpjWrite.buildScript(MPJ_LogicTreeInversionRunner.class.getName(), argz);
		
		int mins = (int)(nodeRounds*invMins + 60d);
		System.out.println("Total job time: "+mins+" mins = "+(float)((double)mins/60d)+" hours");
		
		pbsWrite.writeScript(new File(localDir, "batch_inversion.slurm"), script, mins, nodes, remoteToalThreads, queue);
	}

}
