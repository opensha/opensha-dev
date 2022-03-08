package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree.FileBuilder;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_EpistemicLogicTreeNode;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_LogicTreeNode;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class GMM_LogicTreeBranchAttach {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/");
		
		File mainDir = new File(invDir, "2021_11_30-u3_branches-orig_calcs-5h");
		int samples = 0;
		boolean writeHazardScript = true;
		IncludeBackgroundOption hazBgkOp = IncludeBackgroundOption.INCLUDE;
		File resultsFile = new File(mainDir, "results.zip");
		
		String destName = mainDir.getName()+"-with_gmm_branches";
		if (samples > 0)
			destName += "-"+samples+"_samples";
		File destDir = new File(invDir, destName);
		
		List<LogicTreeLevel<?>> addLevels = new ArrayList<>();
		addLevels.add(LogicTreeLevel.forEnum(NGAW2_LogicTreeNode.class, "NGA-West2 GMM", "NGAW2 GMM"));
		addLevels.add(LogicTreeLevel.forEnum(NGAW2_EpistemicLogicTreeNode.class,
				"NGA-West2 Additional Epistemic Uncertainty Option", "NGAW2 Epi"));
		
		File resultsModFile = new File(destDir, resultsFile.getName());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		Preconditions.checkState(destDir.exists() || destDir.mkdir());
		
		LogicTree<?> origTree = slt.getLogicTree();
		if (samples > 0)
			origTree = origTree.sample(samples, true, new Random((long)origTree.size()*(long)samples));
		
		FileBuilder builder = new SolutionLogicTree.FileBuilder(slt.getProcessor(), resultsModFile);
		
		for (LogicTreeBranch<?> branch : origTree) {
			System.out.println("Processing branch: "+branch);
			
			FaultSystemSolution sol = slt.forBranch(branch, false);
			
			for (LogicTreeBranch<LogicTreeNode> subBranch : buildSubBranches(branch, addLevels))
				builder.solution(sol, subBranch);
		}
		
		SolutionLogicTree finalTree = builder.build();
		
		System.out.println("Built "+finalTree.getLogicTree().size()+" branches");
		
		if (writeHazardScript) {
			File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
			int remoteToalThreads = 20;
			int remoteTotalMemGB = 53;
			String queue = "scec";
			int nodes = 38;
			
//			int minsEach = 10;
//			double gridSpacing = 1d;
			
			int minsEach = 45;
			double gridSpacing = 0.1d;
			
//			JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//					USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
			JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
					USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
			BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
			
			File remoteJobDir = new File(remoteMainDir, destDir.getName());
			File remoteResultsDir = new File(remoteJobDir, "results_hazard");
			File remoteResultsFile = new File(remoteJobDir, resultsModFile.getName());
			
			List<File> classpath = new ArrayList<>();
			classpath.add(new File(remoteJobDir, "opensha-dev-all.jar"));
			
			mpjWrite.setClasspath(classpath);
			if (mpjWrite instanceof MPJExpressShellScriptWriter)
				((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			else if (mpjWrite instanceof FastMPJShellScriptWriter)
				((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			
			String argz = "--input-file "+remoteResultsFile.getAbsolutePath();
			argz += " --output-dir "+remoteResultsDir.getAbsolutePath();
			if (hazBgkOp != null) {
				argz += " --gridded-seis "+hazBgkOp.name();
				argz += " --max-distance 200"; // reduce max distance when we have gridded seis
			}
			argz += " --grid-spacing "+gridSpacing;
			argz += " --periods 0,1";
			argz += " "+MPJTaskCalculator.argumentBuilder().exactDispatch(1).threads(remoteToalThreads).build();
			List<String> script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
			
			int nodeRounds = (int)Math.ceil((double)finalTree.getLogicTree().size()/(double)(nodes));
			int mins = minsEach*nodeRounds + 60;
			nodes = Integer.min(40, nodes);
			pbsWrite.writeScript(new File(destDir, "batch_hazard.slurm"), script, mins, nodes, remoteToalThreads, queue);
		}
	}
	
	private static List<LogicTreeBranch<LogicTreeNode>> buildSubBranches(LogicTreeBranch<?> curBranch,
			List<LogicTreeLevel<?>> addLevels) {
		List<LogicTreeBranch<LogicTreeNode>> ret = new ArrayList<>();
		
		Preconditions.checkState(!addLevels.isEmpty());
		LogicTreeLevel<?> curLevel = addLevels.get(0);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> modLevels = new ArrayList<>();
		modLevels.addAll(curBranch.getLevels());
		modLevels.add(curLevel);
		
		List<LogicTreeLevel<?>> nextLevels = addLevels.size() == 1 ? null : addLevels.subList(1, addLevels.size());
		
		for (LogicTreeNode choice : curLevel.getNodes()) {
			LogicTreeBranch<LogicTreeNode> modBranch = new LogicTreeBranch<>(modLevels);
			
			for (int i=0; i<curBranch.size(); i++)
				modBranch.setValue(i, curBranch.getValue(i));
			
			modBranch.setValue(curBranch.size(), choice);
			
			if (nextLevels == null)
				ret.add(modBranch);
			else
				ret.addAll(buildSubBranches(modBranch, nextLevels));
		}
		
		return ret;
	}

}
