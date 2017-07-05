package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

import scratch.UCERF3.analysis.MPJDistributedCompoundFSSPlots;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter.RunSites;

public class BulkCompoundScriptWrite {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		List<LogicTreeBranchNode<?>> fmBranches =
				LogicTreePBSWriter.getNonZeroChoices(FaultModels.class, InversionModels.CHAR_CONSTRAINED);
		List<LogicTreeBranchNode<?>> dmBranches =
				LogicTreePBSWriter.getNonZeroChoices(DeformationModels.class, InversionModels.CHAR_CONSTRAINED);
//		List<LogicTreeBranchNode<?>> scaleBranches =
//				LogicTreePBSWriter.getNonZeroChoices(ScalingRelationships.class, InversionModels.CHAR_CONSTRAINED);
		List<LogicTreeBranchNode<?>> scaleBranches = Lists.newArrayList();
		scaleBranches.add(null);
		
//		File remoteMainDir = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_sub_plots");
//		
//		File compoundFile = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_combined/" +
//				"2013_01_14-stampede_3p2_production_runs_combined_COMPOUND_SOL.zip");
		
//		File remoteMainDir = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_sub_plots");
		
//		File compoundFile = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_combined/" +
//				"2013_01_14-stampede_3p2_production_runs_combined_COMPOUND_SOL.zip");
		
//		File remoteMainDir = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_dm_means");
//		
//		File compoundFile = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
//				"2013_01_14-stampede_3p2_production_runs_combined/" +
//				"2013_01_14-stampede_3p2_production_runs_combined_COMPOUND_SOL.zip");
		
//		File remoteMainDir = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/" +
//				"2013_05_10-ucerf3p3-production-10runs_fm_dm_scale_sub_plots");
//		
//		File compoundFile = new File("/work/00950/kevinm/ucerf3/inversion/" +
//				"2013_05_10-ucerf3p3-production-10runs/" +
//				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		
		File remoteMainDir = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
				"2013_05_10-ucerf3p3-production-10runs_fm_dm_sub_plots");
		
		File compoundFile = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/" +
				"2013_05_10-ucerf3p3-production-10runs/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF3/comp_bulk");
		if (!localMainDir.exists())
			localMainDir.mkdir();
		File writeDir = new File(localMainDir, remoteMainDir.getName());
		if (!writeDir.exists())
			writeDir.mkdir();
		
		RunSites site = RunSites.HPCC;
//		RunSites site = RunSites.STAMPEDE;
		int nodes = 30;
//		int bundleSize = 30; // TODO, must be >0
//		int jobMins = 6*60; // TODO
		int jobMins = 2*60; // TODO
		
		List<File> classpath = LogicTreePBSWriter.getClasspath(remoteMainDir, remoteMainDir);
		
		JavaShellScriptWriter mpjWrite;
		if (site.isFastMPJ())
			mpjWrite = new FastMPJShellScriptWriter(site.getJAVA_BIN(), site.getMaxHeapSizeMB(null),
					classpath, site.getMPJ_HOME());
		else
			mpjWrite = new MPJExpressShellScriptWriter(site.getJAVA_BIN(), site.getMaxHeapSizeMB(null),
					classpath, site.getMPJ_HOME());
		mpjWrite.setHeadless(true);
		
		BatchScriptWriter batchWrite = site.forBranch(null);
		if (batchWrite instanceof USC_HPCC_ScriptWriter)
			((USC_HPCC_ScriptWriter)batchWrite).setNodesAddition(null);
		
		for (LogicTreeBranchNode<?> fm : fmBranches) {
			for (LogicTreeBranchNode<?> dm : dmBranches) {
				for (LogicTreeBranchNode<?> scale : scaleBranches) {
					String dirName = fm.getShortName()+"_"+dm.getShortName();
					if (scale != null)
						dirName += "_"+scale.getShortName();
					File remoteJobDir = new File(remoteMainDir, dirName);
					String argVal = fm.getShortName()+"_,_"+dm.getShortName()+"_";
					if (scale != null)
						argVal += ",_"+scale.getShortName()+"_";
					
//					String argss = "--threads 2 --min-dispatch 2 --plot-paleo-faults --plot-parent-mfds" +
//							" --plot-slip-misfits --plot-participations --plot-mini-sect-ris --plot-ave-slips" +
//							" --name-grep "+argVal+" "+compoundFile.getAbsolutePath()+" "+remoteJobDir.getAbsolutePath();
					
					String argss = "--threads 4 --min-dispatch 4 --plot-slip-rates" +
//					String argss = "--threads 4 --min-dispatch 4 --plot-all --no-erf-plots" +
//					String argss = "--threads 4 --min-dispatch 4 --build-mean --plot-parent-mfds --plot-paleo-faults" +
							" --name-grep "+argVal+" "+compoundFile.getAbsolutePath()+" "+remoteJobDir.getAbsolutePath();
					
					List<String> script = mpjWrite.buildScript(MPJDistributedCompoundFSSPlots.class.getName(), argss);
					
					String scriptName = remoteJobDir.getName()+".pbs";
					
					batchWrite.writeScript(new File(writeDir, scriptName), script, jobMins, nodes, site.getPPN(null), null);
				}
			}
		}
	}

}
