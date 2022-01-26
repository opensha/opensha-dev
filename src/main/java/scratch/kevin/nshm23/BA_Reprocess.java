package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class BA_Reprocess {

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteToalThreads = 20;
		int remoteTotalMemGB = 55;
		String queue = "scec";
		int nodes = 2;
		int mins = 300;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteMainDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, classpath, USC_CARC_ScriptWriter.FMPJ_HOME);
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();

//		String dirName = "2021_11_23-u3_branches-FM3_1-5h";
//		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.class;
//		String dirName = "2021_11_24-nshm23_draft_branches-FM3_1";
//		Class<? extends InversionConfigurationFactory> factoryClass = DraftNSHM23InvConfigFactory.class;
//		String dirName = "2021_11_30-nshm23_draft_branches-FM3_1-FaultSpec";
//		Class<? extends InversionConfigurationFactory> factoryClass = DraftNSHM23InvConfigFactory.class;
//		String dirName = "2021_12_01-nshm23_draft_branches-no_paleo-no_parkfield-FM3_1-SysAvg";
//		Class<? extends InversionConfigurationFactory> factoryClass = DraftNSHM23InvConfigFactory.NoPaleoParkfield.class;
//		String dirName = "2021_12_01-u3_branches-no_paleo-no_parkfield-single_mfd_reg-FM3_1-5h";
//		Class<? extends InversionConfigurationFactory> factoryClass = U3InversionConfigFactory.NoPaleoParkfieldSingleReg.class;
//		String dirName = "2021_12_03-nshm23_draft_branches-no_paleo-no_parkfield-FM3_1-FaultSpec";
//		Class<? extends InversionConfigurationFactory> factoryClass = DraftNSHM23InvConfigFactory.NoPaleoParkfield.class;
		String dirName = "2021_12_08-nshm23_draft_branches-FM3_1-TotNuclRate-SubB1-2h";
		Class<? extends InversionConfigurationFactory> factoryClass = NSHM23_InvConfigFactory.class;
		
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists());
		File remoteDir = new File(remoteMainDir, dirName);
		
		File remoteLogicTree = new File(remoteDir, "logic_tree.json");
		
		File resultsDir = new File(remoteDir, "results");
		String argz = "--logic-tree "+remoteLogicTree.getAbsolutePath();
		argz += " --output-dir "+resultsDir.getAbsolutePath();
		argz += " --inversion-factory '"+factoryClass.getName()+"'"; // surround in single quotes to escape $'s
		argz += " --annealing-threads 1";
		List<String> script = mpjWrite.buildScript(MPJ_LogicTreeInversionRunner.class.getName(), argz);
		
		pbsWrite.writeScript(new File(localDir, "ba_regen.slurm"), script, mins, nodes, remoteToalThreads, queue);
	}

}
