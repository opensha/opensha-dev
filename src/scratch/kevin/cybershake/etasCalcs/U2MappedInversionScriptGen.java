package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter.RunSites;

public class U2MappedInversionScriptGen {

	public static void main(String[] args) throws IOException {
		int numRuns = 10;
		RunSites site = RunSites.HPCC;
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		int annealHours = 5;
		String jobName = "2015_02_18-ucerf2-rups-inversion";
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/u2_inverted/"+jobName);
		if (!localDir.exists())
			localDir.mkdir();
		File remoteDir = new File("/home/scec-02/kmilner/ucerf3/inversions/"+jobName);
		
		int ppn = 8;
		int mins = (annealHours+1)*60;
		int memMB = 8*1024;
		
		JavaShellScriptWriter javaWriter = new JavaShellScriptWriter(site.getJAVA_BIN(), memMB,
				LogicTreePBSWriter.getClasspath(site, remoteDir));
		javaWriter.setHeadless(true);
		if (site.getFM_STORE() != null) {
			javaWriter.setProperty(FaultModels.FAULT_MODEL_STORE_PROPERTY_NAME, site.getFM_STORE());
		}
		
		for (int r=0; r<numRuns; r++) {
			File pbsFile = new File(localDir, "inversion_"+r+".pbs");
			String branchStr = "FM2_1_UC2ALL_AveU2_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU2_run"+r;
			String classArgs = "--completion-time "+annealHours+"h --sub-completion 1s --cool FAST_SA --nonneg LIMIT_ZERO_RATES "
					+ "--num-threads 100% --branch-prefix "+branchStr+ " --u2-rups-only --directory "+remoteDir.getAbsolutePath();
			List<String> script = javaWriter.buildScript(CommandLineInversionRunner.class.getName(), classArgs);
			pbsWrite.writeScript(pbsFile, script, mins, 1, ppn, null);
		}
	}

}
