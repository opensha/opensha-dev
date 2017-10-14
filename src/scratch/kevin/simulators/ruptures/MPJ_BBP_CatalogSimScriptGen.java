package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;

public class MPJ_BBP_CatalogSimScriptGen {

	public static void main(String[] args) throws IOException {
		// REMOTE paths
		File myHPCDir = new File("/auto/scec-02/kmilner/simulators/catalogs/");
		File jacquiCSDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/cybershake/");
//		File catalogDir = new File(jacquiCSDir, "UCERF3_millionElement");
//		File catalogDir = new File(jacquiCSDir, "rundir2194_long");
//		File catalogDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/cybershake/rundir2194_long");
		File catalogDir = new File(myHPCDir, "rundir2273");
		
		double slipVel = 1d;
		double minMag = 6;
		int skipYears = 5000;
		
		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 24*60;
		int heapSizeMB = 60*1024;
		String bbpDataDir = "${TMPDIR}";
		
		File localDir = new File("/home/kevin/bbp/parallel");
		File remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-"+catalogDir.getName()+"-all-m"+(float)minMag+"-skipYears"+skipYears;
		if (!RSQSimBBP_Config.DO_HF)
			jobName += "-noHF";
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		// copy sites file
		List<BBP_Site> sites = RSQSimBBP_Config.getStandardSites();
		File sitesFile = new File(localJobDir, "sites.stl");
		BBP_Site.writeToFile(sitesFile, sites);
		File remoteSitesFile = new File(remoteJobDir, sitesFile.getName());
		
		String argz = "--max-dispatch 1000 --threads "+threads;
		argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --sites-file "+remoteSitesFile.getAbsolutePath();
		argz += " --catalog-dir "+catalogDir.getAbsolutePath();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
		argz += " --slip-velocity "+(float)slipVel;
		argz += " --min-mag "+(float)minMag+" --skip-years "+skipYears;
		if (RSQSimBBP_Config.DO_HF)
			argz += " --no-hf";
		if (bbpDataDir != null && !bbpDataDir.isEmpty())
			argz += " --bbp-data-dir "+bbpDataDir;
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		MPJExpressShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		List<String> script = mpjWrite.buildScript(MPJ_BBP_CatalogSim.class.getName(), argz);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "cat_bbp_parallel.pbs"), script);
	}

}
