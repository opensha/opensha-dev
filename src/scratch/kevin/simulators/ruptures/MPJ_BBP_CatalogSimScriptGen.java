package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter.Device;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;

class MPJ_BBP_CatalogSimScriptGen {

	public static void main(String[] args) throws IOException {
		// REMOTE paths
		@SuppressWarnings("unused")
		File myHPCDir = new File("/auto/scec-02/kmilner/simulators/catalogs/");
		File stampedeCatalogDir = new File("/work/00950/kevinm/stampede2/simulators/catalogs");
		File jacquiCSDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/cybershake/");
//		File catalogDir = new File(jacquiCSDir, "UCERF3_millionElement");
//		File catalogDir = new File(jacquiCSDir, "rundir2194_long");
//		File catalogDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/cybershake/rundir2194_long");
//		File catalogDir = new File(myHPCDir, "rundir2342");
//		File catalogDir = new File(jacquiCSDir, "rundir2194_K2");
//		File catalogDir = new File(jacquiCSDir, "modLoad_testB");
//		File catalogDir = new File(jacquiCSDir, "tunedBase1m_ddotEQmod");
		File catalogDir = new File(stampedeCatalogDir, "rundir2616");
		
		boolean standardSites = true;
		boolean griddedSites = false;
		double griddedSpacing = 1d;
		
		double minMag = 6;
//		double minMag = 7;
		int numRG = 0;
//		double minMag = 7;
//		int numRG = 20;
		
		int skipYears = 5000;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
//		int threads = 20;
//		int nodes = 25;
//		String queue = "scec";
//		int mins = 24*60;
//		int heapSizeMB = 45*1024;
//		String bbpDataDir = "${TMPDIR}";
//		File remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
//		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
//		List<File> classpath = new ArrayList<>();
//		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
//		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
//		Preconditions.checkState(catalogDir.getAbsolutePath().contains("scec-"),
//				"You forgot the catalog dir on HPC, dummy");
		
		int threads = 96;
		int nodes = 10;
		String queue = "skx-normal";
		int mins = 24*60;
		int heapSizeMB = 100*1024;
		String bbpDataDir = "/tmp";
		File remoteDir = new File("/work/00950/kevinm/stampede2/bbp/parallel");
		BatchScriptWriter pbsWrite = new StampedeScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, heapSizeMB, classpath,
				StampedeScriptWriter.FMPJ_HOME, Device.NIODEV);
		((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		Preconditions.checkState(catalogDir.getAbsolutePath().contains("kevinm"),
				"You forgot the catalog dir on Stampede, dummy");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-"+catalogDir.getName()+"-all-m"+(float)minMag+"-skipYears"+skipYears;
		if (!RSQSimBBP_Config.DO_HF)
			jobName += "-noHF";
		if (numRG > 0)
			jobName += "-rg"+numRG;
		if (standardSites)
			jobName += "-standardSites";
		if (griddedSites)
			jobName += "-griddedSites";
		Preconditions.checkState(standardSites || griddedSites);
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		// copy sites file
		List<BBP_Site> sites = new ArrayList<>();
		if (standardSites)
			sites.addAll(RSQSimBBP_Config.getStandardSites());
		if (griddedSites)
			sites.addAll(RSQSimBBP_Config.getCAGriddedSites(griddedSpacing));
		File sitesFile = new File(localJobDir, "sites.stl");
		System.out.println("Writing "+sites.size()+" sites to "+sitesFile.getAbsolutePath());
		BBP_Site.writeToFile(sitesFile, sites);
		File remoteSitesFile = new File(remoteJobDir, sitesFile.getName());
		
		String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).maxDispatch(500).threads(threads).endTimeSlurm().build();
		argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --sites-file "+remoteSitesFile.getAbsolutePath();
		argz += " --catalog-dir "+catalogDir.getAbsolutePath();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
		argz += " --min-mag "+(float)minMag+" --skip-years "+skipYears;
		if (!RSQSimBBP_Config.DO_HF)
			argz += " --no-hf";
		if (bbpDataDir != null && !bbpDataDir.isEmpty())
			argz += " --bbp-data-dir "+bbpDataDir;
		if (numRG > 0)
			argz += " --rup-gen-sims "+numRG;
		
		List<String> script = mpjWrite.buildScript(MPJ_BBP_CatalogSim.class.getName(), argz);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "cat_bbp_parallel.pbs"), script);
	}

}
