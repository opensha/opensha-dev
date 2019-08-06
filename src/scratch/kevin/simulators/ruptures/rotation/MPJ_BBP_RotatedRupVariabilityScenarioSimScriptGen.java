package scratch.kevin.simulators.ruptures.rotation;

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
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.MPJ_BBP_CatalogSimScriptGen;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

class MPJ_BBP_RotatedRupVariabilityScenarioSimScriptGen {

	public static void main(String[] args) throws IOException {
		// REMOTE paths
//		String catalogDirName = "rundir2585_1myrs";
		String catalogDirName = "rundir4317";
		
		int skipYears = 2000;

		Scenario[] scenarios = Scenario.values();
//		Scenario[] scenarios = {Scenario.M6p6_VERT_SS_SURFACE};
//		double[] distances = BBP_PartBValidationConfig.OFFICIAL_DISTANCES;
		double[] distances = { 20d, 50d, 100d };
		int numSourceAz = 18;
//		int numSiteToSourceAz = 36;
		int numSiteToSourceAz = 1;
		int maxRuptures = 400;
		
//		RSQSimBBP_Config.VM = VelocityModel.LA_BASIN_863;
		VelocityModel vm = RSQSimBBP_Config.VM;
		
//		List<BBP_Site> sites = RSQSimBBP_Config.getCyberShakeInitialLASites();
//		String stitesStr = "csLASites";
		
		List<BBP_Site> sites = RSQSimBBP_Config.getCyberShakeInitialLASites().subList(0, 1);
		String stitesStr = "1site";
		
		System.out.println("Expected num: "+(numSourceAz*numSiteToSourceAz*distances.length*scenarios.length*maxRuptures*sites.size()));
		
		double timeScalar = 1d;
		boolean scaleVelocities = true;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 48*60;
		int heapSizeMB = 45*1024;
		String bbpDataDir = "${TMPDIR}";
		String nodeScratchDir = null;
		String bbpCopyParentDir = "/staging/pjm/kmilner";
		File bbpEnvFile = new File("/auto/scec-02/kmilner/bbp/bbp_env.sh");
		String sharedScratchDir = "${SCRATCHDIR}";
		File remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
//		int threads = 96;
//		int nodes = 10;
//		String queue = "skx-normal";
//		int mins = 24*60;
//		int heapSizeMB = 100*1024;
//		String bbpDataDir = "/tmp";
//		String nodeScratchDir = null;
//		String bbpCopyParentDir = "/scratch/00950/kevinm/bbp";
//		File bbpEnvFile = new File("/work/00950/kevinm/stampede2/bbp/bbp_env.sh");
//		String sharedScratchDir = "/scratch/00950/kevinm/";
//		File remoteDir = new File("/work/00950/kevinm/stampede2/bbp/parallel");
//		BatchScriptWriter pbsWrite = new StampedeScriptWriter();
//		List<File> classpath = new ArrayList<>();
//		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
//		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, heapSizeMB, classpath,
//				StampedeScriptWriter.FMPJ_HOME, Device.NIODEV);
//		((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-"+catalogDirName+"-rotatedRups";
		if (scenarios.length == 1)
			jobName += "-"+scenarios[0].getPrefix();
		else
			jobName += "-"+scenarios.length+"scenarios";
		if (distances.length == 1)
			jobName += "-"+(float)distances[0]+"km";
		else
			jobName += "-"+distances.length+"dists";
		jobName += "-"+numSourceAz+"srcAz-"+numSiteToSourceAz+"siteSrcAz";
		if (maxRuptures > 0 && maxRuptures < Integer.MAX_VALUE)
			jobName += "-"+maxRuptures+"rups";
		jobName += "-skipYears"+skipYears;
		jobName += "-vm"+vm.name();
		if (!RSQSimBBP_Config.DO_HF)
			jobName += "-noHF";
		if (timeScalar != 1d) {
			jobName += "-timeScale"+(float)timeScalar;
			if (scaleVelocities)
				jobName += "-velScale";
		}
		jobName += "-"+stitesStr;
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		File sitesFile = new File(localJobDir, "sites.stl");
		System.out.println("Writing "+sites.size()+" sites to "+sitesFile.getAbsolutePath());
		BBP_Site.writeToFile(sitesFile, sites);
		File remoteSitesFile = new File(remoteJobDir, sitesFile.getName());
		
		String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).maxDispatch(500).threads(threads).endTimeSlurm().build();
		argz += " --vm "+vm.name()+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --catalog-dir "+catalogDirName;
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
		argz += " --skip-years "+skipYears;
		argz += " --scenarios "+scenarios[0].name();
		for (int i=1; i<scenarios.length; i++)
			argz += ","+scenarios[i].name();
		argz += " --distances "+commaSeparate(Doubles.asList(distances).toArray());
		argz += " --num-source-az "+numSourceAz;
		argz += " --num-site-source-az "+numSiteToSourceAz;
		argz += " --sites "+remoteSitesFile.getAbsolutePath();
		if (maxRuptures > 0 && maxRuptures < Integer.MAX_VALUE)
			argz += " --max-ruptures "+maxRuptures;
		if (!RSQSimBBP_Config.DO_HF)
			argz += " --no-hf";
		if (bbpDataDir != null && !bbpDataDir.isEmpty())
			argz += " --bbp-data-dir "+bbpDataDir;
		if (timeScalar != 1d) {
			argz += " --time-scalar "+(float)timeScalar;
			if (scaleVelocities)
				argz += " --velocity-scale";
		}
		List<String> addLines = new ArrayList<>();
		argz = MPJ_BBP_CatalogSimScriptGen.addBBP_EnvArgs(argz, addLines, remoteJobDir, nodeScratchDir,
				sharedScratchDir, bbpCopyParentDir, bbpEnvFile);
		
		List<String> script = mpjWrite.buildScript(MPJ_BBP_RotatedRupVariabilityScenarioSim.class.getName(), argz);
		
		if (!addLines.isEmpty())
			script.addAll(2, addLines);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "cat_bbp_rotated.slurm"), script);
	}
	
	private static String commaSeparate(Object... vals) {
		String[] strs = new String[vals.length];
		for (int i=0; i<vals.length; i++)
			strs[i] = vals[i].toString();
		return commaSeparate(strs);
	}
	
	private static String commaSeparate(String... vals) {
		Preconditions.checkState(vals.length > 0);
		String str = vals[0].toString();
		for (int i=1; i<vals.length; i++)
			str += ","+vals[i];
		return str;
	}

}
