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

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.bbp.BBP_Site;

class MPJ_BBP_PartBSimScriptGen {

	public static void main(String[] args) throws IOException {
		// REMOTE paths
		@SuppressWarnings("unused")
		String catalogDirName = "rundir4314";
		
		int numSites = 100;
		boolean randomAz = false;
		int skipYears = 2000;
		int maxRuptures = 500;
		
		double timeScalar = 1d;
		boolean scaleVelocities = true;
		
		VelocityModel vm = RSQSimBBP_Config.VM;
		float[] distances = { 20, 50, 100 };
		FilterMethod filter = BBP_PartBValidationConfig.FILTER_METHOD_DEFAULT;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 24*60;
		int heapSizeMB = 45*1024;
		String bbpDataDir = USC_HPCC_ScriptWriter.NODE_TEMP_DIR;
		String nodeScratchDir = null;
		String bbpCopyParentDir = USC_HPCC_ScriptWriter.SHARED_SCRATCH_DIR+"/kmilner";
		File bbpEnvFile = new File("/auto/scec-02/kmilner/bbp/bbp_env.sh");
//		String sharedScratchDir = "${SCRATCHDIR}";
		String sharedScratchDir = null;
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
		jobName += "-"+catalogDirName+"-partB-skipYears"+skipYears+"-"+numSites+"sites-vm"+vm.name()+"-"+distances.length+"dists";
		if (randomAz)
			jobName += "-randomAz";
		if (!RSQSimBBP_Config.DO_HF)
			jobName += "-noHF";
		if (timeScalar != 1d) {
			jobName += "-timeScale"+(float)timeScalar;
			if (scaleVelocities)
				jobName += "-velScale";
		}
		jobName += "-filter_"+filter.getPrefix();
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).maxDispatch(500).threads(threads).endTimeSlurm().build();
		argz += " --vm "+vm.name()+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --catalog-dir "+catalogDirName;
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
		argz += " --skip-years "+skipYears;
		argz += " --num-sites "+numSites;
		argz += " --distances "+Joiner.on(",").join(Floats.asList(distances));
		if (maxRuptures > 0)
			argz += " --max-ruptures "+maxRuptures;
		argz += " --filter "+filter.name();
		if (randomAz)
			argz += " --random-azimuth";
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
		
		List<String> script = mpjWrite.buildScript(MPJ_BBP_PartBSim.class.getName(), argz);
		
		if (!addLines.isEmpty())
			script.addAll(2, addLines);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "cat_bbp_partb.slurm"), script);
	}

}
