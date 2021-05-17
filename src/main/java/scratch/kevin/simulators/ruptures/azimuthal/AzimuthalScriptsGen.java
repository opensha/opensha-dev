package scratch.kevin.simulators.ruptures.azimuthal;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter.Device;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.MPJ_BBP_CatalogSimScriptGen;

public class AzimuthalScriptsGen {

	public static void main(String[] args) throws IOException {
//		String catalogDirName = "rundir2585_1myrs";
		String catalogDirName = "rundir5046";
//		String catalogDirName = "rundir4860_multi_combine";
//		String catalogDirName = "rundir4983_stitched";
		boolean gp = false;

//		String catalogDirName = null;
//		boolean gp = true;
		
		int skipYears = 20000;
		double gpPatchArea = 1d;
		
		int numRuptures = 100;
		double spacing = 5d;
//		double spacing = 10d;
		double buffer = 100d;
		boolean positiveX = false;

//		Scenario[] scenarios = Scenario.values();
//		Scenario[] scenarios = {Scenario.M6p6_VERT_SS_SURFACE_RELAXED, Scenario.M7p2_VERT_SS_SURFACE_RELAXED};
//		Scenario[] scenarios = {Scenario.M7p2_VERT_SS_SURFACE, Scenario.M7p2_VERT_SS_SURFACE_RANDMAG};
//		Scenario[] scenarios = {Scenario.M7p2_VERT_SS_SURFACE_RANDMAG_0p15, Scenario.M7p2_VERT_SS_SURFACE_RANDMAG_0p2,
//				Scenario.M7p2_VERT_SS_SURFACE_RANDMAG_0p25, Scenario.M7p2_VERT_SS_SURFACE_RANDMAG_0p3};
//		Scenario[] scenarios = {Scenario.M7p2_VERT_SS_SURFACE};
//		Scenario[] scenarios = {Scenario.M6p6_REVERSE};
		Scenario[] scenarios = {Scenario.M6p6_VERT_SS_SURFACE, Scenario.M6p6_REVERSE,
				Scenario.M7p2_VERT_SS_SURFACE};
//		Scenario[] scenarios = {Scenario.M6p6_VERT_SS_SURFACE,
//				Scenario.M7p2_VERT_SS_SURFACE};
////		double[] distances = BBP_PartBValidationConfig.OFFICIAL_DISTANCES;
		FilterMethod filter = FilterMethod.MEDIAN_LENGTH;
//		FilterMethod filter = BBP_PartBValidationConfig.FILTER_METHOD_DEFAULT;
		
//		RSQSimBBP_Config.VM = VelocityModel.LA_BASIN_863;
		VelocityModel vm = RSQSimBBP_Config.VM;
		
		double timeScalar = 1d;
		boolean scaleVelocities = false;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads = 20;
		int nodes = 30;
		String queue = "scec";
		int mins = 48*60;
		int heapSizeMB = 45*1024;
		String bbpDataDir = USC_HPCC_ScriptWriter.NODE_TEMP_DIR;
		String nodeScratchDir = null;
		String bbpCopyParentDir = USC_HPCC_ScriptWriter.SHARED_SCRATCH_DIR+"/kmilner";
		String nodeGFDir = USC_HPCC_ScriptWriter.NODE_TEMP_DIR+"/gfs";
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
		
//		int threads = 48;
//		int nodes = 10;
//		String queue = "skx-normal";
//		int mins = 24*60;
//		int heapSizeMB = 100*1024;
//		String bbpDataDir = "/tmp";
//		String nodeScratchDir = null;
//		String nodeGFDir = "/tmp/gfs";
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
		if (gp) {
			Preconditions.checkState(catalogDirName == null, "can't supply catalog with GP");
			jobName += "-gp";
		} else {
			Preconditions.checkNotNull(catalogDirName, "GP is false but no catalog specified");
			jobName += "-"+catalogDirName;
		}
		jobName += "-azimuthal";
		if (scenarios.length == 1)
			jobName += "-"+scenarios[0].getPrefix();
		else
			jobName += "-"+scenarios.length+"scenarios";
		if (!gp)
			jobName += "-filter_"+filter.getPrefix();
		jobName += "-"+numRuptures+"rups";
		jobName += "-buffer"+(int)buffer;
		jobName += "-spacing"+(int)spacing;
		if (positiveX)
			jobName += "-positiveX";
		if (gp)
			jobName += "-patchArea"+(float)gpPatchArea;
		else
			jobName += "-skipYears"+skipYears;
		jobName += "-vm"+vm.name();
		if (!RSQSimBBP_Config.DO_HF)
			jobName += "-noHF";
		if (timeScalar != 1d) {
			jobName += "-timeScale"+(float)timeScalar;
			if (scaleVelocities)
				jobName += "-velScale";
		}
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).maxDispatch(500).threads(threads).endTimeSlurm().build();
		argz += " --vm "+vm.name()+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		if (!gp) {
			argz += " --catalog-dir "+catalogDirName;
			argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
			argz += " --skip-years "+skipYears;
			if (timeScalar != 1d) {
				argz += " --time-scalar "+(float)timeScalar;
				if (scaleVelocities)
					argz += " --velocity-scale";
			}
		} else {
			argz += " --patch-area "+(float)gpPatchArea;
		}
		argz += " --num-ruptures "+numRuptures;
		argz += " --scenarios "+scenarios[0].name();
		for (int i=1; i<scenarios.length; i++)
			argz += ","+scenarios[i].name();
		if (!gp)
			argz += " --filter "+filter.name();
		argz += " --buffer "+(float)buffer;
		argz += " --spacing "+(float)spacing;
		if (positiveX)
			argz += " --positive-x";
		if (!RSQSimBBP_Config.DO_HF)
			argz += " --no-hf";
		if (bbpDataDir != null && !bbpDataDir.isEmpty())
			argz += " --bbp-data-dir "+bbpDataDir;
		if (nodeGFDir != null && !nodeGFDir.isEmpty())
			argz += " --node-gf-dir "+nodeGFDir;
		List<String> addLines = new ArrayList<>();
		argz = MPJ_BBP_CatalogSimScriptGen.addBBP_EnvArgs(argz, addLines, remoteJobDir, nodeScratchDir,
				sharedScratchDir, bbpCopyParentDir, bbpEnvFile);
		
		String clazz = gp ? MPJ_GP_AzimuthalSim.class.getName()
				: MPJ_RSQSimAzimuthalSim.class.getName();
		
		List<String> script = mpjWrite.buildScript(clazz, argz);
		
		if (!addLines.isEmpty())
			script.addAll(2, addLines);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		String scriptName = gp ? "gp_bbp_azimuthal.slurm" : "cat_bbp_azimuthal.slurm";
		pbsWrite.writeScript(new File(localJobDir, scriptName), script);
	}

}
