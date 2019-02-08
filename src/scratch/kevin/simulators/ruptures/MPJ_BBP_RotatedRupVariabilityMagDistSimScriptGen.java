package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityMagDistPageGen.RuptureType;

class MPJ_BBP_RotatedRupVariabilityMagDistSimScriptGen {

	public static void main(String[] args) throws IOException {
		// REMOTE paths
		@SuppressWarnings("unused")
		File myHPCDir = new File("/auto/scec-02/kmilner/simulators/catalogs/");
		File stampedeCatalogDir = new File("/work/00950/kevinm/stampede2/simulators/catalogs");
		File jacquiCSDir = new File("/home/scec-00/gilchrij/RSQSim/CISM/cybershake/");
//		File catalogDir = new File(myHPCDir, "rundir2585_1myrs");
		File catalogDir = new File(myHPCDir, "rundir2740");
//		File catalogDir = new File(stampedeCatalogDir, "rundir2829");
		
		int skipYears = 5000;

//		RuptureType[] rupTypes = RuptureType.values();
		RuptureType[] rupTypes = { RuptureType.REVERSE };
		boolean writeIndividual = true;
		
		VelocityModel vm = RSQSimBBP_Config.VM;
		
		double minDist = 20d;
		int numDist = 8;
		double deltaDist = 20;
		System.out.println("Dist range: "+(float)minDist+" "+(float)(minDist + (numDist-1)*deltaDist));
		
		double minMag = 6.5d;
		int numMag = 11;
		double deltaMag = 0.1;
		System.out.println("Mag range: "+(float)minMag+" "+(float)(minMag + (numMag-1)*deltaMag));
		
		int numSourceAz = 36;
		int numSiteToSourceAz = 1;
		int minRuptures = 10;
		int maxRuptures = 100;

		System.out.println("Min num: "+(numSourceAz*numSiteToSourceAz*numDist*numMag*rupTypes.length*minRuptures));
		System.out.println("Max num: "+(numSourceAz*numSiteToSourceAz*numDist*numMag*rupTypes.length*maxRuptures));
		
		double timeScalar = 1d;
		boolean scaleVelocities = true;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 7*24*60;
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
		Preconditions.checkState(catalogDir.getAbsolutePath().contains("scec-"),
				"You forgot the catalog dir on HPC, dummy");
		
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
//		Preconditions.checkState(catalogDir.getAbsolutePath().contains("kevinm"),
//				"You forgot the catalog dir on Stampede, dummy");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-"+catalogDir.getName()+"-rotatedRupsMagDist";
		if (rupTypes.length == 1)
			jobName += "-"+rupTypes[0].getPrefix();
		else
			jobName += "-"+rupTypes.length+"types";
		jobName += "-"+numDist+"dists";
		jobName += "-"+numMag+"mags";
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
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		writeScript(catalogDir, skipYears, rupTypes, minDist, numDist, deltaDist, minMag, numMag, deltaMag, numSourceAz,
				numSiteToSourceAz, minRuptures, maxRuptures, timeScalar, scaleVelocities, threads, nodes, queue, mins,
				bbpDataDir, nodeScratchDir, bbpCopyParentDir, bbpEnvFile, sharedScratchDir, pbsWrite, mpjWrite,
				localJobDir, remoteJobDir, "cat_bbp_rotated.slurm", vm);
		if (writeIndividual && rupTypes.length > 1) {
			for (RuptureType rupType : rupTypes) {
				RuptureType[] myTypes = { rupType };
				File myRemoteDir = new File(remoteJobDir, rupType.getPrefix());
				new File(localJobDir, rupType.getPrefix()).mkdir();
				writeScript(catalogDir, skipYears, myTypes, minDist, numDist, deltaDist, minMag, numMag, deltaMag, numSourceAz,
						numSiteToSourceAz, minRuptures, maxRuptures, timeScalar, scaleVelocities, threads, nodes, queue, mins,
						bbpDataDir, nodeScratchDir, bbpCopyParentDir, bbpEnvFile, sharedScratchDir, pbsWrite, mpjWrite,
						localJobDir, myRemoteDir, rupType.getPrefix()+".slurm", vm);
			}
		}
	}

	private static void writeScript(File catalogDir, int skipYears, RuptureType[] rupTypes, double minDist, int numDist,
			double deltaDist, double minMag, int numMag, double deltaMag, int numSourceAz, int numSiteToSourceAz,
			int minRuptures, int maxRuptures, double timeScalar, boolean scaleVelocities, int threads, int nodes,
			String queue, int mins, String bbpDataDir, String nodeScratchDir, String bbpCopyParentDir, File bbpEnvFile,
			String sharedScratchDir, BatchScriptWriter pbsWrite, JavaShellScriptWriter mpjWrite, File localJobDir,
			File remoteJobDir, String scriptFileName, VelocityModel vm) throws IOException {
		String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).maxDispatch(500).threads(threads).endTimeSlurm().build();
		argz += " --vm "+vm+" --method "+RSQSimBBP_Config.METHOD.name();
		argz += " --catalog-dir "+catalogDir.getAbsolutePath();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --time-step "+(float)RSQSimBBP_Config.SRF_DT+" --srf-interp "+RSQSimBBP_Config.SRF_INTERP_MODE.name();
		argz += " --skip-years "+skipYears;
		argz += " --rupture-type "+rupTypes[0].name();
		for (int i=1; i<rupTypes.length; i++)
			argz += ","+rupTypes[i].name();
		argz += " --min-distance "+minDist+" --num-distance "+numDist+" --delta-distance "+deltaDist;
		argz += " --min-mag "+minMag+" --num-mag "+numMag+" --delta-mag "+deltaMag;
		argz += " --num-source-az "+numSourceAz;
		argz += " --num-site-source-az "+numSiteToSourceAz;
		if (minRuptures > 0 && minRuptures < Integer.MAX_VALUE)
			argz += " --min-ruptures "+minRuptures;
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
		
		List<String> script = mpjWrite.buildScript(MPJ_BBP_RotatedRupVariabilityMagDistSim.class.getName(), argz);
		
		if (!addLines.isEmpty())
			script.addAll(2, addLines);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, scriptFileName), script);
	}

}
