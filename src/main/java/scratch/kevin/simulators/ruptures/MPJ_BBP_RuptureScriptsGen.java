package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Somerville_2006_MagAreaRel;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.SRF_PointData;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.MPJ_BBP_RupGenSim;
import scratch.kevin.bbp.MPJ_BBP_ShakeMapSim;
import scratch.kevin.bbp.MPJ_BBP_Utils;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

class MPJ_BBP_RuptureScriptsGen {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int eventID = 136704;
//		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
//		int eventID = 4099020;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2273.instance(baseDir);
//		int eventID = 412778;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2310.instance(baseDir);
//		int eventID = 3802809;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2320.instance(baseDir);
//		int eventID = 6195527;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2336.instance(baseDir);
//		int eventID = 131670;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2337.instance(baseDir);
//		int eventID = 203769;
//		RSQSimCatalog catalog = Catalogs.JG_2194_K2.instance(baseDir);
//		int eventID = 18840012;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int eventID = 526885;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
////		int eventID = 81854;
////		int eventID = 2637969;
//		int eventID = 1670183;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
//		int eventID = 385955;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
//		int eventID = 5304;
//		int eventID = 31324;
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		int eventID = 9955310;
//		RSQSimCatalog catalog = Catalogs.BRUCE_4322.instance(baseDir);
//		int eventID = 40636;
////		int eventID = 92236;
//		RSQSimCatalog catalog = Catalogs.BRUCE_4827.instance(baseDir);
//		int eventID = 195167;
////		int eventID = 128149;
//		RSQSimCatalog catalog = Catalogs.BRUCE_4841.instance(baseDir);
//		int eventID = 755070;
////		int eventID = 2441060;
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance(baseDir);
////		int eventID = 51863;
////		int eventID = 39055;
//		int eventID = 12581;
////		int eventID = 77272;
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance(baseDir);
//		int eventID = 1499589;
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(baseDir);
//		int eventID = 7119753;
//		int eventID = 8242900;
//		int eventID = 1499589;
//		int eventID = 7028377;
//		int eventID = 13383629;
//		int eventID = 6553169;
//		int eventID = 13272163;
//		int eventID = 1651575;
//		int eventID = 3015199;
//		int eventID = 1879413;
//		int eventID = 10096082;
//		int eventID = 7992279;
//		int eventID = 7748150;
//		int eventID = 3012841;
//		int eventID = 2614773;
		int eventID = 1809975;
		
		double timeScalar = 1d;
		boolean scaleVelocities = false;
		
		File catalogDir = catalog.getCatalogDir();
		
		boolean stampede = false;
		
		boolean doGP = true;
		boolean doShakeMap = false;
		boolean doGPShakeMaps = false;

		boolean csSites = false;
		boolean cs500Sites = true;
		
		int numGP = 400;
		int numShakeMapGP = 5;
		double mapSpacing = 0.05;
//		double mapSpacing = 0.02;
//		int maxNodes = 16;
		int maxNodes = 8;
//		int maxNodes = 10;
		boolean gpAdjustDDW = false;
		
		File srcFile = RSQSimBBP_Config.getEventSrcFile(catalog, eventID);
		File srfFile = RSQSimBBP_Config.getEventSRFFile(catalog, eventID, RSQSimBBP_Config.SRF_INTERP_MODE,
				RSQSimBBP_Config.SRF_DT, timeScalar, scaleVelocities);
		if (!srcFile.exists() || (doShakeMap && !srfFile.exists())) {
			System.out.println("need to load event "+eventID);
			RSQSimEvent event = catalog.loader().byID(eventID);
			System.out.println("done loading");
			
			RSQSimBBP_Config.generateBBP_Inputs(catalog, event, false, timeScalar, scaleVelocities);
			Preconditions.checkState(srcFile.exists(), "SRC file not generated?");
			Preconditions.checkState(srfFile.exists(), "SRF file not generated?");
			
			System.out.println("done writing event inputs");
		}
		Preconditions.checkState(srcFile.exists(), "Source file doesn't exist: %s", srcFile.getAbsolutePath());
		
		int gpMins = 24*60;
		int mapMins = 24*60;
		boolean gpSplitSites = true;
		int heapSizeMB = 40*1024;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads;
		String queue;
		File remoteDir;
		String bbpDataDir;
		String nodeScratchDir;
		String bbpCopyParentDir;
		String nodeGFDir;
		File bbpEnvFile;
		String sharedScratchDir;
		BatchScriptWriter pbsWrite;
		
		JavaShellScriptWriter mpjWrite;
		if (stampede) {
			threads = 96;
			queue = "skx-normal";
			remoteDir = new File("/work/00950/kevinm/stampede2/bbp/parallel");
			bbpDataDir = "/tmp";
			nodeScratchDir = null;
			nodeGFDir = "/tmp/gfs";
			bbpCopyParentDir = "/scratch/00950/kevinm/bbp";
			bbpEnvFile = new File("/work/00950/kevinm/stampede2/bbp/bbp_env.sh");
//			sharedScratchDir = "/scratch/00950/kevinm/";
			sharedScratchDir = null;
			pbsWrite = new StampedeScriptWriter(true);
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, heapSizeMB, null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		} else {
//			threads = 20;
//			queue = "scec";
//			remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
//			bbpDataDir = USC_HPCC_ScriptWriter.NODE_TEMP_DIR;
//			nodeScratchDir = null;
//			nodeGFDir = USC_HPCC_ScriptWriter.NODE_TEMP_DIR+"/gfs";
//			bbpCopyParentDir = USC_HPCC_ScriptWriter.SHARED_SCRATCH_DIR+"/kmilner";
//			bbpEnvFile = new File("/auto/scec-02/kmilner/bbp/bbp_env.sh");
////			sharedScratchDir = "${SCRATCHDIR}";
//			sharedScratchDir = null;
//			pbsWrite = new USC_HPCC_ScriptWriter();
//			mpjWrite = new MPJExpressShellScriptWriter(
//					USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, null, USC_HPCC_ScriptWriter.MPJ_HOME);
			threads = 20;
			queue = "scec";
			bbpDataDir = USC_CARC_ScriptWriter.NODE_TEMP_DIR;
			nodeScratchDir = null;
//			String bbpCopyParentDir = USC_CARC_ScriptWriter.SHARED_SCRATCH_DIR+"/kmilner";
			bbpCopyParentDir = null;
			nodeGFDir = USC_CARC_ScriptWriter.NODE_TEMP_DIR+"/gfs";
			bbpEnvFile = new File("/project/scec_608/kmilner/bbp/bbp_env.sh");
//			String sharedScratchDir = "${SCRATCHDIR}";
			sharedScratchDir = null;
			remoteDir = new File("/project/scec_608/kmilner/bbp/parallel");
			pbsWrite = new USC_CARC_ScriptWriter();
//			List<File> classpath = new ArrayList<>();
//			classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
//			JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//					USC_CARC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_CARC_ScriptWriter.MPJ_HOME);
//			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			mpjWrite = new FastMPJShellScriptWriter(
					USC_CARC_ScriptWriter.JAVA_BIN, heapSizeMB, null, USC_CARC_ScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		}
		
		String dateStr = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		mpjWrite.setClasspath(classpath);
		
		BBP_SourceFile src = BBP_SourceFile.readFile(srcFile);
		BBP_SourceFile gpSRC = null;
		if (gpAdjustDDW && (doGP || doGPShakeMaps)) {
			Somerville_2006_MagAreaRel ma = new Somerville_2006_MagAreaRel();
			double newArea = ma.getMedianArea(src.getMag());
			BBP_PlanarSurface origSurf = src.getSurface();
			double newWidth = newArea/origSurf.getLength();
			BBP_PlanarSurface newSurf = new BBP_PlanarSurface(origSurf.getTopCenter(), origSurf.getLength(),
					newWidth, origSurf.getFocalMechanism());
			gpSRC = new BBP_SourceFile(newSurf, src.getMag(), src.getHypoAlongStrike(), src.getHypoDownDip(),
					src.getdWid(), src.getdLen(), src.getCornerFreq(), src.getSeed());
		} else {
			gpSRC = src;
		}
		
		if (doGP) {
			String jobName = dateStr;
			jobName += "-"+catalogDir.getName()+"-event"+eventID+"-gp";
			
			String len = gpSRC.getdLen()+"";
			if (len != null)
				jobName += "-dx"+len;
			if (!RSQSimBBP_Config.DO_HF)
				jobName += "-noHF";
			if (stampede)
				jobName += "-stampede";
			if (csSites)
				jobName += "-csLASites";
			if (cs500Sites)
				jobName += "-cs500Sites";
			if (gpAdjustDDW)
				jobName += "-adjustDDW";
			
			System.out.println("Writing GP sim to "+jobName);
			
			List<BBP_Site> sites;
			if (csSites)
				sites = RSQSimBBP_Config.getCyberShakeInitialLASites();
			else if (cs500Sites)
				sites = RSQSimBBP_Config.getCyberShakeVs500LASites();
			else
				sites = RSQSimBBP_Config.getStandardSites(gpSRC);
			
			Preconditions.checkState(!sites.isEmpty(), "No sites!");
			int totalSims;
			if (gpSplitSites)
				totalSims = numGP*sites.size();
			else
				totalSims = numGP;
			int nodes = totalSims/threads;
			while (nodes > maxNodes)
				nodes /= 2;
			
			System.out.println("Doing "+totalSims+" total sims on "+nodes+" nodes with "+threads+" threads");
			
			File localJobDir = new File(localDir, jobName);
			System.out.println(localJobDir.getAbsolutePath());
			Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
			File remoteJobDir = new File(remoteDir, jobName);
			
			// write sites file
			File sitesFile = new File(localJobDir, "sites.stl");
			BBP_Site.writeToFile(sitesFile, sites);
			File remoteSitesFile = new File(remoteJobDir, sitesFile.getName());
			
			// write src file
			String srcName = srcFile.getName();
			if (gpAdjustDDW)
				srcName = srcName.substring(0, srcName.indexOf(".src"))+"_adjust_ddw.src";
			gpSRC.writeToFile(new File(localJobDir, srcName));
			File remoteSrcFile = new File(remoteJobDir, srcName);
			
			String argz = MPJTaskCalculator.argumentBuilder().exactDispatch(threads).threads(threads).endTimeSlurm().build();
			argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
			argz += " --sites-file "+remoteSitesFile.getAbsolutePath();
			argz += " --src-file "+remoteSrcFile.getAbsolutePath();
			argz += " --output-dir "+remoteJobDir.getAbsolutePath();
			argz += " --num-sims "+numGP;
			if (gpSplitSites)
				argz += " --split-sites";
			if (!RSQSimBBP_Config.DO_HF)
				argz += " --no-hf";
			if (bbpDataDir != null && !bbpDataDir.isEmpty())
				argz += " --bbp-data-dir "+bbpDataDir;
			if (nodeGFDir != null && !nodeGFDir.isEmpty())
				argz += " --node-gf-dir "+nodeGFDir;

			List<String> addLines = new ArrayList<>();
			argz = MPJ_BBP_CatalogSimScriptGen.addBBP_EnvArgs(argz, addLines, remoteJobDir, nodeScratchDir,
					sharedScratchDir, bbpCopyParentDir, bbpEnvFile);
			
			List<String> script = mpjWrite.buildScript(MPJ_BBP_RupGenSim.class.getName(), argz);
			
			if (!addLines.isEmpty())
				script.addAll(2, addLines);
			
			script = pbsWrite.buildScript(script, gpMins, nodes, threads, queue);
			pbsWrite.writeScript(new File(localJobDir, "gp_bbp_parallel.slurm"), script);
		}
		if (doShakeMap) {
			String jobName = dateStr;
			jobName += "-"+catalogDir.getName()+"-event"+eventID+"-shakemap";
			if (timeScalar != 1d) {
				jobName += "-timeScale"+(float)timeScalar;
				if (scaleVelocities)
					jobName += "-velScale";
			}
			if (!RSQSimBBP_Config.DO_HF)
				jobName += "-noHF";
			if (stampede)
				jobName += "-stampede";
			
			System.out.println("Writing ShakeMap sim to "+jobName);
			
			int nodes = maxNodes;
			
			File localJobDir = new File(localDir, jobName);
			System.out.println(localJobDir.getAbsolutePath());
			Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
			File remoteJobDir = new File(remoteDir, jobName);
			
			// copy src file
			Files.copy(srcFile, new File(localJobDir, srcFile.getName()));
			File remoteSrcFile = new File(remoteJobDir, srcFile.getName());
			
			// copy srf file
			Files.copy(srfFile, new File(localJobDir, srfFile.getName()));
			File remoteSrfFile = new File(remoteJobDir, srfFile.getName());
			String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).threads(threads).endTimeSlurm().build();
			argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
			argz += " --src-file "+remoteSrcFile.getAbsolutePath();
			argz += " --srf-file "+remoteSrfFile.getAbsolutePath();
			argz += " --output-dir "+remoteJobDir.getAbsolutePath();
			argz += " --buffer-dist "+RSQSimBBP_Config.MAX_DIST+" --spacing "+mapSpacing;
			if (!RSQSimBBP_Config.DO_HF)
				argz += " --no-hf";
			if (bbpDataDir != null && !bbpDataDir.isEmpty())
				argz += " --bbp-data-dir "+bbpDataDir;
			
			List<String> addLines = new ArrayList<>();
			argz = MPJ_BBP_CatalogSimScriptGen.addBBP_EnvArgs(argz, addLines, remoteJobDir, nodeScratchDir,
					sharedScratchDir, bbpCopyParentDir, bbpEnvFile);
			
			List<String> script = mpjWrite.buildScript(MPJ_BBP_ShakeMapSim.class.getName(), argz);
			
			if (!addLines.isEmpty())
				script.addAll(2, addLines);
			
			script = pbsWrite.buildScript(script, mapMins, nodes, threads, queue);
			pbsWrite.writeScript(new File(localJobDir, "map_bbp_parallel.slurm"), script);
		}
		if (doGPShakeMaps) {
			String jobPrefix = dateStr;
			jobPrefix += "-"+catalogDir.getName()+"-event"+eventID+"-shakemap-gp";
			if (!RSQSimBBP_Config.DO_HF)
				jobPrefix += "-noHF";
			if (stampede)
				jobPrefix += "-stampede";
			if (gpAdjustDDW)
				jobPrefix += "-adjustDDW";
			
			for (int i=0; i<numShakeMapGP; i++) {
				int seed = 1000 + i;
				String jobName = jobPrefix+"-seed"+seed;
				System.out.println("Writing ShakeMap sim to "+jobName);
				
				int nodes = maxNodes;
				
				File localJobDir = new File(localDir, jobName);
				System.out.println(localJobDir.getAbsolutePath());
				Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
				File remoteJobDir = new File(remoteDir, jobName);
				
				// write src file with unique seed
				gpSRC.setSeed(seed);
				File localSrcFile = new File(localJobDir, srcFile.getName());
				gpSRC.writeToFile(localSrcFile);
				File remoteSrcFile = new File(remoteJobDir, srcFile.getName());
				
				File runDir = new File(localJobDir, "srf_build");
				MPJ_BBP_Utils.waitOnDir(runDir, 5, 1000);
				BBP_Wrapper wrapper = new BBP_Wrapper(RSQSimBBP_Config.VM, RSQSimBBP_Config.METHOD, localSrcFile,
						null, null, null, runDir);
				wrapper.setSRFGenOnly(true);
				wrapper.run();
				
				File newSRF = null;
				for (File file : runDir.listFiles())
					if (file.getName().endsWith(".srf"))
						newSRF = file;
				Preconditions.checkNotNull(newSRF, "Couldn't file new SRF file in "+runDir.getAbsolutePath());
				Files.copy(newSRF, srfFile);
				
				Files.copy(newSRF, new File(localJobDir, newSRF.getName()));
				File remoteSrfFile = new File(remoteJobDir, newSRF.getName());
				String argz = MPJTaskCalculator.argumentBuilder().minDispatch(threads).threads(threads).endTimeSlurm().build();
				argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
				argz += " --src-file "+remoteSrcFile.getAbsolutePath();
				argz += " --srf-file "+remoteSrfFile.getAbsolutePath();
				argz += " --output-dir "+remoteJobDir.getAbsolutePath();
				argz += " --buffer-dist "+RSQSimBBP_Config.MAX_DIST+" --spacing "+mapSpacing;
				if (!RSQSimBBP_Config.DO_HF)
					argz += " --no-hf";
				if (bbpDataDir != null && !bbpDataDir.isEmpty())
					argz += " --bbp-data-dir "+bbpDataDir;
				
				List<String> addLines = new ArrayList<>();
				argz = MPJ_BBP_CatalogSimScriptGen.addBBP_EnvArgs(argz, addLines, remoteJobDir, nodeScratchDir,
						sharedScratchDir, bbpCopyParentDir, bbpEnvFile);
				
				List<String> script = mpjWrite.buildScript(MPJ_BBP_ShakeMapSim.class.getName(), argz);
				
				if (!addLines.isEmpty())
					script.addAll(2, addLines);
				
				script = pbsWrite.buildScript(script, mapMins, nodes, threads, queue);
				pbsWrite.writeScript(new File(localJobDir, "map_bbp_parallel.slurm"), script);
			}
		}
	}

}
