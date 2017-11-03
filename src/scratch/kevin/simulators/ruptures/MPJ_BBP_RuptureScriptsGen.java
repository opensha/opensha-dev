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
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.SRF_PointData;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.MPJ_BBP_RupGenSim;
import scratch.kevin.bbp.MPJ_BBP_ShakeMapSim;
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
		RSQSimCatalog catalog = Catalogs.JG_2194_K2.instance(baseDir);
		int eventID = 18840012;
		
		File catalogDir = catalog.getCatalogDir();
		
		boolean stampede = false;
		
		boolean doGP = true;
		boolean doShakeMap = true;
		
		int numGP = 200;
		double mapSpacing = 0.05;
		int maxNodes = 36;
		
		File srcFile = RSQSimBBP_Config.getEventSrcFile(catalog, eventID);
		File srfFile = RSQSimBBP_Config.getEventSRFFile(catalog, eventID, RSQSimBBP_Config.SRF_INTERP_MODE, RSQSimBBP_Config.SRF_DT);
		BBP_SourceFile src = null;
		if (!srcFile.exists() || (doShakeMap && !srfFile.exists())) {
			System.out.println("need to load event "+eventID);
			RSQSimEvent event = catalog.loader().byID(eventID);
			System.out.println("done loading");
			
			RSQSimBBP_Config.generateBBP_Inputs(catalog, event, false);
			Preconditions.checkState(srcFile.exists(), "SRC file not generated?");
			Preconditions.checkState(srfFile.exists(), "SRF file not generated?");
			
			System.out.println("done writing event inputs");
		}
		Preconditions.checkState(srcFile.exists(), "Source file doesn't exist: %s", srcFile.getAbsolutePath());
		
		int gpMins = 10*60;
		int mapMins = 24*60;
		boolean gpSplitSites = true;
		int heapSizeMB = 40*1024;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		
		int threads;
		String queue;
		File remoteDir;
		String bbpDataDir;
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		JavaShellScriptWriter mpjWrite;
		if (stampede) {
			threads = 68;
			queue = "normal";
			remoteDir = new File("/work/00950/kevinm/stampede2/bbp/parallel");
			bbpDataDir = "/tmp";
			pbsWrite = new StampedeScriptWriter(true);
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, heapSizeMB, null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		} else {
			threads = 20;
			queue = "scec";
			remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
			bbpDataDir = "${TMPDIR}";
			pbsWrite = new USC_HPCC_ScriptWriter();
			mpjWrite = new MPJExpressShellScriptWriter(
					USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, null, USC_HPCC_ScriptWriter.MPJ_HOME);
		}
		
		String dateStr = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		mpjWrite.setClasspath(classpath);
		
		if (doGP) {
			String jobName = dateStr;
			jobName += "-"+catalogDir.getName()+"-event"+eventID+"-gp";
			
			if (src == null)
				src = BBP_SourceFile.readFile(srcFile);
			String len = src.getdLen()+"";
			if (len != null)
				jobName += "-dx"+len;
			if (!RSQSimBBP_Config.DO_HF)
				jobName += "-noHF";
			if (stampede)
				jobName += "-stampede";
			
			System.out.println("Writing GP sim to "+jobName);
			
			List<BBP_Site> sites = RSQSimBBP_Config.getStandardSites(src);
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
			
			// copy src file
			Files.copy(srcFile, new File(localJobDir, srcFile.getName()));
			File remoteSrcFile = new File(remoteJobDir, srcFile.getName());
			
			String argz = "--exact-dispatch "+threads+" --threads "+threads;
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
			
			List<String> script = mpjWrite.buildScript(MPJ_BBP_RupGenSim.class.getName(), argz);
			
			script = pbsWrite.buildScript(script, gpMins, nodes, threads, queue);
			pbsWrite.writeScript(new File(localJobDir, "gp_bbp_parallel.pbs"), script);
		}
		if (doShakeMap) {
			String jobName = dateStr;
			jobName += "-"+catalogDir.getName()+"-event"+eventID+"-shakemap";
			
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
			
			String argz = "--min-dispatch "+threads+" --threads "+threads;
			argz += " --vm "+RSQSimBBP_Config.VM.name()+" --method "+RSQSimBBP_Config.METHOD.name();
			argz += " --src-file "+remoteSrcFile.getAbsolutePath();
			argz += " --srf-file "+remoteSrfFile.getAbsolutePath();
			argz += " --output-dir "+remoteJobDir.getAbsolutePath();
			argz += " --buffer-dist "+RSQSimBBP_Config.MAX_DIST+" --spacing "+mapSpacing;
			if (!RSQSimBBP_Config.DO_HF)
				argz += " --no-hf";
			if (bbpDataDir != null && !bbpDataDir.isEmpty())
				argz += " --bbp-data-dir "+bbpDataDir;
			
			List<String> script = mpjWrite.buildScript(MPJ_BBP_ShakeMapSim.class.getName(), argz);
			
			script = pbsWrite.buildScript(script, mapMins, nodes, threads, queue);
			pbsWrite.writeScript(new File(localJobDir, "map_bbp_parallel.pbs"), script);
		}
	}

}
