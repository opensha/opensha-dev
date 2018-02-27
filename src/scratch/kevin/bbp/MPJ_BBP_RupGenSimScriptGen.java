package scratch.kevin.bbp;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class MPJ_BBP_RupGenSimScriptGen {

	public static void main(String[] args) throws IOException {
		VelocityModel vm = VelocityModel.LA_BASIN; // 1 on HPC
		Method method = Method.GP;
		boolean noHF = true;
//		File catalogDir = new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement");
//		int eventID = 4099020;
		File catalogDir = new File("/data/kevin/simulators/catalogs/rundir2194_long");
		int eventID = 136704;
		
		int numSites = 2;
		File sitesFile = new File("/home/kevin/bbp/bbp_data/run/stations_cs_sites.stl");
		File eventsDir = new File(catalogDir, "event_srfs");
		File srcFile = new File(eventsDir, "event_"+eventID+".src");
		Preconditions.checkState(srcFile.exists(), "Source file doesn't exist: %s", srcFile.getAbsolutePath());
		
		int threads = 20;
		int nodes = 20;
		String queue = "scec";
		int mins = 10*60;
		boolean splitSites = true;
		String bbpDataDir = "${TMPDIR}";
		
		int numSims = threads*nodes;
		if (splitSites)
			numSims /= numSites;
		
		File localDir = new File("/home/kevin/bbp/parallel");
		File remoteDir = new File("/auto/scec-02/kmilner/bbp/parallel");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-"+catalogDir.getName()+"-event"+eventID;
		String len = detectLenScale(srcFile);
		if (len != null)
			jobName += "-dx"+len;
		if (noHF)
			jobName += "-noHF";
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		// copy sites file
		Files.copy(sitesFile, new File(localJobDir, sitesFile.getName()));
		File remoteSitesFile = new File(remoteJobDir, sitesFile.getName());
		
		// copy src file
		Files.copy(srcFile, new File(localJobDir, srcFile.getName()));
		File remoteSrcFile = new File(remoteJobDir, srcFile.getName());
		
		String argz = MPJTaskCalculator.argumentBuilder().exactDispatch(threads).threads(threads).endTimeSlurm().build();
		argz += " --vm "+vm.name()+" --method "+method.name();
		argz += " --sites-file "+remoteSitesFile.getAbsolutePath();
		argz += " --src-file "+remoteSrcFile.getAbsolutePath();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --num-sims "+numSims;
		if (splitSites)
			argz += " --split-sites";
		if (noHF)
			argz += " --no-hf";
		if (bbpDataDir != null && !bbpDataDir.isEmpty())
			argz += " --bbp-data-dir "+bbpDataDir;
		
		int heapSizeMB = 10*1024;
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		MPJExpressShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		List<String> script = mpjWrite.buildScript(MPJ_BBP_RupGenSim.class.getName(), argz);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "bbp_parallel.pbs"), script);
	}
	
	private static String detectLenScale(File srcFile) throws IOException {
		for (String line : Files.readLines(srcFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("DLEN")) {
				String[] split = line.split("\\s+");
				return split[split.length-1];
			}
		}
		return null;
	}

}
