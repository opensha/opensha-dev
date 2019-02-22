package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.earthquake.param.BackgroundRupType;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_UCERF3_EAL_CombinerScriptGen {

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		File remoteDir = new File("/home/scec-02/kmilner/ucerf3/eal");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		jobName += "-ucerf3-ngaw2-cea-consolidate";
		
		File trueMeanFile = new File(remoteDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		File cfssFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/2013_05_10-ucerf3p3-production-10runs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		
		File willsDir = new File(remoteDir, "2017_05_24-ucerf3-ngaw2-cea-proxy-wills2015");
		File waldDir = new File(remoteDir, "2017_05_26-ucerf3-ngaw2-cea-proxy-wald");
		BackgroundRupType bgType = BackgroundRupType.CROSSHAIR;
		// only used for tracts, can be from either directory
		File portfolioFile = new File(remoteDir, "Porter-24May2017-CA-RES1-2017-Wills2015.csv");
		
		File erfProbsDir = new File("/home/scec-02/kmilner/ucerf3/erf_probs/2014_10_07-ucerf3-erf-probs");
		double erfProbsDuration = 1d;
		
//		Location tractLoc = null;
//		String cityPrefix = null;
//		double tractRadius = 0;
		
//		Location tractLoc = new Location(34.108300, -117.289646);
//		String cityPrefix = "san-bernardino";
//		double tractRadius = 5;
		
		Location tractLoc = new Location(38.2494, -122.0400);
		String cityPrefix = "fairfield";
		double tractRadius = 15;
		
		if (tractLoc != null)
			jobName += "-"+cityPrefix;
		
		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 24*60;
		int heapSizeMB = 45*1024;
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		String argz = MPJTaskCalculator.argumentBuilder().threads(threads).endTimeSlurm().minDispatch(threads).maxDispatch(threads*10).build();
		argz += " --wills-dir "+willsDir.getAbsolutePath();
		argz += " --wald-dir "+waldDir.getAbsolutePath();
		argz += " --erf-probs-dir "+erfProbsDir.getAbsolutePath();
		argz += " --erf-probs-duration "+(float)erfProbsDuration;
		argz += " --true-mean-sol "+trueMeanFile.getAbsolutePath();
		argz += " --compound-sol "+cfssFile.getAbsolutePath();
		if (tractLoc != null) {
			argz += " --tract-location "+tractLoc.getLatitude()+","+tractLoc.getLongitude();
			if (tractRadius > 0)
				argz += " --tract-radius "+tractRadius;
			argz += " --background-type "+bgType.name();
			argz += " --portfolio "+portfolioFile.getAbsolutePath();
		}
		argz += " "+remoteJobDir.getAbsolutePath();
		
		List<String> script = mpjWrite.buildScript(MPJ_UCERF3_EAL_Combiner.class.getName(), argz);
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "eal_consolidate.slurm"), script);
	}

}
