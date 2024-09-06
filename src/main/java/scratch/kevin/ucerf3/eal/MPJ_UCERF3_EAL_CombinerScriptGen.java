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
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.earthquake.param.BackgroundRupType;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_UCERF3_EAL_CombinerScriptGen {

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		File remoteDir = new File("/project/scec_608/kmilner/ucerf3/eal");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		File trueMeanFile = new File(remoteDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		File cfssFile = new File("/project/scec_608/kmilner/ucerf3/inversion_compound_plots/2013_05_10-ucerf3p3-production-10runs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		
//		File willsDir = new File(remoteDir, "2017_05_24-ucerf3-ngaw2-cea-proxy-wills2015");
//		File waldDir = new File(remoteDir, "2017_05_26-ucerf3-ngaw2-cea-proxy-wald");
//		jobName += "-ucerf3-ngaw2-cea-consolidate";
		// only used for tracts, can be from either directory
//		File portfolioFile = new File(remoteDir, "Porter-24May2017-CA-RES1-2017-Wills2015.csv");
		
		File willsDir = new File(remoteDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wills2015");
		File waldDir = new File(remoteDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wald2007");
		jobName += "-ucerf3-ngaw2-cea-100pct-consolidate";
		// only used for tracts, can be from either directory
		File portfolioFile = new File(remoteDir, "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wills2015.csv");
		
		BackgroundRupType bgType = BackgroundRupType.CROSSHAIR;
		
		File erfProbsDir = new File("/project/scec_608/kmilner/ucerf3/erf_probs/2014_10_07-ucerf3-erf-probs");
		double erfProbsDuration = 1d;
		
		Location tractLoc = null;
		String cityPrefix = null;
		double tractRadius = 0;
		boolean tractIndividual = false;
		boolean lec = true;
//		LossCOV_Model covModel = LossCOV_Model.PORTER_POWER_LAW_2020_09_01;
		LossCOV_Model covModel = LossCOV_Model.PORTER_POWER_LAW_2020_09_01_fixed;
		
//		Location tractLoc = new Location(34.108300, -117.289646);
//		String cityPrefix = "san-bernardino";
//		double tractRadius = 5;
		
//		Location tractLoc = new Location(38.2494, -122.0400);
//		String cityPrefix = "fairfield";
//		double tractRadius = 15;
		
//		Location tractLoc = new Location(34.0407, -118.2468);
//		String cityPrefix = "los-angeles";
//		double tractRadius = 100;
		
//		Location tractLoc = new Location(37.7946, -122.3999);
//		String cityPrefix = "san-francisco";
//		double tractRadius = 100;
		
		if (tractLoc != null)
			jobName += "-"+cityPrefix;
		
		if (tractIndividual)
			jobName += "-tractEALs";
		
		if (lec) {
			jobName += "-calcLEC";
			if (covModel != null)
				jobName += "-covModel";
		}
		
		int threads = 20;
		int nodes = 32;
		String queue = "scec";
		int mins = 24*60;
		if (tractIndividual)
			mins = mins*7;
		int heapSizeMB = 55*1024;
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_CARC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		int minDispatch = tractIndividual ? 1 : threads;
		int maxDispatch = tractIndividual ? 10 : threads*10;
		
		String argz = MPJTaskCalculator.argumentBuilder().threads(threads).endTimeSlurm()
					.minDispatch(minDispatch).maxDispatch(maxDispatch).build();
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
			argz += " --portfolio "+portfolioFile.getAbsolutePath();
			argz += " --node-temp-dir ${TMPDIR}";
		}
		if (tractIndividual) {
			argz += " --background-type "+bgType.name();
			argz += " --portfolio "+portfolioFile.getAbsolutePath();
			argz += " --tract-branch-eals";
		} else if (lec) {
			argz += " --calc-lec";
			if (covModel != null)
				argz += " --lec-cov "+covModel.name();
		}
		argz += " "+remoteJobDir.getAbsolutePath();
		
		List<String> script = mpjWrite.buildScript(MPJ_UCERF3_EAL_Combiner.class.getName(), argz);
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "eal_consolidate.slurm"), script);
	}

}
