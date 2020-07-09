package scratch.kevin.ucerf3.eal.spatialCorr;

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

public class MPJ_SpatiallyCorrelatedLossCalcScriptGen {

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		File remoteDir = new File("/home/scec-02/kmilner/ucerf3/eal");
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		File trueMeanFile = new File(remoteDir, "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		
//		File willsDir = new File(remoteDir, "2017_05_24-ucerf3-ngaw2-cea-proxy-wills2015");
//		File waldDir = new File(remoteDir, "2017_05_26-ucerf3-ngaw2-cea-proxy-wald");
//		jobName += "-ucerf3-ngaw2-cea-consolidate";
		// only used for tracts, can be from either directory
//		File portfolioFile = new File(remoteDir, "Porter-24May2017-CA-RES1-2017-Wills2015.csv");
		
		File willsDir = new File(remoteDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wills2015");
		File waldDir = new File(remoteDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wald2007");
		jobName += "-ucerf3-ngaw2-cea-100pct";
		// only used for tracts, can be from either directory
		File portfolioFile = new File(remoteDir, "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wills2015.csv");
		File vulnFile = new File(remoteDir, "2014_05_16_VUL06.txt");
		
		String taus = "-2,-1,0,1,2";
		File randFieldsDir = new File("/home/scec-02/kmilner/ucerf3/eal/random_fields/sa10_1km_800x800");
		double fieldSpacing = 1d;
		
		int threads = 20;
		int nodes = 10;
		String queue = "scec";
		int mins = 24*60;
		int heapSizeMB = 55*1024;
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		jobName += "-spatial";
		
		File localJobDir = new File(localDir, jobName);
		System.out.println(localJobDir.getAbsolutePath());
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		int minDispatch = 1;
		int maxDispatch = 5;
		
		String argz = MPJTaskCalculator.argumentBuilder().threads(threads).endTimeSlurm()
					.minDispatch(minDispatch).maxDispatch(maxDispatch).build();
		argz += " --wills-dir "+willsDir.getAbsolutePath();
		argz += " --wald-dir "+waldDir.getAbsolutePath();
		argz += " --true-mean-sol "+trueMeanFile.getAbsolutePath();
		argz += " --portfolio "+portfolioFile.getAbsolutePath();
		argz += " --vuln-file "+vulnFile.getAbsolutePath();
		argz += " --taus "+taus;
		argz += " --fields-dir "+randFieldsDir.getAbsolutePath();
		argz += " --field-spacing "+fieldSpacing;
		
		argz += " "+remoteJobDir.getAbsolutePath();
		
		List<String> script = mpjWrite.buildScript(MPJ_SpatiallyCorrelatedLossCalc.class.getName(), argz);
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(new File(localJobDir, "spatial_calc.slurm"), script);
	}

}
