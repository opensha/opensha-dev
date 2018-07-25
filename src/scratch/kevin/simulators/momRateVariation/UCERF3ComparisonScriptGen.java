package scratch.kevin.simulators.momRateVariation;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class UCERF3ComparisonScriptGen {

	public static void main(String[] args) throws IOException {
		int numBatches =  1;
		int threadsPerBatch = 20;
		int catalogDuration = 1000000;
		int memGB = 46;
		String queue = "scec";
		int timeHours = 72;
		MagDependentAperiodicityOptions cov = MagDependentAperiodicityOptions.MID_VALUES;
		String nameAdd = "fm3_1-geol";
		
		String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date())+"-"+cov.name();
		if (nameAdd != null)
			jobName += "-"+nameAdd;
//		File localDir = new File("/home/kevin/Simulators/time_series/ucerf3_compare");
//		File remoteDir = new File("/auto/scec-02/kmilner/ucerf3/synth_catalog/mom_rate_tests");
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_catalogs");
		File remoteDir = new File("/auto/scec-02/kmilner/ucerf3/time_dep_catalogs");
//		File fssFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/"
//				+ "2013_05_10-ucerf3p3-production-10runs/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_"
//				+ "FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		File fssFile = new File("/auto/scec-02/kmilner/ucerf3/inversion_compound_plots/"
				+ "2013_05_10-ucerf3p3-production-10runs_fm_dm_sub_plots/FM3_1_GEOL/"
				+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
		
		
		File localJobDir = new File(localDir, jobName);
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		List<File> classpath = Lists.newArrayList(new File(remoteJobDir, "opensha-dev-all.jar"));
		
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, memGB*1024, classpath);
		
		int mins = timeHours*60;
		
		for (int batch=0; batch<numBatches; batch++) {
			String dirName = "batch"+batch;
			File localOutputDir = new File(localJobDir, dirName);
			Preconditions.checkState(localOutputDir.exists() || localOutputDir.mkdir());
			File remoteOutputDir = new File(remoteJobDir, dirName);
			
			String argsStr = remoteOutputDir.getAbsolutePath()+" "+fssFile.getAbsolutePath()+" "+threadsPerBatch
					+" "+catalogDuration+" "+cov.name();
			
			List<String> script = javaWrite.buildScript(UCERF3ComparisonCalc.class.getName(), argsStr);
			File outputFile = new File(localOutputDir, dirName+".pbs");
			pbsWrite.writeScript(outputFile, script, mins, 1, threadsPerBatch, queue);
		}
	}

}
