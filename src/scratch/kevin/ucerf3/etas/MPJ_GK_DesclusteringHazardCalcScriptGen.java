package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import avro.shaded.com.google.common.base.Preconditions;
import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_GK_DesclusteringHazardCalcScriptGen {

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF3/etas/etas_decluster");
		Preconditions.checkState(localMainDir.exists() || localMainDir.mkdir());
		File remoteMainDir = new File("/home/scec-06/kmilner/ucerf3/etas_decluster");
		
		File catalogFile = new File("/home/scec-06/kmilner/ucerf3/etas_sim/"
				+ "2019_11_05-Start2012_500yr_kCOV1p5_Spontaneous_HistoricalCatalog/results_m5_preserve_chain.bin");
		File fssFile = new File("/home/scec-06/kmilner/ucerf3/ucerf3-etas-launcher/inputs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		double spacing = 0.1d;
		
		String dirName = "2020_12_11-decluster-full-td-kCOV1.5-scale1.14-spacing"+(float)spacing;

		int threads = 20;
		int nodes = 36;
		String queue = "scec";
		int mins = 20*60;
		int heapSizeMB = 45*1024;
		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(remoteMainDir, "opensha-dev-all.jar"));
		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, heapSizeMB, classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		File localJobDir = new File(localMainDir, dirName);
		if (!localJobDir.exists())
			localJobDir.mkdir();
		File remoteJobDir = new File(remoteMainDir, dirName);
		
		File pbsFile = new File(localJobDir, "declustered_hazard.slurm");
		
		String argz = MPJTaskCalculator.argumentBuilder().threads(threads).minDispatch(threads).maxDispatch(10*threads).build(" ");
		argz += " --fss-file "+fssFile.getAbsolutePath();
		argz += " --catalog-file "+catalogFile.getAbsolutePath();
		argz += " --output-dir "+new File(remoteJobDir, "results");
		argz += " --grid-spacing "+spacing;
		
		List<String> script = mpjWrite.buildScript(MPJ_GK_DesclusteringHazardCalc.class.getName(), argz);
		
		script = pbsWrite.buildScript(script, mins, nodes, threads, queue);
		pbsWrite.writeScript(pbsFile, script);
	}

}
