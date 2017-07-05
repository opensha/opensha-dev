package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;

import com.google.common.collect.Lists;

public class TimeDepFSS_ERF_Simulator_ScriptGen {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/scec-02/kmilner/simulators");
		List<File> classpath = Lists.newArrayList();
		classpath.add(new File(dir, "OpenSHA_complete.jar"));
		
		MagDependentAperiodicityOptions cov = MagDependentAperiodicityOptions.LOW_VALUES;
		double duration = 50;
		String prefix_add = "cov"+cov.name()+"_dur"+(int)duration+"_";
		
		String outputDirName = "2013_11_20-erf-audit-cov-LOW-dur50";
		File localOutputDir = new File("/tmp", outputDirName);
		if (!localOutputDir.exists())
			localOutputDir.mkdir();
		File remoteOutputDir = new File(dir, outputDirName);
		
//		int numJobs = 500;
		int numJobs = 100;
//		int trialsPerJob = 5000;
//		int trialsPerJob = 2500;
		int trialsPerJob = 1000;
		int secsPerTrial = 12; // very conservative
		int mins = secsPerTrial*trialsPerJob/60;
		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		int maxHeapSizeMB = 8000;
		float hours = (float)mins/60f;
		System.out.println("Run time: "+hours);
		
		JavaShellScriptWriter writer = new JavaShellScriptWriter(javaBin, maxHeapSizeMB, classpath);
		
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		int jobDigits = (""+(numJobs-1)).length();
		
		for (int i=0; i<numJobs; i++) {
			String prefix = i+"";
			while (prefix.length() < jobDigits)
				prefix = "0"+prefix;
			prefix = "erf_audit_"+prefix_add+prefix;
			
			File pbsFile = new File(localOutputDir, prefix+".pbs");
			
			String scriptArgs = remoteOutputDir.getAbsolutePath()+" "+prefix+" "+trialsPerJob+" "+cov.name()+" "+duration;
			List<String> script = writer.buildScript(TimeDepFSS_ERF_Simulator_Test.class.getName(), scriptArgs);
			pbsWrite.writeScript(pbsFile, script, mins, 1, 8, null);
		}
	}

}
