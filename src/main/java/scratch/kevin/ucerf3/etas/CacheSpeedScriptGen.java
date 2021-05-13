package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

public class CacheSpeedScriptGen {

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/etas/cache_tests/soft_10Gjvm");
		File remoteDir = new File("/auto/scec-02/kmilner/ucerf3/etas_sim/cache_test/soft_10Gjvm");
		
		double[] sizes = { 0d, 0.5d, 1d, 2d, 4d, 8d, 16d };
		int numRuns = 5;
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(new File(remoteDir, "OpenSHA_complete.jar"));
		classpath.add(new File(remoteDir.getParentFile().getParentFile(), "commons-cli-1.2.jar"));
		
		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, 10*1024, classpath);
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter("dodecacore");
		
		for (double size : sizes) {
			for (int run=0; run<numRuns; run++) {
				String jobName = "cache_"+(float)size+"gb_run"+run+".pbs";
				
				List<String> script = javaWrite.buildScript(CacheSpeedTester.class.getName(),
						remoteDir.getAbsolutePath()+" "+(float)size+" "+run);
				script.add(script.size()-2, "cd "+remoteDir.getParentFile().getAbsolutePath());
				pbsWrite.writeScript(new File(localDir, jobName), script, 90, 1, 24, null);
			}
		}
	}

}
