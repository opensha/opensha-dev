package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

public class MPJ_HPC_DebugScriptGen {

	public static void main(String[] args) throws IOException {
		int numTasks = 1000;
		int sleepSeconds = 1;
		
		int memGigs = 24;
		
		int mins = 60;
		int ppn = 8;
		String queue = "scec";
		
		int[] nodeCounts = { 2, 4, 8, 16, 24, 36 };
		boolean[] fmpjs = { false, true };
		
		File remoteDir = new File("/home/scec-02/kmilner/mpj_debug");
		File localDir = new File("/tmp");
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteDir, "commons-cli-1.2.jar"));
		classpath.add(new File(remoteDir, "OpenSHA_complete.jar"));
		FastMPJShellScriptWriter fmpjWrite = new FastMPJShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
				classpath, USC_HPCC_ScriptWriter.FMPJ_HOME);
		fmpjWrite.setUseLaunchWrapper(true);
		MPJExpressShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
				USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
				classpath, USC_HPCC_ScriptWriter.MPJ_HOME);
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		String className = MPJ_DummyTest.class.getName();
		String argz = "--num "+numTasks+" --time "+sleepSeconds;
		
		int nodeDigits = 0;
		for (int nodes : nodeCounts) {
			int len = (nodes+"").length();
			if (len > nodeDigits)
				nodeDigits = len;
		}
		
		for (int nodes : nodeCounts) {
			for (boolean fmpj : fmpjs) {
				String nodeStr = nodes+"";
				while (nodeStr.length() < nodeDigits)
					nodeStr = "0"+nodeStr;
				String prefix;
				if (fmpj)
					prefix = "fmpj_"+nodeStr+"nodes";
				else
					prefix = "mpj_"+nodeStr+"nodes";
				
				System.out.println("Writing for "+prefix);
				
				File statusFile = new File(remoteDir, prefix+"_status.txt");
				String myArgs = argz+" --file "+statusFile.getAbsolutePath();
				
				List<String> script;
				if (fmpj)
					script = fmpjWrite.buildScript(className, myArgs);
				else
					script = mpjWrite.buildScript(className, myArgs);
				script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
				
				pbsWrite.writeScript(new File(localDir, prefix+".pbs"), script);
			}
		}
	}

}
