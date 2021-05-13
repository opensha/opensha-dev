package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

public class HazardBranchesPostProcessPBSWriter {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File(args[0]);
		File outDir = new File(dir, "assemble");
		if (!outDir.exists())
			outDir.mkdir();
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(new File("/home/scec-02/kmilner/hazMaps/svn/dist/OpenSHA_complete.jar"));
		
		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, 7000, classpath);
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter("dodecacore");
		
		int mins = 60;
		int nodes = 1;
		int ppn = 8;
		String queue = "nbns";
		
		for (File curveDir : dir.listFiles()) {
			if (!curveDir.isDirectory())
				continue;
			String name = curveDir.getName();
			if (!name.startsWith("curves"))
				continue;
			
			String cliargs = curveDir.getAbsolutePath();
			
			List<String> script = javaWrite.buildScript(HazardBranchesPostProcess.class.getName(), cliargs);
			File pbsFile = new File(outDir, "assemble_"+name+".pbs");
			pbsWrite.writeScript(pbsFile, script, mins, nodes, ppn, queue);
		}
	}

}
