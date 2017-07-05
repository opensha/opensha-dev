package scratch.peter.ucerf3.scripts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.IOUtils;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.util.ClassUtils;

import scratch.peter.ucerf3.calc.UC3_CalcUtils;
import scratch.peter.ucerf3.calc.UC3_CalcMPJ_MapCompound;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class MapsFromCompound {

	private static final String NEWLINE = IOUtils.LINE_SEPARATOR;
	private static final Joiner J_NL = Joiner.on(NEWLINE);
	private static final Joiner J_SPACE = Joiner.on(" ");
	private static final File MPJ_HOME;
	private static final File JAVA_BIN;
	
	static {
		MPJ_HOME = new File("/home/rcf-40/pmpowers/FastMPJ");
		JAVA_BIN = new File("/usr/usc/jdk/default/jre/bin/java");
	}

	public static void main(String[] args) throws IOException {
		if (args.length != 11) {
			System.out
				.println("USAGE: " +
					ClassUtils.getClassNameWithoutPackage(MapsFromCompound.class) +
					" <queue> <nodes> <hours> <libDir> <scriptPath> <branchfile>" +
					" <solPath> <grid> <spacing> <period> <outDir>");
			System.exit(1);
		}

		String queue = args[0];
		int nodes = Integer.parseInt(args[1]);
		int hours = Integer.parseInt(args[2]);
		String libDir = args[3];
		String scriptPath = args[4];
		List<String>  branches = UC3_CalcUtils.readBranchFile(args[5]);
		
		List<String> otherArgs = Lists.newArrayList(Iterables.skip(
			Arrays.asList(args), 6));

		writeScript(libDir, hours, nodes, queue, branches, scriptPath, otherArgs);
	}

	private static void writeScript(String libDir, int hrs,
			int nodes, String queue, List<String> branches, String scriptPath, 
			List<String> args) {
		try {
			File shaJAR = new File(libDir, "OpenSHA_complete.jar");
			File cliJAR = new File(libDir, "commons-cli-1.2.jar");
			ArrayList<File> classpath = Lists.newArrayList(shaJAR, cliJAR);
			FastMPJShellScriptWriter mpj = new FastMPJShellScriptWriter(
				JAVA_BIN, 4096, classpath, MPJ_HOME);

			List<String> classNames = Lists.newArrayList();
			List<String> argss = Lists.newArrayList();
			for (String branchID : branches) {
				classNames.add(UC3_CalcMPJ_MapCompound.class.getName());
				// insert branchID as second arg
				args.add(1, branchID);
				String scriptArgs = J_SPACE.join(args);
				argss.add(scriptArgs);
				args.remove(1);
			}
			List<String> script = mpj.buildScript(classNames, argss);
			script.add(NEWLINE);
			HPCC_ScriptWriter writer = new HPCC_ScriptWriter();
			script = writer.buildScript(script, hrs, nodes, 8, queue);

			File pbsFile = new File(scriptPath);
			String scriptStr = J_NL.join(script);
			Files.write(scriptStr, pbsFile, Charsets.US_ASCII);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

}
