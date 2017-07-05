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

import scratch.peter.ucerf3.calc.UC3_CalcMPJ_Map;
import scratch.peter.ucerf3.calc.UC3_CalcMPJ_MapAverage;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;
import scratch.peter.ucerf3.calc.UC3_CalcMPJ_MapCompound;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class MapsFromSolution {

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
		if (args.length != 13) {
			System.out
				.println("USAGE: " +
					ClassUtils.getClassNameWithoutPackage(MapsFromSolution.class) +
					" <queue> <nodes> <hours> <libDir> <scriptPath>" +
					" <solPath> <imr> <grid> <spacing> <period> <bgOption>" +
					" <outDir> <nshmp>");
			System.exit(1);
		}

		String queue = args[0];
		int nodes = Integer.parseInt(args[1]);
		int hours = Integer.parseInt(args[2]);
		String libDir = args[3];
		String scriptPath = args[4];
		
		List<String> otherArgs = Lists.newArrayList(Iterables.skip(
			Arrays.asList(args), 5));

		writeScript(libDir, hours, nodes, queue, scriptPath, otherArgs);
	}

	public static void writeScript(String libDir, int hrs,
			int nodes, String queue, String scriptPath, 
			List<String> args) {
		try {
			File shaJAR = new File(libDir, "OpenSHA_complete.jar");
			File cliJAR = new File(libDir, "commons-cli-1.2.jar");
			ArrayList<File> classpath = Lists.newArrayList(shaJAR, cliJAR);
			FastMPJShellScriptWriter mpj = new FastMPJShellScriptWriter(
				JAVA_BIN, 4096, classpath, MPJ_HOME);

			String scriptArgs = J_SPACE.join(args);
			List<String> script = mpj.buildScript(
				UC3_CalcMPJ_Map.class.getName(), scriptArgs);
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
