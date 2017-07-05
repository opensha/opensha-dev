package scratch.peter.ucerf3.scripts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.IOUtils;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.util.ClassUtils;
import org.opensha.nshmp2.util.Period;

import scratch.peter.ucerf3.calc.UC3_CalcCurve;

import com.google.common.base.Charsets;
import com.google.common.base.Enums;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class CurvesFromSolution {

	private static final String NEWLINE = IOUtils.LINE_SEPARATOR;
	private static final Joiner J = Joiner.on(NEWLINE);
	private static final Splitter S = Splitter.on(',');
	private static final File JAVA_BIN;

	static {
		JAVA_BIN = new File("/usr/usc/jdk/default/jre/bin/java");
	}

	public static void main(String[] args) throws IOException {
		if (args.length != 6) {
			System.out.println("USAGE: " +
				ClassUtils.getClassNameWithoutPackage(CurvesFromSolution.class) +
				" <javaLib> <script> <filepath> <sitefile> <periods> <outDir>");
			System.exit(1);
		}

		String libDir = args[0];
		String scriptpath = args[1];
		String filepath = args[2];
		String sitefile = args[3];
		String periods = args[4];
		String outDir = args[5];

		int hours = 1;
		int nodes = 1;
		String queue = "nbns";

		writeScript(scriptpath, filepath, sitefile, periods, libDir, outDir,
			hours, nodes, queue);
	}

	private static void writeScript(String scriptpath, String filepath,
			String sitefile, String periods, String libDir, String outDir, int hrs,
			int nodes, String queue) {
		try {
			File shaJAR = new File(libDir, "OpenSHA_complete.jar");
			File cliJAR = new File(libDir, "commons-cli-1.2.jar");
			ArrayList<File> classpath = Lists.newArrayList(shaJAR, cliJAR);
			JavaShellScriptWriter jssw = new JavaShellScriptWriter(JAVA_BIN,
				5120, classpath);

			String cliArgs = filepath + " " + sitefile + " " + periods + " " + outDir;
			List<String> script = jssw.buildScript(
				UC3_CalcCurve.class.getName(), cliArgs);
			script.add(NEWLINE);
			HPCC_ScriptWriter writer = new HPCC_ScriptWriter();
			script = writer.buildScript(script, hrs, nodes, 8, queue);

			File pbsFile = new File(scriptpath);
			String scriptStr = J.join(script);
			Files.write(scriptStr, pbsFile, Charsets.US_ASCII);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

}
