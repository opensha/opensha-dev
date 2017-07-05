package scratch.peter.timedep;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.nshmp2.calc.HazardResultWriter;
import org.opensha.nshmp2.calc.HazardResultWriterLocal;
import org.opensha.nshmp2.calc.HazardResultWriterMPJ;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import scratch.peter.ucerf3.calc.UC3_CalcMPJ_Map;

import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.io.Closeables;
import com.google.common.io.Files;
import com.google.common.io.Flushables;

public class CalcMPJ_TimeDepMap extends MPJTaskCalculator {
	
	private static final String S = File.separator;
	
	private ThreadedTimeDepCalc calc;
	private LocationList locs;
	
	// these will only end up getting used by the dispatch (root)
	// node during doFinalAssembly(); ignored on other nodes
	private File outDir;
	private Period period;
	
	public CalcMPJ_TimeDepMap(CommandLine cmd, String[] args)
			throws IOException, InvocationTargetException, FileNotFoundException {
		
		super(cmd);
		if (args.length != 9) {
			System.err.println("USAGE: CalcMPJ_TimeDepMap [<options>] " +
					"<erfStr> <imr> <grid> <spacing> <period> <bgInclude> " +
					"<duration> <timeDep> <outPath>");
			abortAndExit(2);
		}

		Preconditions.checkArgument(getNumThreads() >= 1, 
				"threads must be >= 1. you supplied: "+getNumThreads());
		debug(rank, null, "setup for "+getNumThreads()+" threads");
		
		String erfStr = args[0]; // "UC3" | "UC2"
		SourceIMR imr = SourceIMR.valueOf(args[1]);
		TestGrid grid = TestGrid.valueOf(args[2]);
		double spacing = Double.parseDouble(args[3]);
		locs = grid.grid(spacing).getNodeList();
		period = Period.valueOf(args[4]);
		IncludeBackgroundOption bg = IncludeBackgroundOption.valueOf(args[5]);
		double duration = Double.parseDouble(args[6]);
		boolean timeDep = Boolean.parseBoolean(args[7]);
		String outPath = args[8];

		outDir = new File(outPath + S + erfStr + S + grid + S + period);
		
		HazardResultWriter writer = new HazardResultWriterMPJ(outDir);
		calc = new ThreadedTimeDepCalc(erfStr, imr, locs, period, bg, timeDep,
			duration, writer);
	}
	
	@Override
	public int getNumTasks() {
		return locs.size();
	}
	
	@Override
	public void calculateBatch(int[] batch) throws Exception, InterruptedException {
		calc.calculate(batch);
		System.out.println("Batch complete");
	}
	

	@Override
	protected void doFinalAssembly() throws Exception {
		if (rank == 0) aggregateResults(outDir, period, false);
	}
	
	
	// overridden for testing
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option erfOp = new Option("e", "mult-erfs", false,
			"If set, a copy of the ERF will be instantiated for each thread.");
		erfOp.setRequired(false);
		ops.addOption(erfOp);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		try {
			Options options = createOptions();
			CommandLine cmd = parse(options, args, CalcMPJ_TimeDepMap.class);
			args = cmd.getArgs();
			CalcMPJ_TimeDepMap driver = new CalcMPJ_TimeDepMap(cmd, args);
			driver.run();
			finalizeMPJ();
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}
	
	/**
	 * Utility method to aggregate hazard curves stored in indivudual files
	 * with names lat_lon.txt.
	 * 
	 * @param dir containing curve files
	 * @param period for which curves were calculated
	 */
	public static void aggregateResults(File dir, Period period, boolean determ) {
		if (determ) {
			try {
				File[] pFiles = dir.listFiles(new FileFilter() {
					@Override public boolean accept(File f) {
						return f.getName().endsWith("prob.txt");
					}
				});
				File curves = new File(dir, "curves.csv");
				BufferedWriter brP = Files.newWriter(curves, Charsets.US_ASCII);
				HazardResultWriterLocal.writeCurveHeader(brP, period);
				for (File file : pFiles) {
					StringBuilder sb = new StringBuilder();
					String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
						file.getName(), "_"), '_', ',');
					sb.append(latlon).append(",");
					Files.copy(file, Charsets.US_ASCII, sb);
					brP.write(sb.toString());
					brP.newLine();
					file.delete();
				}
				Flushables.flushQuietly(brP);
				Closeables.close(brP, true);

				File[] dFiles = dir.listFiles(new FileFilter() {
					@Override public boolean accept(File f) {
						return f.getName().endsWith("det.txt");
					}
				});
				File detVals = new File(dir, "determ.txt");
				BufferedWriter brD = Files.newWriter(detVals, Charsets.US_ASCII);
				for (File file : dFiles) {
					StringBuilder sb = new StringBuilder();
					String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
						file.getName(), "_"), '_', '\t');
					sb.append(latlon).append("\t");
					Files.copy(file, Charsets.US_ASCII, sb);
					brD.write(sb.toString());
					brD.newLine();
					file.delete();
				}
				Flushables.flushQuietly(brD);
				Closeables.close(brD, true);

			} catch (Exception e) {
				e.printStackTrace();
			}
			
		} else {
			String[] exts = {"txt"};
			try {
				Collection<File> files = FileUtils.listFiles(dir, exts, false);
				File curves = new File(dir, "curves.csv");
				BufferedWriter br = Files.newWriter(curves, Charsets.US_ASCII);
				HazardResultWriterLocal.writeCurveHeader(br, period);
				for (File file : files) {
					StringBuilder sb = new StringBuilder();
					String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
						file.getName(), "."), '_', ',');
					sb.append(latlon).append(",");
					Files.copy(file, Charsets.US_ASCII, sb);
					br.write(sb.toString());
					br.newLine();
					file.delete();
				}
				Flushables.flushQuietly(br);
				Closeables.close(br, true);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
}
