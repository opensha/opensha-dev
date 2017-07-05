package scratch.peter.ucerf3.calc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardResultWriter;
import org.opensha.nshmp2.calc.ThreadedHazardCalc;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EpistemicListERF;

import scratch.peter.nshmp.HazardResultWriterMPJ_NSHMP_Det;
import scratch.peter.ucerf3.NSHMP13_DeterminisiticERF;

import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.io.Closeables;
import com.google.common.io.Files;
import com.google.common.io.Flushables;

public class UC3_CalcMPJ_AltDetermMap extends MPJTaskCalculator {
	
	private static final String S = File.separator;
	private ThreadedHazardCalc calc;
	private LocationList locs;
	private HazardResultWriterMPJ_NSHMP_Det writer;
	
	// these will only end up getting used by the dispatch (root)
	// node during doFinalAssembly(); ignored on other nodes
	private File outDir;
	private Period period;
	
	public UC3_CalcMPJ_AltDetermMap(CommandLine cmd, String[] args)
			throws IOException, InvocationTargetException, FileNotFoundException {
		
		super(cmd);
		if (args.length != 6) {
			System.err.println("ARGS: " + Arrays.toString(args));
			System.err.println("USAGE: UC3_CalcMPJ_AltDetermMap [<options>] " +
					"<imr> <grid> <spacing> <period> <outPath> <alea>");
			abortAndExit(2);
		}

		Preconditions.checkArgument(getNumThreads() >= 1, 
				"threads must be >= 1. you supplied: "+getNumThreads());
		debug(rank, null, "setup for "+getNumThreads()+" threads");
		
		SourceIMR imr = SourceIMR.valueOf(args[0]);
		TestGrid grid = TestGrid.valueOf(args[1]);
		double spacing = Double.parseDouble(args[2]);
		locs = grid.grid(spacing).getNodeList();
		period = Period.valueOf(args[3]);
		String outPath = args[4];
		boolean alea = Boolean.parseBoolean(args[5]);

		outDir = new File(outPath + S + grid + S + period);
		outDir.mkdirs();
		
		// only output deterministic
		writer = new HazardResultWriterMPJ_NSHMP_Det();
		AbstractERF erf = NSHMP13_DeterminisiticERF.create(alea);
		EpistemicListERF wrapped = ERF_ID.wrapInList(erf);
		wrapped.updateForecast();
		
		calc = new ThreadedHazardCalc(wrapped, imr, locs, period, false, writer, true);
	}
	
	@Override
	public int getNumTasks() {
		return locs.size();
	}
	
	@Override
	public void calculateBatch(int[] batch) throws Exception, InterruptedException {
		calc.calculate(batch);
	}
	

	@Override
	protected void doFinalAssembly() throws Exception {
//		if (rank == 0) aggregateResults(outDir, period, true);
		writer.toFile(outDir, rank);
		// will combine per-node files independently
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
//		TestGrid gr = TestGrid.CA_RELM;
//		gr.grid(0.1);
//		RegionUtils.locListToKML(gr.grid(0.1).getBorder(), "CA_RELM", Color.magenta);
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			CommandLine cmd = parse(options, args, UC3_CalcMPJ_AltDetermMap.class);
			args = cmd.getArgs();
			UC3_CalcMPJ_AltDetermMap driver = new UC3_CalcMPJ_AltDetermMap(cmd, args);
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
//		if (determ) {
			try {
//				File[] pFiles = dir.listFiles(new FileFilter() {
//					@Override public boolean accept(File f) {
//						return f.getName().endsWith("prob.txt");
//					}
//				});
//				File curves = new File(dir, "curves.csv");
//				BufferedWriter brP = Files.newWriter(curves, Charsets.US_ASCII);
//				HazardResultWriterLocal.writeCurveHeader(brP, period);
//				for (File file : pFiles) {
//					StringBuilder sb = new StringBuilder();
//					String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
//						file.getName(), "_"), '_', ',');
//					sb.append(latlon).append(",");
//					Files.copy(file, Charsets.US_ASCII, sb);
//					brP.write(sb.toString());
//					brP.newLine();
//					file.delete();
//				}
//				Flushables.flushQuietly(brP);
//				Closeables.closeQuietly(brP);

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
			
//		} else {
//			String[] exts = {"txt"};
//			try {
//				Collection<File> files = FileUtils.listFiles(dir, exts, false);
//				File curves = new File(dir, "curves.csv");
//				BufferedWriter br = Files.newWriter(curves, Charsets.US_ASCII);
//				HazardResultWriterLocal.writeCurveHeader(br, period);
//				for (File file : files) {
//					StringBuilder sb = new StringBuilder();
//					String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
//						file.getName(), "."), '_', ',');
//					sb.append(latlon).append(",");
//					Files.copy(file, Charsets.US_ASCII, sb);
//					br.write(sb.toString());
//					br.newLine();
//					file.delete();
//				}
//				Flushables.flushQuietly(br);
//				Closeables.closeQuietly(br);
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//		}
	}
	
	private static String nameFromPath(String solPath) {
		int ssIdx1 = StringUtils.lastIndexOf(solPath, "/") + 1;
		int ssIdx2 = StringUtils.lastIndexOf(solPath, ".");
		return solPath.substring(ssIdx1, ssIdx2);
	}

}
