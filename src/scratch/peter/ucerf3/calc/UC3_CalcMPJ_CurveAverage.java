package scratch.peter.ucerf3.calc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.calc.HazardResultWriterSites;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

import com.google.common.base.Enums;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Curves at sites for an average solution.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC3_CalcMPJ_CurveAverage extends MPJTaskCalculator {

	private static final String S = File.separator;
	private static final Splitter SPLIT = Splitter.on(',');

	private String solPath;
	private Map<String, Location> locMap;
	private LocationList locs;
	private int solCount;
	private List<Period> periods;
	private IncludeBackgroundOption bg;
	private String outDir;

	private boolean epiUncert = false;

	private ExecutorService ex;
	private ExecutorCompletionService<HazardResult> ecs;

	public UC3_CalcMPJ_CurveAverage(CommandLine cmd, String[] args)
			throws IOException, InvocationTargetException,
			FileNotFoundException {

		super(cmd);
		if (args.length < 6) {
			System.err.println("USAGE: UC3_CalcMPJ_CurveAverage [<options>] "
				+ "<solfile> <sitefile> <solCount> <periods> <bgOption> <outDir>");
			abortAndExit(2);
		}

		Preconditions.checkArgument(getNumThreads() >= 1,
			"threads must be >= 1. you supplied: " + getNumThreads());
		debug(rank, null, "setup for " + getNumThreads() + " threads");

		// read args
		solPath = args[0];
		locMap = UC3_CalcUtils.readSiteFile(args[1]);
		locs = new LocationList();
		for (Location loc : locMap.values()) {
			locs.add(loc);
		}
		solCount = Integer.parseInt(args[2]);
		periods = readArgAsList(args[3], Period.class);
		bg = IncludeBackgroundOption.valueOf(args[4]);
		outDir = args[5];

		// init executor
		int numProc = Runtime.getRuntime().availableProcessors();
		ex = Executors.newFixedThreadPool(numProc);
		ecs = new ExecutorCompletionService<HazardResult>(ex);
	}

	@Override
	public int getNumTasks() {
		return solCount;
	}

	@Override
	public void calculateBatch(int[] indices) throws InterruptedException,
			ExecutionException, IOException {

		for (int idx : indices) {

			// init erf for branch
			FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solPath,
				idx, bg, false,
				true, 1.0);
			erf.updateForecast();
			EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);

			String erfOutDir = outDir + S + erf.getName();
			HazardResultWriterSites writer = new HazardResultWriterSites(
				erfOutDir, locMap);

			// loop periods and locations
			for (Period period : periods) {
				writer.writeHeader(period);
				for (Location loc : locs) {
					Site site = new Site(loc);
					HazardCalc hc = HazardCalc.create(wrappedERF, site, period,
						epiUncert);
					ecs.submit(hc);
				}
			}

			// collect period-location pair results
			int resultCount = locs.size() * periods.size();
			for (int j = 0; j < resultCount; j++) {
				writer.write(ecs.take().get());
			}

		}
		System.out.println("Batch complete");
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		if (ex != null) ex.shutdown();
	}

	private static <T extends Enum<T>> List<T> readArgAsList(String arg,
			Class<T> clazz) {
		Iterable<T> it = Iterables.transform(SPLIT.split(arg),
			Enums.stringConverter(clazz));
		return Lists.newArrayList(it);
	}

	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);

		try {
			Options options = createOptions();
			CommandLine cmd = parse(options, args,
				UC3_CalcMPJ_CurveAverage.class);
			args = cmd.getArgs();
			UC3_CalcMPJ_CurveAverage driver = new UC3_CalcMPJ_CurveAverage(cmd,
				args);
			driver.run();
			finalizeMPJ();
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
