package scratch.kevin.ucerf3.etas;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CurveMetadata;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import mpi.MPI;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.etas.ETAS_HazardMapCalc.Duration;
import scratch.kevin.ucerf3.etas.ETAS_HazardMapCalc.DurationConstants;
import scratch.kevin.ucerf3.etas.ETAS_HazardMapCalc.MapType;

public class MPJ_ETAS_HazardMapCalc extends MPJTaskCalculator {
	
	private List<List<ETAS_EqkRupture>> catalogs;
	
	// for precomputed shakemaps
	private RandomAccessFile raFile;
	private long[] filePositions;
	private int[] fileLengths;
	
	private GriddedRegion region;
	
	private String imt;
	private double period;
	
	private ETAS_HazardMapCalc mapCalc;
	
	private AttenRelRef gmpeRef;
	
	private List<Site> sites;
	
	private boolean calcGridded;
	private boolean calcFault;
	private boolean calcLongTerm;
	
	private MapType[] mapTypes;
	private Duration[] durations;
	private static final Duration[] DURATIONS_DEFAULT = ETAS_HazardMapCalc.getDurationDefaults();
	
	private final boolean griddedConditional = true;
	private ETAS_CatalogGridSourceProvider griddedSources;
	
	private ExecutorService executor;
	
	private String imtName;
	private BinaryCurveArchiver archiver;
	private DiscretizedFunc xVals;
	
	private boolean printEach;

	public MPJ_ETAS_HazardMapCalc(CommandLine cmd) throws IOException, DocumentException {
		super(cmd);
		
		File catalogsFile = new File(cmd.getOptionValue("catalogs"));
		catalogs = ETAS_CatalogIO.loadCatalogsBinary(catalogsFile);
		
		FaultSystemSolution sol = null;
		if (cmd.hasOption("fault-data-file")) {
			// precalc mode
			File faultDataFile = new File(cmd.getOptionValue("fault-data-file"));
			Preconditions.checkState(faultDataFile.exists());
			raFile = new RandomAccessFile(faultDataFile, "r");
		} else {
			Preconditions.checkArgument(cmd.hasOption("solution-file"),
					"Must supply fault system solution file if no fault data precalc file");
			File solFile = new File(cmd.getOptionValue("solution-file"));
			Preconditions.checkState(solFile.exists());
			sol = FaultSystemIO.loadSol(solFile);
		}
		
		double spacing = Double.parseDouble(cmd.getOptionValue("spacing"));
		region = new CaliforniaRegions.RELM_TESTING_GRIDDED(spacing);
		
		calcGridded = !cmd.hasOption("no-gridded");
		calcFault = !cmd.hasOption("no-fault");
		calcLongTerm = cmd.hasOption("calc-long-term");
		
		Preconditions.checkArgument(calcGridded || calcFault || calcLongTerm);
		
		imt = cmd.getOptionValue("imt").toUpperCase();
		if (imt.equals(SA_Param.NAME)) {
			Preconditions.checkArgument(cmd.hasOption("period"), "Must supply period if Sa");
			period = Double.parseDouble(cmd.getOptionValue("period"));
			if (period == Math.floor(period))
				imtName = "sa_"+(int)period+"s";
			else
				imtName = "sa_"+(float)period+"s";
		} else {
			imtName = imt.toLowerCase();
		}
		
		File outputDir = new File(cmd.getOptionValue("output-dir"));
		if (rank == 0)
			Preconditions.checkState(outputDir.exists() && outputDir.isDirectory() || outputDir.mkdir(),
				"Output directory doesn't exist or couldn't be created: %s", outputDir.getAbsoluteFile());
		
		gmpeRef = AttenRelRef.valueOf(cmd.getOptionValue("gmpe"));
		
		xVals = new IMT_Info().getDefaultHazardCurve(imt);
		
		debug("Loading sites");
		File sitesFile = new File(cmd.getOptionValue("sites-file"));
		sites = ETAS_HazardMapCalc.loadSitesFile(gmpeRef, sitesFile);
		Preconditions.checkState(sites.size() == region.getNodeCount(), "Supplied sites file is wrong size");
		
		double griddedSpacing = Double.parseDouble(cmd.getOptionValue("gridded-spacing"));
		
		if (calcGridded)
			griddedSources = new ETAS_CatalogGridSourceProvider(catalogs, griddedSpacing, griddedConditional);
		
		if (cmd.hasOption("durations")) {
			String durStr = cmd.getOptionValue("durations");
			durations = parseDurations(durStr);
		} else {
			durations = DURATIONS_DEFAULT;
		}
		
		// don't give it the fault file, we're reading externally
		mapCalc = new ETAS_HazardMapCalc(catalogs, region, xVals, null, sol, griddedSources, gmpeRef, imt, spacing, sites, durations);
		mapCalc.setCalcFaults(calcFault);
		mapCalc.setCalcGridded(calcGridded);
		mapCalc.setCalcLongTerm(calcLongTerm);
		
		if (calcLongTerm && cmd.hasOption("long-term-durations")) {
			String durStr = cmd.getOptionValue("long-term-durations");
			Duration[] longTermDurations = parseDurations(durStr);
			mapCalc.setLongTermCalcDurations(longTermDurations);
		}
		
		if (calcLongTerm && cmd.hasOption("elastic-rebound")) {
			TestScenario scenario = TestScenario.valueOf(cmd.getOptionValue("elastic-rebound"));
			if (rank == 0)
				debug("Applying elastic rebound for scenario: "+scenario.name());
			mapCalc.setScenarioForElasticRebound(scenario);
		}
		
		executor = mapCalc.createExecutor(getNumThreads());
		
		if (calcFault && raFile != null)
			loadFilePositions();
		
		List<MapType> mapTypesList = Lists.newArrayList();
		if (calcFault)
			mapTypesList.add(MapType.FAULT_ONLY);
		if (calcGridded)
			mapTypesList.add(MapType.GRIDDED_ONLY);
		if (calcFault && calcGridded)
			mapTypesList.add(MapType.COMBINED);
		if (calcLongTerm) {
			mapTypesList.add(MapType.U3TD);
			mapTypesList.add(MapType.U3TI);
		}
		Preconditions.checkState(!mapTypesList.isEmpty(), "Must calculate something...");
		mapTypes = mapTypesList.toArray(new MapType[0]);
		
		Map<String, DiscretizedFunc> xValsMap = Maps.newHashMap();
		for (MapType type : mapTypes)
			for (Duration duration : durationsForType(type))
				xValsMap.put(getArchiverName(type, duration), xVals);
		
		archiver = new BinaryCurveArchiver(outputDir, getNumTasks(), xValsMap);
		if (rank == 0)
			archiver.initialize();
		
		printEach = getNumTasks() < 50000;
	}
	
	private Duration[] parseDurations(String durStr) {
		String[] split = durStr.trim().split(",");
		Duration[] durations = new Duration[split.length];
		for (int i=0; i<split.length; i++) {
			if (Character.isDigit(split[i].charAt(0))) {
				// it's a day range
				Preconditions.checkState(Character.isDigit(split[i].charAt(split[i].length()-1)));
				Preconditions.checkState(split[i].contains("-"));
				int dashIndex = split[i].indexOf("-");
				String first = split[i].substring(0, dashIndex);
				int startDay = Integer.parseInt(first);
				String last = split[i].substring(dashIndex+1);
				int endDay = Integer.parseInt(last);
				durations[i] = ETAS_HazardMapCalc.durationForDayRange(startDay, endDay);
			} else {
				durations[i] = DurationConstants.valueOf(split[i]).duration;
			}
			if (rank == 0)
				debug("Parsed duration '"+split[i]+"' to: "+durations[i]);
		}
		return durations;
	}
	
	private String getArchiverName(MapType type, Duration duration) {
		return "results_"+imtName+"_"+type.fileName+"_"+duration.fileName;
	}
	
	private void loadFilePositions() throws IOException {
		filePositions = new long[getNumTasks()];
		fileLengths = new int[getNumTasks()];
		if (rank == 0) {
			debug("Loading file positions");
			
			byte[] buf = new byte[4];
			DataInputStream in = new DataInputStream(new ByteArrayInputStream(buf));
			
			raFile.seek(0l);
			raFile.readFully(buf);
			long pos = 4;
			int count = in.readInt();
			Preconditions.checkState(count == getNumTasks(), "Bad count in file! Expected %s, got %s", getNumTasks(), count);
			
			for (int index=0; index<filePositions.length; index++) {
				filePositions[index] = pos;
				raFile.readFully(buf);
				pos += 4;
				in.reset();
				int checkIndex = in.readInt();
				Preconditions.checkState(index == checkIndex, "Bad index. Expected %s, got %s", index, checkIndex);
				// skip 2 double vals for lat/lon
				pos += 16;
				raFile.seek(pos);
				// read in number of ruptures
				raFile.readFully(buf);
				pos += 4;
				in.reset();
				int numRups = in.readInt();
				Preconditions.checkState(numRups >= 0);
				// now skip ahead to next index
				// 20 bytes per rupture: index (int=4), mean (double=8), stdDev (double=8)
				pos = pos + numRups*20l;
				fileLengths[index] = (int)(pos - filePositions[index]);
				raFile.seek(pos);
			}
			
			debug("Distributing file positions");
		}
		
		MPI.COMM_WORLD.Bcast(filePositions, 0, filePositions.length, MPI.LONG, 0);
		if (rank == 0)
			debug("Distributing file lengths");
		MPI.COMM_WORLD.Bcast(fileLengths, 0, fileLengths.length, MPI.INT, 0);
	}

	@Override
	protected int getNumTasks() {
		return region.getNodeCount();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		List<Future<Integer>> futures = Lists.newArrayList();
		
		for (int index : batch)
			futures.add(executor.submit(new CalcRunnable(index), index));
		
		for (Future<Integer> future : futures)
			future.get();
	}
	
	private class CalcRunnable implements Runnable {
		
		private int index;
		
		public CalcRunnable(int index) {
			this.index = index;
		}

		@Override
		public void run() {
			try {
				calculate(index);
				if (printEach) debug("done with "+index);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
	}
	
	private Duration[] durationsForType(MapType type) {
		if (type.isETAS())
			return durations;
		Preconditions.checkState(mapCalc.isCalcLongTerm());
		return mapCalc.getLongTermCalcDurations();
	}
	
	private void calculate(int index) throws IOException {
		Table<Duration, MapType, CurveMetadata> curveMetas = HashBasedTable.create();
		
		boolean allDone = true;
		for (MapType type : mapTypes) {
			for (Duration duration : durationsForType(type)) {
				CurveMetadata meta = new CurveMetadata(sites.get(index), index, null, getArchiverName(type, duration));
				curveMetas.put(duration, type, meta);
				allDone = allDone && archiver.isCurveCalculated(meta, xVals);
			}
		}
		
		if (allDone) {
			debug(index+" already done, skipping");
			return;
		}
		
		Map<Integer, double[]> precomputedFaultVals = null;
		if (raFile != null) {
			// load precomputed fault shakemaps
			long pos = filePositions[index];
			int len = fileLengths[index];
			if (printEach) debug("calculating fault "+index+", pos="+pos+", len="+len);

			DataInputStream in;
			synchronized (raFile) {
				byte[] buf = new byte[len];
				in = new DataInputStream(new ByteArrayInputStream(buf));

				raFile.seek(pos);
				raFile.readFully(buf);
			}

			precomputedFaultVals = mapCalc.loadSiteFromInputStream(in, index);
		}
		
		if (printEach)
			debug("Calculating "+index);
		Table<Duration, MapType, DiscretizedFunc> curves = mapCalc.calculateCurves(sites.get(index), precomputedFaultVals);
		
		for (Cell<Duration, MapType, CurveMetadata> cell : curveMetas.cellSet()) {
			Duration duration = cell.getRowKey();
			MapType type = cell.getColumnKey();
			DiscretizedFunc curve = curves.get(duration, type);
			Preconditions.checkState(curve != null, "Curve not calculated for %s, %s! Size=%s", duration, type, curves.size());
			archiver.archiveCurve(curve, cell.getValue());
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option catalogs = new Option("c", "catalogs", true, "ETAS Catalogs Binary File");
		catalogs.setRequired(true);
		ops.addOption(catalogs);
		
		Option faultFile = new Option("f", "fault-data-file", true,
				"Fault shakemap precalc data file. Must supply this, or --solution-file.");
		faultFile.setRequired(false);
		ops.addOption(faultFile);
		
		Option solFile = new Option("sol", "solution-file", true,
				"FaultSystemSolution file for calculating fault IMs on the fly. Must supply this, or --fault-data-file.");
		solFile.setRequired(false);
		ops.addOption(solFile);
		
		Option spacing = new Option("s", "spacing", true, "Grid spacing in degrees");
		spacing.setRequired(true);
		ops.addOption(spacing);
		
		Option imt = new Option("i", "imt", true, "IMT. One of 'SA','PGA','PGV'");
		imt.setRequired(true);
		ops.addOption(imt);
		
		Option period = new Option("p", "period", true, "Period, required if IMT is SA");
		period.setRequired(false);
		ops.addOption(period);
		
		Option outputDir = new Option("o", "output-dir", true, "Output directory");
		outputDir.setRequired(true);
		ops.addOption(outputDir);
		
		Option gmpeFile = new Option("g", "gmpe", true, "GMPE reference name");
		gmpeFile.setRequired(true);
		ops.addOption(gmpeFile);
		
		Option siteFile = new Option("sites", "sites-file", true, "Sites XML file");
		siteFile.setRequired(true);
		ops.addOption(siteFile);
		
		Option griddedSpacing = new Option("gs", "gridded-spacing", true, "Spacing for gridded ruptures in degrees");
		griddedSpacing.setRequired(true);
		ops.addOption(griddedSpacing);
		
		Option noGridded = new Option("ng", "no-gridded", false, "Flag to disable gridded calculation");
		noGridded.setRequired(false);
		ops.addOption(noGridded);
		
		Option noFault = new Option("nf", "no-fault", false, "Flag to disable fault calculation");
		noFault.setRequired(false);
		ops.addOption(noFault);
		
		Option durations = new Option("d", "durations", true,
				"Durations to calculate (comma separated)");
		durations.setRequired(false);
		ops.addOption(durations);
		
		Option longTerm = new Option("lt", "calc-long-term", false, "Flag enable long term comparisons");
		longTerm.setRequired(false);
		ops.addOption(longTerm);
		
		Option longTermDurations = new Option("ltd", "long-term-durations", true,
				"Durations to calculate for long term (comma separated)");
		longTermDurations.setRequired(false);
		ops.addOption(longTermDurations);
		
		Option elastic = new Option("el", "elastic-rebound", true,
				"Applies elastic rebound for UCERF3-TD for the given ETAS scenario name");
		elastic.setRequired(false);
		ops.addOption(elastic);
		
		return ops;
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		executor.shutdown();
		archiver.close();
		if (raFile != null)
			raFile.close();
	}
	
	public static void main(String args[]) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_ETAS_HazardMapCalc.class);
			
			MPJ_ETAS_HazardMapCalc driver = new MPJ_ETAS_HazardMapCalc(cmd);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
