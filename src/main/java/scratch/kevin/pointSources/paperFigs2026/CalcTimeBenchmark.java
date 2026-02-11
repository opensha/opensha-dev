package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.FIGURES_DIR;
import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.ORIG_SOL_FILE;
import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.PROPOSED_DIST_CORR_MODEL;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy.CacheTypes;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.NEHRP_TestCity;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.kevin.latex.LaTeXUtils;
import scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.Models;

public class CalcTimeBenchmark {
	
	// times all in thread-seconds (not wall time)
	
	
	/*
	 * for quicker and dirtier stats
	 */
//	private static final double MIN_SPINUP_TIME = 60d;
//	private static final int MIN_SPINUP_ROUNDS = 1;
//	
//	private static final double MIN_BENCHMARK_TIME = 60d;
//	private static final int MIN_BENCHMARK_ROUNDS = 3;
//	private static final int THREADS = 10;
	/*
	 * for production stats
	 */
	private static final double MIN_SPINUP_TIME = 60d;
	private static final int MIN_SPINUP_ROUNDS = 2;
	
	private static final double MIN_BENCHMARK_TIME = 600d; // 10 minutes
	private static final int MIN_BENCHMARK_ROUNDS = 5;
	private static final int THREADS = 1;

	private static final boolean ADD_ON_FAULT = true;
	private static final boolean CACHE_GRID_SOURCES = true;

	public static void main(String[] args) throws IOException {
		SurfaceCachingPolicy.force(CacheTypes.THREAD_LOCAL);
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		File outputDir, solFile;
		boolean optimize;
		
		boolean loadExisting;
		
		if (args.length == 0) {
			outputDir = FIGURES_DIR;
			solFile = ORIG_SOL_FILE;
			optimize = true;
			loadExisting = false;
		} else {
			Preconditions.checkState(args.length == 3, "Usage: <output-dir> <sol-file> <optimize>");
			outputDir = new File(args[0]);
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
					"Output dir doesn't exist and couldn't be created: %s", outputDir.getAbsolutePath());
			solFile = new File(args[1]);
			Preconditions.checkState(solFile.exists(), "Solution file doesn't exist: %s", solFile.getAbsolutePath());
			optimize = Boolean.parseBoolean(args[2]);
			loadExisting = false;
		}
		
		System.out.println("Loading solution and instantiating ERF");
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(FaultSystemSolution.load(solFile));
		
//		// test to see the memory requirements of 100x
//		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
//		erf.setGriddedSeismicitySettings(Models.FINITE_100X_UNCENTERED.getGridProps());
//		erf.setCacheGridSources(CACHE_GRID_SOURCES);
//		System.out.println("Updating forecast");
//		erf.updateForecast();
//		System.out.println("Done, iterating");
//		int sourceCount = 0;
//		int rupCount = 0;
//		for (ProbEqkSource source : erf) {
//			sourceCount++;
//			for (ProbEqkRupture rup : source)
//				rupCount++;
//		}
//		System.out.println("Iterated over "+groupedDF.format(sourceCount)+" sources and "+groupedDF.format(rupCount)+" ruptures");
//		sourceCount = 0;
//		rupCount = 0;
//		for (ProbEqkSource source : erf) {
//			sourceCount++;
//			for (ProbEqkRupture rup : source)
//				rupCount++;
//		}
//		System.out.println("Iterated over "+groupedDF.format(sourceCount)+" sources and "+groupedDF.format(rupCount)+" ruptures");
//		System.exit(0);
		
		String outputPrefix, texPrefix;
		if (optimize) {
			outputPrefix = "benchmark_stats_optimized";
			texPrefix = "BenchmarkOptimized";
		} else {
			outputPrefix = "benchmark_stats";
			texPrefix = "Benchmark";
		}
		
		System.out.println("***************************************");
		System.out.println("Benchmarking with ADD_ON_FAULT="+ADD_ON_FAULT
				+", CACHE_GRID_SOURCES="+CACHE_GRID_SOURCES+", optimize="+optimize+", THREADS="+THREADS);
		System.out.println("\tSpinup settings:\t"+MIN_SPINUP_ROUNDS+" rounds in "+(float)MIN_SPINUP_TIME+"s");
		System.out.println("\tBenchmark settings:\t"+MIN_BENCHMARK_ROUNDS+" rounds in "+(float)MIN_BENCHMARK_TIME+"s");
		System.out.println("\tOutput dir:\t"+outputDir.getAbsolutePath());
		System.out.println("\tOutput prefix:\t"+outputPrefix);
		System.out.println("***************************************");
		
//		boolean optimize = true;
//		String outputPrefix = "benchmark_stats_optimized";
//		String texPrefix = "BenchmarkOptimized";
		
//		boolean optimize = false;
//		String outputPrefix = "benchmark_stats";
//		String texPrefix = "Benchmark";
		
		List<Location> siteLocs = new ArrayList<>();
		Region siteReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
//		GriddedRegion siteGridded = new GriddedRegion(siteReg, 2d, GriddedRegion.ANCHOR_0_0);
//		System.out.println("Gridded region has "+siteGridded.getNodeCount()+" nodes");
//		for (int i=0; i<siteGridded.getNodeCount(); i++)
//			siteLocs.add(siteGridded.getLocation(i));
		
		for (NEHRP_TestCity site : NEHRP_TestCity.values()) {
			Location loc = site.location();
			if (siteReg.contains(loc))
				siteLocs.add(loc);
		}
		
		File csvFile = new File(outputDir, outputPrefix+".csv");
		if (loadExisting)
			Preconditions.checkState(csvFile.exists(), "loadExisting=true but CSV doesn't exist: %s", csvFile.getAbsolutePath());
		
		System.out.println("Using "+siteLocs.size()+" site locations");
//		System.exit(0);

		ArrayDeque<ScalarIMR> gmms = new ArrayDeque<>();
		ArrayDeque<HazardCurveCalculator> calcs = new ArrayDeque<>();
		for (int i=0; i<THREADS; i++) {
			ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
//			ScalarIMR gmm = new NSHMP_GMM_Wrapper.Single(Gmm.COMBINED_ACTIVE_CRUST_2023);
//			ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
			gmm.setIntensityMeasure(PGA_Param.NAME);
			gmms.push(gmm);
			HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
			calc.setPointSourceOptimizationsEnabled(optimize);
			calcs.push(calc);
		}
		ScalarIMR gmm0 = gmms.peek();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmm0.getIntensityMeasure());
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		List<Site> sites = new ArrayList<>();
		
		for (Location loc : siteLocs) {
			Site site = new Site(loc);
			site.addParameterList(gmm0.getSiteParams());
			sites.add(site);
		}
		
		ExecutorService exec = THREADS == 1 ? null : Executors.newFixedThreadPool(THREADS);
		
		// benchmark off fault
		BenchmarkResult onFaultResult = !loadExisting && ADD_ON_FAULT ? benchmark(erf, gmms, sites, calcs, logXVals, exec, null) : null;
		
		List<Models> models = new ArrayList<>();
		models.add(Models.TRUE_POINT);
		models.add(Models.AS_PUBLISHED);
		models.add(PROPOSED_DIST_CORR_MODEL);
		models.add(Models.FINITE_1X_UNCENTERED);
		models.add(Models.FINITE_2X_UNCENTERED);
		models.add(Models.FINITE_5X_UNCENTERED);
		models.add(Models.FINITE_10X_UNCENTERED);
		models.add(Models.FINITE_20X_UNCENTERED);
		models.add(Models.FINITE_50X_UNCENTERED);
		models.add(Models.FINITE_100X_UNCENTERED);

		Map<Models, List<Models>> comparisons = new HashMap<>();
		for (Models model : models)
			comparisons.put(model, new ArrayList<>());
		
		// compare all to as published
		for (Models model : models)
			if (model != Models.AS_PUBLISHED)
				comparisons.get(model).add(Models.AS_PUBLISHED);
		// compare all to as proposed
		for (Models model : models)
			if (model != PROPOSED_DIST_CORR_MODEL)
				comparisons.get(model).add(PROPOSED_DIST_CORR_MODEL);
		
		Map<Models, BenchmarkResult> results = new HashMap<>();
		
		if (loadExisting) {
			CSVFile<String> csv = CSVFile.readFile(csvFile, true);
			int startRow=1;
			if (ADD_ON_FAULT) {
				Preconditions.checkState(csv.get(startRow, 0).equals("On-Fault"), "Expected on-fault to be first");
				int rounds = csv.getInt(startRow, 1);
				int curves = csv.getInt(startRow, 2);
				double median = csv.getDouble(startRow, 4);
				List<Double> batchTimes = new ArrayList<>(rounds);
				for (int i=0; i<rounds; i++)
					batchTimes.add(median);
				onFaultResult = new BenchmarkResult(curves/rounds, batchTimes);
				startRow++;
			}
			for (int row=startRow; row<csv.getNumRows(); row++) {
				String name = csv.get(row, 0);
				Models model = null;
				for (Models candidate : models) {
					if (candidate.name.equals(name)) {
						model = candidate;
						break;
					}
				}
				Preconditions.checkNotNull(model, "Didn't find a match for %s", name);
				int rounds = csv.getInt(row, 1);
				int curves = csv.getInt(row, 2);
				double median = csv.getDouble(row, 4);
				List<Double> batchTimes = new ArrayList<>(rounds);
				for (int i=0; i<rounds; i++)
					batchTimes.add(median);
				BenchmarkResult result = new BenchmarkResult(curves/rounds, batchTimes);
				if (ADD_ON_FAULT)
					result = addOnFault(result, onFaultResult);
				results.put(model, result);
			}
		} else {
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Name", "Rounds Calculated", "Curves Calculated", "Total Time (s)", "Median Round Time (s)", "Median Rate (/s)");
			if (ADD_ON_FAULT)
				csv.addLine("On-Fault", onFaultResult.numRounds+"", onFaultResult.curvesCalculated+"",
						secDF.format(onFaultResult.totalTimeSec), secDF.format(onFaultResult.medianBatchTime), (float)onFaultResult.curvesPerSec+"");
			for (Models model : models) {
				BenchmarkResult result = benchmark(erf, gmms, sites, calcs, logXVals, exec, model);
				csv.addLine(model.getName(), result.numRounds+"", result.curvesCalculated+"",
						secDF.format(result.totalTimeSec), secDF.format(result.medianBatchTime), (float)result.curvesPerSec+"");
				if (ADD_ON_FAULT)
					result = addOnFault(result, onFaultResult);
				results.put(model, result);
			}
			csv.writeToFile(csvFile);
		}
		
		if (exec != null)
			exec.shutdown();
		
		File texFile = new File(outputDir, outputPrefix+".tex");
		File texTemp = new File(outputDir, texFile.getName()+".tmp");
		FileWriter texFW = new FileWriter(texTemp);

		DecimalFormat oneDF = new DecimalFormat("0.0");
		DecimalFormat intDF = new DecimalFormat("0");
		
		for (Models model : models) {
			System.out.println("**************************************");
			System.out.println(model.getName());
			BenchmarkResult result = results.get(model);
			Preconditions.checkNotNull(result, "No result calculated for %s", model);
			System.out.println(result);
			
			for (Models comp : comparisons.get(model)) {
				System.out.println("Compared against: "+comp.getName());
				BenchmarkResult compResult = results.get(comp);
				Preconditions.checkNotNull(comp, "No result calculated for %s", comp);
				double speedup = result.curvesPerSec/compResult.curvesPerSec;
				double invSpeedup = 1d/speedup;
				double pDiff = 100d*(result.curvesPerSec - compResult.curvesPerSec)/compResult.curvesPerSec;
				double invPDiff = 100d*(compResult.curvesPerSec - result.curvesPerSec)/result.curvesPerSec;
				texFW.write("% "+model.getName()+" vs "+comp.getName()+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+model.texName+"Vs"+comp.texName, oneDF.format(speedup))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+""+comp.texName+"Vs"+model.texName, oneDF.format(invSpeedup))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Percent"+model.texName+"Vs"+comp.texName, getPDiffStr(pDiff))+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"Percent"+comp.texName+"Vs"+model.texName, getPDiffStr(invPDiff))+"\n");
				System.out.println("\tSpeedup:\t"+(float)speedup+"\t("+getPDiffStr(pDiff)+")");
				System.out.println("\tInv. Speedup:\t"+(float)(1d/speedup)+"\t("+getPDiffStr(invPDiff)+")");
				texFW.write("\n");
			}
			
			System.out.println("**************************************\n");
		}
		
		texFW.close();
		Files.move(texTemp, texFile);
	}
	
	private static String getPDiffStr(double pDiff) {
		String ret = LaTeXUtils.groupedIntNumber(pDiff)+"%";
		if (pDiff > 0d)
			ret = "+"+ret;
		return ret;
	}
	
	private static BenchmarkResult benchmark(BaseFaultSystemSolutionERF erf, Deque<ScalarIMR> gmms, List<Site> sites,
			Deque<HazardCurveCalculator> calcs, DiscretizedFunc logXVals, ExecutorService exec, Models model) {
		String modelName;
		Runnable cleanup = null;
		if (model == null) {
			modelName = "Fault-Only";
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			if (erf.getSolution().hasModule(RupMFDsModule.class))
				erf.setParameter(UseRupMFDsParam.NAME, false);
			erf.updateForecast();
		} else {
			modelName = model.getName();
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
			erf.setCacheGridSources(CACHE_GRID_SOURCES);
			erf.setGriddedSeismicitySettings(model.getGridProps());
			Function<GridSourceList, GridSourceList> modFunc = model.getGridModFunction();
			FaultSystemSolution origSol = erf.getSolution();
			if (modFunc != null) {
				FaultSystemSolution modSol = new FaultSystemSolution(origSol.getRupSet(), origSol.getRateForAllRups());
				modSol.setGridSourceProvider(modFunc.apply(origSol.requireModule(GridSourceList.class)));
				erf.setSolution(modSol);
				cleanup = () -> erf.setSolution(origSol);
			}
			erf.updateForecast();
		}
		
		System.out.println("Warming up for "+modelName);
		BenchmarkResult result = benchmark(erf, gmms, sites, calcs, logXVals, exec, MIN_SPINUP_TIME, MIN_SPINUP_ROUNDS);
		System.out.println("\tDONE; "+result);
		
		System.out.println("Running "+modelName);
		result = benchmark(erf, gmms, sites, calcs, logXVals, exec, MIN_BENCHMARK_TIME, MIN_BENCHMARK_ROUNDS);
		System.out.println("\tDONE; "+result);
		
		if (cleanup != null)
			cleanup.run();
		return result;
	}
	
	public static volatile double BLACKHOLE;
	
	private static BenchmarkResult benchmark(BaseFaultSystemSolutionERF erf, Deque<ScalarIMR> gmms, List<Site> sites,
			Deque<HazardCurveCalculator> calcs, DiscretizedFunc logXVals, ExecutorService exec, double minTime, int minRounds) {
		
		double totalTime = 0d;
		List<Double> batchTimes = new ArrayList<>();
		while (batchTimes.size() < minRounds || totalTime < minTime) {
			
			System.gc();
			try {
				Thread.sleep(500l);
				double myTime = 0d;
				Stopwatch watch = Stopwatch.createStarted();
				boolean first = true;
				if (exec == null) {
					for (Site site : sites) {
						myTime += new BenchmarkCalc(erf, gmms, site, calcs, logXVals).call();
						if (first) {
							System.out.print("\t\t");
							first = false;
						}
						System.out.print(".");
					}
				} else {
					List<Future<Double>> futures = new ArrayList<>();
					for (Site site : sites)
						futures.add(exec.submit(new BenchmarkCalc(erf, gmms, site, calcs, logXVals)));
					
					for (Future<Double> future : futures) {
						myTime += future.get();
						if (first) {
							System.out.print("\t\t");
							first = false;
						}
						System.out.print(".");
					}
				}
				double myWallTime = watch.stop().elapsed(TimeUnit.MILLISECONDS)/1000d;
				totalTime += myTime;
				batchTimes.add(myTime);
				System.out.println(" done round "+batchTimes.size()+"/"+minRounds+" in "+secDF.format(myTime)
						+"s, wallClock="+secDF.format(myWallTime)+" s, "+secDF.format(totalTime)+"/"+(float)minTime+"s");
			} catch (Exception e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		return new BenchmarkResult(sites.size(), batchTimes);
	}
	
	private static class BenchmarkCalc implements Callable<Double> {

		private BaseFaultSystemSolutionERF erf;
		private Deque<ScalarIMR> gmms;
		private Site site;
		private Deque<HazardCurveCalculator> calcs;
		private DiscretizedFunc logXVals;

		public BenchmarkCalc(BaseFaultSystemSolutionERF erf, Deque<ScalarIMR> gmms, Site site,
				Deque<HazardCurveCalculator> calcs, DiscretizedFunc logXVals) {
			this.erf = erf;
			this.gmms = gmms;
			this.site = site;
			this.calcs = calcs;
			this.logXVals = logXVals;
		}

		@Override
		public Double call() throws Exception {
			HazardCurveCalculator calc;
			ScalarIMR gmm;
			synchronized (BenchmarkCalc.class) {
				calc = calcs.pop();
				gmm = gmms.pop();
			}
			DiscretizedFunc curve = logXVals.deepClone();
			
			Stopwatch watch = Stopwatch.createStarted();
			calc.getHazardCurve(curve, site, gmm, erf);
			watch.stop();
			double sum = 0d;
			for (int i=0; i<curve.size(); i++)
				sum += curve.getY(i);
			
			synchronized (BenchmarkCalc.class) {
				calcs.push(calc);
				gmms.push(gmm);
				BLACKHOLE += sum;
			}
			return watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		}
		
	}
	
	public static final DecimalFormat secDF = new DecimalFormat("0.0");
	public static final DecimalFormat groupedDF = new DecimalFormat("0");
	static {
		groupedDF.setGroupingSize(3);
		groupedDF.setGroupingUsed(true);
	}
	
	private static BenchmarkResult addOnFault(BenchmarkResult result, BenchmarkResult onFaultResult) {
		List<Double> modBatchTimes = new ArrayList<>();
		for (double batchTime : result.batchTimes)
			modBatchTimes.add(batchTime + onFaultResult.medianBatchTime);
		return new BenchmarkResult(result.curvesPerRound, modBatchTimes);
	}
	
	private static class BenchmarkResult {
		public final int numRounds;
		public final int curvesPerRound;
		public final int curvesCalculated;
		public final double[] batchTimes;
		public final double totalTimeSec;
		public final double medianBatchTime;
		public final double curvesPerSec;
		
		private BenchmarkResult(int curvesPerRound, List<Double> batchTimes) {
			this.numRounds = batchTimes.size();
			this.curvesPerRound = curvesPerRound;
			this.curvesCalculated = numRounds * curvesPerRound;
			this.batchTimes = Doubles.toArray(batchTimes);
			this.totalTimeSec = StatUtils.sum(this.batchTimes);
			this.medianBatchTime = DataUtils.median(this.batchTimes);
			this.curvesPerSec = (double)curvesPerRound/medianBatchTime;
		}
		@Override
		public String toString() {
			return numRounds+" rounds ("+groupedDF.format(curvesCalculated)+" curves) in "
					+secDF.format(totalTimeSec)+"s;\tmedian="+secDF.format(medianBatchTime)+"s, "+(float)curvesPerSec+" curves/s";
		}
	}

}
