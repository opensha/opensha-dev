package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class SuperSamplingResolutionTests {

	public static void main(String[] args) throws IOException {
		Region testReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		double testSpacing = 0.33d;
//		double testSpacing = 3d;
		GriddedRegion testGridded = new GriddedRegion(testReg, testSpacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Have "+testGridded.getNodeCount()+" test sites");
		
		LocationList siteLocs = new LocationList(testGridded.getNodeCount());
		// randomly shift so that they're not perfectly co-located
		Random rand = new Random(testGridded.getNodeCount());
		for (Location loc : testGridded.getNodeList())
			siteLocs.add(new Location(loc.lat + testSpacing*(rand.nextDouble() - 0.5d),
					loc.lon + testSpacing*(rand.nextDouble() - 0.5d)));
		
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
		
		File outputDir = new File("/tmp");
		String csvPrefix = "super_sampling_results_"+siteLocs.size()+"sites";
		
		GridSourceList gridSources = sol.requireModule(GridSourceList.class);
		
		ReturnPeriods[] testRPs = ReturnPeriods.values();
		
//		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		AttenRelRef gmmRef = AttenRelRef.WRAPPED_ASK_2014;
		
		double period = 0d;
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		
//		double period = 1d;
//		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
		
		boolean cacheGridSources = true;
		
		int threads = FaultSysTools.defaultNumThreads();
		
		SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS);
		
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>(threads);
		ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>(threads);
		
		ExecutorService exec = Executors.newFixedThreadPool(threads);
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.setCacheGridSources(cacheGridSources);
		erf.updateForecast();
		
		System.out.println("Will calculate "+siteLocs.size()+" curves per test");
		
		System.out.println("Calculating fualt-only curves");
		DiscretizedFunc[] faultCurves = calcCurves(siteLocs, erf, period, gmmRef, gmmDeque, calcDeque, xVals, exec, sourceFilters, null);
		System.out.println();
		
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.updateForecast();
		
		List<double[]> samplingParams = new ArrayList<>();
		
		int refIndex = 0; // whatever is first is the reference for comparing hazard; put your best foot forward

		samplingParams.add(new double[] {0.5d, 300d, 0d, 0d}); // double resolution to 300 km
		samplingParams.add(new double[] {1d, 300d, 0d, 0d}); // full to 300 km
		
		int benchmarkIndex = samplingParams.size(); // no supersampling is time benchmark
		samplingParams.add(null); // no supersampling
		
		double[] fixedDists = { 20, 30, 40, 60 };
		double[] borderDists = { 0d, 40d, 60d, 80d, 100d };
		double[] cornerDists = { 0d, 60d, 100d, 200d, 300d };
		for (double fixedDist : fixedDists) {
			for (double borderDist : borderDists) {
				if (borderDist > 0d && borderDist <= fixedDist)
					continue;
				for (double cornerDist : cornerDists) {
					if (cornerDist > 0d && (cornerDist <= fixedDist || cornerDist <= borderDist))
						continue;
					samplingParams.add(new double[] {1d, fixedDist, borderDist, cornerDist});
				}
			}
		}
//		samplingParams.add(new double[] {1d, 20d, 0d, 0d}); // full to 20 km
//		samplingParams.add(new double[] {1d, 30d, 0d, 0d}); // full to 30 km
//		samplingParams.add(new double[] {1d, 40d, 0d, 0d}); // full to 40 km
//		samplingParams.add(new double[] {1d, 40d, 0d, 120d}); // variable
//		samplingParams.add(new double[] {1d, 30d, 60d, 120d}); // variable
//		samplingParams.add(new double[] {0.5d, 30d, 60d, 120d}); // variable
//		samplingParams.add(new double[] {1d, 30d, 100d, 200d}); // variable
//		samplingParams.add(new double[] {1d, 40d, 100d, 200d}); // variable
//		samplingParams.add(new double[] {1d, 20d, 40d, 120d}); // variable
//		samplingParams.add(new double[] {1d, 20d, 40d, 0d}); // variable
//		samplingParams.add(new double[] {1d, 20d, 60d, 0d}); // variable
//		samplingParams.add(new double[] {1d, 30d, 60d, 0d}); // variable
//		samplingParams.add(new double[] {1d, 40d, 80d, 0d}); // variable
		
		List<DiscretizedFunc[]> samplingCurves = new ArrayList<>(samplingParams.size());
		List<Double> times = new ArrayList<>(samplingParams.size());
		List<Integer> indexes = new ArrayList<>(samplingCurves.size());
		
		DiscretizedFunc[] refCurves = null;
		
		for (int s=0; s<samplingParams.size(); s++) {
			double[] params = samplingParams.get(s);
			
			System.out.println(s+"/"+samplingParams.size()+". "+samplingName(params));
			
			DiscretizedFunc[] curves = null;
			double secs = 0;
			
			File cacheFile = null;
			if (s == refIndex && params != null) {
				cacheFile = new File(outputDir, "cache_"+oDF.format(params[0])+"_"+oDF.format(params[1])
						+"_"+oDF.format(params[2])+"_"+oDF.format(params[3])+"_s"+siteLocs.size()+"_p"+oDF.format(period)+".csv");
				if (cacheFile.exists()) {
					System.out.println("Reading cache from: "+cacheFile.getAbsolutePath());
					curves = SolHazardMapCalc.loadCurvesCSV(CSVFile.readFile(cacheFile, true), null);
					Preconditions.checkState(curves.length == siteLocs.size());
					secs = Double.NaN;
				} else {
					System.out.println("Will cache to: "+cacheFile.getAbsolutePath());
				}
			}
			
			if (curves == null) {
				if (params == null)
					gridSources.setSupersamplingParams(0d, 0d, 0d, 0d);
				else
					gridSources.setSupersamplingParams(params[0], params[1], params[2], params[3]);
				erf.clearCachedGridSources();
				Stopwatch watch = Stopwatch.createStarted();
				curves = calcCurves(siteLocs, erf, period, gmmRef, gmmDeque, calcDeque, xVals, exec, sourceFilters, faultCurves);
				watch.stop();
				secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
				
				if (cacheFile != null) {
					System.out.println("Caching to "+cacheFile.getAbsolutePath());
					SolHazardMapCalc.writeCurvesCSV(cacheFile, curves, siteLocs);
				}
			}
			
			indexes.add(s);
			times.add(secs);
			System.out.println("Took "+(float)secs+" s");
			System.out.println();
			if (s == 0)
				refCurves = curves;
			
			samplingCurves.add(curves);
		}
		
		exec.shutdown();
		
		// sort by time, decreasing
		indexes = ComparablePairing.getSortedData(times, indexes);
		Collections.reverse(indexes);
		
		List<CSVFile<String>> csvs = new ArrayList<>(testRPs.length+1);
		for (int r=0; r<testRPs.length+1; r++) {
			CSVFile<String> csv = new CSVFile<>(true);
			
			csv.addLine("Full Sample Distance (km)", "Perimeter Sample Distance (km)", "Corner Distance (km)",
					"Time (s)", "Time Factor", "Minimum", "Maximum", "Average", "MAD", "SRSS", "|Min|+|Max|+MAD");
			
			// placeholder
			List<String> empty = new ArrayList<>(csv.getNumCols());
			for (int i=0; i<csv.getNumCols(); i++)
				empty.add(csvPrefix);
			for (int i=0; i<indexes.size()-1; i++) // -1 here because we skip the first (reference)
				csv.addLine(empty);
			
			csvs.add(csv);
		}
		
		double benchmarkSecs = times.get(benchmarkIndex);
		for (int s : indexes) {
			double[] params = samplingParams.get(s);
			String label = samplingName(params);
			double secs = times.get(s);
			double timeFactor = secs/benchmarkSecs;
			if (s == benchmarkIndex)
				label += "\t"+(float)secs+" s (benchmark)";
			else
				label += "\t"+(float)secs+" s ("+twoDF.format(timeFactor)+"x slower)";
			System.out.println(label);
			if (s == refIndex) {
				System.out.println("\tREFERENCE");
			} else {
				DiscretizedFunc[] curves = samplingCurves.get(s);
				MinMaxAveTracker pDiffTrackAll = new MinMaxAveTracker();
				MinMaxAveTracker absPDiffTrackAll = new MinMaxAveTracker();
				double sumSqAll = 0d;
				
				for (int r=0; r<csvs.size(); r++) {
					MinMaxAveTracker pDiffTrack;
					MinMaxAveTracker absPDiffTrack;;
					double sumSq;
					String rpLabel;
					
					if (r < testRPs.length) {
						rpLabel = testRPs[r].label;
						double[] map = calcMap(curves, testRPs[r]);
						double[] refMap = calcMap(refCurves, testRPs[r]);

						pDiffTrack = new MinMaxAveTracker();
						absPDiffTrack = new MinMaxAveTracker();
						sumSq = 0d;
						for (int i=0; i<map.length; i++) {
							if (refMap[i] == 0d)
								continue;
							double pDiff = 100d*(map[i] - refMap[i])/refMap[i];
							sumSq += pDiff*pDiff;
							sumSqAll += pDiff*pDiff;
							pDiffTrack.addValue(pDiff);
							pDiffTrackAll.addValue(pDiff);
							absPDiffTrack.addValue(Math.abs(pDiff));
							absPDiffTrackAll.addValue(Math.abs(pDiff));
						}
					} else {
						rpLabel = "All RPs";
						pDiffTrack = pDiffTrackAll;
						absPDiffTrack = absPDiffTrackAll;
						sumSq = sumSqAll;
					}
					double minMaxAvgSum = Math.abs(pDiffTrack.getMin()) + Math.abs(pDiffTrack.getMax()) + absPDiffTrack.getAverage();
					System.out.println("\t"+rpLabel+":\t["+pFormat(pDiffTrack.getMin())
							+", "+pFormat(pDiffTrack.getMax())+"], avg="+pFormat(pDiffTrack.getAverage())
							+", mad="+pFormat(absPDiffTrack.getAverage(), false)
							+"\t\tsrss="+twoDF.format(Math.sqrt(sumSq))
							+", |min|+|max|+mad="+twoDF.format(minMaxAvgSum));
					
					CSVFile<String> csv = csvs.get(r);
					
					List<String> line = List.of(
							params == null ? "0.0" : params[1]+"",
							params == null ? "0.0" : params[2]+"",
							params == null ? "0.0" : params[3]+"",
							secs+"",
							timeFactor+"",
							pDiffTrack.getMin()+"",
							pDiffTrack.getMax()+"",
							pDiffTrack.getAverage()+"",
							absPDiffTrack.getAverage()+"",
							Math.sqrt(sumSq)+"",
							minMaxAvgSum+"");
					csv.setLine(s, line); // don't add 1, we're skipping the first index (reference) and the header takes its place
				}
			}
		}
		
		for (int r=0; r<csvs.size(); r++) {
			CSVFile<String> csv = csvs.get(r);
			String rpPrefix = r < testRPs.length ? testRPs[r].name() : "COMBINED";
			File outputFile = new File(outputDir, csvPrefix+"_"+rpPrefix+".csv");
			csv.writeToFile(outputFile);
		}
	}
	
	private static final DecimalFormat twoDF = new DecimalFormat("0.00");
	
	private static String pFormat(double percent) {
		return pFormat(percent, false);
	}
	
	private static String pFormat(double percent, boolean plus) {
		String ret = twoDF.format(percent)+"%";
		if (plus && !ret.startsWith("-"))
			ret = "+"+ret;
		return ret;
	}
	
	private static double[] calcMap(DiscretizedFunc[] curves, ReturnPeriods rp) {
		double[] ret = new double[curves.length];
		
		double level = rp.oneYearProb;
		
		for (int i=0; i<curves.length; i++) {
			// curveLevel is a probability, return the IML at that probability
			double val;
			if (level > curves[i].getMaxY())
				val = 0d;
			else if (level < curves[i].getMinY())
				// saturated
				val = curves[i].getMaxX();
			else
				val = curves[i].getFirstInterpolatedX_inLogXLogYDomain(level);
			ret[i] = val;
		}
		
		return ret;
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static String samplingName(double[] params) {
		if (params == null)
			return "No supersampling";
		if (params[2] == 0d && params[3] == 0d)
			return "Full up to "+oDF.format(params[1])+" km ("+oDF.format(params[0])+" km spacing)";
		return "Full up to "+oDF.format(params[1])+"-"+oDF.format(params[2])+"-"+oDF.format(params[3])+" km ("+oDF.format(params[0])+" km spacing)";
	}
	
	private static DiscretizedFunc[] calcCurves(LocationList siteLocs, FaultSystemSolutionERF erf, double period,
			AttenRelRef gmmRef, ArrayDeque<ScalarIMR> gmmDeque, ArrayDeque<HazardCurveCalculator> calcDeque,
			DiscretizedFunc xVals, ExecutorService exec, SourceFilterManager sourceFilters,
			DiscretizedFunc[] combineWith) {
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logXVals.set(Math.log(xVals.getX(i)), 1d);
		
		DiscretizedFunc[] ret = new DiscretizedFunc[siteLocs.size()];
		
		List<Future<DiscretizedFunc>> futures = new ArrayList<>(siteLocs.size());
		
//		System.out.println("Calculating "+siteLocs.size()+" hazard curves");

		int printMod;
		if (siteLocs.size() > 500)
			printMod = 100;
		else if (siteLocs.size() > 250)
			printMod = 50;
		else
			printMod = 20;
		
		for (Location siteLoc : siteLocs) {
			futures.add(exec.submit(new Callable<DiscretizedFunc>() {

				@Override
				public DiscretizedFunc call() throws Exception {
					HazardCurveCalculator calc = null;
					ScalarIMR gmm = null;
					synchronized (SuperSamplingResolutionTests.class) {
						if (!calcDeque.isEmpty())
							calc = calcDeque.pop();
						if (!gmmDeque.isEmpty())
							gmm = gmmDeque.pop();
					}
					if (calc == null) {
						calc = new HazardCurveCalculator(sourceFilters);
						calc.setTrackProgress(false);
					}
					if (gmm == null)
						gmm = gmmRef.get();
					else
						gmm.setParamDefaults();
					
					if (period == 0d) {
						gmm.setIntensityMeasure(PGA_Param.NAME);
					} else {
						Preconditions.checkState(period > 0d);
						gmm.setIntensityMeasure(SA_Param.NAME);
						SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
					}
					
					DiscretizedFunc logCurve = new LightFixedXFunc(logXVals);
					
					Site site = new Site(siteLoc);
					site.addParameterList(gmm.getSiteParams());
					
					calc.getHazardCurve(logCurve, site, gmm, erf);
					
					synchronized (SuperSamplingResolutionTests.class) {
						calcDeque.push(calc);
						gmmDeque.push(gmm);
					}
					
					DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
					for (int i=0; i<logCurve.size(); i++)
						ret.set(xVals.getX(i), logCurve.getY(i));
					
					return ret;
				}
			}));
		}
		
		for (int i=0; i<ret.length; i++) {
			DiscretizedFunc curve;
			try {
				curve = futures.get(i).get();
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (combineWith != null) {
				Preconditions.checkState(combineWith[i].size() == curve.size());
				DiscretizedFunc combined = new ArbitrarilyDiscretizedFunc();
				for (int j=0; j<curve.size(); j++)
					combined.set(curve.getX(j), 1d - (1d - curve.getY(j)*(1d - combineWith[i].getY(j))));
				curve = combined;
			}
			ret[i] = curve;
			System.out.print(".");
			if ((i+1) % printMod == 0 && i < (ret.length-1))
				System.out.println(" ("+(i+1)+")");
		}
		System.out.println(" DONE");
		
		return ret;
	}

}
