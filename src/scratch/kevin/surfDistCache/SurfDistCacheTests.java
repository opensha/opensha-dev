package scratch.kevin.surfDistCache;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.FaultRuptureSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.MultiDistanceCache;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class SurfDistCacheTests {
	
	static enum TestType {
		DIRECT_CALC,
		HAZARD_CURVE;
	}

	public static void main(String[] args) throws IOException, DocumentException, InterruptedException {
//		int numThreads=4, numSites=24;
		int numThreads=1, numSites=8;
		TestType type = TestType.HAZARD_CURVE;
		File baSolFile = new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		if (args.length == 4) {
			type = TestType.valueOf(args[0]);
			baSolFile = new File(args[1]);
			numThreads = Integer.parseInt(args[2]);
			numSites = Integer.parseInt(args[3]);
		} else if (args.length != 0) {
			System.out.println("USAGE: <sol file> <num threads> <num sites>");
			System.exit(2);
		}
		
//		System.setProperty(SurfaceCachingPolicy.SIZE_PROP, "8");
		final FaultSystemSolutionERF erf = new FaultSystemSolutionERF(FaultSystemIO.loadSol(baSolFile));
		erf.updateForecast();
		
//		final DiscretizedFunc xVals = IMT_Info.getUSGS_PGA_Function();
		final DiscretizedFunc xVals = new LightFixedXFunc(IMT_Info.getUSGS_PGA_Function());
		
		final TestType myType = type;
		
		System.out.println("Caching policy: "+SurfaceCachingPolicy.getPolicyStr());
		System.out.println("Starting for "+numSites+" sites, "+numThreads+" threads...");
		System.out.println("Test type: "+type);
		
		List<Task> tasks = Lists.newArrayList();
		
		System.gc();
		Thread.sleep(1000);
		
		final Stopwatch watch = Stopwatch.createStarted();
		
		for (int i=0; i<numSites; i++) {
			final int j = i;
			tasks.add(new Task() {
				
				@Override
				public void compute() {
					ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
					imr.setParamDefaults();
					imr.setIntensityMeasure(PGA_Param.NAME);
					
					Location loc = new Location(35+Math.random()*0.01, -118+Math.random()*0.01);
					
					switch (myType) {
					case DIRECT_CALC:
						for (ProbEqkSource source : erf) {
							if (!(source instanceof FaultRuptureSource))
								continue;
							for (ProbEqkRupture rup : source) {
								RuptureSurface surf = rup.getRuptureSurface();
								surf.getDistanceJB(loc);
								surf.getDistanceRup(loc);
								surf.getDistanceX(loc);
								surf.getDistanceSeis(loc);
							}
						}
						break;
					case HAZARD_CURVE:
						HazardCurveCalculator calc = new HazardCurveCalculator();
						DiscretizedFunc func = xVals.deepClone();
						Site site = new Site(loc);
						site.addParameterList(imr.getSiteParams());
						calc.getHazardCurve(func, site, imr, erf);
						break;

					default:
						throw new IllegalStateException("Unknown test type: "+myType);
					}
					System.out.println("Done "+j+" (elapsed: "+watch.elapsed(TimeUnit.SECONDS)+")");
				}
			});
		}
		
		new ThreadedTaskComputer(tasks).computeThreaded(numThreads);
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" s");
		System.gc();
		Thread.sleep(1000);
		System.gc();
		Runtime rt = Runtime.getRuntime();
		long totalMB = rt.totalMemory() / 1024 / 1024;
		long freeMB = rt.freeMemory() / 1024 / 1024;
		long usedMB = totalMB - freeMB;
		
		System.out.println("mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
		
		MultiDistanceCache.printDebugStats();
	}

}
