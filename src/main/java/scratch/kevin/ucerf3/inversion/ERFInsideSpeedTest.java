package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipException;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RupInRegionCache;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.CompoundFSSPlots.RegionalMFDPlot;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;

public class ERFInsideSpeedTest {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws ZipException 
	 */
	public static void main(String[] args) throws ZipException, IOException {
		File dir = new File("/tmp");
		File file = new File(dir, "2012_10_29-logic-tree-fm3_1_x7-fm3_2_x1_COMPOUND_SOL.zip");
		FaultSystemSolutionFetcher fetch = CompoundFaultSystemSolution.fromZipFile(file);
		
		InversionFaultSystemSolution sol = fetch.iterator().next();
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
		erf.updateForecast();
		
//		GenericSmartSurfaceFilter filter = new GenericSmartSurfaceFilter(3, 5, 50d);
//		GenericSmartSurfaceFilter filter = null;
		
		RupInRegionCache cache = new RupInRegionCache() {
			
			private ConcurrentMap<Region, ConcurrentMap<String, Boolean>> map = Maps.newConcurrentMap();

			@Override
			public boolean isRupInRegion(ERF erf, ProbEqkSource src, EqkRupture rup, int srcIndex, int rupIndex,
					Region region) {
				RuptureSurface surf = rup.getRuptureSurface();
				if (surf instanceof CompoundSurface) {
					ConcurrentMap<String, Boolean> regMap = map.get(region);
					if (regMap == null) {
						regMap = Maps.newConcurrentMap();
						map.putIfAbsent(region, regMap);
						// in case another thread put it in first
						regMap = map.get(region);
					}
					String key = srcIndex+"_"+rupIndex;
					Boolean inside = regMap.get(key);
					if (inside == null) {
						inside = false;
						for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface())
							if (region.contains(loc)) {
								inside = true;
								break;
							}
						regMap.putIfAbsent(key, inside);
					}
					return inside;
				}
				for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface())
					if (region.contains(loc))
						return true;
				return false;
			}
			
		};
		
		System.out.println("Calculating!");
		
		for (Region r : RegionalMFDPlot.getDefaultRegions()) {
			Stopwatch watch = Stopwatch.createStarted();
			ERF_Calculator.getParticipationMagFreqDistInRegion(erf, r, 5.05, 40, 0.1, true, cache);
			watch.stop();
			System.out.println("Took "+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d)+" secs for "+r.getName());
		}
		
		for (Region r : RegionalMFDPlot.getDefaultRegions()) {
			Stopwatch watch = Stopwatch.createStarted();
			ERF_Calculator.getParticipationMagFreqDistInRegion(erf, r, 5.05, 40, 0.1, true, cache);
			watch.stop();
			System.out.println("Took "+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d)+" secs for "+r.getName());
		}
	}

}
