package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.EnumMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class SubductionHazardBenchmark {

	public static void main(String[] args) throws IOException {
		Region reg = PRVI25_RegionLoader.loadPRVI_ModelBroad();
		double spacing = 0.1;
		GriddedRegion gridReg = new GriddedRegion(reg, spacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Will calculate for "+gridReg.getNodeCount()+" sites");
		
		File solFile = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_08_16-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.ONLY;
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, bgOp);
		erf.updateForecast();
		
		Map<TectonicRegionType, ScalarIMR> gmmMap = new EnumMap<>(TectonicRegionType.class);
		gmmMap.put(TectonicRegionType.SUBDUCTION_INTERFACE, AttenRelRef.USGS_PRVI_INTERFACE.get());
		gmmMap.put(TectonicRegionType.SUBDUCTION_SLAB, AttenRelRef.USGS_PRVI_SLAB.get());
		
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logXVals.set(Math.log(xVals.getX(i)), 0d);
		
		int printDelta = 10;
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Site site = new Site(gridReg.getLocation(i));
			site.addParameterList(gmmMap.values().iterator().next().getSiteParams());
			calc.getHazardCurve(logXVals, site, gmmMap, erf);
			System.out.print(".");
			if ((i+1) % printDelta == 0) {
				double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
				double secsEach = secs/(i+1);
				System.out.println(" "+(i+1)+" DONE in "+(float)secs+" s, "+(float)secsEach+" s/calc");
			}
		}
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		double secsEach = secs/(gridReg.getNodeCount());
		System.out.println(" "+(gridReg.getNodeCount())+" DONE in "+(float)secs+" s, "+(float)secsEach+" s/calc");
	}

}
