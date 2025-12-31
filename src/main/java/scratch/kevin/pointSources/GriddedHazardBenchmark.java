package scratch.kevin.pointSources;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.Random;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Stopwatch;

import gov.usgs.earthquake.nshmp.gmm.Gmm;

public class GriddedHazardBenchmark {

	public static void main(String[] args) throws IOException {
		NSHM23_WUS_BranchAveragedERF erf = new NSHM23_WUS_BranchAveragedERF();
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		
		GriddedSeismicitySettings settings = erf.getGriddedSeismicitySettings();
		
//		settings = settings.forSurfaceType(BackgroundRupType.FINITE);
//		settings = settings.forPointSourceMagCutoff(5d);
//		settings = settings.forFiniteRuptureSettings(settings.finiteRuptureSettings.forNumSurfaces(12));
//		settings = settings.forSupersamplingSettings(GridCellSupersamplingSettings.DEFAULT.forApplyToFinite(true));
////		settings = settings.forSupersamplingSettings(null);
		
		settings = settings.forSurfaceType(BackgroundRupType.POINT);
//		settings = settings.forDistanceCorrection(PointSourceDistanceCorrections.NONE.get());
//		settings = settings.forDistanceCorrection(PointSourceDistanceCorrections.MEDIAN_RJB.get());
		settings = settings.forDistanceCorrection(PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_CENTERED.get());
//		settings = settings.forDistanceCorrection(PointSourceDistanceCorrections.SUPERSAMPLING_0p1_FIVE_POINT_RJB_DIST.get());
		settings = settings.forPointSourceMagCutoff(5d);
		settings = settings.forSupersamplingSettings(GridCellSupersamplingSettings.DEFAULT);
//		settings = settings.forSupersamplingSettings(null);
		
		erf.setCacheGridSources(true);
		erf.setGriddedSeismicitySettings(settings);
		
		erf.updateForecast();
		
		int numToCalc = 1000;
		ScalarIMR gmm = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.get();
//		ScalarIMR gmm = new NSHMP_GMM_Wrapper.Single(Gmm.COMBINED_ACTIVE_CRUST_2023);
//		ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
		gmm.setIntensityMeasure(PGA_Param.NAME);
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		
		Region siteReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		GriddedRegion siteGridded = new GriddedRegion(siteReg, 0.05, GriddedRegion.ANCHOR_0_0);
		
		Site site = new Site(siteGridded.getLocation(0));
		site.addParameterList(gmm.getSiteParams());
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmm.getIntensityMeasure());
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		Random r = new Random(siteGridded.getNodeCount());
		
		Stopwatch watch = Stopwatch.createStarted();
		
		for (int n=0; n<numToCalc; n++) {
			site.setLocation(siteGridded.getLocation(r.nextInt(siteGridded.getNodeCount())));
			
			double secsStart = watch.elapsed().toMillis()/1000d;
			calc.getHazardCurve(logXVals, site, gmm, erf);
			double secsEnd = watch.elapsed().toMillis()/1000d;
			
			double myDuration = secsEnd - secsStart;
			double avgDuration = secsEnd/(n+1d);
			System.out.println("Finished curve "+(n+1)+" in "+(float)myDuration+" s;\taverage duration: "+(float)avgDuration+" s");
		}
		System.out.println("DONE with "+numToCalc);
	}

}
