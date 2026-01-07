package scratch.kevin.pointSources;

import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.calc.PointSourceOptimizedExceedProbCalc;
import org.opensha.sha.calc.PointSourceOptimizedSpectraCalc;
import org.opensha.sha.calc.RuptureExceedProbCalculator;
import org.opensha.sha.calc.RuptureSpectraCalculator;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.PointSource.PoissonPointSource;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

public class CalculationSetupExample {

	public static void main(String[] args) {
		// this shows how to create stubs for testing
		
		Location gridLoc = new Location(0d, 0d);
		
		// should test with and without a distance correction
		PointSourceDistanceCorrection distCorr = null;
//		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING.get();
		
		// should test that the cache isn't stale after changing focal mechanisms
		FocalMech mech = FocalMech.STRIKE_SLIP;
//		FocalMech mech = FocalMech.REVERSE;
		
		// create a PointSource instance
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(5.05, 29, 0.1);
		mfd.setAllButTotMoRate(mfd.getMinX(), mfd.getMaxX(), 1d, 1d); // total rate=1/yr, b=1
		PoissonPointSource pointSource = PointSource.poissonBuilder(gridLoc)
				.surfaceBuilder(PointSourceNshm.SURF_BUILDER_DEFAULT)
				.forMFDAndFocalMech(mfd, mech.mechanism, TectonicRegionType.ACTIVE_SHALLOW)
				.duration(1d) // 1 year
				.build();
		
		// here are some IMRs to test with
		// we will need to pass in multiple IMRs and make sure the cache isn't shared between them
		AttenRelRef[] gmmRefs = {
			AttenRelRef.WRAPPED_ASK_2014,
			AttenRelRef.WRAPPED_BSSA_2014,
		};
		
		// here's how you get an actual instance
		ScalarIMR imr = gmmRefs[0].get();
		
		// set an IMT; will need to pass in multiple IMTs and make sure the cache isn't shared between them
		// PGA
		imr.setIntensityMeasure(PGA_Param.NAME);
		// 1s SA
//		imr.setIntensityMeasure(SA_Param.NAME);
//		SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), 1d);
		
		// IMLs at which to compute exceedance probabilities, in ln units
		EvenlyDiscretizedFunc logXValues = new EvenlyDiscretizedFunc(-5d, 1d, 50);
		
		// create a site x distance away from the grid location
		double azimuth = 0d; // radians
		double siteDistance = 100d; // km
		Site site = new Site(LocationUtils.location(gridLoc, azimuth, siteDistance));
		// need to attach site params from the IMR to the site
		site.addParameterList(imr.getSiteParams());
		// set the site in the IMR
		imr.setSite(site);
		
		// instantiate our optimized rupture exceedance calculator
		RuptureExceedProbCalculator basicExceedCalc = RuptureExceedProbCalculator.BASIC_IMPLEMENTATION;
		RuptureExceedProbCalculator optimizedExceedCalc = new PointSourceOptimizedExceedProbCalc();
		
		for (ProbEqkRupture rup : pointSource) {
			EvenlyDiscretizedFunc exceedProbs = logXValues.deepClone();
			optimizedExceedCalc.getExceedProbabilities(imr, rup, exceedProbs);
		}
		
		// now the same, but for spectra
		imr.setIntensityMeasure(SA_Param.NAME);
		RuptureSpectraCalculator basicSpectraCalc = RuptureSpectraCalculator.BASIC_IMPLEMENTATION;
		RuptureSpectraCalculator optimizedSpectraCalc = new PointSourceOptimizedSpectraCalc();
		
		PeriodParam periodParam = (PeriodParam)imr.getIntensityMeasure().getIndependentParameter(PeriodParam.NAME);
		List<Double> periods = periodParam.getAllowedDoubles();
		
		for (ProbEqkRupture rup : pointSource) {
			// probability spectrum of exceeding fixed ln(IML) = -1
			DiscretizedFunc probAtIML = optimizedSpectraCalc.getSA_ExceedProbSpectrum(imr, rup, -1);
			
			// calculate exceedance probabilities at our IMLs for each period
			// this is used when calculating spectrum of IMLs at a fixed exceedance probability
			EvenlyDiscretizedFunc[] exceedProbs = new EvenlyDiscretizedFunc[periods.size()];
			for (int p=0; p<exceedProbs.length; p++)
				exceedProbs[p] = logXValues.deepClone();
			optimizedSpectraCalc.getMultiPeriodExceedProbabilities(imr, periodParam, periods, rup, exceedProbs);
		}
	}

}
