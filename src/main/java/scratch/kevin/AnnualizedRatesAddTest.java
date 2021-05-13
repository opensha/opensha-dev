package scratch.kevin;

import org.apache.commons.lang3.time.StopWatch;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

public class AnnualizedRatesAddTest {

	/**
	 * @param args
	 * @throws RemoteException 
	 */
	public static void main(String[] args) {
		MeanUCERF2 erf = new MeanUCERF2();
		
		ScalarIMR imr = new CB_2008_AttenRel(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(SA_Param.NAME);
		
		Site site = new Site(new Location(34, -118));
		for (Parameter<?> param : imr.getSiteParams())
			site.addParameter(param);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		DiscretizedFunc func = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
		
		erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
		erf.updateForecast();
		
		StopWatch watch = new StopWatch();
		
		watch.start();
		calc.getHazardCurve(func, site, imr, erf);
		watch.stop();
		System.out.println("Include: "+(watch.getTime() / 1000d)+" s");
		DiscretizedFunc include_rates = calc.getAnnualizedRates(func, 30);
		
		erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
		erf.updateForecast();
		
		watch.reset();
		watch.start();
		calc.getHazardCurve(func, site, imr, erf);
		watch.stop();
		System.out.println("EXclude: "+(watch.getTime() / 1000d)+" s");
		DiscretizedFunc exclude_rates = calc.getAnnualizedRates(func, 30);
		
		erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_ONLY);
		erf.updateForecast();

		watch.reset();
		watch.start();
		calc.getHazardCurve(func, site, imr, erf);
		watch.stop();
		System.out.println("Only: "+(watch.getTime() / 1000d)+" s");
		DiscretizedFunc only_rates = calc.getAnnualizedRates(func, 30);
		
		DiscretizedFunc comb_rates = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<exclude_rates.size(); i++) {
			comb_rates.set(exclude_rates.getX(i), exclude_rates.getY(i)+only_rates.getY(i));
		}
		
		System.out.println(include_rates);
		System.out.println(comb_rates);
	}

}
