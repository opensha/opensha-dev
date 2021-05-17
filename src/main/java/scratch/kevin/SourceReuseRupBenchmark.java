package scratch.kevin;

import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.Frankel96.Frankel96_AdjustableEqkRupForecast;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Stopwatch;

public class SourceReuseRupBenchmark {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		Frankel96_AdjustableEqkRupForecast erf = new Frankel96_AdjustableEqkRupForecast();
//		erf.setParameter(Frankel96_AdjustableEqkRupForecast.BACK_SEIS_NAME,
//				Frankel96_AdjustableEqkRupForecast.BACK_SEIS_ONLY);
		MeanUCERF2 erf = new MeanUCERF2();
		erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_ONLY);
		erf.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_CROSSHAIR);
		erf.updateForecast();
		
		ScalarIMR imr = AttenRelRef.CB_2008.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGA_Param.NAME);
		
		Location loc = new Location(34, -119);
		Site site = new Site(loc);
		
		for (Parameter<?> param : imr.getSiteParams())
			site.addParameter(param);
		
		int num = 10;
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		ArbitrarilyDiscretizedFunc func = IMT_Info.getUSGS_PGA_Function();
		
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<num; i++) {
			System.out.print("  "+i);
			calc.getHazardCurve(func, site, imr, erf);
		}
		watch.stop();
		System.out.println();
		
		double secs = (double)watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
		System.out.println("Time: "+secs+" secs");
	}

}
