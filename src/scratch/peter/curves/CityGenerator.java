package scratch.peter.curves;

import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.nshmp.NEHRP_TestCity.*;

import java.util.Collection;
import java.util.EnumSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.exceptions.ConstraintException;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.UCERF3.utils.ModUCERF2.ModMeanUCERF2;
import scratch.UCERF3.utils.UpdatedUCERF2.GridSources;
import scratch.UCERF3.utils.UpdatedUCERF2.MeanUCERF2update;
import scratch.UCERF3.utils.UpdatedUCERF2.MeanUCERF2update_FM2p1;
import scratch.UCERF3.utils.UpdatedUCERF2.ModMeanUCERF2update_FM2p1;
import scratch.UCERF3.utils.UpdatedUCERF2.UCERF2_FM2pt1_FSS_ERFupdate;

/**
 * Utility class to generate hazard curves for NEHRP test cities.
 */
class CityGenerator {

	private static final String OUT_DIR = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/ModMeanUC2";
//	private static final String OUT_DIR = "?/Volumes/Scratch/rtgm/MeanUCERF2";
//	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/MeanUCERF2update";
//	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/MeanUCERF2update_FM2P1";
//	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/ModMeanUCERF2update_FM2P1";
//	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/FSS_UC2map";
	private static Period[] periods = { GM0P00, GM0P20, GM1P00};
//	private static Period[] periods = { GM0P00, GM0P20 };
//	private static Period[] periods = { GM0P00 };
	private static Collection<NEHRP_TestCity> cities;

	static {
//		cities = EnumSet.of(LOS_ANGELES, VENTURA); 
		cities = NEHRP_TestCity.getCA();
//		cities = EnumSet.of(NEHRP_TestCity.VENTURA);
	}
	
	public static void main(String[] args) {
		new CityGenerator();
	}

	private CityGenerator() {
		try {
			int threadCt = Math.min(Runtime.getRuntime().availableProcessors(),
				periods.length);
			ExecutorService ex = Executors.newFixedThreadPool(threadCt);
			ERF erf = newERF(); // assumes ERF is threadsafe
			for (Period period : periods) {
				ScalarIMR imr = newIMR(period);
//				ERF erf = newERF();
				CityProcessor proc = new CityProcessor(imr, erf, cities, period, OUT_DIR);
//				proc.run();
				ex.submit(proc);
			}
			ex.shutdown();
			ex.awaitTermination(48, TimeUnit.HOURS);
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}
	}

	static ERF newERF() {
//		MeanUCERF2 erf = new MeanUCERF2();
		ModMeanUCERF2 erf = new ModMeanUCERF2();
//		MeanUCERF2 erf = new MeanUCERF2update(GridSources.ALL);
//		MeanUCERF2 erf = new MeanUCERF2update_FM2p1();
//		ModMeanUCERF2 erf = new ModMeanUCERF2update_FM2p1();
		setParams(erf); // UC2 erfs
		
//		FaultSystemSolutionPoissonERF erf = new UCERF2_FM2pt1_FSS_ERFupdate();
//		erf.getTimeSpan().setDuration(1.0);
		
//		erf.setParameter(MeanUCERF2.RUP_OFFSET_PARAM_NAME, 1.0);
//		erf.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
//		erf.setParameter(UCERF2.FLOATER_TYPE_PARAM_NAME, UCERF2.FULL_DDW_FLOATER);
//		erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, UCERF2.PROB_MODEL_POISSON);
		
//		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
//		ts.setDuration(1);
//		erf.setTimeSpan(ts);
		
		erf.updateForecast();
		
		return erf;
	}

	static ScalarIMR newIMR(Period period) {
//		ScalarIMR imr = SourceIMR.WUS_FAULT.instance(period); //new NSHMP08_WUS();
		
		ScalarIMR imr = new NSHMP08_WUS();
		imr.setIntensityMeasure((period == GM0P00) ? PGA_Param.NAME : SA_Param.NAME);
		try {
			imr.getParameter(PeriodParam.NAME).setValue(period.getValue());
		}  catch (ConstraintException ce) { /* do nothing */ }

		imr.getParameter(NSHMP08_WUS.IMR_UNCERT_PARAM_NAME).setValue(false);
		return imr;
	}

	private static void setParams(ERF uc2) {
		uc2.setParameter(MeanUCERF2.RUP_OFFSET_PARAM_NAME, 1.0);
		uc2.setParameter(UCERF2.PROB_MODEL_PARAM_NAME,
			UCERF2.PROB_MODEL_POISSON);
		uc2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
		uc2.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
		uc2.setParameter(UCERF2.FLOATER_TYPE_PARAM_NAME,
			UCERF2.FULL_DDW_FLOATER);
		uc2.getTimeSpan().setDuration(1.0);
	}

//	static ScalarIMR newIMR(AttenRelRef imrRef, Period period) {
//		ScalarIMR imr = imrRef.instance(null); 
//		imr.setParamDefaults();
//		if (period == Period.GM0P00) {
//			imr.setIntensityMeasure("PGA");
//		} else {
//			imr.setIntensityMeasure("SA");
//			imr.getParameter(PeriodParam.NAME).setValue(period.getValue());
//		}
//		return imr;
//	}


}
