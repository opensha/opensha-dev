package scratch.peter.curves;

import static org.opensha.nshmp2.util.Period.*;

import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.exceptions.ConstraintException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.AbstractEpistemicListERF;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2_TimeIndependentEpistemicList;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import scratch.peter.ucerf3.calc.UC3_CalcUtils;

/**
 * Used to generate hazard curves for epistemic list ERFs.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
class RTGM_ListGenerator {

	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/UCERF2-TimeIndepSRP";
//	private static AttenRelRef[] imrRefs = { AttenRelRef.NSHMP_2008 };
//	private static Period[] periods = { Period.GM0P20, Period.GM1P00 };
	private static Period[] periods = { Period.GM0P00, Period.GM0P20, Period.GM1P00};
	private static Map<String, Location> locMap;
	private static String sitePath;
	
	static {
		// sitePath = "tmp/UC3sites/NEHRPsites.txt";
		// sitePath = "tmp/UC3sites/PBRsites.txt";
		sitePath = "tmp/UC3sites/SRPsites.txt";
		try {
			locMap = UC3_CalcUtils.readSiteFile(sitePath);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		new RTGM_ListGenerator();
//		System.out.println(cities);
//		RTGM_ListProcessor proc = new RTGM_ListProcessor(
//			newIMR(Period.GM1P00), newERF(), NEHRP_TestCity.LOS_ANGELES, Period.GM1P00, OUT_DIR);
//		proc.run();
//		
	}

	private RTGM_ListGenerator() {
		try {
			int threadCt = Runtime.getRuntime().availableProcessors();
			ExecutorService ex = Executors.newFixedThreadPool(threadCt);
			for (Period period : periods) {
				for (String locName : locMap.keySet()) {
					ScalarIMR imr = newIMR(period);
					EpistemicListERF erfs = newERF();
					RTGM_ListProcessor proc = new RTGM_ListProcessor(imr, erfs,
						locName, locMap.get(locName), period, OUT_DIR);
					ex.submit(proc);
//					proc.run();
				}
			}
			ex.shutdown();
			ex.awaitTermination(48, TimeUnit.HOURS);
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}
	}

	static EpistemicListERF newERF() {
		AbstractEpistemicListERF erf = new UCERF2_TimeIndependentEpistemicList();
//		AbstractEpistemicListERF erf = new UCERF2_TimeDependentEpistemicList();
		
		Parameter bgSrcParam = erf.getParameter(UCERF2.BACK_SEIS_RUP_NAME);
		bgSrcParam.setValue(UCERF2.BACK_SEIS_RUP_POINT);
		Parameter floatParam = erf.getParameter(UCERF2.FLOATER_TYPE_PARAM_NAME);
		floatParam.setValue(UCERF2.FULL_DDW_FLOATER);

		// prob model is set to poisson by default
		
		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
		ts.setDuration(1);
		erf.setTimeSpan(ts);
		return erf;
	}

	static ScalarIMR newIMR(Period period) {
		
//		ScalarIMR imr = SourceIMR.WUS_FAULT.instance(period); //new NSHMP08_WUS();
//		imr.getParameter(NSHMP08_WUS.IMR_UNCERT_PARAM_NAME).setValue(false);

		ScalarIMR imr = new NSHMP08_WUS();
		imr.setIntensityMeasure((period == GM0P00) ? PGA_Param.NAME : SA_Param.NAME);
		try {
			imr.getParameter(PeriodParam.NAME).setValue(period.getValue());
		}  catch (ConstraintException ce) { /* do nothing */ }

		imr.getParameter(NSHMP08_WUS.IMR_UNCERT_PARAM_NAME).setValue(false);
		return imr;

//		ScalarIMR imr = AttenRelRef.NSHMP_2008.instance(null);
//		ScalarIMR imr = AttenRelRef.BA_2008.instance(null);
		
//		imr.setParamDefaults();
//		if (period == Period.GM0P00) {
//			imr.setIntensityMeasure("PGA");
//		} else {
//			imr.setIntensityMeasure("SA");
//			imr.getParameter(PeriodParam.NAME).setValue(period.getValue());
//		}
		
//		return imr;
	}


}
