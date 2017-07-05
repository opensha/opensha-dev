package scratch.peter.curves;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.StringUtils;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;

import com.google.common.collect.Lists;

/**
 * USed to generate hazard curves for UCERF2 Time Independent List
 * 
 * @author Peter Powers
 * @version $Id:$
 */
class UCERF2_RTGM_Generator {

//	private static final String IN_DIR = "/Volumes/Scratch/UCERF3invSols/refCH";
	private static final String OUT_DIR = "/Volumes/Scratch/rtgm/UCERF3-UCERF2map";
	private static AttenRelRef[] imrRefs = { AttenRelRef.NSHMP_2008 };
	private static Period[] periods = { Period.GM0P20, Period.GM1P00 };
//	private static Period[] periods = { Period.GM0P20};
	private static Collection<NEHRP_TestCity> cities;
//	private static File zip;
	private static List<File> zips;
	
	private List<Future<?>> futures;

	static {
		cities = NEHRP_TestCity.getCA();
//		cities = EnumSet.of(LOS_ANGELES, SAN_FRANCISCO, SACRAMENTO);
//		cities = EnumSet.of(NEHRP_TestCity.LOS_ANGELES);
		
		zips = Lists.newArrayList();

		// UCERF2 mapped
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/ucerf2-noInv/FM3_1_UCERF2_COMPARISON_sol.zip"));
		zips.add(new File("/Volumes/Scratch/UCERF3invSols/ucerf2-noInv/FM2_1_UCERF2_COMPARISON_sol.zip"));

		// UCERF2 as inv
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/ucerf2/FM2_1_UC2ALL_AveU2_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2_sol.zip"));

		// UCERF3 ref
//		zip = new File("/Volumes/Scratch/UCERF3invSols/refCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
		
		// UCERF3 var
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_ShConStrDrp_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_2_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_Shaw09Mod_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_ZENG_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllBsqrtLen_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_RelaxMFD_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff8.0_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_HB08_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_ApplyCC_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.2_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate7.0_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate10.0_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_GEOL_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		zips.add(new File("/Volumes/Scratch/UCERF3invSols/varCH/FM3_1_ABM_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip"));
//		File zipDir = new File(IN_DIR);
//		Collection<File> files = FileUtils.listFiles(
//			zipDir, new String[] {"zip"}, false);
//		erfZips = Lists.newArrayList(files);
	}
	
	public static void main(String[] args) {
		for (File zip : zips) {
			new UCERF2_RTGM_Generator(zip);
		}
	}

	private UCERF2_RTGM_Generator(File zip) {
		try {
			int numProc = Runtime.getRuntime().availableProcessors();
			ExecutorService ex = Executors.newFixedThreadPool(3);
			System.out.println("NumProc: " + numProc);
			futures = Lists.newArrayList();
			for (Period period : periods) {
				for (AttenRelRef imrRef : imrRefs) {
					ScalarIMR imr = newIMR(imrRef, period);
//					UCERF2_RTGM_Processor proc = new UCERF2_RTGM_Processor(imr,
//						newU3_ERF(zip), cities, period, OUT_DIR+File.separator+StringUtils.removeEnd(zip.getName(), ".zip"));
					UCERF2_RTGM_Processor proc = new UCERF2_RTGM_Processor(imr,
						newU2_ERF(zip), cities, period, OUT_DIR+File.separator+StringUtils.removeEnd(zip.getName(), ".zip"));
					futures.add(ex.submit(proc));
				}
			}
			ex.shutdown();
			ex.awaitTermination(48, TimeUnit.HOURS);
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}
	}

	
//	private ERF newU3_ERF(File zip) {
//		UCERF3_FaultSysSol_ERF erf = new UCERF3_FaultSysSol_ERF();
//		erf.getParameter(FaultSystemSolutionPoissonERF.FILE_PARAM_NAME).setValue(zip);
//		erf.getParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME).setValue(true);
//		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
//		ts.setDuration(1);
//		erf.setTimeSpan(ts);
//		erf.updateForecast();
//		return erf;
//	}
	
	private ERF newU2_ERF(File zip) {
//		UCERF2_FaultSysSol_ERF erf = new UCERF2_FaultSysSol_ERF(zip);
//		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
//		ts.setDuration(1);
//		erf.setTimeSpan(ts);
//		erf.updateForecast();
//		return erf;
		return null;
	}
	

	// UCERF2
	//	/Volumes/Scratch/UCERF3invSols/ucerf2map/ALLCAL_UCERF2_rakesfixed.zip
	
	// CH vars
	//	FM3_1_NEOK_ShConStrDrp_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_2_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_Shaw09Mod_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_ZENG_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllBsqrtLen_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_RelaxMFD_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff8.0_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_HB08_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_ApplyCC_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.2_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate7.0_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_NEOK_EllB_DsrUni_CharConst_M5Rate10.0_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_GEOL_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip
	//	FM3_1_ABM_EllB_DsrUni_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_sol.zip

//	static ERF newERFflt() {
//		FaultSystemSolutionPoissonERF erf = new FaultSystemSolutionPoissonERF(
//			"/Users/pmpowers/projects/OpenSHA/tmp/invSols/ucerf2map/ALLCAL_UCERF2_rakesfixed.zip");
//		
//		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
//		ts.setDuration(1);
//		erf.setTimeSpan(ts);
//		
//		erf.updateForecast();
//		
//		return erf;
//	}
//	
//	static ERF newERFbg() {
//		MeanUCERF2 erf = new MeanUCERF2();
//		
//		Parameter bgInclude = erf.getParameter(UCERF2.BACK_SEIS_NAME);
//		bgInclude.setValue(UCERF2.BACK_SEIS_ONLY);
//		Parameter bgSrcParam = erf.getParameter(UCERF2.BACK_SEIS_RUP_NAME);
//		bgSrcParam.setValue(UCERF2.BACK_SEIS_RUP_POINT);
//		Parameter floatParam = erf.getParameter(UCERF2.FLOATER_TYPE_PARAM_NAME);
//		floatParam.setValue(UCERF2.FULL_DDW_FLOATER);
//		Parameter probParam = erf.getParameter(UCERF2.PROB_MODEL_PARAM_NAME);
//		probParam.setValue(UCERF2.PROB_MODEL_POISSON);
//		
//		TimeSpan ts = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
//		ts.setDuration(1);
//		erf.setTimeSpan(ts);
//		
//		erf.updateForecast();
//		
//		return erf;
//		
//	}

	static ScalarIMR newIMR(AttenRelRef imrRef, Period period) {
		ScalarIMR imr = imrRef.instance(null); 
		imr.setParamDefaults();
		if (period == Period.GM0P00) {
			imr.setIntensityMeasure("PGA");
		} else {
			imr.setIntensityMeasure("SA");
			imr.getParameter(PeriodParam.NAME).setValue(period.getValue());
		}
		return imr;
	}


}
