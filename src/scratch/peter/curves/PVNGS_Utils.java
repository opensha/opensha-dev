package scratch.peter.curves;

import static org.opensha.nshmp2.util.Period.*;
import static org.opensha.sha.earthquake.param.IncludeBackgroundOption.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.ClassUtils;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.calc.HazardResultWriter;
import org.opensha.nshmp2.calc.HazardResultWriterLocal;
import org.opensha.nshmp2.calc.HazardResultWriterSites;
import org.opensha.nshmp2.calc.ThreadedHazardCalc;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.UpdatedUCERF2.GridSources;
import scratch.UCERF3.utils.UpdatedUCERF2.MeanUCERF2update;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class PVNGS_Utils {

	private static final String S = File.separator;

	private static String brAvgSolPath = "tmp/invSols/tree/2013_01_14-UC32-MEAN_BRANCH_AVG_SOL_FM31.zip";
	private static String UC32solPath = "tmp/invSols/tree/2013_01_14-UC32-COMPOUND_SOL.zip";
	private static String outPath = "tmp/hazard/PALO_VERDE";
	private static String sitePath = "tmp/curves/sites/palo-verde.txt";
	private static Location PALO_VERDE = new Location(33.40,-112.85);
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// run UC32 branch averaged solution with 400 and 600 km cutoffs
		// for PGA, 5Hz, 1Hz, and 4sec and bg, flt, and all
//		runPVNGS();
		
		// run MeanUC2 for same as above

//		runBranches();
		calcRatios();
		
	}
	
	// run the 20 fault branches for 
	private static String solList = "tmp/invSolSets/UC32-DM-MS-U2.txt";
	
	public static void runBranches() throws IOException {
		
		File brFile = new File(solList);
		List<String> branchList = Files.readLines(brFile, Charsets.US_ASCII);
		Period[] periods = { GM0P00, GM1P00};
		boolean epiUnc = false;
		Map<String, Location> siteMap = UC3_CalcUtils.readSiteFile(sitePath);
		LocationList locs = new LocationList();
		for (Location loc : siteMap.values()) {
			locs.add(loc);
		}

		for (String branch : branchList) {
			FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(UC32solPath, branch,
				IncludeBackgroundOption.EXCLUDE, false, true, 1.0);
			erf.updateForecast();
			EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);

			for (Period period : periods) {
				String outDir = outPath + S + "minitree" + S + erf.getName() ;
				HazardResultWriterSites writer = new HazardResultWriterSites(outDir,
					 siteMap);
				writer.writeHeader(period);
				Site site = new Site(PALO_VERDE);
				HazardCalc hc = HazardCalc.create(wrappedERF, site, period, epiUnc);
//				hc.distanceCutoff = 400.0; // re-privatized in calculator
				HazardResult hr = hc.call();
				writer.write(hr);
				writer.close();
			}
		}
	}

	
	private static void runPVNGS() throws IOException {
		String outDir = "MUC2"; //"UC32brAvg";
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00, GM4P00);
		List<String> srcTypes = Lists.newArrayList("all", "bg", "flt");
		List<Double> distances = Lists.newArrayList(400.0, 800.0);

		for (String srcType : srcTypes) {

			// UCERF3
//			IncludeBackgroundOption bgOption = srcType.equals("bg") ? ONLY : 
//				srcType.equals("flt") ? EXCLUDE : INCLUDE;
//			String fName = srcType + "_curves.csv";
//			UCERF3_FaultSysSol_ERF erf = UC3_CalcUtils.getUC3_ERF(brAvgSolPath,
//				bgOption, false, true, 1.0);
//			erf.updateForecast();
//			EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
		
			// UCERF2
			String bgOption = srcType.equals("bg") ? UCERF2.BACK_SEIS_ONLY : 
				srcType.equals("flt") ? UCERF2.BACK_SEIS_EXCLUDE : 
					UCERF2.BACK_SEIS_INCLUDE;
			String fName = srcType + "_curves.csv";
			MeanUCERF2 erf = new MeanUCERF2();
			erf.setParameter(MeanUCERF2.RUP_OFFSET_PARAM_NAME, 5.0);
			erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, UCERF2.PROB_MODEL_POISSON);
			erf.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
			erf.setParameter(UCERF2.FLOATER_TYPE_PARAM_NAME, UCERF2.FULL_DDW_FLOATER);
			erf.getTimeSpan().setDuration(1.0);
			erf.setParameter(UCERF2.BACK_SEIS_NAME, bgOption);
			erf.updateForecast();
			EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);


			for (Period period : periods) {
				for (Double distance : distances) {
					String out = outDir + distance.intValue();
					File outFile = new File(outPath + S + out + S + period + S + fName);
					HazardResultWriter writer = new HazardResultWriterLocal(outFile, period);
					Site site = new Site(PALO_VERDE);
					HazardCalc hc = HazardCalc.create(wrappedERF, site, period, false);
//					hc.distanceCutoff = distance;
					HazardResult hr = hc.call();
					System.out.println(hr.curve());
					writer.write(hr);
					writer.close();
				}
			}
		}
	}
	
	private static void calcRatios() {
		double[] X_PGA = new double[] {0.0050,0.0070,0.0098,0.0137,0.0192,0.0269,0.0376,0.0527,0.0738,0.103,0.145};
		double[] X_1Hz = new double[] {0.0025,0.00375,0.00563,0.00844,0.0127,0.019,0.0285,0.0427,0.0641,0.0961,0.144,0.216,0.324};
		
		double[] UC2_PGA = new double[] {0.03022074256228455,0.02277368494979446,0.01587147481120066,0.009721634468341478,0.004895448629758334,0.001939595153247317,5.832576691547278E-4,1.1805545232053196E-4,1.3588836830571048E-5,7.508749962013619E-7,5.551585206964458E-9};
		double[] UC2_1Hz = new double[] {0.06100000566108958,0.056547013946159484,0.04813566359053752,0.03596234247790725,0.022460001690390752,0.011485152521903852,0.004641611563084904,0.0014544756737834867,3.3132113471247466E-4,5.260794317056016E-5,5.677790463080918E-6,2.665584075685033E-7,1.9621702060201156E-9};
		
		double[] UC3_PGA = new double[] {0.030200566516563644,0.02224955679031878,0.014681187843002775,0.008518981423068544,0.004231364787153559,0.0017591225610641633,5.981816604790419E-4,1.555834080017771E-4,2.8597849477981996E-5,3.0861801963463746E-6,6.10349547322092E-8};
		double[] UC3_1Hz = new double[] {0.06792774289862909,0.057069463816862794,0.043513322133364216,0.029756562708278048,0.017910500838443318,0.009453894219212687,0.0042431035505638625,0.0015950566265212141,4.8172893583496096E-4,1.1136897654695494E-4,1.8202997724564194E-5,1.5840410175064124E-6,5.675379515797179E-9};

		System.out.println(X_PGA.length + " " + UC2_PGA.length + " " + UC3_PGA.length);
		System.out.println(X_1Hz.length + " " + UC2_1Hz.length + " " + UC3_1Hz.length);
		
		DiscretizedFunc f_PGA_U2 = new ArbitrarilyDiscretizedFunc();
		for (int i=1; i<X_PGA.length; i++) {
			f_PGA_U2.set(X_PGA[i], UC2_PGA[i]);
		}

		DiscretizedFunc f_PGA_U3 = new ArbitrarilyDiscretizedFunc();
		for (int i=1; i<X_PGA.length; i++) {
			f_PGA_U3.set(X_PGA[i], UC3_PGA[i]);
		}

		DiscretizedFunc f_1Hz_U2 = new ArbitrarilyDiscretizedFunc();
		for (int i=1; i<X_1Hz.length; i++) {
			f_1Hz_U2.set(X_1Hz[i], UC2_1Hz[i]);
		}

		DiscretizedFunc f_1Hz_U3 = new ArbitrarilyDiscretizedFunc();
		for (int i=1; i<X_1Hz.length; i++) {
			f_1Hz_U3.set(X_1Hz[i], UC3_1Hz[i]);
		}
		
		// 1e-6 ground motions
		double gmU2pga = ProbOfExceed.get(f_PGA_U2, ProbOfExceed.PE1IN10000);
		double gmU3pga = ProbOfExceed.get(f_PGA_U3, ProbOfExceed.PE1IN10000);
		double gmU21hz = ProbOfExceed.get(f_1Hz_U2, ProbOfExceed.PE1IN10000);
		double gmU31hz = ProbOfExceed.get(f_1Hz_U3, ProbOfExceed.PE1IN10000);

		System.out.println("U2 PGA 1e-6 g: " + gmU2pga);
		System.out.println("U3 PGA 1e-6 g: " + gmU3pga);
		System.out.println("U2 1Hz 1e-6 g: " + gmU21hz);
		System.out.println("U3 1Hz 1e-6 g: " + gmU31hz);

	}

}
