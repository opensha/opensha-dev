package scratch.peter.ucerf3;

import static org.opensha.nshmp2.util.Period.GM0P20;

import java.io.File;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

/**
 * Add comments here
 * 
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class DistCalcTest {

	private static final String S = File.separator;

	DistCalcTest(String solSetPath, String branchID, List<Location> locs,
		Period period, boolean epi) {
		FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solSetPath,
			branchID, IncludeBackgroundOption.EXCLUDE, false, true, 1.0);
		erf.updateForecast();
		EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
		Stopwatch sw = Stopwatch.createUnstarted();
		for (Location loc : locs) {
			System.out.println("Starting calc for " + loc);
			sw.reset().start();
			Site site = new Site(loc);
			HazardCalc hc = HazardCalc.create(wrappedERF, site, period, epi);
			hc.call();
			System.out.println("Compute time: " + sw.stop().elapsed(TimeUnit.SECONDS));
		}
	}

	DistCalcTest(String solSetPath, List<Location> locs,
		Period period, boolean epi) {
		FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solSetPath,
			IncludeBackgroundOption.INCLUDE, false, true, 1.0);
		erf.updateForecast();
		EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
		Stopwatch sw = Stopwatch.createUnstarted();
		for (Location loc : locs) {
			System.out.println("Starting calc for " + loc);
			sw.reset().start();
			Site site = new Site(loc);
			HazardCalc hc = HazardCalc.create(wrappedERF, site, period, epi);
			hc.call();
			System.out.println("Compute time: " + sw.stop().elapsed(TimeUnit.SECONDS));
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Period period = GM0P20;
//		String solSetPath = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/tree/"
//			+ "2013_01_14-UC32-COMPOUND_SOL.zip";
		String solSetPath = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/tree/"
		+ "2013_01_14-UC32-MEAN_BRANCH_AVG_SOL_FM31.zip";

//		String solSetPath = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/coulombTest/"
//				+ "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_VarCoulomb0_mean_sol.zip";
//		String refBranch = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3";
		boolean epi = false;

		List<Location> locs = Lists.newArrayList(
//			NEHRP_TestCity.SAN_FRANCISCO.location(),
//			NEHRP_TestCity.OAKLAND.location(),
			NEHRP_TestCity.LOS_ANGELES.location() //,
//			NEHRP_TestCity.RIVERSIDE.location(),
//			NEHRP_TestCity.SAN_DIEGO.location()
			);

		try {
			new DistCalcTest(solSetPath, locs, period, epi);
//			new DistCalcTest(solSetPath, refBranch, locs, period, epi);
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}
	}

}
