package scratch.peter.ucerf3.calc;

import static org.opensha.nshmp2.util.Period.GM0P00;

import java.io.File;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.calc.ERF_ID;
import org.opensha.nshmp2.calc.HazardCalc;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class TestCalc {

	private static final String S = File.separator;

	TestCalc(String solSetPath, List<Location> locs,
		Period period, boolean epi) {
		FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(solSetPath,
			IncludeBackgroundOption.EXCLUDE, false, true, 1.0);
		erf.updateForecast();
		EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
		Stopwatch sw = Stopwatch.createUnstarted();
		for (Location loc : locs) {
			System.out.println("Starting calc for " + loc);
			sw.reset().start();
			Site site = new Site(loc);
			HazardCalc hc = HazardCalc.create(wrappedERF, site, period, epi);
			HazardResult result = hc.call();
			System.out.println(result.curve());
			System.out.println("Compute time: " + sw.stop().elapsed(TimeUnit.SECONDS));
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Period period = GM0P00;
		String solSetPath = "/Users/pmpowers/projects/OpenSHA/tmp/UC33/src/bravg/2013_05_03-ucerf3p3-production-first-five_MEAN_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";

		boolean epi = false;

		List<Location> locs = Lists.newArrayList(
			NEHRP_TestCity.LOS_ANGELES.location()
			);

		try {
			
			String path = "tmp/UC33/src/mean/mean_ucerf3_sol.zip";
			FaultSystemSolutionERF erf = UC3_CalcUtils.getNSHMP_UC3_ERF(path,
				IncludeBackgroundOption.INCLUDE, false,
				true, 1.0);
			erf.updateForecast();
			EpistemicListERF wrappedERF = ERF_ID.wrapInList(erf);
			Site site = new Site(NEHRP_TestCity.LOS_ANGELES.location());
			HazardCalc hc = HazardCalc.create(wrappedERF, site, period, epi);
			HazardResult result = hc.call();
			System.out.println(result.curve());

////			String path = "tmp/UC33/src/vars/2013_05_09-ucerf3p3-branch-wt-test_COMPOUND_SOL.zip";
////			String branch = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarSlipWt0.1_VarSlipWtUnNorm10.0_VarPaleo1.2_VarMFDWt10.0_VarSectNuclMFDWt0.01_VarSmoothPaleoSect1000";
//			
//			String path = "tmp/UC33/src/bravg/FM-DM-MS/UC33brAvg_FM31_ABM_ELLB.zip";
////			String branch="FM3_1_ABM_Shaw09Mod_DsrUni_CharConst_M5Rate7.6_MMaxOff7.2_NoFix_SpatSeisU2";
//			
//			// init erf for branch
//			FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(path,
//				IncludeBackgroundOption.EXCLUDE, false,
//				true, 1.0);
//			erf.updateForecast();
//			
//			List<FaultSectionPrefData> sectData = erf.getSolution().getRupSet().getFaultSectionDataList();
//			Map<String, Integer> parentSectMap = Maps.newTreeMap();
//			for (FaultSectionPrefData fspd : sectData) {
//				String name = fspd.getParentSectionName();
//				if (parentSectMap.containsKey(name)) continue;
//				parentSectMap.put(name, fspd.getParentSectionId());
//			}
//			for (String key : parentSectMap.keySet()) {
//				System.out.println(parentSectMap.get(key) + " " + key);
//			}
//			
//			FaultSystemSolution fss = erf.getSolution();
//			FaultSystemRupSet rupSet = fss.getRupSet();
//			double[] rates = fss.getRateForAllRups();
//			
//			List<Integer> IDs = Lists.newArrayList(719, 720); //, 721);
//			for (int id : IDs) {
//				System.out.println("======= " + id + " =======");
//				List<Integer> rupIDs = rupSet.getRupturesForParentSection(id);
//				for (int rupID : rupIDs) {
//					List<FaultSectionPrefData> data = rupSet.getFaultSectionDataForRupture(rupID);
//					rates[rupID] = 0.0;
//					String name = data.size()+" SECTIONS BETWEEN "+data.get(0).getName()+" AND "+data.get(data.size()-1).getName();
//					
//					System.out.println(rupID + " " + erf.getSolution().getRateForRup(rupID) + " " + name);
//				}
//			}
//			List<FaultSectionPrefData> sectData = erf.getSolution().getRupSet().

			
//			new TestCalc(solSetPath, locs, period, epi);
//			mendoTest(solSetPath);
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}
	}
	
	private static void mendoTest(String path) throws Exception {
		
		File file = new File(path);
		FaultSystemSolution fss = FaultSystemIO.loadSol(file);
		for(int s=0; s<fss.getRupSet().getNumSections(); s++) {
			System.out.println(fss.getRupSet().getFaultSectionData(s).getName());
		}
	}

}
