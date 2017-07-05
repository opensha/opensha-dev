package scratch.peter.ucerf3.calc;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Add comments here
 *
 * @author Peter Powers
 */
public class KlamathTest {
	/**
	 * @param args
	 */
	public static void main(String[] args) {


		try {
			String path = "tmp/UC33/src/bravg/FM-DM-MS/UC33brAvg_FM31_ABM_ELLB.zip";
			
			// init erf for branch
			FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(path,
				IncludeBackgroundOption.EXCLUDE, false,
				true, 1.0);
			erf.updateForecast();
			
			List<FaultSectionPrefData> sectData = erf.getSolution().getRupSet().getFaultSectionDataList();
			Map<String, Integer> parentSectMap = Maps.newTreeMap();
			for (FaultSectionPrefData fspd : sectData) {
				String name = fspd.getParentSectionName();
				if (parentSectMap.containsKey(name)) continue;
				parentSectMap.put(name, fspd.getParentSectionId());
			}
			for (String key : parentSectMap.keySet()) {
				System.out.println(parentSectMap.get(key) + " " + key);
			}
			
			FaultSystemSolution fss = erf.getSolution();
			FaultSystemRupSet rupSet = fss.getRupSet();
			double[] rates = fss.getRateForAllRups();
			
			System.out.println("======= Klamath East Total =======");
			SummedMagFreqDist mfdTotal = fss.calcNucleationMFD_forParentSect(719, 5.0, 7.6, 17);
			System.out.println(mfdTotal);
			
			System.out.println("======= Klamath East Alone =======");
			// remove those ruptures that span klamath east
			List<Integer> rupIDs = rupSet.getRupturesForParentSection(720);
			for (int rupID : rupIDs) {
				List<FaultSectionPrefData> data = rupSet
					.getFaultSectionDataForRupture(rupID);
				rates[rupID] = 0.0;
//				String name = data.size() + " SECTIONS BETWEEN " +
//					data.get(0).getName() + " AND " +
//					data.get(data.size() - 1).getName();
//
//				System.out.println(rupID + " " +
//					erf.getSolution().getRateForRup(rupID) + " " + name);
			}
			SummedMagFreqDist mfdAlone = fss.calcNucleationMFD_forParentSect(719, 5.0, 7.6, 17);
			System.out.println(mfdAlone);
			
			ArrayList funcs = Lists.newArrayList();
			funcs.add(mfdTotal);
			funcs.add(mfdAlone);
			
			ArrayList<PlotCurveCharacterstics> plotCharsB = Lists.newArrayList(
				getLineChar(new Color(0, 84, 166)),
				getLineChar(new Color(230, 20, 100)));

			GraphWindow graph = new GraphWindow(funcs,
					"Klamath Graben East", plotCharsB);
				graph.setX_AxisLabel("Magnitude");
				graph.setY_AxisLabel("Incremental Rate");
				graph.setYLog(true);
				graph.setX_AxisRange(5.0, 7.6);
				graph.setY_AxisRange(1e-6, 1e-3);

//			new TestCalc(solSetPath, locs, period, epi);
//			mendoTest(solSetPath);
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}
	}
	
	private static PlotCurveCharacterstics getLineChar(Color c) {
		return new PlotCurveCharacterstics(PlotLineType.SOLID,2f, c);
	}


}
