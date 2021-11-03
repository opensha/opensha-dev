package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;
import scratch.kevin.nshm23.InversionsCLI;

public class WaterLevelReproduce {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/OpenSHA/UCERF3/waterlevel_reproduce");
		
		String prefix = "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_run000";
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File(baseDir, prefix+"_sol.zip"));
		
		double[] noMinRates = MatrixIO.doubleArrayFromFile(new File(baseDir, prefix+"_noMinRates.bin"));
		
		double[] waterLevel = new double[noMinRates.length];
		for (int i=0; i<waterLevel.length; i++)
			waterLevel[i] = sol.getRateForRup(i)-noMinRates[i];

		List<double[]> comps = new ArrayList<>();
		List<SectSlipRates> compSlips = new ArrayList<>();
		List<String> compNames = new ArrayList<>();
		List<double[]> compOrigSlips = new ArrayList<>();
		
		// build it using the old solution
		FaultSystemRupSet rupSet = sol.getRupSet();
		SectSlipRates origTargetSlips = rupSet.requireModule(SectSlipRates.class);
		U3LogicTreeBranch branch = rupSet.requireModule(U3LogicTreeBranch.class);
		InversionModels model = branch.getValue(InversionModels.class);
		FaultModels fm = branch.getValue(FaultModels.class);
		InversionTargetMFDs targetMFDs = rupSet.requireModule(InversionTargetMFDs.class);
		
		UCERF3InversionConfiguration config = UCERF3InversionConfiguration.forModel(model, rupSet, fm, targetMFDs);
		double[] newBasis = config.getMinimumRuptureRateBasis();
		
		double fract = config.getMinimumRuptureRateFraction();
		
		double[] calcWaterLevel = new double[newBasis.length];
		for (int i=0; i<calcWaterLevel.length; i++)
			calcWaterLevel[i] = newBasis[i]*fract;
		
		comps.add(calcWaterLevel);
		compNames.add("New from old solution");
		compSlips.add(rupSet.requireModule(SectSlipRates.class));
		compOrigSlips.add(rupSet.getSlipRateForAllSections());
		
//		// build it using the new pipelines with old rup set
//		rupSet = FaultSystemRupSet.buildFromExisting(rupSet, false).forU3Branch(branch).build();
//		comps.add(InversionsCLI.getU3Generator(rupSet).getWaterLevelRates());
//		compNames.add("New method, old sects");
//		compSlips.add(rupSet.requireModule(SectSlipRates.class));
//		compOrigSlips.add(rupSet.getSlipRateForAllSections());
		
		List<? extends FaultSection> subSects = RuptureSets.getU3SubSects(fm);
		rupSet = FaultSystemRupSet.builder(subSects, rupSet.getSectionIndicesForAllRups()).forU3Branch(branch).build();
		comps.add(InversionsCLI.getU3Generator(rupSet).getWaterLevelRates());
		compNames.add("Fully New Method");
		compSlips.add(rupSet.requireModule(SectSlipRates.class));
		compOrigSlips.add(rupSet.getSlipRateForAllSections());
		
		Range wlRange = new Range(1e-12, 1e-4);
		Range slipRange = new Range(1e-6, 1);
		
		for (int c=0; c<comps.size(); c++) {
			String compName = compNames.get(c);
			double[] newWaterLevel = comps.get(c);
			SectSlipRates newSlips = compSlips.get(c);

			System.out.println("*** "+compName+" ***");
			scatter(waterLevel, newWaterLevel, "Water Level Compare", compName, wlRange);
			System.out.println("Slip Rates: "+newSlips.getClass().getName());
			scatter(sol.getRupSet().getSlipRateForAllSections(), compOrigSlips.get(c), "DM Orig Slip Rate Compare", compName, slipRange);
			System.out.println("Slip Rates: "+newSlips.getClass().getName());
			scatter(origTargetSlips.getSlipRates(), newSlips.getSlipRates(), "Target Slip Rate Compare", compName, slipRange);
		}
	}
	
	private static void scatter(double[] orig, double[] newVals, String compName, String title, Range range) {
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		MinMaxAveTracker diffTrack = new MinMaxAveTracker();
		MinMaxAveTracker ratioTrack = new MinMaxAveTracker();
		
		int zeroDiffs = 0;
		for (int i=0; i<orig.length; i++) {
			scatter.set(orig[i], newVals[i]);
			diffTrack.addValue(newVals[i]-orig[i]);
			if (orig[i] == 0d || newVals[i] == 0d) {
				if (orig[i] != 0d || newVals[i] != 0d)
					zeroDiffs++;
			} else {
				ratioTrack.addValue(newVals[i]/orig[i]);
			}
		}
		
		System.out.println("Orig Range: \t"+StatUtils.min(orig)+", "+StatUtils.max(orig));
		
		System.out.println("New Range:\t"+StatUtils.min(newVals)+", "+StatUtils.max(newVals));
		
		System.out.println("Diffs:\t"+diffTrack);
		System.out.println("Ratios:\t"+ratioTrack);
		System.out.println("# 0-diffs:\t"+zeroDiffs);
		
		GraphWindow gw = new GraphWindow(scatter, title, new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		gw.setX_AxisLabel("Original");
		gw.setY_AxisLabel(compName);
		gw.setAxisRange(range, range);
		gw.setXLog(true);
		gw.setYLog(true);
		gw.setVisible(true);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		System.out.println();
	}

}
