package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;

import com.google.common.base.Preconditions;

public class SegAdjustMisfitsPlot {

	public static void main(String[] args) throws IllegalStateException, IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
//		File mainDir = new File(invDir, "2022_01_25-coulomb-fm31-ref_branch-dist_cutoff-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate");
////		File mainDir = new File(invDir, "2022_01_25-coulomb-fm31-ref_branch-dist_cutoff-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-NuclMFD");
//		
//		String[] prefixes = {
//			"MaxDist1km",
//			"MaxDist5km",
//			"MaxDist12km"
//		};
//		double[] xVals = {
//			1d,
//			5d,
//			12d,
//		};
//		String xLabel = "Max Jump Dist";

		File mainDir = new File(invDir, "2022_01_28-coulomb-fm31-ref_branch-seg_model-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate");
//		File mainDir = new File(invDir, "2022_01_28-coulomb-fm31-ref_branch-seg_model-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-NuclMFD");
		String[] prefixes = {
			"ShawR0_1",
			"ShawR0_2",
			"ShawR0_3",
			"ShawR0_4",
			"ShawR0_5",
		};
		double[] xVals = {
			1d,
			2d,
			3d,
			4d,
			5d,
		};
		String xLabel = "Shaw & Dieterich (2007) R0";
		
		DiscretizedFunc adjustedFunc = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc unadjustedFunc = new ArbitrarilyDiscretizedFunc();
		
		for (File subDir : mainDir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			File solFile = new File(subDir, "solution.zip");
			
			if (!solFile.exists())
				continue;
			
			InversionMisfitStats stats = FaultSystemSolution.load(solFile).requireModule(InversionMisfitStats.class);
			
			double avgMAD = 0d;
			int numAVG = 0;
			for (MisfitStats misfits : stats.getStats()) {
				if (misfits.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
					numAVG++;
					avgMAD += misfits.absMean;
				}
			}
			avgMAD /= numAVG;
			
			System.out.println(subDir.getName()+": "+avgMAD);
			
			double x = Double.NaN;
			for (int i=0; i<prefixes.length; i++) {
				String prefix = prefixes[i];
				if (subDir.getName().startsWith(prefix)) {
					x = xVals[i];
					break;
				}
			}
			Preconditions.checkState(Double.isFinite(x));
			
			boolean adjusted = subDir.getName().endsWith("adj_targets");
			if (adjusted)
				adjustedFunc.set(x, avgMAD);
			else
				unadjustedFunc.set(x, avgMAD);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		unadjustedFunc.setName("Original Targets");
		funcs.add(unadjustedFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		adjustedFunc.setName("Adjusted Targets");
		funcs.add(adjustedFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Segmentation/MFD Target Adjustment Misfits",
				xLabel, "Average Constraint MAD");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, new Range(StatUtils.min(xVals), StatUtils.max(xVals)), null);
		
		PlotUtils.writePlots(mainDir, "constraint_misfits", gp, 800, 650, true, false, false);
	}

}
