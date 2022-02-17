package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;

import com.google.common.base.Preconditions;

public class LogicTreeMisfitsCompare {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File primaryDir = new File(invDir, "2022_02_16-nshm23_u3_hybrid_branches-no_reweight_use_prev-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip");
//		File primaryDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-CappedRdst-2000ip");
		
		File refDir = new File(invDir, "2022_02_15-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-JumpProb-2000ip");
		
		Quantity quantity = Quantity.MAD;
		
		Map<LogicTreeBranch<?>, InversionMisfitStats> primaryBranchMisfits =
				LogicTreeMisfitPageGen.loadBranchMisfits(new File(primaryDir, "results.zip"));
		Map<LogicTreeBranch<?>, InversionMisfitStats> refBranchMisfits =
				LogicTreeMisfitPageGen.loadBranchMisfits(new File(refDir, "results.zip"));
		
		Preconditions.checkState(primaryBranchMisfits.size() == refBranchMisfits.size());
		
		DefaultXY_DataSet avgScatter = new DefaultXY_DataSet();
		Map<String, DefaultXY_DataSet> constraintScatters = new HashMap<>();
		
		int maxMinAway = 0;
		
		for (LogicTreeBranch<?> branch : primaryBranchMisfits.keySet()) {
			InversionMisfitStats stats1 = primaryBranchMisfits.get(branch);
			
			LogicTreeBranch<?> refBranch;
			if (refBranchMisfits.containsKey(branch)) {
				refBranch = branch;
			} else {
				// find most similar branch, with might be slightly different if one run uses a different branch choice
				int minAway = Integer.MAX_VALUE;
				int numAtMin = 0;
				LogicTreeBranch<?> closest = null;
				for (LogicTreeBranch<?> oBranch : refBranchMisfits.keySet()) {
					int numAway = oBranch.getNumAwayFrom(branch);
					if (numAway < minAway) {
						minAway = numAway;
						numAtMin = 1;
						closest = oBranch;
					} else if (numAway == minAway) {
						numAtMin++;
					}
				}
				Preconditions.checkState(numAtMin == 1, "Multiple branches are equally close to target branch");
				maxMinAway = Integer.max(maxMinAway, minAway);
				refBranch = closest;
			}
			InversionMisfitStats stats2 = refBranchMisfits.get(refBranch);
			
			Map<String, Double> vals1 = loadUncertWeighted(stats1, quantity);
			Map<String, Double> vals2 = loadUncertWeighted(stats2, quantity);
			
			avgScatter.set(avg(vals2), avg(vals1));
			
			Preconditions.checkState(vals1.size() == vals2.size());
			for (String constrName : vals1.keySet()) {
				DefaultXY_DataSet scatter = constraintScatters.get(constrName);
				if (scatter == null) {
					scatter = new DefaultXY_DataSet();
					constraintScatters.put(constrName, scatter);
				}
				scatter.set(vals2.get(constrName), vals1.get(constrName));
			}
		}
		if (maxMinAway > 0)
			System.out.println("Had to find closest branches, furthest away was "+maxMinAway);
		
		String prefix = "misfits_compare";
		plotScatter(primaryDir, prefix, quantity, "All Constraint Average", avgScatter);
		
		for (String constrName : constraintScatters.keySet()) {
			String myPrefix = prefix+"_"+constrName.replaceAll("\\W+", "_");
			plotScatter(primaryDir, myPrefix, quantity, constrName, constraintScatters.get(constrName));
		}
	}
	
	private static Map<String, Double> loadUncertWeighted(InversionMisfitStats stats, Quantity quantity) {
		Map<String, Double> ret = new HashMap<>();
		for (MisfitStats misfits : stats.getStats())
			if (misfits.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
				ret.put(misfits.range.shortName, misfits.get(quantity));
		return ret;
	}
	
	private static double avg(Map<String, Double> stats) {
		double avg = 0d;
		for (Double val : stats.values())
			avg += val;
		avg /= stats.size();
		return avg;
	}
	
	private static void plotScatter(File outputDir, String prefix, Quantity quantity, String title, XY_DataSet scatter)
			throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		double max = Math.max(scatter.getMaxX(), scatter.getMaxY());
		max = 0.5*Math.ceil(2d*max);
		
		Range range = new Range(0d, max);
		
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.set(0d, 0d);
		line.set(max, max);
		
		funcs.add(line);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Reference "+quantity, "Comparison "+quantity);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, range, range);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, false, true, false, false);
	}

}
