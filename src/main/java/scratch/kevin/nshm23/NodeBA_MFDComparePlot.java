package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.MaxJumpDistModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;

public class NodeBA_MFDComparePlot {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/");
		
		LogicTreeLevel<?> level = NSHM23_LogicTreeBranch.SUPRA_B;
		File invDir = new File(mainDir, 
//				"2022_05_09-nshm23_u3_hybrid_branches-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvg");
				"2022_05_27-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg");
		
//		LogicTreeLevel<?> level = LogicTreeLevel.forEnum(MaxJumpDistModels.class, "Max Dist Segmentation", "MaxDist");
////		LogicTreeLevel<?> level = NSHM23_LogicTreeBranch.MAX_DIST;
//		File invDir = new File(mainDir, 
//				"2022_05_09-nshm23_u3_hybrid_branches-strict_cutoff_seg-shift_seg_1km-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1");
		
		File baSolFile = new File(invDir, "results_FM3_1_CoulombRupSet_branch_averaged.zip");
		
		File nodeBAdir = new File(invDir, "node_branch_averaged");
		
		FaultSystemSolution baSol = FaultSystemSolution.load(baSolFile);
		
		Range xRange = new Range(5.8d, 8.8d);
		Range yRange = new Range(1e-5, 1e-1);
		Range yRangeCml = new Range(1e-5, 1e0);
		IncrementalMagFreqDist[] dataMFDs = {
				new GutenbergRichterMagFreqDist(1.0, TotalMag5Rate.RATE_6p5.getRateMag5(), 5.05, 9.95, 50),
				new GutenbergRichterMagFreqDist(1.0, TotalMag5Rate.RATE_7p9.getRateMag5(), 5.05, 9.95, 50),
				new GutenbergRichterMagFreqDist(1.0, TotalMag5Rate.RATE_9p6.getRateMag5(), 5.05, 9.95, 50)
		};
		
		PlotCurveCharacterstics dataChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.GREEN.darker());
		
		Comparator<LogicTreeNode> nodeComp = new Comparator<LogicTreeNode>() {
			
			@Override
			public int compare(LogicTreeNode o1, LogicTreeNode o2) {
				if (o1.getClass().isEnum())
					return Integer.compare(((Enum<?>)o1).ordinal(), ((Enum<?>)o2).ordinal());
				return o1.getShortName().compareTo(o2.getShortName());
			}
		};
		
		List<FaultSystemSolution> nodeSols = new ArrayList<>();
		List<LogicTreeNode> nodeVals = new ArrayList<>();
		String solPrefix = level.getShortName().replaceAll("\\W+", "_");
		for (File file : nodeBAdir.listFiles()) {
			String name = file.getName();
			if (file.isFile() && name.startsWith(solPrefix) && name.endsWith(".zip")) {
				FaultSystemSolution sol = FaultSystemSolution.load(file);
				LogicTreeBranch<?> branch = sol.requireModule(LogicTreeBranch.class);
				LogicTreeNode value = branch.getValue(level.getType());
				if (nodeSols.isEmpty()) {
					nodeSols.add(sol);
					nodeVals.add(value);
				} else {
					// sort
					for (int i=0; i<=nodeVals.size(); i++) {
						if (i == nodeVals.size()) {
							nodeSols.add(sol);
							nodeVals.add(value);
							break;
						} else {
							LogicTreeNode oValue = nodeVals.get(i);
							int cmp = nodeComp.compare(value, oValue);
							if (cmp < 0) {
								nodeSols.add(i, sol);
								nodeVals.add(i, value);
								break;
							}
						}
					}
				}
			}
		}
		
		double weight0 = nodeVals.get(0).getNodeWeight(null);
		boolean allSameWeight = true;
		for (LogicTreeNode node : nodeVals)
			allSameWeight = allSameWeight && (float)weight0 == (float)node.getNodeWeight(null);
		
		Preconditions.checkState(!nodeVals.isEmpty(), "No matches found starting with %s", solPrefix);
		
		System.out.println("Loaded "+nodeVals.size()+" node BA solutions:");
		for (LogicTreeNode value : nodeVals)
			System.out.println("\t"+value.getName());
		
		CPT bCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, nodeVals.size()-1);
		
		FaultSystemSolution straightAvg = AverageSolutionCreator.buildAverage(nodeSols.toArray(new FaultSystemSolution[0]));
		Preconditions.checkState(straightAvg.hasModule(RupMFDsModule.class));
		
		List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (IncrementalMagFreqDist dataMFD : dataMFDs) {
			if (incrFuncs.isEmpty())
				dataMFD.setName("UCERF3 Data Constraints");
			else
				dataMFD.setName(null);
			incrFuncs.add(dataMFD);
			chars.add(dataChar);
		}
		
		IncrementalMagFreqDist baMFD = loadMFD(baSol, xRange);
		baMFD.setName("Full Branch-Average");
		
		incrFuncs.add(baMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		if (!allSameWeight) {
			IncrementalMagFreqDist avgMFD = loadMFD(straightAvg, xRange);
			avgMFD.setName("Even-Weight Average");
			
			incrFuncs.add(avgMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		}
		
		for (int i=0; i<nodeVals.size(); i++) {
			LogicTreeNode value = nodeVals.get(i);
			IncrementalMagFreqDist mfd = loadMFD(nodeSols.get(i), xRange);
			
			mfd.setName(value.getShortName());
			incrFuncs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, bCPT.getColor((float)i)));
		}
		
//		if (baSol.getRupSet().hasModule(InversionTargetMFDs.class)) {
//			InversionTargetMFDs targets = baSol.getRupSet().getModule(InversionTargetMFDs.class);
//			
//			IncrementalMagFreqDist target = targets.getTotalOnFaultSupraSeisMFD();
//			if (target != null) {
//				target.setName("BA Target");
//				Color color = SolMFDPlot.SUPRA_SEIS_TARGET_COLOR;
//				incrFuncs.add(target);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
//				
//				if (target instanceof UncertainIncrMagFreqDist) {
//					UncertainBoundedIncrMagFreqDist sigmaIncrBounds =
//							((UncertainIncrMagFreqDist)target).estimateBounds(UncertaintyBoundType.ONE_SIGMA);
//					sigmaIncrBounds.setName("± σ");
//					
//					incrFuncs.add(sigmaIncrBounds);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
//							new Color(color.getRed(), color.getGreen(), color.getBlue(), 60)));
//				}
//			}
//		}
		
		PlotSpec spec = new PlotSpec(incrFuncs, chars, level.getName()+" Sweep", "Magnitude", "Incremental Rate (/yr)");
		spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		String prefix = solPrefix+"_mfds";
		PlotUtils.writePlots(nodeBAdir, prefix, gp, 900, 800, true, true, false);
		
		prefix += "_cml";
		
		List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
		for (int j=0; j<incrFuncs.size(); j++) {
			IncrementalMagFreqDist incr = incrFuncs.get(j);
			PlotCurveCharacterstics pChar = chars.get(j);
			//				if (incr instancoe)
			if (pChar.getLineType() == PlotLineType.SHADED_UNCERTAIN) {
				Preconditions.checkState(incr instanceof UncertainBoundedIncrMagFreqDist);
				UncertainBoundedIncrMagFreqDist bounded = (UncertainBoundedIncrMagFreqDist)incr;
				
				EvenlyDiscretizedFunc cumulative = bounded.getCumRateDistWithOffset();
				
				EvenlyDiscretizedFunc upperCumulative = bounded.getUpper().getCumRateDistWithOffset();
				EvenlyDiscretizedFunc lowerCumulative = bounded.getLower().getCumRateDistWithOffset();
				Preconditions.checkState(cumulative.size() == upperCumulative.size());
				for (int k=0; k<cumulative.size(); k++) {
					upperCumulative.set(k, Math.max(cumulative.getY(k), upperCumulative.getY(k)));
					lowerCumulative.set(k, Math.max(0, Math.min(cumulative.getY(k), lowerCumulative.getY(k))));
				}
				
				UncertainArbDiscFunc cmlBounded = new UncertainArbDiscFunc(cumulative, lowerCumulative, upperCumulative);
				cmlBounded.setName("± σ");
				cmlFuncs.add(cmlBounded);
			} else {
				cmlFuncs.add(incr.getCumRateDistWithOffset());
			}
		}
		spec = new PlotSpec(cmlFuncs, chars, level.getName()+" Sweep", "Magnitude", "Cumulative Rate (/yr)");
		spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRangeCml);
		
		PlotUtils.writePlots(nodeBAdir, prefix, gp, 900, 800, true, true, false);
	}
	
	private static IncrementalMagFreqDist loadMFD(FaultSystemSolution sol, Range xRange) {
		IncrementalMagFreqDist defaultMFD = SolMFDPlot.initDefaultMFD(xRange.getLowerBound()+0.01, xRange.getUpperBound()-0.01);
		
		return sol.calcTotalNucleationMFD(defaultMFD.getMinX(), defaultMFD.getMaxX(), defaultMFD.getDelta());
	}

}
