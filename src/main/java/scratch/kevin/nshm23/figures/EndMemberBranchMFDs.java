package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class EndMemberBranchMFDs {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File modelDir = new File(invsDir,
				"2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		FaultSystemSolution baSol = FaultSystemSolution.load(new File(modelDir,
				"results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		LogicTree<?> tree = LogicTree.read(new File(modelDir, "logic_tree.json"));
		
		double minMag = 6d;
		
		BranchAveragingOrder order = baSol.requireModule(BranchAveragingOrder.class);
		BranchRegionalMFDs regMFDs = baSol.requireModule(BranchRegionalMFDs.class);
		RegionsOfInterest roi = baSol.getRupSet().getModule(RegionsOfInterest.class);
		
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8.45);
		
		Region fullRegion = NSHM23_RegionLoader.loadFullConterminousWUS();
		IncrementalMagFreqDist obsMFD = NSHM23_RegionalSeismicity.getRemapped(fullRegion,
				NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, 8.55);
		
		List<String> branchNames = new ArrayList<>();
		List<LogicTreeNode[]> branchNodes = new ArrayList<>();
		List<Color> branchColors = new ArrayList<>();
		
		branchNames.add("Classic (all)");
		branchColors.add(Color.RED);
		branchNodes.add(new LogicTreeNode[] {
				NSHM23_SegmentationModels.CLASSIC,
		});
		
		branchNames.add("Classic (b=0)");
		branchColors.add(Color.MAGENTA);
		branchNodes.add(new LogicTreeNode[] {
				NSHM23_SegmentationModels.CLASSIC,
				SupraSeisBValues.B_0p0
		});
		
		branchNames.add("No Segmentaion");
		branchColors.add(Color.BLUE);
		branchNodes.add(new LogicTreeNode[] {
				NSHM23_SegmentationModels.NONE
		});
		
		branchNames.add("No Segmentaion (b=1)");
		branchColors.add(Color.CYAN);
		branchNodes.add(new LogicTreeNode[] {
				NSHM23_SegmentationModels.NONE,
				SupraSeisBValues.B_1p0
		});
		
		SummedMagFreqDist[] sumMFDs = new SummedMagFreqDist[tree.size()];
		
		Region[] sumRegions = {
				AnalysisRegions.CONUS_U3_RELM.load(),
				AnalysisRegions.CONUS_IMW.load(),
				AnalysisRegions.CONUS_PNW.load()
		};
		
		for (Region region : sumRegions) {
			int regionIndex = -1;
			List<Region> allRegions = roi.getRegions();
			for (int i=0; i<allRegions.size(); i++) {
				Region testReg = allRegions.get(i);
				if (testReg.equalsRegion(region)) {
					regionIndex = i;
					break;
				}
			}
			System.out.println(region.getName()+" is index "+regionIndex);
			Preconditions.checkState(regionIndex >= 0);
			IncrementalMagFreqDist[] mfds = regMFDs.getRegionalBranchMFDs(MFDType.SUPRA_ONLY, regionIndex);
			
			for (int i=0; i<mfds.length; i++) {
				if (sumMFDs[i] == null) {
					sumMFDs[i] = new SummedMagFreqDist(mfds[i].getMinX(), mfds[i].getMaxX(), mfds[i].size());
				} else {
					Preconditions.checkState(sumMFDs[i].getMinX() == mfds[i].getMinX());
					Preconditions.checkState(sumMFDs[i].size() == mfds[i].size());
				}
				sumMFDs[i].addIncrementalMagFreqDist(mfds[i]);
			}
		}
		
		refMFD = sumMFDs[0];
		
		SummedMagFreqDist avgMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		SummedMagFreqDist[] branchMFDs = new SummedMagFreqDist[branchNames.size()];
		double[] branchSumWeights = new double[branchNames.size()];
		for (int i=0; i<branchNames.size(); i++)
			branchMFDs[i] = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		double sumWeight = 0d;
		
		for (LogicTreeBranch<?> branch : tree) {
			double weight = tree.getBranchWeight(branch);
			
			int index = order.getBranchAveragingIndex(branch);
			IncrementalMagFreqDist mfd = sumMFDs[index];
			
			mfd = mfd.deepClone();
			mfd.scale(weight);
			
			avgMFD.addIncrementalMagFreqDist(mfd);
			sumWeight += weight;
			
			for (int i=0; i<branchMFDs.length; i++) {
				if (matches(branch, branchNodes.get(i))) {
					branchMFDs[i].addIncrementalMagFreqDist(mfd);
					branchSumWeights[i] += weight;
				}
			}
		}
		
		avgMFD.scale(1d/sumWeight);
		for (int i=0; i<branchMFDs.length; i++)
			branchMFDs[i].scale(1d/branchSumWeights[i]);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		obsMFD.setName("Observed (extrapolated)");
		funcs.add(obsMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		avgMFD.setName("Full Branch Average");
		funcs.add(avgMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		for (int i=0; i<branchNames.size(); i++) {
			IncrementalMagFreqDist mfd = branchMFDs[i];
			mfd.setName(branchNames.get(i));
			
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, branchColors.get(i)));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "End Member Comparison", "Magnitude", "Incremental Rate (1/yr)");
		spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		Range magRage = new Range(6d, 8.5d);
		Range yRange = new Range(1e-6, 2e0);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, magRage, yRange);
		
		PlotUtils.writePlots(outputDir, "end_member_mfds", gp, 850, 800, true, true, false);
	}
	
	private static boolean matches(LogicTreeBranch<?> branch, LogicTreeNode... nodes) {
		for (LogicTreeNode node : nodes)
			if (!branch.hasValue(node))
				return false;
		return true;
	}

}
