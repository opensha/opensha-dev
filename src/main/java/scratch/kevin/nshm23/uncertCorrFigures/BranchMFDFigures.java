package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet.RuptureProperties;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class BranchMFDFigures {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/data/kevin/nshm23/batch_inversions/");
		
		File corrSolDir = new File(invsDir, "2024_02_02-nshm23_branches-WUS_FM_v3");
		File corrBAFile = new File(corrSolDir, "results_WUS_FM_v3_branch_averaged.zip"); // can't be gridded, need order
		
		File avgSolDir = new File(invsDir, "2024_05_07-nshm23_branches-WUS_FM_v3-AvgSupraB-AvgSeg");
		File avgBAFile = new File(avgSolDir, "results_WUS_FM_v3_branch_averaged.zip"); // can't be gridded, need order
		
		File randSolDir = new File(invsDir, "2023_11_16-nshm23_branches-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File randBAFile = new File(randSolDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"); // can't be gridded, need order
		
		File outputDir = new File("/home/kevin/Documents/papers/2024_nshm23_uncert_correlation/figures/branch_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		LogicTree<?> corrTree = LogicTree.read(new File(corrSolDir, "logic_tree.json"));
		FaultSystemSolution corrBASol = FaultSystemSolution.load(corrBAFile);
		LogicTree<?> avgTree = LogicTree.read(new File(avgSolDir, "logic_tree.json"));
		FaultSystemSolution avgBASol = FaultSystemSolution.load(avgBAFile);
		LogicTree<?> randTree = LogicTree.read(new File(randSolDir, "logic_tree.json"));
		FaultSystemSolution randBASol = FaultSystemSolution.load(randBAFile);

		List<EvenlyDiscretizedFunc> corrMFDs = loadCmlMFDs(corrBASol, corrTree);
		
		// load average and random in the background while we plot correlated
		CompletableFuture<List<EvenlyDiscretizedFunc>> avgMFDsFuture = CompletableFuture.supplyAsync(new Supplier<List<EvenlyDiscretizedFunc>>() {

			@Override
			public List<EvenlyDiscretizedFunc> get() {
				try {
					return loadCmlMFDs(avgBASol, avgTree);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		});
		CompletableFuture<List<EvenlyDiscretizedFunc>> randMFDsFuture = CompletableFuture.supplyAsync(new Supplier<List<EvenlyDiscretizedFunc>>() {

			@Override
			public List<EvenlyDiscretizedFunc> get() {
				try {
					return loadCmlMFDs(randBASol, randTree);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		});
		
		EvenlyDiscretizedFunc corrBA = branchAvgMFD(corrMFDs, corrTree);
		
		NSHM23_SegmentationModels[] segModels = {
				NSHM23_SegmentationModels.NONE,
				NSHM23_SegmentationModels.LOW,
				NSHM23_SegmentationModels.MID,
				NSHM23_SegmentationModels.HIGH,
				NSHM23_SegmentationModels.CLASSIC
		};
		SupraSeisBValues[] bVals = {
				SupraSeisBValues.B_0p0,
				SupraSeisBValues.B_0p25,
				SupraSeisBValues.B_0p5,
				SupraSeisBValues.B_0p75,
				SupraSeisBValues.B_1p0
		};
		
		CPT segCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().trim(0d, 0.9d);
		CPT bCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().trim(0d, 0.95d);
		
		double obsRateM5 = Double.NaN;
		double obsBVal = Double.NaN;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.95);
		UncertainBoundedIncrMagFreqDist observedIncr = NSHM23_RegionalSeismicity.getRemapped(NSHM23_RegionLoader.loadFullConterminousWUS(),
					NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
		observedIncr.setName("Observed");
		observedIncr.setBoundName("95% Bounds");
			
		obsRateM5 = observedIncr.getCumRate(observedIncr.getClosestXIndex(5.01));
		obsBVal = ((GutenbergRichterMagFreqDist)NSHM23_RegionalSeismicity.PREFFERRED.build(
					NSHM23_RegionLoader.SeismicityRegions.CONUS_WEST, refMFD, 8.95)).get_bValue();
		
		new Plot(corrMFDs, corrBA).plotObsMFD(obsRateM5, obsBVal).plot(outputDir, "correlated", "Fully Correlated");
		
		List<EvenlyDiscretizedFunc> segMFDs = nodeBranchAvgMFDs(corrMFDs, corrTree, segModels);
		new Plot(corrMFDs, corrBA).plotColoredAverages(segMFDs, segModels, segCPT)
			.plot(outputDir, "correlated_seg", "Segmentation Models");
		List<EvenlyDiscretizedFunc> bMFDs = nodeBranchAvgMFDs(corrMFDs, corrTree, bVals);
		new Plot(corrMFDs, corrBA).plotColoredAverages(bMFDs, bVals, bCPT)
			.plot(outputDir, "correlated_bval", "b-values");
		List<EvenlyDiscretizedFunc> bSegCombos = nodePairwiseBranchAvgMFDs(corrMFDs, corrTree, segModels, bVals);
		new Plot(corrMFDs, corrBA).plotThinAverages(bSegCombos, "Seg & b-Value Combinations")
			.plot(outputDir, "correlated_seg_bval", "Branch Combinations");
		
		CPT endMemberCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().trim(0.15, 0.85);
		EvenlyDiscretizedFunc noneB0 = branchAvgMFD(corrMFDs, corrTree, NSHM23_SegmentationModels.NONE, SupraSeisBValues.B_0p0);
		EvenlyDiscretizedFunc classicB1 = branchAvgMFD(corrMFDs, corrTree, NSHM23_SegmentationModels.CLASSIC, SupraSeisBValues.B_1p0);
		new Plot(corrMFDs, corrBA).plotThinAverages(bSegCombos, "Seg & b-Value Combinations")
			.plotHighlight(noneB0, "None, b=0.0", endMemberCPT.getMinColor())
			.plotHighlight(classicB1, "Classic, b=1.0", endMemberCPT.getMaxColor())
			.plot(outputDir, "correlated_seg_bval_highlight", "Branch Combinations");
		
		List<EvenlyDiscretizedFunc> avgMFDs = avgMFDsFuture.join();
		EvenlyDiscretizedFunc avgBA = branchAvgMFD(avgMFDs, avgTree);
		
		new Plot(avgMFDs, avgBA).plotObsMFD(obsRateM5, obsBVal).plotComparison(corrBA, "Fully Correlated")
			.plot(outputDir, "average", "Average Seg & b-Value");
		
		List<EvenlyDiscretizedFunc> randMFDs = randMFDsFuture.join();
		EvenlyDiscretizedFunc randBA = branchAvgMFD(randMFDs, randTree);
		
		new Plot(randMFDs, randBA).plotObsMFD(obsRateM5, obsBVal).plotComparison(corrBA, "Fully Correlated")
			.plot(outputDir, "random", "Randomly Sampled");
	}
	
	private static List<EvenlyDiscretizedFunc> loadCmlMFDs(SolutionLogicTree slt, IncrementalMagFreqDist refIncrMFD) throws IOException {
		LogicTree<?> tree = slt.getLogicTree();
		List<EvenlyDiscretizedFunc> ret = new ArrayList<>(tree.size());
		for (int i=0; i<tree.size(); i++) {
			LogicTreeBranch<?> branch = tree.getBranch(i);
			System.out.println("Branch "+i+"/"+tree.size());
			RuptureProperties props = slt.loadPropsForBranch(branch);
			double[] rates = slt.loadRatesForBranch(branch);
			IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refIncrMFD.getMinX(), refIncrMFD.size(), refIncrMFD.getDelta());
			for (int r=0; r<rates.length; r++)
				if (rates[r] > 0)
					mfd.add(mfd.getClosestXIndex(props.mags[r]), rates[r]);
			ret.add(mfd.getCumRateDistWithOffset());
		}
		return ret;
	}
	
	private static List<EvenlyDiscretizedFunc> loadCmlMFDs(FaultSystemSolution baSol, LogicTree<?> tree) throws IOException {
		BranchRegionalMFDs branchMFDs = baSol.requireModule(BranchRegionalMFDs.class);
		BranchAveragingOrder branchOrder = baSol.requireModule(BranchAveragingOrder.class);
		
		IncrementalMagFreqDist[] allBranchMFDs = branchMFDs.getSupraTotalBranchMFDs();
		List<EvenlyDiscretizedFunc> ret = new ArrayList<>(tree.size());
		for (LogicTreeBranch<?> branch : tree) {
			int index = branchOrder.getBranchAveragingIndex(branch);
			ret.add(allBranchMFDs[index].getCumRateDistWithOffset());
		}
		return ret;
	}
	
	private static List<EvenlyDiscretizedFunc> nodeBranchAvgMFDs(List<EvenlyDiscretizedFunc> mfds, LogicTree<?> tree,
			LogicTreeNode[] nodes) {
		List<EvenlyDiscretizedFunc> ret = new ArrayList<>();
		for (LogicTreeNode node : nodes)
			ret.add(branchAvgMFD(mfds, tree, node));
		return ret;
	}
	
	private static List<EvenlyDiscretizedFunc> nodePairwiseBranchAvgMFDs(List<EvenlyDiscretizedFunc> mfds,
			LogicTree<?> tree, LogicTreeNode[] nodes1, LogicTreeNode[] nodes2) {
		List<EvenlyDiscretizedFunc> ret = new ArrayList<>();
		for (LogicTreeNode node1 : nodes1)
			for (LogicTreeNode node2 : nodes2)
				ret.add(branchAvgMFD(mfds, tree, node1, node2));
		return ret;
	}
	
	private static EvenlyDiscretizedFunc branchAvgMFD(List<EvenlyDiscretizedFunc> mfds, LogicTree<?> tree,
			LogicTreeNode... requiredNodes) {
		EvenlyDiscretizedFunc ret = null;
		double weightSum = 0d;
		Preconditions.checkState(mfds.size() == tree.size());
		
		for (int i=0; i<tree.size(); i++) {
			LogicTreeBranch<?> branch = tree.getBranch(i);
			if (requiredNodes != null && requiredNodes.length > 0) {
				boolean skip = false;
				for (LogicTreeNode required : requiredNodes) {
					if (!branch.hasValue(required)) {
						skip = true;
						break;
					}
				}
				if (skip)
					continue;
			}
			double weight = tree.getBranchWeight(i);
			weightSum += weight;
			
			EvenlyDiscretizedFunc mfd = mfds.get(i);
			
			if (ret == null)
				ret = new EvenlyDiscretizedFunc(mfd.getMinX(), mfd.size(), mfd.getDelta());
			else
				Preconditions.checkState(mfd.getMinX() == ret.getMinX() && mfd.size() == ret.size());
			
			for (int j=0; j<mfd.size(); j++) {
				double val = mfd.getY(j);
				if (val > 0d)
					ret.add(j, val*weight);
			}
		}
		Preconditions.checkState(weightSum > 0d);
		ret.scale(1d/weightSum);
		return ret;
	}
	
	private static class Plot {
		
		private List<EvenlyDiscretizedFunc> mfds;
		private EvenlyDiscretizedFunc baMFD;
		
		private List<EvenlyDiscretizedFunc> thinAverages;
		private String thinAveragesLabel;
		
		private List<EvenlyDiscretizedFunc> coloredAverages;
		private LogicTreeNode[] coloredNodes;
		private CPT cpt;
		
		private List<EvenlyDiscretizedFunc> highlights;
		private List<String> highlightNames;
		private List<Color> highlightColors;
		
		private EvenlyDiscretizedFunc compMFD;
		private String compName;
		
		private double obsN5;
		private double obsB;

		public Plot(List<EvenlyDiscretizedFunc> mfds, EvenlyDiscretizedFunc baMFD) {
			this.mfds = mfds;
			this.baMFD = baMFD;
		}
		
		public Plot plotThinAverages(List<EvenlyDiscretizedFunc> averages, String label) {
			this.thinAverages = averages;
			this.thinAveragesLabel = label;
			return this;
		}
		
		public Plot plotColoredAverages(List<EvenlyDiscretizedFunc> averages, LogicTreeNode[] nodes, CPT cpt) {
			this.coloredAverages = averages;
			this.coloredNodes = nodes;
			this.cpt = cpt;
			return this;
		}
		
		public Plot plotComparison(EvenlyDiscretizedFunc compMFD, String compName) {
			this.compMFD = compMFD;
			this.compName = compName;
			return this;
		}
		
		public Plot plotHighlight(EvenlyDiscretizedFunc mfd, String name, Color color) {
			if (highlights == null) {
				highlights = new ArrayList<>();
				highlightNames = new ArrayList<>();
				highlightColors = new ArrayList<>();
			}
			highlights.add(mfd);
			highlightNames.add(name);
			highlightColors.add(color);
			return this;
		}
		
		public Plot plotObsMFD(double n5, double b) {
			this.obsN5 = n5;
			this.obsB = b;
			return this;
		}
		
		public void plot(File outputDir, String prefix, String title) throws IOException {
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			PlotCurveCharacterstics indvChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(215, 215, 215));
			for (int i=0; i<mfds.size(); i++) {
				EvenlyDiscretizedFunc mfd = mfds.get(i);
				mfd.setName(null);
				if (i == 0) {
					mfd = mfd.deepClone();
					mfd.setName("Individual Branches ("+mfds.size()+")");
				}
				funcs.add(mfd);
				chars.add(indvChar);
			}
			
			if (thinAverages != null) {
				PlotCurveCharacterstics thinChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY);
				for (int i=0; i<thinAverages.size(); i++) {
					EvenlyDiscretizedFunc mfd = thinAverages.get(i);
					if (i == 0) {
						mfd = mfd.deepClone();
						mfd.setName(thinAveragesLabel+" ("+thinAverages.size()+")");
					}
					funcs.add(mfd);
					chars.add(thinChar);
				}				
			}
			
			if (coloredAverages != null) {
				Preconditions.checkState(coloredAverages.size() > 1);
				cpt = cpt.rescale(0d, coloredAverages.size()-1d);
				for (int i=0; i<coloredAverages.size(); i++) {
					EvenlyDiscretizedFunc mfd = coloredAverages.get(i);
					String name = coloredNodes[i].getShortName();
					if (name.endsWith("Seg"))
						name = name.substring(0, name.length()-3);
					mfd.setName(name);
					Color color = cpt.getColor((float)i);
					funcs.add(mfd);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
				}				
			}
			
			if (highlights != null) {
				for (int i=0; i<highlights.size(); i++) {
					EvenlyDiscretizedFunc mfd = highlights.get(i);
					String name = highlightNames.get(i);
					if (name.endsWith("Seg"))
						name = name.substring(0, name.length()-3);
					mfd.setName(name);
					Color color = highlightColors.get(i);
					funcs.add(mfd);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
				}
			}
			
			if (compMFD != null) {
				compMFD = compMFD.deepClone();
				compMFD.setName(compName);;
				
				funcs.add(compMFD);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.DARK_GRAY));				
			}
			
			Range xRange = new Range(6d, 8.5);
			Range yRange = new Range(1e-5, 2e0);
			
			if (obsN5 > 0d) {
				Color obsColor = new Color(125, 80, 145); // "indigo"
				EvenlyDiscretizedFunc refFunc = FaultSysTools.initEmptyMFD(xRange.getUpperBound()+0.05);
				GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refFunc.getMinX(), refFunc.size(), refFunc.getDelta());
				gr.setAllButTotMoRate(5.05, gr.getMaxX(), obsN5, obsB);
				// now make it cumulative, but without the roll off
				gr.scaleToIncrRate(5.05d, obsN5);
				
				gr.setName("Observed (Extrapolation)");
				funcs.add(gr);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, obsColor));
				
			}
			
			baMFD.setName("Branch Averaged");
			funcs.add(baMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Magnitude", "Cumulative Rate (1/yr)");
			spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.getPlotPrefs().scaleFontSizes(1.25);
			
			gp.drawGraphPanel(spec, false, true, xRange, yRange);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
		}
	}

}
