package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;

import net.mahdilamb.colormap.Colors;

public class SlipRateFigures {
	
	static Region CRUSTAL_FAULT_MAP_REG = new Region(new Location(16.4, -70.2), new Location(20.2, -61.7));;

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/data/kevin/nshm23/batch_inversions");
		File figsDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures");
		File crustalOutputDir = new File(figsDir, "crustal_dm");
		
		File crustalSolDir = new File(invsDir, "2024_12_05-prvi25_crustal_branches-dmSample5x/");
		FaultSystemSolution crustalSol = FaultSystemSolution.load(
				new File(crustalSolDir, "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		FaultSystemSolution crustalNoClassicSol = buildNoClassic(new File(crustalSolDir, "node_branch_averaged"));
		
		plotCrustal(crustalOutputDir, crustalSol, crustalNoClassicSol);
//		plotCrustal(outputDir, null);
	}
	
	private static void plotCrustal(File outputDir, FaultSystemSolution sol, FaultSystemSolution solNoClassic) throws IOException {
		PRVI25_CrustalFaultModels fm = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1;
		
		int numTotSects = 0;
		int numProxySects = 0;
		for (FaultSection sect : fm.getFaultSections()) {
			numTotSects++;
			if (sect.isProxyFault())
				numProxySects++;
		}
		
		List<? extends FaultSection> origSubSects = PRVI25_CrustalDeformationModels.GEOLOGIC.build(fm);
		List<? extends FaultSection> targetSubSects = PRVI25_CrustalDeformationModels.GEOLOGIC_DIST_AVG.build(fm);
		System.out.println("Model has "+numTotSects+" sections ("+(numTotSects-numProxySects)+" regular, "
				+numProxySects+" proxy) and "+origSubSects.size()+" subsections");
		GeographicMapMaker mapMaker = new GeographicMapMaker(CRUSTAL_FAULT_MAP_REG, origSubSects);
		mapMaker.setWriteGeoJSON(false);
		
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 11d);
		slipCPT.setPreferredTickInterval(2d);
		
		for (boolean target : new boolean[] {false,true}) {
			String prefix = target ? "crustal_target_slip_rates" : "crustal_original_slip_rates";
			List<? extends FaultSection> subSects = target ? targetSubSects : origSubSects;
			List<Double> slipRates = getDMSlipRates(subSects, false);
			System.out.println("Max crustal is "+slipRates.stream().mapToDouble(D->D).max().getAsDouble());
			
			String label = (target ? "Target " : "")+"Slip Rate (mm/yr)";
			mapMaker.plotSectScalars(slipRates, slipCPT, label);
			
			mapMaker.plot(outputDir, prefix, " ");
			mapMaker.setCPTLocation(RectangleEdge.TOP);
			PlotSpec linearPlot = mapMaker.buildPlot(" ");
			mapMaker.setCPTLocation(RectangleEdge.BOTTOM);
			
			CPT logSlipCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(-1, 1);
//			slipCPT.setPreferredTickInterval(2d);
			
			List<Double> logSlipRates = getDMSlipRates(subSects, true);
			
			mapMaker.plotSectScalars(logSlipRates, logSlipCPT, "Log10 "+label);
			
			mapMaker.plot(outputDir, prefix+"_log", " ");
			PlotSpec logPlot = mapMaker.buildPlot(" ");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(List.of(linearPlot, logPlot), false, false, List.of(mapMaker.getXRange()),
					List.of(mapMaker.getYRange(), mapMaker.getYRange()));
			
			PlotUtils.writePlots(outputDir, prefix+"_combined", gp, mapMaker.getDefaultPlotWidth(), true, true, true, false);
		}

		FaultSystemRupSet rupSet = sol.getRupSet();
		List<Double> targetSlips = new ArrayList<>();
		Map<Integer, List<Double>> targetSlipsByParent = new HashMap<>();
		SectSlipRates rupSetTargets = rupSet.getSectSlipRates();
		for (int i=0; i<rupSetTargets.size(); i++) {
			double slip = rupSetTargets.getSlipRate(i)*1e3;
			targetSlips.add(slip);
			int parentID = rupSet.getFaultSectionData(i).getParentSectionId();
			if (!targetSlipsByParent.containsKey(parentID))
				targetSlipsByParent.put(parentID, new ArrayList<>());
			targetSlipsByParent.get(parentID).add(slip);
		}
		Map<Integer, Double> targetAverageSlipsByParent = targetSlipsByParent.entrySet().stream()
				.collect(Collectors.toMap(Map.Entry::getKey,
						entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).average().orElse(0.0)
						));
		
		for (boolean includeClassic : new boolean[] {false,true}) {
			String prefix = "crustal_solution_slip_rates";
			String solName;
			SolutionSlipRates solSlipsModule;
			if (includeClassic) {
				solSlipsModule = sol.requireModule(SolutionSlipRates.class);
				solName = "Solution";
			} else {
				solSlipsModule = solNoClassic.requireModule(SolutionSlipRates.class);
				prefix += "_no_classic";
				solName = "Solution (excl. classic)";
			}
			List<Double> solSlips = new ArrayList<>();
			for (int i=0; i<targetSubSects.size(); i++)
				solSlips.add(solSlipsModule.get(i)*1e3);
			Map<Integer, List<Double>> solSlipsByParent = new HashMap<>();
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
				int parentID = sect.getParentSectionId();
				if (!solSlipsByParent.containsKey(parentID))
					solSlipsByParent.put(parentID, new ArrayList<>());
				solSlipsByParent.get(parentID).add(solSlipsModule.get(sect.getSectionId())*1e3);
			}
			Map<Integer, Double> solAverageSlipsByParent = solSlipsByParent.entrySet().stream()
					.collect(Collectors.toMap(Map.Entry::getKey,
							entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).average().orElse(0.0)
							));
//			List<Double> solSlips = new ArrayList<>();
//			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
//				int parentID = sect.getParentSectionId();
//				solSlips.add(solAverageSlipsByParent.get(parentID));
//				if (sect.getSectionId() > 0 && parentID != rupSet.getFaultSectionData(sect.getSectionId()-1).getParentSectionId())
//					System.out.println(sect.getParentSectionName()+": "+solAverageSlipsByParent.get(parentID));
//			}
			
			mapMaker.plotSectScalars(solSlips, slipCPT, solName+" Slip Rates (mm/yr)");
			mapMaker.plot(outputDir, prefix, " ");
			
			CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1d, 1d);
			CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-20d, 20d);
			
			mapMaker.plotSectScalars(difference(targetSlips, solSlips), diffCPT, solName+" - Target Slip Rates (mm/yr)");
			mapMaker.plot(outputDir, prefix+"_diff", " ");
			
			mapMaker.plotSectScalars(percentDifference(targetSlips, solSlips), pDiffCPT, solName+" vs Target Slip Rates (% Difference)");
			mapMaker.plot(outputDir, prefix+"_pDiff", " ");
			
			plotScatter(outputDir, prefix+"_scatter", targetSlips, solSlips, targetAverageSlipsByParent, solAverageSlipsByParent);
		}
	}
	
	private static FaultSystemSolution buildNoClassic(File nodeSolDir) throws IOException {
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.add(PRVI25_LogicTreeBranch.SEG);
		for (NSHM23_SegmentationModels segModel : NSHM23_SegmentationModels.values()) {
			if (segModel.getNodeWeight(null) == 0d || segModel == NSHM23_SegmentationModels.CLASSIC)
				continue;
			File solFile = new File(nodeSolDir, "SegModel_"+segModel.getFilePrefix()+".zip");
			LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
			branch.setValue(segModel);
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			baCreator.addSolution(sol, branch);
		}
		return baCreator.build();
	}
	
	private static List<Double> getDMSlipRates(List<? extends FaultSection> subSects, boolean log) {
		List<Double> slipRates = new ArrayList<>();
		
		for (FaultSection sect : subSects) {
			double slip = sect.getReducedAveSlipRate();
			if (log)
				slip = Math.log10(slip);
			slipRates.add(slip);
		}
		
		return slipRates;
	}
	
	private static List<Double> difference(List<Double> target, List<Double> solution) {
		List<Double> ret = new ArrayList<>();
		for (int i=0; i<target.size(); i++)
			ret.add(solution.get(i) - target.get(i));
		return ret;
	}
	
	private static List<Double> percentDifference(List<Double> target, List<Double> solution) {
		List<Double> ret = new ArrayList<>();
		for (int i=0; i<target.size(); i++)
			ret.add(100d*(solution.get(i) - target.get(i))/target.get(i));
		return ret;
	}
	
	private static void plotScatter(File outputDir, String prefix, List<Double> targetSlips, List<Double> solSlips,
			Map<Integer, Double> targetParentSlips, Map<Integer, Double> solParentSlips) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range linearRange = new Range(0d, 13d);
		Range logRange = new Range(1e-1, 2e1);
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(linearRange.getLowerBound(), linearRange.getLowerBound());
		oneToOne.set(logRange.getLowerBound(), logRange.getLowerBound());
		oneToOne.set(linearRange.getUpperBound(), linearRange.getUpperBound());
		oneToOne.set(logRange.getUpperBound(), logRange.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		if (solSlips != null) {
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			for (int i=0; i<solSlips.size(); i++)
				scatter.set(targetSlips.get(i), solSlips.get(i));
			
			funcs.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 4f, Color.BLACK));
		}
		
		if (solParentSlips != null) {
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			for (int id : solParentSlips.keySet())
				scatter.set(targetParentSlips.get(id), solParentSlips.get(id));
			
			funcs.add(scatter);
			Color color = Colors.tab_blue;
//			color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 180);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, color));
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, false, linearRange, linearRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, false, true, true, false);
		
		gp.drawGraphPanel(plot, true, true, logRange, logRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_log", gp, 800, false, true, true, false);
	}

}
