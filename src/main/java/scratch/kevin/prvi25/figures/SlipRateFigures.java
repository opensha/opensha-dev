package scratch.kevin.prvi25.figures;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
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
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class SlipRateFigures {
	
	public static Region CRUSTAL_FAULT_MAP_REG = new Region(new Location(16.4, -70.2), new Location(20.2, -61.7));

	public static void main(String[] args) throws IOException {
		File crustalDMOutputDir = new File(FIGURES_DIR, "crustal_dm");
		File crustalSolOutputDir = new File(FIGURES_DIR, "crustal_sol");
		
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_SUPRA_ONLY);
		FaultSystemSolution crustalNoClassicSol = buildNoClassic(new File(CRUSTAL_DIR, "node_branch_averaged"));
		
		plotCrustal(crustalDMOutputDir, crustalSolOutputDir, crustalSol, crustalNoClassicSol);
//		plotCrustal(outputDir, null);

		File subSolOutputDir = new File(FIGURES_DIR, "sub_sol");
		FaultSystemSolution subLargeSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution subSmallSol = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		plotSubduction(subSolOutputDir, subSmallSol, subLargeSol);
	}
	
	private static void plotCrustal(File dmOutputDir, File solOutputDir, FaultSystemSolution sol, FaultSystemSolution solNoClassic) throws IOException {
		Preconditions.checkState(dmOutputDir.exists() || dmOutputDir.mkdir());
		Preconditions.checkState(solOutputDir.exists() || solOutputDir.mkdir());
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
		
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 10d);
		slipCPT.setPreferredTickInterval(2d);
		
		for (boolean target : new boolean[] {false,true}) {
			String prefix = target ? "crustal_target_slip_rates" : "crustal_original_slip_rates";
			List<? extends FaultSection> subSects = target ? targetSubSects : origSubSects;
			List<Double> slipRates = getDMSlipRates(subSects);
			System.out.println("Max crustal is "+slipRates.stream().mapToDouble(D->D).max().getAsDouble());
			
			String label = (target ? "Target " : "")+"slip rate (mm/yr)";
			mapMaker.plotSectScalars(slipRates, slipCPT, label);
			
			mapMaker.plot(dmOutputDir, prefix, " ");
			mapMaker.setCPTLocation(RectangleEdge.TOP);
			PlotSpec linearPlot = mapMaker.buildPlot(" ");
			mapMaker.setCPTLocation(RectangleEdge.BOTTOM);
			
//			CPT logSlipCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(-1, 1);
			CPT logSlipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-1, 1);
			logSlipCPT.setLog10(true);
			
			mapMaker.plotSectScalars(slipRates, logSlipCPT, label);
			
			mapMaker.plot(dmOutputDir, prefix+"_log", " ");
			PlotSpec logPlot = mapMaker.buildPlot(" ");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(List.of(linearPlot, logPlot), false, false, List.of(mapMaker.getXRange()),
					List.of(mapMaker.getYRange(), mapMaker.getYRange()));
			
			PlotUtils.writePlots(dmOutputDir, prefix+"_combined", gp, mapMaker.getDefaultPlotWidth(), true, true, true, false);
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
		List<? extends FaultSection> lowSubSects = PRVI25_CrustalDeformationModels.GEOLOGIC_LOW.build(fm);
		List<? extends FaultSection> highSubSects = PRVI25_CrustalDeformationModels.GEOLOGIC_HIGH.build(fm);
		List<Double> binWidths = new ArrayList<>();
		for (int i=0; i<lowSubSects.size(); i++)
			binWidths.add(highSubSects.get(i).getOrigAveSlipRate() - lowSubSects.get(i).getOrigAveSlipRate());
		Map<Integer, List<Double>> binWidthsByParent = new HashMap<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			int parentID = sect.getParentSectionId();
			if (!binWidthsByParent.containsKey(parentID))
				binWidthsByParent.put(parentID, new ArrayList<>());
			binWidthsByParent.get(parentID).add(binWidths.get(sect.getSectionId()));
		}
		Map<Integer, Double> averageBinWidthByParent = binWidthsByParent.entrySet().stream()
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
			Map<Integer, String> parentNames = new HashMap<>();
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
				int parentID = sect.getParentSectionId();
				if (!solSlipsByParent.containsKey(parentID)) {
					solSlipsByParent.put(parentID, new ArrayList<>());
					parentNames.put(parentID, sect.getParentSectionName());
				}
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
			mapMaker.plot(solOutputDir, prefix, " ");
			
			CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1d, 1d);
			CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-20d, 20d);
			
			List<? extends FaultSection> sects = sol.getRupSet().getFaultSectionDataList();
			List<XYTextAnnotation> anns = new ArrayList<>();
			List<LocationList> arrows = new ArrayList<>();
			Font labelFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
			CrustalFaultNamesFigure.buildLabelsAndArrows(sects, "Bouillante Montserrat", false, new Location(16.55, -62.4),
					labelFont, TextAnchor.CENTER_RIGHT, 0, anns, arrows);
			CrustalFaultNamesFigure.buildLabelsAndArrows(sects, "Main Ridge 2", false, new Location(19.2, -65.9),
					labelFont, TextAnchor.TOP_CENTER, 0d, anns, arrows);
			CrustalFaultNamesFigure.buildLabelsAndArrows(sects, "Main Ridge 1", false, new Location(19.8, -65.1),
					labelFont, TextAnchor.BOTTOM_LEFT, 0d, anns, arrows);
			CrustalFaultNamesFigure.buildLabelsAndArrows(sects, "Bunce 6", false, new Location(19.4, -64.6),
					labelFont, TextAnchor.TOP_CENTER, 0d, anns, arrows);
			
			mapMaker.addAnnotations(anns);
			mapMaker.plotArrows(arrows, 6d, new Color(0, 0, 0, 180), 1f);
			mapMaker.setFillArrowheads(true);
			
			mapMaker.plotSectScalars(difference(targetSlips, solSlips), diffCPT, solName+" - target slip rates (mm/yr)");
			mapMaker.plot(solOutputDir, prefix+"_diff", " ");
			
			mapMaker.plotSectScalars(percentDifference(targetSlips, solSlips), pDiffCPT, solName+" vs target slip rates (% Difference)");
			mapMaker.plot(solOutputDir, prefix+"_pDiff", " ");
			
			List<XYAnnotation> linearAnns = circleFaultAnns(targetSubSects, "Main Ridge", solSlips, targetSlips, TextAnchor.TOP_LEFT);
			linearAnns.addAll(circleFaultAnns(targetSubSects, "Bouillante Montserrat", solSlips, targetSlips, TextAnchor.BASELINE_RIGHT));
			linearAnns.addAll(circleFaultAnns(targetSubSects, "Septentrional 1", solSlips, targetSlips, TextAnchor.BASELINE_RIGHT));
//			linearAnns.addAll(circleFaultAnns(targetSubSects, "Septentrional 2", solSlips, targetSlips, TextAnchor.BASELINE_RIGHT));
			
			Range linearRange = new Range(0d, 10d);
			Range logRange = new Range(1e-1, 2e1);
			
			plotScatter(solOutputDir, prefix+"_scatter", targetSlips, solSlips,
					targetAverageSlipsByParent, solAverageSlipsByParent,
					linearRange, logRange, linearAnns, null);
			
			CSVFile<String> slipsCSV = new CSVFile<>(true);
			
			slipsCSV.addLine("", "Absolute slip rate misfit (mm/yr)", "Misfit percentage of rate", "Misfit percentage of bound width");
			DecimalFormat pDF = new DecimalFormat("0.#%");
			DecimalFormat twoDigits = new DecimalFormat("0.00");
			
			FileWriter slipTEX = new FileWriter(new File(solOutputDir, prefix+"_slip_misfits.tex"));
			for (boolean subsection : new boolean[] {true,false}) {
				List<Double> absValues = new ArrayList<>();
				List<Double> scaledToRates = new ArrayList<>();
				List<Double> scaledToBins = new ArrayList<>();
				if (subsection) {
					for (int i=0; i<solSlips.size(); i++) {
						double solVal = solSlips.get(i);
						double targetVal = targetSlips.get(i);
						double binWidth = binWidths.get(i);
						double absDiff = Math.abs(solVal - targetVal);
						absValues.add(absDiff);
						scaledToRates.add(absDiff/targetVal);
						scaledToBins.add(absDiff/binWidth);
					}
				} else {
					CSVFile<String> parentSlipsCSV = new CSVFile<>(true);
					parentSlipsCSV.addLine("Parent Section ID", "Parent Section Name", "Target Slip Rate (mm/yr)",
							"Solution Slip Rate (mm/yr)", "Difference (mm/yr)", "Fractional Mismit", "Bin Fractional Misfit");
					for (int parentID : solAverageSlipsByParent.keySet()) {
						double solVal = solAverageSlipsByParent.get(parentID);
						double targetVal = targetAverageSlipsByParent.get(parentID);
						double binWidth = averageBinWidthByParent.get(parentID);
						double absDiff = Math.abs(solVal - targetVal);
						absValues.add(absDiff);
						scaledToRates.add(absDiff/targetVal);
						scaledToBins.add(absDiff/binWidth);
						
						String parentName = parentNames.get(parentID);
						parentSlipsCSV.addLine(
								parentID+"",
								parentName,
								(float)targetVal+"",
								(float)solVal+"",
								(float)(solVal - targetVal)+"",
								(float)(absDiff/targetVal)+"",
								(float)(absDiff/binWidth)+"");
						
						if (parentName.equals("Septentrional 1")) {
							slipTEX.write(LaTeXUtils.defineValueCommand("SeptentrionalSlipMisfitPct", pDF.format(absDiff/targetVal))+"\n");
						}
					}
					parentSlipsCSV.writeToFile(new File(solOutputDir, prefix+"_parent_slips.csv"));
				}
				for (boolean max : new boolean[] {false,true}) {
					String label = subsection ? "Subsection" : "Fault section aggregated";
					if (max)
						label += ", maximum";
					else
						label += ", average";
					List<String> line = new ArrayList<>();
					line.add(label);
					double abs, scaledToRate, scaledToBin;
					String texPrefix = subsection ? "CrustalSubsect" : "CrustalSect";
					if (!includeClassic)
						texPrefix += "ExclClassic";
					if (max) {
						texPrefix += "Max";
						abs = absValues.stream().mapToDouble(D->D).max().getAsDouble();
						scaledToRate = scaledToRates.stream().mapToDouble(D->D).max().getAsDouble();
						scaledToBin = scaledToBins.stream().mapToDouble(D->D).max().getAsDouble();
					} else {
						texPrefix += "Avg";
						abs = absValues.stream().mapToDouble(D->D).average().getAsDouble();
						scaledToRate = scaledToRates.stream().mapToDouble(D->D).average().getAsDouble();
						scaledToBin = scaledToBins.stream().mapToDouble(D->D).average().getAsDouble();
					}
					texPrefix += "SlipMsft";
					slipTEX.write(LaTeXUtils.defineValueCommand(texPrefix, twoDigits.format(abs))+"\n");
					slipTEX.write(LaTeXUtils.defineValueCommand(texPrefix+"Pct", pDF.format(scaledToRate))+"\n");
					slipTEX.write(LaTeXUtils.defineValueCommand(texPrefix+"BinPct", pDF.format(scaledToBin))+"\n");
					line.add(twoDigits.format(abs));
					line.add(pDF.format(scaledToRate));
					line.add(pDF.format(scaledToBin));
					slipsCSV.addLine(line);
				}
			}
			slipTEX.close();
			slipsCSV.writeToFile(new File(solOutputDir, prefix+"_slip_misfits.csv"));
		}
	}
	
	private static void plotSubduction(File outputDir, FaultSystemSolution smallSol, FaultSystemSolution largeSol) throws IOException {
		double weightLarge = PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE.getNodeWeight(null);
		double weightSmall = PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL.getNodeWeight(null);
		
		Preconditions.checkState(smallSol.getRupSet().getNumSections() == largeSol.getRupSet().getNumSections());
		
		int[] commonParents = {
				FaultSectionUtils.findParentSectionID(smallSol.getRupSet().getFaultSectionDataList(), "Hispaniola"),
				FaultSectionUtils.findParentSectionID(smallSol.getRupSet().getFaultSectionDataList(), "Muertos")
		};
		
		List<FaultSection> combSects = new ArrayList<>();
		List<Double> combSectSolSlipRates = new ArrayList<>();
		List<Double> combSectTargetSlipRates = new ArrayList<>();
		Map<Integer, Double> combAverageSlipsByParent = new HashMap<>();
		Map<Integer, Double> combAverageTargetByParent = new HashMap<>();
		Map<Integer, Integer> idToCombMap = new HashMap<>();
		
		List<Double> absValues = new ArrayList<>();
		List<Double> scaledToRates = new ArrayList<>();
		
		for (boolean small : new boolean[] {false,true}) {
			List<? extends FaultSection> sects;
			SolutionSlipRates solSlipsModule;
			SectSlipRates targets;
			double weight;
			if (small) {
				targets = smallSol.getRupSet().getSectSlipRates();
				sects = smallSol.getRupSet().getFaultSectionDataList();
				solSlipsModule = smallSol.requireModule(SolutionSlipRates.class);
				weight = weightSmall/(weightSmall+weightLarge);
			} else {
				targets = largeSol.getRupSet().getSectSlipRates();
				sects = largeSol.getRupSet().getFaultSectionDataList();
				solSlipsModule = largeSol.requireModule(SolutionSlipRates.class);
				weight = weightLarge/(weightSmall+weightLarge);
			}
			
			List<Double> solSlips = new ArrayList<>();
			for (int i=0; i<sects.size(); i++)
				solSlips.add(solSlipsModule.get(i)*1e3);
			Map<Integer, List<Double>> solSlipsByParent = new HashMap<>();
			for (FaultSection sect : sects) {
				int parentID = sect.getParentSectionId();
				if (!solSlipsByParent.containsKey(parentID))
					solSlipsByParent.put(parentID, new ArrayList<>());
				solSlipsByParent.get(parentID).add(solSlipsModule.get(sect.getSectionId())*1e3);
			}
			Map<Integer, List<Double>> targetSlipsByParent = new HashMap<>();
			for (int i=0; i<sects.size(); i++) {
				double slip = targets.getSlipRate(i)*1e3;
				int parentID = sects.get(i).getParentSectionId();
				if (!targetSlipsByParent.containsKey(parentID))
					targetSlipsByParent.put(parentID, new ArrayList<>());
				targetSlipsByParent.get(parentID).add(slip);
			}
			for (int s=0; s<sects.size(); s++) {
				FaultSection sect = sects.get(s);
				int origID = sect.getSectionId();
				double solSlip = solSlipsModule.get(s)*1e3;
				double targetSlip = targets.getSlipRate(s)*1e3;
				if (Ints.contains(commonParents, sect.getParentSectionId())) {
					Integer index = idToCombMap.get(origID);
					if (index == null) {
						// first time
						sect = sect.clone();
						index = combSects.size();
						sect.setSectionId(index);
						combSects.add(sect);
						combSectSolSlipRates.add(0d);
						combSectTargetSlipRates.add(0d);
						idToCombMap.put(origID, index);
					}
					combSectSolSlipRates.set(index, combSectSolSlipRates.get(index) + weight*solSlip);
					combSectTargetSlipRates.set(index, combSectTargetSlipRates.get(index) + weight*targetSlip);
				} else {
					sect = sect.clone();
					int index = combSects.size();
					idToCombMap.put(sect.getSectionId(), index);
					sect.setSectionId(index);
					combSects.add(sect);
					combSectSolSlipRates.add(solSlip);
					combSectTargetSlipRates.add(targetSlip);
				}
				double absDiff = Math.abs(solSlip - targetSlip);
				absValues.add(absDiff);
				scaledToRates.add(absDiff/targetSlip);
			}
			Map<Integer, Double> solAverageSlipsByParent = solSlipsByParent.entrySet().stream()
					.collect(Collectors.toMap(Map.Entry::getKey,
							entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).average().orElse(0.0)
							));
			Map<Integer, Double> targetAverageSlipsByParent = targetSlipsByParent.entrySet().stream()
					.collect(Collectors.toMap(Map.Entry::getKey,
							entry -> entry.getValue().stream().mapToDouble(Double::doubleValue).average().orElse(0.0)
							));
			for (int parent : solAverageSlipsByParent.keySet()) {
				double parentSolSlip = solAverageSlipsByParent.get(parent);
				double parentTargetSlip = targetAverageSlipsByParent.get(parent);
				if (Ints.contains(commonParents, parent)) {
					double prev = combAverageSlipsByParent.containsKey(parent) ? combAverageSlipsByParent.get(parent) : 0d;
					combAverageSlipsByParent.put(parent, prev + weight*parentSolSlip);
					combAverageTargetByParent.put(parent, prev + weight*parentTargetSlip);
				} else {
					// unique fake parent ID
					int fakeParent = -solAverageSlipsByParent.size();
					combAverageSlipsByParent.put(fakeParent, parentSolSlip);
					combAverageTargetByParent.put(fakeParent, parentTargetSlip);
				}
			}
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(PRVI_SubductionSubSectPlots.plotReg, combSects);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setFillSurfaces(true);
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
		mapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		
		CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1d, 1d);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-20d, 20d);
		
		List<Double> sortables = new ArrayList<>(combSectSolSlipRates.size());
		for (int i=0; i<combSectSolSlipRates.size(); i++)
			sortables.add((double)i);
		
		mapMaker.plotSectScalars(difference(combSectTargetSlipRates, combSectSolSlipRates),
				sortables, diffCPT, "Solution - target slip deficit rates (mm/yr)");
		mapMaker.plot(outputDir, "sub_slip_diff", " ");
		
		mapMaker.plotSectScalars(percentDifference(combSectTargetSlipRates, combSectSolSlipRates),
				sortables, pDiffCPT, "Solution vs target slip deficit rates (% difference)");
		mapMaker.plot(outputDir, "sub_slip_pDiff", " ");
		
		mapMaker.addAnnotations(PRVI_SubductionSubSectPlots.getLabelAnns(mapMaker));
		
		mapMaker.plotSectScalars(difference(combSectTargetSlipRates, combSectSolSlipRates),
				sortables, diffCPT, "Solution - target slip deficit rates (mm/yr)");
		mapMaker.plot(outputDir, "sub_slip_diff_labeled", " ");
		
		mapMaker.plotSectScalars(percentDifference(combSectTargetSlipRates, combSectSolSlipRates),
				sortables, pDiffCPT, "Solution vs target slip deficit rates (% difference)");
		mapMaker.plot(outputDir, "sub_slip_pDiff_labeled", " ");
		
		Range linearRange = new Range(0d, 4d);
		Range logRange = new Range(1e-1, 1e1);
		
		plotScatter(outputDir, "sub_slip_scatter", combSectTargetSlipRates, combSectSolSlipRates,
				combAverageTargetByParent, combAverageSlipsByParent, linearRange, logRange);
		
		CSVFile<String> slipsCSV = new CSVFile<>(true);
		
		slipsCSV.addLine("", "Absolute slip rate misfit (mm/yr)", "Misfit percentage of rate");
		DecimalFormat pDF = new DecimalFormat("0.#%");
		DecimalFormat twoDigits = new DecimalFormat("0.00");
		
		FileWriter slipTEX = new FileWriter(new File(outputDir, "sub_slip_misfits.tex"));
		
		for (boolean max : new boolean[] {false,true}) {
			String label = "Subsection";
			if (max)
				label += ", maximum";
			else
				label += ", average";
			List<String> line = new ArrayList<>();
			line.add(label);
			double abs, scaledToRate;
			String texPrefix = "SubductionSubsect";
			if (max) {
				texPrefix += "Max";
				abs = absValues.stream().mapToDouble(D->D).max().getAsDouble();
				scaledToRate = scaledToRates.stream().mapToDouble(D->D).max().getAsDouble();
			} else {
				texPrefix += "Avg";
				abs = absValues.stream().mapToDouble(D->D).average().getAsDouble();
				scaledToRate = scaledToRates.stream().mapToDouble(D->D).average().getAsDouble();
			}
			texPrefix += "SlipMsft";
			slipTEX.write(LaTeXUtils.defineValueCommand(texPrefix, twoDigits.format(abs))+"\n");
			slipTEX.write(LaTeXUtils.defineValueCommand(texPrefix+"Pct", pDF.format(scaledToRate))+"\n");
			line.add(twoDigits.format(abs));
			line.add(pDF.format(scaledToRate));
			slipsCSV.addLine(line);
		}
		
		slipTEX.close();
		slipsCSV.writeToFile(new File(outputDir, "sub_slip_misfits.csv"));
	}
	
	private static List<XYAnnotation> circleFaultAnns(List<? extends FaultSection> sects, String name,
			List<Double> solSlips, List<Double> targetSlips, TextAnchor anchor) {
		List<XYAnnotation> anns = new ArrayList<>();
		
		MinMaxAveTracker solTrack = new MinMaxAveTracker();
		MinMaxAveTracker targetTrack = new MinMaxAveTracker();
		for (FaultSection sect : sects) {
			if (sect.getSectionName().contains(name)) {
				solTrack.addValue(solSlips.get(sect.getSectionId()));
				targetTrack.addValue(targetSlips.get(sect.getSectionId()));
			}
		}
		Preconditions.checkState(solTrack.getNum() > 0);
		double w = 1 + targetTrack.getLength();
		double h = 1 + solTrack.getLength();
		double centerX = targetTrack.getCenter();
		double centerY = solTrack.getCenter();
		double lowerX = centerX - 0.5*w;
		double lowerY = centerY - 0.5*h;
		Color color = Color.DARK_GRAY;
		XYShapeAnnotation oval = new XYShapeAnnotation(new Ellipse2D.Double(lowerX, lowerY, w, h), new BasicStroke(3f), color);
		anns.add(oval);
		
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		double labelX, labelY;
		if (anchor.isBaseline())
			// baseline means the anchor is on the bottom which means we're above it
			labelY = lowerY+h;
		else if (anchor.isVerticalCenter())
			labelY = centerY;
		else
			labelY = lowerY;
		if (anchor.isLeft())
			// anchor is left which means we're on the right
			labelX = lowerX+w;
		else if (anchor.isHorizontalCenter())
			labelX = centerX;
		else
			labelX = lowerX;
		XYTextAnnotation label = new XYTextAnnotation(name, labelX, labelY);
		label.setFont(annFont);
//		label.setPaint(color);
		label.setTextAnchor(anchor);
		anns.add(label);
		
		return anns;
	}
	
	private static FaultSystemSolution buildNoClassic(File nodeSolDir) throws IOException {
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.add(PRVI25_LogicTree.SEG);
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
	
	private static List<Double> getDMSlipRates(List<? extends FaultSection> subSects) {
		List<Double> slipRates = new ArrayList<>();
		
		for (FaultSection sect : subSects) {
			double slip = sect.getReducedAveSlipRate();
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
			Map<Integer, Double> targetParentSlips, Map<Integer, Double> solParentSlips,
			Range linearRange, Range logRange) throws IOException {
		plotScatter(outputDir, prefix, targetSlips, solSlips, targetParentSlips, solParentSlips, linearRange, logRange, null, null);
	}
	
	private static void plotScatter(File outputDir, String prefix, List<Double> targetSlips, List<Double> solSlips,
			Map<Integer, Double> targetParentSlips, Map<Integer, Double> solParentSlips,
			Range linearRange, Range logRange,
			List<? extends XYAnnotation> linearAnns, List<? extends XYAnnotation> logAnns) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
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
			scatter.setName("Subsections");
			
			funcs.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 4f, Color.BLACK));
		}
		
		if (solParentSlips != null) {
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			for (int id : solParentSlips.keySet())
				scatter.set(targetParentSlips.get(id), solParentSlips.get(id));
			scatter.setName("Fault section averages");
			
			funcs.add(scatter);
			Color color = Colors.tab_blue;
//			color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 180);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, color));
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Target slip rate (mm/yr)", "Solution slip rate (mm/yr)");
		if (solSlips != null && solParentSlips != null)
			plot.setLegendInset(RectangleAnchor.TOP_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		plot.setPlotAnnotations(linearAnns);
		gp.drawGraphPanel(plot, false, false, linearRange, linearRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, false, true, true, false);
		
		plot.setPlotAnnotations(logAnns);
		gp.drawGraphPanel(plot, true, true, logRange, logRange);
		
		PlotUtils.writePlots(outputDir, prefix+"_log", gp, 800, false, true, true, false);
	}

}
