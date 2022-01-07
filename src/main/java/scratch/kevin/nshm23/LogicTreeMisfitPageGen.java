package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

public class LogicTreeMisfitPageGen {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");

//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
		File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//		File mainDir = new File(invDir, "");
		File resultsFile = new File(mainDir, "results.zip");
		
		File outputDir = new File(mainDir, "logic_tree_misfits");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
		Map<LogicTreeLevel<?>, HashSet<LogicTreeNode>> levelNodes = new HashMap<>();
		List<? extends LogicTreeLevel<?>> levels = tree.getLevels();
		for (LogicTreeLevel<?> level : levels)
			levelNodes.put(level, new HashSet<>());
		
		List<LogicTreeBranch<?>> branches = new ArrayList<>();
		List<InversionMisfitStats> branchStats = new ArrayList<>();
		
		Quantity[] quantities = { Quantity.STD_DEV, Quantity.MEAN, Quantity.MAD, Quantity.RMSE};
		
		Quantity[] summaryQuantities = { Quantity.STD_DEV };
		
		int numSummaryBranches = 5;
		
		ZipFile	zip = new ZipFile(resultsFile);
		
		AveragingAccumulator<InversionMisfitStats> fullAccumulator = null;
		for (LogicTreeBranch<?> branch : tree) {
			branches.add(branch);
			
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(InversionMisfitStats.MISFIT_STATS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += InversionMisfitStats.MISFIT_STATS_FILE_NAME;
			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			InversionMisfitStats stats = new InversionMisfitStats(null);
			stats.initFromCSV(csv);
			branchStats.add(stats);
			
			if (fullAccumulator == null)
				fullAccumulator = stats.averagingAccumulator();
			fullAccumulator.process(stats, branch.getBranchWeight());
		}
		
		zip.close();
		
		InversionMisfitStats fullAvgStats = fullAccumulator.getAverage();
		
		List<String> constraintNames = new ArrayList<>();
		Map<String, List<MisfitStats>> constraintStatsMap = new HashMap<>();
		
		for (InversionMisfitStats bstats : branchStats) {
			for (MisfitStats stats : bstats.getStats()) {
				List<MisfitStats> list = constraintStatsMap.get(stats.range.name);
				if (list == null) {
					list = new ArrayList<>();
					constraintStatsMap.put(stats.range.name, list);
					constraintNames.add(stats.range.name);
				}
				list.add(stats);
			}
		}
		
		Map<String, HistogramFunction[]> constraintFullHistsMap = new HashMap<>();
		for (String constraintName : constraintStatsMap.keySet()) {
			List<MisfitStats> constraintStats = constraintStatsMap.get(constraintName);
			HistogramFunction[] fullDistHists = new HistogramFunction[quantities.length];
			for (int i=0; i<quantities.length; i++) {
				double min = Double.POSITIVE_INFINITY;
				double max = Double.NEGATIVE_INFINITY;
				List<Double> values = new ArrayList<>();
				for (MisfitStats stats : constraintStats) {
					double val = stats.get(quantities[i]);
					if (Double.isFinite(val)) {
						values.add(val);
						min = Math.min(min, val);
						max = Math.max(max, val);
					}
				}
				if (!values.isEmpty()) {
					if (max <= min)
						max = 1d;
					double span = max - min;
					Preconditions.checkState(span > 0d);
					double histDelta = Math.max(1e-6, Math.pow(10, Math.floor(Math.log10(span))-1)/2);
					if ((float)max < (float)(min+histDelta))
						max = min+histDelta*1.1;
					fullDistHists[i] = HistogramFunction.getEncompassingHistogram(min, max, histDelta);
					for (double val : values)
						fullDistHists[i].add(fullDistHists[i].getClosestXIndex(val), 1d);
				}
			}
			constraintFullHistsMap.put(constraintName, fullDistHists);
		}
		
		Table<String, LogicTreeNode, List<MisfitStats>> constraintChoiceAllVals = HashBasedTable.create();
		Table<String, LogicTreeBranch<?>, MisfitStats> constraintBranchVals = HashBasedTable.create();
		Map<LogicTreeNode, AveragingAccumulator<InversionMisfitStats>> choiceAccumulators = new HashMap<>();
		for (int i=0; i<branches.size(); i++) {
			LogicTreeBranch<?> branch = branches.get(i);
			InversionMisfitStats bstats = branchStats.get(i);
			
			for (MisfitStats stats : bstats.getStats())
				constraintBranchVals.put(stats.range.name, branch, stats);
			
			for (LogicTreeNode node : branch) {
				AveragingAccumulator<InversionMisfitStats> accumulator = choiceAccumulators.get(node);
				if (accumulator == null) {
					accumulator = bstats.averagingAccumulator();
					choiceAccumulators.put(node, accumulator);
				}
				accumulator.process(bstats, branch.getBranchWeight());
				for (MisfitStats stats : bstats.getStats()) {
					List<MisfitStats> statList = constraintChoiceAllVals.get(stats.range.name, node);
					if (statList == null) {
						statList = new ArrayList<>();
						constraintChoiceAllVals.put(stats.range.name, node, statList);
					}
					statList.add(stats);
				}
			}
		}
		Table<String, LogicTreeNode, MisfitStats> constraintChoiceAvgVals = HashBasedTable.create();
		for (LogicTreeNode node : choiceAccumulators.keySet()) {
			InversionMisfitStats astats = choiceAccumulators.get(node).getAverage();
			for (MisfitStats stats : astats.getStats())
				constraintChoiceAvgVals.put(stats.range.name, node, stats);
		}
		
		List<String> lines = new ArrayList<>();
		
		Map<LogicTreeNode, LogicTreeLevel<?>> nodeLevels = new HashMap<>();
		
		for (LogicTreeBranch<?> branch : tree) {
			for (int i=0; i<levels.size(); i++) {
				LogicTreeNode node = branch.getValue(i);
				levelNodes.get(levels.get(i)).add(node);
				if (!nodeLevels.containsKey(node))
					nodeLevels.put(node, levels.get(i));
			}
		}
		
		CPT colorCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
		
		lines.add("# Logic Tree Misfits");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		
		if (summaryQuantities != null && summaryQuantities.length > 0) {
			for (Quantity quantity : summaryQuantities) {
				lines.add("## "+quantity+" Summary");
				lines.add(topLink); lines.add("");
				
				for (String constraintName : constraintNames) {
					lines.add("### "+constraintName+", "+quantity+" Summary");
					lines.add(topLink); lines.add("");
					
					if (numSummaryBranches > 0) {
						List<LogicTreeBranch<?>> branchesWith = new ArrayList<>();
						List<Double> branchVals = new ArrayList<>();
						for (LogicTreeBranch<?> branch : branches) {
							MisfitStats stats = constraintBranchVals.get(constraintName, branch);
							if (stats != null) {
								double val = Math.abs(stats.get(quantity));
								branchesWith.add(branch);
								branchVals.add(val);
							}
						}
						
						List<ComparablePairing<Double, LogicTreeBranch<?>>> pairs = ComparablePairing.build(branchVals, branchesWith);
						Collections.sort(pairs);
						
						TableBuilder table = MarkdownUtils.tableBuilder();
						
						table.initNewLine();
						table.addColumn("Rank").addColumn("Value");
						for (LogicTreeLevel<?> level : levels)
							if (levelNodes.get(level).size() > 1)
								table.addColumn(level.getName());
						table.finalizeLine();
						
						for (int i=0; i<pairs.size(); i++) {
							if (i < numSummaryBranches || i >= pairs.size()-numSummaryBranches) {
								ComparablePairing<Double, LogicTreeBranch<?>> pair = pairs.get(i);
								table.initNewLine();
								table.addColumn(i+1).addColumn(pair.getComparable().floatValue());
								LogicTreeBranch<?> branch = pair.getData();
								for (int j=0; j<levels.size(); j++)
									if (levelNodes.get(levels.get(j)).size() > 1)
										table.addColumn(branch.getValue(j).getShortName());
								table.finalizeLine();
							} else if (i == numSummaryBranches) {
								table.initNewLine();
								table.addColumn("...").addColumn("...");
								for (LogicTreeLevel<?> level : levels)
									if (levelNodes.get(level).size() > 1)
										table.addColumn("...");
								table.finalizeLine();
							}
						}
						lines.add("**Sorted branch values:**");
						lines.add("");
						lines.addAll(table.build());
						lines.add("");
					}
					
					TableBuilder table = MarkdownUtils.tableBuilder();
					
					table.addLine("Choice", "Mean", "Min", "Max");
					
					for (LogicTreeLevel<?> level : levels) {
						HashSet<LogicTreeNode> nodes = levelNodes.get(level);
						if (nodes.size() < 2)
							continue;
						
						List<LogicTreeNode> sortedNodes = new ArrayList<>(nodes);
						sortedNodes.sort(new Comparator<LogicTreeNode>() {

							@Override
							public int compare(LogicTreeNode o1, LogicTreeNode o2) {
								return o1.getName().compareTo(o2.getName());
							}
						});
						table.addLine("**"+level.getName()+"**", "", "", "");
						for (LogicTreeNode node : sortedNodes) {
							MisfitStats stats = constraintChoiceAvgVals.get(constraintName, node);
							table.initNewLine();
							table.addColumn(node.getName());
							if (stats == null) {
								table.addColumn("_(N/A)_").addColumn("_(N/A)_").addColumn("_(N/A)_");
							} else {
								double min = Double.POSITIVE_INFINITY;
								double max = Double.NEGATIVE_INFINITY;
								for (MisfitStats oStats : constraintChoiceAllVals.get(constraintName, node)) {
									double val = oStats.get(quantity);
									if (Double.isFinite(val)) {
										min = Math.min(val, min);
										max = Math.max(val, max);
									}
								}
								table.addColumn((float)stats.get(quantity)).addColumn((float)min).addColumn((float)max);
							}
							table.finalizeLine();
						}
					}
					lines.addAll(table.build());
					lines.add("");
				}
			}
		}
		
		for (LogicTreeLevel<?> level : levels) {
			HashSet<LogicTreeNode> nodes = levelNodes.get(level);
			if (nodes.size() < 2)
				continue;
			
			List<LogicTreeNode> sortedNodes = new ArrayList<>(nodes);
			sortedNodes.sort(new Comparator<LogicTreeNode>() {

				@Override
				public int compare(LogicTreeNode o1, LogicTreeNode o2) {
					return o1.getName().compareTo(o2.getName());
				}
			});
			
			CPT levelCPT = colorCPT.rescale(0d, sortedNodes.size()-1d);
			
			lines.add("## "+level.getName());
			lines.add(topLink); lines.add("");
			
			for (String constraintName : constraintNames) {
				lines.add("### "+level.getName()+", "+constraintName+" Misfits");
				lines.add(topLink); lines.add("");
				
				for (int q=0; q<quantities.length; q++) {
					Quantity quantity = quantities[q];
					HistogramFunction fullHist = constraintFullHistsMap.get(constraintName)[q];
					if (fullHist == null)
						continue;
					
					lines.add("**"+quantity+"**");
					lines.add("");
					TableBuilder table = MarkdownUtils.tableBuilder();
					
					for (int i=0; i<sortedNodes.size(); i++) {
						LogicTreeNode node = sortedNodes.get(i);
						
						MisfitStats avgMisfit = constraintChoiceAvgVals.get(constraintName, node);
						if (avgMisfit == null)
							continue;
						
						MisfitStats fullAvgMisfit = null;
						for (MisfitStats misfits : fullAvgStats.getStats()) {
							if (misfits.range.name.equals(constraintName)) {
								fullAvgMisfit = misfits;
								break;
							}
						}
						Preconditions.checkNotNull(fullAvgMisfit);
						
						List<MisfitStats> allMisfits = constraintChoiceAllVals.get(constraintName, node);
						
						table.initNewLine();
						table.addColumn(node.getName());
						
						double min = Double.POSITIVE_INFINITY;
						double max = Double.NEGATIVE_INFINITY;
						HistogramFunction subset = new HistogramFunction(fullHist.getMinX(), fullHist.getMaxX(), fullHist.size());
						for (MisfitStats mstats : allMisfits) {
							double val = mstats.get(quantity);
							if (Double.isFinite(val)) {
								min = Math.min(min, val);
								max = Math.max(max, val);
								subset.add(subset.getClosestXIndex(val), 1d);
							}
						}
						
						table.addColumn("Average "+quantity+": "+(float)avgMisfit.get(quantity));
						table.addColumn(quantity+" Range: ["+(float)min+", "+(float)max+"]");
						
						// build plot
						List<XY_DataSet> funcs = new ArrayList<>();
						List<PlotCurveCharacterstics> chars = new ArrayList<>();
						
						Range yRange = new Range(0, 1.1*fullHist.getMaxY());
						Range xRagne = new Range(fullHist.getMinX()-0.5*fullHist.getDelta(), fullHist.getMaxX()+0.5*fullHist.getDelta());
						
						funcs.add(fullHist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
						
						DefaultXY_DataSet fullMean = new DefaultXY_DataSet();
						fullMean.set(fullAvgMisfit.get(quantity), 0d);
						fullMean.set(fullAvgMisfit.get(quantity), yRange.getUpperBound());
						
						funcs.add(fullMean);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
						
						Color color = levelCPT.getColor((float)i);
						
						funcs.add(subset);
						chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
						
						DefaultXY_DataSet nodeMean = new DefaultXY_DataSet();
						nodeMean.set(avgMisfit.get(quantity), 0d);
						nodeMean.set(avgMisfit.get(quantity), yRange.getUpperBound());
						
						funcs.add(nodeMean);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, color.darker().darker()));
						
						HeadlessGraphPanel gp = PlotUtils.initHeadless();
						
						PlotSpec spec = new PlotSpec(funcs, chars, constraintName+" Misfits", quantity.toString(), "Count");
						
						gp.drawGraphPanel(spec, false, false, xRagne, yRange);
						
						String prefix = level.getShortName()+"_"+node.getFilePrefix()+"_"+avgMisfit.range.shortName+"_"+quantity.name()+"_hist";
						prefix = prefix.replaceAll("\\W+", "_");
						PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 650, true, true, false);
						
						table.addColumn("![histogram]("+resourcesDir.getName()+"/"+prefix+".png)");
						
						table.finalizeLine();
					}
					
					table = table.invert();
					
					lines.addAll(table.build());
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}

}
