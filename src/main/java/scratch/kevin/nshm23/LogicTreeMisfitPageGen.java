package scratch.kevin.nshm23;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
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
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.io.archive.ArchiveInput;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.helpers.FileBackedModule;
import org.opensha.sha.earthquake.faultSysSolution.hazard.LogicTreeCurveAverager;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class LogicTreeMisfitPageGen {
	
	public static void main(String[] args) throws IOException {
		String usage = "USAGE: <results.zip> <output-dir> [<filter-prefix1> ... <filter-prefixN>]";
		usage += "\n\tFor filter prefixes, prepend - for exclusion";
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		File resultsFile;
		SolutionLogicTree slt;
		LogicTree<?> tree;
		File outputDir;
		boolean currentWeights = false;
		if (args.length == 0 && invDir.exists()) {
			System.out.println("Assuming hardcoded. Otherwise, usage is:\n"+usage);
			
//			File mainDir = new File(invDir, "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//			File mainDir = new File(invDir, "2023_05_15-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-NoAdj");
//			File mainDir = new File(invDir, "2023_03_23-nshm23_branches-10000ip-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//			File mainDir = new File(invDir, "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
			File mainDir = new File(invDir, "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
			
			resultsFile = new File(mainDir, "results.zip");
			
			slt = SolutionLogicTree.load(resultsFile);
			tree = slt.getLogicTree();
			
//			outputDir = new File(mainDir, "logic_tree_misfits");
			tree = tree.matchingNone(NSHM23_SegmentationModels.CLASSIC);
			outputDir = new File(mainDir, "logic_tree_misfits_no_classic");
		} else if (args.length < 2) {
			System.err.println(usage);
			System.exit(1);
			throw new IllegalArgumentException(usage); // won't actually get here
		} else {
			resultsFile = new File(args[0]);
			
			slt = SolutionLogicTree.load(resultsFile);
			tree = slt.getLogicTree();
			
			outputDir = new File(args[1]);
			
			if (args.length > 2) {
				// filtered
				for (int i=2; i<args.length; i++) {
					String prefix = args[i];
					boolean exclude = false;
					if (prefix.startsWith("-")) {
						exclude = true;
						prefix = prefix.substring(1);
					}
					System.out.println("Filtering tree for prefix '"+prefix+"', exclusion="+exclude);
					LogicTreeNode matchingNode = null;
					for (LogicTreeBranch<?> branch : tree) {
						for (LogicTreeNode node : branch) {
							if (node.getFilePrefix().equals(prefix)) {
								if (matchingNode == null) {
									// first occurrence
									matchingNode = node;
								} else {
									// multiple occurrences
									Preconditions.checkState(matchingNode.equals(node),
											"Can't filter with prefix '%s' as there are multiple matching nodes:"
											+ "\n\tNode 1:\tname='%s'\tshortName=\t'%s'"
											+ "\n\tNode 2:\tname='%s'\tshortName=\t'%s'",
											prefix, matchingNode.getName(), matchingNode.getShortName(),
											node.getName(), node.getShortName());
								}
							}
						}
					}
					Preconditions.checkNotNull(matchingNode, "No branches found with a node matching '%s'", prefix);
					LogicTree<?> filtered;
					if (exclude)
						filtered = tree.matchingNone(matchingNode);
					else
						filtered = tree.matchingAll(matchingNode);
					System.out.println("Filtered from "+tree.size()+" to "+filtered.size()+" branches");
					tree = filtered;
				}
			}
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		if (currentWeights)
			tree.setWeightProvider(new BranchWeightProvider.CurrentWeights());
		Map<LogicTreeLevel<?>, HashSet<LogicTreeNode>> levelNodes = new HashMap<>();
		List<? extends LogicTreeLevel<?>> levels = tree.getLevels();
		for (LogicTreeLevel<?> level : levels)
			levelNodes.put(level, new HashSet<>());
		
		List<LogicTreeBranch<?>> branches = new ArrayList<>();
		List<InversionMisfitStats> branchStats = new ArrayList<>();
		
		Quantity[] quantities = { Quantity.STD_DEV, Quantity.MEAN, Quantity.MAD, Quantity.RMSE};
		
		Quantity[] summaryQuantities = { Quantity.MAD, Quantity.STD_DEV };
		
		int numSummaryBranches = 5;
		
		Map<LogicTreeBranch<?>, InversionMisfitStats> branchMisfits = loadBranchMisfits(resultsFile, tree);
		System.out.println("Loaded misfit stats for "+branchMisfits.size()+" branches");
		
		AveragingAccumulator<InversionMisfitStats> fullAccumulator = null;
		for (LogicTreeBranch<?> branch : tree) {
			branches.add(branch);
			
			InversionMisfitStats stats = branchMisfits.get(branch);
			Preconditions.checkNotNull(stats, "Stats is null? branch: %s", branch);
			branchStats.add(stats);
			
			if (fullAccumulator == null)
				fullAccumulator = stats.averagingAccumulator();
			fullAccumulator.process(stats, tree.getBranchWeight(branch));
		}
		
		InversionMisfitStats fullAvgStats = fullAccumulator.getAverage();
		
		List<String> constraintNames = new ArrayList<>();
		Map<String, List<MisfitStats>> constraintStatsMap = new HashMap<>();
		Map<String, String> constraintShortNames = new HashMap<>();
		
		int numUncertConstraints = 0;
		
		for (InversionMisfitStats bstats : branchStats) {
			for (MisfitStats stats : bstats.getStats()) {
				List<MisfitStats> list = constraintStatsMap.get(stats.range.name);
				if (list == null) {
					list = new ArrayList<>();
					constraintStatsMap.put(stats.range.name, list);
					constraintNames.add(stats.range.name);
					if (stats.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
						numUncertConstraints++;
					constraintShortNames.put(stats.range.name, stats.range.shortName);
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
				accumulator.process(bstats, tree.getBranchWeight(branch));
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
		Map<LogicTreeNode, InversionMisfitStats> choiceMisfitStats = new HashMap<>();
		for (LogicTreeNode node : choiceAccumulators.keySet()) {
			InversionMisfitStats astats = choiceAccumulators.get(node).getAverage();
			choiceMisfitStats.put(node, astats);
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
		
		levels = new ArrayList<>(levels);
		// see if we should filter any out
		for (int i=levels.size(); --i>=0;) {
			LogicTreeLevel<?> level = levels.get(i);
			int numNodes = levelNodes.get(level).size();
			if (LogicTreeCurveAverager.shouldSkipLevel(level, numNodes)) {
				System.out.println("Will skip level: "+level);
				levels.remove(i);
				levelNodes.remove(level);
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
				
				if (numUncertConstraints > 1) {
					List<Double> averages = new ArrayList<>();
					List<Double> averageWeights = new ArrayList<>();
					List<LogicTreeBranch<?>> branchesWithAvg = new ArrayList<>();
					WeightedTrack avgTrack = new WeightedTrack();
					WeightedTrack avgRatioTrack = new WeightedTrack();

					Map<String, WeightedTrack> constraintVals = new HashMap<>();
					Map<String, WeightedTrack> constraintRatios = new HashMap<>();
					
					for (int i=0; i < branches.size(); i++) {
						LogicTreeBranch<?> branch = branches.get(i);
						InversionMisfitStats bStats = branchStats.get(i);
						int numUncert = 0;
						double mySum = 0d;
						double myMax = 0d;
						double myMin = Double.POSITIVE_INFINITY;
						double weight = tree.getBranchWeight(branch);
						for (MisfitStats stats : bStats.getStats()) {
							if (stats.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
								continue;
							
							double val = stats.get(quantity);
							mySum += val;
							myMax = Math.max(myMax, val);
							myMin = Math.min(myMin, val);
							numUncert++;
						}
						if (numUncert < 2)
							continue;
						double myAvg = mySum/(double)numUncert;
						for (MisfitStats stats : bStats.getStats()) {
							if (stats.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
								continue;
							WeightedTrack cVals = constraintVals.get(stats.range.name);
							if (cVals == null) {
								cVals = new WeightedTrack();
								constraintVals.put(stats.range.name, cVals);
								constraintRatios.put(stats.range.name, new WeightedTrack());
							}
							double val = stats.get(quantity);
							double ratio = val/myAvg;
//							avgRatioTrack.add(Math.max(ratio, 1d/ratio), weight);
							avgRatioTrack.add(ratio, weight);
							cVals.add(val, weight);
							constraintRatios.get(stats.range.name).add(ratio, weight);
						}
						averages.add(myAvg);
						averageWeights.add(weight);
						branchesWithAvg.add(branch);
						avgTrack.add(myAvg, weight);
					}
					System.out.println("Building average uncert-weighted summary for "+averages.size()+" branches");
					if (averages.size() > 1) {
						lines.add("### "+quantity+" Average Uncertainty-Weighted Constraint Fits");
						lines.add(topLink); lines.add("");
						
						TableBuilder table = MarkdownUtils.tableBuilder();
						
						table.addLine("Constraint", "Average "+quantity,
								"Best Fitting", "Worst Fitting", "Ratio to All Uncert-Weighted Constraints");
						table.addLine("All Uncertainty-Weighted Average", "**"+(float)avgTrack.average()+"**",
								"**"+(float)avgTrack.min+"**", "**"+(float)avgTrack.max+"**",
								"**"+(float)avgRatioTrack.average()+"**");
						for (String constraintName : constraintNames) {
							WeightedTrack cVals = constraintVals.get(constraintName);
							if (cVals == null)
								continue;
							WeightedTrack cRatios = constraintRatios.get(constraintName);
							table.addLine(constraintName, (float)cVals.average(),
									(float)cVals.min, (float)cVals.max, (float)cRatios.average());
						}
						
						lines.addAll(table.build());
						
						if (numSummaryBranches > 0) {
							List<ComparablePairing<Double, LogicTreeBranch<?>>> pairs = ComparablePairing.build(averages, branchesWithAvg);
							Collections.sort(pairs);
							
							table = MarkdownUtils.tableBuilder();
							
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
							lines.add("**Sorted branch values:** A list of best and worst-fitting branches, sorted by how "
									+ "well they fit all uncertainthy-weighted constraints on average.");
							lines.add("");
							lines.addAll(table.build());
							lines.add("");
						}
						
						// plot distributions here?
						double min = avgTrack.min;
						if (min > 0)
							min = 0;
						double max = avgTrack.max;
						double span = max - min;
						Preconditions.checkState(span > 0d, "max == min: %s %s", max, min);
						double histDelta = (max - min)/40d;
						double roundOrder = Math.pow(10,  Math.floor(Math.log10(0.5*(max - min))-2));
						histDelta = roundOrder * Math.floor(histDelta/roundOrder);
						if ((float)max < (float)(min+histDelta))
							max = min+histDelta*1.1;
						
						HistogramFunction hist = HistogramFunction.getEncompassingHistogram(min, max, histDelta);
						for (double avg : averages)
							hist.add(hist.getClosestXIndex(avg), 1d);
						
						List<XY_DataSet> funcs = new ArrayList<>();
						List<PlotCurveCharacterstics> chars = new ArrayList<>();
						
						Range yRange = new Range(0, 1.1*hist.getMaxY());
						Range xRange = new Range(hist.getMinX()-0.5*hist.getDelta(), hist.getMaxX()+0.5*hist.getDelta());
						
						hist.setName("Distribution");
						funcs.add(hist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
						
						double avg = avgTrack.average();
						DefaultXY_DataSet fullMean = new DefaultXY_DataSet();
						fullMean.set(avg, 0d);
						fullMean.set(avg, yRange.getUpperBound());
						
						fullMean.setName("Average: "+new DecimalFormat("0.##").format(avg));
						funcs.add(fullMean);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
						
						HeadlessGraphPanel gp = PlotUtils.initHeadless();
						
						PlotSpec spec = new PlotSpec(funcs, chars, " ",
								"Constraint "+quantity, "Count");
						spec.setLegendInset(true);
						
						gp.drawGraphPanel(spec, false, false, xRange, yRange);
						
						String prefix = quantity.name()+"_dist_hist";
						prefix = prefix.replaceAll("\\W+", "_");
						PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 650, true, true, false);
						
						lines.add("**Distribution of branch "+quantity+" values, averaged accross constraints:**");
						lines.add("");
						lines.add("![histogram]("+resourcesDir.getName()+"/"+prefix+".png)");
						lines.add("");
					}
					
					// average values
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
							InversionMisfitStats choiceStats = choiceMisfitStats.get(node);
							table.initNewLine();
							table.addColumn(node.getName());
							if (choiceStats == null) {
								table.addColumn("_(N/A)_").addColumn("_(N/A)_").addColumn("_(N/A)_");
							} else {
								double min = Double.POSITIVE_INFINITY;
								double max = Double.NEGATIVE_INFINITY;
								double avg = 0d;
								int num = 0;
								for (MisfitStats stats : choiceStats.getStats()) {
									if (stats.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
										continue;
									double val = stats.get(quantity);
									min = Math.min(min, val);
									max = Math.max(max, val);
									avg += val;
									num++;
								}
								if (num == 0) {
									table.addColumn("_(N/A)_").addColumn("_(N/A)_").addColumn("_(N/A)_");
								} else {
									avg /= num;
									table.addColumn((float)avg).addColumn((float)min).addColumn((float)max);
								}
							}
							table.finalizeLine();
						}
					}
					lines.add("**Average branch choice values:** A list of average fits across all uncertainthy-weighted "
							+ "constraints for each logic tree branch choice.");
					lines.add("");
					lines.addAll(table.build());
					lines.add("");
				}
				
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
						lines.add("**Sorted branch values:** A list of best and worst-fitting branches, sorted by how "
								+ "well they fit the "+constraintName+" constraint.");
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
					lines.add(constraintName+" misfit "+quantity+" for each logic tree branch choice");
					lines.add("");
					lines.addAll(table.build());
					lines.add("");
				}
			}
		}
		
		lines.add("## Constraint Summaries");
		lines.add(topLink); lines.add("");
		
		// see if we have variable weights
		Map<String, List<Double>> constraintWeights = new HashMap<>();
		
		boolean evenWeighting = true;
		for (String constraintName : constraintNames) {
			Double firstWeight = null;
			List<Double> weights = new ArrayList<>();
			for (MisfitStats stats : constraintBranchVals.row(constraintName).values()) {
				if (firstWeight == null)
					firstWeight = stats.range.weight;
				evenWeighting = evenWeighting && (float)stats.range.weight == firstWeight.floatValue();
				weights.add(stats.range.weight);
			}
			constraintWeights.put(constraintName, weights);
		}
		
		lines.add("### Constraint Weights");
		lines.add(topLink); lines.add("");
		
		TableBuilder wTable = MarkdownUtils.tableBuilder();
		
		if (evenWeighting)
			wTable.addLine("Constraint Name", "Weight");
		else
			wTable.addLine("Constraint Name", "Average Weight", "Range", "Weight Distribution");
		for (String constraintName : constraintNames) {
			MinMaxAveTracker track = new MinMaxAveTracker();
			List<Double> weights = constraintWeights.get(constraintName);
			for (double weight : weights)
				track.addValue(weight);
			
			wTable.initNewLine();
			wTable.addColumn(constraintName);
			wTable.addColumn((float)track.getAverage()+"");
			if (evenWeighting) {
				wTable.finalizeLine();
				continue;
			}
			wTable.addColumn("["+(float)track.getMin()+", "+(float)track.getMax()+"]");
			
			double max = track.getMax();
			double min = track.getMin();
			if ((float)min == (float)max) {
				wTable.addColumn("_(N/A)_");
			} else {
				double span = max - min;
				Preconditions.checkState(span > 0d, "max == min: %s %s", max, min);
				double histDelta = (max - min)/40d;
				double roundOrder = Math.pow(10,  Math.floor(Math.log10(0.5*(max - min))-2));
				histDelta = roundOrder * Math.floor(histDelta/roundOrder);
				if ((float)max < (float)(min+histDelta))
					max = min+histDelta*1.1;
				HistogramFunction hist = HistogramFunction.getEncompassingHistogram(min, max, histDelta);
				for (double val : weights)
					hist.add(hist.getClosestXIndex(val), 1d);
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				Range yRange = new Range(0, 1.1*hist.getMaxY());
				Range xRagne = new Range(hist.getMinX()-0.5*hist.getDelta(), hist.getMaxX()+0.5*hist.getDelta());
				
				hist.setName("Weight Distribution");
				funcs.add(hist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
				
				DefaultXY_DataSet fullMean = new DefaultXY_DataSet();
				fullMean.set(track.getAverage(), 0d);
				fullMean.set(track.getAverage(), yRange.getUpperBound());
				
				fullMean.setName("Average Weight: "+new DecimalFormat("0.##").format(track.getAverage()));
				funcs.add(fullMean);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.DARK_GRAY));
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				String shortName = constraintShortNames.get(constraintName);
				String title = constraintName;
				if (title.startsWith("Uncertain "))
					title = title.substring(10);
				title = title.trim();
				PlotSpec spec = new PlotSpec(funcs, chars, title, "Constraint Weight", "Count");
				spec.setLegendInset(true);
				
				gp.drawGraphPanel(spec, false, false, xRagne, yRange);
				
				String prefix = shortName+"_weight_hist";
				prefix = prefix.replaceAll("\\W+", "_");
				PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 650, true, true, false);
				
				wTable.addColumn("![histogram]("+resourcesDir.getName()+"/"+prefix+".png)");
				wTable.finalizeLine();
			}
		}
		
		lines.addAll(wTable.build());
		
		for (String constraintName : constraintNames) {
			lines.add("### "+constraintName+" Summary");
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			MisfitStats fullAvgMisfit = null;
			for (MisfitStats misfits : fullAvgStats.getStats()) {
				if (misfits.range.name.equals(constraintName)) {
					fullAvgMisfit = misfits;
					break;
				}
			}
			Preconditions.checkNotNull(fullAvgMisfit);
			
			table.initNewLine();
			for (Quantity quantity : quantities)
				table.addColumn(quantity+"");
			table.finalizeLine();
			table.initNewLine();
			for (Quantity quantity : quantities)
				table.addColumn("Average: "+(float)fullAvgMisfit.get(quantity));
			table.finalizeLine();
			for (Quantity quantity : quantities) {
				double min = Double.POSITIVE_INFINITY;
				double max = Double.NEGATIVE_INFINITY;
				for (MisfitStats stats : constraintBranchVals.row(constraintName).values()) {
					double val = stats.get(quantity);
					min = Math.min(min, val);
					max = Math.max(max, val);
				}
				table.addColumn("Range: ["+(float)min+", "+(float)max+"]");
			}
			table.finalizeLine();
			table.initNewLine();
			for (int q=0; q<quantities.length; q++) {
				Quantity quantity = quantities[q];
				HistogramFunction fullHist = constraintFullHistsMap.get(constraintName)[q];
				
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
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				PlotSpec spec = new PlotSpec(funcs, chars, constraintName+" Misfits", quantity.toString(), "Count");
				
				gp.drawGraphPanel(spec, false, false, xRagne, yRange);
				
				String prefix = fullAvgMisfit.range.shortName+"_"+quantity.name()+"_full_hist";
				prefix = prefix.replaceAll("\\W+", "_");
				PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 650, true, true, false);
				
				table.addColumn("![histogram]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine();
			
			if (numUncertConstraints > 1 &&
					constraintBranchVals.row(constraintName).values().iterator().next().range.weightingType
					== ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				table.initNewLine();
				
				for (Quantity quantity : quantities) {
					// scatter plot of how well this constraint was fit relative to others
					DefaultXY_DataSet scatter = new DefaultXY_DataSet();
					
					String shortName = null;
					
					for (int i=0; i<branches.size(); i++) {
						LogicTreeBranch<?> branch = branches.get(i);
						InversionMisfitStats allStats = branchStats.get(i);
						MisfitStats myConstrStats = constraintBranchVals.get(constraintName, branch);
						if (myConstrStats == null || myConstrStats.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
							// this branch doesn't have this constraint
							continue;
						
						if (shortName == null)
							shortName = myConstrStats.range.shortName;
						
						double myVal = myConstrStats.get(quantity);
						double avgOther = 0d;
						int numOther = 0;
						for (MisfitStats oStats : allStats.getStats()) {
							if (oStats.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY
									&& !oStats.range.name.equals(constraintName)) {
								numOther++;
								avgOther += oStats.get(quantity);
							}
						}
						if (numOther == 0)
							continue;
						avgOther /= (double)numOther;
						scatter.set(avgOther, myVal);
					}
					
					if (scatter.size() == 0)
						continue;
					
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					funcs.add(scatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, Color.BLACK));
					
					double max = Math.ceil(Math.max(scatter.getMaxX(), scatter.getMaxY()));
					double min = Math.floor(Math.min(scatter.getMinX(), scatter.getMinY()));
					Preconditions.checkState(max > min, "%s %s", max, min);
					
					DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
					oneToOne.set(min, min);
					oneToOne.set(max, max);
					
					funcs.add(oneToOne);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
					
					PlotSpec spec = new PlotSpec(funcs, chars, shortName+" Relative Fits",
							"Average Other Constraint "+quantity, shortName+" "+quantity);
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					Range range = new Range(min, max);
					gp.drawGraphPanel(spec, false, false, range, range);
					
					String prefix = shortName+"_"+quantity.name()+"_misfit_scatter";
					prefix = prefix.replaceAll("\\W+", "_");
					
					PlotUtils.writePlots(resourcesDir, prefix, gp, 800, -1, true, false, false);
					
					table.addColumn("!["+quantity+" Scatter]("+resourcesDir.getName()+"/"+prefix+".png)");
				}
				table.finalizeLine();
			}
			
			lines.addAll(table.build());
			lines.add("");
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
			
			if (numUncertConstraints > 0) {
				lines.add("### "+level.getName()+" Average Constraint Fits");
				lines.add(topLink); lines.add("");
				
				lines.add("Average misfits across all constraints as a function of "+level.getName()+" choice");
				lines.add("");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				table.addLine("Choice", "Quantity", "Average Constraint Misfit");
				for (int q=0; q<quantities.length; q++) {
					Quantity quantity = quantities[q];
					
					for (int i=0; i<sortedNodes.size(); i++) {
						LogicTreeNode node = sortedNodes.get(i);
						
						double avgVal = 0d;
						int avgNum = 0;
						for (MisfitStats stats : constraintChoiceAvgVals.column(node).values()) {
							if (stats.range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
								avgVal += stats.get(quantity);
								avgNum++;
							}
						}
						if (avgNum == 0)
							continue;
						
						avgVal /= (double)avgNum;
						
						table.addLine(node.getShortName(), quantity+"", (float)avgVal);
					}
				}
				
				lines.addAll(table.build());
				lines.add("");
			}
			
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
	
	public static Map<LogicTreeBranch<?>, InversionMisfitStats> loadBranchMisfits(File resultsFile) throws IOException {
		return loadBranchMisfits(resultsFile, null);
	}
	
	public static Map<LogicTreeBranch<?>, InversionMisfitStats> loadBranchMisfits(File resultsFile, LogicTree<?> tree)
			throws IOException {
		ArchiveInput input = new ArchiveInput.ZipFileInput(resultsFile);
		
		if (tree == null) {
			BufferedInputStream logicTreeIS = FileBackedModule.getInputStream(input, "solution_logic_tree/", "logic_tree.json");
			Gson gson = new GsonBuilder().registerTypeAdapter(LogicTree.class, new LogicTree.Adapter<>()).create();
			InputStreamReader reader = new InputStreamReader(logicTreeIS);
			tree = gson.fromJson(reader, LogicTree.class);
		}
		Map<LogicTreeBranch<?>, InversionMisfitStats> ret = new HashMap<>();
		for (LogicTreeBranch<?> branch : tree) {
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(InversionMisfitStats.MISFIT_STATS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += InversionMisfitStats.MISFIT_STATS_FILE_NAME;
//			System.out.println("Loading "+entryName);
			Preconditions.checkNotNull(input.hasEntry(entryName), "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(input.getInputStream(entryName), true);
			InversionMisfitStats stats = new InversionMisfitStats(null);
			stats.initFromCSV(csv);
			ret.put(branch, stats);
		}
		
		input.close();
		return ret;
	}
	
	private static class WeightedTrack {
		double sumWeighted = 0d;
		double sumWeights = 0d;
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		
		public void add(double val, double weight) {
			sumWeighted += val*weight;
			sumWeights += weight;
			min = Math.min(min, val);
			max = Math.max(max, val);
		}
		
		public double average() {
			return sumWeighted/sumWeights;
		}
	}

}
