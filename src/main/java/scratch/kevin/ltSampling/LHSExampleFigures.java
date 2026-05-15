package scratch.kevin.ltSampling;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.TruncatedNormalDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.ContinuousDistributionSampledLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.SimpleValuedNode;
import org.opensha.commons.logicTree.LogicTreeNode.ValuedLogicTreeNode;
import org.opensha.commons.logicTree.LogicTreePairwiseLHSIteration;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class LHSExampleFigures {

	public static void main(String[] args) throws IOException {
		int samples = 30;
//		int samples = 100;
//		int samples = 1000;
//		int samples = 2000;

		int dpi = 300;
		
		File outputDir = new File("/tmp/lhs_samples");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(), "Can't create output dir: %s", outputDir.getAbsolutePath());
		
		outputDir = new File(outputDir, samples+"_samples");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(), "Can't create output dir: %s", outputDir.getAbsolutePath());
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.add(LogicTreeLevel.forEnum(DefModelEnum.class, "Deformation Model", "Deformation Model"));
		levels.add(LogicTreeLevel.forEnum(ScaleEnum.class, "Scaling Relationship", "Scaling Relationship"));
		levels.add(new ContinuousDistributionSampledLevel(
				"GR b-value", "GR b-value", UniformContinuousDistribution.of(0d, 1d), "Sample ", "Sample", "Sample"));
		levels.add(new ContinuousDistributionSampledLevel(
				"Off-fault Mmax", "Off-fault Mmax", TruncatedNormalDistribution.of(7.6, 0.2, 7.15, 8.05), -1, "Sample ", "Sample", "Sample"));
		
		long seed = 123456789l;
		
		CPT tab10cpt = GMT_CPT_Files.CATEGORICAL_TAB10.instance();
		Color[] tab10 = new Color[tab10cpt.size()];
		for (int i=0; i<tab10.length; i++)
			tab10[i] = tab10cpt.get(i).minColor;
		
		Map<LogicTreeNode, Color> nodeColors = new HashMap<>();
		int colorI = 0;
		for (LogicTreeLevel<? extends LogicTreeNode> level : levels) {
			if (!(level instanceof ContinuousDistributionSampledLevel)) {
				for (LogicTreeNode node : level.getNodes())
					nodeColors.put(node, tab10[colorI++ % tab10.length]);
			}
		}

		Color[] distColors = new Color[levels.size()];
		CPT[] distCPTs = new CPT[levels.size()];
		boolean useDistCPTforPDF = true;
		
		for (int l=0; l<levels.size(); l++) {
			if (levels.get(l) instanceof ContinuousDistributionSampledLevel) {
				distColors[l] = tab10[colorI++];
				distCPTs[l] = new CPT(0d, 1d, distColors[l].brighter().brighter(), distColors[l], distColors[l].darker().darker());
			}
		}
		
		for (SamplingMethod sm : SamplingMethod.values()) {
			LogicTree<LogicTreeNode> tree = LogicTree.buildSampled(levels, samples, seed, sm);
			if (sm == SamplingMethod.MONTE_CARLO) {
				// write tree plot
				LogicTreeFigureWriter tfw = new LogicTreeFigureWriter(tree, false, true);
				tfw.write(outputDir, "logic_tree", true, true);
			}
			
			for (int l=0; l<levels.size(); l++) {
				LogicTreeLevel<? extends LogicTreeNode> level = levels.get(l);
				
				String samplePrefix = level.getFilePrefix()+"_samples_"+sm.name();
				
				if (level instanceof ContinuousDistributionSampledLevel) {
					ContinuousDistributionSampledLevel distLevel = (ContinuousDistributionSampledLevel)level;
					ContinuousDistribution dist = distLevel.getDistribution();
					
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					Range xRange = new Range(dist.getSupportLowerBound(), dist.getSupportUpperBound());
					EvenlyDiscretizedFunc densityFunc = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), 1000);
					for (int i=0; i<densityFunc.size(); i++)
						densityFunc.set(i, dist.density(densityFunc.getX(i)));
					
					double xRangeDelta = 0.003*xRange.getLength();
					xRange = new Range(xRange.getLowerBound()-xRangeDelta, xRange.getUpperBound()+xRangeDelta);
					
					List<? extends SimpleValuedNode<Double>> nodes = distLevel.getNodes();
					
					double maxY;
					if (samples < 200) {
						maxY = Math.max(2d, densityFunc.getMaxY()*1.1);
						
						funcs.add(densityFunc);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.GRAY));
						
						if (sm.isLHS()) {
							// draw ticks
							double tickDelta = maxY*0.05;
							EvenlyDiscretizedFunc binEdges = new EvenlyDiscretizedFunc(
									0d, 1d, nodes.size()+1);
							PlotCurveCharacterstics tickChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 0.7f,	Color.BLACK);
							for (int i=0; i<binEdges.size(); i++) {
								double p = binEdges.getX(i);
								double x = dist.inverseCumulativeProbability(p);
								double yCenter = dist.density(x);
								
								funcs.add(new DefaultXY_DataSet(x, yCenter-tickDelta, x, yCenter+tickDelta));
								chars.add(tickChar);
							}
						}
						
						float symbolWidth = 5f;
						for (SimpleValuedNode<Double> node : nodes) {
							DefaultXY_DataSet xy = new DefaultXY_DataSet();
							double x = node.getValue();
							double y = dist.density(x);
							xy.set(x, y);
							
							Color color = useDistCPTforPDF ? distCPTs[l].getColor(dist.cumulativeProbability(x)) : distColors[l];
							
							funcs.add(xy);
							chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, symbolWidth, color));
							
							funcs.add(xy);
							chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, symbolWidth, distColors[l].darker().darker()));
						}
					} else {
						maxY = Math.max(2d, densityFunc.getMaxY()*1.5);
						double length = dist.getSupportUpperBound() - dist.getSupportLowerBound();
						int bins = 100;
						double binWidth = length / bins;
						EvenlyDiscretizedFunc hist = new EvenlyDiscretizedFunc(dist.getSupportLowerBound()+0.5*binWidth, bins, binWidth);
						for (SimpleValuedNode<Double> node : nodes)
							hist.add(hist.getClosestXIndex(node.getValue()), 1d);
						
						hist.scale(1d/(samples*binWidth));
						
						funcs.add(hist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, distColors[l]));
						
						funcs.add(densityFunc);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
					}
					
					PlotSpec plot = new PlotSpec(funcs, chars, sm.getName(), level.getName(), "PDF Density");
					
					Range yRange = new Range(0d, maxY);
					
					HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
					gp.drawGraphPanel(plot, false, false, xRange, yRange);
					PlotUtils.setYTick(gp, 0.5);
					
					PlotUtils.writePrintPlots(outputDir, samplePrefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*0.5, 3d, dpi, true, true, false);
				} else {
					List<? extends LogicTreeNode> nodes = level.getNodes();
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					List<XYTextAnnotation> anns = new ArrayList<>();
					
					int maxCount = 0;

					Font nameFont = new Font(Font.SANS_SERIF, Font.BOLD, 12);
					Font countFont = new Font(Font.SANS_SERIF, Font.BOLD, 12);
					Font subNameFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
					Font subCountFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
					
					double offset = samples*0.01;
					
					for (int i=0; i<nodes.size(); i++) {
						LogicTreeNode node = nodes.get(i);
						int count = 0;
						for (LogicTreeBranch<LogicTreeNode> branch : tree)
							if (branch.getValue(l).equals(node))
								count++;
						maxCount = Integer.max(count, maxCount);
						EvenlyDiscretizedFunc hist = new EvenlyDiscretizedFunc(0.5, nodes.size(), 1d);
						hist.set(i, (double)count);
						
						funcs.add(hist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, nodeColors.get(node)));
						
						XYTextAnnotation ann = new XYTextAnnotation(node.getShortName(), hist.getX(i), count+offset);
						ann.setFont(nameFont);
						ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
						anns.add(ann);
						
						ann = new XYTextAnnotation(count+"", hist.getX(i), offset);
						ann.setFont(countFont);
						ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
						anns.add(ann);
					}
					
					PlotSpec plot = new PlotSpec(funcs, chars, sm.getName(), level.getName(), "Sample Count");
					plot.setPlotAnnotations(anns);
					
					HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
					double maxY = switch (nodes.size()){
					case 2:
						yield 0.65*samples;
					case 3:
						yield 0.45*samples;
					default:
						yield maxCount*1.2d;
					};
					gp.drawGraphPanel(plot, false, false, new Range(0d, nodes.size()), new Range(0d, maxY));
					PlotUtils.setAxisVisible(gp, false, true);
					PlotUtils.setGridLinesVisible(gp, false, true);
					
					PlotUtils.writePrintPlots(outputDir, samplePrefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*0.5, 3d, dpi, true, true, false);
					
					for (int m=0; m<levels.size(); m++) {
						LogicTreeLevel<? extends LogicTreeNode> oLevel = levels.get(m);
						if (l == m || oLevel instanceof ContinuousDistributionSampledLevel)
							continue;
						
						String pairPrefix = level.getFilePrefix()+"_and_"+oLevel.getFilePrefix()+"_samples_"+sm.name();
						
						funcs.clear();
						chars.clear();
						anns.clear();
						
						List<? extends LogicTreeNode> oNodes = levels.get(m).getNodes();
						
						for (int i=0; i<nodes.size(); i++) {
							LogicTreeNode node = nodes.get(i);
							int count = 0;
							double x = 0.5 + i;
							for (int j=0; j<oNodes.size(); j++) {
								LogicTreeNode oNode = oNodes.get(j);
								int subCount = 0;
								for (LogicTreeBranch<LogicTreeNode> branch : tree)
									if (branch.getValue(l).equals(node) && branch.getValue(m).equals(oNode))
										subCount++;
								
								int countStart = count;
								int countEnd = count+subCount;
								EvenlyDiscretizedFunc hist = new EvenlyDiscretizedFunc(0.5, nodes.size(), 1d);
								hist.set(i, (double)countEnd);
								
								funcs.add(hist);
								chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, nodeColors.get(oNode)));
								
								XYTextAnnotation ann = new XYTextAnnotation(oNode.getShortName(), x, countEnd);
								ann.setFont(subNameFont);
								ann.setTextAnchor(TextAnchor.TOP_CENTER);
								anns.add(ann);
								
								ann = new XYTextAnnotation(subCount+"", x, countStart+0.5*offset);
								ann.setFont(subCountFont);
								ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
								anns.add(ann);
								
								count = countEnd;
							}
							
							XYTextAnnotation ann = new XYTextAnnotation(node.getShortName(), x, count+offset);
							ann.setFont(nameFont);
							ann.setTextAnchor(TextAnchor.BASELINE_CENTER);
							anns.add(ann);
						}
						
						plot = new PlotSpec(funcs, chars, sm.getName(), level.getName(), "Sample Count");
						plot.setPlotAnnotations(anns);
						
						gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
						gp.drawGraphPanel(plot, false, false, new Range(0d, nodes.size()), new Range(0d, maxY));
						PlotUtils.setAxisVisible(gp, false, true);
						PlotUtils.setGridLinesVisible(gp, false, true);
						
						PlotUtils.writePrintPlots(outputDir, pairPrefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*0.5, 3d, dpi, true, true, false);
					}
				}
			}
			
			if (samples <= 30) {
				// build branch vector plot
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				List<XYAnnotation> anns = new ArrayList<>();
				
				BasicStroke outlineStroke = new BasicStroke(1f);
//				Font indexOrigFont = new Font(Font.SANS_SERIF, Font.PLAIN, 10);
//				Font indexSwappedFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
				

				Font indexOrigFont = null;
				Font indexSwappedFont = new Font(Font.SANS_SERIF, Font.PLAIN, 10);
				
				List<int[]> origBranchIndexes = null;
				if (sm ==SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE) {
					LogicTree<LogicTreeNode> tree2 = LogicTree.buildSampled(levels, samples, seed, SamplingMethod.LATIN_HYPERCUBE);
					List<double[]> levelFixedWeights = new ArrayList<>();
					for (int l=0; l<levels.size(); l++) {
						LogicTreeLevel<? extends LogicTreeNode> level = levels.get(l);
						if (level instanceof ContinuousDistributionSampledLevel) {
							levelFixedWeights.add(null);
						} else {
							List<? extends LogicTreeNode> nodes = level.getNodes();
							double[] weights = new double[nodes.size()];
							for (int n=0; n<nodes.size(); n++)
								weights[n] = nodes.get(n).getNodeWeight(null);
							levelFixedWeights.add(weights);
						}
					}
					LogicTreePairwiseLHSIteration<LogicTreeNode> iter = new LogicTreePairwiseLHSIteration<>(
							levels, tree2.getBranches(), levelFixedWeights);
					
					iter.setTrackSwaps(true);
					
					iter.iterate(Integer.max(10000, samples*100), new Random(seed), false);
					
					origBranchIndexes = iter.getOriginalBranchIndexes();
					tree = tree2;
				}
				
				for (int n=0; n<samples; n++) {
					double x = n;
					
					LogicTreeBranch<LogicTreeNode> branch = tree.getBranch(n);
					
					double x0 = n-0.4;
					double x1 = n+0.4;
					
					for (int l=0; l<levels.size(); l++) {
						LogicTreeLevel<? extends LogicTreeNode> level = levels.get(l);
						LogicTreeNode node = branch.getValue(l);
						
						Color color;
						if (level instanceof ContinuousDistributionSampledLevel) {
							ContinuousDistribution dist = ((ContinuousDistributionSampledLevel)level).getDistribution();
							double value = ((ValuedLogicTreeNode<Double>)node).getValue();
							color = distCPTs[l].getColor(dist.cumulativeProbability(value));
						} else {
							color = nodeColors.get(node);
						}
						
						double y1 = (levels.size()-l);
						double y0 = y1-1;
						double y = y0+0.5;
//						double y0 = l;
//						double y1 = l+1;
//						double y = l+0.5;
						
						if (origBranchIndexes == null) {
							anns.add(new XYBoxAnnotation(x0, y0, x1, y1, outlineStroke, Color.BLACK, color));
						} else {
							int index = origBranchIndexes.get(n)[l];
							XYTextAnnotation indexAnn = new XYTextAnnotation(index+"", x, y);
							indexAnn.setTextAnchor(TextAnchor.CENTER);
							if (index == n) {
								if (indexOrigFont != null) {
									indexAnn.setFont(indexOrigFont);
									anns.add(indexAnn);
								}
								anns.add(new XYBoxAnnotation(x0, y0, x1, y1, outlineStroke, Color.BLACK,
										new Color(color.getRed(), color.getGreen(), color.getBlue(), 60)));
							} else {
								anns.add(new XYBoxAnnotation(x0, y0, x1, y1, outlineStroke, Color.BLACK, color));
								indexAnn.setFont(indexSwappedFont);
								anns.add(indexAnn);
							}
						}
						
						
					}
				}
				
				PlotSpec plot = new PlotSpec(funcs, chars, sm.getName(), "Branch Index", " ");
				plot.setPlotAnnotations(anns);
				
				HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
				gp.drawGraphPanel(plot, false, false, new Range(-0.5d, samples-0.5), new Range(-0.05d, levels.size()+0.05));
				PlotUtils.setAxisVisible(gp, true, false);
				PlotUtils.setGridLinesVisible(gp, false, false);
				
				PlotUtils.writePrintPlots(outputDir, "branches_"+sm.name(), gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH, 3d, dpi, true, true, false);
			}
		}
	}
	
	private enum DefModelEnum implements LogicTreeNode.FixedWeightNode {
		GEOLOGIC("Geologic", 0.5),
		GEODETIC("Geodetic", 0.5);
		
		private String name;
		private double weight;

		private DefModelEnum(String name, double weight) {
			this.name = name;
			this.weight = weight;
		}

		@Override
		public String getFilePrefix() {
			return name();
		}

		@Override
		public String getShortName() {
			return name;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public double getNodeWeight() {
			return weight;
		}
	}
	
	private enum ScaleEnum implements LogicTreeNode.FixedWeightNode {
		LOGA_4p1("LogA+4.1", 1d/3d),
		LOGA_4p2("LogA+4.2", 1d/3d),
		LOGA_4p3("LogA+4.3", 1d/3d);
		
		private String name;
		private double weight;

		private ScaleEnum(String name, double weight) {
			this.name = name;
			this.weight = weight;
		}

		@Override
		public String getFilePrefix() {
			return name();
		}

		@Override
		public String getShortName() {
			return name;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public double getNodeWeight() {
			return weight;
		}
	}

}
