package scratch.kevin.sampling;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;
import java.util.function.IntToDoubleFunction;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.sampling.CategoricalSamplingDimension;
import org.opensha.commons.data.sampling.ContinuousSamplingDimension;
import org.opensha.commons.data.sampling.PointSet;
import org.opensha.commons.data.sampling.SamplingDimension;
import org.opensha.commons.data.sampling.scoring.CenteredDiscrepancy;
import org.opensha.commons.data.sampling.scoring.ProjectionDiscrepancyScore;
import org.opensha.commons.data.sampling.scoring.ProjectionDiscrepancyScorer;
import org.opensha.commons.data.sampling.scoring.ProjectionDiscrepancyScore.ProjectionResult;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.sampling.LogicTreePointSetMapper;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.RandomSeedUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_ModelRegimeNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import net.mahdilamb.colormap.Colors;

public class SamplingScoreFigures {
	
	static List<SamplingDimension> getDimsNSHM23() {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>(NSHM23_LogicTreeBranch.levelsCombined);
		// remove the first fault model (fixed) level
		levels.remove(0);
		LogicTreePointSetMapper<LogicTreeNode> mapper = new LogicTreePointSetMapper<>(levels);
		List<SamplingDimension> dims = mapper.getSamplingDimensions();
		System.out.println("NSHM23 levels:");
		for (int l=0; l<levels.size(); l++)
			System.out.println("\t"+l+". "+levels.get(l).getShortName()+":\t"+dims.get(l));
		return dims;
	}
	
	static List<SamplingDimension> getDimsNSHM27_AmSam() {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.addAll(NSHM27_LogicTree.buildLevels(NSHM27_SeismicityRegions.AMSAM, TectonicRegionType.SUBDUCTION_INTERFACE, true, true, true, true));
		levels.addAll(NSHM27_LogicTree.buildLevels(NSHM27_SeismicityRegions.AMSAM, TectonicRegionType.SUBDUCTION_SLAB, true, true, true, false));
		levels.addAll(NSHM27_LogicTree.buildLevels(NSHM27_SeismicityRegions.AMSAM, TectonicRegionType.ACTIVE_SHALLOW, true, true, true, false));
		// remove the model/regime (fixed) level
		for (int l=levels.size(); --l>=0;)
			if (levels.get(l) instanceof NSHM27_ModelRegimeNode.Level)
				levels.remove(l);
		LogicTreePointSetMapper<LogicTreeNode> mapper = new LogicTreePointSetMapper<>(levels);
		List<SamplingDimension> dims = mapper.getSamplingDimensions();
		System.out.println("NSHM27 levels:");
		for (int l=0; l<levels.size(); l++)
			System.out.println("\t"+l+". "+levels.get(l).getShortName()+":\t"+dims.get(l));
		return dims;
	}

	public static void main(String[] args) throws IOException {
		File mainDir = new File(PaperPaths.FIGURES_DIR, "scores");
		Preconditions.checkState(mainDir.exists() || mainDir.mkdir());
		int scoreOrders = 4;
		int[] sampleCounts = {256, 512, 1024, 2048, 4096, 8192, 16384};
//		int[] sampleCounts = {256, 512, 1024, 2048, 4096, 8192};
//		int[] sampleCounts = {256, 512, 1024, 2048, 4096};
//		int[] sampleCounts = {256, 512, 1024};
//		int[] sampleCounts = {256, 512};
		SamplingMethod[] methods = {
				SamplingMethod.MONTE_CARLO,
				SamplingMethod.LATIN_HYPERCUBE,
				SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE,
				SamplingMethod.SOBOL,
				SamplingMethod.OWEN_SCRAMBLED_SOBOL
		};
		Map<SamplingMethod, PlotCurveCharacterstics> combPlotChars = Map.of(
				SamplingMethod.MONTE_CARLO,
						new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Colors.tab_red),
				SamplingMethod.LATIN_HYPERCUBE,
						new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Colors.tab_green),
				SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE,
						new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Colors.tab_orange),
				SamplingMethod.OWEN_SCRAMBLED_SOBOL,
						new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_blue));
//		int numPlotTrials = 10;
		int numPlotTrials = 0;
//		int numAvgTrials = 20;
//		int numAvgTrials = 50;
		int numAvgTrials = 100;
//		int numAvgTrials = 500;
		
		boolean redoNormScores = false;
		boolean redoCenteredDiscrepancies = false;
		boolean replotIndvSamples = false;
		
		String treeName = null;
		List<SamplingDimension> samplingDimensions = new ArrayList<>();
		for (int i=0; i<10; i++)
			samplingDimensions.add(ContinuousSamplingDimension.INSTANCE);
		String samplingPrefix = "continuous_"+samplingDimensions.size()+"d";
		
//		String treeName = "NSHM23-WUS";
//		List<SamplingDimension> samplingDimensions = getDimsNSHM23();
//		String samplingPrefix = "nshm23_"+samplingDimensions.size()+"d";
		
//		String treeName = "NSHM27-AmSam";
//		List<SamplingDimension> samplingDimensions = getDimsNSHM27_AmSam();
//		String samplingPrefix = "nshm27_amsam_"+samplingDimensions.size()+"d";
		
		final int dimensions = samplingDimensions.size();
		int numContinuous = 0;
		int numCategorical = 0;
		for (SamplingDimension dim : samplingDimensions) {
			if (dim instanceof CategoricalSamplingDimension)
				numCategorical++;
			else if (dim instanceof ContinuousSamplingDimension)
				numContinuous++;
		}
		
		if (treeName == null)
			treeName = "";
		else
			treeName += ": ";
		if (dimensions == numContinuous)
			treeName += dimensions+"D, all continuous";
		else if (dimensions == numCategorical)
			treeName += dimensions+"D, all categorical";
		else
			treeName += numContinuous+" continuous, "+numCategorical+" categorical";
		
		System.out.println(treeName);
		System.out.println();
		
		File outputDir = new File(mainDir, samplingPrefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File scoresCSVFile = new File(outputDir, "combined_scores.csv");
		File centeredDiscrepancyCSVFile = new File(outputDir, "combined_centered_discrepancies.csv");
		
		redoNormScores |= !scoresCSVFile.exists();
		redoCenteredDiscrepancies |= !centeredDiscrepancyCSVFile.exists();
		
		ProjectionDiscrepancyScorer serialScorer = ProjectionDiscrepancyScorer.exact(1);
//		ProjectionDiscrepancyScorer serialScorer = ProjectionDiscrepancyScorer.exact(4);
		ProjectionDiscrepancyScorer parallelScorer = ProjectionDiscrepancyScorer.exact(16);
		
		Color[] orderColors = new Color[scoreOrders];
		Color[] oderLightColors = new Color[scoreOrders];
		CPT catCPT = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
		CPT catLightCPT = GMT_CPT_Files.CATEGORICAL_TAB10_LIGHT_NOGRAY.instance();
		for (int d=0; d<scoreOrders; d++) {
			orderColors[d] = catCPT.get(d % catCPT.size()).minColor;
			oderLightColors[d] = catLightCPT.get(d % catLightCPT.size()).minColor;
		}
		
//		IntToDoubleFunction orderThicknessFunc = (order)->1d+(scoreOrders-order)/2d;
		
		// this results in the following, and lower orders don't change as more are added:
		// 1D:	3.375
		// 2D:	2.25
		// 3D:	1.5
		// 4D:	1.0
		// 5D:	0.6666667
		// 6D:	0.44444445
		IntToDoubleFunction orderThicknessFunc = (order)->Math.pow(1.5, 4-order);
//		for (int order=1; order<=scoreOrders; order++)
//			System.out.println(order+"D:\t"+(float)orderThicknessFunc.applyAsDouble(order));
//		System.exit(0);
		
		Range dimXRange = new Range(1d, dimensions);
		Range logYRange = new Range(1e-4, 2e0);
		Range equivYRange = new Range(1e2, sampleCounts[sampleCounts.length-1] > 3000 ? 1e8 : 1e7);
		Range centeredYRange = new Range(1e-6, 1e-1);
		
		List<List<List<ProjectionDiscrepancyScore>>> methodScores = new ArrayList<>();
		List<List<List<Double>>> methodCenteredDiscrepancyScores = new ArrayList<>();
		for (int m=0; m<methods.length; m++) {
			methodScores.add(new ArrayList<>());
			methodCenteredDiscrepancyScores.add(new ArrayList<>());
		}
		
		
		if (!redoNormScores && !redoCenteredDiscrepancies) {
			System.out.println("Replotting combined results only");
		} else {
			Stopwatch totalWatch = Stopwatch.createStarted();
			for (int sampleCount : sampleCounts) {
				Stopwatch sampleWatch = Stopwatch.createStarted();
				System.out.println("Doing "+sampleCount+" samples");
				File subDir = new File(outputDir, sampleCount+"_samples");
				Preconditions.checkState(subDir.exists() || subDir.mkdir());
				PointSet[] firstPointSets = new PointSet[methods.length];
				for (int m=0; m<methods.length; m++) {
					Stopwatch methodWatch = Stopwatch.createStarted();				
					SamplingMethod method = methods[m];
					String prefix = method.name();
					
					// use the same seed for pairwise and regular LHS
					String seedName = method == SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE ? SamplingMethod.LATIN_HYPERCUBE.name() : method.name();
					Random baseRand = new Random(RandomSeedUtils.seedForStrings(seedName));
					
					int myTrials = method == SamplingMethod.SOBOL ? 1 : numAvgTrials;
					System.out.println("Doing "+method+", "+sampleCount+" samples x "+myTrials+" trials");
					
					LinkedList<CompletableFuture<PointSet>> sampleFutures = new LinkedList<>();
					
					for (int i=0; i<myTrials; i++) {
						long seed = baseRand.nextLong();
						sampleFutures.add(CompletableFuture.supplyAsync(()->method.prepare(sampleCount, samplingDimensions, seed)));
					}
					
					List<PointSet> samples = new ArrayList<>(myTrials);
					
					while (!sampleFutures.isEmpty())
						samples.add(sampleFutures.removeFirst().join());
					
					if (redoNormScores) {
						// if we only have 1 trial, do that one in parallel
						// if we have many, rely on across-trial parallelism instead
						ProjectionDiscrepancyScorer scorer = myTrials == 1 ? parallelScorer : serialScorer;
						
						List<CompletableFuture<ProjectionDiscrepancyScore>> scoreFutures = new ArrayList<>();
						for (PointSet sample : samples) {
							if (firstPointSets[m] == null)
								firstPointSets[m] = sample;
							scoreFutures.add(CompletableFuture.supplyAsync(()->scorer.score(sample, scoreOrders)));
						}
						
						List<ProjectionDiscrepancyScore> scores = scoreFutures.stream().map(F->F.join()).toList();
						
						methodScores.get(m).add(scores);
						
						if (method == SamplingMethod.OWEN_SCRAMBLED_SOBOL) {
							// rebuild it to remove the row scrambling
							firstPointSets[m] = method.createGenerator(baseRand.nextLong()).generate(sampleCount, dimensions);
						}
						
						List<XY_DataSet> funcs = new ArrayList<>();
						List<PlotCurveCharacterstics> chars = new ArrayList<>();

						double[][][] scores2D = new double[dimensions][dimensions][myTrials];
						for (int i=0; i<dimensions; i++)
							for (int j=0; j<dimensions; j++)
								Arrays.fill(scores2D[i][j], Double.NaN);
						double avgScore2D = Double.NaN;
						
						List<UncertainArbDiscFunc> shadedFuncs = new ArrayList<>();
						List<PlotCurveCharacterstics> shadedChars = new ArrayList<>();
						
						for (int order=1; order<=scoreOrders; order++) {
							double overallAverage = 0d;
							double[][] dimScores = new double[dimensions][scores.size()];
							double[] dimAverages = new double[dimensions];
							
							for (int s=0; s<scores.size(); s++) {
								ProjectionDiscrepancyScore score = scores.get(s);
								overallAverage += score.getOrderMeanScore(order);
								int[] dimCounts = new int[dimensions];
								for (ProjectionResult proj : score.getProjectionResults()) {
									if (proj.getProjection().order() != order)
										continue;
									double dimScore = proj.getNormalizedScore();
									int[] dims = proj.getProjection().getDimensions();
									Preconditions.checkState(dims.length == order);
									for (int d : dims) {
										dimScores[d][s] += dimScore;
										dimCounts[d]++;
									}
									
									if (order == 2) {
										Preconditions.checkState(Double.isNaN(scores2D[dims[0]][dims[1]][s]));
										Preconditions.checkState(Double.isNaN(scores2D[dims[1]][dims[0]][s]));
										scores2D[dims[0]][dims[1]][s] = dimScore;
										scores2D[dims[1]][dims[0]][s] = dimScore;
									}
								}
								for (int d=0; d<dimensions; d++) {
									Preconditions.checkState(dimCounts[d] >= 1);
									dimScores[d][s] /= dimCounts[d];
									dimAverages[d] += dimScores[d][s];
								}
								if (s < numPlotTrials) {
									// plot it
									EvenlyDiscretizedFunc dimFunc = new EvenlyDiscretizedFunc(1d, dimensions, 1d);
									for (int d=0; d<dimensions; d++)
										dimFunc.set(d, dimScores[d][s]);
//									if (s == 0)
//										dimFunc.setName("Individual samples");
									funcs.add(dimFunc);
									chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, oderLightColors[order-1]));
								}
							}
							if (scores.size() > 1) {
								EvenlyDiscretizedFunc upperDimFunc = new EvenlyDiscretizedFunc(1d, dimensions, 1d);
								EvenlyDiscretizedFunc middleDimFunc = new EvenlyDiscretizedFunc(1d, dimensions, 1d);
								EvenlyDiscretizedFunc lowerDimFunc = new EvenlyDiscretizedFunc(1d, dimensions, 1d);
								for (int d=0; d<dimensions; d++) {
									double lower = StatUtils.percentile(dimScores[d], 2.5d);
									lowerDimFunc.set(d, lower);
									double upper = StatUtils.percentile(dimScores[d], 97.5d);
									upperDimFunc.set(d, upper);
									middleDimFunc.set(d, 0.5*(upper+lower));
								}
								shadedFuncs.add(new UncertainArbDiscFunc(middleDimFunc, lowerDimFunc, upperDimFunc));
								shadedChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, oderLightColors[order-1]));
							}
							
							overallAverage /= scores.size();
							DefaultXY_DataSet overallFunc = new DefaultXY_DataSet(1d, overallAverage, dimensions, overallAverage);
//							overallFunc.setName("Overall average");
							overallFunc.setName(order+"D");
							funcs.add(overallFunc);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, orderColors[order-1]));
							
							for (int d=0; d<dimensions; d++)
								dimAverages[d] /= scores.size();
							EvenlyDiscretizedFunc dimFunc = new EvenlyDiscretizedFunc(1d, dimensions, 1d);
							for (int d=0; d<dimensions; d++)
								dimFunc.set(d, dimAverages[d]);
//							dimFunc.setName("Single-dimension averages");
							funcs.add(dimFunc);
							chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, oderLightColors[order-1]));
							
							if (order == 2)
								avgScore2D = overallAverage;
						}
						
						funcs.addAll(shadedFuncs);
						chars.addAll(shadedChars);
						
						PlotSpec plot = new PlotSpec(funcs, chars, method.getShortName()+" (N="+sampleCount+")", "Dimension #", "Normalized projection score");
						plot.setLegendVisible(true);
						double xTick = dimensions > 20 ? 2d : 1d;
						
						HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
						
						gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
						
						gp.drawGraphPanel(plot, false, true, dimXRange, logYRange);
						PlotUtils.setXTick(gp, xTick);
						
						if (replotIndvSamples || !new File(subDir, "scores_"+prefix+".png").exists())
							PlotUtils.writePrintPlots(subDir, "scores_"+prefix, gp,
									PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 3d, 300, true, true, false);
						
						EvenlyDiscrXYZ_DataSet avgXYZ = new EvenlyDiscrXYZ_DataSet(dimensions, dimensions, 1d, 1d, 1d);
						EvenlyDiscrXYZ_DataSet avgAbsXYZ = new EvenlyDiscrXYZ_DataSet(dimensions, dimensions, 1d, 1d, 1d);
						CPT logRatioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
						logRatioCPT.setLog10(true);
						CPT logAbsCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d).trim(0d, 1d);
//						CPT logAbsCPT = logRatioCPT.trim(0d, 1d);
						logAbsCPT.setLog10(true);
						for (int i=0; i<dimensions; i++) {
							for (int j=0; j<dimensions; j++) {
								if (i == j) {
									avgXYZ.set(i, j, Double.NaN);
									avgAbsXYZ.set(i, j, Double.NaN);
									continue;
								}
								double[] myScores = scores2D[i][j];
								double sum = 0d;
								double absSum = 0d;
								for (double score : myScores) {
									double ratio = score / avgScore2D;
									sum += ratio;
									absSum += Math.max(ratio, 1/ratio);
								}
								avgXYZ.set(i, j, sum/scores.size());
								avgAbsXYZ.set(i, j, absSum/scores.size());
							}
						}
						
						XYZPlotSpec xyzPlot = new XYZPlotSpec(avgXYZ, logRatioCPT, plot.getTitle(),
								"Dimension #", "Dimension #", "Average pair score / overall 2D score");
						
						Range xyzRange = new Range(0.5, dimensions+0.5);
						gp.drawGraphPanel(xyzPlot, false, false, xyzRange, xyzRange);
						
						if (replotIndvSamples || !new File(subDir, "scores_2D_"+prefix+".png").exists())
							PlotUtils.writePrintPlots(subDir, "scores_2D_"+prefix, gp,
								PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, false, 300, true, true, false);
						
						xyzPlot = new XYZPlotSpec(avgAbsXYZ, logAbsCPT, plot.getTitle(),
								"Dimension #", "Dimension #", "Average realization pair score factor");
						
						gp.drawGraphPanel(xyzPlot, false, false, xyzRange, xyzRange);
						
						if (replotIndvSamples || !new File(subDir, "scores_2D_"+prefix+"_deviation.png").exists())
							PlotUtils.writePrintPlots(subDir, "scores_2D_"+prefix+"_deviation", gp,
								PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, false, 300, true, true, false);
					}
					if (redoCenteredDiscrepancies) {
						System.out.println("Calculating centered discrepancy scores");
						List<CompletableFuture<Double>> scoreFutures = new ArrayList<>();
						for (PointSet sample : samples) {
							scoreFutures.add(CompletableFuture.supplyAsync(()->{
								double score = CenteredDiscrepancy.score(sample);
//								System.out.println("Score: "+(float)score);
								return score;
							}));
						}
						
						List<Double> scores = scoreFutures.stream().map(F->F.join()).toList();
						
						methodCenteredDiscrepancyScores.get(m).add(scores);
					}
					
					methodWatch.stop();
					System.out.println("\tDONE in "+timeStr(methodWatch));
				}
				
				if (sampleCount <= 1024 && redoNormScores) {
					// now plot 2D scatters
					int[][] plotDims = {
							{0, 1},
							{2, 3},
							{4, 5}
					};
					for (int p=0; p<plotDims.length; p++) {
						int dim1 = plotDims[p][0];
						int dim2 = plotDims[p][1];
						
						List<PlotSpec> plots = new ArrayList<>();
						Range range = new Range(0d, 1d);
						List<Range> xRanges = new ArrayList<>();
						List<Range> yRanges = List.of(range);
						List<String> subtitles = new ArrayList<>();
						for (int m=0; m<methods.length; m++) {
							List<XY_DataSet> funcs = new ArrayList<>();
							List<PlotCurveCharacterstics> chars = new ArrayList<>();
							PointSet sample = firstPointSets[m];
							DefaultXY_DataSet xy = new DefaultXY_DataSet();
							for (int i=0; i<sampleCount; i++)
								xy.set(sample.get(i, dim1), sample.get(i, dim2));
							funcs.add(xy);
							chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1f, Color.BLACK));
							
							PlotSpec plot = new PlotSpec(funcs, chars, null, null, null);
							plots.add(plot);
							xRanges.add(range);
							subtitles.add(methods[m].getShortName());
						}

						
						HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
						PlotPreferences prefs = gp.getPlotPrefs();
						prefs.setPlotPadding(new RectangleInsets(0, 5, 5, 5));
						
						gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
						PlotUtils.setAxisVisible(gp, false, false);
						
						Font subtitleFont = new Font(Font.SANS_SERIF, Font.BOLD, 10);
						
						PlotUtils.addSubplotTitles(gp, subtitles, subtitleFont);
						
						String prefix = "scatters_2D";
						if (p > 0)
							prefix += "_"+dim1+"_"+dim2;
						
						if (replotIndvSamples || !new File(subDir, prefix+".png").exists())
							PlotUtils.writePrintPlots(subDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH, false, 300, true, true, false);
					}
				}
				sampleWatch.stop();
				System.out.println("DONE with "+sampleCount+" in "+timeStr(sampleWatch)+"\n");
			}
			
			totalWatch.stop();
			System.out.println("DONE with all calculations in "+timeStr(totalWatch));
		}
		
		// now combined plots
		String prefix = "combined_scores";
		
		List<XY_DataSet> scoreFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> scoreChars = new ArrayList<>();
		List<XY_DataSet> equivCountFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> equivCountChars = new ArrayList<>();
		List<XY_DataSet> centeredFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> centeredChars = new ArrayList<>();
		CSVFile<String> scoresCSV = redoNormScores ? new CSVFile<>(true) : CSVFile.readFile(scoresCSVFile, true);
		CSVFile<String> centeredScoresCSV = redoCenteredDiscrepancies ? new CSVFile<>(true) : CSVFile.readFile(centeredDiscrepancyCSVFile, true);
		if (redoNormScores) {
			List<String> header = new ArrayList<>();
			header.add("");
			for (int order=1; order<=scoreOrders; order++)
				for (int s=0; s<sampleCounts.length; s++)
					header.add(order+"D "+sampleCounts[s]);
			scoresCSV.addLine(header);
		}
		if (redoCenteredDiscrepancies) {
			List<String> header = new ArrayList<>();
			header.add("");
			for (int s=0; s<sampleCounts.length; s++)
				header.add(sampleCounts[s]+"");
			centeredScoresCSV.addLine(header);
		}
		int rowIndex = 1;
		for (int m=0; m<methods.length; m++) {
			SamplingMethod method = methods[m];
			PlotCurveCharacterstics methodChar = combPlotChars.get(method);
			if (methodChar == null)
				continue;

			List<String> scoreLine = new ArrayList<>();
			scoreLine.add(method.getShortName());
			List<String> centeredLine = new ArrayList<>();
			centeredLine.add(method.getShortName());
			
			int colIndex = 1;
			
			for (int order=1; order<=scoreOrders; order++) {
				EvenlyDiscretizedFunc scoreFunc = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
//				EvenlyDiscretizedFunc scoreLowerFunc = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
//				EvenlyDiscretizedFunc scoreUpperFunc = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
				EvenlyDiscretizedFunc equivFunc = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
				
				for (int s=0; s<sampleCounts.length; s++) {
					int sampleCount = sampleCounts[s];
					
					double avg;
					if (redoNormScores) {
						List<ProjectionDiscrepancyScore> scores = methodScores.get(m).get(s);
						double sum = 0d;
//						double min = Double.POSITIVE_INFINITY;
//						double max = 0d;
						for (ProjectionDiscrepancyScore score : scores) {
							double orderScore = score.getOrderMeanScore(order);
							sum += orderScore;
//							min = Math.min(min, orderScore);
//							max = Math.max(max, orderScore);
						}
						avg = sum / scores.size();
					} else {
						avg = scoresCSV.getDouble(rowIndex, colIndex++);
					}
					scoreLine.add((float)avg+"");

					scoreFunc.set(s, avg);
//					scoreLowerFunc.set(s, min);
//					scoreUpperFunc.set(s, max);
					double equivCount = (double)sampleCount / avg;
					equivFunc.set(s, equivCount);
				}
				
				if (order == 1) {
					scoreFunc.setName(method.getShortName());
					equivFunc.setName(method.getShortName());
				} else if (method == SamplingMethod.MONTE_CARLO) {
					// they all overlap, cleaner to just show 1D
					continue;
				}
				
				double thickness = orderThicknessFunc.applyAsDouble(order);
				scoreFuncs.add(scoreFunc);
				scoreChars.add(getForThickness(methodChar, thickness));

//				UncertainArbDiscFunc rangeFunc = new UncertainArbDiscFunc(scoreFunc, scoreLowerFunc, scoreUpperFunc);
//				rangeFunc.setName(null);
//				scoreFuncs.add(0, rangeFunc);
//				scoreChars.add(0, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, methodTransColor));
				
				equivCountFuncs.add(equivFunc);
				equivCountChars.add(getForThickness(methodChar, thickness));
			}
			
			EvenlyDiscretizedFunc centeredFunc = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
			centeredFunc.setName(method.getShortName());
			for (int s=0; s<sampleCounts.length; s++) {
				double centered;
				if (redoCenteredDiscrepancies) {
					centered = methodCenteredDiscrepancyScores.get(m).get(s).stream().mapToDouble(d->d).average().getAsDouble();
					centeredLine.add((float)centered+"");
				} else {
					centered = centeredScoresCSV.getDouble(rowIndex, s+1);
				}
				centeredFunc.set(s, centered);
			}
			centeredFuncs.add(centeredFunc);
			centeredChars.add(getForThickness(methodChar, 3f));
			
			rowIndex++;
			if (redoNormScores)
				scoresCSV.addLine(scoreLine);
			if (redoCenteredDiscrepancies)
				centeredScoresCSV.addLine(centeredLine);
		}
		List<XY_DataSet> orderTicknessFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> orderThicknessChars = new ArrayList<>();
		for (int order=1; order<=scoreOrders; order++) {
			double thickness = orderThicknessFunc.applyAsDouble(order);
			
			XY_DataSet fakeXY = new DefaultXY_DataSet(-100d, 1d);
			fakeXY.setName(order+"D");
			PlotCurveCharacterstics orderChar = new PlotCurveCharacterstics(PlotLineType.SOLID, (float)thickness, Color.GRAY);
			orderTicknessFuncs.add(fakeXY);
			orderThicknessChars.add(orderChar);
		}
		PlotSpec orderPlot = new PlotSpec(orderTicknessFuncs, orderThicknessChars, null, null, null);
		
		HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
		PlotPreferences prefs = gp.getPlotPrefs();
		prefs.setPlotLabelFontSize(10);
		prefs.setLegendFontSize(8);
		prefs.setLegendLineLength(8d);
		prefs.getPlotPadding();
		prefs.setPlotPadding(new RectangleInsets(4, 0, 0, 12));
		
		gp.drawGraphPanel(orderPlot, false, false);
		LegendItemCollection orderLegendItems = gp.getPlot().getLegendItems();
		
//		EvenlyDiscretizedFunc equivLinear = new EvenlyDiscretizedFunc(0d, sampleCounts.length, 1d);
//		for (int s=0; s<sampleCounts.length; s++)
//			equivLinear.set(s, sampleCounts[s]);
//		equivLinear.setName("Linear");
//		equivCountFuncs.add(0, equivLinear);
//		equivCountChars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		Range categoricalXRange = new Range(0d, sampleCounts.length-1);
		
		NumberFormat categoryFormat = new NumberFormat() {
		    @Override
		    public StringBuffer format(double value, StringBuffer buffer, FieldPosition pos) {
		        int index = (int)Math.round(value);
		        if (index >= 0 && index < sampleCounts.length
		                && Math.abs(value - index) < 1e-6)
		            buffer.append(sampleCounts[index]);
		        return buffer;
		    }

		    @Override
		    public StringBuffer format(long value, StringBuffer buffer, FieldPosition pos) {
		        return format((double)value, buffer, pos);
		    }

		    @Override
		    public Number parse(String source, ParsePosition pos) {
		        pos.setErrorIndex(pos.getIndex());
		        return null;
		    }
		};
		
		PlotSpec plot = new PlotSpec(scoreFuncs, scoreChars, treeName, "Sample count", "Normalized projection score");
//		plot.setLegendInset(true);
		plot.setLegendVisible(true);
		
		orderPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		plot.addPlotAnnotation(orderPlot.buildInsetLegend(orderLegendItems, prefs, false, true, categoricalXRange, logYRange));
		
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(plot, false, true, categoricalXRange, logYRange);
		PlotUtils.setXTick(gp, 1);
		((NumberAxis)gp.getXAxis()).setNumberFormatOverride(categoryFormat);
		
		PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 4, 300, true, true, false);
		if (redoNormScores)
			scoresCSV.writeToFile(scoresCSVFile);
		
		plot = new PlotSpec(equivCountFuncs, equivCountChars, treeName, "Sample count", "Equivalent MCS count");
//		plot.setLegendInset(true);
		plot.setLegendVisible(true);
		
		orderPlot.setLegendInset(RectangleAnchor.TOP_LEFT);
		plot.addPlotAnnotation(orderPlot.buildInsetLegend(orderLegendItems, prefs, false, true, categoricalXRange, equivYRange));
		
		gp.drawGraphPanel(plot, false, true, categoricalXRange, equivYRange);
		PlotUtils.setXTick(gp, 1);
		((NumberAxis)gp.getXAxis()).setNumberFormatOverride(categoryFormat);
		
		prefix = "combined_equivs";
		PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 4, 300, true, true, false);
		
		// now centered
		plot = new PlotSpec(centeredFuncs, centeredChars, treeName, "Sample count", "Squared centered discrepancy");
//		plot.setLegendInset(true);
		plot.setLegendVisible(true);
		
		gp.drawGraphPanel(plot, false, true, categoricalXRange, centeredYRange);
		PlotUtils.setXTick(gp, 1);
		((NumberAxis)gp.getXAxis()).setNumberFormatOverride(categoryFormat);
		
		prefix = "combined_centered_discrepancies";
		PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 4, 300, true, true, false);
		
		if (redoCenteredDiscrepancies)
			centeredScoresCSV.writeToFile(centeredDiscrepancyCSVFile);
	}
	
	private static PlotCurveCharacterstics getForThickness(PlotCurveCharacterstics pChar, double thickness) {
		PlotCurveCharacterstics copy = (PlotCurveCharacterstics)pChar.clone();
		copy.setLineWidth((float)thickness);
		return copy;
	}
	
	private static final DecimalFormat timeDF = new DecimalFormat("0.0");
	private static String timeStr(Stopwatch watch) {
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		if (secs < 90d)
			return timeDF.format(secs)+" s";
		double mins = secs/60d;
		if (mins < 90d)
			return timeDF.format(mins)+" m";
		double hours = mins / 60d;
		return timeDF.format(hours)+" h";
	}

}
