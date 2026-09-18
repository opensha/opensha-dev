package scratch.kevin.nshm27.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateModel;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateFileLoader.Exact;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateFileLoader.PureGR;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateFileLoader.RateRecord;
import org.opensha.sha.earthquake.nshmp.seismicity.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_SeisClassificationMethod;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_SeisRateModelBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_SeisRateModelSamples;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

public class ObsUncertaintyBoundsFigure {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "obs_mfd_bounds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		boolean incremental = false;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(3.05, 10.95);
		
		Color[] colors = {
				Colors.tab_blue,
				Colors.tab_orange,
				Colors.tab_green
		};
		
		PlotLineType[] lineTypes = {
				PlotLineType.DASHED,
				PlotLineType.SOLID,
				PlotLineType.DOTTED
		};
		
		PlotSymbol[] symbolTypes = {
				PlotSymbol.FILLED_SQUARE,
				PlotSymbol.FILLED_CIRCLE,
				PlotSymbol.FILLED_INV_TRIANGLE
		};
		
		RateType[] types = {
			RateType.M1,
			RateType.M1_TO_MMAX,
			RateType.EXACT
		};
		
		boolean includeWtMean = false;
		
		NSHM27_SeismicityRegions[] seisRegions = NSHM27_SeismicityRegions.values();
		TectonicRegionType[] trts = {TectonicRegionType.SUBDUCTION_INTERFACE, TectonicRegionType.SUBDUCTION_SLAB, TectonicRegionType.ACTIVE_SHALLOW};
		
		double weightLow = NSHM27_SeisRateModelBranch.LOW.getNodeWeight();
		double weightPref = NSHM27_SeisRateModelBranch.PREFFERRED.getNodeWeight();
		double weightHigh = NSHM27_SeisRateModelBranch.HIGH.getNodeWeight();
		
		for (NSHM27_SeismicityRegions seisReg : seisRegions) {
			for (NSHM27_SeisClassificationMethod classification : NSHM27_SeisClassificationMethod.values()) {
				for (TectonicRegionType trt : trts) {
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					String prefix = seisReg.name()+"_"+trt.name()+"_"+classification.name();
					String title = seisReg.getShortName()+" ("+NSHM27_RegionLoader.getNameForTRT(trt)+", "+classification.getShortName()+")";
					
					DecimalFormat oDF = new DecimalFormat("0.#");

					Double m1 = null;
					Double mMax = null;
					Map<RateType, EvenlyDiscretizedFunc[]> typeMagFuncs = new HashMap<>();
					EvenlyDiscretizedFunc overallMean = null;
					SeismicityRateModel[] rateModels = new SeismicityRateModel[types.length];
					for (int t=0; t<types.length; t++)
						rateModels[t] = NSHM27_SeisRateModelBranch.loadRateModel(seisReg, classification, trt, types[t]);
					for (int t=0; t<types.length; t++) {
						RateType type = types[t];
						SeismicityRateModel rateModel = rateModels[t];
						
						RateRecord meanRec = rateModel.getMeanRecord();
						if (m1 == null)
							m1 = meanRec.M1;
						if (mMax == null && meanRec instanceof PureGR)
							mMax = ((PureGR)meanRec).Mmax;
						
						EvenlyDiscretizedFunc meanMFD;
						if (incremental)
							meanMFD = SeismicityRateFileLoader.buildIncrementalMFD(meanRec, refMFD, refMFD.getMaxX());
						else
							meanMFD = cmlMFD(meanRec, refMFD);
						
						if (funcs.isEmpty()) {
							meanMFD.setName("Mean");
							funcs.add(meanMFD);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
							overallMean = meanMFD;
						}
						
						Color color = colors[type.ordinal()];
						PlotLineType plt = lineTypes[type.ordinal()];
						RateRecord low = rateModel.getLowerRecord();
						RateRecord high = rateModel.getUpperRecord();

						EvenlyDiscretizedFunc lowMFD;
						EvenlyDiscretizedFunc highMFD;
						if (incremental) {
							lowMFD = SeismicityRateFileLoader.buildIncrementalMFD(low, refMFD, refMFD.getMaxX());
							highMFD = SeismicityRateFileLoader.buildIncrementalMFD(high, refMFD, refMFD.getMaxX());
						} else {
							lowMFD = cmlMFD(low, refMFD);
							highMFD = cmlMFD(high, refMFD);
						}
						
						typeMagFuncs.put(type, new EvenlyDiscretizedFunc[] {lowMFD, highMFD});
						
						lowMFD.setName(typeName(type));
						funcs.add(lowMFD);
						chars.add(new PlotCurveCharacterstics(plt, 3f, color));
						highMFD.setName(null);
						funcs.add(highMFD);
						chars.add(new PlotCurveCharacterstics(plt, 3f, color));
						
						if (includeWtMean && type == RateType.M1_TO_MMAX) {
							EvenlyDiscretizedFunc weightAvg = new EvenlyDiscretizedFunc(meanMFD.getMinX(), meanMFD.size(), meanMFD.getDelta());
							Preconditions.checkState((float)meanMFD.getMinX() == (float)lowMFD.getMinX());
							Preconditions.checkState((float)meanMFD.getMinX() == (float)highMFD.getMinX());
							for (int i=0; i<weightAvg.size(); i++)
								weightAvg.set(i, weightLow*lowMFD.getY(i) + weightPref*meanMFD.getY(i) + weightHigh*highMFD.getY(i));
							
							if (!incremental) {
								System.out.println(title);
								System.out.println("\tM>5: "+(float)weightAvg.getY(weightAvg.getClosestXIndex(5.01)));
								System.out.println("\tM>6: "+(float)weightAvg.getY(weightAvg.getClosestXIndex(6.01)));
								System.out.println("\tM>6 snapped: "+(float)weightAvg.getX(weightAvg.getClosestXIndex(6.01)));
								System.out.println("\tM1="+m1.floatValue());
								System.out.println("\tMmax="+mMax.floatValue());
							}
							
							weightAvg.setName(typeName(type)+" Average");
							funcs.add(weightAvg);
							chars.add(new PlotCurveCharacterstics(plt, 3f, Color.DARK_GRAY));
						}
					}
					
					for (XY_DataSet func : funcs) {
						if (func.getName() != null && func.getName().contains("M1"))
							func.setName(func.getName().replace("M1", "M₁"));
						if (func.getName() != null && func.getName().contains("Mmax"))
							func.setName(func.getName().replace("Mmax", "Mₘₐₓ"));
					}
					
					Range xRange = new Range(4d, 8d);
					Range yRange = incremental ? new Range(1e-4, 1e2) : new Range(1e-3, 1e3);
					
					List<XYTextAnnotation> anns = new ArrayList<>();
					Font annFont = new Font(Font.SANS_SERIF, Font.PLAIN, 22);
					
					DefaultXY_DataSet m1Line = new DefaultXY_DataSet();
					m1Line.set(m1, yRange.getLowerBound());
					m1Line.set(m1, yRange.getUpperBound());
					funcs.add(m1Line);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
					
					DefaultXY_DataSet mMaxLine = new DefaultXY_DataSet();
					mMaxLine.set(mMax, yRange.getLowerBound());
					mMaxLine.set(mMax, yRange.getUpperBound());
					funcs.add(mMaxLine);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
					
					XYTextAnnotation m1Ann = new XYTextAnnotation(" M₁=5", m1, yRange.getUpperBound());
					m1Ann.setFont(annFont);
					m1Ann.setTextAnchor(TextAnchor.TOP_LEFT);
					anns.add(m1Ann);
					
					XYTextAnnotation mMaxAnn = new XYTextAnnotation("Mₘₐₓ="+mMax.floatValue()+" ", mMax, yRange.getLowerBound());
					mMaxAnn.setFont(annFont);
					mMaxAnn.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
					anns.add(mMaxAnn);
					
					PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", incremental ? "Incremental Rate (1/yr)" : "Cumulative Rate (1/yr)");
					plot.setLegendInset(true);
					plot.setPlotAnnotations(anns);
					
					HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
					
					gp.drawGraphPanel(plot, false, true, xRange, yRange);
					
					PlotUtils.writePlots(outputDir, prefix, gp, 700, 650, true, true, false);
					
					if (!incremental) {
						List<PureGR> samples = new NSHM27_SeisRateModelSamples(seisReg, trt).loadOrigSamples(classification);
						Collections.shuffle(samples, new Random(samples.size()));
						
//						int c = 200;
//						int a = 127;
//						int c = 150;
//						int a = 80;
						int c = 180;
						int a = 60;
						PlotCurveCharacterstics indvChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(c, c, c, a));
						int maxNumRates = 1000;
						int numRates = Integer.min(maxNumRates, samples.size());
						GutenbergRichterMagFreqDist[] rateMFDs = new GutenbergRichterMagFreqDist[samples.size()];
						for (int i=0; i<rateMFDs.length; i++) {
							PureGR pair = samples.get(i);
							GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX()-0.5*refMFD.getDelta(), refMFD.size(), refMFD.getDelta());
							// this sets shape, min/max
							// subtract a tiny amount from mMax so that if it's exactly at a bin edge, e.g. 7.9, it rounds down, e.g. to 7.85
							gr.setAllButTotCumRate(gr.getX(0), gr.getMaxX(), 1e16, pair.b);
							// this scales it to match
							// similarly, add a tiny amount to M1 so that if it's exactly at a bin edge (which it should be as it's determined
							// using cumulative binning), it rounds up to the incremental bin for that cumulative edge
//							gr.scaleToCumRate(refMFD.getClosestXIndex(m1+0.001), pair[0]);
							// this is actually cumulative, but the GR dist object works in incremental space
							gr.scaleToIncrRate(gr.getClosestXIndex(m1+0.001), pair.rateAboveM1);
							gr.setName(null);
							
							rateMFDs[i] = gr;
							
							if (i < numRates) {
								funcs.add(0, gr);
								chars.add(0, indvChar);
								
								if (i == 0) {
									DiscretizedFunc forLegend = gr.deepClone();
									forLegend.scale(1e-10);
									forLegend.setName("Individual samples");
									funcs.add(0, forLegend);
									chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(c, c, c))); // alpha disabled
								}
							}
						}
						
						gp.drawGraphPanel(plot, false, true, xRange, yRange);
						
						PlotUtils.writePlots(outputDir, prefix+"_with_indv", gp, 700, 650, true, true, false);
						
						// write histogram
						double[] histMags = {5d, 6d, 7d};
						
						for (double histMag : histMags) {
							double minRate = Double.POSITIVE_INFINITY;
							double maxRate = Double.NEGATIVE_INFINITY;
							
							double[] typeLowers = new double[types.length];
							double[] typeUppers = new double[types.length];
							
							for (int t=0; t<types.length; t++) {
								RateType type = types[t];
								EvenlyDiscretizedFunc[] typeFuncs = typeMagFuncs.get(type);
								int index = typeFuncs[0].getClosestXIndex(histMag);
								Preconditions.checkState(Precision.equals(histMag, typeFuncs[0].getX(index), 1e-4));
								typeLowers[t] = typeFuncs[0].getY(index);
								typeUppers[t] = typeFuncs[1].getY(index);
								minRate = Math.min(minRate, typeLowers[t]);
								maxRate = Math.max(maxRate, typeUppers[t]);
							}
							
							double[] histValues = new double[samples.size()];
							for (int r=0; r<histValues.length; r++) {
								int index = rateMFDs[r].getClosestXIndex(histMag);
								Preconditions.checkState(Precision.equals(histMag, rateMFDs[r].getX(index), 1e-4),
										"Bad sample mag[%s]=%s for histMag=%s", index, rateMFDs[r].getX(index), histMag);
								histValues[r] = rateMFDs[r].getY(index);
								minRate = Math.min(minRate, histValues[r]);
								maxRate = Math.max(maxRate, histValues[r]);
							}
							
							int index = overallMean.getClosestXIndex(histMag);
							Preconditions.checkState(Precision.equals(histMag, overallMean.getX(index), 1e-4));
							double meanValue = overallMean.getY(index);
							
							
							double logMinRate = Math.log10(minRate) - 0.1;
							double logMaxRate = Math.log10(maxRate) + 0.1;
							double logMeanValue = Math.log10(meanValue);
							// center around mean value
							double maxDiff = Math.max(Math.abs(logMeanValue - logMaxRate), Math.abs(logMeanValue - logMinRate));
							logMinRate = logMeanValue - maxDiff;
							logMaxRate = logMeanValue + maxDiff;
							
							HistogramFunction logHist = HistogramFunction.getEncompassingHistogram(logMinRate, logMaxRate, 0.05);
							
							for (double histValue : histValues)
								logHist.add(logHist.getClosestXIndex(Math.log10(histValue)), 1d);
							
							double maxY = logHist.getMaxY()*1.1d;
							DiscretizedFunc linearHist = new ArbitrarilyDiscretizedFunc();
							for (Point2D pt : logHist)
								linearHist.set(Math.pow(10, pt.getX()), pt.getY());
							
							funcs = new ArrayList<>();
							chars = new ArrayList<>();
							
							linearHist.setName("Sampled Distribution");
							funcs.add(linearHist);
							chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
							
							
							DefaultXY_DataSet meanXY = new DefaultXY_DataSet();
							meanXY.set(meanValue, 0d);
							meanXY.set(meanValue, maxY);
							meanXY.setName("Mean");
							funcs.add(meanXY);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));

							for (int t=0; t<types.length; t++) {
								RateType type = types[t];
								Color color = colors[type.ordinal()];
								PlotLineType plt = lineTypes[type.ordinal()];
								
								DefaultXY_DataSet lowerXY = new DefaultXY_DataSet();
								lowerXY.set(typeLowers[t], 0d);
								lowerXY.set(typeLowers[t], maxY);
								lowerXY.setName(typeName(type));
								
								DefaultXY_DataSet upperXY = new DefaultXY_DataSet();
								upperXY.set(typeUppers[t], 0d);
								upperXY.set(typeUppers[t], maxY);
								
								funcs.add(lowerXY);
								chars.add(new PlotCurveCharacterstics(plt, 3f, color));
								
								funcs.add(upperXY);
								chars.add(new PlotCurveCharacterstics(plt, 3f, color));
							}
							
							plot = new PlotSpec(funcs, chars, " ", "M>"+oDF.format(histMag)+" Rate", "Sample Count");
							plot.setLegendInset(true);
							
							gp.drawGraphPanel(plot, true, false, new Range(Math.pow(10, logMinRate), Math.pow(10, logMaxRate)), new Range(0d, maxY));
							
							PlotUtils.writePlots(outputDir, prefix+"_hist_m"+oDF.format(histMag), gp, 700, 650, true, true, false);
						}
						
						// write rate vs b
						MinMaxAveTracker rateTrack = new MinMaxAveTracker();
						MinMaxAveTracker bTrack = new MinMaxAveTracker();
						double[] bValues = new double[samples.size()];
						double[] rates = new double[samples.size()];
						for (int s=0; s<samples.size(); s++) {
							PureGR gr = samples.get(s);
							rateTrack.addValue(gr.rateAboveM1);
							bTrack.addValue(gr.b);
							bValues[s] = gr.b;
							rates[s] = gr.rateAboveM1;
							
						}
						int numBinsEach = 20;
						double rateLength = rateTrack.getLength();
						double bLenth = bTrack.getLength();
						double rateDelta = rateLength/(numBinsEach-1);
						double bDelta = bLenth/(numBinsEach-1);
						
//						double bMin = bTrack.getMin()-0.5*bDelta;
//						double rateMin = Math.max(0, rateTrack.getMin()-0.5*rateDelta);
						double bMin = bTrack.getMin();
						double rateMin = rateTrack.getMin();
						
						System.out.println("b range: "+bTrack);
						System.out.println("\tbinMin="+bMin);
						System.out.println("\tbinDelta="+bDelta);
						System.out.println("rate range: "+rateTrack);
						System.out.println("\tbinMin="+rateMin);
						System.out.println("\tbinDelta="+rateDelta);
						
						float symbolWidth = 7f;
						
						EvenlyDiscrXYZ_DataSet rateVsB = new EvenlyDiscrXYZ_DataSet(numBinsEach, numBinsEach,
								bMin, rateMin, bDelta, rateDelta);
						
						for (PureGR gr : samples)
							rateVsB.add(gr.b, gr.rateAboveM1, 1d);
						
						funcs = new ArrayList<>();
						chars = new ArrayList<>();
						
						DefaultXY_DataSet mean = new DefaultXY_DataSet();
						mean.setName("Mean");
						
						mean.set(StatUtils.mean(bValues), StatUtils.mean(rates));
						
						funcs.add(mean);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, symbolWidth, Color.BLACK));
						
						XY_DataSet meanCopy = mean.deepClone();
						meanCopy.setName(null);
						funcs.add(meanCopy);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, symbolWidth, Color.WHITE));
						
						for (int t=0; t<types.length; t++) {
							RateType type = types[t];
							SeismicityRateModel rateModel = rateModels[t];
							PlotSymbol symbol = symbolTypes[type.ordinal()];
							Color color = colors[type.ordinal()];
							
							RateRecord lower = rateModel.getLowerRecord();
							RateRecord upper = rateModel.getUpperRecord();
							if (lower instanceof PureGR) {
								DefaultXY_DataSet xy = new DefaultXY_DataSet();
								xy.setName(typeName(type));
								
								PureGR gr = (PureGR)lower;
								xy.set(gr.b, gr.rateAboveM1);
								gr = (PureGR)upper;
								xy.set(gr.b, gr.rateAboveM1);
								
								funcs.add(xy);
								chars.add(new PlotCurveCharacterstics(symbol, symbolWidth, color));
							}
						}
						
						double lowerB = StatUtils.percentile(bValues, 2.5d);
						double upperB = StatUtils.percentile(bValues, 97.5d);
						double lowerRate = StatUtils.percentile(rates, 2.5d);
						double upperRate = StatUtils.percentile(rates, 97.5d);
						
						DefaultXY_DataSet xy = new DefaultXY_DataSet();
						xy.setName("95% marginal bounds");
						
						xy.set(lowerB, lowerRate);
						xy.set(upperB, upperRate);
						
						funcs.add(xy);
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_INV_TRIANGLE, symbolWidth, Colors.tab_red));
						
						gp = PlotUtils.initPrintHeadless();
						int tickLabelFontSize = gp.getPlotPrefs().getTickLabelFontSize();
						gp.getPlotPrefs().setPlotPadding(new RectangleInsets(tickLabelFontSize/2d, 5, 0, tickLabelFontSize+4));
//						gp.getPlotPrefs().setCptPadding(10);
						
						rateVsB.scale(1d/samples.size());
						CPT cpt = GMT_CPT_Files.SEQUENTIAL_OSLO_UNIFORM.instance().reverse().rescale(0d, rateVsB.getMaxZ());
						XYZPlotSpec xyzPlot = new XYZPlotSpec(rateVsB, funcs, chars, cpt, title,
								"b-value", "M>"+oDF.format(m1)+" Rate", "Fraction of "+samples.size()+" samples");
						xyzPlot.setLegendInset(RectangleAnchor.TOP_LEFT);
						xyzPlot.setIncludeZlabelInLegend(false);
						
						double ratePlotMax = Math.min(rateVsB.getMaxY()-0.5*rateDelta, StatUtils.percentile(rates, 99.99));
						double ratePlotMin = Math.max(rateVsB.getMinY()+0.5*rateDelta, StatUtils.percentile(rates, 0.01));
						double bPlotMax = Math.min(rateVsB.getMaxX()-0.5*bDelta, StatUtils.percentile(bValues, 99.99));
						double bPlotMin = Math.max(rateVsB.getMinX()+0.5*bDelta, StatUtils.percentile(bValues, 0.01));
						
						gp.drawGraphPanel(xyzPlot, false, false, new Range(bPlotMin, bPlotMax), new Range(ratePlotMin, ratePlotMax));
//								new Range(rateVsB.getMinX()-0.5*rateVsB.getGridSpacingX(), rateVsB.getMaxX()+0.5*rateVsB.getGridSpacingX()),
//								new Range(rateVsB.getMinY()-0.5*rateVsB.getGridSpacingY(), rateVsB.getMaxY()+0.5*rateVsB.getGridSpacingY()));
						
						PlotUtils.writePrintPlots(outputDir, prefix+"_rate_vs_b", gp,
								PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*2d/3d, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*2d/3d, 
								150, true, true, false);
					}
				}
			}
		}
	}
	
	private static String typeName(RateType type) {
		String name = type.toString();
		return name.replace("Branches", "branches");
	}
	
	private static EvenlyDiscretizedFunc cmlMFD(RateRecord record, EvenlyDiscretizedFunc refMFD) {
		if (record.type == RateType.EXACT)
			return ((Exact)record).cumulativeDist;
		Preconditions.checkState(record instanceof PureGR);
		PureGR grRec = (PureGR)record;
		// fake a cml GR
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(
				grRec.b, 1d, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		grMFD.scaleToIncrRate(grMFD.getX(grMFD.getClosestXIndex(grRec.M1+0.01)), grRec.rateAboveM1);

		EvenlyDiscretizedFunc cmlGR = new EvenlyDiscretizedFunc(
				refMFD.getMinX()-0.5*refMFD.getDelta(), refMFD.size(), refMFD.getDelta());
		for (int i=0; i<cmlGR.size(); i++)
			cmlGR.set(i, grMFD.getY(i));
		return cmlGR;
	}

}
