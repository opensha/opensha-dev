package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.MathArrays;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

public class ZScoreHistPlot {
	
	private static final boolean rate_weighted = true;
	
	public static class ZScoreResult {
		
		public final double mean;
		public final double stdDevFract;
		public final HistogramFunction hist;
		public final Map<String, HistogramFunction> sourceHists;
		
		private ZScoreResult(double mean, double stdDevFract, HistogramFunction hist, Map<String, HistogramFunction> sourceHists) {
			super();
			this.mean = mean;
			this.stdDevFract = stdDevFract;
			this.hist = hist;
			this.sourceHists = sourceHists;
		}
	}
	
	public static <E> ZScoreResult[] calcZScores(SimulationRotDProvider<E> simProv,
			Collection<? extends RuptureComparison<E>> eventComps, List<Site> sites, IMT[] imts,
					RuptureComparisonFilter<E> filter) throws IOException {
		return calcZScores(simProv, eventComps, sites, imts, filter, null);
	}
	
	public static <E> ZScoreResult[] calcZScores(SimulationRotDProvider<E> simProv,
			Collection<? extends RuptureComparison<E>> eventComps, List<Site> sites, IMT[] imts,
					RuptureComparisonFilter<E> filter, Table<String, E, Double> sourceRupContribFracts) throws IOException {
		int numComputed = 0;
		int numMatches = 0;
		for (Site site : sites) {
			for (RuptureComparison<E> comp : eventComps) {
				if (!comp.isComputed(site, imts[0]))
					continue;
				numComputed++;
				if ((filter != null && !filter.matches(comp, site)))
					continue;
				numMatches += simProv.getNumSimulations(site, comp.getRupture());
			}
		}
		System.out.println(numMatches+" matches (of "+numComputed+" computed)");
		int numBins;
		if (numMatches < 100)
			numBins = 10;
		else if (numMatches < 500)
			numBins = 40;
		else
			numBins = 100;
		System.out.println("Binning with "+numBins+" bins");
		
		ZScoreResult[] ret = new ZScoreResult[imts.length];
		
		for (int p=0; p<imts.length; p++) {
			IMT imt = imts[p];
			HistogramFunction hist = new HistogramFunction(-numStdDev, numStdDev, numBins);
			
			Map<String, HistogramFunction> sourceHists = sourceRupContribFracts == null ? null : new HashMap<>();
			
			List<Double> allVals = new ArrayList<>();
			List<Double> allWeights = rate_weighted ? new ArrayList<>() : null;
			int count = 0;
			for (Site site : sites) {
				for (RuptureComparison<E> comp : eventComps) {
					if (!comp.hasSite(site))
						continue;
					Preconditions.checkState(comp.isComputed(site, imt),
							"Computed for %s but not %s", imts[0], imt);
					if ((filter != null && !filter.matches(comp, site)))
						continue;
					double gmpeVal = comp.getLogMean(site, imt);
					Map<String, Double> sourceFracts = sourceHists == null ?
							null : sourceRupContribFracts.column(comp.getRupture());
					
					List<Double> simVals = simProv.getValues(site, comp.getRupture(), imt);
					
					double rateEach = comp.getAnnualRate()/(double)simVals.size();
					for (double simVal : simVals) {
						// in log space
						simVal = Math.log(simVal);

						double val = (simVal - gmpeVal)/comp.getStdDev(site, imt);
						
						int ind = hist.getClosestXIndex(val);
						if (sourceFracts != null) {
							double sum = 0d;
							for (String sourceName : sourceFracts.keySet()) {
								HistogramFunction sourceHist = sourceHists.get(sourceName);
								if (sourceHist == null) {
									sourceHist = new HistogramFunction(hist.getMinX(), hist.getMaxX(), hist.size());
									sourceHist.setName(sourceName);
									sourceHists.put(sourceName, sourceHist);
								}
								double fract = sourceFracts.get(sourceName);
								sourceHist.add(ind, rateEach*fract);
								sum += fract;
							}
							Preconditions.checkState((float)sum == 1f, "Bad sum of contribution fracts: %s", sum);
						}
						
						hist.add(ind, rateEach);
						allVals.add(val);
						if (rate_weighted)
							allWeights.add(rateEach);
						count++;
					}
				}
			}
			if (count == 0)
				return null;
			double[] valsArray = Doubles.toArray(allVals);
			double mean, stdDev;
			if (rate_weighted) {
				double[] weightsArray = Doubles.toArray(allWeights);
				weightsArray = MathArrays.normalizeArray(weightsArray, weightsArray.length);
				mean = new Mean().evaluate(valsArray, weightsArray);
				double var = new Variance().evaluate(valsArray, weightsArray, mean);
				stdDev = Math.sqrt(var);
				System.out.println(imt.getDisplayName()+" mean="+(float)mean
						+"\tvar="+(float)var+"\tsd="+(float)stdDev);
			} else {
				mean = new Mean().evaluate(valsArray);
				stdDev = Math.sqrt(new Variance().evaluate(valsArray, mean));
			}
			
			double area = calcArea(hist);
			hist.scale(1d/area);
			
			if (sourceHists != null)
				for (HistogramFunction sourceHist : sourceHists.values())
					sourceHist.scale(1d/area);
			
			ret[p] = new ZScoreResult(mean, stdDev, hist, sourceHists);
		}
		
		return ret;
	}
	
	public static <E> boolean plotStandardNormal(SimulationRotDProvider<E> simProv, Collection<? extends RuptureComparison<E>> eventComps,
			List<Site> sites, IMT[] imts, AttenRelRef gmpe, RuptureComparisonFilter<E> filter, List<String> binDescriptions,
			File outputDir, String prefix) throws IOException {
		return plotStandardNormal(simProv, eventComps, sites, imts, gmpe, filter, binDescriptions, outputDir, prefix, null, 0);
	}
	
	private static final double maxY = 0.7d;
	private static final double numStdDev = 3.75;
	
	public static <E> boolean plotStandardNormal(SimulationRotDProvider<E> simProv, Collection<? extends RuptureComparison<E>> eventComps,
			List<Site> sites, IMT[] imts, AttenRelRef gmpe, RuptureComparisonFilter<E> filter, List<String> binDescriptions,
			File outputDir, String prefix, Table<String, E, Double> sourceRupContribFracts, int maxNumSourceContribs) throws IOException {
		
		List<PlotSpec> specs = new ArrayList<>();
		
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(new Range(-numStdDev, numStdDev));
		
		List<Double> means = new ArrayList<>();
		List<Double> stdDevs = new ArrayList<>();
		
		Color stdDevColor = new Color(0, 150, 0);
		
		DatasetRenderingOrder order = DatasetRenderingOrder.FORWARD;
		
		ZScoreResult[] scores = calcZScores(simProv, eventComps, sites, imts, filter, sourceRupContribFracts);
		if (scores == null)
			return false;
		
		for (int p=0; p<imts.length; p++) {
			IMT imt = imts[p];
			ZScoreResult score = scores[p];
			means.add(score.mean);
			stdDevs.add(score.stdDevFract);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			EvenlyDiscretizedFunc stdNormal = new EvenlyDiscretizedFunc(score.hist.getMinX(), score.hist.getMaxX(), 1000);
			double scalar = 1d/Math.sqrt(2d*Math.PI);
			for (int i=0; i<stdNormal.size(); i++) {
				double x = stdNormal.getX(i);
				double y = scalar*Math.exp(-0.5*x*x);
				stdNormal.set(i, y);
			}
			
			score.hist.setName(simProv.getName());
			stdNormal.setName("Standard Normal");
			
//			maxY = Math.max(maxY, Math.max(stdNormal.getMaxY(), hist.getMaxY()));
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			meanLine.set(score.mean, 0);
			meanLine.set(score.mean, maxY-0.1);
			meanLine.setName("Mean");
			
			if (score.sourceHists != null && !score.sourceHists.isEmpty()) {
				order = DatasetRenderingOrder.REVERSE;
				List<HistogramFunction> sourceHistList = new ArrayList<>();
				for (String sourceName : score.sourceHists.keySet())
					sourceHistList.add(score.sourceHists.get(sourceName));
				// sort by area, decreasing
				Collections.sort(sourceHistList, histComparator);
				if (sourceHistList.size() > maxNumSourceContribs) {
					HistogramFunction otherHist = new HistogramFunction(
							score.hist.getMinX(), score.hist.getMaxX(), score.hist.size());
					otherHist.setName("Other");
					for (int i=maxNumSourceContribs-1; i<sourceHistList.size(); i++) {
						HistogramFunction sourceHist = sourceHistList.get(i);
						for (int j=0; j<otherHist.size(); j++)
							otherHist.add(j, sourceHist.getY(j));
					}
					sourceHistList = sourceHistList.subList(0, maxNumSourceContribs-1);
					sourceHistList.add(otherHist);
				}
				CPT colorCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, Integer.max(sourceHistList.size()-1, 1));
				colorCPT = colorCPT.reverse();
				colorCPT.setBelowMinColor(colorCPT.getMinColor());
				colorCPT.setAboveMaxColor(colorCPT.getMaxColor());
				
				funcs.add(meanLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
				funcs.add(stdNormal);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				
				HistogramFunction runningTotal = new HistogramFunction(score.hist.getMinX(), score.hist.getMaxX(), score.hist.size());
				for (int i=0; i<sourceHistList.size(); i++) {
					HistogramFunction sourceHist = sourceHistList.get(i);
					
					for (int j=0; j<score.hist.size(); j++)
						runningTotal.add(j, sourceHist.getY(j));
					
					EvenlyDiscretizedFunc clone = runningTotal.deepClone();
					clone.setName(sourceHist.getName());
					
					// this will stagger it
					float cptVal = i % 2 == 0 ? (float)i*0.5f :
						(float)(i-1)*0.5f+colorCPT.getMaxValue()*0.5f + (sourceHistList.size() % 2 == 0 ? 0.5f : 1f);
//					System.out.println(i+" => "+cptVal);
					
					funcs.add(clone);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colorCPT.getColor(cptVal)));
				}
				
				funcs.add(score.hist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			} else {
				funcs.add(score.hist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
				funcs.add(stdNormal);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				funcs.add(meanLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
			}
			
			for (double sigma=Math.ceil(-numStdDev); sigma<=numStdDev; sigma++) {
				DefaultXY_DataSet sigmaLine = new DefaultXY_DataSet();
				sigmaLine.set(sigma, 0);
				sigmaLine.set(sigma, maxY-0.075);
				funcs.add(sigmaLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, stdDevColor));
			}
			
			String title = imts.length == 1 ? " " : gmpe.getShortName()+" Log-Normal Comparision";
			String xAxisLabel = "z-score (Standard Deviations)";
			String yAxisLabel = "Density";
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			spec.setLegendVisible(imt == imts[imts.length-1]);
			
			specs.add(spec);
		}
		
		DecimalFormat twoSigFig = new DecimalFormat("0.00");
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<imts.length; i++) {
			List<String> labels = new ArrayList<>(binDescriptions);
			labels.add(0, imts[i].getShortName());
			
			double yEach = maxY/8d;
			double x = -numStdDev + 0.2;
			double y = maxY - yEach*1.2;
			
			Font bigFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
			Font smallFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
			
			List<XYAnnotation> anns = new ArrayList<>();
			XYTextAnnotation meanAnn = new XYTextAnnotation(
					"Mean = "+twoSigFig.format(means.get(i)), -x, y);
			meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			meanAnn.setFont(bigFont);
			anns.add(meanAnn);
			XYTextAnnotation stdDevAnn = new XYTextAnnotation(
					"σ-fract = "+twoSigFig.format(stdDevs.get(i)), -x, y-yEach);
			stdDevAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			stdDevAnn.setFont(bigFont);
			anns.add(stdDevAnn);
			
			for (int j=0; j<labels.size(); j++) {
				String label = labels.get(j);
				XYTextAnnotation ann = new XYTextAnnotation(label, x, y);
				y -= yEach;
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
				if (j == 0) {
					ann.setFont(bigFont);
					yEach *= (double)smallFont.getSize()/(double)bigFont.getSize();
				} else {
					ann.setFont(smallFont);
				}
				anns.add(ann);
			}
			
			for (double sigma=Math.ceil(-numStdDev); sigma<=numStdDev; sigma++) {
				int s = (int)Math.round(Math.abs(sigma));
				String label;
				if (sigma < -0.1)
					label = "-"+s+" σ";
				else if (s == 0)
					label = s+" σ";
				else
					label = "+"+s+" σ";
				XYTextAnnotation ann = new XYTextAnnotation(label, sigma, maxY);
				ann.setTextAnchor(TextAnchor.TOP_CENTER);
				ann.setFont(bigFont);
				ann.setPaint(stdDevColor);
				anns.add(ann);
			}
			
			specs.get(i).setPlotAnnotations(anns);
			yRanges.add(new Range(0, maxY));
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(order);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(outputDir, prefix);
		int width, height;
		if (specs.size() > 1) {
			width = 800;
			height = 200 + 300*specs.size();
		} else {
			width = 600;
			height = 450;
		}
		gp.getChartPanel().setSize(width, height);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		return true;
	}
	
	private static double calcArea(HistogramFunction hist) {
		double area = 0d;
		for (Point2D pt : hist)
			area += hist.getDelta()*pt.getY();
		return area;
	}
	
	private static final Comparator<HistogramFunction> histComparator = new Comparator<HistogramFunction>() {

		@Override
		public int compare(HistogramFunction o1, HistogramFunction o2) {
			double a1 = calcArea(o1);
			double a2 = calcArea(o2);
			return -Double.compare(a1, a2);
		}
	};

}
