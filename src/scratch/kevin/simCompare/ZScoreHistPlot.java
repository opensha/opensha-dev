package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
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
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public class ZScoreHistPlot {
	
	public static <E> boolean plotStandardNormal(SimulationRotDProvider<E> simProv, Collection<? extends RuptureComparison<E>> eventComps,
			List<Site> sites, double[] periods, AttenRelRef gmpe, RuptureComparisonFilter<E> filter, List<String> binDescriptions,
			File outputDir, String prefix) throws IOException {
		
		List<PlotSpec> specs = new ArrayList<>();
		double maxY = 0.7d;
		double numStdDev = 3.75;
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(new Range(-numStdDev, numStdDev));
		
		List<Double> means = new ArrayList<>();
		List<Double> stdDevs = new ArrayList<>();
		int numComputed = 0;
		int numMatches = 0;
		for (Site site : sites) {
			for (RuptureComparison<E> comp : eventComps) {
				if (!comp.isComputed(site, periods[0]))
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
		
		Color stdDevColor = new Color(0, 150, 0);
		
		for (double period : periods) {
			HistogramFunction hist = new HistogramFunction(-numStdDev, numStdDev, numBins);
			
			double mean = 0d;
			List<Double> allVals = new ArrayList<>();
			int count = 0;
			for (Site site : sites) {
				for (RuptureComparison<E> comp : eventComps) {
					if (!comp.hasSite(site))
						continue;
					Preconditions.checkState(comp.isComputed(site, period),
							"Computed for %ss but not %ss", periods[0], period);
					if ((filter != null && !filter.matches(comp, site)))
						continue;
					double gmpeVal = comp.getLogMean(site, period);
					for (DiscretizedFunc spectra : simProv.getRotD50s(site, comp.getRupture())) {
						// in log space
						double simVal = Math.log(spectra.getY(period));

						double val = (simVal - gmpeVal)/comp.getStdDev(site, period);
						hist.add(hist.getClosestXIndex(val), 1d);
						mean += val;
						allVals.add(val);
						count++;
					}
				}
			}
			if (count == 0)
				return false;
			mean /= count;
			means.add(mean);
			double stdDev = Math.sqrt(StatUtils.variance(Doubles.toArray(allVals), mean));
			stdDevs.add(stdDev);
			
			double area = 0d;
			for (Point2D pt : hist)
				area += hist.getDelta()*pt.getY();
			hist.scale(1d/area);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			EvenlyDiscretizedFunc stdNormal = new EvenlyDiscretizedFunc(hist.getMinX(), hist.getMaxX(), 1000);
			double scalar = 1d/Math.sqrt(2d*Math.PI);
			for (int i=0; i<stdNormal.size(); i++) {
				double x = stdNormal.getX(i);
				double y = scalar*Math.exp(-0.5*x*x);
				stdNormal.set(i, y);
			}
			
			funcs.add(hist);
			hist.setName(simProv.getName());
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			funcs.add(stdNormal);
			stdNormal.setName("Standard Normal");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
//			maxY = Math.max(maxY, Math.max(stdNormal.getMaxY(), hist.getMaxY()));
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			meanLine.set(mean, 0);
			meanLine.set(mean, maxY-0.1);
			meanLine.setName("Mean");
			funcs.add(meanLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
			
			for (double sigma=Math.ceil(-numStdDev); sigma<=numStdDev; sigma++) {
				DefaultXY_DataSet sigmaLine = new DefaultXY_DataSet();
				sigmaLine.set(sigma, 0);
				sigmaLine.set(sigma, maxY-0.075);
				funcs.add(sigmaLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, stdDevColor));
			}
			
			String title = gmpe.getShortName()+" Log-Normal Comparision";
			String xAxisLabel = "z-score (Standard Deviations)";
			String yAxisLabel = "Density";
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			spec.setLegendVisible(period == periods[periods.length-1]);
			
			specs.add(spec);
		}
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<periods.length; i++) {
			List<String> labels = new ArrayList<>(binDescriptions);
			labels.add(0, MultiRupGMPE_ComparePageGen.optionalDigitDF.format(periods[i])+"s SA");
			
			double yEach = maxY/8d;
			double x = -numStdDev + 0.2;
			double y = maxY - yEach*1.2;
			
			Font bigFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
			Font smallFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
			
			List<XYAnnotation> anns = new ArrayList<>();
			XYTextAnnotation meanAnn = new XYTextAnnotation(
					"Mean = "+MultiRupGMPE_ComparePageGen.optionalDigitDF.format(means.get(i)), -x, y);
			meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			meanAnn.setFont(bigFont);
			anns.add(meanAnn);
			XYTextAnnotation stdDevAnn = new XYTextAnnotation(
					"σ = "+MultiRupGMPE_ComparePageGen.optionalDigitDF.format(stdDevs.get(i)), -x, y-yEach);
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
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 200 + 300*specs.size());
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		return true;
	}

}
