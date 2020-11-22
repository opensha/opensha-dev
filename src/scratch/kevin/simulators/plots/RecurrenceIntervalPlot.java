package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class RecurrenceIntervalPlot extends AbstractPlot {
	
	private Collection<SimulatorElement> elementsToInclude;
	private double[] minMags;
	private double[] prevTimes;
	private List<MinMaxAveTracker> riTrackers;
	private List<List<Double>> intervals;
	
	private String title = "Interevent Times";
	
	private double[] compVals;
	private String compName;
	
	public RecurrenceIntervalPlot(double... minMags) {
		this(null, minMags);
	}

	public RecurrenceIntervalPlot(Collection<SimulatorElement> elementsToInclude, double... minMags) {
		Preconditions.checkArgument(minMags.length > 0, "Must supply at least one min mag");
		this.elementsToInclude = elementsToInclude;
		this.minMags = minMags;
		prevTimes = new double[minMags.length];
		riTrackers = Lists.newArrayList();
		intervals = Lists.newArrayList();
		for (int i=0; i<minMags.length; i++) {
			prevTimes[i] = Double.NaN;
			riTrackers.add(new MinMaxAveTracker());
			intervals.add(new ArrayList<Double>());
		}
	}
	
	public void setComparison(double[] compVals, String compName) {
		Preconditions.checkArgument(compVals.length == minMags.length);
		this.compVals = compVals;
		this.compName = compName;
	}
	
	public void setPlotTitle(String title) {
		this.title = title;
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return elementsToInclude;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double time = e.getTimeInYears();
		for (int i=0; i<minMags.length; i++) {
			if (e.getMagnitude() < minMags[i])
				continue;
			
			if (!Double.isNaN(prevTimes[i])) {
				double interval = time - prevTimes[i];
				riTrackers.get(i).addValue(interval);
				intervals.get(i).add(interval);
			}
			prevTimes[i] = time;
		}
	}
	
	private static int preferred_num_bins = 25;
	
	private HistogramFunction buildLogHist(List<Double> intervals) {
		if (intervals.isEmpty())
			return null;
		double minNonZero = Double.POSITIVE_INFINITY;
		double max = 0d;
		for (double interval : intervals) {
			if (interval > 0) {
				minNonZero = Math.min(minNonZero, interval);
				max = Math.max(max, interval);
			}
		}
		if (minNonZero == max) {
			minNonZero = minNonZero*0.9;
			max = max*1.1;
		}
		minNonZero = Math.max(minNonZero, 1e-4);
		double diff = max - minNonZero;
		Preconditions.checkState(diff >= 0);
		
//		double calcDelta = diff/preferred_num_bins;
//		double delta;
//		if (calcDelta <= 2.5d)
//			delta = 1d;
//		else if (calcDelta < 7.5d)
//			delta = 5d;
//		else if (calcDelta < 25d)
//			delta = 10d;
//		else if (calcDelta < 75d)
//			delta = 50d;
//		else if (calcDelta < 250d)
//			delta = 100d;
//		else if (calcDelta < 750d)
//			delta = 500d;
//		else if (calcDelta < 2500d)
//			delta = 1000d;
//		else if (calcDelta < 7500d)
//			delta = 5000d;
//		else
//			delta = 10000d;
		
		
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(Math.log10(minNonZero), Math.log10(max), 0.1);
		double globalMin = hist.getMinX() - 0.5*hist.getDelta();
		for (double interval : intervals) {
			double logInterval = Math.log10(interval);
			if (logInterval < globalMin)
				continue;
			int ind = hist.getClosestXIndex(Math.log10(interval));
			hist.add(ind, 1d);
		}
		
		hist.normalizeBySumOfY_Vals();
		
		// now convert to density
		for (int i=0; i<hist.size(); i++) {
			double middle = hist.getX(i);
			double left = middle - 0.5*hist.getDelta();
			double right = middle + 0.5*hist.getDelta();
			double width = Math.pow(10, right) - Math.pow(10, left);
			hist.set(i, hist.getY(i)/width);
		}
		
		return hist;
	}

	@Override
	public void finalizePlot() throws IOException {
		for (int i=0; i<minMags.length; i++) {
			String magLabel = getCleanMagLabel(minMags[i]);
			String myPrefix = getOutputPrefix()+"_m"+magLabel;
			String myTitle = title+", M≥"+magLabel;
			
			MinMaxAveTracker tracker = riTrackers.get(i);
			
			HistogramFunction hist = buildLogHist(intervals.get(i));
			if (hist == null)
				continue;
			double mean = tracker.getAverage();
			
			double minNonZeroY = Double.POSITIVE_INFINITY;
			for (Point2D pt : hist)
				if (pt.getY() > 0)
					minNonZeroY = Math.min(minNonZeroY, pt.getY());
			Range yRange = calcEncompassingLog10Range(minNonZeroY, hist.getMaxY());
			double minY = yRange.getLowerBound();
			double maxY = yRange.getUpperBound();
			
			double minX = hist.getMinX()-0.5*hist.getDelta();
			double maxX = hist.getMaxX()+0.5*hist.getDelta();
			
			Range xRange = new Range(minX, maxX);
			
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			hist.setName("Simulated");
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, getSecondaryColor()));
			
			funcs.add(getLine("Mean="+yearDF.format(mean), Math.log10(mean), minY, Math.log10(mean), maxY));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getPrimaryColor()));
			
			if (compVals != null) {
				funcs.add(getLine(compName+" Mean="+yearDF.format(compVals[i]), Math.log10(compVals[i]), minY,
						Math.log10(compVals[i]), maxY));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getComparableColor()));
				
				minX = Math.min(minX, 0.9*compVals[i]);
				maxX = Math.max(maxX, 1.1*compVals[i]);
			}
			
			double poissonRate = (intervals.get(i).size()+1)/getCurrentDurationYears();
			EvenlyDiscretizedFunc poissonFunc = new EvenlyDiscretizedFunc(hist.getMinX(), hist.getMaxX(), hist.size());
			for (int j=0; j<poissonFunc.size(); j++) {
				double x = poissonFunc.getX(j);
				double binStart = Math.pow(10, x-hist.getDelta());
				double binEnd = Math.pow(10, x+hist.getDelta());
				double probRupBeforeEnd = 1d - Math.exp(-poissonRate*binEnd);
				double probRupBeforeStart = 1d - Math.exp(-poissonRate*binStart);
				double probInBin = probRupBeforeEnd - probRupBeforeStart;
				poissonFunc.set(j, probInBin/(binEnd-binStart)); // normalized to prob density
			}
			poissonFunc.setName("Poisson");
			funcs.add(poissonFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			
			SimpleRegression regression = new SimpleRegression();
			double maxPoiss = poissonFunc.getMaxY();
			for (Point2D pt : hist) {
				if (pt.getY() < maxPoiss)
					// only regress over power law section
					continue;
				regression.addData(pt.getX(), Math.log10(pt.getY()));
			}
			List<XYAnnotation> anns = null;
			if (regression.getN() > 0) {
				ArbitrarilyDiscretizedFunc fitAboveFunc = new ArbitrarilyDiscretizedFunc("Power-Law Fit");
				ArbitrarilyDiscretizedFunc fitBelowFunc = new ArbitrarilyDiscretizedFunc();
				double slope = regression.getSlope();
				double intercept = regression.getIntercept();
				for (int j=0; j<hist.size(); j++) {
					double x = hist.getX(j);
					double y = Math.pow(10, slope*x + intercept);
					if (y > maxPoiss)
						fitAboveFunc.set(x, y);
					else
						fitBelowFunc.set(x, y);
				}
				double xAtMaxPoiss = (Math.log10(maxPoiss) - intercept)/slope;
				fitAboveFunc.set(xAtMaxPoiss, maxPoiss);
				fitBelowFunc.set(xAtMaxPoiss, maxPoiss);
				funcs.add(fitAboveFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
				funcs.add(fitBelowFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));
				
				anns = new ArrayList<>();
				double midAboveX = 0.5*(fitAboveFunc.getMinX()+fitAboveFunc.getMaxX());
				double midAboveY = fitAboveFunc.getInterpolatedY(midAboveX);
				XYTextAnnotation ann = new XYTextAnnotation(" Slope = "+optionalDigitDF.format(slope), midAboveX, midAboveY);
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
				ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
				anns.add(ann);
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, myTitle, "Log10 Interevent Time (years)", "Probability Density (1/yr)");
			plot.setLegendVisible(true);
			plot.setPlotAnnotations(anns);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(plot, false, true, xRange, yRange);
			gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
			gp.saveAsTXT(new File(getOutputDir(), myPrefix+".txt").getAbsolutePath());
			gp.saveAsPNG(new File(getOutputDir(), myPrefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(getOutputDir(), myPrefix+".pdf").getAbsolutePath());
		}
	}

}
