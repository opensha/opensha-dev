package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
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
	
	private HistogramFunction buildHist(MinMaxAveTracker tracker, List<Double> intervals) {
		if (tracker.getNum() == 0)
			return null;
		double min = tracker.getMin();
		double max = tracker.getMax();
		if (min == max) {
			min = min*0.9;
			max = max*1.1;
		}
		double diff = max - min;
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
		double mean = tracker.getAverage();
		double delta;
		if (mean < 25d)
			delta = 1d;
		else if (mean < 60d)
			delta = 5d;
		else if (mean < 250d)
			delta = 10d;
		else if (mean < 500d)
			delta = 20d;
		else if (mean < 1000d)
			delta = 50d;
		else if (mean < 2500d)
			delta = 100d;
		else if (mean < 5000d)
			delta = 200d;
		else if (mean < 10000d)
			delta = 500d;
		else
			delta = 1000d;
			
		
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(min, max, delta);
//		System.out.println(tracker);
//		System.out.println("input min="+min+", max="+max+", delta="+delta);
//		System.out.println("output min="+hist.getMinX()+", max="+hist.getMaxX());
		for (double interval : intervals)
			hist.add(interval, 1d);
		
		return hist;
	}

	@Override
	protected void finalize() throws IOException {
		for (int i=0; i<minMags.length; i++) {
			String magLabel = getCleanMagLabel(minMags[i]);
			String myPrefix = getOutputPrefix()+"_m"+magLabel;
			String myTitle = title+", M≥"+magLabel;
			
			MinMaxAveTracker tracker = riTrackers.get(i);
			
			HistogramFunction hist = buildHist(tracker, intervals.get(i));
			if (hist == null)
				continue;
			double mean = tracker.getAverage();
			
			double minY = 0;
			double maxY = hist.getMaxY()*1.2;
			
			double minX = hist.getMinX()-0.5*hist.getDelta();
			double maxX = hist.getMaxX()+0.5*hist.getDelta();
			
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			hist.setName(getCatalogName()+" Histogram");
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, getSecondaryColor()));
			
			funcs.add(getLine(getCatalogName()+" Mean="+yearDF.format(mean), mean, minY, mean, maxY));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getPrimaryColor()));
			
			if (compVals != null) {
				funcs.add(getLine(compName+" Mean="+yearDF.format(compVals[i]), compVals[i], minY, compVals[i], maxY));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getComparableColor()));
				
				minX = Math.min(minX, 0.9*compVals[i]);
				maxX = Math.max(maxX, 1.1*compVals[i]);
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, myTitle, "Interevent Time (years)", "Count");
			plot.setLegendVisible(true);
			
			Range xRange = new Range(minX, maxX);
			Range yRange = new Range(minY, maxY);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(plot, false, false, xRange, yRange);
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsTXT(new File(getOutputDir(), myPrefix+".txt").getAbsolutePath());
			gp.saveAsPNG(new File(getOutputDir(), myPrefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(getOutputDir(), myPrefix+".pdf").getAbsolutePath());
		}
	}

}
