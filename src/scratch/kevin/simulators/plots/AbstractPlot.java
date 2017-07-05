package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Random;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public abstract class AbstractPlot {
	
	private int eventCount = 0;
	private double minTime, maxTime;
	
	private static final int plot_width_default = 1000;
	private static final int plot_height_default = 800;
	
	private HeadlessGraphPanel gp;
	private XYZGraphPanel xyzGP;
	private String catalogName;
	private File outputDir;
	private String outputPrefix;
	
	public synchronized void processEvent(SimulatorEvent e) {
		double time = e.getTime();
		if (eventCount == 0) {
			minTime = time;
			maxTime = time;
		} else {
			Preconditions.checkState(time >= maxTime, "Events out of order! %s < %s", time, maxTime);
			maxTime = time;
		}
		eventCount++;
		doProcessEvent(e);
	}
	
	public void initialize(String catalogName, File outputDir, String outputPrefix) {
		this.catalogName = catalogName;
		this.outputDir = outputDir;
		this.outputPrefix = outputPrefix;
	}
	
	public String getCatalogName() {
		return catalogName;
	}

	public File getOutputDir() {
		return outputDir;
	}

	public String getOutputPrefix() {
		return outputPrefix;
	}

	protected abstract void doProcessEvent(SimulatorEvent e);
	
	protected abstract void finalize() throws IOException;
	
	/**
	 * This allows a plot to specify which elements it is applicable to for efficiency. If null, processEvent(e)
	 * will be called on all events. If non-null, processEvent(e) will only be called on events which include at least
	 * one of the supplied elements.
	 * @return null if this plot is applicable to all events, or a collection of elements for which
	 * this plot is applicable
	 */
	public abstract Collection<SimulatorElement> getApplicableElements();
	
	/**
	 * 
	 * @return current event count, including the current event if called during processing of an event
	 */
	protected int getCurrentEventCount() {
		return eventCount;
	}
	
	/**
	 * @return minimum event time in seconds
	 */
	protected double getMinTime() {
		return minTime;
	}
	
	/**
	 * @return current maximum event time processed, including current event if called during processing of an event 
	 */
	protected double getCurrentMaxTime() {
		return maxTime;
	}
	
	/**
	 * @return current duration in seconds, including current event if called during processing of an event
	 */
	protected double getCurrentDuration() {
		double duration = maxTime - minTime;
		Preconditions.checkState(duration >= 0, "Negative duration: %s", duration);
		return duration;
	}
	
	/**
	 * @return current duration in years, including current event if called during processing of an event
	 */
	protected double getCurrentDurationYears() {
		return getCurrentDuration()/General_EQSIM_Tools.SECONDS_PER_YEAR;
	}
	
	protected static Range calcEncompassingLog10Range(double min, double max) {
		Preconditions.checkState(min > 0, "Min must be positive for log plot! %s", min);
		Preconditions.checkState(min < max, "Min must be < max: %s >= %s", min, max);
		double logMin = Math.floor(Math.log10(min));
		double logMax = Math.ceil(Math.log10(max));
		
		return new Range(Math.pow(10, logMin), Math.pow(10, logMax));
	}
	
	protected static double minNonZero(XY_DataSet func) {
		return minNonZero(func, false);
	}
	
	protected static XY_DataSet downsampleByMag(XY_DataSet scatter, boolean magX, int targetNum) {
		double minMag, maxMag;
		if (magX) {
			minMag = scatter.getMinX();
			maxMag = scatter.getMaxX();
		} else {
			minMag = scatter.getMinY();
			maxMag = scatter.getMaxY();
		}
		EvenlyDiscretizedFunc magFunc = new EvenlyDiscretizedFunc(minMag, maxMag, 100);
		for (Point2D pt : scatter) {
			double mag = magX ? pt.getX() : pt.getY();
			magFunc.add(magFunc.getClosestXIndex(mag), 1);
		}
		double numEachBin = (double)targetNum / (double)magFunc.size();
		// now convert to prob of inclusion
		for (int i=0; i<magFunc.size(); i++) {
			double prob = numEachBin / magFunc.getY(i);
			magFunc.set(i, prob);
		}
		Random r = new Random();
		DefaultXY_DataSet sampled = new DefaultXY_DataSet();
		for (Point2D pt : scatter) {
			double mag;
			if (magX)
				mag = pt.getX();
			else
				mag = pt.getY();
			double prob = magFunc.getY(magFunc.getClosestXIndex(mag));
			Preconditions.checkState(prob > 0);
			if (r.nextDouble() < prob)
				sampled.set(pt);
		}
		return sampled;
	}
	
	protected static double minNonZero(XY_DataSet func, boolean x) {
		double min = Double.POSITIVE_INFINITY;
		for (Point2D pt : func) {
			double val = x ? pt.getX() : pt.getY();
			if (val > 0 && Doubles.isFinite(val))
				min = Math.min(val, min);
		}
		return min;
	}
	
	protected Color getPrimaryColor() {
		return Color.BLACK;
	}
	
	protected Color getSecondaryColor() {
		return Color.GRAY;
	}
	
	protected Color getComparableColor() {
		return Color.RED;
	}
	
	protected synchronized HeadlessGraphPanel getGraphPanel() {
		if (gp == null) {
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);
			gp = new HeadlessGraphPanel(plotPrefs);
		}
		
		return gp;
	}
	
	public synchronized XYZGraphPanel getXYZGraphPanel() {
		if (xyzGP == null)
			xyzGP = new XYZGraphPanel(getGraphPanel().getPlotPrefs());
		return xyzGP;
	}
	
	protected int getPlotHeight() {
		return plot_height_default;
	}
	
	protected int getPlotWidth() {
		return plot_width_default;
	}
	
	protected XY_DataSet getLine(String name, double x1, double y1, double x2, double y2) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		xy.setName(name);
		xy.set(x1, y1);
		xy.set(x2, y2);
		return xy;
	}
	
	protected String getCleanMagLabel(double mag) {
		if (mag == Math.floor(mag))
			return (int)mag+"";
		else return (float)mag+"";
	}
	
	protected static DecimalFormat yearDF = new DecimalFormat("0.00");

}
