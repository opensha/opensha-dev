package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class MFDPlot extends AbstractPlot {
	
	private IncrementalMagFreqDist mfd;
	private Collection<SimulatorElement> elementsToInclude;
	
	private static double min_mag_default = 4d;
	private static double max_mag_default = 9d;
	private static double delta_default = 0.1;
	private Range xRange;
	private Range yRange;
	
	private IncrementalMagFreqDist comparableMFD;
	private EvenlyDiscretizedFunc comparableCumulativeMFD;
	private String comparableName;
	
	private String plotTitle = "Magnitude Frequency Distribution";
	
	private boolean plotIncremental = true;
	private boolean plotCumulative = true;
	private boolean plotCombined = true;
	
	private boolean plotGR = true;
	private UncertainArbDiscFunc compRange;
	
	public MFDPlot() {
		this(min_mag_default);
	}
	
	public MFDPlot(Double minMag) {
		this(minMag, null);
	}
	
	public MFDPlot(Double minMag, Collection<SimulatorElement> elementsToInclude) {
		this(minMag, calcNum(minMag), delta_default, elementsToInclude);
	}
	
	private static int calcNum(Double minMag) {
		if (minMag == null)
			minMag = min_mag_default;
		return (int)((max_mag_default - minMag)/delta_default);
	}
	
	static IncrementalMagFreqDist buildIncrementalFunc(double minMag, double delta) {
		minMag = Math.floor(minMag/delta)*delta;
		int num = calcNum(minMag);
		return new IncrementalMagFreqDist(minMag+0.5*delta, num, delta);
	}
	
	public MFDPlot(Double minMag, int num, double delta, Collection<SimulatorElement> elementsToInclude) {
		if (minMag == null)
			minMag = min_mag_default;
		// shift by a half bin, cumulative will land on minMag
		mfd = new IncrementalMagFreqDist(minMag+0.5*delta, num, delta);
		this.elementsToInclude = elementsToInclude;
		xRange = new Range(minMag, mfd.getMaxX()+0.5*delta);
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		// if elementsToInclude is non null, then this will only be called if it's already a match
		double mag = e.getMagnitude();
		if (!xRange.contains(mag))
			return;
		int ind = mfd.getClosestXIndex(mag);
		mfd.set(ind, mfd.getY(ind)+1d);
	}
	
	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return elementsToInclude;
	}

	public void setComparableMFD(IncrementalMagFreqDist comparableMFD, String comparableName) {
		this.comparableMFD = comparableMFD;
		if (comparableMFD ==  null)
			comparableCumulativeMFD = null;
		else
			comparableCumulativeMFD = comparableMFD.getCumRateDistWithOffset();
		this.comparableName = comparableName;
	}
	
	public void setComparableRange(UncertainArbDiscFunc compRange) {
		this.compRange = compRange;
	}

	public void setPlotIncremental(boolean plotIncremental) {
		this.plotIncremental = plotIncremental;
	}

	public void setPlotCumulative(boolean plotCumulative) {
		this.plotCumulative = plotCumulative;
	}

	public void setPlotCombined(boolean plotCombined) {
		this.plotCombined = plotCombined;
	}
	
	public void setPlotTitle(String plotTitle) {
		this.plotTitle = plotTitle;
	}
	
	public void setPlotGR(boolean plotGR) {
		this.plotGR = plotGR;
	}

	@Override
	public void finalizePlot() throws IOException {
		mfd = mfd.deepClone();
		double durationYears = getCurrentDurationYears();
		Preconditions.checkState(durationYears >= 0d);
		// annualize
		mfd.scale(1d/durationYears);
		EvenlyDiscretizedFunc cumulative = mfd.getCumRateDistWithOffset();
		
		String prefix = getOutputPrefix();
		if (plotCombined) {
			makePlot(getOutputDir(), getCatalogName(), prefix, mfd, cumulative);
		} else {
			if (plotIncremental)
				makePlot(getOutputDir(), getCatalogName(), prefix+"_incremental", mfd, null);
			if (plotCumulative)
				makePlot(getOutputDir(), getCatalogName(), prefix+"_cumulative", null, cumulative);
		}
		
		// now write CSV
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Magnitude", "Incremental Rate (1/yr)", "Cumulative Rate (1/yr)");
		if (comparableMFD != null) {
			header.add(comparableName+" Incremental Rate (1/yr)");
			header.add(comparableName+" Cumulative Rate (1/yr)");
		}
		csv.addLine(header);
		
		for (int i=0; i<mfd.size(); i++) {
			// ues cumulative mags
			List<String> line = Lists.newArrayList((float)cumulative.getX(i)+"", mfd.getY(i)+"", cumulative.getY(i)+"");
			if (comparableMFD != null) {
				int index = comparableCumulativeMFD.getClosestXIndex(cumulative.getX(i));
				double closestX = comparableCumulativeMFD.getX(index);
				if (Math.abs(closestX-cumulative.getX(i)) < 0.51*comparableCumulativeMFD.getDelta()) {
					line.add(comparableMFD.getY(index)+"");
					line.add(comparableCumulativeMFD.getY(index)+"");
				} else {
					line.add("");
					line.add("");
				}
			}
			csv.addLine(line);
		}
		csv.writeToFile(new File(getOutputDir(), prefix+".csv"));
	}
	
	public void setYRange(Range yRange) {
		this.yRange = yRange;
	}
	
	private void makePlot(File outputDir, String catalogName, String prefix, EvenlyDiscretizedFunc mfd, EvenlyDiscretizedFunc cumulative)
			throws IOException {
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		boolean hasBoth = mfd != null && cumulative != null;
		
		if (comparableMFD != null) {
			if (hasBoth || cumulative != null) {
				// add cumulative plot
				if (hasBoth)
					comparableCumulativeMFD.setName(comparableName+" Cumulative");
				else
					comparableCumulativeMFD.setName(comparableName);
				funcs.add(comparableCumulativeMFD);
			} else {
				funcs.add(comparableMFD);
				comparableMFD.setName(comparableName);
			}
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getComparableColor()));
		}
		
		if (mfd != null) {
			if (hasBoth) {
				// histogram, on bottom
				funcs.add(0, mfd);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, getSecondaryColor()));
				mfd.setName(catalogName);
			} else {
				// line
				funcs.add(mfd);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getPrimaryColor()));
				mfd.setName(catalogName);
			}
		}
		
		if (cumulative != null) {
			funcs.add(cumulative);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getPrimaryColor()));
			if (hasBoth)
				cumulative.setName("Cumulative");
			else
				cumulative.setName(catalogName);
		}
		
		if (plotGR) {
			double totCmlRate;
			if (cumulative != null)
				totCmlRate = cumulative.getY(0);
			else
				totCmlRate = mfd.getY(0);
			GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(1d, 1d,
					xRange.getLowerBound(), xRange.getUpperBound(), 100);
			grMFD.scaleToIncrRate(0, totCmlRate);
			grMFD.setName("G-R B=1");
			funcs.add(grMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
		}
		
		if (compRange != null) {
			funcs.add(compRange);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getComparableColor()));

			funcs.add(compRange.getUpper());
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, getComparableColor()));

			funcs.add(compRange.getLower());
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, getComparableColor()));
		}
		
		String xAxisLabel = "Magnitude";
		String yAxisLabel;
		if (hasBoth)
			yAxisLabel = "Annual Rate";
		else if (mfd != null)
			yAxisLabel = "Incremental Rate (1/yr)";
		else
			yAxisLabel = "Cumulative Rate (1/yr)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, plotTitle, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		for (DiscretizedFunc func : funcs) {
			minY = Math.min(minY, minNonZero(func));
			maxY = Math.max(maxY, func.getMaxY());
		}
		Range yRange = this.yRange;
		if (yRange == null) {
			if (!Doubles.isFinite(minY))
				yRange = new Range(1d, 10d);
			else
				yRange = calcEncompassingLog10Range(minY, maxY);
		}
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
	}
	
	public static void plotMultiMFDs(Collection<DiscretizedFunc> mfds, DiscretizedFunc baselineMFD, File outputDir, String prefix)
			throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (baselineMFD != null) {
			funcs.add(baselineMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		}
		
		for (DiscretizedFunc mfd : mfds) {
			if (mfd == baselineMFD)
				continue;
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		}
		
		if (baselineMFD != null && baselineMFD.getName() != null)
			funcs.get(funcs.size()-1).setName((funcs.size()-1)+" variations");
		
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Cumulative Rate (1/yr)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, "MFD Comparison", xAxisLabel, yAxisLabel);
		plot.setLegendVisible(baselineMFD != null && baselineMFD.getName() != null);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		double minX = Double.POSITIVE_INFINITY;
		double maxX = 0d;
		for (DiscretizedFunc func : funcs) {
			minY = Math.min(minY, AbstractPlot.minNonZero(func));
			maxY = Math.max(maxY, func.getMaxY());
			minX = Math.min(minX, func.getMinX());
			for (Point2D pt : func)
				if (pt.getY() > 0)
					maxX = Math.max(maxX, pt.getX());
		}
		Range yRange;
		if (!Doubles.isFinite(minY))
			yRange = new Range(1d, 10d);
		else
			yRange = AbstractPlot.calcEncompassingLog10Range(minY, maxY);
		Range xRange = new Range(minX, maxX+0.1);
		
		HeadlessGraphPanel gp = buildGraphPanel();
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		gp.getChartPanel().setSize(650, 600);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4665.instance();
		double minMag = 5d;
		
		int[] skipYears = { 5000, 5000, 315000 };
		double[] maxDurs = { 0d, 310000d, 0d };
		
		File outputDir = new File(catalog.getCatalogDir(), "sub_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (int i=0; i<skipYears.length; i++) {
			String prefix = "mfd";
			if (skipYears[i] < 50000 && maxDurs[i] == 0d)
				prefix += "_full";
			else if (skipYears[i] < 50000 && maxDurs[i] > 0d)
				prefix += "_first_half";
			else if (skipYears[i] > 50000)
				prefix += "_second_half";
			Loader loader = catalog.loader().skipYears(skipYears[i]);
			if (maxDurs[i] > 0)
				loader.maxDuration(maxDurs[i]);
			System.out.println("Processing: "+prefix);
			MFDPlot plot = new MFDPlot(minMag);
			plot.setYRange(new Range(1e-6, 1e2));
			plot.initialize(catalog.getName(), outputDir, prefix);
			for (RSQSimEvent e : loader.iterable())
				plot.processEvent(e);
			plot.finalizePlot();
		}
	}

}
