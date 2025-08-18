package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public class StationarityPlot extends AbstractPlot {
	
	// for mag-vs-time
	private ArbitrarilyDiscretizedFunc[] timeCountFuncs;
	private double[] plotMags;
	
	// for MFDs
	private EvenlyDiscretizedFunc refMFD;
	private List<double[]> mfdsForTimes;
	
	private static final double storeDeltaYears = 1;
	private static final double rateEach = 1d/storeDeltaYears;
	private static final int numPlotBins = 10;
	
	public StationarityPlot() {
		this(5d, 8d);
	}
	
	public StationarityPlot(double minMag, double maxMag) {
		Preconditions.checkState(minMag <= maxMag);
		Preconditions.checkState(Double.isFinite(minMag));
		Preconditions.checkState(Double.isFinite(maxMag));
		List<Double> magsList = new ArrayList<>();
		magsList.add(minMag);
		if (minMag != Math.round(minMag))
			minMag = Math.ceil(minMag);
		else
			minMag += 1d;
		while (minMag <= maxMag) {
			magsList.add(minMag);
			minMag += 1;
		}
		Preconditions.checkState(!magsList.isEmpty());
		plotMags = Doubles.toArray(magsList);
		
		refMFD = FaultSysTools.initEmptyMFD(5.01, 9.41);
	}

	@Override
	protected synchronized void doProcessEvent(SimulatorEvent e) {
		double t = e.getTimeInYears();
		if (timeCountFuncs == null) {
			timeCountFuncs = new ArbitrarilyDiscretizedFunc[plotMags.length];
			for (int i = 0; i < timeCountFuncs.length; i++) {
				timeCountFuncs[i] = new ArbitrarilyDiscretizedFunc();
				timeCountFuncs[i].set(t, 0);
			}
			mfdsForTimes = new ArrayList<>();
			mfdsForTimes.add(new double[refMFD.size()]);
		}
		
		double mag = e.getMagnitude();
		if (mag < plotMags[0])
			return;
		
		double curMinTime = timeCountFuncs[0].getMaxX();
		double curMaxTime = curMinTime + storeDeltaYears;
		Preconditions.checkState(t >= curMinTime);
		
		while (t >= curMaxTime) {
			for (ArbitrarilyDiscretizedFunc func : timeCountFuncs)
				func.set(curMaxTime, 0);
			mfdsForTimes.add(new double[refMFD.size()]);
			curMinTime = curMaxTime;
			curMaxTime += storeDeltaYears;
		}
		Preconditions.checkState(mfdsForTimes.size() == timeCountFuncs[0].size());
		
		for (int i=0; i<plotMags.length; i++) {
			if (mag >= plotMags[i]) {
				int index = timeCountFuncs[i].size()-1;
				timeCountFuncs[i].set(index, timeCountFuncs[i].getY(index) + rateEach);
			}
		}
		
		int magIndex = refMFD.getClosestXIndex(mag);
		mfdsForTimes.get(mfdsForTimes.size()-1)[magIndex] += rateEach;
	}

	@Override
	public void finalizePlot() throws IOException {
		double minMag = plotMags[0];
		double maxMag;
		if (plotMags.length > 1)
			maxMag = plotMags[plotMags.length-1];
		else
			maxMag = minMag+1;
		CPT cpt = new CPT(minMag, maxMag, Color.BLUE, Color.RED);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		int bins;
		int numXeach;
		if (timeCountFuncs[0].size() < numPlotBins) {
			bins = timeCountFuncs[0].size();
			numXeach = 1;
		} else {
			bins = numPlotBins;
			numXeach = timeCountFuncs[0].size()/numPlotBins;
		}
		
		double minX = timeCountFuncs[0].getMinX();
		double binWidth = timeCountFuncs[0].getX(numXeach) - minX;
		double maxX = minX + bins*binWidth;
		String binString;
		if (binWidth < 10)
			binString = bins+" bins, "+(float)binWidth+" yr each";
		else
			binString = bins+" bins, "+(int)Math.round(binWidth)+" yr each";
		
		int[] startIndexes = new int[bins];
		double[] startTimes = new double[bins];
		int[] endIndexes = new int[bins];
		double[] endTimes = new double[bins];
		for (int i=0; i<bins; i++) {
			int start = i*numXeach;
			int end = start + numXeach;
			Preconditions.checkState(end <= timeCountFuncs[0].size());
			startIndexes[i] = start;
			startTimes[i] = timeCountFuncs[0].getX(start);
			endIndexes[i] = end;
			endTimes[i] = startTimes[i] + binWidth;
		}
		
		MinMaxAveTracker yTrack = new MinMaxAveTracker();
		for (int m=0; m<plotMags.length; m++) {
			ArbitrarilyDiscretizedFunc inputFunc = timeCountFuncs[m];
			if (inputFunc.getMaxY() == 0)
				continue;
			
			Color c = cpt.getColor((float)plotMags[m]);
			PlotCurveCharacterstics binnedChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, c);
			double meanForMag = 0;
			
			boolean first = true;
			
			for (int i=0; i<bins; i++) {
				int start = startIndexes[i];
				int end = endIndexes[i];
				Preconditions.checkState(end <= inputFunc.size());
				double mean = 0d;
				for (int index=start; index<end; index++)
					mean += inputFunc.getY(index);
				mean /= numXeach;
				if (mean > 0) {
					yTrack.addValue(mean);
					
					XY_DataSet binnedFunc = new DefaultXY_DataSet();
					if (first)
						binnedFunc.setName("M≥"+(float)plotMags[m]);
					first = false;
					binnedFunc.set(startTimes[i], mean);
					binnedFunc.set(endTimes[i], mean);
					
					funcs.add(binnedFunc);
					chars.add(binnedChar);
				}
				meanForMag += mean;
			}
			meanForMag /= bins;
			
			DefaultXY_DataSet straightLine = new DefaultXY_DataSet();
			straightLine.set(minX, meanForMag);
			straightLine.set(maxX, meanForMag);
			funcs.add(straightLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f,
					new Color(c.getRed(), c.getGreen(), c.getBlue(), 100)));
		}
		
		String title = getCatalogName()+" Stationarity";
		String xAxisLabel = "Years ("+binString+")";
		String yAxisLabel = "Annual Rate";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		Range yRange = calcEncompassingLog10Range(yTrack.getMin(), yTrack.getMax());
		Range xRange = new Range(minX, maxX);
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(getOutputDir(), getOutputPrefix(), gp, getPlotWidth(), (int)(getPlotHeight()*0.6), true, true, false);
		
		List<IncrementalMagFreqDist> mfds = new ArrayList<>();
		List<PlotCurveCharacterstics> mfdChars = new ArrayList<>();
		
		CPT startYearCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(startTimes[0], startTimes[startTimes.length-1]+1d);
		minMag = Double.POSITIVE_INFINITY;
		maxMag = 0d;
		double minRate = rateEach / (double)numXeach;
		double maxRate = 0d;
		for (int b=0; b<bins; b++) {
			IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			
			for (int i=startIndexes[b]; i<endIndexes[b]; i++) {
				double[] yValues = mfdsForTimes.get(i);
				for (int m=0; m<yValues.length; m++) {
					if (yValues[m] > 0) {
						mfd.add(m, yValues[m]);
						minMag = Math.min(minMag, mfd.getX(m));
						maxMag = Math.max(maxMag, mfd.getX(m));
					}
				}
			}
			
			mfd.scale(1d/(double)numXeach);
			
			maxRate = Math.max(maxRate, mfd.getMaxY());
			
			mfd.setName(null);
			mfds.add(mfd);
			mfdChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, startYearCPT.getColor(startTimes[b])));
		}
		
		minMag = 0.5*Math.floor(minMag*2d);
		maxMag = 0.5*Math.ceil(maxMag*2d);
		
		if (maxMag <= minMag) {
			minMag = 5d;
			maxMag = 9d;
		}
		
		if (minRate >= maxRate)
			maxRate = minRate * 10d;
		
		minRate = Math.pow(10, Math.floor(Math.log10(minRate)));
		maxRate = Math.pow(10, Math.ceil(Math.log10(maxRate)));
		
		plot = new PlotSpec(mfds, mfdChars, title, "Magnitude", "Incremental Rate (1/yr)");
		plot.setLegendVisible(true);
		PlotPreferences prefs = gp.getPlotPrefs();
		plot.addSubtitle(GraphPanel.getLegendForCPT(startYearCPT, "Start Year ("+binString+")",
				prefs.getAxisLabelFontSize(), prefs.getTickLabelFontSize(), Double.NaN, RectangleEdge.BOTTOM));
		
		xRange = new Range(minMag, maxMag);
		yRange = new Range(minRate, maxRate);
		
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(getOutputDir(), getOutputPrefix()+"_mfds", gp, getPlotWidth(), getPlotHeight(), true, true, false);
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
