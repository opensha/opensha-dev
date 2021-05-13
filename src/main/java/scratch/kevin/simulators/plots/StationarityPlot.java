package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public class StationarityPlot extends AbstractPlot {
	
	private ArbitrarilyDiscretizedFunc[] timeCountFuncs;
	private double[] mags;
	
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
		mags = Doubles.toArray(magsList);
	}

	@Override
	protected synchronized void doProcessEvent(SimulatorEvent e) {
		double t = e.getTimeInYears();
		if (timeCountFuncs == null) {
			timeCountFuncs = new ArbitrarilyDiscretizedFunc[mags.length];
			for (int i = 0; i < timeCountFuncs.length; i++) {
				timeCountFuncs[i] = new ArbitrarilyDiscretizedFunc();
				timeCountFuncs[i].set(t, 0);
			}
		}
		
		double mag = e.getMagnitude();
		if (mag < mags[0])
			return;
		
		double curMinTime = timeCountFuncs[0].getMaxX();
		double curMaxTime = curMinTime + storeDeltaYears;
		Preconditions.checkState(t >= curMinTime);
		
		while (t >= curMaxTime) {
			for (ArbitrarilyDiscretizedFunc func : timeCountFuncs)
				func.set(curMaxTime, 0);
			curMinTime = curMaxTime;
			curMaxTime += storeDeltaYears;
		}
		
		for (int i=0; i<mags.length; i++) {
			if (mag >= mags[i]) {
				int index = timeCountFuncs[i].size()-1;
				timeCountFuncs[i].set(index, timeCountFuncs[i].getY(index) + rateEach);
			}
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		double minMag = mags[0];
		double maxMag;
		if (mags.length > 1)
			maxMag = mags[mags.length-1];
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
		
		MinMaxAveTracker yTrack = new MinMaxAveTracker();
		for (int m=0; m<mags.length; m++) {
			ArbitrarilyDiscretizedFunc inputFunc = timeCountFuncs[m];
			if (inputFunc.getMaxY() == 0)
				continue;
			
			Color c = cpt.getColor((float)mags[m]);
			PlotCurveCharacterstics binnedChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, c);
			double meanForMag = 0;
			
			boolean first = true;
			
			for (int i=0; i<bins; i++) {
				int start = i*numXeach;
				int end = start + numXeach;
				Preconditions.checkState(end <= inputFunc.size());
				double mean = 0d;
				for (int index=start; index<end; index++)
					mean += inputFunc.getY(index);
				mean /= numXeach;
				if (mean > 0) {
					yTrack.addValue(mean);
					
					XY_DataSet binnedFunc = new DefaultXY_DataSet();
					if (first)
						binnedFunc.setName("M≥"+(float)mags[m]);
					first = false;
					double startX = inputFunc.getX(start);
					double endX = startX + binWidth;
					binnedFunc.set(startX, mean);
					binnedFunc.set(endX, mean);
					
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
		gp.getChartPanel().setSize(getPlotWidth(), (int)(getPlotHeight()*0.6));
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
