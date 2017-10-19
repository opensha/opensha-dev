package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

import scratch.UCERF3.utils.IDPairing;

public class RuptureVelocityPlot extends AbstractPlot {
	
	private LoadingCache<IDPairing, Double> elemDistCache;
	
	private Map<Integer, SimulatorElement> elemsMap;
	private double minMag;
	
	private EvenlyDiscretizedFunc distFunc;
	private HistogramFunction[] distBins;
	private double[] binnedAverages;
	private double[] binnedMins;
	private double[] binnedMaxs;
	
	private XY_DataSet magVelFunc;
	private XY_DataSet distVelFunc;
	
	private static Map<Integer, SimulatorElement> buildMap(List<SimulatorElement> elems) {
		Map<Integer, SimulatorElement> elemsMap = new HashMap<>();
		for (SimulatorElement elem : elems)
			elemsMap.put(elem.getID(), elem);
		return elemsMap;
	}
	
	public RuptureVelocityPlot(List<SimulatorElement> elems, double minMag) {
		this(buildMap(elems), minMag);
	}
	
	public RuptureVelocityPlot(Map<Integer, SimulatorElement> elemsMap, double minMag) {
		this(elemsMap, minMag, 200d, 20d, 0, 6, 0.5);
	}
	
	public RuptureVelocityPlot(Map<Integer, SimulatorElement> elemsMap, double minMag,
			double maxDist, double deltaDist, double minVelocity, double maxVelocity, double deltaVelocity) {
		this.elemsMap = elemsMap;
		this.minMag = minMag;
		
		CacheBuilder<Object, Object> build = CacheBuilder.newBuilder();
		build.concurrencyLevel(Integer.max(4, Runtime.getRuntime().availableProcessors()));
		elemDistCache = build.build(new Loader());
		
		distFunc = HistogramFunction.getEncompassingHistogram(0d, maxDist, deltaDist);
		distBins = new HistogramFunction[distFunc.size()];
		for (int i=0; i<distBins.length; i++)
			distBins[i] = HistogramFunction.getEncompassingHistogram(minVelocity, maxVelocity, deltaVelocity);
		binnedAverages = new double[distBins.length];
		binnedMins = new double[distBins.length];
		for (int i=0; i<binnedMins.length; i++)
			binnedMins[i] = Double.POSITIVE_INFINITY;
		binnedMaxs = new double[distBins.length];
		
		magVelFunc = new DefaultXY_DataSet();
		distVelFunc = new DefaultXY_DataSet();
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		if (e.getMagnitude() < minMag)
			return;
		double minTime = Double.POSITIVE_INFINITY;
		int hypoID = -1;
		for (EventRecord rec : e) {
			double[] times = rec.getElementTimeFirstSlips();
			Preconditions.checkNotNull(times, "Event doesn't have timing information");
			int[] ids = rec.getElementIDs();
			
			for (int i=0; i<ids.length; i++) {
				if (times[i] < minTime) {
					minTime = times[i];
					hypoID = ids[i];
				}
			}
		}
		
		SimulatorElement hypo = elemsMap.get(hypoID);
		Preconditions.checkNotNull(hypo);
		
		double aveVel = 0d;
		int num = 0;
		
		if (distVelFunc != null && distVelFunc.size() > 0)
			// only for single events
			distVelFunc = null;
		
		for (EventRecord rec : e) {
			double[] times = rec.getElementTimeFirstSlips();
			int[] ids = rec.getElementIDs();
			
			for (int i=0; i<ids.length; i++) {
				if (ids[i] != hypoID) {
					try {
						double dist = elemDistCache.get(new IDPairing(hypoID, ids[i]));
						double tDelta = times[i] - minTime;
						if (tDelta == 0)
							continue;
						double vel = dist/(tDelta);
						Preconditions.checkState(Double.isFinite(vel) && vel > 0,
								"Bad velocity! vel = %s / %s = %s", dist, tDelta, vel);
						int distBin = distFunc.getClosestXIndex(dist);
						int velBin = distBins[distBin].getClosestXIndex(vel);
						distBins[distBin].add(velBin, 1d);
						binnedAverages[distBin] += vel;
						binnedMins[distBin] = Math.min(binnedMins[distBin], vel);
						binnedMaxs[distBin] = Math.max(binnedMaxs[distBin], vel);
						aveVel += vel;
						num++;
						if (distVelFunc != null)
							distVelFunc.set(dist, vel);
					} catch (ExecutionException e1) {
						throw ExceptionUtils.asRuntimeException(e1);
					}
				}
			}
		}
		
		if (num > 0)
			aveVel /= num;
		magVelFunc.set(e.getMagnitude(), aveVel);
	}

	@Override
	public void finalizePlot() throws IOException {
		if (distVelFunc == null)
			plotScatter(magVelFunc, true);
		else
			plotScatter(distVelFunc, false);
		plotDistScaling();
	}
	
	private static final int max_scatter_points = 100000;
	private static final double plot_max_vel = 8;
	
	private void plotScatter(XY_DataSet scatter, boolean mag) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		XY_DataSet plotScatter = scatter;
		if (mag && scatter.size() > max_scatter_points) {
			System.out.println("Filtering Mag/Vel scatter from "+scatter.size()+" to ~"+max_scatter_points+" points");
			plotScatter = downsampleByMag(scatter, true, max_scatter_points);
			System.out.println("Filter done (random mag-dependent sample): "+plotScatter.size());
		}
		funcs.add(plotScatter);
		plotScatter.setName(getCatalogName());
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		double minX = Math.floor(scatter.getMinX());
		double maxX = 0.5*Math.ceil(2d*scatter.getMaxX());
		Range xRange = new Range(minX, maxX);
		double minY = Math.floor(scatter.getMinY());
		double maxY = Math.min(plot_max_vel, Math.ceil(scatter.getMaxY()));
		Range yRange = new Range(minY, maxY);
		
		if (mag) {
			EvenlyDiscretizedFunc meanFunc = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), 20);
			int[] counts = new int[meanFunc.size()];
			double[] means = new double[meanFunc.size()];
			for (Point2D pt : scatter) {
				int xInd = meanFunc.getClosestXIndex(pt.getX());
				means[xInd] += pt.getY();
				counts[xInd]++;
			}
			for (int i=0; i<means.length; i++)
				meanFunc.set(i, means[i]/counts[i]);
			meanFunc.setName("Mean");
			funcs.add(meanFunc);
		} else {
			SimpleRegression regression = new SimpleRegression();
			for (Point2D pt : scatter)
				regression.addData(pt.getX(), pt.getY());
			double b = regression.getIntercept();
			double m = regression.getSlope();
			EvenlyDiscretizedFunc fit = new EvenlyDiscretizedFunc(minX, maxX, 10);
			// use one to one for x values
			for (int i=0; i<fit.size(); i++) {
				double x = fit.getX(i);
				double y = m*x + b;
				fit.set(x, y);
			}
			fit.setName("Best Fit Line");
			funcs.add(fit);
		}
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		String title = "Rupture Velocity Scatter";
		String xAxisLabel, yAxisLabel;
		if (mag) {
			xAxisLabel = "Magnitude";
			yAxisLabel = "Average Rupture Velocity (km/s)";
		} else {
			xAxisLabel = "Element Hypocentral Distance (km)";
			yAxisLabel = "Element Rupture Velocity (km/s)";
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.drawGraphPanel(plot, false, false, xRange, yRange);
		gp.getChartPanel().setSize(getPlotWidth(), (int)(getPlotHeight()*0.6));
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix()+"_scatter";
		if (!mag)
			prefix += "_dist";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private void plotDistScaling() throws IOException {
		double totalMean = 0d;
		ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
//		ArbitrarilyDiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
//		ArbitrarilyDiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc minusSigmaFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc plusSigmaFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<distBins.length; i++) {
			double numInBin = distBins[i].calcSumOfY_Vals();
			if (numInBin == 0)
				continue;
			totalMean += binnedAverages[i]/numInBin;
			binnedAverages[i]/=numInBin;
			double stdDev = 0;
			stdDev = distBins[i].computeStdDev();
			
			double x = distFunc.getX(i);
			meanFunc.set(x, binnedAverages[i]);
//			minFunc.set(x, binnedMins[i]);
//			maxFunc.set(x, binnedMaxs[i]);
			minusSigmaFunc.set(x, binnedAverages[i] - stdDev);
			plusSigmaFunc.set(x, binnedAverages[i] + stdDev);
			
//			System.out.println(x+" "+binnedAverages[i]+" "+stdDev);
		}
		
		double maxY = Math.min(plot_max_vel, plusSigmaFunc.getMaxY()*1.2);
//		double maxY = Math.min(plot_max_vel, maxFunc.getMaxY()*1.2);
		
		UncertainArbDiscDataset stdDevData = new UncertainArbDiscDataset(meanFunc, minusSigmaFunc, plusSigmaFunc);
//		UncertainArbDiscDataset minMaxData = new UncertainArbDiscDataset(meanFunc, minFunc, maxFunc);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
//		minMaxData.setName("Min/Max Range");
//		funcs.add(minMaxData);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, Color.LIGHT_GRAY));
		
		stdDevData.setName("±σ Range");
		funcs.add(stdDevData);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 255, 0, 100)));
		
		meanFunc.setName("Mean");
		funcs.add(meanFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		String title = "Rupture Velocity vs Distance";
		String xAxisLabel = "Hypocentral Distance (km)";
		String yAxisLabel = "Rupture Velocity (km/s)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.drawGraphPanel(plot, false, false, null, new Range(0, maxY));
		gp.getChartPanel().setSize(getPlotWidth(), (int)(getPlotHeight()*0.6));
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix()+"_vs_dist";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	private class Loader extends CacheLoader<IDPairing, Double> {

		@Override
		public Double load(IDPairing key) throws Exception {
			Location p1 = elemsMap.get(key.getID1()).getCenterLocation();
			Location p2 = elemsMap.get(key.getID2()).getCenterLocation();
			return LocationUtils.linearDistanceFast(p1, p2);
		}
		
	}

}
