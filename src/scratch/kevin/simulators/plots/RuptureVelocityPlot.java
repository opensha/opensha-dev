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

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
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
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RuptureVelocityPlot extends AbstractPlot {
	
	private LoadingCache<IDPairing, Double> elemDistCache;
	
	private Map<Integer, SimulatorElement> elemsMap;
	private double minMag;
	private List<Range> subMagRanges;
	
	private EvenlyDiscretizedFunc distFunc;
	private HistogramFunction[][] distBins;
	private double[][] binnedAverages;
	private double[][] binnedMins;
	private double[][] binnedMaxs;
	
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
		this(elemsMap, 200d, 20d, 0, 6, 0.5, minMag);
	}
	
	public RuptureVelocityPlot(Map<Integer, SimulatorElement> elemsMap, double maxDist, double deltaDist,
			double minVelocity, double maxVelocity, double deltaVelocity, double minMag) {
		this.elemsMap = elemsMap;
		this.minMag = minMag;
		
		subMagRanges = new ArrayList<>();
		subMagRanges.add(new Range(minMag, 10d));
		double lowerMag = Math.floor(minMag*2)/2;
		double deltaMag = 0.5;
		while ((float)lowerMag < (float)8.5) {
			subMagRanges.add(new Range(lowerMag, lowerMag + deltaMag));
			lowerMag += deltaMag;
		}
		// this makes the max bin unbounded
		if (subMagRanges.size() > 2) {
			int maxIndex = subMagRanges.size()-1;
			Range maxRange = subMagRanges.get(maxIndex);
			subMagRanges.set(maxIndex, new Range(maxRange.getLowerBound(), 10d));
		}
		
		CacheBuilder<Object, Object> build = CacheBuilder.newBuilder();
		build.concurrencyLevel(Integer.max(4, Runtime.getRuntime().availableProcessors()));
		elemDistCache = build.build(new Loader());
		
		distFunc = HistogramFunction.getEncompassingHistogram(0d, maxDist, deltaDist);
		distBins = new HistogramFunction[subMagRanges.size()][distFunc.size()];
		binnedAverages = new double[subMagRanges.size()][distFunc.size()];
		binnedMins = new double[subMagRanges.size()][distFunc.size()];
		binnedMaxs = new double[subMagRanges.size()][distFunc.size()];
		for (int m=0; m<subMagRanges.size(); m++) {
			for (int i=0; i<distBins[m].length; i++)
				distBins[m][i] = HistogramFunction.getEncompassingHistogram(minVelocity, maxVelocity, deltaVelocity);
			for (int i=0; i<binnedMins[m].length; i++)
				binnedMins[m][i] = Double.POSITIVE_INFINITY;
		}

		magVelFunc = new DefaultXY_DataSet();
		// single event scatter
		distVelFunc = new DefaultXY_DataSet();
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < minMag)
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
						for (int m=0; m<subMagRanges.size(); m++) {
							if (subMagRanges.get(m).contains(mag)) {
								int velBin = distBins[m][distBin].getClosestXIndex(vel);
								distBins[m][distBin].add(velBin, 1d);
								binnedAverages[m][distBin] += vel;
								binnedMins[m][distBin] = Math.min(binnedMins[m][distBin], vel);
								binnedMaxs[m][distBin] = Math.max(binnedMaxs[m][distBin], vel);
							}
						}
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
		magVelFunc.set(mag, aveVel);
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
	
	private List<DiscretizedFunc> calcDistScalingFuncs(ArbitrarilyDiscretizedFunc minusSigmaFunc, ArbitrarilyDiscretizedFunc plusSigmaFunc) {
		List<DiscretizedFunc> subMagFuncs = new ArrayList<>();
		for (int m=0; m<subMagRanges.size(); m++) {
			ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
			subMagFuncs.add(meanFunc);
			for (int i=0; i<distFunc.size(); i++) {
				double numInBin = distBins[m][i].calcSumOfY_Vals();
//				System.out.println("Mag "+m+" ("+subMagRanges.get(m).getLowerBound()+" "
//						+subMagRanges.get(m).getUpperBound()+"), dist "+i+" has "+numInBin);
				if (numInBin == 0)
					continue;
				binnedAverages[m][i]/=numInBin;
				
				double x = distFunc.getX(i);
				meanFunc.set(x, binnedAverages[m][i]);
				if (m == 0 && minusSigmaFunc != null && plusSigmaFunc != null) {
					double stdDev = 0;
					stdDev = distBins[m][i].computeStdDev();
					minusSigmaFunc.set(x, binnedAverages[m][i] - stdDev);
					plusSigmaFunc.set(x, binnedAverages[m][i] + stdDev);
				}
				
//				System.out.println(x+" "+binnedAverages[i]+" "+stdDev);
			}
		}
		return subMagFuncs;
	}
	
	private void plotDistScaling() throws IOException {
		ArbitrarilyDiscretizedFunc minusSigmaFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc plusSigmaFunc = new ArbitrarilyDiscretizedFunc();
		
		List<DiscretizedFunc> subMagFuncs = calcDistScalingFuncs(minusSigmaFunc, plusSigmaFunc);
		
		double maxY = Math.min(plot_max_vel, plusSigmaFunc.getMaxY()*1.4);
//		double maxY = Math.min(plot_max_vel, maxFunc.getMaxY()*1.2);
		
		UncertainArbDiscDataset stdDevData = new UncertainArbDiscDataset(subMagFuncs.get(0), minusSigmaFunc, plusSigmaFunc);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
//		minMaxData.setName("Min/Max Range");
//		funcs.add(minMaxData);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, Color.LIGHT_GRAY));
		
		stdDevData.setName("±σ Range");
		funcs.add(stdDevData);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 255, 0, 100)));
		
		subMagFuncs.get(0).setName("Mean");
		funcs.add(subMagFuncs.get(0));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		double maxMin = minMag;
		for (Range range : subMagRanges)
			maxMin = Math.max(range.getLowerBound(), maxMin);
		CPT subCPT = new CPT(minMag, maxMin, Color.LIGHT_GRAY, Color.DARK_GRAY);
		for (int m=1; m<subMagRanges.size() && distVelFunc == null; m++) {
			Range range = subMagRanges.get(m);
			Color c = subCPT.getColor((float)range.getLowerBound());
			DiscretizedFunc func = subMagFuncs.get(m);
			if (func.calcSumOfY_Vals() == 0d)
				continue;
			if (range.getUpperBound() > 9)
				func.setName("M"+optionalDigitDF.format(range.getLowerBound())+"+");
			else
				func.setName("M"+optionalDigitDF.format(range.getLowerBound())+"-"+optionalDigitDF.format(range.getUpperBound()));
			funcs.add(func);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, c));
		}
		
		String title = "Rupture Velocity vs Distance";
		String xAxisLabel = "Hypocentral Distance (km)";
		String yAxisLabel = "Rupture Velocity (km/s)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.setLegendFontSize(20);
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
	
	public static void plotDistanceScalingComparison(RuptureVelocityPlot plot1, String name1, RuptureVelocityPlot plot2, String name2,
			File outputDir, String prefix) throws IOException {
		// first mag range func is all mags
		DiscretizedFunc distFunc1 = plot1.calcDistScalingFuncs(null, null).get(0);
		DiscretizedFunc distFunc2 = plot2.calcDistScalingFuncs(null, null).get(0);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		distFunc1.setName(name1);
		funcs.add(distFunc1);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		distFunc2.setName(name2);
		funcs.add(distFunc2);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		String title = "Rupture Velocity vs Distance";
		String xAxisLabel = "Hypocentral Distance (km)";
		String yAxisLabel = "Mean Rupture Velocity (km/s)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(true);
		
		double maxY = Math.min(plot_max_vel, Math.max(distFunc1.getMaxY(), distFunc2.getMaxY())*1.2);
		
		HeadlessGraphPanel gp = plot1.getGraphPanel();
		gp.setLegendFontSize(20);
		gp.drawGraphPanel(plot, false, false, null, new Range(0, maxY));
		gp.getChartPanel().setSize(plot1.getPlotWidth(), (int)(plot1.getPlotHeight()*0.6));
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
		double minMag = 6.5;
		double skipYears = 5000;
		
		RSQSimCatalog catalog1 = Catalogs.BRUCE_2585.instance(baseDir);
		String name1 = "Original";
		RSQSimCatalog catalog2 = Catalogs.BRUCE_2740.instance(baseDir);
		String name2 = "Finite Receiver Modifications";
		
		File outputDir = new File("/tmp");
		String prefix = "rupture_velocity_comparison";
		
		RuptureVelocityPlot plot1 = new RuptureVelocityPlot(catalog1.getElements(), minMag);
		RuptureVelocityPlot plot2 = new RuptureVelocityPlot(catalog2.getElements(), minMag);

		for (RSQSimEvent e : catalog1.loader().minMag(minMag).skipYears(skipYears).load())
			plot1.processEvent(e);
		for (RSQSimEvent e : catalog2.loader().minMag(minMag).skipYears(skipYears).load())
			plot2.processEvent(e);
		
		plotDistanceScalingComparison(plot1, name1, plot2, name2, outputDir, prefix);
	}

}
