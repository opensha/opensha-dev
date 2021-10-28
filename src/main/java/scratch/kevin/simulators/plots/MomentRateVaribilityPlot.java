package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import flanagan.math.FourierTransform;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class MomentRateVaribilityPlot extends AbstractPlot {

	private int windowLen;
	private double lenToPlot;
	private double halfTaper;
	
	private double[] taper;
	
	private List<Double> eventMoments;
	private List<Double> eventTimes;
	private double minTime = Double.POSITIVE_INFINITY;
	private double maxTime = Double.NEGATIVE_INFINITY;
	
	private static int welch_window_len = 4096;
	public static double MIN_CAT_LEN_YEARS = welch_window_len*3;

	public MomentRateVaribilityPlot(int windowLen, double lenToPlot) {
		this.windowLen = windowLen;
		this.lenToPlot = lenToPlot;
		this.halfTaper = windowLen * 0.5;
		
		eventMoments = new ArrayList<>();
		eventTimes = new ArrayList<>();
		
		taper = buildHanningTaper();
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double time = e.getTimeInYears();
		minTime = Math.min(time, minTime);
		maxTime = Math.max(time, maxTime);
		
		double moment = 0d;
		for (EventRecord r : e)
			moment += r.getMoment();
		
		eventMoments.add(moment);
		eventTimes.add(time);
	}
	
	private double[] calcMomRateTimeSeries(List<Double> eventMoments, List<Double> eventTimes) {
		double startYear = minTime - halfTaper;
		double endYear = maxTime + halfTaper;
		double lenDouble = endYear - startYear;
		int len = (int)Math.ceil(lenDouble);
		// nexpow2
		len = (int)Math.pow(2, Math.ceil(Math.log(len)/Math.log(2)));
		double[] data = new double[len];
		
		for (int i=0; i<eventMoments.size(); i++) {
			double moment = eventMoments.get(i);
			double time = eventTimes.get(i);
			
			double taperStart = time - halfTaper;
			for (int t=0; t<taper.length; t++) {
				double taperYear = taperStart + t;
				int index = indexForYear(startYear, taperYear);
				if (index < 0)
					continue;
				if (index >= len)
					break;
				data[index] += moment*taper[t];
			}
		}
		
		return data;
	}
	
	private static int indexForYear(double startYear, double year) {
		double indexDouble = (year - startYear);
		return (int)Math.round(indexDouble);
	}
	
	private ArbitrarilyDiscretizedFunc calcPSD(double[] data) {
		if (data.length < 2*welch_window_len) {
			System.out.println("Not enough data for moment variation calculation with only "+data.length+" years");
			return new ArbitrarilyDiscretizedFunc();
		}
		
		FourierTransform trans = new FourierTransform(data);
		trans.setSegmentLength(welch_window_len);
		int numSeg = data.length/welch_window_len - 1;
//		System.out.println("Num segments: "+numSeg);
		trans.setSegmentNumber(numSeg);
		trans.setWelch();
		trans.setOverlapOption(true);
		trans.setDeltaT(1d);
		double[][] psd = trans.powerSpectrum();
		
		ArbitrarilyDiscretizedFunc psdFunc = new ArbitrarilyDiscretizedFunc();
		
	
		for (int i=0; i<psd[0].length; i++) {
			double freq = psd[0][i];
			double power = psd[1][i];
			double period = 1d/freq;
			psdFunc.set(period, power);
		}
		
		return psdFunc;
	}
	
	public static final int num_poisson = 200;

	@Override
	public void finalizePlot() throws IOException {
		double[] data = calcMomRateTimeSeries(eventMoments, eventTimes);
		ArbitrarilyDiscretizedFunc dataPSD = calcPSD(data);
		double duration = maxTime - minTime;
		Random r = new Random(eventMoments.size());
		List<Double> randomTimes = new ArrayList<>();
		for (int i=0; i<eventTimes.size(); i++)
			randomTimes.add(r.nextDouble()*duration + minTime);
		double[] poissonData = calcMomRateTimeSeries(eventMoments, randomTimes);
		ArbitrarilyDiscretizedFunc poissonPSD = calcPSD(poissonData);
		
		// now poisson conf
		XY_DataSetList poissonFuncs = new XY_DataSetList();
		List<Double> poissonWeights = new ArrayList<>();
		poissonFuncs.add(poissonPSD);
		poissonWeights.add(1d);
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		System.out.print("Genering "+num_poisson+" random PSDs...");
		List<Future<ArbitrarilyDiscretizedFunc>> futures = new ArrayList<>();
		for (int i=1; i<num_poisson; i++) {
			final long seed = eventMoments.size()*i;
			futures.add(exec.submit(new Callable<ArbitrarilyDiscretizedFunc>() {
				
				@Override
				public ArbitrarilyDiscretizedFunc call() {
					Random r = new Random(seed);
					List<Double> myRandomTimes = new ArrayList<>();
					for (int j=0; j<eventTimes.size(); j++)
						myRandomTimes.add(r.nextDouble()*duration + minTime);
					double[] poissonData = calcMomRateTimeSeries(eventMoments, myRandomTimes);
					return calcPSD(poissonData);
				}
			}));
		}
		try {
			for (Future<ArbitrarilyDiscretizedFunc> future : futures) {
				ArbitrarilyDiscretizedFunc psd;
				try {
					psd = future.get();
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				if (psd.size() != poissonPSD.size())
					throw new IllegalStateException("Poisson func size mismatch. First is "+poissonPSD.size()+", encountered "+psd.size()+"."
							+" Origin x range: ["+poissonPSD.getMinX()+" "+poissonPSD.getMaxX()+"], "
									+ "current x range: ["+psd.getMinX()+" "+psd.getMaxX()+"]");
				poissonFuncs.add(psd);
				poissonWeights.add(1d);
			}
		} catch (RuntimeException e) {
			throw e;
		} finally {
			exec.shutdown();
		}
		System.out.println("DONE.");
		FractileCurveCalculator poissonFractiles = new FractileCurveCalculator(poissonFuncs, poissonWeights);
		ArbitrarilyDiscretizedFunc poissonLower = (ArbitrarilyDiscretizedFunc)poissonFractiles.getFractile(0.025);
		ArbitrarilyDiscretizedFunc poissonUpper = (ArbitrarilyDiscretizedFunc)poissonFractiles.getFractile(0.975);
		ArbitrarilyDiscretizedFunc poissonMedian = (ArbitrarilyDiscretizedFunc)poissonFractiles.getFractile(0.5);
		UncertainArbDiscFunc poissonRange = new UncertainArbDiscFunc(poissonMedian, poissonLower, poissonUpper);
		
		Range xRange = new Range(30, 3000);
		double minY = Double.POSITIVE_INFINITY, maxY = 0d;
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		dataPSD.setName(getCatalogName());
		funcs.add(dataPSD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		poissonPSD.setName("Poisson");
		funcs.add(poissonPSD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.GRAY));
		
		poissonRange.setName("95% Conf");
		funcs.add(poissonRange);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1.5f, Color.LIGHT_GRAY));
		
		for (XY_DataSet func : funcs) {
			for (Point2D pt : func) {
				if (xRange.contains(pt.getX())) {
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		
		Range yRange = calcEncompassingLog10Range(minY, maxY);
		yRange = new Range(Math.min(yRange.getLowerBound(), 1e35), Math.max(yRange.getUpperBound(), 1e38));
		
		
		PlotSpec psdSpec = new PlotSpec(funcs, chars, "Welch PSD of Tapered Moment Time Series", "Period (years)", "|Amplitude|");
		psdSpec.setLegendVisible(true);
		
		String prefix = getOutputPrefix()+"_welch";
		HeadlessGraphPanel gp = getGraphPanel();
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(psdSpec, true, true, xRange, yRange);
		gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		gp.saveAsTXT(new File(getOutputDir(), prefix+".txt").getAbsolutePath());
		gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
		
		if (lenToPlot > 0) {
			// now plot time series
			int midIndex = indexForYear(minTime, minTime + 0.5*(maxTime - minTime));
			int startIndex = midIndex + (int)Math.round(lenToPlot/2d);
			EvenlyDiscretizedFunc timeSeries = new EvenlyDiscretizedFunc(0d, (int)lenToPlot, 1d);
			for (int i=0; i<timeSeries.size(); i++)
				timeSeries.set(i, data[startIndex+i]);
			
			funcs.clear();
			chars.clear();
			
			funcs.add(timeSeries);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			double minNonZero = Double.POSITIVE_INFINITY;
			for (Point2D pt : timeSeries)
				if (pt.getY() > 0)
					minNonZero = Math.min(minNonZero, pt.getY());
			
			xRange = new Range(0, timeSeries.getMaxX());
			yRange = calcEncompassingLog10Range(minNonZero, timeSeries.getMaxY());
			
			PlotSpec timeSpec = new PlotSpec(funcs, chars, "Moment Release Time Series", "Time (yr)", "Moment (Nm)");
			timeSpec.setLegendVisible(true);
			
			prefix = getOutputPrefix()+"_time_series";
			gp = getGraphPanel();
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			gp.drawGraphPanel(timeSpec, false, true, xRange, yRange);
			gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight()*2/3);
			gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
		}
	}
	
	public double[] buildHanningTaper() {
		double[] ret = new double[windowLen];
		
		double twoPi = Math.PI*2;
		double innerCosMult = twoPi/((double)windowLen - 1d);
		
		for (int i=0; i<windowLen; i++)
			ret[i] = 0.5*(1d-Math.cos((double)i*innerCosMult));
		
		// now normalize
		double sum = 0;
		for (double v : ret)
			sum += v;
		for (int i=0; i<windowLen; i++)
			ret[i] /= sum;
		
		return ret;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
		double skipYears = 5000;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		MomentRateVaribilityPlot plot = new MomentRateVaribilityPlot(25, 2000);
		plot.initialize(catalog.getName(), outputDir, "mom_rate");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).load())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
			
//		RuptureVelocityPlot plot2 = new RuptureVelocityPlot(catalog2.getElements(), minMag);
//
//		for (RSQSimEvent e : catalog1.loader().minMag(minMag).skipYears(skipYears).load())
//			plot1.processEvent(e);
//		for (RSQSimEvent e : catalog2.loader().minMag(minMag).skipYears(skipYears).load())
//			plot2.processEvent(e);
//		
//		plotDistanceScalingComparison(plot1, name1, plot2, name2, outputDir, prefix);
	}

}
