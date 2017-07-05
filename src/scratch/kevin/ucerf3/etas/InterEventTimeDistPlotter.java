package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.kevin.simulators.SimAnalysisCatLoader;

public class InterEventTimeDistPlotter {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		
//		String simName = "UCERF3-ETAS";
//		File resultFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin");
//		List<List<ETAS_EqkRupture>> etasCatalogs = ETAS_CatalogIO.loadCatalogsBinary(resultFile);
//		List<SimulatorEvent> simCatalog = null;
		
		String simName = "RSQSim-ALLCAL2";
		List<List<ETAS_EqkRupture>> etasCatalogs = null;
		List<? extends SimulatorEvent> simCatalog = new SimAnalysisCatLoader(true, null, false).getEvents();
		
		List<SimulatorEvent> simCatalogPoisson = null;
		List<List<ETAS_EqkRupture>> etasCatalogsPoisson = null;
		
		boolean doPoisson = true;
		
		// in seconds
		double minTime = 0.1d; // 0.1 second
		double maxTime = 1d*60d*60d*24*365.25*1000; // 1000 years
		
		List<double[]> magBins = Lists.newArrayList();
		magBins.add(new double[] {5d, 6d});
		magBins.add(new double[] {6d, 7d});
		magBins.add(new double[] {7d, 10d});
		
		double logMinTime = Math.log10(minTime);
		double logMaxTime = Math.log10(maxTime);
		double logDelta = 0.05;
		
		if (doPoisson) {
			if (etasCatalogs == null) {
				// simulator
				simCatalogPoisson = Lists.newArrayList();
				double startSecs = simCatalog.get(0).getTime();
				double endSecs = simCatalog.get(simCatalog.size()-1).getTime();
				double durationSecs = endSecs - startSecs;
				for (SimulatorEvent e : simCatalog) {
					double timeSeconds = startSecs + Math.random()*(durationSecs);
					simCatalogPoisson.add(e.cloneNewTime(timeSeconds, e.getID()));
				}
				Collections.sort(simCatalogPoisson);
			} else {
				// ETAS
				Comparator<ETAS_EqkRupture> comp = new Comparator<ETAS_EqkRupture>() {

					@Override
					public int compare(ETAS_EqkRupture o1, ETAS_EqkRupture o2) {
						return new Long(o1.getOriginTime()).compareTo(o2.getOriginTime());
					}
				};
				etasCatalogsPoisson = Lists.newArrayList();
				for (List<ETAS_EqkRupture> catalog : etasCatalogs) {
					long startMillis = catalog.get(0).getOriginTime();
					long endMillis = catalog.get(catalog.size() - 1).getOriginTime();
					long durationMillis = endMillis - startMillis;
					
					List<ETAS_EqkRupture> catalogPoisson = Lists.newArrayList();
					
					for (ETAS_EqkRupture rup : catalog) {
						long time = startMillis + (long)(Math.random()*durationMillis);
						ETAS_EqkRupture newRup = new ETAS_EqkRupture((ObsEqkRupture)rup.clone());
						newRup.setOriginTime(time);
						catalogPoisson.add(newRup);
					}
					
					Collections.sort(catalogPoisson, comp);
					etasCatalogsPoisson.add(catalogPoisson);
				}
			}
		}
		
		for (double[] magBin : magBins) {
			double minMag = magBin[0];
			double maxMag = magBin[1];
			
			System.out.println("Processing mag range "+minMag+" to "+maxMag);
			
			List<Double> returnPeriods;
			List<Double> returnPeriodsPoisson = null;
			
			if (etasCatalogs == null) {
				// Simulators
				returnPeriods = calcReturnPeriodsSimulators(simCatalog, minMag, maxMag);
				if (simCatalogPoisson != null)
					returnPeriodsPoisson = calcReturnPeriodsSimulators(simCatalogPoisson, minMag, maxMag);
			} else {
				// ETAS
				returnPeriods = calcReturnPeriodsETAS(etasCatalogs, minMag, maxMag);
				if (etasCatalogsPoisson != null)
					returnPeriodsPoisson = calcReturnPeriodsETAS(etasCatalogsPoisson, minMag, maxMag);
			}
			
			ArbitrarilyDiscretizedFunc densityHist = calcDensity(returnPeriods, logMinTime, logMaxTime, logDelta);
			ArbitrarilyDiscretizedFunc densityHistPoisson = null;
			if (returnPeriodsPoisson != null)
				densityHistPoisson = calcDensity(returnPeriodsPoisson, logMinTime, logMaxTime, logDelta);
			
//			PoissonDistribution p = new PoissonDistribution(mean);
//			HistogramFunction poissonLogTimeHist = HistogramFunction.getEncompassingHistogram(logMinTime, logMaxTime, logDelta);
//			for (int i=0; i<poissonLogTimeHist.size(); i++) {
//				double logX = poissonLogTimeHist.getX(i);
//				double logBinStart = logX - logDelta;
//				double logBinEnd = logX + logDelta;
//				double binStart = Math.pow(10, logBinStart);
//				double binEnd = Math.pow(10, logBinEnd);
//				
//				double binProb = p.cum
//			}
			
			String title;
			String prefix;
			if (maxMag > 9) {
				title = simName+" M"+(int)minMag+"+ Interevent Distribution";
				prefix = "interevent_m"+(int)minMag+"+_"+simName;
			} else {
				title = simName+" M"+(int)minMag+"-"+(int)maxMag+" Interevent Distribution";
				prefix = "interevent_m"+(int)minMag+"_"+(int)maxMag+"_"+simName;
			}
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			if (densityHistPoisson != null) {
				densityHistPoisson.setName("Poisson");
				funcs.add(densityHistPoisson);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 2f, Color.BLUE.darker()));
			}
			
			densityHist.setName("Catalog");
			funcs.add(densityHist);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 2f, Color.RED.darker()));
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Interevent Time (s)", "Density");
			spec.setLegendVisible(densityHistPoisson != null);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(1e-1, 1e10, 1e-12, 1e-1);
			
			ETAS_MultiSimAnalysisTools.setFontSizes(gp);
			
			gp.drawGraphPanel(spec, true, true);
			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(new File(outputDir, prefix+"_density.png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, prefix+"_density.pdf").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, prefix+"_density.txt").getAbsolutePath());
		}
	}
	
	private static void normalizeUnityArea(DiscretizedFunc linearFunc, double logDelta) {
		double curArea = calcArea(linearFunc, logDelta);
		double areaScale = 1/curArea;
		linearFunc.scale(areaScale);
		double newArea = calcArea(linearFunc, logDelta);
		Preconditions.checkState((float)newArea == 1f, "Bad area after adjustment. Orig: %s, Mod: %s", curArea, newArea);
	}
	
	private static double calcArea(DiscretizedFunc linearFunc, double logDelta) {
		double curArea = 0d;
		for (Point2D pt : linearFunc) {
			double x = pt.getX();
			double binStart = Math.pow(10, Math.log10(x) - logDelta);
			double binEnd = Math.pow(10, Math.log10(x) + logDelta);
			double binWidth = binEnd - binStart;
			double binHeight = pt.getY();
			curArea += binHeight * binWidth;
		}
		return curArea;
	}
	
	private static List<Double> calcReturnPeriodsSimulators(List<? extends SimulatorEvent> simCatalog, double minMag, double maxMag) {
		List<Double> returnPeriods = Lists.newArrayList();
		double prevTime = Double.NEGATIVE_INFINITY;
		for (SimulatorEvent e : simCatalog) {
			if (e.getMagnitude() < minMag || e.getMagnitude() > maxMag)
				continue;
			double t = e.getTime(); // seconds
			if (prevTime > Double.NEGATIVE_INFINITY) {
				double diffSeconds = t - prevTime;
				returnPeriods.add(diffSeconds);
			}
			prevTime = t;
		}
		return returnPeriods;
	}
	
	private static List<Double> calcReturnPeriodsETAS(List<List<ETAS_EqkRupture>> etasCatalogs, double minMag, double maxMag) {
		List<Double> returnPeriods = Lists.newArrayList();
		for (List<ETAS_EqkRupture> catalog : etasCatalogs) {
			long prevTime = Long.MIN_VALUE;
			for (ETAS_EqkRupture rup : catalog) {
				if (rup.getMag() < minMag || rup.getMag() > maxMag)
					continue;
				long t = rup.getOriginTime(); // millis
				if (prevTime > Long.MIN_VALUE) {
					long diff = t - prevTime;
					double diffSeconds = (double)diff/1000d;
					returnPeriods.add(diffSeconds);
				}
				prevTime = t;
			}
		}
		return returnPeriods;
	}
	
	private static ArbitrarilyDiscretizedFunc calcDensity(List<Double> returnPeriods,
			double logMinTime, double logMaxTime, double logDelta) {
		HistogramFunction logTimeHist = HistogramFunction.getEncompassingHistogram(logMinTime, logMaxTime, logDelta);
		double firstBinStart = logTimeHist.getMinX() - 0.5*logDelta;
		double lastBinEnd = logTimeHist.getMaxX() + 0.5*logDelta;
		int numBelow = 0;
		int numAbove = 0;
		double total = 0d;
		for (double returnPeriod : returnPeriods) {
			total += returnPeriod;
			double logDiffSeconds = Math.log10(returnPeriod);
			if (logDiffSeconds < firstBinStart) {
				numBelow++;
			} else if (logDiffSeconds > lastBinEnd) {
				numAbove++;
			} else {
				logTimeHist.add(logDiffSeconds, 1d);
			}
		}
		
		System.out.println("New below: "+numBelow);
		System.out.println("New above: "+numAbove);
		
		double mean = total/(double)returnPeriods.size();
		System.out.println("Mean: "+mean);
		
		// now convert to linear space
		ArbitrarilyDiscretizedFunc linearTimeHist = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logTimeHist)
			if (pt.getY() > 0)
				linearTimeHist.set(Math.pow(10, pt.getX()), pt.getY());
		
		// now density
		ArbitrarilyDiscretizedFunc densityHist = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : linearTimeHist) {
			double density = pt.getY()/pt.getX();
			densityHist.set(pt.getX(), density);
		}
		
		normalizeUnityArea(densityHist, logDelta);
		
		return densityHist;
	}

}
