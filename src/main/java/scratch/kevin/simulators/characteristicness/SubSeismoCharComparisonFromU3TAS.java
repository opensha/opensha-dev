package scratch.kevin.simulators.characteristicness;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.kevin.simulators.RSQSimCatalog;

public class SubSeismoCharComparisonFromU3TAS {

	public static void main(String[] args) throws IOException {
		String prefix = "u3tas_char_test";
		
		File simDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_05_24-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog");
		File addSimDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_06_06-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog");
		File outputDir = new File(simDir, "char_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		ETAS_Config config = ETAS_Config.readJSON(new File(simDir, "config.json"));
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		GriddedRegion gridReg = launcher.getRegion();
		
		List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(new File(simDir, "results_m5_preserve_chain.bin"), 5d);
		if (addSimDir != null) {
			catalogs = new ArrayList<>(catalogs);
			catalogs.addAll(ETAS_CatalogIO.loadCatalogsBinary(new File(addSimDir, "results_m5_preserve_chain.bin"), 5d));
		}
		
		boolean cumulative = true;
		double obsDuration = 93;
//		double obsDuration = 200;
		
		int numGridNodes = gridReg.getNodeCount();

		double durationEach = config.getDuration();
		double durationTotal = durationEach*catalogs.size();
		System.out.println("Loaded "+catalogs.size()+" catalogs with "+(float)durationTotal+" total years");
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(5.01, 8.51);
		SummedMagFreqDist totalMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		int m5Index = 0;
		int m7Index = refMFD.getClosestXIndex(7.01);
		
		// first calculate M7 rates from the long full simulation
		IncrementalMagFreqDist[] longTermMFDs = new IncrementalMagFreqDist[numGridNodes];
		for (int i=0; i<numGridNodes; i++)
			longTermMFDs[i] = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (ETAS_Catalog catalog : catalogs) {
			for (ETAS_EqkRupture rup : catalog) {
				double mag = rup.getMag();
				Location hypo = rup.getHypocenterLocation();
				Preconditions.checkNotNull(hypo);
				int node = gridReg.indexForLocation(hypo);
				if (node < 0)
					continue;
				longTermMFDs[node].add(refMFD.getClosestXIndex(mag), 1d);
			}
		}
		double[] longTermM5Rates = new double[numGridNodes];
		double[] longTermM7Rates = new double[numGridNodes];
		MinMaxAveTracker m7RatesTrack = new MinMaxAveTracker();
		for (int i=0; i<numGridNodes; i++) {
			longTermMFDs[i].scale(1d/durationEach);
			longTermM5Rates[i] = cumulative ? longTermMFDs[i].getCumRate(m5Index) : longTermMFDs[i].getY(m5Index);
			longTermM7Rates[i] = cumulative ? longTermMFDs[i].getCumRate(m7Index) : longTermMFDs[i].getY(m7Index);
			totalMFD.addIncrementalMagFreqDist(longTermMFDs[i]);
			m7RatesTrack.addValue(longTermM7Rates[i]);
		}
		totalMFD.scale(1d/catalogs.size());
		
		System.out.println("Total MFD:\n"+totalMFD);
		System.out.println("M7 rates:\t"+m7RatesTrack);
		
		double totalM7 = cumulative ? totalMFD.getCumRate(m7Index) : totalMFD.getY(m7Index);
		double totalM5 = cumulative ? totalMFD.getCumRate(m5Index) : totalMFD.getY(m5Index);
		double totalRatio = totalM7/totalM5;
		System.out.println("Total MFD ratio:\t"+(float)totalM7+" / "+(float)totalM5+" = "+(float)totalRatio);
		double minCharRatio = totalRatio + 0.001;
		double maxAnticharRatio = totalRatio - 0.001;
		System.out.println("Thresholds: <"+(float)maxAnticharRatio+", >"+(float)minCharRatio);
		
		Random rand = new Random((long)durationTotal);
		SummedMagFreqDist avgCharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
		SummedMagFreqDist avgAnticharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
		MinMaxAveTracker ratiosTrack = new MinMaxAveTracker();
		MinMaxAveTracker fractCharTrack = new MinMaxAveTracker();
		MinMaxAveTracker fractAnticharTrack = new MinMaxAveTracker();
		MinMaxAveTracker windowEventCountTrack = new MinMaxAveTracker();
		Range ratioRange = new Range(0d, Math.ceil(totalRatio*10d)/5d);
//		double ratioDelta = 0.001;
		double ratioDelta = 0.0005;
		EvenlyDiscretizedFunc ratioFuncLT = HistogramFunction.getEncompassingHistogram(ratioRange.getLowerBound(), ratioRange.getUpperBound(), ratioDelta);
		EvenlyDiscretizedFunc charRatioFuncLT = new EvenlyDiscretizedFunc(ratioFuncLT.getMinX(), ratioFuncLT.size(), ratioFuncLT.getDelta());
		EvenlyDiscretizedFunc anticharRatioFuncLT = new EvenlyDiscretizedFunc(ratioFuncLT.getMinX(), ratioFuncLT.size(), ratioFuncLT.getDelta());
		boolean[] charNodes = new boolean[numGridNodes];
		boolean[] anticharNodes = new boolean[numGridNodes];
		for (int i=0; i<numGridNodes; i++) {
			double ratio = longTermM7Rates[i]/longTermM5Rates[i];
			ratioFuncLT.add(ratioFuncLT.getClosestXIndex(ratio), 1d);
			ratiosTrack.addValue(ratio);
			if (ratio > minCharRatio) {
				charNodes[i] = true;
				charRatioFuncLT.add(ratioFuncLT.getClosestXIndex(ratio), 1d);
			} else if (ratio < maxAnticharRatio) {
				anticharNodes[i] = true;
				anticharRatioFuncLT.add(ratioFuncLT.getClosestXIndex(ratio), 1d);
			}
		}
		double numRatiosLT = ratioFuncLT.calcSumOfY_Vals();
		ratioFuncLT.scale(1d/numRatiosLT);
		charRatioFuncLT.scale(1d/numRatiosLT);
		anticharRatioFuncLT.scale(1d/numRatiosLT);
		EvenlyDiscretizedFunc ratioFunc = HistogramFunction.getEncompassingHistogram(ratioRange.getLowerBound(), ratioRange.getUpperBound(), ratioDelta);
		EvenlyDiscretizedFunc charRatioFunc = new EvenlyDiscretizedFunc(ratioFunc.getMinX(), ratioFunc.size(), ratioFunc.getDelta());
		EvenlyDiscretizedFunc anticharRatioFunc = new EvenlyDiscretizedFunc(ratioFunc.getMinX(), ratioFunc.size(), ratioFunc.getDelta());
		
		List<DiscretizedFunc> allSubCharMFDs = new ArrayList<>();
		List<DiscretizedFunc> allSubAnticharMFDs = new ArrayList<>();
		
		long obsDurationMillis = (long)(obsDuration * ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		long timeDeltaMillis = (long)(obsDuration * ProbabilityModelsCalc.MILLISEC_PER_YEAR * 0.1);

		int numWindows = 0;
		for (ETAS_Catalog catalog : catalogs) {
			int startIndex = 0;
			long startTimeMillis = config.getSimulationStartTimeMillis()
					+ (long)(obsDurationMillis*rand.nextDouble());
			long lastEventTimeMillis = startTimeMillis + (long)(config.getDuration()*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
			
			while (startTimeMillis + obsDurationMillis < lastEventTimeMillis && startIndex < catalog.size()) {
				while (startIndex < catalog.size() && catalog.get(startIndex).getOriginTime() < startTimeMillis)
					startIndex++;
				long endTimeMillis = startTimeMillis + obsDurationMillis;
				Preconditions.checkState(endTimeMillis <= lastEventTimeMillis);
				List<ETAS_EqkRupture> windowedEvents = new ArrayList<>();
				for (int i=startIndex; i<catalog.size(); i++) {
					ETAS_EqkRupture event = catalog.get(i);
					Preconditions.checkState(event.getOriginTime() >= startTimeMillis, "StartIndex problem");
					if (event.getOriginTime() > endTimeMillis)
						break;
					windowedEvents.add(event);
				}
				Preconditions.checkState(!windowedEvents.isEmpty(), "No events in %s yr window %s starting in %s?",
						obsDuration, numWindows, startIndex);
				windowEventCountTrack.addValue(windowedEvents.size());
				
				// now figure out M5 rates in this window
				IncrementalMagFreqDist[] windowMFDs = new IncrementalMagFreqDist[numGridNodes];
				for (int s=0; s<numGridNodes; s++)
					windowMFDs[s] = new IncrementalMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
				for (ETAS_EqkRupture e : windowedEvents) {
					int magIndex = totalMFD.getClosestXIndex(e.getMag());
					int node = gridReg.indexForLocation(e.getHypocenterLocation());
					if (node < 0)
						continue;
					windowMFDs[node].add(magIndex, 1d);
				}
				for (int s=0; s<numGridNodes; s++)
					windowMFDs[s].scale(1d/obsDuration);
				
				SummedMagFreqDist charMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
				SummedMagFreqDist anticharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
				
				int numChar = 0;
				int numAntichar = 0;
				MinMaxAveTracker m5RatesTrack = new MinMaxAveTracker();
				double[] ratios = new double[numGridNodes];
				for (int s=0; s<numGridNodes; s++) {
					IncrementalMagFreqDist  mfd = windowMFDs[s];
//					if (numWindows == 0 && s % 100 == 0)
//						System.out.println("1st window sect "+s+" MFD:\n"+mfd);
					double m5Rate = cumulative ? mfd.getTotalIncrRate() : mfd.getY(0);
					m5RatesTrack.addValue(m5Rate);
					double ratio = longTermM7Rates[s] / m5Rate;
					if (!Double.isFinite(ratio)) {
						ratios[s] = Double.NaN;
						continue;
					}
					ratios[s] = ratio;
					ratioFunc.add(ratioFunc.getClosestXIndex(ratio), 1d);
					ratiosTrack.addValue(ratio);
					if (ratio > minCharRatio) {
						numChar++;
						charMFD.addIncrementalMagFreqDist(mfd);
						charRatioFunc.add(ratioFunc.getClosestXIndex(ratio), 1d);
					} else if (ratio < maxAnticharRatio) {
						numAntichar++;
						anticharMFD.addIncrementalMagFreqDist(mfd);
						anticharRatioFunc.add(ratioFunc.getClosestXIndex(ratio), 1d);
					}
				}
				if (numWindows % 100 == 0) {
					System.out.println("Window "+numWindows+" M5 rates:\t"+m5RatesTrack);
					System.out.println("\t"+numChar+" characteristic");
					System.out.println("\t"+numAntichar+" anti-characteristic");
				}
				fractCharTrack.addValue((double)numChar/(double)numGridNodes);
				fractAnticharTrack.addValue((double)numAntichar/(double)numGridNodes);
				
				avgCharMFD.addIncrementalMagFreqDist(charMFD);
				avgAnticharMFD.addIncrementalMagFreqDist(anticharMFD);
				
				allSubCharMFDs.add(charMFD);
				allSubAnticharMFDs.add(anticharMFD);
				
				numWindows++;
				startTimeMillis += timeDeltaMillis;
			}
		}
		
		System.out.println("Processed "+numWindows+" time windows");
		System.out.println("Window event counts:\t"+windowEventCountTrack);
		System.out.println("Ratios:\t"+ratiosTrack);
		System.out.println("Fract chars:\t"+fractCharTrack);
		System.out.println("Fract antichars:\t"+fractAnticharTrack);
		
		avgCharMFD.scale(1d/numWindows);
		avgAnticharMFD.scale(1d/numWindows);
		double numRatios = ratioFunc.calcSumOfY_Vals();
		ratioFunc.scale(1d/numRatios);
		charRatioFunc.scale(1d/numRatios);
		anticharRatioFunc.scale(1d/numRatios);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		ratioFunc.setName("Total");
		funcs.add(ratioFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));

		charRatioFunc.setName("Most characteristic");
		funcs.add(charRatioFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_red));

		anticharRatioFunc.setName("Least characteristic");
		funcs.add(anticharRatioFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_blue));
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "M7 rate / Windowed M5 rate", "Fraction");
//		plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		gp.drawGraphPanel(plot, false, true, ratioRange, new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(outputDir, prefix+"_ratios", gp, 800, 800, true, false, false);
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();

		ratioFuncLT.setName("Total");
		funcs.add(ratioFuncLT);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));

		charRatioFuncLT.setName("Most characteristic");
		funcs.add(charRatioFuncLT);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_red));

		anticharRatioFuncLT.setName("Least characteristic");
		funcs.add(anticharRatioFuncLT);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_blue));
		
		plot = new PlotSpec(funcs, chars, " ", "M7 rate / M5 rate", "Fraction");
//		plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		gp.drawGraphPanel(plot, false, true, ratioRange, new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(outputDir, prefix+"_ratios_lt", gp, 800, 800, true, false, false);
		
		if (numWindows > 100) {
			Collections.shuffle(allSubCharMFDs, rand);
			Collections.shuffle(allSubAnticharMFDs, rand);
			allSubCharMFDs = allSubCharMFDs.subList(0, 100);
			allSubAnticharMFDs = allSubAnticharMFDs.subList(0, 100);
		}
		
		for (int i=0; i<3; i++) {
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			String myPrefix = prefix;
			if (i == 1) {
				// add individual chars
				myPrefix += "_char";
				for (DiscretizedFunc mfd : allSubCharMFDs) {
					mfd.setName(funcs.isEmpty() ? "Individual windows (most characteristic)" : null);
					funcs.add(mfd);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_lightred));
				}
			} else if (i == 2) {
				// add individual chars
				myPrefix += "_antichar";
				for (DiscretizedFunc mfd : allSubAnticharMFDs) {
					mfd.setName(funcs.isEmpty() ? "Individual windows (least characteristic)" : null);
					funcs.add(mfd);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_lightblue));
				}
			}
			
			totalMFD.setName("Total");
			funcs.add(totalMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
			
			avgCharMFD.setName("Most characteristic");
			funcs.add(avgCharMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_red));
			
			avgAnticharMFD.setName("Least characteristic");
			funcs.add(avgAnticharMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
			
			plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
//			plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			gp.drawGraphPanel(plot, false, true, new Range(5d,  8.5d), new Range(1e-6, 1e0));
			
			PlotUtils.writePlots(outputDir, myPrefix, gp, 800, 800, true, false, false);
		}
	}

}
