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
import scratch.kevin.simulators.RSQSimCatalog;

public class SubSeismoCharComparison {

	public static void main(String[] args) throws IOException {
		String prefix = "rsqsim_char_test";
		
		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_5892.instance(); prefix += "_wus";
//		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_2585_1MYR.instance(); prefix += "_ca";
		
		List<RSQSimEvent> events = catalog.loader().skipSlipsAndTimes().skipYears(20000).minMag(5d).load();
		
		boolean cumulative = true;
		double obsDuration = 93;
//		double obsDuration = 200;
		double timeDelta = obsDuration;
		
		catalog.setFractForInclusion(0d);
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		
		int numSubSects = mapper.getSubSections().size();

		System.out.println("Loaded "+events.size()+" events");
		
		IncrementalMagFreqDist totalMFD = FaultSysTools.initEmptyMFD(5.01, 8.51);
		
		// first calculate M7 rates from the long term catalog
		double[] longTermM7Rates = new double[numSubSects];
		for (RSQSimEvent e : events) {
			totalMFD.add(totalMFD.getClosestXIndex(e.getMagnitude()), 1d);
			if (e.getMagnitude() >= 7d) {
				if (!cumulative && e.getMagnitude() >= 7.1)
					continue;
				int totalCount = 0;
				List<List<SubSectionMapping>> allMappings = mapper.getAllSubSectionMappings(e);
				for (List<SubSectionMapping> mappings : allMappings)
					totalCount += mappings.size();
				double fractEach = 1d/(double)totalCount;
				for (List<SubSectionMapping> mappings : allMappings)
					for (SubSectionMapping mapping : mappings)
						longTermM7Rates[mapping.getSubSect().getSectionId()] += fractEach;
			}
		}
		double startTime = events.get(0).getTimeInYears();
		double lastEventTime = events.get(events.size()-1).getTimeInYears();
		double catLength = lastEventTime - startTime;
		MinMaxAveTracker m7RatesTrack = new MinMaxAveTracker();
		for (int s=0; s<numSubSects; s++) {
			longTermM7Rates[s] /= catLength;
			m7RatesTrack.addValue(longTermM7Rates[s]);
		}
		totalMFD.scale(1d/catLength);
		
		System.out.println("Total MFD:\n"+totalMFD);
		
		System.out.println("Catalog spans "+catLength+" years");
		System.out.println("M7 rates:\t"+m7RatesTrack);
		
		double totalM7 = cumulative ? totalMFD.getCumRate(totalMFD.getClosestXIndex(7.01)) : totalMFD.getY(totalMFD.getClosestXIndex(7.01));
		double totalM5 = cumulative ? totalMFD.getCumRate(0) : totalMFD.getY(0);
		double totalRatio = totalM7/totalM5;
		System.out.println("Total MFD ratio:\t"+(float)totalM7+" / "+(float)totalM5+" = "+(float)totalRatio);
		double minCharRatio = totalRatio + 0.01;
		double maxAnticharRatio = totalRatio - 0.01;
		System.out.println("Thresholds: <"+(float)maxAnticharRatio+", >"+(float)minCharRatio);
		
		Random rand = new Random(events.size());
		startTime += rand.nextDouble()*timeDelta;
		
		int startIndex = 0;
		int numWindows = 0;
		SummedMagFreqDist avgCharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
		SummedMagFreqDist avgAnticharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
		MinMaxAveTracker ratiosTrack = new MinMaxAveTracker();
		MinMaxAveTracker fractCharTrack = new MinMaxAveTracker();
		MinMaxAveTracker fractAnticharTrack = new MinMaxAveTracker();
		MinMaxAveTracker windowEventCountTrack = new MinMaxAveTracker();
		Range ratioRange = new Range(0d, Math.ceil(totalRatio*10d)/5d);
//		double ratioDelta = 0.001;
		double ratioDelta = 0.0005;
		EvenlyDiscretizedFunc ratioFunc = HistogramFunction.getEncompassingHistogram(ratioRange.getLowerBound(), ratioRange.getUpperBound(), ratioDelta);
		EvenlyDiscretizedFunc charRatioFunc = new EvenlyDiscretizedFunc(ratioFunc.getMinX(), ratioFunc.size(), ratioFunc.getDelta());
		EvenlyDiscretizedFunc anticharRatioFunc = new EvenlyDiscretizedFunc(ratioFunc.getMinX(), ratioFunc.size(), ratioFunc.getDelta());
		
		List<DiscretizedFunc> allSubCharMFDs = new ArrayList<>();
		List<DiscretizedFunc> allSubAnticharMFDs = new ArrayList<>();
		
		while (startTime + obsDuration < lastEventTime && startIndex < events.size()) {
			while (startIndex < events.size() && events.get(startIndex).getTimeInYears() < startTime)
				startIndex++;
			double endTime = startTime + obsDuration;
			Preconditions.checkState(endTime <= lastEventTime);
			List<RSQSimEvent> windowedEvents = new ArrayList<>();
			for (int i=startIndex; i<events.size(); i++) {
				RSQSimEvent event = events.get(i);
				Preconditions.checkState(event.getTimeInYears() >= startTime, "StartIndex problem");
				if (event.getTimeInYears() > endTime)
					break;
				windowedEvents.add(event);
			}
			Preconditions.checkState(!windowedEvents.isEmpty(), "No events in %s yr window %s starting in %s?",
					obsDuration, numWindows, startIndex);
			windowEventCountTrack.addValue(windowedEvents.size());
			
			// now figure out M5 rates in this window
			IncrementalMagFreqDist[] windowMFDs = new IncrementalMagFreqDist[numSubSects];
			for (int s=0; s<numSubSects; s++)
				windowMFDs[s] = new IncrementalMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
			for (RSQSimEvent e : windowedEvents) {
				int magIndex = totalMFD.getClosestXIndex(e.getMagnitude());
				int totalCount = 0;
				List<List<SubSectionMapping>> allMappings = mapper.getAllSubSectionMappings(e);
				for (List<SubSectionMapping> mappings : allMappings)
					totalCount += mappings.size();
				Preconditions.checkState(totalCount > 0, "No mappings for event?");
				double fractEach = 1d/(double)totalCount;
				for (List<SubSectionMapping> mappings : allMappings)
					for (SubSectionMapping mapping : mappings)
						windowMFDs[mapping.getSubSect().getSectionId()].add(magIndex, fractEach);
			}
			for (int s=0; s<numSubSects; s++)
				windowMFDs[s].scale(1d/obsDuration);
			
			SummedMagFreqDist charMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
			SummedMagFreqDist anticharMFD = new SummedMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta());
			
			int numChar = 0;
			int numAntichar = 0;
			MinMaxAveTracker m5RatesTrack = new MinMaxAveTracker();
			double[] ratios = new double[numSubSects];
			for (int s=0; s<numSubSects; s++) {
				IncrementalMagFreqDist  mfd = windowMFDs[s];
//				if (numWindows == 0 && s % 100 == 0)
//					System.out.println("1st window sect "+s+" MFD:\n"+mfd);
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
			fractCharTrack.addValue((double)numChar/(double)numSubSects);
			fractAnticharTrack.addValue((double)numAntichar/(double)numSubSects);
			
			avgCharMFD.addIncrementalMagFreqDist(charMFD);
			avgAnticharMFD.addIncrementalMagFreqDist(anticharMFD);
			
			allSubCharMFDs.add(charMFD);
			allSubAnticharMFDs.add(anticharMFD);
			
			numWindows++;
			startTime += timeDelta;
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
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "M7 rate / M5 rate", "Fraction");
//		plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
		
		gp.drawGraphPanel(plot, false, true, ratioRange, new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(new File("/tmp"), prefix+"_ratios", gp, 800, 800, true, false, false);
		
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
			
			PlotUtils.writePlots(new File("/tmp"), myPrefix, gp, 800, 800, true, false, false);
		}
	}

}
