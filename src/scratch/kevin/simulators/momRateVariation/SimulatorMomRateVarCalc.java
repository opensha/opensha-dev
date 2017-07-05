package scratch.kevin.simulators.momRateVariation;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.EQSIM_EventRecord;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.utils.MatrixIO;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class SimulatorMomRateVarCalc {

	public static void main(String[] args) throws IOException {
		int[] windowLens = { 10, 25, 50, 75, 100, 150, 200 };
//		int[] windowLens = { 75 };
		
		// plot windows
//		for (int windowLen : windowLens) {
//			double[] taper = buildHanningTaper(windowLen);
//			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
//			for (int i=0; i<taper.length; i++)
//				func.set((double)i, taper[i]);
//			DefaultXY_DataSet vert = new DefaultXY_DataSet();
//			double x = (windowLen-1d)*0.5;
//			System.out.println("x="+x);
//			vert.set(x, 0);
//			vert.set(x, StatUtils.max(taper));
//			List<XY_DataSet> funcs = Lists.newArrayList();
//			funcs.add(func);
//			funcs.add(vert);
//			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
//			new GraphWindow(funcs, "Taper", chars);
//		}
		
//		File oDir = new File("/tmp/mom_rate");
//		Preconditions.checkState(oDir.exists() || oDir.mkdir());
		
//		List<RuptureIdentifier> idens = SynchIdens.getStandardSoCal();
		
		// don't use idens in loading, but rather
		List<RuptureIdentifier> loadIdens = Lists.newArrayList();
		// only so cal ruptures
		loadIdens.add(new RegionIden(new CaliforniaRegions.RELM_SOCAL()));
//		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, loadIdens, false).getEvents();
//		File mainDir = new File("/home/kevin/Simulators/UCERF3_interns");
////		String catName = "base";
////		String catName = "sigmahigh";
//		String catName = "sigmalow";
//		File dir = new File(mainDir, catName);
//		File dir = new File("/home/kevin/Simulators/UCERF3_JG_supraSeisGeo2");
		File dir = new File("/home/kevin/Simulators/UCERF3_interns/sigmahigh");
		File oDir = new File(dir, "mom_rate_time_series");
		Preconditions.checkState(oDir.exists() || oDir.mkdir());
		File geomFile = new File(dir, "UCERF3.D3.1.1km.tri.2.flt");
		
		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		loadIdens.add(new MagRangeRuptureIdentifier(5d, 10d));
		loadIdens = LogicalAndRupIden.buildAsList(loadIdens);
		List<? extends SimulatorEvent> events = RSQSimFileReader.readEventsFile(dir, elements, loadIdens);
		
		// randomize
		List<SimulatorEvent> trueRandom = Lists.newArrayList();
		double start = events.get(0).getTime();
		double end = events.get(events.size()-1).getTime();
		int id = 0;
		for (SimulatorEvent e : events)
			trueRandom.add(e.cloneNewTime((end-start)*Math.random()+start, id++));
		Collections.sort(trueRandom);
		
//		int windowLen = 50;
		int plotStartYear = 10000;
		int plotEndYear = 14000;
//		boolean hanningTaper = true;
		boolean[] hanningTapers = { true, false} ;
		boolean twoWay = false;
		
		for (boolean hanningTaper : hanningTapers) {
			String taperFileStr = "";
			if (hanningTaper)
				taperFileStr = "_taper";
			
			// plot regular
			plotMomRateVar(events, windowLens, "RSQSim Actual", plotStartYear, plotEndYear, hanningTaper, twoWay,
					new File(oDir, "actual"+taperFileStr+".png"));

			// now straight random
			plotMomRateVar(trueRandom, windowLens, "RSQSim Random", plotStartYear, plotEndYear, hanningTaper, twoWay,
					new File(oDir, "random"+taperFileStr+".png"));
		}
		
		// now write out time series for Matlab
		for (int windowLen : windowLens) {
			List<List<? extends SimulatorEvent>> eventLists = Lists.newArrayList();
			List<String> prefixes = Lists.newArrayList();
			eventLists.add(events);
			prefixes.add("actual");
			eventLists.add(trueRandom);
			prefixes.add("random");
			
			for (int n=0; n<eventLists.size(); n++) {
				List<? extends SimulatorEvent> myEvents = eventLists.get(n);
				String prefix = prefixes.get(n);
				
				System.out.println("Generating timeseries for "+prefix+", "+windowLen+" yr");
				
				File outputFile = new File(oDir, prefix+"_"+windowLen+"yr.bin");
				
				writeMomRateTimeSeries(windowLen, myEvents, outputFile);
			}
		}
		
//		List<EQSIM_Event> matchedEvents = Lists.newArrayList();
//		List<EQSIM_Event> unmatchedEvents = Lists.newArrayList();
//		eventLoop:
//		for (EQSIM_Event e : events) {
//			for (RuptureIdentifier iden : idens) {
//				if (iden.isMatch(e)) {
//					matchedEvents.add(e);
//					continue eventLoop;
//				}
//			}
//			// no match
//			unmatchedEvents.add(e);
//		}
//		
//		// plot only matched
//		plotMomRateVar2(matchedEvents, windowLen, "RSQSim Matched", plotStartYear, plotEndYear, logSmooth,
//				new File(oDir, "matched.png"));
//		// plot only unmatched
//		plotMomRateVar2(unmatchedEvents, windowLen, "RSQSim Unmatched", plotStartYear, plotEndYear, logSmooth,
//				new File(oDir, "unmatched.png"));
//		
//		// now randomize poisson
//		List<EQSIM_Event> poissonEvents = RandomCatalogBuilder.getRandomResampledCatalog(
//				events, idens, RandomDistType.POISSON, false, 1);
//		plotMomRateVar2(poissonEvents, windowLen, "RSQSim Poisson", plotStartYear, plotEndYear, logSmooth,
//				new File(oDir, "poisson.png"));
//		
//		// now unsynchronized
//		List<EQSIM_Event> unsynchEvents = RandomCatalogBuilder.getRandomResampledCatalog(
//				events, idens, RandomDistType.ACTUAL, false, 1);
//		plotMomRateVar2(unsynchEvents, windowLen, "RSQSim Unsynchronized", plotStartYear, plotEndYear, logSmooth,
//				new File(oDir, "unsynchronized.png"));
	}

	static void writeMomRateTimeSeries(int windowLen,
			List<? extends SimulatorEvent> myEvents, File outputFile) throws IOException {
		List<Double> yearsList = Lists.newArrayList();
		double startYear = myEvents.get(0).getTimeInYears() + windowLen;
		double endYear = myEvents.get(myEvents.size()-1).getTimeInYears() - windowLen;
		
		for (double year=startYear; year<endYear; year++)
			yearsList.add(year);
		
		double[] years = Doubles.toArray(yearsList);
		
		double[] momRates = calcTaperedMomRates(myEvents, years, buildHanningTaper(windowLen));
		
		MatrixIO.doubleArrayToFile(momRates, outputFile);
	}
	
	private static void plotMomRateVar(List<SimulatorEvent> events, int windowLen, String name,
			int plotStartYears, int plotEndYears, boolean logSmooth, File outputFile) throws IOException {
		int years = (int)General_EQSIM_Tools.getSimulationDurationYears(events);
		double minTime = events.get(0).getTimeInYears();
		
		double[] yearlyMoRates = new double[years];
		
		// seconds
		double startTime = events.get(0).getTime();
		
		for (SimulatorEvent e : events) {
			double secsFromStart = e.getTime()-startTime;
			int year = (int)(secsFromStart / General_EQSIM_Tools.SECONDS_PER_YEAR);
			if (year == years)
				break;
			double moment = 0;
			for (EventRecord r : e)
				moment += r.getMoment();
			yearlyMoRates[year] = yearlyMoRates[year] + moment;
		}
		
		double meanMoRate = StatUtils.mean(yearlyMoRates);
		double maxMoRate = StatUtils.max(yearlyMoRates);
		double minMoRate = StatUtils.min(yearlyMoRates);
		
		System.out.println("Long term: mean="+meanMoRate+"\tmax="+maxMoRate+"\tmin="+minMoRate);
		
		double[] windows = new double[years-windowLen];
		
		for (int i=0; (i+windowLen)<yearlyMoRates.length; i++) {
			double tot = 0;
			if (logSmooth)
				for (int j=i; j<i+windowLen; j++)
					tot += zeroIfNotFinite(Math.log(yearlyMoRates[j]));
			else
				for (int j=i; j<i+windowLen; j++)
					tot += yearlyMoRates[j];
			double avg = tot / (double)windowLen;
			windows[i] = avg;
			if (logSmooth)
				windows[i] = Math.exp(windows[i]);
		}
		
		System.out.println("Windows: mean="+StatUtils.mean(windows)+"\tmax="+StatUtils.max(windows)
				+"\tmin="+StatUtils.min(windows));
		
		System.out.println("Max window / mean = "+(StatUtils.max(windows) / meanMoRate));
		System.out.println("mean / min window = "+(meanMoRate / StatUtils.min(windows)));
		
//		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(
//				minTime+((double)windowLen*0.5), windows.length, 1d);
//		for (int i=0; i<windows.length; i++)
//			func.set(i, windows[i]);
		
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		for (int i=plotStartYears; i<plotEndYears; i++)
			func.set((double)i, windows[i]);
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		String title = name+" "+windowLen+"yr Moving Average of Seismic Moment Release";
		
		GraphWindow gw = new GraphWindow(funcs, title, chars);
		gw.setX_AxisLabel("Years");
		gw.setY_AxisLabel("Moment Rate (N-m/yr)");
		gw.setYLog(true);
		gw.setY_AxisRange(1e18, 1e20);
		gw.setPlotLabelFontSize(24);
		gw.setAxisLabelFontSize(18);
		gw.setTickLabelFontSize(16);
		gw.setSize(1200, 700);
		
		if (outputFile != null)
			gw.saveAsPNG(outputFile.getAbsolutePath());
	}
	
	static int findFirstEventIndexAfter(double timeSecs, List<? extends SimulatorEvent> events) {
		SimulatorEvent fakeEvent = new SimulatorEvent(new EQSIM_EventRecord(null));
		fakeEvent.setTime(timeSecs);
		int index = Collections.binarySearch(events, fakeEvent);
		if (index < 0) {
			// index = -ins_pt - 1
			// index + 1 = -ins_pt
			// ins_pt = -(index + 1)
			int ins_pt = -(index + 1);
			return ins_pt;
		}
		return index;
	}
	
	public static double calcWindowedMomentRate(List<? extends SimulatorEvent> events, int windowLen,
			double timeYears, double[] taper, boolean twoWay) {
		double windowSecs = General_EQSIM_Tools.SECONDS_PER_YEAR*(double)windowLen;
		double halfWindowSecs = windowSecs*0.5;
		
		double windowStart;
		if (twoWay)
			windowStart = timeYears*General_EQSIM_Tools.SECONDS_PER_YEAR-halfWindowSecs;
		else
			windowStart = timeYears*General_EQSIM_Tools.SECONDS_PER_YEAR-windowSecs;
		double windowEnd = windowStart + windowSecs;
		
		double taperWindowStart, taperWindowEnd;
		if (taper != null) {
			Preconditions.checkState(taper.length == windowLen);
			taperWindowStart = windowStart - halfWindowSecs;
			taperWindowEnd = windowEnd + halfWindowSecs;
		} else {
			taperWindowStart = windowStart;
			taperWindowEnd = windowEnd;
		}
		
		int index = findFirstEventIndexAfter(taperWindowStart, events);
		double moment = 0;
		for (int j=index; j<events.size(); j++) {
			SimulatorEvent e = events.get(j);
			double t = e.getTime();
			
			Preconditions.checkState(t >= taperWindowStart);
			if (t > taperWindowEnd)
				break;
			double myMoment = 0d;
			for (EventRecord r : e)
				myMoment += r.getMoment();
			if (taper != null) {
				// apply portion of taper within window
				for (int n=0; n<taper.length; n++) {
					// calculate time for this bin
					double relativeTime = (double)n/(double)(taper.length-1);
					double delta = windowSecs*relativeTime - halfWindowSecs;
					double absoluteTime = t + delta;
					if (absoluteTime >= windowStart && absoluteTime <= windowEnd)
						moment += myMoment*taper[n];
				}
			} else {
				moment += myMoment;
			}
		}
		
		moment /= (double)windowLen;
		
		return moment;
	}
	
	public static double[] calcTaperedMomRates(List<? extends SimulatorEvent> events, double[] years, double[] taper) {
		double startYear = years[0] - taper.length;
		
		int startIndex = findFirstEventIndexAfter(startYear, events);
		
		double[] momRates = new double[years.length];
		
		int taperHalfWidth = taper.length/2;
		
		for (int i=startIndex; i<events.size(); i++) {
			// find closest year
			SimulatorEvent e = events.get(i);
			
			double t = e.getTimeInYears();
			
			// center index of the taper
			int index = Arrays.binarySearch(years, t);
			if (index < 0)
				index = -(index + 1);
			// edge case treatment
			if (index == 0) {
				// it's before
				double delta = years[0] - t;
				Preconditions.checkState(delta >= 0, "Negative delta? "+delta);
				index = -(int)(delta + 0.5);
			} else if (index == years.length) {
				// after the end
				double delta = t - years[years.length - 1];
				Preconditions.checkState(delta > 0);
				index = years.length + (int)(delta + 0.5);
			}
			
			int taperStart = index - taperHalfWidth;
			if (taperStart >= years.length)
				break;
			
			double m = 0;
			for (EventRecord r : e)
				m += r.getMoment();
			
//			if (i == 1000) {
//				System.out.println("Debug for event "+i+" with moment "+m);
//				System.out.println("Event time: "+t);
//				System.out.println("Mapped to year "+years[index]+" at index "+index);
//			}
			
			for (int j=0; j<taper.length; j++) {
				int yearIndex = taperStart + j;
				if (yearIndex < 0)
					continue;
				if (yearIndex >= years.length)
					break;
				
				momRates[yearIndex] += m*taper[j];
//				if (i == 1000) {
//					System.out.println("\tAdding "+m*taper[j]+" moment to year "
//							+years[yearIndex]+" at index "+yearIndex);
//				}
			}
		}
		
		return momRates;
	}
	
	static void plotMomRateVar(List<? extends SimulatorEvent> events, int[] windowLens, String name,
			int plotStartYears, int plotEndYears, boolean hanningTaper, boolean twoWay, File outputFile)
					throws IOException {
//		int years = (int)General_EQSIM_Tools.getSimulationDurationYears(events);
		
		double secsPerYear = General_EQSIM_Tools.SECONDS_PER_YEAR;
		
		List<PlotSpec> specs = Lists.newArrayList();
		
		Range xRange = null;
		Range yRange = new Range(1e18, 1e20);
		
		for (int windowLen : windowLens) {
			double[] taper = null;
			if (hanningTaper)
				taper = buildHanningTaper(windowLen);
			
			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
			if (taper != null) {
				double[] years = new double[plotEndYears - plotStartYears];
				for (int i=0; i<years.length; i++)
					years[i] = plotStartYears+i;
				double[] momRates = calcTaperedMomRates(events, years, taper);
				for (int i=0; i<momRates.length; i++)
					func.set((double)i, momRates[i]);
			} else {
				for (int i=plotStartYears; i<plotEndYears; i++) {
					double moment = calcWindowedMomentRate(events, windowLen, (double)i, taper, twoWay);
					
					func.set((double)i, moment);
				}
			}
			
			if (xRange == null)
				xRange = new Range(func.getMinX(), func.getMaxX());
			
			ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(func);
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			
			String title = name+" Seismic Moment Release";
			if (hanningTaper)
				title += ", Tapered";
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Years", "Moment Rate (N-m/yr)");
			
			double x = xRange.getLowerBound()+(xRange.getUpperBound()-xRange.getLowerBound())*0.1;
			double y = yRange.getLowerBound()+(yRange.getUpperBound()-yRange.getLowerBound())*0.9;
			String label = windowLen+"yr";
			if (hanningTaper)
				label += " Hanning Taper";
			else
				label += " Moving Average";
			XYTextAnnotation ann = new XYTextAnnotation(label, x, y);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
			ann.setTextAnchor(TextAnchor.TOP_LEFT);
			List<XYTextAnnotation> anns = Lists.newArrayList(ann);
			spec.setPlotAnnotations(anns);
			
			specs.add(spec);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		List<Range> xRanges = Lists.newArrayList(xRange);
		List<Range> yRanges = Lists.newArrayList();
		for (int i=0; i<specs.size(); i++)
			yRanges.add(yRange);
		gp.drawGraphPanel(specs, false, true, xRanges, yRanges);
		gp.getChartPanel().setSize(1200, 60+260*windowLens.length);
		gp.saveAsPNG(outputFile.getAbsolutePath());
		
//		GraphWindow gw = new GraphWindow(funcs, title, chars);
//		gw.setX_AxisLabel("Years");
//		gw.setY_AxisLabel("Moment Rate (N-m/yr)");
//		gw.setYLog(true);
//		gw.setY_AxisRange(1e18, 1e20);
//		gw.setPlotLabelFontSize(24);
//		gw.setAxisLabelFontSize(18);
//		gw.setTickLabelFontSize(16);
//		gw.setSize(1200, 700);
//		
//		if (outputFile != null)
//			gw.saveAsPNG(outputFile.getAbsolutePath());
	}
	
	private static double zeroIfNotFinite(double val) {
		if (!Doubles.isFinite(val))
			return 0d;
		return val;
	}
	
	public static double[] buildHanningTaper(int windowLen) {
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
		
//		System.out.println("Window has sum="+StatUtils.sum(ret));
		
		return ret;
	}
	
	public static double[] getOnlyBeforeWindowTaper(int windowLenBefore) {
		double[] ret = new double[windowLenBefore*2+1];
//		for (int i=0; i<=windowLenBefore; i++)
		for (int i=windowLenBefore+1; i<ret.length; i++)
			ret[i] = 1d;
		
		double sum = StatUtils.sum(ret);
		for (int i=0; i<ret.length; i++)
			ret[i] /= sum;
		
		System.out.println(Joiner.on(",").join(Doubles.asList(ret)));
		
		return ret;
	}

}
