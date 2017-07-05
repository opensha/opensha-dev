package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.EventsInWindowsMatcher;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.inversion.CommandLineInversionRunner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class BatchPlotGen {
	
	static final double DAYS_PER_YEAR = 365.242;
	private static final DecimalFormat df_hundreths = new DecimalFormat("0.00");
	private static final boolean PLOT_MFD_FAULT_COMPONENTS = false;

	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File plotBaseDir = new File(dir, "plots_short");
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File plotBaseDir = new File(dir, "plots_long");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		List<double[]> magRanges= Lists.newArrayList();
		magRanges.add(toArray(7, 10));
		magRanges.add(toArray(7, 7.5));
		magRanges.add(toArray(7.5, 10));
		
		Map<Integer, String> elementsOfInterest = Maps.newHashMap();
		elementsOfInterest.put(ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, "SAF Mojave");
		elementsOfInterest.put(ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, "SAF Carrizo");
		elementsOfInterest.put(ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, "SAF Coachella");
		elementsOfInterest.put(ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, "San Jacinto");
		elementsOfInterest.put(ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, "Garlock West");
		
		Region reg = new CaliforniaRegions.RELM_SOCAL();
		HashSet<Integer> elementsInRegion = MFDCalc.getElementsInsideRegion(tools.getElementsList(), reg);
		
		// These are special faults whose MFDs we want to plot separately
		List<List<Integer>> mfdSpecialFaultParents = Lists.newArrayList();
		List<String> mfdSpecialFaultNames = Lists.newArrayList();
		List<Color> mfdSpecialFaultColors = Lists.newArrayList();
		List<HashSet<Integer>> mfdSpecialFaultSections = Lists.newArrayList();
		
		if (PLOT_MFD_FAULT_COMPONENTS) {
			mfdSpecialFaultNames.add("San Andreas");
			List<Integer> safParents = Lists.newArrayList();
			// sections 1-18
			for (int i=1; i<=18; i++)
				safParents.add(i);
			mfdSpecialFaultParents.add(safParents);
			mfdSpecialFaultColors.add(Color.GREEN);

			mfdSpecialFaultNames.add("San Jacinto");
			List<Integer> sjParents = Lists.newArrayList();
			// sections 23-28
			for (int i=23; i<=28; i++)
				sjParents.add(i);
			mfdSpecialFaultParents.add(sjParents);
			mfdSpecialFaultColors.add(Color.MAGENTA);

			mfdSpecialFaultNames.add("Garlock");
			// sections 71, 72
			mfdSpecialFaultParents.add(Lists.newArrayList(71, 72));
			mfdSpecialFaultColors.add(Color.CYAN);

			for (List<Integer> parents : mfdSpecialFaultParents) {
				HashSet<Integer> sections = new HashSet<Integer>();
				for (SimulatorElement elem : tools.getElementsList()) {
					if (parents.contains(elem.getSectionID()) && elementsInRegion.contains(elem.getID()))
						sections.add(elem.getID());
				}
				mfdSpecialFaultSections.add(sections);
			}

			// now add other
			HashSet<Integer> otherElems = new HashSet<Integer>(elementsInRegion);
			for (HashSet<Integer> elems : mfdSpecialFaultSections)
				otherElems.removeAll(elems);

			mfdSpecialFaultNames.add("Other So Cal");
			mfdSpecialFaultColors.add(Color.GRAY);
			mfdSpecialFaultSections.add(otherElems);
		}
		
		double mfdPlotMinMag = 5.5;
		double mfdPlotMaxMag = 8.0;
		double delta = 0.1;
		int mfdPlotNumMag = (int)((mfdPlotMaxMag - mfdPlotMinMag) / delta) + 1;
		
		double[] durations = { 1d/DAYS_PER_YEAR, 7d/DAYS_PER_YEAR, 1d, 5d, 10d, 29d, 30d};
		double minWindowDurationYears = 0d;
		boolean randomizeEventTimes = false;
		
		List<double[]> particRanges = Lists.newArrayList();
		particRanges.add(toArray(6.5, 7));
		particRanges.add(toArray(7, 10));
		particRanges.add(toArray(6.5, 10));
		
		Map<Integer, List<Integer>> elemsByFaultMap = Maps.newHashMap();
		Map<Integer, String> faultNameMap = Maps.newHashMap();
		for (SimulatorElement element : tools.getElementsList()) {
			Integer faultID = element.getSectionID();
			List<Integer> elemsByFault = elemsByFaultMap.get(faultID);
			if (elemsByFault == null) {
				elemsByFault = Lists.newArrayList();
				elemsByFaultMap.put(faultID, elemsByFault);
				faultNameMap.put(faultID, element.getName());
			}
			elemsByFault.add(element.getID());
		}
		
		MinMaxAveTracker eventMagTrack = new MinMaxAveTracker();
		for (SimulatorEvent e : events)
			eventMagTrack.addValue(e.getMagnitude());
		double totalEventDuration = General_EQSIM_Tools.getSimulationDurationYears(events);
		
		List<Map<Integer, Double>> indepParticRatesList = Lists.newArrayList();
		for (double[] particRange : particRanges)
			indepParticRatesList.add(calcParticRates(events, totalEventDuration, particRange[0], particRange[1]));
		
		List<Map<Integer, double[]>> indepFaultParticRatesList = Lists.newArrayList();
		for (int i=0; i<indepParticRatesList.size(); i++) {
			Map<Integer, Double> indepParticRates = indepParticRatesList.get(i);
			Map<Integer, double[]> indepFaultParticRates = Maps.newHashMap();
			indepFaultParticRatesList.add(indepFaultParticRates);
			for (int faultID : elemsByFaultMap.keySet()) {
				List<Integer> elemIDs = elemsByFaultMap.get(faultID);
				double[] faultParticRates = new double[elemIDs.size()];
				for (int j=0; j<elemIDs.size(); j++) {
					Integer elemID = elemIDs.get(j);
					Double particRate = indepParticRates.get(elemID);
					if (particRate != null)
						faultParticRates[j] = particRate;
				}
				Arrays.sort(faultParticRates);
				indepFaultParticRates.put(faultID, faultParticRates);
			}
		}
		
		IncrementalMagFreqDist eventMFD =
				MFDCalc.calcMFD(events, elementsInRegion, totalEventDuration*DAYS_PER_YEAR, mfdPlotMinMag, mfdPlotNumMag, delta);
		EvenlyDiscretizedFunc eventCmlMFD = eventMFD.getCumRateDist();
		
		ArrayList<EvenlyDiscretizedFunc> eventMFDs = Lists.newArrayList(eventMFD, eventCmlMFD);
		ArrayList<PlotCurveCharacterstics> eventMFDChars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE.darker()));
		
		ArrayList<RangePlotter> tasks = Lists.newArrayList();
		
		for (double[] magRange : magRanges) {
			double minMag = magRange[0];
			double maxMag = magRange[1];
			
			String magStr;
			String magFileStr;
			if (maxMag < 10) {
				magStr = (float)minMag+"=>"+(float)maxMag;
				magFileStr = (float)minMag+"_"+(float)maxMag;
			} else {
				magStr = (float)minMag+"+";
				magFileStr = magStr;
			}
			
			List<RuptureIdentifier> rupIdens = Lists.newArrayList();
			// just mag range, no specific fault section
			rupIdens.add(new MagRangeRuptureIdentifier(minMag, maxMag, elementsInRegion));
			// now each fault section
			for (Integer elementID : elementsOfInterest.keySet()) {
				rupIdens.add(new ElementMagRangeDescription("", elementID, minMag, maxMag));
			}
			
			for (RuptureIdentifier rupIden : rupIdens) {
				String plotDirName;
				String plotTitle;
				if (rupIden instanceof ElementMagRangeDescription) {
					int elementID = ((ElementMagRangeDescription)rupIden).getElementIDs().get(0);
					String name = elementsOfInterest.get(elementID);
					String nameFileSafe = name.replaceAll(" ", "_").toLowerCase();
					
					plotDirName = nameFileSafe+"_";
					plotTitle = name+" ";
				} else {
					plotDirName = "mag_only_";
					plotTitle = "";
				}
				plotDirName += magFileStr;
				plotTitle += "Mag "+magStr;
				
				File plotDir = new File(plotBaseDir, plotDirName);
				if (!plotDir.exists())
					plotDir.mkdir();
				
				List<? extends SimulatorEvent> matches = rupIden.getMatches(events);
				MinMaxAveTracker matchMagTrack = new MinMaxAveTracker();
				for (SimulatorEvent e : matches)
					matchMagTrack.addValue(e.getMagnitude());
				
				System.out.println("Num Matches: "+matches.size());
				
				IncrementalMagFreqDist matchMFD =
						MFDCalc.calcMFD(matches, null, -1, mfdPlotMinMag, mfdPlotNumMag, delta);
				EvenlyDiscretizedFunc matchCmlMFD = matchMFD.getCumRateDist();
				PlotSpec matchSpec = new PlotSpec(
						Lists.newArrayList(matchMFD, matchCmlMFD),
						Lists.newArrayList(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED),
								new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED.darker())),
								plotTitle+" Matched Event MFD", "Magnitude", "Count");
				
				writeMFDPlot(matchSpec, plotDir, plotDirName+"_match_mfd", mfdPlotMinMag, mfdPlotMaxMag, 1e1, 1e4);
				
				for (double windowDurationYears : durations) {
					RangePlotter rp = new RangePlotter(mfdPlotMinMag, delta, mfdPlotNumMag,
							events, rupIden, plotDirName, plotTitle, plotDir,
							minWindowDurationYears, randomizeEventTimes,
							matches, elementsInRegion, eventMFDs,
							eventMFDChars, windowDurationYears,
							particRanges, elemsByFaultMap, faultNameMap, indepFaultParticRatesList,
							mfdSpecialFaultSections, mfdSpecialFaultNames, mfdSpecialFaultColors);
					
					tasks.add(rp);
				}
			}
		}
		
		ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks);
		comp.computeThreaded();
	}
	
	private static String getParticRangeStr(double[] particRange) {
		if (particRange[1] > 9)
			return (float)particRange[0]+"+";
		else
			return (float)particRange[0]+"=>"+(float)particRange[1];
	}
	
	private static class RangePlotter implements Task {
		
		private double mfdPlotMinMag;
		private double delta;
		private int mfdPlotNumMag;
		private List<? extends SimulatorEvent> events;
		private RuptureIdentifier rupIden;
		private String plotDirName;
		private String plotTitle;
		private File plotDir;
		private double minWindowDurationYears;
		private boolean randomizeEventTimes;
		private List<? extends SimulatorEvent> matches;
		private HashSet<Integer> elementsInRegion;
		private ArrayList<EvenlyDiscretizedFunc> eventMFDs;
		private ArrayList<PlotCurveCharacterstics> eventMFDChars;
		private double windowDurationYears;
		private List<double[]> particRanges;
		private Map<Integer, List<Integer>> elemsByFaultMap;
		private Map<Integer, String> faultNameMap;
		private List<Map<Integer, double[]>> indepFaultParticRatesList;
		private List<HashSet<Integer>> mfdSpecialFaultSections;
		private List<String> mfdSpecialFaultNames;
		private List<Color> mfdSpecialFaultColors;
		
		public RangePlotter(double mfdPlotMinMag, double delta,
				int mfdPlotNumMag, List<? extends SimulatorEvent> events,
				RuptureIdentifier rupIden, String plotDirName, String plotTitle,
				File plotDir, double minWindowDurationYears,
				boolean randomizeEventTimes, List<? extends SimulatorEvent> matches,
				HashSet<Integer> elementsInRegion,
				ArrayList<EvenlyDiscretizedFunc> eventMFDs,
				ArrayList<PlotCurveCharacterstics> eventMFDChars,
				double windowDurationYears,
				List<double[]> particRanges,
				Map<Integer, List<Integer>> elemsByFaultMap,
				Map<Integer, String> faultNameMap,
				List<Map<Integer, double[]>> indepFaultParticRatesList,
				List<HashSet<Integer>> mfdSpecialFaultSections,
				List<String> mfdSpecialFaultNames,
				List<Color> mfdSpecialFaultColors) {
			this.mfdPlotMinMag = mfdPlotMinMag;
			this.delta = delta;
			this.mfdPlotNumMag = mfdPlotNumMag;
			this.events = events;
			this.rupIden = rupIden;
			this.plotDirName = plotDirName;
			this.plotTitle = plotTitle;
			this.plotDir = plotDir;
			this.minWindowDurationYears = minWindowDurationYears;
			this.randomizeEventTimes = randomizeEventTimes;
			this.matches = matches;
			this.elementsInRegion = elementsInRegion;
			this.eventMFDs = eventMFDs;
			this.eventMFDChars = eventMFDChars;
			this.windowDurationYears = windowDurationYears;
			this.particRanges = particRanges;
			this.elemsByFaultMap = elemsByFaultMap;
			this.faultNameMap = faultNameMap;
			this.indepFaultParticRatesList = indepFaultParticRatesList;
			this.mfdSpecialFaultSections = mfdSpecialFaultSections;
			this.mfdSpecialFaultNames = mfdSpecialFaultNames;
			this.mfdSpecialFaultColors = mfdSpecialFaultColors;
		}

		@Override
		public void compute() {
			EventsInWindowsMatcher eventsInWindows = new EventsInWindowsMatcher(
					events, rupIden, minWindowDurationYears, windowDurationYears,
					randomizeEventTimes);
			
			double timeDepDuration = eventsInWindows.getTotalWindowDurationYears();
			
			IncrementalMagFreqDist timeDepMFD =
					MFDCalc.calcMFD(eventsInWindows.getEventsInWindows(), elementsInRegion, timeDepDuration*DAYS_PER_YEAR, mfdPlotMinMag, mfdPlotNumMag, delta);
			EvenlyDiscretizedFunc timeDepCmlMFD = timeDepMFD.getCumRateDist();
			
			ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			funcs.addAll(eventMFDs);
			chars.addAll(eventMFDChars);
			
			funcs.add(timeDepMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			funcs.add(timeDepCmlMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED.darker()));
			
			// now special faults
			for (int i=0; i<mfdSpecialFaultSections.size(); i++) {
				HashSet<Integer> elements = mfdSpecialFaultSections.get(i);
				Color color = mfdSpecialFaultColors.get(i);
				String faultName = mfdSpecialFaultNames.get(i);
				
				IncrementalMagFreqDist specialMFD =
						MFDCalc.calcMFD(eventsInWindows.getEventsInWindows(), elements, timeDepDuration*DAYS_PER_YEAR, mfdPlotMinMag, mfdPlotNumMag, delta);
				specialMFD.setName("Time Dep MFD for: "+faultName);
				funcs.add(0, specialMFD);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, color));
			}
			
			String durationStr;
			if (windowDurationYears < 1d) {
				durationStr = (int)(windowDurationYears*DAYS_PER_YEAR+0.5)+" Day";
			} else {
				durationStr = (int)windowDurationYears+" Year";
				if (windowDurationYears < 10)
					durationStr = "0"+durationStr;
			}
			
			String title = durationStr+" "+plotTitle+" Post Event MFDs ("+matches.size()+" matches, "
					+df_hundreths.format(timeDepDuration)+" yr total window length)";
			String xAxisLabel = "Magnitude";
			String yAxisLabel = "Rate (1/day)";
			PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			
			String durationStrFName = durationStr.replaceAll(" ", "_").toLowerCase();
			
			double mfdPlotMaxMag = mfdPlotMinMag + delta * (mfdPlotNumMag-1);
			
			try {
				writeMFDPlot(spec, plotDir, plotDirName+"_mfds_"+durationStrFName, mfdPlotMinMag, mfdPlotMaxMag, 1e-7, 1e1);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}

			for (int i=0; i<particRanges.size(); i++) {
				double[] particRange = particRanges.get(i);
				CSVFile<String> csv = new CSVFile<String>(true);

				List<String> header = Lists.newArrayList("Fault ID", "Fault Name");

				String rangeStr = getParticRangeStr(particRange);

				header.add("Median Indep Expected ("+rangeStr+")");
				header.add("Median Time Dep Expected ("+rangeStr+")");
				header.add("Median Indep Rate ("+rangeStr+")");
				header.add("Median Time Dep Rate ("+rangeStr+")");
				header.add("Median Ratio ("+rangeStr+")");
				for (int j=2; j<7; j++)
					header.add(header.get(j).replace("Median", "Max"));

				csv.addLine(header);
				
				Map<Integer, Double> rates = calcParticRates(eventsInWindows.getEventsInWindows(),
						timeDepDuration, particRange[0], particRange[1]);

				for (Integer faultID : elemsByFaultMap.keySet()) {
					List<Integer> elements = elemsByFaultMap.get(faultID);
					double[] faultRates = new double[elements.size()];
					for (int j=0; j<elements.size(); j++) {
						int id = elements.get(j);
						Double rate = rates.get(id);
						if (rate == null)
							rate = 0d;
						faultRates[j] = rate;
					}
					// sort
					Arrays.sort(faultRates);
					
					double[] indepFaultRates = indepFaultParticRatesList.get(i).get(faultID);
					
					double medianIndep = median(indepFaultRates);
					double medianDep = median(faultRates);
					double maxIndep = indepFaultRates[indepFaultRates.length-1];
					double maxDep = faultRates[faultRates.length-1];
					
					List<String> line = Lists.newArrayList();
					line.add(faultID+"");
					line.add(faultNameMap.get(faultID));
					double durationDays = windowDurationYears * DAYS_PER_YEAR;
					line.addAll(getComparisonLines(medianDep, medianIndep, durationDays));
					line.addAll(getComparisonLines(maxDep, maxIndep, durationDays));
					
					csv.addLine(line);
				}
				
				// sort by num expected
				csv.sort(csv.getNumCols()-4, 1, new Comparator<String>() {
					
					@Override
					public int compare(String o1, String o2) {
						Double d1 = Double.parseDouble(o1);
						Double d2 = Double.parseDouble(o2);
						if (d1.isNaN())
							d1 = -1d;
						if (d2.isNaN())
							d2 = -1d;
						return -d1.compareTo(d2);
					}
				});
				String csvName = plotDirName+"_fault_partic_rates_"+durationStrFName+"_"+rangeStr.replaceAll("=>", "_")+".csv";
				try {
					csv.writeToFile(new File(plotDir, csvName));
				} catch (IOException e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
		}
		
	}
	
	private static List<String> getComparisonLines(double timeDep, double timeIndep, double duration) {
		List<String> lines = Lists.newArrayList();
		
		// indep expected
		lines.add(""+(timeIndep * duration));
		// dep expected
		lines.add(""+(timeDep * duration));
		// indep rate
		lines.add(""+timeIndep);
		// dep rate
		lines.add(""+timeDep);
		// ratio
		lines.add(""+(timeDep/timeIndep));
		
		return lines;
	}
	
	/**
	 * returns the median...array must be sorted
	 * 
	 * @param m
	 * @return
	 */
	static double median(double[] m) {
		int middle = (m.length)/2;  // subscript of middle element
		if (m.length%2 == 1) {
			// Odd number of elements -- return the middle one.
			return m[middle];
		} else {
			// Even number -- return average of middle two
			// Must cast the numbers to double before dividing.
			return (m[middle-1] + m[middle]) / 2.0;
		}
	}
	
	private static void writeMFDPlot(PlotSpec spec, File dir, String prefix, double minX, double maxX, double minY, double maxY) throws IOException {
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		gp.setYLog(true);
		
//		MinMaxAveTracker xTrack = new MinMaxAveTracker();
//		MinMaxAveTracker yTrack = new MinMaxAveTracker();
//		
//		for (DiscretizedFunc func : spec.getFuncs()) {
//			for (Point2D pt : func) {
//				xTrack.addValue(pt.getX());
//				if (pt.getY() > 0)
//					yTrack.addValue(pt.getY());
//			}
//		}
//		
//		gp.setUserBounds(xTrack.getMin(), xTrack.getMax(), yTrack.getMin()*0.9d, yTrack.getMax()*1.1d);
		
//		gp.setUserBounds(5.5d, 8d, 1e-7, 1e1);
		gp.setUserBounds(minX, maxX, minY, maxY);
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(18);
		
		gp.drawGraphPanel(spec);
		
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPDF(new File(dir, prefix+".pdf").getAbsolutePath());
		gp.saveAsPNG(new File(dir, prefix+".png").getAbsolutePath());
	}
	
	private static Map<Integer, Double> calcParticRates(
			List<? extends SimulatorEvent> events, double duration, double minMag, double maxMag) {
		Map<Integer, Double> rates = Maps.newHashMap();
		
		double eventRate = 1d/(duration*DAYS_PER_YEAR);
		
		for (SimulatorEvent e : events) {
			double mag = e.getMagnitude();
			if (mag < minMag || mag >= maxMag)
				continue;
			
			for (Integer elementID : e.getAllElementIDs()) {
				Double rate = rates.get(elementID);
				if (rate == null)
					rate = 0d;
				rate += eventRate;
				rates.put(elementID, rate);
			}
		}
		
		return rates;
	}
	
	private static double[] toArray(double... vals) {
		return vals;
	}

}
