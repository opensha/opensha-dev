package scratch.kevin.simulators.erf;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.BPT_AperiodicityParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.EventsInWindowsMatcher;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SupraSeisRupIden;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;

public class TimeDepFSS_ERF_Simulator_Test {

	public static void main(String[] args) throws IOException {
		boolean interactive;
		String outputPrefix;
		File outputDir;
		File dataDir;
		int numTrials;
		final MagDependentAperiodicityOptions default_cov = MagDependentAperiodicityOptions.LOW_VALUES;
		MagDependentAperiodicityOptions cov;
		boolean test_time_indep = false;
		double defaultForecastDuration = 50d;
		double forecastDuration;
		if (args.length > 0) {
			interactive = false;
			outputDir = new File(args[0]);
			dataDir = outputDir.getParentFile();
			outputPrefix = args[1];
			numTrials = Integer.parseInt(args[2]);
			if (args.length > 3)
				cov = MagDependentAperiodicityOptions.valueOf(args[3]);
			else
				cov = default_cov;
			if (args.length > 4)
				forecastDuration = Double.parseDouble(args[4]);
			else
				forecastDuration = defaultForecastDuration;
			
		} else {
			interactive = true;
			dataDir = new File("/home/kevin/Simulators");
			outputDir = dataDir;
			numTrials = 100000;
			outputPrefix = "erf_audit_"+numTrials;
			cov = default_cov;
			forecastDuration = defaultForecastDuration;
		}
		if (test_time_indep)
			outputPrefix += "_INDEP";
		File geomFile = new File(dataDir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry from "+geomFile.getAbsolutePath());
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dataDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		boolean supraSeismogenic = true;
		double minMag = 6.5;
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		if (supraSeismogenic)
			rupIdens.add(new SupraSeisRupIden(tools));
		if (minMag > 0)
			rupIdens.add(new MagRangeRuptureIdentifier(minMag, 10d));
		LogicalAndRupIden andIden = new LogicalAndRupIden(rupIdens);
		rupIdens = Lists.newArrayList();
		rupIdens.add(andIden);
		
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList(), rupIdens);
		
//		Region region = null;
		Region region = new CaliforniaRegions.RELM_SOCAL();
		
		if (region != null) {
			Map<Integer, Boolean> elementsInRegionsCache = Maps.newHashMap();
			
			// just uese centers since they're small enough elements
			for (SimulatorElement elem : tools.getElementsList())
				elementsInRegionsCache.put(elem.getID(), region.contains(elem.getCenterLocation()));
			
			List<SimulatorEvent> eventsInRegion = Lists.newArrayList();
			for (SimulatorEvent e : events) {
				for (int elemID : e.getAllElementIDs()) {
					if (elementsInRegionsCache.get(elemID)) {
						eventsInRegion.add(e);
						break;
					}
				}
			}
			
			System.out.println(eventsInRegion.size()+"/"+events.size()+" events in region: "+region.getName());
			events = eventsInRegion;
		}
		
//		if (supraSeismogenic) {
//			List<EQSIM_Event> supraSeisEvents = Lists.newArrayList();
//			for (EQSIM_Event e : events)
//				if (tools.isEventSupraSeismogenic(e, Double.NaN))
//					supraSeisEvents.add(e);
//			events = supraSeisEvents;
//		}
		
		double durationYears = General_EQSIM_Tools.getSimulationDurationYears(events);
		
		List<SimulatorElement> elements = tools.getElementsList();
		SubSectionBiulder subSectBuilder = new SubSectionBiulder(elements);
		FaultSystemRupSet rupSet = SimulatorFaultSystemSolution.buildRupSet(elements, events, durationYears, subSectBuilder);
		FaultSystemSolution sol = new SimulatorFaultSystemSolution(rupSet, subSectBuilder, events, durationYears);
		System.out.println("Precombine sol has "+rupSet.getNumRuptures()+" rups");
		sol = combineIdenticalRups(sol);
		rupSet = sol.getRupSet();
		System.out.println("Combined sol has "+rupSet.getNumRuptures()+" rups");
		
		// now map each EQSIM event to a rup index
		Map<Integer, Integer> simToRupIndex = Maps.newHashMap();
		Map<UniqueIndexSet, Integer> subSectsToRupMap = Maps.newHashMap();
		for (int r=0; r<rupSet.getNumRuptures(); r++)
			subSectsToRupMap.put(new UniqueIndexSet(new HashSet<Integer>(rupSet.getSectionsIndicesForRup(r))), r);
		Map<Integer, Integer> elemsIDsToSubSects = subSectBuilder.getElemIDToSubSectsMap();
		for (SimulatorEvent e: events) {
			HashSet<Integer> subSects = new HashSet<Integer>();
			for (int elemID : e.getAllElementIDs())
				subSects.add(elemsIDsToSubSects.get(elemID));
			if (!subSectsToRupMap.containsKey(new UniqueIndexSet(subSects))) {
				System.out.println("Hmm, couldn't find rup: "+Joiner.on(",").join(subSects));
				System.out.println("Event mag: "+e.getMagnitude());
				for (int r=0; r<rupSet.getNumRuptures(); r++) {
					List<Integer> subSectsForRup = rupSet.getSectionsIndicesForRup(r);
					if (subSectsForRup.size() == subSects.size()) {
						boolean match = true;
						for (int s : subSectsForRup) {
							if (!subSects.contains(s)) {
								match = false;
								break;
							}
						}
						if (match) {
							System.out.println("But there is a match! rup "+r+": "+Joiner.on(",").join(subSectsForRup));
							break;
						}
					}
				}
			}
			int rupIndex = subSectsToRupMap.get(new UniqueIndexSet(subSects));
			simToRupIndex.put(e.getID(), rupIndex);
		}
//		System.out.println("First 10 mappings:");
//		for (int i=0; i<10; i++) {
//			int simID = events.get(i).getID();
//			System.out.println("\t"+simID+": "+simToRupIndex.get(simID));
//		}
//		System.exit(0);
		
		double minDuration = 1000;
		double maxDuration = durationYears-forecastDuration;
		double durationDelta = maxDuration - minDuration;
		double firstEventTimeYears = events.get(0).getTimeInYears();
		
		int timeSpanStartYear = 2013;
		
		double[] forecastOccurences = new double[rupSet.getNumRuptures()];
		double[] simulatorOccurances = new double[rupSet.getNumRuptures()];
		
		FaultSystemSolutionERF indepERF = null;
		if (test_time_indep) {
			indepERF = new FaultSystemSolutionERF(sol);
			indepERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			indepERF.getTimeSpan().setDuration(forecastDuration);
			// this will make it use the simulator dates of last events
			indepERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			indepERF.updateForecast();
		}
		
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<numTrials; i++) {
			System.out.println("Performing trial "+i);
			double randDuration = durationDelta*Math.random();
			
			double randTrainEndTime = firstEventTimeYears+minDuration+randDuration;
			double forecastEndTime = randTrainEndTime + forecastDuration;
			List<SimulatorEvent> eventsBefore = Lists.newArrayList();
			List<SimulatorEvent> eventsDuring = Lists.newArrayList();
			for (SimulatorEvent e : events) {
				if (e.getTimeInYears() <= randTrainEndTime)
					eventsBefore.add(e);
				else if (e.getTimeInYears() <= forecastEndTime)
					eventsDuring.add(e);
				else
					break;
			}
			System.out.println("Events during: "+eventsDuring.size());
//			try {
//				Thread.sleep(1000);
//			} catch (InterruptedException e1) {}
			
			FaultSystemSolutionERF erf;
			if (test_time_indep) {
				erf = indepERF;
			} else {
				populateOpenIntervals(subSectBuilder, eventsBefore, randTrainEndTime, timeSpanStartYear);
				
				erf = new FaultSystemSolutionERF(sol);
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
				erf.setParameter(MagDependentAperiodicityParam.NAME, cov);
				erf.getTimeSpan().setStartTime(timeSpanStartYear);
				erf.getTimeSpan().setDuration(forecastDuration);
				// make sure this didn't clear the start year
				Preconditions.checkState(erf.getTimeSpan().getStartTimeYear() == timeSpanStartYear);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
				erf.updateForecast();
			}
			
			if (!test_time_indep || i == 0) {
				// don't bother if this is a time independent erf test
				for (ProbEqkSource source : erf) {
					String srcName = source.getName();
					int rupIndex = Integer.parseInt(srcName.substring(srcName.indexOf("#")+1, srcName.indexOf(";")));
					double prob = source.computeTotalProb();
					forecastOccurences[rupIndex] += prob;
				}
			}
			// now log what actually happened
			// do it this way to only count rups once even if they happen multiple times
			HashSet<Integer> rupsOccurredThisTime = new HashSet<Integer>();
			for (SimulatorEvent e : eventsDuring) {
				int rupIndex = simToRupIndex.get(e.getID());
				rupsOccurredThisTime.add(rupIndex);
			}
			for (int rupIndex : rupsOccurredThisTime) {
				simulatorOccurances[rupIndex] += 1d;
			}
		}
		if (test_time_indep) {
			for (int i=0; i<forecastOccurences.length; i++)
				forecastOccurences[i] *= (double)numTrials;
		}
		watch.stop();
		long secs = watch.elapsed(TimeUnit.SECONDS);
		System.out.println("Time per: "+(float)secs/(float)numTrials+" secs");
		
		DefaultXY_DataSet ratioData = new DefaultXY_DataSet(simulatorOccurances, forecastOccurences);
		double max = StatUtils.max(simulatorOccurances);
		max = Math.max(max, StatUtils.max(forecastOccurences));
		ArbitrarilyDiscretizedFunc eventRatio = new ArbitrarilyDiscretizedFunc();
		eventRatio.set(0d, 0d);
		eventRatio.set(1e-1, 1e-1);
		eventRatio.set(max, max);
		
		List<XY_DataSet> elems = Lists.newArrayList();
		elems.add(ratioData);
		elems.add(eventRatio);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		double numOccur = StatUtils.sum(simulatorOccurances);
		double numForecast = StatUtils.sum(forecastOccurences);
		System.out.println("Tot sim occurences: "+numOccur);
		System.out.println("Tot predicted: "+numForecast);
		
		if (interactive) {
			GraphWindow gw = new GraphWindow(elems, "Time Dependent Simulator Audit ("+numTrials+" sims)", chars);
			gw.setX_AxisLabel("Actual RSQSim Occurances (sum="+(float)numOccur+")");
			gw.setY_AxisLabel("FSS ERF Predicted Occurances (sum="+(float)numForecast+")");
			Range range = new Range(0, max);
			gw.setAxisRange(range, range);
			gw.setSize(1000, 800);
			gw.saveAsPNG("/tmp/fss_erf_audit.png");
			gw.setXLog(true);
			gw.setYLog(true);
			range = new Range(1d, max);
			gw.setAxisRange(range, range);
			gw.saveAsPNG("/tmp/fss_erf_audit_log.png");
		}
		if (outputPrefix != null) {
			File simOccurFile = new File(outputDir, outputPrefix+"_sim_occurances.bin");
			File erfOccurFile = new File(outputDir, outputPrefix+"_erf_occurances.bin");
			MatrixIO.doubleArrayToFile(simulatorOccurances, simOccurFile);
			MatrixIO.doubleArrayToFile(forecastOccurences, erfOccurFile);
		}
	}
	
	static class UniqueIndexSet {
		private int[] sortedIndexes;
		
		public UniqueIndexSet(HashSet<Integer> indexes) {
			List<Integer> list = Lists.newArrayList(indexes);
			Collections.sort(list);
			sortedIndexes = Ints.toArray(list);
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(sortedIndexes);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			UniqueIndexSet other = (UniqueIndexSet) obj;
			if (!Arrays.equals(sortedIndexes, other.sortedIndexes))
				return false;
			return true;
		}
	}
	
	static FaultSystemSolution combineIdenticalRups(FaultSystemSolution orig) {
		FaultSystemRupSet origRupSet = orig.getRupSet();
		
		List<FaultSectionPrefData> faultSectionData = origRupSet.getFaultSectionDataList();
		double[] sectSlipRates = origRupSet.getSlipRateForAllSections();
		double[] sectSlipRateStdDevs = origRupSet.getSlipRateStdDevForAllSections();
		double[] sectAreas = origRupSet.getAreaForAllSections();
		String info = origRupSet.getInfoString();
		
		int origNumRups = origRupSet.getNumRuptures();
		int curRupCount = 0;
		Map<UniqueIndexSet, Integer> rupIndexMap = Maps.newHashMap();
		
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		double[] mags = new double[origNumRups];
		double[] rakes = new double[origNumRups];
		double[] rupAreas = new double[origNumRups];
		double[] rupLengths = new double[origNumRups];
		
		int[] duplicateCounts = new int[origNumRups];
		
		double[] rates = new double[origNumRups];
		
		for (int r=0; r<origRupSet.getNumRuptures(); r++) {
			UniqueIndexSet rupIndexes = new UniqueIndexSet(new HashSet<Integer>(origRupSet.getSectionsIndicesForRup(r)));
			Integer newIndex = rupIndexMap.get(rupIndexes);
			if (newIndex == null) {
				newIndex = curRupCount++;
				rupIndexMap.put(rupIndexes, newIndex);
				sectionForRups.add(origRupSet.getSectionsIndicesForRup(r));
			}
			
			mags[newIndex] += origRupSet.getMagForRup(r);
			rakes[newIndex] = origRupSet.getAveRakeForRup(r);
			rupAreas[newIndex] += origRupSet.getAreaForRup(r);
			rupLengths[newIndex] += origRupSet.getLengthForRup(r);
			duplicateCounts[newIndex]++;
			rates[newIndex] += orig.getRateForRup(r);
		}
		
		mags = Arrays.copyOf(mags, curRupCount);
		rakes = Arrays.copyOf(rakes, curRupCount);
		rupAreas = Arrays.copyOf(rupAreas, curRupCount);
		rupLengths = Arrays.copyOf(rupLengths, curRupCount);
		rates = Arrays.copyOf(rates, curRupCount);
		
		for (int i=0; i<curRupCount; i++) {
			double cnt = (double)duplicateCounts[i];
			mags[i] /= cnt;
			rupAreas[i] /= cnt;
			rupLengths[i] /= cnt;
		}
		
		FaultSystemRupSet rupSet = new FaultSystemRupSet(faultSectionData, sectSlipRates,
				sectSlipRateStdDevs, sectAreas, sectionForRups, mags, rakes, rupAreas, rupLengths, info);
		
		return new FaultSystemSolution(rupSet, rates);
	}
	
	private static void populateOpenIntervals(SubSectionBiulder builder, List<SimulatorEvent> events, double catalogTime,
			int timeSpanStartYear) {
		List<FaultSectionPrefData> subSects = builder.getSubSectsList();
		// clear any open intervals
		for (FaultSectionPrefData subSect : subSects)
			subSect.setDateOfLastEvent(Long.MIN_VALUE);
		
		GregorianCalendar intervalBasis = new GregorianCalendar(timeSpanStartYear, 0, 0);
		
		Map<Integer, Integer> elemToSubSectMap = builder.getElemIDToSubSectsMap();
		
		for (int i=events.size(); --i>=0;) {
			SimulatorEvent e = events.get(i);
			double openInterval = catalogTime - e.getTimeInYears();
			Preconditions.checkState(openInterval >= 0);
			
			long relativeEventTime = LastEventData.calcDate(intervalBasis, openInterval).getTimeInMillis();
			
			for (int elemID : e.getAllElementIDs()) {
				FaultSectionPrefData subSect = subSects.get(elemToSubSectMap.get(elemID));
				if (subSect.getDateOfLastEvent() == Long.MIN_VALUE)
					subSect.setDateOfLastEvent(relativeEventTime);
			}
		}
	}

}
