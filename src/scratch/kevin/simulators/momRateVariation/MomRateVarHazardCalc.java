package scratch.kevin.simulators.momRateVariation;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.region.CaliforniaRegions.RELM_SOCAL;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.NGAWest_2014_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.erf.SimulatorFaultSystemSolution;
import scratch.kevin.simulators.erf.SubSectionBiulder;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class MomRateVarHazardCalc {
	
	private static enum CatalogTypes {
		RSQSIM,
		UCERF3_TD,
		UCERF3_ETAS
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/Simulators/mom_rate_hazard");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		// don't use idens in loading, but rather
		List<RuptureIdentifier> loadIdens = Lists.newArrayList();
		// only so cal ruptures
		Region region = new CaliforniaRegions.RELM_SOCAL();
//		Region region = new CaliforniaRegions.RELM_TESTING();
		loadIdens.add(new RegionIden(region));
		SimAnalysisCatLoader loader = new SimAnalysisCatLoader(true, loadIdens, false);
		List<? extends SimulatorEvent> events = loader.getEvents();
		List<SimulatorElement> elements = loader.getElements();
		
		// UCERF3
		File u3MainDir = new File("/home/kevin/Simulators/time_series/ucerf3_compare/2015_07_30-MID_VALUES");
		FaultSystemSolution u3Sol = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		File u3EtasCatalogs = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_04-spontaneous-10000yr-full_td-subSeisSupraNucl-gridSeisCorr/results_m4.bin");
		Map<Integer, SimulatorElement> u3Elems = UCERF3ComparisonAnalysis.loadElements(u3Sol.getRupSet());
		
		CatalogTypes[] types = CatalogTypes.values();
		boolean[] poissons = { false, true }; // must always be false first, will shuffle in place
		
		double[] durations = { 1d, 5d, 10d, 15d, 30d, 50d, 100d };
		
//		int windowLen = 75;
//		int windowLen = 25;
//		boolean before = false;
		
		int windowLen = 100;
		boolean before = true;
		
		outputDir = new File(outputDir, windowLen+"yr");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double[] taper;
		if (before)
			taper = SimulatorMomRateVarCalc.getOnlyBeforeWindowTaper(windowLen);
		else
			taper = SimulatorMomRateVarCalc.buildHanningTaper(windowLen);
		
		for (CatalogTypes type : types) {
			List<? extends SimulatorEvent> myEvents;
			List<List<SimulatorEvent>> eventLists;
			switch (type) {
			case RSQSIM:
				myEvents = events;
				break;
			case UCERF3_TD:
				eventLists = UCERF3ComparisonAnalysis.loadUCERF3Catalogs(
						u3MainDir, u3Sol, region, u3Elems, 0);
				myEvents = UCERF3ComparisonAnalysis.stitch(eventLists);
				double maxDelta = 0d;
				int maxDeltaIndex = -1;
				for (int i=1; i<myEvents.size(); i++) {
					SimulatorEvent e0 = myEvents.get(i-1);
					SimulatorEvent e1 = myEvents.get(i);
					double delta = e1.getTimeInYears() - e0.getTimeInYears();
					if (delta > maxDelta) {
						maxDelta = delta;
						maxDeltaIndex = i;
					}
					if (delta > 100) {
						System.out.println("100 yr delta found");
						System.out.println("\tDelta="+delta);
						System.out.println("\tindexBefore="+e0.getID());
						System.out.println("\tindexAfter="+e1.getID());
					}
				}
				
				System.out.println("Max U3 delta of "+maxDelta+" found at index="+maxDeltaIndex
						+", t="+myEvents.get(maxDeltaIndex).getTimeInYears());
				break;
			case UCERF3_ETAS:
				List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(u3EtasCatalogs, 5);
				eventLists = UCERF3_ETASComparisons.loadUCERF3EtasCatalogs(catalogs, u3Sol, region, u3Elems);
				myEvents = UCERF3ComparisonAnalysis.stitch(eventLists);
				break;

			default:
				throw new IllegalStateException();
			}
			for (boolean poisson : poissons) {
				double startYear = myEvents.get(0).getTimeInYears() + windowLen;
				double endYear = myEvents.get(myEvents.size()-1).getTimeInYears();
				
				String prefix = type.name();
				if (poisson)
					prefix += "_POISSON";
				
				System.out.println("Working on "+prefix);
				
				if (poisson) {
					List<SimulatorEvent> poissonEvents = Lists.newArrayList();
					double startSecs = myEvents.get(0).getTime();
					double endSecs = myEvents.get(myEvents.size()-1).getTime();
					double durationSecs = endSecs - startSecs;
					for (SimulatorEvent e : myEvents) {
						double timeSeconds = startSecs + Math.random()*(durationSecs);
						poissonEvents.add(e.cloneNewTime(timeSeconds, e.getID()));
					}
					Collections.sort(poissonEvents);
					myEvents = poissonEvents;
				}
				
				File myDir = new File(outputDir, prefix);
				Preconditions.checkState(myDir.exists() || myDir.mkdir());
				
				List<Double> yearsList = Lists.newArrayList();
				for (double year=startYear+windowLen; year<endYear-windowLen; year+=1d)
					yearsList.add(year);
				double[] years = Doubles.toArray(yearsList);
				double[] momRates;
				
//				String beforeStr = "";
//				if (before)
//					beforeStr = "_before";
//				File momRateFile = new File(myDir,
//						"mom_rates"+"_taper"+windowLen+"yr"+beforeStr+"_"+years.length+"yrs.bin");
//				if (momRateFile.exists()) {
//					momRates = MatrixIO.doubleArrayFromFile(momRateFile);
//					Preconditions.checkState(momRates.length == years.length);
//				} else {
					momRates = SimulatorMomRateVarCalc.calcTaperedMomRates(myEvents, years, taper);
//					MatrixIO.doubleArrayToFile(momRates, momRateFile);
//				}
				doMomRateCalc(myDir, events, years, momRates, durations, windowLen);
				doEventRateCalc(myDir, myEvents, years, momRates, 7d, durations, windowLen);
			}
		}
		
		
//		if (before) {
//			int startIndex = 0;
//			ArbitrarilyDiscretizedFunc origFunc = new ArbitrarilyDiscretizedFunc();
//			ArbitrarilyDiscretizedFunc newFunc = new ArbitrarilyDiscretizedFunc();
//			for (int i=0; i<years.length; i++) {
//				double endTime = years[i];
//				double startTime = endTime-windowLen;
//				
//				if (i < 1000)
//					origFunc.set(endTime, momRates[i]);
//				
//				momRates[i] = 0;
//				
//				for (int j=startIndex; j<events.size(); j++) {
//					EQSIM_Event e = events.get(j);
//					double t = e.getTimeInYears();
//					if (t < startTime) {
//						startIndex = j;
//						continue;
//					}
//					if (t > endTime)
//						break;
//					for (EventRecord rec : e)
//						momRates[i] += rec.getMoment();
//				}
//				
//				momRates[i] /= windowLen;
//				if (i < 1000)
//					newFunc.set(endTime, momRates[i]);
//			}
//			List<PlotElement> funcs = Lists.newArrayList();
//			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
//			funcs.add(origFunc);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//			funcs.add(newFunc);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//			DefaultXY_DataSet xy = new DefaultXY_DataSet();
//			for (EQSIM_Event e : events) {
//				if (e.getMagnitude() < 7d)
//					continue;
//				double t = e.getTimeInYears();
//				if (t < years[0])
//					continue;
//				if (t > years[1000])
//					break;
//				xy.set(t, (e.getMagnitude()-6)*1e19);
//			}
//			funcs.add(xy);
//			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.RED));
//			new GraphWindow(funcs, "Test", chars);
//		}
		
//		double[] hazard_durations = { 5d, 30d, 50d, 75d };
//		double hazardMinMag = 5.5d;
//		doHazardCalc(outputDir, events, elements, hazard_durations,
//				hazardMinMag, years, momRates);
	
//		doEventRateCalc(events, years, momRates, 7d, durations);
		
	}

	private static void doHazardCalc(File outputDir, List<SimulatorEvent> events,
			List<SimulatorElement> elements, double[] hazard_durations,
			double hazardMinMag, double[] years, double[] momRates)
			throws IOException {
		double upperMomRate = 1.5e19;
		double lowerMomRate = 9e18;
		int numAbove = 0;
		int numBelow = 0;
		int numBetween = 0;
		int numYears = 0;
		
		State[] states = new State[years.length];
		
		for (int i=0; i<years.length; i++) {
			if (momRates[i] > upperMomRate) {
				numAbove++;
				states[i] = State.ABOVE;
			} else if (momRates[i] < lowerMomRate) {
				numBelow++;
				states[i] = State.BELOW;
			} else {
				numBetween++;
				states[i] = State.BETWEEN;
			}
			numYears++;
		}
		System.out.println(numAbove+"/"+numYears+" ("+(float)(100d*numAbove/(double)numYears)+" %) above");
		System.out.println(numBelow+"/"+numYears+" ("+(float)(100d*numBelow/(double)numYears)+" %) below");
		System.out.println(numBetween+"/"+numYears+" ("+(float)(100d*numBetween/(double)numYears)+" %) between");
		
		// now build sub catalogs, use the midpoint in each span below/above/between
		int curStartIndex = 0;
		
		Map<State, List<Integer>> windowCenters = Maps.newHashMap();
		windowCenters.put(State.ABOVE, new ArrayList<Integer>());
		windowCenters.put(State.BELOW, new ArrayList<Integer>());
		windowCenters.put(State.BETWEEN, new ArrayList<Integer>());
		
		State curState = states[0];
		List<Double> curVals = Lists.newArrayList();
		
		for (int i=0; i<states.length; i++) {
			if (states[i] != curState) {
//				int endIndex = i-1;
//				int center = (endIndex+curStartIndex)/2;
				int center = curState.getCentralIndex(curVals)+curStartIndex;
				windowCenters.get(curState).add(center);
				
				curState = states[i];
				curStartIndex = i;
				curVals.clear();;
			}
			curVals.add(momRates[i]);
		}
		
		SubSectionBiulder subSectBuilder = new SubSectionBiulder(elements);
		
		Site site = new Site(new Location(34.055, -118.2467)); // LA Civic Center
		ScalarIMR imr = new NGAWest_2014_Averaged_AttenRel(null, false);
		imr.setParamDefaults();
		imr.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), 1d);
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
		SiteTranslator trans = new SiteTranslator();
		for (Parameter<?> param : imr.getSiteParams()) {
			trans.setParameterValue(param, datas);
			site.addParameter(param);
		}
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(imr.getIntensityMeasure());
		
		List<SimulatorEvent> hazardEvents = Lists.newArrayList();
		for (SimulatorEvent event : events)
			if (event.getMagnitude() >= hazardMinMag)
				hazardEvents.add(event);
		
		for (double duration : hazard_durations) {
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			for (State state : State.values()) {
				System.out.println("Calculating "+(int)duration+"yr, "+state);
				SimulatorFaultSystemSolution fss = buildFSS(
						events, subSectBuilder, hazardMinMag, windowCenters.get(state), years, duration);
				DiscretizedFunc func = calcHazardCurve(fss, duration, site, imr, xVals);
				func.setName(state.name());
				funcs.add(func);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, state.c));
			}
//			// now full catalog
//			SimulatorFaultSystemSolution fss = SimulatorFaultSystemSolution.build(
//					subSectBuilder, hazardEvents, totalDurationYears);
//			DiscretizedFunc func = calcHazardCurve(fss, duration, site, imr, xVals);
//			func.setName("Full Catalog");
//			funcs.add(0, func);
//			chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, (int)duration+"yr Hazard Curves",
					"1s SA", "Exceed Prob");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setAxisRange(new Range(1e-2, 3), new Range(1e-4, 1));
			gw.saveAsPNG(new File(outputDir, "curves_"+(int)duration+"yr.png").getAbsolutePath());
		}
		
//		FaultSystemRupSet rupSet = buildRupSet(elements, events, totalDurationYears, subSectBuilder);
	}
	
	private static enum State {
		ABOVE(Color.RED),
		BELOW(Color.BLUE),
		BETWEEN(Color.GREEN);
		
		private Color c;
		
		private State(Color c) {
			this.c = c;
		}
		
		private int getCentralIndex(List<Double> values) {
			Preconditions.checkArgument(!values.isEmpty());
			
			switch (this) {
			case ABOVE:
				// find max value
				int maxIndex = -1;
				double maxVal = 0;
				for (int i=0; i<values.size(); i++) {
					double v = values.get(i);
					if (v > maxVal) {
						maxVal = v;
						maxIndex = i;
					}
				}
				return maxIndex;
			case BELOW:
				// find min value
				int minIndex = -1;
				double minVal = Double.POSITIVE_INFINITY;
				for (int i=0; i<values.size(); i++) {
					double v = values.get(i);
					if (v < minVal) {
						minVal = v;
						minIndex = i;
					}
				}
				return minIndex;
			case BETWEEN:
				// find median
				double median = DataUtils.median(Doubles.toArray(values));
				int closestIndex = 0-1;
				double closestDelta = Double.POSITIVE_INFINITY;
				for (int i=0; i<values.size(); i++) {
					double delta = Math.abs(values.get(i)-median);
					if (delta < closestDelta) {
						closestDelta = delta;
						closestIndex = i;
					}
				}
				return closestIndex;

			default:
				throw new IllegalStateException("Unknown Satate");
			}
		}
	}
	
	private static SimulatorFaultSystemSolution buildFSS(
			List<SimulatorEvent> allEvents, SubSectionBiulder subSectBuilder, double minMag,
			List<Integer> windowCenters, double[] years, double forecastDuration) {
		
		List<SimulatorEvent> includedEvents = Lists.newArrayList();
		double includedDuration = 0;
		
		Preconditions.checkArgument(!windowCenters.isEmpty());
		
		for (int center : windowCenters) {
			double year = years[center];
			double yearSecs = year*General_EQSIM_Tools.SECONDS_PER_YEAR;
			double endYear = year + forecastDuration;
			
			int startIndex = SimulatorMomRateVarCalc.findFirstEventIndexAfter(yearSecs, allEvents);
			
			for (int i=startIndex; i<allEvents.size(); i++) {
				SimulatorEvent e = allEvents.get(i);
				double t = e.getTimeInYears();
				Preconditions.checkState(t >= year);
				if (t > endYear)
					break;
				if (e.getMagnitude() >= minMag)
					includedEvents.add(e);
			}
			includedDuration += forecastDuration;
		}
		
		System.out.println("Built SimFSS from "+includedEvents.size()+" events, dur="+includedDuration);
		
		return SimulatorFaultSystemSolution.build(subSectBuilder, includedEvents, includedDuration);
	}
	
	private static DiscretizedFunc calcHazardCurve(FaultSystemSolution fss, double duration,
			Site site, ScalarIMR imr, DiscretizedFunc xVals) {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(fss);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(duration);
		erf.updateForecast();
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		DiscretizedFunc func = xVals.deepClone();
		func = HazardCurveSetCalculator.getLogFunction(func);
		
		calc.getHazardCurve(func, site, imr, erf);
		
		func = HazardCurveSetCalculator.unLogFunction(xVals, func);
		
		return func;
	}
	
	private static void doEventRateCalc(File outputDir, List<? extends SimulatorEvent> events, double[] years,
			double[] moRates, double minMag, double[] durations, int durationBefore) throws IOException {
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		int colorIndex = 0;
		
		for (double duration : durations) {
			DiscretizedFunc hist = calcEventRatesForMomRates(events, years, moRates, minMag, duration);
			hist.setName((int)duration+"yr");
			funcs.add(hist);
			if (colorIndex == colors.size())
				colorIndex = 0;
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(colorIndex++)));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Event Rates vs Moment Rates",
				"Moment Rate for "+durationBefore+"yrs Before", "Rate M>="+(double)minMag+" Following");
		spec.setLegendVisible(true);
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(22);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, true, false);
		gp.getChartPanel().setSize(1000, 800);
		File outputFile = new File(outputDir, "event_rate_following");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		
		// now gain func
		int countAbove = 0;
		for (SimulatorEvent e : events)
			if (e.getMagnitude() >= minMag)
				countAbove++;
		double simDuration = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
		double rateAbove = (double)countAbove/simDuration;
		
		Range xRange = new Range(1e18, 1e20);
		
		List<DiscretizedFunc> gainFuncs = Lists.newArrayList();
		for (int i=0; i<durations.length; i++) {
			double duration = durations[i];
			DiscretizedFunc countFunc = funcs.get(i);
			DiscretizedFunc gainFunc = new ArbitrarilyDiscretizedFunc();
			gainFunc.setName(countFunc.getName());
			
			for (Point2D pt : countFunc) {
				double count = pt.getY();
				double myRate = count/duration;
				double gain = myRate/rateAbove;
				gainFunc.set(pt.getX(), gain);
			}
			while (gainFunc.getMinX() < xRange.getLowerBound())
				xRange = new Range(xRange.getLowerBound()/10, xRange.getUpperBound());
			while (gainFunc.getMaxX() > xRange.getUpperBound())
				xRange = new Range(xRange.getLowerBound()/10, xRange.getUpperBound());
			
			gainFuncs.add(gainFunc);
		}
		
		ArbitrarilyDiscretizedFunc flatGainFunc = new ArbitrarilyDiscretizedFunc();
		flatGainFunc.set(xRange.getLowerBound(), 1d);
		flatGainFunc.set(xRange.getUpperBound(), 1d);
		gainFuncs.add(0, flatGainFunc);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		spec = new PlotSpec(gainFuncs, chars, "Event Rate Gain vs Moment Rates",
				"Moment Rate for "+durationBefore+"yrs Before", "Rate Gain M>="+(double)minMag+" Following");
		spec.setLegendVisible(true);
		gp.drawGraphPanel(spec, true, false, xRange, new Range(0d, 2.5d));
		gp.getChartPanel().setSize(1000, 800);
		outputFile = new File(outputDir, "event_rate_gain");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
	}
	
	private static ArbitrarilyDiscretizedFunc calcEventRatesForMomRates(
			List<? extends SimulatorEvent> events, double[] years, double[] moRates, double minMag, double durationYears) {
		Preconditions.checkArgument(moRates.length == years.length);
		
		double maxMoRate = StatUtils.max(moRates);
		double minMoRate = StatUtils.min(moRates);
		Preconditions.checkState(minMoRate > 0d, "Cannot have mo rate of zero");
		
		// bin in Log10 space
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(
				Math.log10(minMoRate), Math.log10(maxMoRate), 0.1);
		int[] binCounts = new int[hist.size()];
		
		int startEventIndex = 0;
		
		for (int i=0; i<years.length; i++) {
			double moRate = moRates[i];
			int xIndex = hist.getClosestXIndex(Math.log10(moRate));
			binCounts[xIndex]++;
			
			double startYear = years[i];
			double endYear = startYear + durationYears;
			
			for (int j=startEventIndex; j<events.size(); j++) {
				SimulatorEvent e = events.get(j);
				double t = e.getTimeInYears();
				if (t < startYear) {
					startEventIndex = j;
					continue;
				}
				if (t > endYear)
					break;
				if (e.getMagnitude() >= minMag) {
					hist.add(xIndex, 1d);
				}
			}
		}
		
		for (int i=0; i<hist.size(); i++)
			if (binCounts[i] > 0)
				hist.set(i, hist.getY(i)/(double)binCounts[i]);
		
		// now go back to linear space
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<hist.size(); i++)
			if (binCounts[i] >= 10)
				ret.set(Math.pow(10d, hist.getX(i)), hist.getY(i));
		
		return ret;
	}
	
	private static void doMomRateCalc(File outputDir, List<? extends SimulatorEvent> events, double[] years,
			double[] moRates, double[] durations, int durationBefore) throws IOException {
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		int colorIndex = 0;
		
		for (double duration : durations) {
			DiscretizedFunc hist = calcMomRatesForMomRates(events, years, moRates, duration);
			hist.setName((int)duration+"yr");
			funcs.add(hist);
			if (colorIndex == colors.size())
				colorIndex = 0;
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(colorIndex++)));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 2f, colors.get(colorIndex++)));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Moment Rate Before/After",
				"Moment Rate for "+durationBefore+"yrs Before", "Moment Rate Following");
		spec.setLegendVisible(true);
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(22);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, true, false);
		gp.getChartPanel().setSize(1000, 800);
		File outputFile = new File(outputDir, "mom_rate_following");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		
		double meanMoRate = StatUtils.mean(moRates);
		
		Range xRange = new Range(1e18, 1e20);
		List<DiscretizedFunc> gainFuncs = Lists.newArrayList();
		for (int i=0; i<durations.length; i++) {
			DiscretizedFunc countFunc = funcs.get(i);
			DiscretizedFunc gainFunc = new ArbitrarilyDiscretizedFunc();
			gainFunc.setName(countFunc.getName());
			
			for (Point2D pt : countFunc) {
//				double gain = pt.getY()/pt.getX();
				double gain = pt.getY()/meanMoRate;
				gainFunc.set(pt.getX(), gain);
			}
			while (gainFunc.getMinX() < xRange.getLowerBound())
				xRange = new Range(xRange.getLowerBound()/10, xRange.getUpperBound());
			while (gainFunc.getMaxX() > xRange.getUpperBound())
				xRange = new Range(xRange.getLowerBound()/10, xRange.getUpperBound());
			
			gainFuncs.add(gainFunc);
		}
		
		ArbitrarilyDiscretizedFunc flatGainFunc = new ArbitrarilyDiscretizedFunc();
		flatGainFunc.set(xRange.getLowerBound(), 1d);
		flatGainFunc.set(xRange.getUpperBound(), 1d);
		gainFuncs.add(0, flatGainFunc);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		spec = new PlotSpec(gainFuncs, chars, "Moment Rate Gain vs Moment Rates",
				"Moment Rate for "+durationBefore+"yrs Before", "Moment Rate Gain Following");
		spec.setLegendVisible(true);
		gp.drawGraphPanel(spec, true, false, xRange, new Range(0d, 2.5d));
		gp.getChartPanel().setSize(1000, 800);
		outputFile = new File(outputDir, "mom_rate_gain");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		
		// now histogram
		double minLogMoRate = Math.log10(StatUtils.min(moRates));
		double maxLogMoRate = Math.log10(StatUtils.max(moRates));
		HistogramFunction momRateHist = HistogramFunction.getEncompassingHistogram(minLogMoRate, maxLogMoRate, 0.1);
		HistogramFunction momRateAfterLowHist =
				new HistogramFunction(momRateHist.getMinX(), momRateHist.getMaxX(), momRateHist.size());
		double lowerLimit = StatUtils.percentile(moRates, 25);
		
		double overallMean = 0d;
		double afterLowMean = 0d;
		for (int i=0; i<moRates.length; i++) {
			momRateHist.add(Math.log10(moRates[i]), 1d);
			overallMean += moRates[i];
			if (moRates[i] < lowerLimit && (i+100)<moRates.length) {
				momRateAfterLowHist.add(Math.log10(moRates[i+100]), 1d);
				afterLowMean += moRates[i+100];
			}
		}
		overallMean /= momRateHist.calcSumOfY_Vals();
		afterLowMean /= momRateAfterLowHist.calcSumOfY_Vals();
		double fract = momRateAfterLowHist.calcSumOfY_Vals()/momRateHist.calcSumOfY_Vals();
		momRateAfterLowHist.normalizeBySumOfY_Vals();
		System.out.println("Fraction below "+lowerLimit+": "+fract);
		momRateAfterLowHist.scale(fract);
		momRateHist.normalizeBySumOfY_Vals();
		DecimalFormat meanDF = new DecimalFormat("0.##E0");
		momRateHist.setName("All Years (mean="+meanDF.format(overallMean)+")");
		double gain = afterLowMean / overallMean;
		momRateAfterLowHist.setName("100yrs After Bottom 25% (gain="+new DecimalFormat("0.00").format(gain)+")");
		
		funcs = Lists.newArrayList();
		chars = Lists.newArrayList();
		
		funcs.add(momRateHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		spec = new PlotSpec(funcs, chars, "Moment Rate Histogram",
				"Log10(Moment Rate)", "Density");
		spec.setLegendVisible(false);
		gp.drawGraphPanel(spec, false, false);
		gp.getChartPanel().setSize(1000, 800);
		outputFile = new File(outputDir, "mom_rate_hist");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		
		funcs.add(momRateAfterLowHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.RED));
		
		spec = new PlotSpec(funcs, chars, "Moment Rate Histogram",
				"Log10(Moment Rate)", "Density");
		spec.setLegendVisible(true);
		gp.drawGraphPanel(spec, false, false);
		gp.getChartPanel().setSize(1000, 800);
		outputFile = new File(outputDir, "mom_rate_hist_after_low");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
	}
	
	private static UncertainArbDiscDataset calcMomRatesForMomRates(
			List<? extends SimulatorEvent> events, double[] years, double[] moRates, double durationYears) {
		Preconditions.checkArgument(moRates.length == years.length);
		
		double maxMoRate = StatUtils.max(moRates);
		double minMoRate = StatUtils.min(moRates);
//		for (int i=0; i<years.length; i++) {
//			if (moRates[i] == minMoRate) {
//				System.out.println("Minimum of "+minMoRate+" at bin "+i+"/"+years.length+", t="+years[i]);
//				break;
//			}
//		}
		Preconditions.checkState(minMoRate > 0d, "Cannot have mo rate of zero");
		
		// bin in Log10 space
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(
				Math.log10(minMoRate), Math.log10(maxMoRate), 0.1);
		List<List<Double>> momRateAfters = Lists.newArrayList();
		for (int i=0; i<hist.size(); i++)
			momRateAfters.add(new ArrayList<Double>());
		
		int startEventIndex = 0;
		
		for (int i=0; i<years.length; i++) {
			double moRate = moRates[i];
			int xIndex = hist.getClosestXIndex(Math.log10(moRate));
			
			double startYear = years[i];
			double endYear = startYear + durationYears;
			
			double momentAfter = 0d;
			
			for (int j=startEventIndex; j<events.size(); j++) {
				SimulatorEvent e = events.get(j);
				double t = e.getTimeInYears();
				if (t < startYear) {
					startEventIndex = j;
					continue;
				}
				if (t > endYear)
					break;
				for (EventRecord rec : e)
					momentAfter += rec.getMoment();
			}
//			hist.add(xIndex, momentAfter/durationYears);
			momRateAfters.get(xIndex).add(momentAfter/durationYears);
		}
//		
//		for (int i=0; i<hist.size(); i++)
//			if (binCounts[i] > 0)
//				hist.set(i, hist.getY(i)/(double)binCounts[i]);
//		
//		// now go back to linear space
//		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
//		for (int i=0; i<hist.size(); i++)
//			if (binCounts[i] >= 10)
//				ret.set(Math.pow(10d, hist.getX(i)), hist.getY(i));
//		
//		return ret;
		ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<hist.size(); i++) {
			double x = Math.pow(10d, hist.getX(i));
			List<Double> vals = momRateAfters.get(i);
			if (vals.size() < 100)
				continue;
			
			double[] valsArray = Doubles.toArray(vals);
			
			double mean = StatUtils.mean(valsArray);
			double stdDev = Math.sqrt(StatUtils.variance(valsArray));
			
			meanFunc.set(x, mean);
			lowerFunc.set(x, mean - stdDev);
			upperFunc.set(x, mean + stdDev);
		}
		
		return new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
	}

}
