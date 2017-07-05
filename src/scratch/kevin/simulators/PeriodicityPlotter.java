package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.pdfbox.exceptions.COSVisitorException;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.QuietPeriodIdenMatcher;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.RandomDistType;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class PeriodicityPlotter {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		boolean doRandom = true;
		boolean display = false;
		boolean displayEventTimes = false;
		RandomDistType randDistType = RandomDistType.ACTUAL;
		boolean randSplitMults = true;
		
		File writeDir = new File(dir, "period_plots");
		if (!writeDir.exists())
			writeDir.mkdir();
		
		double[] interevent_mags = { 5d, 5.5d, 6d, 6.5d, 7d, 7.5d };
		for (double mag : interevent_mags) {
			List<RuptureIdentifier> rupIdens = Lists.newArrayList();
			
			rupIdens.add(new ElementMagRangeDescription("Cholame",
					ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, mag, 10d));
			
			rupIdens.add(new ElementMagRangeDescription("Carrizo",
					ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, mag, 10d));
			
			makeMultiRecurrPlots(writeDir, display, mag, events, rupIdens);
		}
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		rupIdens.add(new ElementMagRangeDescription("SAF Cholame 7+",
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, 7d, 10d));
		colors.add(Color.RED);
		
		rupIdens.add(new ElementMagRangeDescription("SAF Carrizo 7+",
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d));
		colors.add(Color.BLUE);
		
		rupIdens.add(new ElementMagRangeDescription("Garlock 7+",
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, 7d, 10d));
		colors.add(Color.GREEN);
		
		rupIdens.add(new ElementMagRangeDescription("SAF Mojave 7+",
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d));
		colors.add(Color.BLACK);
		
		rupIdens.add(new ElementMagRangeDescription("SAF Coachella 7+",
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d));
		colors.add(Color.RED);
		
		rupIdens.add(new ElementMagRangeDescription("San Jacinto 7+",
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, 7d, 10d));
		colors.add(Color.CYAN);
		
		ElementMagRangeDescription mojaveCoachellCorupture = new ElementMagRangeDescription("SAF Mojave/Coachella Corupture",
				6d, 10d, ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID);
		
		double quietAftershockBuffer = 1;
		
//		rupIdens.add(new MagRangeRuptureIdentifier(7d, 10d));
//		rupIdenNames.add("All 7+");
//		colors.add(Color.GRAY);
		
		int cholameIndex = 0;
		int carrizoIndex = 1;
		int garlockIndex = 2;
		int mojaveIndex = 3;
		int coachellaIndex = 4;
		int sanJacintoIndex = 5;
		
		List<RuptureIdentifier> rupIdensSubset = Lists.newArrayList();
		List<Color> colorsSubset = Lists.newArrayList();
		rupIdensSubset.add(rupIdens.get(mojaveIndex));
		colorsSubset.add(colors.get(mojaveIndex));
		rupIdensSubset.add(rupIdens.get(coachellaIndex));
		colorsSubset.add(colors.get(coachellaIndex));
		
		List<RuptureIdentifier> rupIdensNoCholame = Lists.newArrayList(rupIdens);
		rupIdensNoCholame.remove(cholameIndex);
		List<Color> colorsNoCholame = Lists.newArrayList(colors);
		colorsNoCholame.remove(cholameIndex);
		
		List<RuptureIdentifier> elemRupIdens = getOnlyElemMagDescriptions(rupIdens);
		
		boolean[] randoms;
		if (doRandom) {
			randoms = new boolean[2];
			randoms[1] = true;
		} else {
			randoms = new boolean [1];
		}
		
		ArrayList<EvenlyDiscretizedFunc> slidingWindows = null;
		ArrayList<EvenlyDiscretizedFunc> randomizedSlidingWindows = null;
		
		List<SimulatorEvent> randomResampledCatalog = null;
		for (boolean randomized : randoms)
			if (randomized)
				randomResampledCatalog = RandomCatalogBuilder.getRandomResampledCatalog(events, elemRupIdens, randDistType, randSplitMults);
		
		for (boolean randomized : randoms) {
			if (randomized)
				events = randomResampledCatalog;
			
			File myWriteDir;
			if (randomized)
				myWriteDir = new File(writeDir, "randomized");
			else
				myWriteDir = writeDir;
			if (!myWriteDir.exists())
				myWriteDir.mkdir();
			
			plotPeriodsAndEvents(events, display, displayEventTimes, myWriteDir,
					rupIdensSubset, colorsSubset, randomized);
			plotPeriodsAndEvents(events, display, displayEventTimes, myWriteDir,
					rupIdensNoCholame, colorsNoCholame, randomized);
			
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(coachellaIndex),
					rupIdens.get(carrizoIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(coachellaIndex),
					rupIdens.get(mojaveIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(carrizoIndex),
					rupIdens.get(coachellaIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(carrizoIndex),
					rupIdens.get(mojaveIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(mojaveIndex),
					rupIdens.get(carrizoIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(cholameIndex),
					rupIdens.get(carrizoIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(carrizoIndex),
					rupIdens.get(cholameIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(mojaveIndex),
					rupIdens.get(garlockIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(garlockIndex),
					rupIdens.get(mojaveIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(sanJacintoIndex),
					rupIdens.get(mojaveIndex));
			plotTimeBetweenIdens(myWriteDir, display, randomized, events, rupIdens.get(mojaveIndex),
					rupIdens.get(coachellaIndex));
			
			if (!randomized) {
				Map<IDPairing, HistogramFunction[]> origFuncs =
						plotACDF_CCDFs(myWriteDir, events, rupIdens, null, null, 2000d, 10d);
				for (RandomDistType randDist : RandomDistType.values())
					plotACDF_CCDFs(myWriteDir, events, rupIdens, randDist, origFuncs, 2000d, 10d);
				// zoomed in
				File zoomWriteDir = new File(myWriteDir, "corr_zoomed");
				if (!zoomWriteDir.exists())
					zoomWriteDir.mkdir();
				plotACDF_CCDFs(zoomWriteDir, events, rupIdens, null, null, 20d, 1d);
				// zoomed out
				zoomWriteDir = new File(myWriteDir, "corr_wide");
				if (!zoomWriteDir.exists())
					zoomWriteDir.mkdir();
				double totYears = events.get(events.size()-1).getTimeInYears()-events.get(0).getTimeInYears();
				// make it round
				totYears = Math.ceil(totYears / 1000d) * 1000d;
				plotACDF_CCDFs(zoomWriteDir, events, rupIdens, null, null, totYears, 10d);
			}
			
//			double[] windowLengths = { 5d, 10d, 25d, 50d, 100d };
			double[] windowLengths = new double[30];
			for (int i=1; i<=windowLengths.length; i++)
				windowLengths[i-1] = 5d*(double)i;
			double xInc = 1d;
			
			if (randomized)
				randomizedSlidingWindows =
					plotSlidingWindowCounts(myWriteDir, display, randomized, windowLengths, xInc, events, elemRupIdens);
			else
				slidingWindows =
					plotSlidingWindowCounts(myWriteDir, display, randomized, windowLengths, xInc, events, elemRupIdens);
			
			plotInterEventBetweenAllDist(myWriteDir, display, randomized, events, elemRupIdens);
			
			boolean[] initials = {true, false};
			double cumulativePlotYears = 1000d;
			if (!randomized) {
				if (randomResampledCatalog == null)
					randomResampledCatalog = RandomCatalogBuilder.getRandomResampledCatalog(events, elemRupIdens, randDistType, randSplitMults);
				
				for (boolean includeInitialCorupture : initials) {
					
					double[] quiets;
					if (includeInitialCorupture)
						quiets = toDoubleArray(-1d);
					else
						quiets = toDoubleArray(-1d, 156);
					for (double quiet : quiets) {
						if (quiet > 0) {
							RuptureIdentifier quietIden = new QuietPeriodIdenMatcher(rupIdens.get(mojaveIndex), quietAftershockBuffer, quiet, rupIdens);
							plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
									// target
									rupIdens.get(coachellaIndex),
									// given
									quietIden,
									cumulativePlotYears, includeInitialCorupture);
							quietIden = new QuietPeriodIdenMatcher(rupIdens.get(coachellaIndex), quietAftershockBuffer, quiet, rupIdens);
							plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
									// target
									rupIdens.get(mojaveIndex),
									// given
									rupIdens.get(coachellaIndex),
									cumulativePlotYears, includeInitialCorupture);
						} else {
							plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
									// target
									rupIdens.get(coachellaIndex),
									// given
									rupIdens.get(mojaveIndex),
									cumulativePlotYears, includeInitialCorupture);
							plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
									// target
									rupIdens.get(mojaveIndex),
									// given
									rupIdens.get(coachellaIndex),
									cumulativePlotYears, includeInitialCorupture);
						}
						
//						plotConditionalProbs(myWriteDir, display, randomized, events,
//								// target
//								rupIdens.get(coachellaIndex), rupIdenNames.get(coachellaIndex),
//								// given
//								rupIdens.get(mojaveIndex), rupIdenNames.get(mojaveIndex), quiet, 5d, includeInitialCorupture);
					}
				}
				RuptureIdentifier quietIden = new QuietPeriodIdenMatcher(mojaveCoachellCorupture, quietAftershockBuffer, 156, rupIdens);
				plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
						// target
						rupIdens.get(mojaveIndex),
						// given
						quietIden,
						cumulativePlotYears, true);
				plotConditionalProbs(myWriteDir, display, events, randomResampledCatalog,
						// target
						rupIdens.get(coachellaIndex),
						// given
						quietIden,
						cumulativePlotYears, true);
			}
			
			if (!randomized) {
				double[] omoris = { 6d, 6.5d, 7d, 7.5d };
				
				for (int i=0; i<elemRupIdens.size(); i++) {
					RuptureIdentifier rupIden = elemRupIdens.get(i);
					plotOmoriDecay(myWriteDir, display, randomized, events, rupIden, omoris, 365);
				}
			}
		}
		
		if (randomizedSlidingWindows != null && slidingWindows != null) {
			ArrayList<EvenlyDiscretizedFunc> ratios = Lists.newArrayList();
			
			CPT cpt = null;
			try {
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, slidingWindows.size()-1);
			} catch (IOException e1) {
				ExceptionUtils.throwAsRuntimeException(e1);
			}
			
			double min = 0;
			int num = slidingWindows.get(0).size();
			double delta = 1d;
			
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			for (int i=0; i<slidingWindows.size(); i++) {
				EvenlyDiscretizedFunc ratio = new EvenlyDiscretizedFunc(min, num, delta);
				
				EvenlyDiscretizedFunc slidingFunc = slidingWindows.get(i);
				EvenlyDiscretizedFunc randomFunc = randomizedSlidingWindows.get(i);
				
				for (int x=0; x<num; x++) {
					ratio.set(x, slidingFunc.getY(x) / randomFunc.getY(x));
				}
				
				ratios.add(ratio);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cpt.getColor((float)i)));
			}
			
			makePlot(writeDir, "matches_in_windows_ratio", display, false, ratios, chars,
					"Ratio of # Matching Events In Windows (catalog / randomized)",
					"# Events (M7+ on specific faults)",  "Ratio (catalog / randomized)", null, null);
		}
	}
	
	private static Color duplicateColor = Color.ORANGE.darker();

	public static void plotPeriodsAndEvents(List<? extends SimulatorEvent> events,
			boolean display, boolean displayEventTimes, File writeDir,
			List<RuptureIdentifier> rupIdens,
			List<Color> colors, boolean randomized) throws IOException {
		int[] dimensions = { 1000, 500 };
		
		HashSet<Integer> duplicates = new HashSet<Integer>();
		Map<Integer, SimulatorEvent> eventsMap = Maps.newHashMap();
		List<HashSet<Integer>> idsList = Lists.newArrayList();
		
		ArrayList<DiscretizedFunc> allPeriodFuncs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> allPeriodChars = Lists.newArrayList();
		
		ArrayList<ArbitrarilyDiscretizedFunc> labelFuncs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> labelChars = Lists.newArrayList();
		
		for (int i=0; i<rupIdens.size(); i++) {
			RuptureIdentifier rupIden = rupIdens.get(i);
			String name = rupIden.getName();
			
			HistogramFunction hist;
			if (rupIden instanceof ElementMagRangeDescription)
				hist = new HistogramFunction(5d, 100, 10d);
			else
				hist = new HistogramFunction(5d, 100, 1d);
			hist.setName(name);
			
			HashSet<Integer> ids = new HashSet<Integer>();
			idsList.add(ids);
			
			double prevTime = -1;
			
			ArbitrarilyDiscretizedFunc labelFunc = new ArbitrarilyDiscretizedFunc(name+" Rups");
			
			double labelY = (rupIdens.size()-1) - i + 1;
			
			ArrayList<Double> rpList = Lists.newArrayList();
			
			for (SimulatorEvent match : rupIden.getMatches(events)) {
				double time = match.getTimeInYears();
				labelFunc.set(time, labelY);
				Integer id = match.getID();
				if (eventsMap.containsKey(id))
					duplicates.add(id);
				else
					eventsMap.put(id, match);
				ids.add(match.getID());
				
				if (prevTime > 0) {
					double diff = time - prevTime;
					rpList.add(diff);
					int ind = hist.getXIndex(diff);
					if (ind >= 0)
						hist.add(ind, 1d);
				}
				prevTime = time;
			}
			
			double[] rps = Doubles.toArray(rpList);
			Arrays.sort(rps);
			double mean = StatUtils.mean(rps);
			double median = BatchPlotGen.median(rps);
			
			double[] deviations = new double[rps.length];
			for (int j=0; j<rps.length; j++)
				deviations[j] = Math.abs(rps[j] - median);
			Arrays.sort(deviations);
			double mad = BatchPlotGen.median(deviations);
			
			System.out.println(name+": mean="+mean+"\tmedian="+median+"\tmad="+mad);
			
			labelFuncs.add(labelFunc);
			labelChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 8f, colors.get(i)));
			
			ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(hist);
			allPeriodFuncs.add(hist);
//			funcs.add(getCmlGreaterOrEqual(hist));
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, colors.get(i)));
			allPeriodChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i)));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
			
			makePlot(writeDir, "period_"+getFileSafeString(name), display, randomized,
					funcs, chars, name+" Inter-Event Times", "Years", "Number", null, null, dimensions, true);
		}
		
		double[][] periodRanges = new double[2][2];
		periodRanges[0][0] = 0;
		periodRanges[0][1] = 1000;
		periodRanges[1][0] = 0.9;
		periodRanges[1][1] = 2000;
		boolean[] periodLogs = { false, true };
		makePlot(writeDir, "period_all_log", display, randomized,
				allPeriodFuncs, allPeriodChars, " Inter-Event Times", "Years", "Number", periodRanges, periodLogs, dimensions, true);
		
		periodLogs[1] = false;
		periodRanges[1][0] = 0;
		makePlot(writeDir, "period_all", display, randomized,
				allPeriodFuncs, allPeriodChars, " Inter-Event Times", "Years", "Number", periodRanges, periodLogs, dimensions, true);
		
		
		ArrayList<ArbitrarilyDiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		ArbitrarilyDiscretizedFunc duplicateFunc = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<rupIdens.size(); i++) {
			ArbitrarilyDiscretizedFunc overlayFunc = new ArbitrarilyDiscretizedFunc();
			HashSet<Integer> ids = idsList.get(i);
			for (Integer id : ids) {
				SimulatorEvent event = eventsMap.get(id);
				if (duplicates.contains(id))
					duplicates.add(id);
				else
					overlayFunc.set(event.getTimeInYears(), event.getMagnitude());
			}
			overlayFunc.setName(rupIdens.get(i).getName()+" events"); 
			funcs.add(overlayFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, colors.get(i)));
		}
		
		for (Integer id : duplicates) {
			SimulatorEvent event = eventsMap.get(id);
			duplicateFunc.set(event.getTimeInYears(), event.getMagnitude());
		}
		duplicateFunc.setName("Duplicates");
		funcs.add(duplicateFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 3f, duplicateColor));
		
		String prefix = "event_times";
		if (rupIdens.size() < 3)
			prefix += "_subset";
		
		makePlot(writeDir, prefix, displayEventTimes, randomized, funcs, chars, "Event Times");
		
		List<int[]> xRanges = Lists.newArrayList();
		xRanges.add(toIntArray(0, 5000));
		xRanges.add(toIntArray(0, 10000));
		xRanges.add(toIntArray(5000, 10000));
		xRanges.add(toIntArray(200000, 205000));
		xRanges.add(toIntArray(205000, 210000));
		xRanges.add(toIntArray(210000, 215000));
		xRanges.add(toIntArray(215000, 220000));
		xRanges.add(toIntArray(220000, 225000));
		
		colors = Lists.newArrayList(colors);
		colors.add(duplicateColor);
		
		for (int[] xRange : xRanges) {
			double[][] ranges = new double[2][2];
			ranges[0][0] = xRange[0];
			ranges[0][1] = xRange[1];
			ranges[1][0] = 6.75;
			ranges[1][1] = 8;
			
			// need to do this as lines to make it readible, smallest on top
			
			List<HashSet<ArbitrarilyDiscretizedFunc>> funcSets = Lists.newArrayList();
			
			for (DiscretizedFunc func : funcs) {
				HashSet<ArbitrarilyDiscretizedFunc> subFuncs = new HashSet<ArbitrarilyDiscretizedFunc>();
				
				for (Point2D pt : func) {
					if (pt.getX() < ranges[0][0])
						continue;
					if (pt.getX() > ranges[0][1])
						break;
					ArbitrarilyDiscretizedFunc subFunc = new ArbitrarilyDiscretizedFunc();
					subFunc.set(pt.getX(), 6d);
					subFunc.set(pt.getX()+1e-10, pt.getY());
					subFuncs.add(subFunc);
				}
				
				funcSets.add(subFuncs);
			}
			
			ArrayList<ArbitrarilyDiscretizedFunc> subFuncs = Lists.newArrayList();
			ArrayList<PlotCurveCharacterstics> subChars = Lists.newArrayList();
			
			for (HashSet<ArbitrarilyDiscretizedFunc> set : funcSets)
				subFuncs.addAll(set);
			
			Collections.sort(subFuncs, new Comparator<ArbitrarilyDiscretizedFunc>() {

				@Override
				public int compare(ArbitrarilyDiscretizedFunc o1,
						ArbitrarilyDiscretizedFunc o2) {
					double y1 = o1.getY(1);
					double y2 = o2.getY(1);
					return -Double.compare(y1, y2);
				}
			});
			
			funcLoop:
			for (ArbitrarilyDiscretizedFunc func : subFuncs) {
				for (int i=0; i<funcSets.size(); i++) {
					HashSet<ArbitrarilyDiscretizedFunc> set = funcSets.get(i);
					if (set.contains(func)) {
						subChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors.get(i)));
						continue funcLoop;
					}
				}
				throw new IllegalStateException("Func not found in any set!");
			}
			
			makePlot(writeDir, prefix+"_"+xRange[0]+"_"+xRange[1], displayEventTimes, randomized,
					subFuncs, subChars, "Event Times", "Years", "Magnitude", ranges, null, dimensions, true);
			
			double[][] labelRange = new double[2][2];
			labelRange[0][0] = xRange[0];
			labelRange[0][1] = xRange[1];
			labelRange[1][0] = 0.5;
			labelRange[1][1] = labelFuncs.size()+0.5;
			
			int[] labelDims = { 1000, 300 };
			
			makePlot(writeDir, prefix+"_labels_"+xRange[0]+"_"+xRange[1], displayEventTimes, randomized,
					labelFuncs, labelChars, "Event Times", "Years", "Magnitude", labelRange, null, labelDims, true);
		}
	}
	
	private static List<RuptureIdentifier> getOnlyElemMagDescriptions(List<RuptureIdentifier> rupIdens) {
		List<RuptureIdentifier> ret = Lists.newArrayList();
		
		for (RuptureIdentifier rupIden : rupIdens)
			if (rupIden instanceof ElementMagRangeDescription)
				ret.add(rupIden);
		
		return ret;
	}
	
	private static void plotTimeBetweenIdens(File writeDir, boolean display, boolean randomized,
			List<? extends SimulatorEvent> events, RuptureIdentifier iden1,
			RuptureIdentifier iden2) throws IOException {
		String name1 = iden1.getName();
		String name2 = iden2.getName();
		double maxX = 500d;
		HistogramFunction hist = new HistogramFunction(5d, 50, 10d);
		HistogramFunction absHist = new HistogramFunction(5d, 50, 10d);
		
		List<? extends SimulatorEvent> matches1 = iden1.getMatches(events);
		List<? extends SimulatorEvent> matches2 = iden2.getMatches(events);
		HashSet<SimulatorEvent> matches1set = new HashSet<SimulatorEvent>(matches1);
		HashSet<SimulatorEvent> matches2set = new HashSet<SimulatorEvent>(matches2);
		
		int numWithin5years = 0;
		int numAfterWithin5years = 0;
		int totalNum = 0;
		int numCoruptures = 0;
		
		for (SimulatorEvent event1 : matches1) {
			double timeYears1 = event1.getTimeInYears();
			double waitingTime = -1;
			double absMin = Double.MAX_VALUE;
			for (SimulatorEvent event2 : matches2) {
				double timeYears2 = event2.getTimeInYears();
				double absWaitingTime = Math.abs(timeYears2 - timeYears1);
				if (absWaitingTime < absMin)
					absMin = absWaitingTime;
				if (timeYears2 < timeYears1)
					continue;
				waitingTime = timeYears2 - timeYears1;
				break;
			}
			
			if (waitingTime >= 0 && waitingTime <= maxX) {
				hist.add(waitingTime, 1d);
			}
			
			if (absMin >= 0 && absMin <= maxX) {
				absHist.add(absMin, 1d);
			}
			
			totalNum++;
			if (absMin <= 10d)
				numWithin5years++;
			if (absMin == 0d)
				numCoruptures++;
		}
		
		if (!randomized) {
			double percentWithin5 = (double)numWithin5years / (double)totalNum;
			percentWithin5 *= 100;
			System.out.println(name1+" to "+name2+": "+(float)percentWithin5+" % within 5 years");
			double percentCorupture = (double)numCoruptures / (double)totalNum;
			percentCorupture *= 100;
			System.out.println(name1+" to "+name2+": "+(float)percentCorupture+" % co-ruptures");
			double percentWithin5NonCorupture = (double)(numWithin5years - numCoruptures) / (double)totalNum;
			percentWithin5NonCorupture *= 100;
			System.out.println(name1+" to "+name2+": "+(float)percentWithin5NonCorupture+" % within 5 years no co-rupture");
//			double timeIndepExpectedWithin = 100d * ((double)(matches1.size() - numCoruptures) * 10d
//					/ General_EQSIM_Tools.getSimulationDurationYears(events));
//			System.out.println(name1+" to "+name2+": "+(float)timeIndepExpectedWithin+" % within 5 years no co-rupture (RANDOM)");
		}
		
		double[][] ranges = new double[2][2];
		ranges[0][0] = 0;
		ranges[0][1] = maxX;
		ranges[1][0] = 0;
		ranges[1][1] = 3500;
		
		String plotTitle = "Inter-Event Times From "+name1+" to "+name2;
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(hist);
//		funcs.add(getCmlGreaterOrEqual(hist));
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, duplicateColor));
		makePlot(writeDir, "inter_event_"+getFileSafeString(name1)+"_to_"+getFileSafeString(name2),
				display, randomized, funcs, chars, plotTitle, null, null, ranges, null);
		
		plotTitle = "Absolute Inter-Event Times From "+name1+" to "+name2;
		funcs = Lists.newArrayList();
		funcs.add(absHist);
		makePlot(writeDir, "inter_event_abs_"+getFileSafeString(name1)+"_to_"+getFileSafeString(name2),
				display, randomized, funcs, chars, plotTitle, null, null, ranges, null);
		
		HistogramFunction matrixHist = new HistogramFunction(-1000, 200, 10d);
		double matHistMin = matrixHist.getMinX() - 5d;
		double matHistMax = matrixHist.getMaxX() + 5d;
		HistogramFunction absMatrixHist = new HistogramFunction(5d, 100, 10d);
		HistogramFunction absBothMatrixHist = new HistogramFunction(5d, 100, 10d);
		
		for (SimulatorEvent event1 : matches1) {
			double timeYears1 = event1.getTimeInYears();
			for (SimulatorEvent event2 : matches2) {
				double timeYears2 = event2.getTimeInYears();
				
				double timeDelta = timeYears1 - timeYears2;
				double absDelta = Math.abs(timeDelta);
				
				if (timeDelta >= matHistMin && timeDelta <= matHistMax)
					matrixHist.add(timeDelta, 1d);
				if (absDelta <= matHistMax) {
					absMatrixHist.add(absDelta, 1d);
					if (matches1set.contains(event2) && matches2set.contains(event1))
						absBothMatrixHist.add(absDelta, 1d);
				}
			}
		}
		
		double[][] allRanges = new double[2][2];
		
		allRanges[0][0] = 0;
		allRanges[0][1] = 1000;
		allRanges[1][0] = 0;
		allRanges[1][1] = absMatrixHist.getMaxY()*1.1d;
		
		plotTitle = "Inter-Event Times From All "+name1+" to All "+name2;
		funcs = Lists.newArrayList();
		funcs.add(matrixHist);
//		funcs.add(getCmlGreaterOrEqual(hist));
		chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, Color.BLACK));
		makePlot(writeDir, "all_inter_event_"+getFileSafeString(name1)+"_to_"+getFileSafeString(name2),
				display, randomized, funcs, chars, plotTitle, "Years", "Number", null, null);
		
		plotTitle = "Absolute Inter-Event Times From All "+name1+" to All "+name2;
		funcs = Lists.newArrayList();
		funcs.add(absMatrixHist);
		funcs.add(absBothMatrixHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, duplicateColor));
		makePlot(writeDir, "all_inter_event_abs_"+getFileSafeString(name1)+"_to_"+getFileSafeString(name2),
				display, randomized, funcs, chars, plotTitle, "Years", "Number", allRanges, null);
	}
	
	private static PlotSpec getCDF_Func(int ind1, String name1, List<? extends SimulatorEvent> matches1,
			int ind2, String name2, List<? extends SimulatorEvent> matches2, double plotYears, double binWidth) {
		HistogramFunction matrixHist = new HistogramFunction(-plotYears-0.5*binWidth, (int)(plotYears*2d/binWidth)+1, binWidth);
		HistogramFunction corupHist = new HistogramFunction(-plotYears, (int)(plotYears*2d/binWidth), binWidth);
		double matHistMin = matrixHist.getMinX() - 0.5*binWidth;
		double matHistMax = matrixHist.getMaxX() + 0.5*binWidth;
		
		for (SimulatorEvent event1 : matches1) {
			double timeYears1 = event1.getTimeInYears();
			for (SimulatorEvent event2 : matches2) {
				double timeYears2 = event2.getTimeInYears();
				
				double timeDelta = timeYears2 - timeYears1;
				
				if (event1.getID() == event2.getID())
					corupHist.add(0d, 1d);
				else if (timeDelta > matHistMin && timeDelta < matHistMax)
					matrixHist.add(timeDelta, 1d);
			}
		}
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(corupHist);
		funcs.add(matrixHist);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		String title = null;
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Event Time Difference (years)", "Number");
		List<XYTextAnnotation> annotations = Lists.newArrayList();
		double annY;
		if (matrixHist.getMaxY() > corupHist.getMaxY())
			annY = matrixHist.getMaxY()*0.95;
		else
			annY = corupHist.getMaxY()*0.95;
		double annX = matrixHist.getMaxX()*0.9;
		Font font = new Font(Font.SERIF, Font.PLAIN, 14);
		XYTextAnnotation leftAnn = new XYTextAnnotation(name1+" TO EACH "+name2, -annX, annY);
		leftAnn.setFont(font);
		leftAnn.setTextAnchor(TextAnchor.TOP_LEFT);
		annotations.add(leftAnn);
//		XYTextAnnotation rightAnn = new XYTextAnnotation(name2+" => "+name1, annX, annY);
//		rightAnn.setFont(font);
//		rightAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
//		annotations.add(rightAnn);
		spec.setPlotAnnotations(annotations);
		return spec;
	}
	
	public static Map<IDPairing, HistogramFunction[]> plotACDF_CCDFs(File writeDir,
			List<? extends SimulatorEvent> events, List<RuptureIdentifier> idens,
			RandomDistType randDistType, Map<IDPairing, HistogramFunction[]> origCorrHists,
			double plotYears, double binWidth)
					throws IOException {
		if (randDistType != null) {
			events = RandomCatalogBuilder.getRandomResampledCatalog(events, idens, randDistType, true);
			writeDir = new File(writeDir, randDistType.getFNameAdd()+"_corr_plots");
			if (!writeDir.exists())
				writeDir.mkdir();
		}
		
		List<List<? extends SimulatorEvent>> matchesList = Lists.newArrayList();
		for (RuptureIdentifier iden : idens)
			matchesList.add(iden.getMatches(events));
		
		String xAxisLabel = "Event Time Difference (years)";
		String yAxisLabel = "Number";
		
		Map<IDPairing, HistogramFunction> corrHists = Maps.newHashMap();
		Map<IDPairing, HistogramFunction> corupHists = Maps.newHashMap();
		Map<IDPairing, HistogramFunction[]> combHists = Maps.newHashMap();
		
		File acdfDir = new File(writeDir, "acdfs");
		if (!acdfDir.exists())
			acdfDir.mkdir();
		File ccdfDir = new File(writeDir, "ccdfs");
		if (!ccdfDir.exists())
			ccdfDir.mkdir();
		
		List<PlotSpec> specs = Lists.newArrayList();
		for (int i=0; i<idens.size(); i++) {
			String name1 = idens.get(i).getName();
			List<? extends SimulatorEvent> matches1 = matchesList.get(i);
			for (int j=i; j<idens.size(); j++) {
				String name2 = idens.get(j).getName();
				List<? extends SimulatorEvent> matches2 = matchesList.get(j);
				
				PlotSpec spec = getCDF_Func(i, name1, matches1, j, name2, matches2, plotYears, binWidth);
				
				HistogramFunction matrixHist = (HistogramFunction)spec.getPlotElems().get(1);
				HistogramFunction corupHist = (HistogramFunction)spec.getPlotElems().get(0);
				
				IDPairing pair = new IDPairing(i, j);
				
				HistogramFunction[] comb = { corupHist, matrixHist };
				combHists.put(pair, comb);
				
				corrHists.put(pair, matrixHist);
				corupHists.put(pair, corupHist);
				
				specs.add(spec);
				
				HistogramFunction randMatrixHist = null;
				PlotCurveCharacterstics corupChar =
						new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY);
				PlotCurveCharacterstics cdfChar =
						new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK);
				PlotCurveCharacterstics randChar =
						new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLUE);
				
				// write out individual plots
				File myWriteDir;
				String fName;
				String title;
				int ySize;
				if (i == j) {
					myWriteDir = acdfDir;
					fName = "acdf_"+getFileSafeString(name1);
					title = "ACDF: "+name1;
					// skip corup hist
					corupHist = null;
					ySize = 250;
				} else {
					myWriteDir = ccdfDir;
					fName = "ccdf_"+getFileSafeString(name1)+"_"+getFileSafeString(name2);
					title = "CCDF: "+name2+" - "+name1;
					ySize = 400;
				}
				if (origCorrHists != null && randDistType != null) {
					// this is a random one, show actual in bold
					randMatrixHist = matrixHist;
					matrixHist = origCorrHists.get(pair)[1];
					corupHist = origCorrHists.get(pair)[0];
				}
				List<PlotElement> elems = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				if (corupHist != null) {
					elems.add(corupHist);
					chars.add(corupChar);
				}
				if (randMatrixHist != null) {
					elems.add(randMatrixHist);
					chars.add(randChar);
				}
				elems.add(matrixHist);
				chars.add(cdfChar);
				PlotSpec newSpec = new PlotSpec(elems, chars, title, xAxisLabel, yAxisLabel);
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.drawGraphPanel(newSpec, false, false, null, null);
				gp.getChartPanel().setSize(1500, ySize);
				String fileName = myWriteDir.getAbsolutePath()+File.separator+fName;
				gp.saveAsPDF(fileName+".pdf");
			}
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setCombinedOnYAxis(false);
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.drawGraphPanel(specs, false, false, null, null);
		gp.getChartPanel().setSize(2000, 5000);
		String fileName = writeDir.getAbsolutePath()+File.separator+"comb_all_inter_events";
		gp.saveAsPNG(fileName+".png");
		gp.saveAsPDF(fileName+".pdf");
		
		Map<IDPairing, HistogramFunction> corrFirstHists = Maps.newHashMap();
		
		specs = Lists.newArrayList();
		boolean[] swaps = { false, true };
		for (int i=0; i<idens.size(); i++) {
			for (int j=i+1; j<idens.size(); j++) {
				
				for (boolean swap : swaps) {
					String name1, name2;
					List<? extends SimulatorEvent> matches1, matches2;
					IDPairing pair;
					if (swap) {
						name1 = idens.get(j).getName();
						matches1 = matchesList.get(j);
						name2 = idens.get(i).getName();
						matches2 = matchesList.get(i);
						pair = new IDPairing(j, i);
					} else {
						name1 = idens.get(i).getName();
						matches1 = matchesList.get(i);
						name2 = idens.get(j).getName();
						matches2 = matchesList.get(j);
						pair = new IDPairing(i, j);
					}
					
					double halfPlotYears = plotYears*0.5;
					HistogramFunction matrixHist = new HistogramFunction(-halfPlotYears-0.5*binWidth, (int)(halfPlotYears*2d/binWidth)+1, binWidth);
					HistogramFunction corupHist = new HistogramFunction(-halfPlotYears, (int)(halfPlotYears*2d/binWidth), binWidth);
					double matHistMin = matrixHist.getMinX() - 0.5*binWidth;
					double matHistMax = matrixHist.getMaxX() + 0.5*binWidth;
					
					int startK = 0;
					for (SimulatorEvent event1 : matches1) {
						double timeYears1 = event1.getTimeInYears();
						for (int k=startK; k<matches2.size(); k++) {
							SimulatorEvent event2 = matches2.get(k);
							double timeYears2 = event2.getTimeInYears();
							if (timeYears2 < timeYears1)
								continue;
							
//							startK = k;
							double afterDelta = timeYears2 - timeYears1;
							
							if (event1.getID() == event2.getID()) {
								corupHist.add(0d, 1d);
							} else {
								if (afterDelta >= matHistMin && afterDelta <= matHistMax)
									matrixHist.add(afterDelta, 1d);
								if (k > 0) {
									double beforeDelta = matches2.get(k-1).getTimeInYears()-timeYears1;
									if (beforeDelta >= matHistMin && beforeDelta <= matHistMax)
										matrixHist.add(beforeDelta, 1d);
								}
							}
							break;
						}
					}
					
					corrFirstHists.put(pair, matrixHist);
					
					List<DiscretizedFunc> funcs = Lists.newArrayList();
					funcs.add(corupHist);
					funcs.add(matrixHist);
					List<PlotCurveCharacterstics> chars = Lists.newArrayList(
							new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE),
							new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
					String title = null;
					PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
					List<XYTextAnnotation> annotations = Lists.newArrayList();
					double annY = matrixHist.getMaxY()*0.95;
					double annX = matrixHist.getMaxX()*0.9;
					Font font = new Font(Font.SERIF, Font.PLAIN, 14);
					XYTextAnnotation leftAnn = new XYTextAnnotation(name1+" TO CLOSEST "+name2, -annX, annY);
					leftAnn.setFont(font);
					leftAnn.setTextAnchor(TextAnchor.TOP_LEFT);
					annotations.add(leftAnn);
//					XYTextAnnotation rightAnn = new XYTextAnnotation(name1+" => "+name2, annX, annY);
//					rightAnn.setFont(font);
//					rightAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
//					annotations.add(rightAnn);
					spec.setPlotAnnotations(annotations);
					specs.add(spec);
				}
			}
		}
		
		gp = new HeadlessGraphPanel();
		gp.setCombinedOnYAxis(false);
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.drawGraphPanel(specs, false, false, null, null);
		gp.getChartPanel().setSize(1000, 10000);
		fileName = writeDir.getAbsolutePath()+File.separator+"comb_inter_events";
		gp.saveAsPNG(fileName+".png");
		gp.saveAsPDF(fileName+".pdf");
		
		// now assemble pages into PDFs
		File pdfDir = new File(writeDir, "corr_pdfs");
		if (!pdfDir.exists())
			pdfDir.mkdir();
		List<File> pdfFiles = Lists.newArrayList();
		
		for (int i=0; i<idens.size(); i++) {
			for (int j=i+1; j<idens.size(); j++) {
				HistogramFunction ithAutoFunc = corrHists.get(new IDPairing(i, i));
				HistogramFunction ithCorupFunc = corupHists.get(new IDPairing(i, i));
				HistogramFunction jthAutoFunc = corrHists.get(new IDPairing(j, j));
				HistogramFunction jthCorupFunc = corupHists.get(new IDPairing(j, j));
				
				HistogramFunction crossFunc = corrHists.get(new IDPairing(i, j));
				HistogramFunction crossCorupFunc = corupHists.get(new IDPairing(i, j));
				
				HistogramFunction ijFirstFunc = corrFirstHists.get(new IDPairing(i, j));
				HistogramFunction jiFirstFunc = corrFirstHists.get(new IDPairing(j, i));
				
				// M is our j
				// N is our i
				
				String name1 = idens.get(i).getName();
				String name2 = idens.get(j).getName();
				String title = "m: "+name2+"\nn: "+name1;
				if (randDistType != null)
					title += "\nRANDOMIZED: "+randDistType.getName();
				
				List<PlotCurveCharacterstics> chars = Lists.newArrayList(
						new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY),
						new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				if (origCorrHists != null)
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
				int annotationSize = 18;
				
				specs = Lists.newArrayList();
				
				List<HistogramFunction> funcs = Lists.newArrayList();
				funcs.add(jthCorupFunc);
				funcs.add(jthAutoFunc);
				if (origCorrHists != null)
					funcs.add(origCorrHists.get(new IDPairing(j, j))[1]);
				
				PlotSpec mAutoSpec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				mAutoSpec.setPlotAnnotations(Lists.newArrayList(
						getTopLeftAnnotation(funcs, "m ("+name2+") autocorrelation", annotationSize)));
				specs.add(mAutoSpec);
				
				funcs = Lists.newArrayList();
				funcs.add(ithCorupFunc);
				funcs.add(ithAutoFunc);
				if (origCorrHists != null)
					funcs.add(origCorrHists.get(new IDPairing(i, i))[1]);
				
				PlotSpec nAutoSpec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				nAutoSpec.setPlotAnnotations(Lists.newArrayList(
						getTopLeftAnnotation(funcs, "n ("+name1+") autocorrelation", annotationSize)));
				specs.add(nAutoSpec);
				
				funcs = Lists.newArrayList();
				funcs.add(crossCorupFunc);
				funcs.add(crossFunc);
				if (origCorrHists != null)
					funcs.add(origCorrHists.get(new IDPairing(i, j))[1]);
				
				PlotSpec crossSpec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				crossSpec.setPlotAnnotations(Lists.newArrayList(
						getTopLeftAnnotation(funcs, "m ("+name2+"), n ("+name1+")\ncross-correlation", annotationSize)));
				specs.add(crossSpec);
				
				funcs = Lists.newArrayList();
				funcs.add(crossCorupFunc);
				funcs.add(jiFirstFunc);
				
				PlotSpec mnFirstSpec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				mnFirstSpec.setPlotAnnotations(Lists.newArrayList(
						getTopLeftAnnotation(funcs, "m ("+name2+") to first n ("+name1+")", annotationSize)));
				specs.add(mnFirstSpec);
				
				funcs = Lists.newArrayList();
				funcs.add(crossCorupFunc);
				funcs.add(ijFirstFunc);
				
				PlotSpec nmFirstSpec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				nmFirstSpec.setPlotAnnotations(Lists.newArrayList(
						getTopLeftAnnotation(funcs, "n ("+name1+") to first m ("+name2+")", annotationSize)));
				specs.add(nmFirstSpec);
				
				gp = new HeadlessGraphPanel();
				gp.setCombinedOnYAxis(false);
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.drawGraphPanel(specs, false, false, null, null);
				gp.getChartPanel().setSize(2000, 3000);
				File prefixFile = new File(pdfDir, getFileSafeString(name2)+"_corr_"+getFileSafeString(name1));
				fileName = prefixFile.getAbsolutePath();
				gp.saveAsPNG(fileName+".png");
				gp.saveAsPDF(fileName+".pdf");
				
				pdfFiles.add(new File(fileName+".pdf"));
			}
		}
		
		// now combined pdf
		String fName = "corr_combined";
		if (randDistType != null)
			fName += "_"+randDistType.getFNameAdd();
		File outputFile = new File(writeDir, fName+".pdf");
		combinePDFs(pdfFiles, outputFile);
		
		if (randDistType == null) {
			// now make rolling window autocorrelation plots
			double windowLen = 50000;
			double windowStart = events.get(0).getTimeInYears();
			double windowEnd = windowStart+windowLen;
			double lastEvent = events.get(events.size()-1).getTimeInYears();
			
			List<List<PlotSpec>> rollingSpecs = Lists.newArrayList();
			for (int i=0; i<idens.size(); i++)
				rollingSpecs.add(new ArrayList<PlotSpec>());
			
			while (windowEnd <= lastEvent) {
				List<List<SimulatorEvent>> subMatchesList = Lists.newArrayList();
				for (List<? extends SimulatorEvent> matches : matchesList) {
					List<SimulatorEvent> subMatches = Lists.newArrayList();
					for (SimulatorEvent e : matches) {
						double time = e.getTimeInYears();
						if (time < windowStart)
							continue;
						if (time > windowEnd)
							break;
						subMatches.add(e);
					}
					subMatchesList.add(subMatches);
				}
				
				for (int i=0; i<subMatchesList.size(); i++) {
					List<? extends SimulatorEvent> matches = subMatchesList.get(i);
					String name = idens.get(i).getName();
					PlotSpec spec = getCDF_Func(i, name, matches, i, name, matches, plotYears, binWidth);
					spec.setTitle(name+" Rolling Autocorrelations ("+(int)windowLen+" year binning)");
					rollingSpecs.get(i).add(spec);
				}
				
				windowStart += windowLen;
				windowEnd += windowLen;
			}
			
			// now build PDFs
			pdfFiles = Lists.newArrayList();
			
			pdfDir = new File(writeDir, "rolling_autocorr_pdfs");
			if (!pdfDir.exists())
				pdfDir.mkdir();
			
			for (int i=0; i<idens.size(); i++) {
				String name = idens.get(i).getName();
				specs = rollingSpecs.get(i);
				
				gp = new HeadlessGraphPanel();
				gp.setCombinedOnYAxis(false);
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.drawGraphPanel(specs, false, false, null, null);
				gp.getChartPanel().setSize(2000, 3000);
				File prefixFile = new File(pdfDir, getFileSafeString(name)+"_rolling_autocorr");
				fileName = prefixFile.getAbsolutePath();
				gp.saveAsPNG(fileName+".png");
				gp.saveAsPDF(fileName+".pdf");
				
				pdfFiles.add(new File(fileName+".pdf"));
			}
			
			outputFile = new File(writeDir, "rolling_autocorr_combined.pdf");
			combinePDFs(pdfFiles, outputFile);
		}
		
		// now do the progress plots
		double startTime = events.get(0).getTimeInYears();
		double endTime = events.get(events.size()-1).getTimeInYears();
		
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		List<DiscretizedFunc> progFuncs = Lists.newArrayList();
		
		double delta = 100d;
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		
		for (int i=0; i<idens.size(); i++) {
			HistogramFunction hist = new HistogramFunction(0.5*delta, (int)((endTime-startTime)/delta), delta);
			
			List<? extends SimulatorEvent> matches = matchesList.get(i);
			for (SimulatorEvent e : matches) {
				double t = e.getTimeInYears()-startTime;
				if (t < hist.getMaxX()+0.5*delta)
					hist.add(t, 1d);
			}
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(delta, hist.size(), delta);
			double cnt = 0;
			func.setName(idens.get(i).getName());
			for (int j=0; j<hist.size(); j++) {
				cnt += hist.getY(j);
//				if (i == 0)
//					System.out.println("j="+j+", cnt="+cnt+", hist.getY="+hist.getY(j));
				func.set(j, cnt/(double)matches.size());
			}
			
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colors.get(i%colors.size())));
			progFuncs.add(func);
//			progFuncs.add(hist);
		}
		
		PlotSpec spec = new PlotSpec(progFuncs, chars, "Events Over Time", "Years", "Fraction Of Events");
		gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.drawGraphPanel(spec, false, false, null, null);
		gp.getChartPanel().setSize(1000, 800);
		File prefixFile = new File(writeDir, "events_over_time");
		fileName = prefixFile.getAbsolutePath();
		gp.saveAsPNG(fileName+".png");
		gp.saveAsPDF(fileName+".pdf");
//		gp.saveAsTXT(fileName+".txt");
		
		return combHists;
	}
	
	public static void combinePDFs(List<File> pdfFiles, File outputFile) throws IOException {
		PDDocument document = new PDDocument();
		List<PDDocument> subDocs = Lists.newArrayList();
		for (File pdfFile : pdfFiles) {
			PDDocument part = PDDocument.load(pdfFile);
			List<PDPage> list = part.getDocumentCatalog().getAllPages();
			document.addPage(list.get(0));
			subDocs.add(part);
		}
		try {
			document.save(outputFile);
		} catch (COSVisitorException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		document.close();
		for (PDDocument doc : subDocs)
			doc.close();
	}
	
	private static XYTextAnnotation getTopLeftAnnotation(List<? extends DiscretizedFunc> funcs, String text, int fontSize) {
		MinMaxAveTracker xTrack = new MinMaxAveTracker();
		MinMaxAveTracker yTrack = new MinMaxAveTracker();
		for (DiscretizedFunc func : funcs) {
			xTrack.addValue(func.getMinX());
			xTrack.addValue(func.getMaxX());
			yTrack.addValue(func.getMinY());
			yTrack.addValue(func.getMaxY());
		}
		
		double xPos = xTrack.getMin();
		double yPos = yTrack.getMax();
		
		XYTextAnnotation ann = new XYTextAnnotation(text, xPos, yPos);
		ann.setFont(new Font(Font.SERIF, Font.PLAIN, fontSize));
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		
		return ann;
	}
	
	private static EvenlyDiscretizedFunc getCmlGreaterOrEqual(EvenlyDiscretizedFunc func) {
		EvenlyDiscretizedFunc cml = new EvenlyDiscretizedFunc(func.getMinX(), func.size(), func.getDelta());
		
		double tot = 0d;
		for (int i=func.size(); --i>=0;) {
			tot += func.getY(i);
			cml.set(i, tot);
		}
		
		return cml;
	}
	
	private static ArrayList<EvenlyDiscretizedFunc> plotSlidingWindowCounts(File writeDir, boolean display, boolean randomized,
			double[] windowLengths, double xInc,
			List<? extends SimulatorEvent> events, List<RuptureIdentifier> rupIdens) throws IOException {
		List<List<? extends SimulatorEvent>> eventLists = Lists.newArrayList();
		
		double totEventTime = General_EQSIM_Tools.getSimulationDurationYears(events);
		
		Arrays.sort(windowLengths);
		
		for (RuptureIdentifier rupIden : rupIdens)
			eventLists.add(rupIden.getMatches(events));
		
		double maxWindow = StatUtils.max(windowLengths);
		// tim off half of the max window length from each end of the catalog
		double catalogTrim = maxWindow * 0.51;
		double minEventTime = catalogTrim;
		double maxEventTime = totEventTime - catalogTrim;
		
		int numWindows = (int)((totEventTime - 2d*catalogTrim) / xInc);
		double windowWeight = 1d / (double)numWindows;
		
		int numRupIdens = rupIdens.size();
		
		CPT cpt = null;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(windowLengths[0], windowLengths[windowLengths.length-1]);
		} catch (IOException e1) {
			ExceptionUtils.throwAsRuntimeException(e1);
		}
		ArrayList<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		for (double windowLength : windowLengths) {
			EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, 10, 1d);
			double halfLength = windowLength * 0.5d;
			int[] eventIndexesToStart = new int[numRupIdens];
			for (double x=minEventTime; x<maxEventTime; x+=xInc) {
				int numInWindow = 0;
				double minTime = x-halfLength;
				double maxTime = x+halfLength;
				for (int i=0; i<numRupIdens; i++) {
					List<? extends SimulatorEvent> idenEvents = eventLists.get(i);
					for (int e=eventIndexesToStart[i]; e<idenEvents.size(); e++) {
						double t = idenEvents.get(e).getTimeInYears();
						if (t < minTime) {
							eventIndexesToStart[i] = e+1;
							continue;
						}
						if (t > maxTime)
							break;
						numInWindow++;
					}
				}
				if (numInWindow <= func.getMaxX())
					func.add((double)numInWindow, windowWeight);
			}
			func.setName("Window Length: "+windowLength+"\nFractions:\t"+Joiner.on(",\t").join(func.getYValuesIterator()));
			funcs.add(func);
			Color color = cpt.getColor((float)windowLength);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
		}
		
		makePlot(writeDir, "matches_in_windows", display, randomized, funcs, chars,
				"# Matching Events In Windows", "# Events (M7+ on specific faults)",  "Fraction of windows", null, null);
		
		return funcs;
	}
	
	public static double[] getRPs(List<? extends SimulatorEvent> matches) {
		List<Double> rps = Lists.newArrayList();
		
		double prevTime = -1;
		for (SimulatorEvent e : matches) {
			double time = e.getTimeInYears();
			
			if (prevTime >= 0)
				rps.add(time - prevTime);
			
			prevTime = time;
		}
		
		return Doubles.toArray(rps);
	}
	
	private static void plotInterEventBetweenAllDist(File writeDir, boolean display, boolean randomized,
			List<? extends SimulatorEvent> events, List<RuptureIdentifier> rupIdens) throws IOException {
		HashSet<SimulatorEvent> matchesSet = new HashSet<SimulatorEvent>();
		
		for (RuptureIdentifier rupIden : rupIdens)
			for (SimulatorEvent e : rupIden.getMatches(events))
				matchesSet.add(e);
		
		List<SimulatorEvent> matches = Lists.newArrayList(matchesSet);
		Collections.sort(matches);
		
		HistogramFunction hist = new HistogramFunction(5d, 20, 10d);
		
		double prevTime = matches.get(0).getTimeInYears();
		
		for (int i=1; i<matches.size(); i++) {
			SimulatorEvent e = matches.get(i);
			
			double timeDelta = e.getTimeInYears()-prevTime;
			
			if (timeDelta <= hist.getMaxX()+5d) {
				hist.add(timeDelta, 1d);
			}
			
			prevTime = e.getTimeInYears();
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(hist);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, Color.BLACK));
		
		double[][] ranges = new double[2][2];
		ranges[0][0] = 0;
		ranges[0][1] = hist.getMaxX()+15d;
		ranges[1][0] = 0;
		ranges[1][1] = 6000;
		
		makePlot(writeDir, "inter_any_event_dist", display, randomized, funcs, chars, "Inter Event Time Between Any",
				null, null, ranges, null);
	}
	
	private static void plotOmoriDecay(File writeDir, boolean display, boolean randomized,
			List<? extends SimulatorEvent> events, RuptureIdentifier rupIden, double[] minMags, double maxDays)
					throws IOException {
		
		String idenName = rupIden.getName();
		
		if (idenName.contains("7"))
			idenName = idenName.substring(0, idenName.indexOf("7")).trim();
		
		System.out.println("Plotting Omori "+idenName);
		List<? extends SimulatorEvent> matches = rupIden.getMatches(events);
		
		double maxYears = maxDays / BatchPlotGen.DAYS_PER_YEAR;
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		CPT cpt = null;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, minMags.length-1);
		} catch (IOException e1) {
			ExceptionUtils.throwAsRuntimeException(e1);
		}
		
		for (int m=0; m<minMags.length; m++) {
			int startEventI = 0;
			
			double minMag = minMags[m];
			HistogramFunction hist = new HistogramFunction(0.5d, (int)Math.ceil(maxDays), 1d);
			
			for (SimulatorEvent match : matches) {
				double matchTime = match.getTimeInYears();
				double maxTime = matchTime + maxYears;
				
				for (int i=startEventI; i<events.size(); i++) {
					SimulatorEvent e = events.get(i);
					if (e.getID() == match.getID()) {
						startEventI = i;
						continue;
					}
					if (e.getMagnitude() < minMag)
						continue;
					double eventYears = e.getTimeInYears();
					if (eventYears < matchTime)
						continue;
					if (eventYears > maxTime)
						break;
					
					double deltaDays = (eventYears - matchTime) * BatchPlotGen.DAYS_PER_YEAR;
					if (deltaDays <= hist.getMaxX()+0.5d)
						hist.add(deltaDays, 1);
				}
			}
			
			EvenlyDiscretizedFunc cmlWithOffset = new EvenlyDiscretizedFunc(
					hist.getMinX()+0.5d, hist.getMaxX()+0.5d, hist.size());
			double cnt = 0;
			for (int i=0; i<hist.size(); i++) {
				cnt += hist.getY(i);
				cmlWithOffset.set(i, cnt);
			}
			funcs.add(cmlWithOffset);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, cpt.getColor((float)m)));
			
//			if (idenName.startsWith("SAF Mojave")) {
//				System.out.println("Funcs for M >= "+minMag);
//				for (int i=0; i<cmlWithOffset.getNum(); i++)
//					System.out.println((i+1)+". "+hist.getY(i)+"\t"+cmlWithOffset.getY(i));
//			}
		}
		
//		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
//		funcs.add(hist);
//		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
//		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 5f, Color.BLACK));
		
		String prefix = "omori_"+getFileSafeString(idenName);
		
		double minX = 0.5d;
		double maxX = 0d;
		double minY = Double.POSITIVE_INFINITY;
		double maxY = 0d;
		for (DiscretizedFunc func : funcs) {
			if (func.getMaxX() > maxX)
				maxX = func.getMaxX();
			if (func.getMinY() < minY)
				minY = func.getMinY();
			if (func.getMaxY() > maxY)
				maxY = func.getMaxY();
		}
		if (minY == 0)
			minY = 1d;
		
		double[][] ranges = new double[2][2];
		ranges[0][0] = minX * 0.9;
		ranges[0][1] = maxX * 1.1;
		ranges[1][0] = minY * 0.9;
		ranges[1][1] = maxY * 1.1;
		
		System.out.println(minX+", "+maxX+", "+minY+", "+maxY);
		
		boolean[] logs = { true, true };
		
		makePlot(writeDir, prefix, display, randomized, funcs, chars, "Omori Comparison "+idenName,
				"Day", "Cumulative # Events", ranges, logs);
	}
	
	private static void makeMultiRecurrPlots(File dir, boolean display, double mag,
			List<? extends SimulatorEvent> events, List<RuptureIdentifier> rupIdens)
					throws IOException {
		HashSet<SimulatorEvent> matches = new HashSet<SimulatorEvent>();
		
		for (RuptureIdentifier rupIden : rupIdens)
			matches.addAll(rupIden.getMatches(events));
		
		List<SimulatorEvent> matchesList = Lists.newArrayList(matches);
		Collections.sort(matchesList);
		
		double delta = 2.5d;
		double min = delta * 0.5d;
		int num = (int)(500d / delta);
		HistogramFunction hist = new HistogramFunction(min, num, delta);
		
		double maxDelta = hist.getMaxX() + (delta * 0.5d);
		
		double prevTime = matchesList.get(0).getTimeInYears();
		
		double cnt = 0d;
		for (int i=1; i<matchesList.size(); i++) {
			double time = matchesList.get(i).getTimeInYears();
			double timeDelta = time - prevTime;
			
			if (timeDelta <= maxDelta) {
				hist.add(timeDelta, 1d);
				cnt += 1;
			}
			
			prevTime = time;
		}
		
		hist.scale(1d / cnt);
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(hist);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.RED));
		
		List<String> idenNames = Lists.newArrayList();
		for (RuptureIdentifier rupIden : rupIdens)
			idenNames.add(rupIden.getName());
		
		String title = Joiner.on(", ").join(idenNames)+" M"+(float)mag+"+";
		String xAxisLabel = "Interevent Time (years)";
		String yAxisLabel = null;
		
		String prefix = null;
		for (String idenName : idenNames) {
			if (prefix == null)
				prefix = "";
			else
				prefix += "_";
			prefix += getFileSafeString(idenName);
		}
		prefix = "interevent_"+prefix+"_"+(float)mag+"+";
		
		makePlot(dir, prefix, display, false, funcs, chars, title, xAxisLabel, yAxisLabel,
				null, null);
	}
	
	private static void plotConditionalProbs(File dir, boolean display,
			List<? extends SimulatorEvent> events, List<? extends SimulatorEvent> randomizedEvents, RuptureIdentifier targetIden,
			RuptureIdentifier givenIden, double maxTimeYears, boolean includeInitialCorupture)
					throws IOException {
		ArbitrarilyDiscretizedFunc cumulativeFunc = getCumulativeProbDist(
				events, targetIden, givenIden, maxTimeYears,
				includeInitialCorupture);
		cumulativeFunc.setName("Cumulative Probabilities");
		ArbitrarilyDiscretizedFunc randCumulativeFunc = getCumulativeProbDist(
				randomizedEvents, targetIden, givenIden, maxTimeYears,
				includeInitialCorupture);
		randCumulativeFunc.setName("Cumulative Probabilities (Randomized Catalog)");
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(cumulativeFunc);
		funcs.add(randCumulativeFunc);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
		
		String targetName = targetIden.getName();
		String givenName = givenIden.getName();
		
		String title;
		if (includeInitialCorupture)
			title = "Prob("+targetName+"|"+givenName+") incl Initial Co-rupture";
		else
			title = "Prob("+targetName+"|"+givenName+") excl Initial Co-rupture";
		String xAxisLabel = "Time (years)";
		String yAxisLabel = "Cumulative Probability";
		
		String prefix;
		if (includeInitialCorupture)
			prefix = "cumulative_prob_"+getFileSafeString(targetName)+"_given_"+getFileSafeString(givenName);
		else
			prefix = "cumulative_prob_no_initial_"+getFileSafeString(targetName)+"_given_"+getFileSafeString(givenName);
		
		double[][] ranges = new double[2][2];
		ranges[0][0] = 0d;
		ranges[0][1] = maxTimeYears;
		if (maxTimeYears <= 20) {
			ranges[1][0] = 0d;
			ranges[1][1] = 0.25;
		} else {
			ranges[1][0] = 0d;
			ranges[1][1] = 1d;
		}
		
		makePlot(dir, prefix, display, false, funcs, chars, title, xAxisLabel, yAxisLabel,
				ranges, null);
	}

	public static ArbitrarilyDiscretizedFunc getCumulativeProbDist(
			List<? extends SimulatorEvent> events, RuptureIdentifier targetIden,
			RuptureIdentifier givenIden,double maxTimeYears,
			boolean includeInitialCorupture) {
		List<? extends SimulatorEvent> targetMatches = targetIden.getMatches(events);
		List<? extends SimulatorEvent> givenMatches = givenIden.getMatches(events);
		System.out.println("Target matches: "+targetMatches.size());
		System.out.println("Given matches: "+givenMatches.size());
		
		HashSet<Integer> coruptures = null;
		if (!includeInitialCorupture) {
			coruptures = new HashSet<Integer>();
			for (SimulatorEvent e1 : targetMatches)
				for (SimulatorEvent e2 : givenMatches)
					if (e1.getID() == e2.getID())
						coruptures.add(e1.getID());
		}
		
		ArbitrarilyDiscretizedFunc timeFunc = new ArbitrarilyDiscretizedFunc();
		
		int targetStartIndex = 0;
		
		double yVal;
		if (includeInitialCorupture)
			yVal = 1d/(double)givenMatches.size();
		else
			yVal = 1d/(double)(givenMatches.size() - coruptures.size());
		
		for (SimulatorEvent given : givenMatches) {
			double givenTime = given.getTimeInYears();
			double targetMaxTime = givenTime + maxTimeYears;
			if (!includeInitialCorupture && coruptures.contains(given.getID()))
				continue;
			for (int i=targetStartIndex; i<targetMatches.size(); i++) {
				SimulatorEvent target = targetMatches.get(i);
				double targetTime = target.getTimeInYears();
				if (targetTime < givenTime) {
					targetStartIndex = i;
					continue;
				}
				if (targetTime > targetMaxTime)
					break;
				double deltaTime = targetTime - givenTime;
				if (givenIden instanceof QuietPeriodIdenMatcher && deltaTime < ((QuietPeriodIdenMatcher)givenIden).getQuietYears())
					continue;
				int xIndex = timeFunc.getXIndex(deltaTime);
				if (xIndex < 0)
					timeFunc.set(deltaTime, yVal);
				else
					timeFunc.set(xIndex, yVal+timeFunc.getY(xIndex));
//				if (givenIden instanceof QuietPeriodIdenMatcher)
//					System.out.println("Got a hit after "+deltaTime+" years");
				// we only want the first occurrence as we're doing cumulative probabilities
				break;
			}
		}
		
		ArbitrarilyDiscretizedFunc cumulativeFunc = new ArbitrarilyDiscretizedFunc();
		
		double cumulativeRate = 0;
		for (int i=0; i<timeFunc.size(); i++) {
			double x = timeFunc.getX(i);
			cumulativeRate += timeFunc.getY(i);
			cumulativeFunc.set(x, cumulativeRate);
		}
		return cumulativeFunc;
	}
	
	private static void makePlot(File dir, String prefix, boolean display, boolean randomized,
			ArrayList<? extends DiscretizedFunc> funcs, ArrayList<PlotCurveCharacterstics> chars, String plotTitle)
					throws IOException {
		makePlot(dir, prefix, display, randomized, funcs, chars, plotTitle, null, null, null, null);
	}
	
	private static void makePlot(File dir, String prefix, boolean display, boolean randomized,
			ArrayList<? extends DiscretizedFunc> funcs, ArrayList<PlotCurveCharacterstics> chars, String plotTitle,
			String xAxisLabel, String yAxisLabel, double[][] ranges, boolean[] logs) throws IOException {
		makePlot(dir, prefix, display, randomized, funcs, chars, plotTitle,
				xAxisLabel, yAxisLabel, ranges, logs, null, false);
	}
	
	private static void makePlot(File dir, String prefix, boolean display, boolean randomized,
			ArrayList<? extends DiscretizedFunc> funcs, ArrayList<PlotCurveCharacterstics> chars, String plotTitle,
			String xAxisLabel, String yAxisLabel, double[][] ranges, boolean[] logs, int[] dimensions, boolean legend)
					throws IOException {
		if (randomized) {
			plotTitle = "RANDOMIZED "+plotTitle;
			prefix = prefix+"_randomized";
		}
		
		String fileName = new File(dir, prefix).getAbsolutePath();
		
		PlotSpec spec = new PlotSpec(funcs, chars, plotTitle, xAxisLabel, yAxisLabel);
		if (legend)
			spec.setLegendVisible(true);
		
		if (display) {
			GraphWindow gw = new GraphWindow(spec);
			if (dimensions == null)
				gw.setSize(600, 800);
			else
				gw.setSize(dimensions[0], dimensions[1]);
			if (ranges != null) {
				gw.setX_AxisRange(ranges[0][0], ranges[0][1]);
				gw.setY_AxisRange(ranges[1][0], ranges[1][1]);
			}
			if (logs != null) {
				gw.setXLog(logs[0]);
				gw.setYLog(logs[1]);
			}
			
			gw.getGraphWidget().setBackgroundColor(Color.WHITE);
			gw.setTickLabelFontSize(18);
			gw.setAxisLabelFontSize(20);
			gw.setPlotLabelFontSize(21);
			
			gw.saveAsPNG(fileName+".png");
			gw.saveAsPDF(fileName+".pdf");
		} else {
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			if (ranges != null) {
				gp.setUserBounds(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]);
			}
			if (logs != null) {
				gp.setXLog(logs[0]);
				gp.setYLog(logs[1]);
			}
			gp.setBackgroundColor(Color.WHITE);
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.drawGraphPanel(spec);
			if (dimensions == null)
				gp.getChartPanel().setSize(1000, 800);
			else
				gp.getChartPanel().setSize(dimensions[0], dimensions[1]);
			gp.saveAsPNG(fileName+".png");
			gp.saveAsPDF(fileName+".pdf");
		}
		
	}
	
	public static String getFileSafeString(String str) {
		return str.replaceAll("\\W+", "_");
	}
	
	public static int[] toIntArray(int... ints) {
		return ints;
	}
	
	public static double[] toDoubleArray(double... vals) {
		return vals;
	}
	
	public static <E> E[] toArray(E... vals) {
		return vals;
	}

}
