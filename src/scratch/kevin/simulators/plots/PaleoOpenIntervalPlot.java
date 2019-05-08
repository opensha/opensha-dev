package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.sql.Array;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;

public class PaleoOpenIntervalPlot extends AbstractPlot {
	
	private Table<String, Double, Location> paleoSites;
	private Map<Integer, String> elementMappings;
	private HashSet<SimulatorElement> mappedElements;
	
	private Map<String, SiteHistory> siteHistories;
	private SiteHistory allHistory;
	
	private double annotateDuration;
	
	private static final DiscretizedFunc slipProbFunc;
	
	static {
		// From UCERF3 Appendix I Table 1
		// specific to Wrightwood, but here applied to all. slip in meters
		
		slipProbFunc = new ArbitrarilyDiscretizedFunc();
		slipProbFunc.set(0d, 0d);
		slipProbFunc.set(0.10, 0.05);
		slipProbFunc.set(0.20, 0.25);
		slipProbFunc.set(0.30, 0.50);
		slipProbFunc.set(1d, 0.75);
		slipProbFunc.set(2d, 0.95);
		slipProbFunc.set(4d, 0.99);
	}

	public static Table<String, Double, Location> getSetBiasi2019() {
		Table<String, Double, Location> table = HashBasedTable.create();
		table.put("TYS", 1d/329d, new Location(37.5563, -121.9739));
		table.put("SCZ", 1d/106d, new Location(36.9626, -121.6981));
		table.put("FRA", 1d/119d, new Location(34.8122, -118.9034));
		table.put("HOG", 1d/191d, new Location(33.6153, -116.7091));
		table.put("COA", 1d/181d, new Location(33.7274, -116.1701));
		return table;
	}

	public static Table<String, Double, Location> getSetUCERF3() throws IOException {
		Table<String, Double, Location> table = HashBasedTable.create();
		for (PaleoRateConstraint constraint : UCERF3_PaleoRateConstraintFetcher.getConstraints())
			if (!constraint.getPaleoSiteName().contains("Offshore"))
				table.put(constraint.getPaleoSiteName().trim().replaceAll("\\W+", ""), constraint.getMeanRate(), constraint.getPaleoSiteLoction());
		return table;
	}
	
	public PaleoOpenIntervalPlot(Collection<SimulatorElement> elements, Table<String, Double, Location> paleoSites) {
		this(elements, paleoSites, -1);
	}
	
	public PaleoOpenIntervalPlot(Collection<SimulatorElement> elements, Table<String, Double, Location> paleoSites, double annotateDuration) {
		this.annotateDuration = annotateDuration;
		Preconditions.checkArgument(!paleoSites.isEmpty(), "Must supply at least one paleo site");
		this.paleoSites = paleoSites;
		elementMappings = new HashMap<>();
		
		mappedElements = new HashSet<>();
		
		for (Cell<String, Double, Location> cell : paleoSites.cellSet()) {
			String paleoName = cell.getRowKey();
			Preconditions.checkState(paleoSites.row(paleoName).size() == 1);
			Location paleoLoc = cell.getValue();
			double minDist = Double.POSITIVE_INFINITY;
			SimulatorElement closest = null;
			for (SimulatorElement elem : elements) {
				for (Location elemLoc : elem.getVertices()) {
					double dist = LocationUtils.linearDistanceFast(elemLoc, paleoLoc);
					if (dist < minDist) {
						minDist = dist;
						closest = elem;
					}
				}
			}
			System.out.println("Closest element to "+paleoName+" is "+(float)minDist+" km away");
			Preconditions.checkState(minDist < 5d, "No elements within 5 km of paleo site %s with location %s. closest is %s km away.",
					paleoName, paleoLoc, minDist);
			elementMappings.put(closest.getID(), paleoName);
			mappedElements.add(closest);
		}
		
		siteHistories = new HashMap<>();
		for (String name : paleoSites.rowKeySet())
			siteHistories.put(name, new SiteHistory());
		allHistory = new SiteHistory();
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return mappedElements;
	}
	
	private class SiteHistory {
		private List<Double> occurrenceTimes;
		private List<double[]> slips;
		
		public SiteHistory() {
			this(new ArrayList<>(), new ArrayList<>());
		}
		
		public SiteHistory(List<Double> occurrenceTimes, List<double[]> slips) {
			this.occurrenceTimes = occurrenceTimes;
			this.slips = slips;
		}
		
		public void add(double time, double... slips) {
			occurrenceTimes.add(time);
			this.slips.add(slips);
		}
		
		public List<Double> getIntervals() {
			if (occurrenceTimes.size() < 2)
				return null;
			List<Double> intervals = new ArrayList<>();
			for (int i=1; i<occurrenceTimes.size(); i++)
				intervals.add(occurrenceTimes.get(i)-occurrenceTimes.get(i-1));
			return intervals;
		}
		
		public SiteHistory getRandomSampleWithProbModel() {
			List<Double> newTimes = new ArrayList<>();
			List<double[]> newSlips = new ArrayList<>();
			
			for (int i=0; i<occurrenceTimes.size(); i++) {
				double time = occurrenceTimes.get(i);
				double[] eventSlips = slips.get(i);
				boolean match = false;
				for (double slip : eventSlips) {
					double prob = slip >= slipProbFunc.getMaxX() ? slipProbFunc.getMaxY() : slipProbFunc.getInterpolatedY(slip);
					match = match || Math.random() < prob;
				}
				if (match) {
					newTimes.add(time);
					newSlips.add(eventSlips);
				}
			}
			
			return new SiteHistory(newTimes, newSlips);
		}
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		List<Double> matchingSlips = new ArrayList<>();
		double time = e.getTimeInYears();
		int[] ids = e.getAllElementIDs();
		double[] slips = e.getAllElementSlips();
		for (int i=0; i<ids.length; i++) {
			int id = ids[i];
			double slip = slips[i];
			if (slip == 0d)
				continue;
			String paleoMatch = elementMappings.get(id);
			if (paleoMatch != null) {
				matchingSlips.add(slip);
				siteHistories.get(paleoMatch).add(time, slip);
			}
		}
		if (!matchingSlips.isEmpty())
			allHistory.add(time, Doubles.toArray(matchingSlips));
	}
	
	private class SiteResult {
		
		private final ArbitrarilyDiscretizedFunc incrementalFunc;
		private final ArbitrarilyDiscretizedFunc cumulativeCountFunc;
		private final ArbitrarilyDiscretizedFunc cumulativeProbFunc;
		
		private final double numEvents;
		private final double maxInterval;
		private final double meanRI;
		private final double annualRate;
		
		public SiteResult(List<Double> intervals) {
			if (intervals == null) {
				System.out.println("WARNING: Empty interval list encountered");
				incrementalFunc = null;
				cumulativeCountFunc = new ArbitrarilyDiscretizedFunc();
				cumulativeProbFunc = new ArbitrarilyDiscretizedFunc();
				cumulativeCountFunc.set(0d, 0d);
				cumulativeCountFunc.set(1000d, 0d);
				cumulativeProbFunc.set(0d, 0d);
				cumulativeProbFunc.set(1000d, 0d);
				annualRate = 0d;
				meanRI = getCurrentDurationYears();
				maxInterval = meanRI;
				numEvents = 0;
			} else {
				incrementalFunc = new ArbitrarilyDiscretizedFunc();
				double totTime = 0d;
				double maxInterval = 0d;
				for (double interval : intervals) {
					if (incrementalFunc.hasX(interval))
						incrementalFunc.set(interval, incrementalFunc.getY(interval)+1);
					else
						incrementalFunc.set(interval, 1d);
					totTime += interval;
					maxInterval = Math.max(maxInterval, interval);
				}
				
				Preconditions.checkState(incrementalFunc.size() > 0);
				
				numEvents = intervals.size()+1;
				this.maxInterval = maxInterval;
				meanRI = totTime / (double)intervals.size();
				annualRate = 1d/meanRI;
				
				cumulativeCountFunc = new ArbitrarilyDiscretizedFunc();
				cumulativeProbFunc = new ArbitrarilyDiscretizedFunc();
				
				buildCumulativeFuncs(incrementalFunc, cumulativeCountFunc, cumulativeProbFunc);
			}
		}

		private SiteResult(ArbitrarilyDiscretizedFunc incrementalFunc, ArbitrarilyDiscretizedFunc cumulativeCountFunc,
				ArbitrarilyDiscretizedFunc cumulativeProbFunc, double numEvents, double maxInterval, double meanRI, double annualRate) {
			super();
			this.incrementalFunc = incrementalFunc;
			this.cumulativeCountFunc = cumulativeCountFunc;
			this.cumulativeProbFunc = cumulativeProbFunc;
			this.numEvents = numEvents;
			this.maxInterval = maxInterval;
			this.meanRI = meanRI;
			this.annualRate = annualRate;
		}
		
		public double getProbability(double openInterval) {
			if (openInterval > cumulativeProbFunc.getMaxX())
				return 0d;
			if (openInterval <= cumulativeProbFunc.getMinX())
				return 1d;
			return cumulativeProbFunc.getInterpolatedY(openInterval);
		}
		
		public double getCount(double openInterval) {
			if (openInterval > cumulativeCountFunc.getMaxX())
				return 0d;
			if (openInterval <= cumulativeCountFunc.getMinX())
				return cumulativeCountFunc.getMaxY();
			return cumulativeCountFunc.getInterpolatedY(openInterval);
		}
	}
	
	private void buildCumulativeFuncs(ArbitrarilyDiscretizedFunc incrementalFunc, ArbitrarilyDiscretizedFunc cumulativeCountFunc,
			ArbitrarilyDiscretizedFunc cumulativeProbFunc) {
		double totTime = 0;
		for (Point2D pt : incrementalFunc)
			totTime += pt.getX()*pt.getY();
		
		double intervalsLeft = incrementalFunc.calcSumOfY_Vals();
		double sumTimePrev = 0;
		for (int i=0; i<incrementalFunc.size(); i++) {
			double interval = incrementalFunc.getX(i);
			double num = incrementalFunc.getY(i);
			cumulativeCountFunc.set(interval, intervalsLeft);
			double prob = (totTime - sumTimePrev)/totTime;
			cumulativeProbFunc.set(interval, prob);
			
			sumTimePrev += num*interval;
			intervalsLeft -= num;
		}
	}
	
	private ArbitrarilyDiscretizedFunc averageIncrementalFuncs(List<ArbitrarilyDiscretizedFunc> funcs) {
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		double weight = 1d/funcs.size();
		for (ArbitrarilyDiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (ret.hasX(pt.getX()))
					ret.set(pt.getX(), ret.getY(pt.getX())+pt.getY()*weight);
				else
					ret.set(pt.getX(), pt.getY()*weight);
			}
		}
		Preconditions.checkState(ret.size() > 0);
		return ret;
	}
	
	private SiteResult getRandomSampling(SiteHistory fullHistory, int numSamples) {
		List<ArbitrarilyDiscretizedFunc> incrementalFuncs = new ArrayList<>();
		double numEvents = 0d;
		double annualRate = 0d;
		double maxInterval = 0d;
		double weight = 1d/(double)numSamples;
		for (int i=0; i<numSamples; i++) {
			SiteResult result = new SiteResult(fullHistory.getRandomSampleWithProbModel().getIntervals());
			numEvents += weight*result.numEvents;
			annualRate += weight*result.annualRate;
			maxInterval = Math.max(maxInterval, result.maxInterval);
			incrementalFuncs.add(result.incrementalFunc);
		}
		ArbitrarilyDiscretizedFunc incrementalFunc = averageIncrementalFuncs(incrementalFuncs);
		ArbitrarilyDiscretizedFunc countFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc probFunc = new ArbitrarilyDiscretizedFunc();
		buildCumulativeFuncs(incrementalFunc, countFunc, probFunc);
		return new SiteResult(incrementalFunc, countFunc, probFunc, numEvents, maxInterval, 1d/annualRate, annualRate);
	}
	
	private EvenlyDiscretizedFunc getPoissonFunc(double maxDuration, double annualRate) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, maxDuration, 100);
		for (int i=0; i<func.size(); i++) {
			double duration = func.getX(i);
			func.set(i, getPoissonProb(duration, annualRate));
		}
		return func;
	}
	
	private double getPoissonProb(double duration, double annualRate) {
		return Math.exp(-annualRate*duration);
	}

	@Override
	public void finalizePlot() throws IOException {
		List<Double> allIntervals = allHistory.getIntervals();
		if (allIntervals == null) {
			System.out.println("No intervals found!");
			return;
		}
		
		SiteResult totalResult = new SiteResult(allIntervals);
		SiteResult totalProbModelResult = getRandomSampling(allHistory, 100);
		
		double dataRateAny = 0d;
		CSVFile<String> sitesCSV = new CSVFile<>(true);
		sitesCSV.addLine("Site Name", "Data MRI (yr)", "Data Annual Rate", "Catalog MRI (yr)", "Catalog Annual Rate", "Catalog Occurences",
				"Prob Filtered Catalog MRI (yr)", "Prob Filtered Catalog Annual Rate", "Prob Filtered Catalog Occurences");
		for (String paleoName : siteHistories.keySet()) {
			SiteHistory siteHistory = siteHistories.get(paleoName);
			SiteResult result = new SiteResult(siteHistory.getIntervals());
			SiteResult probResult = getRandomSampling(siteHistory, 100);
			double dataRate = paleoSites.row(paleoName).keySet().iterator().next();
			double dataMRI = 1d/dataRate;
			dataRateAny += dataRate;
			sitesCSV.addLine(paleoName, yearDF.format(dataMRI), (float)dataRate+"",
					yearDF.format(result.meanRI), (float)result.annualRate+"", optionalDigitDF.format(result.numEvents),
					yearDF.format(probResult.meanRI), (float)probResult.annualRate+"", optionalDigitDF.format(probResult.numEvents));
		}
		double dataMRI = 1d/dataRateAny;
		sitesCSV.addLine("TOTAL", yearDF.format(dataMRI), (float)dataRateAny+"",
				yearDF.format(totalResult.meanRI), (float)totalResult.annualRate+"", optionalDigitDF.format(totalResult.numEvents),
				yearDF.format(totalProbModelResult.meanRI), (float)totalProbModelResult.annualRate+"",
				optionalDigitDF.format(totalProbModelResult.numEvents));
		sitesCSV.writeToFile(new File(getOutputDir(), getOutputPrefix()+"_sites.csv"));
		
		System.out.println("Rate any: "+totalResult.annualRate);
		System.out.println("Rate any with prob model: "+totalProbModelResult.annualRate);
		System.out.println("Data rate any: "+dataRateAny);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		totalResult.cumulativeCountFunc.setName("Catalog");
		funcs.add(totalResult.cumulativeCountFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		totalProbModelResult.cumulativeCountFunc.setName("Prob. Filtered Catalog");
		funcs.add(totalProbModelResult.cumulativeCountFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE.darker()));
		
		PlotSpec plot = new PlotSpec(funcs, chars, "Paleo Open Interval Count", "Time (years)", "Cumulative Count");
		plot.setLegendVisible(true);
		
		Range xRange = new Range(0, Math.max(totalProbModelResult.maxInterval*1.1, 150));
		Range yRange = new Range(0, Math.max(totalResult.cumulativeCountFunc.getMaxY(), totalProbModelResult.cumulativeCountFunc.getMaxY())*1.1);
		
		DecimalFormat timeDF = new DecimalFormat("0.0");
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		List<String> annLines = new ArrayList<>();
		annLines.add("Unfiltered Catalog");
		annLines.add("Max: "+timeDF.format(totalResult.maxInterval)+" yr");
		if (annotateDuration > 0)
			annLines.add("# ≥ "+timeDF.format(annotateDuration)+" yr: "+timeDF.format(totalResult.getCount(annotateDuration)));
		annLines.add(" ");
		annLines.add("Prob Detection Filtered");
		annLines.add("Max: "+timeDF.format(totalProbModelResult.maxInterval)+" yr");
		if (annotateDuration > 0)
			annLines.add("# ≥ "+timeDF.format(annotateDuration)+" yr: "+timeDF.format(totalProbModelResult.getCount(annotateDuration)));
		plot.setPlotAnnotations(getAnnotations(xRange.getUpperBound(), yRange.getUpperBound(), annFont, annLines));

		EvenlyDiscretizedFunc dataPoisson = getPoissonFunc(xRange.getUpperBound(), dataRateAny);
		EvenlyDiscretizedFunc totalPoisson = getPoissonFunc(xRange.getUpperBound(), totalResult.annualRate);
		EvenlyDiscretizedFunc totalProbModelPoisson = getPoissonFunc(xRange.getUpperBound(), totalProbModelResult.annualRate);
		
		CSVFile<String> probsCSV = new CSVFile<>(true);
		probsCSV.addLine("Open Interval (yr)", "Catalog Probability", "Catalog Poisson Probability",
				"Prob. Filtered Catalog Probability", "Prob. Filtered Catalog Poisson Probability", "Data Poisson Probability");
		double duration = 0d;
		while (duration < totalProbModelResult.maxInterval) {
			duration += 10;
			double catProb = totalResult.getProbability(duration);
			double catPoissonProb = getPoissonProb(duration, totalResult.annualRate);
			double catProbModelProb = totalProbModelResult.getProbability(duration);
			double catProbModelPoissonProb = getPoissonProb(duration, totalProbModelResult.annualRate);
			double dataProb = Math.exp(-dataRateAny*duration);
			probsCSV.addLine(yearDF.format(duration), (float)catProb+"", (float)catPoissonProb+"",
					(float)catProbModelProb+"", (float)catProbModelPoissonProb+"", (float)dataProb+"");
		}
		probsCSV.writeToFile(new File(getOutputDir(), getOutputPrefix()+"_probs.csv"));
		
		String prefix = getOutputPrefix()+"_count";
		HeadlessGraphPanel gp = getGraphPanel();
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(plot, false, false, xRange, yRange);
		gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		gp.saveAsTXT(new File(getOutputDir(), prefix+".txt").getAbsolutePath());
		gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		
		totalResult.cumulativeProbFunc.setName("Catalog");
		funcs.add(totalResult.cumulativeProbFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		totalPoisson.setName("Poisson");
		funcs.add(totalPoisson);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		totalProbModelResult.cumulativeProbFunc.setName("Prob Filtered");
		funcs.add(totalProbModelResult.cumulativeProbFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE.darker()));
		
		totalProbModelPoisson.setName("Poisson");
		funcs.add(totalProbModelPoisson);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE.darker()));
		
		dataPoisson.setName("Data Poisson");
		funcs.add(dataPoisson);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.GRAY));
		
		plot = new PlotSpec(funcs, chars, "Paleo Open Interval Probability", "Time (years)", "Cumulative Probability");
		plot.setLegendVisible(funcs.size() > 1);
		
		yRange = new Range(0, 1);
		
		annLines = new ArrayList<>();
		annLines.add("Unfiltered Catalog");
		annLines.add("Max: "+timeDF.format(totalResult.maxInterval)+" yr");
		if (annotateDuration > 0) {
			annLines.add("Prob ≥ "+timeDF.format(annotateDuration)+" yr: "
					+timeDF.format(totalResult.getProbability(annotateDuration)*100)+" %");
			annLines.add("Poisson Prob ≥ "+timeDF.format(annotateDuration)+" yr: "
					+timeDF.format(getPoissonProb(annotateDuration, totalResult.annualRate)*100)+" %");
		}
		annLines.add(" ");
		annLines.add("Prob. Detection Filtered Catalog");
		annLines.add("Max: "+timeDF.format(totalProbModelResult.maxInterval)+" yr");
		if (annotateDuration > 0) {
			annLines.add("Prob ≥ "+timeDF.format(annotateDuration)+" yr: "
					+timeDF.format(totalProbModelResult.getProbability(annotateDuration)*100)+" %");
			annLines.add("Poisson Prob ≥ "+timeDF.format(annotateDuration)+" yr: "
					+timeDF.format(getPoissonProb(annotateDuration, totalProbModelResult.annualRate)*100)+" %");
		}
		plot.setPlotAnnotations(getAnnotations(xRange.getUpperBound(), yRange.getUpperBound(), annFont, annLines));
		
		prefix = getOutputPrefix()+"_prob";
		gp = getGraphPanel();
		gp.drawGraphPanel(plot, false, false, xRange, yRange);
		gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		gp.saveAsTXT(new File(getOutputDir(), prefix+".txt").getAbsolutePath());
		gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
	}
	
	private List<XYTextAnnotation> getAnnotations(double maxX, double maxY, Font font, List<String> lines) {
		List<XYTextAnnotation> anns = new ArrayList<>();
		for (int i=0; i<lines.size(); i++) {
			XYTextAnnotation ann = new XYTextAnnotation(lines.get(i)+"  ", maxX, maxY*(0.95 - 0.05*i));
			ann.setFont(font);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			anns.add(ann);
		}
		return anns;
	}

}
