package scratch.kevin.simulators.plots;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.calc.recurInterval.BPT_DistCalc;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class U3StyleNormalizedRuptureRecurrenceIntervalPlot extends AbstractPlot {
	
	private Map<Integer, SimulatorElement> elemsMap;
	private BPTAveragingTypeOptions aveType;
	private double[] minMags;
	private double overallMinMag;
	
	private RSQSimSubSectionMapper mapper;
	private List<? extends FaultSection> subSects;
	private double[] subSectPrevTimes;

	private Map<Double, List<RuptureRecord>> magRecords;
	private Map<Double, int[]> sectEventCounts;
	
	public U3StyleNormalizedRuptureRecurrenceIntervalPlot(List<SimulatorElement> elems, RSQSimSubSectionMapper mapper, double... minMags) {
		this(elems, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE, mapper, minMags);
	}
	
	public U3StyleNormalizedRuptureRecurrenceIntervalPlot(List<SimulatorElement> elems, BPTAveragingTypeOptions aveType,
			RSQSimSubSectionMapper mapper, double... minMags) {
		this.aveType = aveType;
		this.mapper = mapper;
		if (minMags == null || minMags.length == 0)
			minMags = new double[] { 0d };
		this.minMags = minMags;
		this.overallMinMag = StatUtils.min(minMags);

		elemsMap = new HashMap<>();
		for (SimulatorElement e : elems)
			elemsMap.put(e.getID(), e);
		
		subSects = mapper.getSubSections();
		subSectPrevTimes = new double[subSects.size()];
		for (int i=0; i<subSectPrevTimes.length; i++)
			subSectPrevTimes[i] = Double.NaN;
		
		magRecords = new HashMap<>();
		sectEventCounts = new HashMap<>();
		for (double mag : minMags) {
			magRecords.put(mag, new ArrayList<>());
			sectEventCounts.put(mag, new int[subSects.size()]);
		}
	}
	
	private class RuptureRecord {
		private int[] sectIDs;
		private double[] sectOIs;
		
		public RuptureRecord(int[] sectIDs, double[] sectOIs) {
			this.sectIDs = sectIDs;
			this.sectOIs = sectOIs;
		}
	}

	public BPTAveragingTypeOptions getAveType() {
		return aveType;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < overallMinMag)
			return;
		
		List<List<SubSectionMapping>> bundled =  mapper.getFilteredSubSectionMappings((RSQSimEvent)e);
		if (bundled.isEmpty())
			return;
		
		List<Integer> sectIDs = new ArrayList<>();
		List<Double> sectOIs = new ArrayList<>();
		boolean preTimesKnown = true;
		
		double eventTime = e.getTimeInYears();
		
		for (List<SubSectionMapping> bundle : bundled) {
			for (SubSectionMapping mapping : bundle) {
				int id = mapping.getSubSect().getSectionId();
				sectIDs.add(id);
				double prevTime = subSectPrevTimes[id];
				preTimesKnown = preTimesKnown && Double.isFinite(prevTime);
				subSectPrevTimes[id] = eventTime;
				sectOIs.add(eventTime - prevTime);
			}
		}
		
		RuptureRecord record = new RuptureRecord(Ints.toArray(sectIDs), Doubles.toArray(sectOIs));
		
		for (double minMag : minMags) {
			if (mag >= minMag) {
				if (preTimesKnown)
					magRecords.get(minMag).add(record);
				int[] sectCounts = sectEventCounts.get(minMag);
				for (int id : sectIDs)
					sectCounts[id]++;
			}
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		double durationYears = getCurrentDurationYears();
		
		String prefix = getOutputPrefix();
		
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Min Mag");
		header.add("Mean");
		header.add("Standard Deviation");
		header.add("COV");
		csv.addLine(header);
		
		for (double minMag : minMags) {
			int[] sectCounts = sectEventCounts.get(minMag);
			double[] sectRates = new double[sectCounts.length];
			for (int i=0; i<sectCounts.length; i++)
				sectRates[i] = sectCounts[i]/durationYears;
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 4d, 0.1);
			
			SummaryStatistics stats = new SummaryStatistics();
			SummaryStatistics logStats = new SummaryStatistics();
			
			for (RuptureRecord rec : magRecords.get(minMag)) {
				double[] rateQuantities = new double[rec.sectIDs.length];
				double[] timeQuantities = new double[rec.sectIDs.length];
				for (int i=0; i<rec.sectIDs.length; i++) {
					int id = rec.sectIDs[i];
					if (aveType.isAveRI())
						rateQuantities[i] = 1d/sectRates[id];
					else
						rateQuantities[i] = sectRates[id];
					if (aveType.isAveNTS())
						timeQuantities[i] = rec.sectOIs[i]/(1d/sectRates[id]);
					else
						timeQuantities[i] = rec.sectOIs[i];
				}
				double aveRate, aveRI;
				if (aveType.isAveRI()) {
					aveRI = StatUtils.mean(rateQuantities);
					aveRate = 1d/aveRI;
				} else {
					aveRate = StatUtils.mean(rateQuantities);
					aveRI = 1d/aveRate;
				}
				double aveTimeSince, aveNormTimeSince;
				if (aveType.isAveNTS()) {
					aveNormTimeSince = StatUtils.mean(timeQuantities);
					aveTimeSince = aveRI*aveNormTimeSince;
				} else {
					aveTimeSince = StatUtils.mean(timeQuantities);
					aveNormTimeSince = aveTimeSince/aveRI;
				}
				int xInd = hist.getClosestXIndex(aveNormTimeSince);
				hist.add(xInd, 1d);
				stats.addValue(aveNormTimeSince);
				logStats.addValue(Math.log(aveNormTimeSince));
			}
			
			hist.normalizeBySumOfY_Vals();
			
			if (stats.getN() == 0)
				continue;
			
			double logMean = logStats.getMean();
			double logSD = logStats.getStandardDeviation();
			double mean = stats.getMean();
			double sd = stats.getStandardDeviation();
			
			System.out.println("Log stats");
			System.out.println(logStats);
			System.out.println();
			System.out.println("Linear stats");
			System.out.println(stats);
			
//			NormalDistribution logNormal = new NormalDistribution(logMean, logSD);
//			double cov = logSD / logMean;
			double cov = sd / mean;
			List<String> line = new ArrayList<>();
			line.add((float)minMag+"");
			line.add((float)mean+"");
			line.add((float)sd+"");
			line.add((float)cov+"");
			csv.addLine(line);
			LogNormalDistribution logNormal = new LogNormalDistribution(Math.log(mean), cov);
			
			System.out.println("COV: "+cov);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			hist.setName(getCatalogName()+" Histogram");
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, getSecondaryColor()));
			
			EvenlyDiscretizedFunc logNormFit = HistogramFunction.getEncompassingHistogram(0d, 4d, 0.01);
			System.out.println("Cum prob at 4: "+logNormal.cumulativeProbability(4d));
			for (int i=0; i<logNormFit.size(); i++) {
				double x = logNormFit.getX(i);
				double x0 = x - 0.5*logNormFit.getDelta();
				double x1 = x + 0.5*logNormFit.getDelta();
				double prob = logNormal.probability(x0, x1) * (hist.getDelta()/logNormFit.getDelta());
				logNormFit.set(i, prob);
			}
			
			logNormFit.setName("Log-Normal Fit");
			funcs.add(logNormFit);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, getPrimaryColor()));
			
			BPT_DistCalc bpt = new BPT_DistCalc();
			bpt.setAll(mean, cov, logNormFit.getDelta(), logNormFit.size());
			EvenlyDiscretizedFunc bptFit = bpt.getPDF();
			bptFit.scale(bptFit.getDelta()/hist.getDelta());
			
			bptFit.setName("BPT Fit");
			funcs.add(bptFit);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, getComparableColor()));
			
//			funcs.add(getLine(getCatalogName()+" Mean="+yearDF.format(mean), mean, minY, mean, maxY));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, getPrimaryColor()));
			
			String magLabel = getCleanMagLabel(minMag);
			String myPrefix = prefix+"_m"+magLabel;
			
			String title = "M≥"+magLabel+" U3 Normalized Interevent Times";
			String xAxisLabel;
			if (aveType.isAveNTS()) {
				title += ", AveNTS";
				xAxisLabel = "Average Normalized Interevent Time";
			} else {
				title += ", AveTS";
				xAxisLabel = "Normalized Average Interevent Time";
			}
			PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, "Probability");
			plot.setLegendVisible(true);
			
			Range xRange = new Range(0, hist.getMaxX() + 0.5*hist.getDelta());
			double maxY = 0d;
			for (XY_DataSet func : funcs)
				maxY = Math.max(maxY, func.getMaxY());
			Range yRange = new Range(0, maxY*1.1);
			
			List<XYTextAnnotation> anns = new ArrayList<>();
			XYTextAnnotation ann = new XYTextAnnotation("COV: "+yearDF.format(cov),
					xRange.getUpperBound()*0.9, yRange.getUpperBound()*0.9);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 30));
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			anns.add(ann);
			plot.setPlotAnnotations(anns);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(plot, false, false, xRange, yRange);
			gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
			gp.saveAsTXT(new File(getOutputDir(), myPrefix+".txt").getAbsolutePath());
			gp.saveAsPNG(new File(getOutputDir(), myPrefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(getOutputDir(), myPrefix+".pdf").getAbsolutePath());
		}
		csv.writeToFile(new File(getOutputDir(), prefix+".csv"));
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
//		D = true;
		
		double skipYears = 5000;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		U3StyleNormalizedRuptureRecurrenceIntervalPlot plot =
				new U3StyleNormalizedRuptureRecurrenceIntervalPlot(
						catalog.getElements(), catalog.getSubSectMapper(), 6d, 6.5d, 7d, 7.5d);
		plot.initialize(catalog.getName(), outputDir, "u3_norm_ts");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
