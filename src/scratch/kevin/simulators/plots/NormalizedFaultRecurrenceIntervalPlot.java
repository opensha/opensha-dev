package scratch.kevin.simulators.plots;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
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
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.plots.SectionRecurrenceComparePlot.SectType;

public class NormalizedFaultRecurrenceIntervalPlot extends AbstractPlot {
	
	private Map<Integer, SimulatorElement> elemsMap;
	private SectType sectType;
	private double minFractForInclusion;
	private double[] minMags;
	private double overallMinMag;
	
	private RSQSimSubSectionMapper mapper;
	
	private Table<Integer, Double, List<Double>> idToTimesTable;
	
	public NormalizedFaultRecurrenceIntervalPlot(List<SimulatorElement> elems, double... minMags) {
		this(elems, SectType.ELEMENT, null, 0d, minMags);
	}
	
	public NormalizedFaultRecurrenceIntervalPlot(List<SimulatorElement> elems, SectType sectType, RSQSimSubSectionMapper mapper,
			double minFractForInclusion, double... minMags) {
		this.sectType = sectType;
		Preconditions.checkArgument(sectType == SectType.ELEMENT || mapper != null, "Must supply mapper if anything but element level selected");
		this.mapper = mapper;
		this.minFractForInclusion = minFractForInclusion;
		if (minMags == null || minMags.length == 0)
			minMags = new double[] { 0d };
		this.minMags = minMags;
		this.overallMinMag = StatUtils.min(minMags);

		elemsMap = new HashMap<>();
		for (SimulatorElement e : elems)
			elemsMap.put(e.getID(), e);
		
		idToTimesTable = HashBasedTable.create();
	}

	public SectType getSectType() {
		return sectType;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < overallMinMag)
			return;
		
		HashSet<Integer> ids = new HashSet<>();
		
		if (sectType == SectType.ELEMENT) {
			for (int id : e.getAllElementIDs())
				ids.add(id);
		} else {
			List<List<SubSectionMapping>> bundled =  mapper.getFilteredSubSectionMappings((RSQSimEvent)e, minFractForInclusion);
			if (minFractForInclusion >= 0 && bundled.isEmpty())
				bundled = mapper.getAllSubSectionMappings((RSQSimEvent)e);
			Preconditions.checkState(!bundled.isEmpty());
			for (List<SubSectionMapping> bundle : bundled) {
				Preconditions.checkState(!bundle.isEmpty());
				for (SubSectionMapping mapping : bundle) {
					FaultSectionPrefData sect = mapping.getSubSect();
					if (sectType == SectType.PARENT)
						ids.add(sect.getParentSectionId());
					else
						ids.add(sect.getSectionId());
				}
			}
		}
		
		double time = e.getTimeInYears();
		
		for (double minMag : minMags) {
			if (mag >= minMag) {
				for (Integer id : ids) {
					List<Double> times = idToTimesTable.get(id, minMag);
					if (times == null) {
						times = new ArrayList<>();
						idToTimesTable.put(id, minMag, times);
					}
					times.add(time);
				}
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
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 4d, 0.1);
			
			SummaryStatistics stats = new SummaryStatistics();
			SummaryStatistics logStats = new SummaryStatistics();
			
			for (int id : idToTimesTable.rowKeySet()) {
				List<Double> times = idToTimesTable.get(id, minMag);
				if (times == null || times.size() == 1)
					continue;
				double meanRI = durationYears/times.size();
				for (int i=1; i<times.size(); i++) {
					double deltaTime = times.get(i) - times.get(i-1);
					Preconditions.checkState(Doubles.isFinite(deltaTime) && deltaTime > 0);
					double normTime = deltaTime / meanRI;
					Preconditions.checkState(Doubles.isFinite(normTime) && normTime > 0);
					hist.add(hist.getClosestXIndex(normTime), 1d);
					stats.addValue(normTime);
					logStats.addValue(Math.log(normTime));
				}
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
			
			String title = "M≥"+magLabel+" "+sectType.getSimType()+" Normalized Interevent Times";
			PlotSpec plot = new PlotSpec(funcs, chars, title, "Normalized Interevent Time", "Probability");
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

}
