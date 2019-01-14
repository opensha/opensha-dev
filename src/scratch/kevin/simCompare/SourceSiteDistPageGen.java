package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

public class SourceSiteDistPageGen<E> {
	
	private SimulationRotDProvider<E> simProv;
	private List<Site> sites;
	
	private static double simLogHistDelta = 0.1;
	private static double gmpeLogHistDelta = 0.01;
	private static double gmpeTruncLevel = 4;

	public SourceSiteDistPageGen(SimulationRotDProvider<E> simProv, List<Site> sites) {
		this.simProv = simProv;
		this.sites = sites;
	}
	
	public void generatePage(Table<AttenRelRef, String, List<RuptureComparison<E>>> sourceComps,
			File outputDir, List<String> headerLines, double[] periods, boolean hypoFilter)
					throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<String> lines = new LinkedList<>();
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		List<String> sourceNames = new ArrayList<>(sourceComps.columnKeySet());
		Collections.sort(sourceNames);
		
		List<AttenRelRef> gmpeRefs = new ArrayList<>(sourceComps.rowKeySet());
		
		List<Range> magRanges = new ArrayList<>();
		magRanges.add(new Range(6, 6.5));
		magRanges.add(new Range(6.5, 7));
		magRanges.add(new Range(7, 7.5));
		magRanges.add(new Range(7.5, 8));
		magRanges.add(new Range(8, 8.5));
		List<String> magStrs = new ArrayList<>();
		List<String> magPrefixes = new ArrayList<>();
		for (Range magRange : magRanges) {
			String magPrefix = "m"+optionalDigitDF.format(magRange.getLowerBound())
				+"_"+optionalDigitDF.format(magRange.getUpperBound());
			String magStr = "M"+optionalDigitDF.format(magRange.getLowerBound())
				+"-"+optionalDigitDF.format(magRange.getUpperBound());
			magStrs.add(magStr);
			magPrefixes.add(magPrefix);
		}
		
		for (String sourceName : sourceNames) {
			Map<AttenRelRef, List<RuptureComparison<E>>> comps = sourceComps.column(sourceName);
			
			lines.add("## "+sourceName);
			lines.add(topLink); lines.add("");
			
			String sourcePrefix = sourceName.replaceAll("\\W+", "_");
			
			for (Site site : sites) {
				Table<Range, Double, File> plotsTable = HashBasedTable.create();
				
				for (int m=0; m<magRanges.size(); m++) {
					Range magRange = magRanges.get(m);
					Map<AttenRelRef, List<RuptureComparison<E>>> magComps = filterByMag(comps, magRange);
					if (magComps.get(gmpeRefs.get(0)).isEmpty())
						continue;
					
					String magStr = magStrs.get(m);
					String magPrefix = magPrefixes.get(m);
					
					int numSims = 0;
					int numRups = 0;
					for (RuptureComparison<E> comp : magComps.get(gmpeRefs.get(0))) {
						if (!comp.isComputed(site, periods[0]))
							continue;
						E rup = comp.getRupture();
						numRups++;
						numSims += simProv.getNumSimulations(site, rup);
					}
					
					System.out.println("Calculating for "+sourceName+", "+site.getName()+", "+magStr);
					
					HistogramFunction[][] simHists = calcSimHist(simProv, site, periods, magComps.get(gmpeRefs.get(0)), hypoFilter);
					if (simHists == null)
						// nothing for this site/mag range
						continue;
					Map<AttenRelRef, HistogramFunction[]> gmpeHists = new HashMap<>();
					for (AttenRelRef gmpeRef : gmpeRefs)
						gmpeHists.put(gmpeRef, calcGMPEHist(site, periods, magComps.get(gmpeRef)));
					
					for (int p=0; p<periods.length; p++) {
						String prefix = sourcePrefix+"_"+site.getName()+"_"+magPrefix+"_"+optionalDigitDF.format(periods[p])+"s";
						File file = new File(resourcesDir, prefix+".png");
						String title = sourceName+", "+site.getName()+", "+magStr;
						Map<AttenRelRef, HistogramFunction> gmpePeriodHists = new HashMap<>();
						for (AttenRelRef gmpeRef : gmpeHists.keySet())
							gmpePeriodHists.put(gmpeRef, gmpeHists.get(gmpeRef)[p]);
						plotSourceSiteHist(resourcesDir, prefix, simHists[p], gmpePeriodHists, title, periods[p], numSims, numRups);
						plotsTable.put(magRange, periods[p], file);
					}
				}
				
				if (plotsTable.isEmpty())
					continue;
				
				lines.add("### "+sourceName+", "+site.getName());
				lines.add(topLink); lines.add("");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				
				table.initNewLine().addColumn("Mag Range");
				for (double period : periods)
					table.addColumn("**"+optionalDigitDF.format(period)+"s**");
				table.finalizeLine();
				
				for (int m=0; m<magRanges.size(); m++) {
					Range magRange = magRanges.get(m);
					if (!plotsTable.containsRow(magRange))
						continue;
					
					table.initNewLine().addColumn("**"+magStrs.get(m)+"**");
					for (double period : periods)
						table.addColumn("![Plot]("+resourcesDir.getName()+"/"+plotsTable.get(magRange, period).getName()+")");
					table.finalizeLine();
				}
				
				lines.addAll(table.build());
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static <E> Map<AttenRelRef, List<RuptureComparison<E>>> filterByMag(
			Map<AttenRelRef, List<RuptureComparison<E>>> sourceComps, Range magRange) {
		Map<AttenRelRef, List<RuptureComparison<E>>> map = new HashMap<>();
		for (AttenRelRef gmpeRef : sourceComps.keySet()) {
			List<RuptureComparison<E>> ret = new ArrayList<>(sourceComps.get(gmpeRef));
			ret.removeIf(rup -> !magRange.contains(rup.getMagnitude()));
			map.put(gmpeRef, ret);
		}
		return map;
	}
	
//	private static final Color[] gmpeColors = { Color.BLUE.darker(), Color.RED.darker(),
//			Color.GREEN.darker(), Color.ORANGE.darker() };
	private static final Color[] gmpeColors = { Color.RED.brighter(), Color.YELLOW.brighter(),
			Color.CYAN.brighter(), Color.MAGENTA.brighter() };
	
	private static final Comparator<AttenRelRef> gmpeRefComparator = new Comparator<AttenRelRef>() {

		@Override
		public int compare(AttenRelRef arg0, AttenRelRef arg1) {
			return arg0.getShortName().compareTo(arg1.getShortName());
		}
	};
	
	private void plotSourceSiteHist(File outputDir, String prefix, HistogramFunction[] simHists,
			Map<AttenRelRef, HistogramFunction> gmpeHists, String title, double period, int numSims, int numRups)
					throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		List<AttenRelRef> gmpes = new ArrayList<>(gmpeHists.keySet());
		gmpes.sort(gmpeRefComparator);
		
		HistogramFunction simHist = simHists[0];
		
		simHist.setName(simProv.getName());
		funcs.add(simHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		if (simHists.length == 3) {
			HistogramFunction nearHist = simHists[1];
			HistogramFunction farHist = simHists[2];
			funcs.add(nearHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, new Color(0, 160, 0, 140)));
			funcs.add(farHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, new Color(160, 0, 0, 140)));
		}
		
		for (int i=0; i<gmpes.size(); i++) {
			AttenRelRef gmpeRef = gmpes.get(i);
			HistogramFunction gmpeHist = gmpeHists.get(gmpeRef);
			gmpeHist.setName(gmpeRef.getShortName());
			funcs.add(gmpeHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, gmpeColors[i % gmpeColors.length]));
		}
		
		double maxY = 0;
		double minX = Double.POSITIVE_INFINITY;
		double maxX = Double.NEGATIVE_INFINITY;
		for (DiscretizedFunc func : funcs) {
			maxY = Math.max(maxY, func.getMaxY());
			minX = Math.min(minX, func.getMinX());
			maxX = Math.max(maxX, func.getMaxX());
		}
		
		maxY *= 1.1;
		minX = Math.floor(minX);
		maxX = Math.ceil(maxX);
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(0, maxY);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Log10 "+optionalDigitDF.format(period)+"s RotD50", "Count");
		spec.setLegendVisible(true);
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		String countText;
		if (numSims == numRups)
			countText = numSims+" Ruptures";
		else
			countText = numSims+" Simulations of "+numRups+" Ruptures";
		XYTextAnnotation countAnn = new XYTextAnnotation("  "+countText, minX, maxY*0.95);
		countAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
		countAnn.setTextAnchor(TextAnchor.TOP_LEFT);
		anns.add(countAnn);
		spec.setPlotAnnotations(anns);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(600, 450);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}
	
	private static class HypoDistResult implements Comparable<HypoDistResult> {
		final Location hypo;
		final double[] vals;
		
		boolean sortLat = true;
		
		public HypoDistResult(Location hypo, double[] vals) {
			super();
			this.hypo = hypo;
			this.vals = vals;
		}
		@Override
		public int compareTo(HypoDistResult o) {
			return sortLat ? Double.compare(hypo.getLatitude(), o.hypo.getLatitude())
					: Double.compare(hypo.getLongitude(), o.hypo.getLongitude());
		}
	}
	
	private static <E> HistogramFunction[][] calcSimHist(SimulationRotDProvider<E> simProv, Site site, double[] periods,
			List<? extends RuptureComparison<E>> sourceComps, boolean hypoFilter) throws IOException {
		List<List<Double>> vals = new ArrayList<>();
		SummaryStatistics[] stats = new SummaryStatistics[periods.length];
		for (int p=0; p<periods.length; p++) {
			vals.add(new ArrayList<>());
			stats[p] = new SummaryStatistics();
		}
		
		List<HypoDistResult> hypoDistList = null;
		if (hypoFilter)
			hypoDistList = new ArrayList<>();
		for (RuptureComparison<E> comp : sourceComps) {
			if (!comp.isComputed(site, periods[0]))
				continue;
			E rup = comp.getRupture();
			int num = simProv.getNumSimulations(site, rup);
			for (int i=0; i<num; i++) {
				DiscretizedFunc func = simProv.getRotD50(site, rup, i);
				double[] myVals = new double[periods.length];
				for (int p=0; p<periods.length; p++)
					myVals[p] = Math.log10(func.getInterpolatedY(periods[p]));
				for (int p=0; p<periods.length; p++) {
					vals.get(p).add(myVals[p]);
					stats[p].addValue(myVals[p]);
				}
				if (hypoFilter) {
					Location hypo = simProv.getHypocenter(rup, i);
					Preconditions.checkNotNull(hypo);
					hypoDistList.add(new HypoDistResult(hypo, myVals));
				}
			}
		}
		
		if (vals.isEmpty())
			return null;
		
		hypoFilter = hypoFilter && hypoDistList.size() > 2;
		
		HistogramFunction[][] ret = hypoFilter ? new HistogramFunction[periods.length][3]
				: new HistogramFunction[periods.length][1];
		
		List<double[]> nearHypoVals = null;
		List<double[]> farHypoVals = null;
		if (hypoFilter) {
			double maxLat = Double.NEGATIVE_INFINITY;
			double minLat = Double.POSITIVE_INFINITY;
			double maxLon = Double.NEGATIVE_INFINITY;
			double minLon = Double.POSITIVE_INFINITY;
			
			for (HypoDistResult result : hypoDistList) {
				double lat = result.hypo.getLatitude();
				double lon = result.hypo.getLongitude();
				maxLat = Double.max(maxLat, lat);
				minLat = Double.min(minLat, lat);
				maxLon = Double.max(maxLon, lon);
				minLon = Double.min(minLon, lon);
			}
			double latSpan = maxLat - minLat;
			double lonSpan = maxLon - minLon;
			boolean sortLat = latSpan > lonSpan;
			for (HypoDistResult result : hypoDistList)
				result.sortLat = sortLat;
			
			Collections.sort(hypoDistList);
			int maxNearIndex = hypoDistList.size()/3;
			int minFarIndex = 2*hypoDistList.size()/3;
			nearHypoVals = new ArrayList<>();
			for (int i=0; i<=maxNearIndex; i++)
				nearHypoVals.add(hypoDistList.get(i).vals);
			farHypoVals = new ArrayList<>();
			for (int i=minFarIndex; i<hypoDistList.size(); i++)
				farHypoVals.add(hypoDistList.get(i).vals);
		}
		
		for (int p=0; p<periods.length; p++) {
			double max = stats[p].getMax();
			double min = stats[p].getMin();
			ret[p][0] = HistogramFunction.getEncompassingHistogram(min, max, simLogHistDelta);
			for (double val : vals.get(p))
				ret[p][0].add(ret[p][0].getClosestXIndex(val), 1d);
			if (hypoFilter) {
				ret[p][1] = HistogramFunction.getEncompassingHistogram(min, max, simLogHistDelta);
				for (double[] val : nearHypoVals)
					ret[p][1].add(ret[p][1].getClosestXIndex(val[p]), 1d);
				ret[p][2] = HistogramFunction.getEncompassingHistogram(min, max, simLogHistDelta);
				for (double[] val : farHypoVals)
					ret[p][2].add(ret[p][2].getClosestXIndex(val[p]), 1d);
				if (hypoDistList.get(0).sortLat) {
					ret[p][1].setName("South Hypos");
					ret[p][2].setName("North Hypos");
				} else {
					ret[p][1].setName("West Hypos");
					ret[p][2].setName("East Hypos");
				}
			}
		}
		
		return ret;
	}
	
	private HistogramFunction[] calcGMPEHist(Site site, double[] periods, List<? extends RuptureComparison<E>> sourceComps)
			throws IOException {
		List<List<Double>> logMeans = new ArrayList<>();
		List<List<Double>> stdDevs = new ArrayList<>();
		List<Integer> simCounts = new ArrayList<>();
		SummaryStatistics[] stats = new SummaryStatistics[periods.length];
		for (int p=0; p<periods.length; p++) {
			logMeans.add(new ArrayList<>());
			stdDevs.add(new ArrayList<>());
			stats[p] = new SummaryStatistics();
		}
		
		for (RuptureComparison<E> comp : sourceComps) {
			if (!comp.isComputed(site, periods[0]))
				continue;
			simCounts.add(simProv.getNumSimulations(site, comp.getRupture()));
			for (int p=0; p<periods.length; p++) {
				double logMean = comp.getLogMean(site, periods[p]);
				double stdDev = comp.getStdDev(site, periods[p]);
				stats[p].addValue(logMean - gmpeTruncLevel*stdDev);
				stats[p].addValue(logMean + gmpeTruncLevel*stdDev);
				logMeans.get(p).add(logMean);
				stdDevs.get(p).add(stdDev);
			}
		}
		
		HistogramFunction[] ret = new HistogramFunction[periods.length];
		
		NormalDistribution stdNorm = new NormalDistribution(0d, 1d);
		
		for (int p=0; p<periods.length; p++) {
			double max = Math.log10(Math.exp(stats[p].getMax()));
			double min = Math.log10(Math.exp(stats[p].getMin()));
			ret[p] = HistogramFunction.getEncompassingHistogram(min, max, gmpeLogHistDelta);
			List<Double> myMeans = logMeans.get(p);
			List<Double> myStdDevs = stdDevs.get(p);
			for (int i=0; i<myMeans.size(); i++) {
				double logMean = myMeans.get(i);
				double stdDev = myStdDevs.get(i);
				double scalar = simCounts.get(i)*simLogHistDelta/gmpeLogHistDelta;
				for (int h=0; h<ret[p].size(); h++) {
					double binCenter = ret[p].getX(h);
					double start = binCenter - 0.5*gmpeLogHistDelta;
					double end = binCenter + 0.5*gmpeLogHistDelta;
					// now from log10 to ln
					start = Math.log(Math.pow(10, start));
					end = Math.log(Math.pow(10, end));
					double startNorm = (start - logMean)/stdDev;
					double endNorm = (end - logMean)/stdDev;
					ret[p].add(h, stdNorm.probability(startNorm, endNorm)*scalar);
				}
			}
		}
		
		return ret;
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

}
