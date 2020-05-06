package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

public class ExtremeRupturePageGen {
	
	private static enum Metric {
		NUM_SUB_SECTS("Subsection Count", "sub_sects", "Subsection[s]",
				"Total count of mapped UCERF3 subsections (e.g. SAF Mojave S Subsection 3), "
				+ "after application of minimum subsection area filter") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				return eventSections.size();
			}
		},
		NUM_PARENT_SECTS("Parent Section Count", "parent_sects", "Section[s]",
				"Total count of mapped UCERF3 subsections (e.g. SAF Mojave S), "
				+ "after application of minimum subsection area filter") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				HashSet<Integer> parentIDs = new HashSet<>();
				for (FaultSectionPrefData sect : eventSections)
					parentIDs.add(sect.getParentSectionId());
				return parentIDs.size();
			}
		},
//		RAW_LENTH_RATIO("Raw Length Ratio", "len_ratio", null, "Ratio of the total rupture length (raw rupture, not UCERF3 mapped)"
//				+ " to the idealized length, defined as the distance between the furthest 2 elements") {
//			@Override
//			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
//					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
//				double totLen = SimulatorUtils.estimateRuptureLength(event);
//				List<SimulatorElement> elems = event.getAllElements();
//				double idealLen = 0d;
//				for (int i=0; i<elems.size(); i++)
//					for (int j=i+1; j<elems.size(); j++)
//						idealLen = Double.max(idealLen, getElementDist2D(elementDistsCache, elems.get(i), elems.get(j)));
//				return totLen / idealLen;
//			}
//		},
		MAPPED_LENTH_RATIO("Mapped Length Ratio", "mapped_len_ratio", null, "Ratio of the total rupture length (UCERF3 mapped subsection rupture)"
				+ " to the idealized length, defined as the straight line distance between the furthest 2 subsections") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				double totLen = 0d;
				for (FaultSectionPrefData sect : eventSections)
					totLen += sect.getTraceLength();
				double idealLen = calcIdealMinLength(eventSections, subSectsDistCache);
				Preconditions.checkState(idealLen > 0);
				double ratio = totLen / idealLen;
//				if (ratio > 50)
//					System.out.println("ratio="+ratio+" = "+totLen+" / "+idealLen);
				return ratio;
			}
		},
		MAPPED_EXCESS_LENTH("Mapped Excess Length", "mapped_len_excess", "km", "Total rupture length (UCERF3 mapped subsection rupture)"
				+ " minus the idealized length, defined as the straight line distance between the furthest 2 subsections") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				double totLen = 0d;
				for (FaultSectionPrefData sect : eventSections)
					totLen += sect.getTraceLength();
				double idealLen = calcIdealMinLength(eventSections, subSectsDistCache);
				Preconditions.checkState(idealLen > 0);
				return totLen - idealLen;
			}
		},
		MOMENT_OFF_MAPPED("Moment Off Mapped Rupture", "moment_off_mapped", "N-m",
				"Moment of simulator elements not included in mapped UCERF3 subsection rupture") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				
				HashSet<Integer> mappedIDs = new HashSet<>();
				for (FaultSectionPrefData sect : eventSections)
					mappedIDs.add(sect.getSectionId());
				
				double momentOff = 0d;
				double[] allSlips = event.getAllElementSlips();
				ArrayList<SimulatorElement> elems = event.getAllElements();
				for (int i=0; i<elems.size(); i++) {
					SimulatorElement elem = elems.get(i);
					int sectIndex = elem.getSectionID() - subSectIDOffset;
					if (!mappedIDs.contains(sectIndex))
						momentOff += FaultMomentCalc.getMoment(elem.getArea(), allSlips[i]);
				}
				return momentOff;
			}
		},
		MAG("Magnitude", "mag", null, "Event Moment Magnitude") {
			@Override
			public double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				return event.getMagnitude();
			}
		},
		MOMENT_FAR_FROM_MAPPED("Moment >100km from Mapped", "mom_far_from_mapped", "N-m",
				"Moment that is at least 100km from the nearest mapped subsection (after application of subsection area threshold)", 1d) {
			private final double buffer = 100;
			private Map<FaultSectionPrefData, Region> sectRegions;
			@Override
			public synchronized double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
					int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache) {
				if (sectRegions == null)
					sectRegions = new HashMap<>();
				List<Region> myRegions = new ArrayList<>(eventSections.size());
				for (FaultSectionPrefData sect : eventSections) {
					Region region = sectRegions.get(sect);
					if (region == null) {
						region = new Region(sect.getFaultTrace(), buffer);
						myRegions.add(region);
					}
				}
				double momentOff = 0d;
				double[] allSlips = event.getAllElementSlips();
				ArrayList<SimulatorElement> elems = event.getAllElements();
				for (int i=0; i<elems.size(); i++) {
					SimulatorElement elem = elems.get(i);
					boolean inside = false;
					for (Region reg : myRegions) {
						if (reg.contains(elem.getCenterLocation())) {
							inside = true;
							break;
						}
					}
					if (!inside)
						momentOff += FaultMomentCalc.getMoment(elem.getArea(), allSlips[i]);
				}
				return momentOff;
			}
		};
		
		private String name;
		private String prefix;
		private String units;
		private String description;
		private double threshold;
		
		private Metric(String name, String prefix, String units, String description) {
			this(name, prefix, units, description, Double.NaN);
		}

		private Metric(String name, String prefix, String units, String description, double threshold) {
			this.name = name;
			this.prefix = prefix;
			this.units = units;
			this.description = description;
			this.threshold = threshold;
		}
		
		public abstract double calculate(RSQSimCatalog catalog, RSQSimEvent event, List<FaultSectionPrefData> eventSections,
				int subSectIDOffset, Map<IDPairing, Double> subSectsDistCache, Map<IDPairing, Double> elementDistsCache);
	}
	
	private static double calcIdealMinLength(List<FaultSectionPrefData> subSects, Map<IDPairing, Double> subSectsDistCache) {
		FaultSectionPrefData farS1 = null;
		FaultSectionPrefData farS2 = null;
		double maxDist = 0d;
		for (int i=0; i<subSects.size(); i++) {
			FaultSectionPrefData s1 = subSects.get(i);
			for (int j=i; j<subSects.size(); j++) {
				FaultSectionPrefData s2 = subSects.get(j);
				double dist = getSubSectDist(subSectsDistCache, s1, s2);
				if (dist >= maxDist) {
					maxDist = dist;
					farS1 = s1;
					farS2 = s2;
				}
			}
		}
		if (farS1 == farS2)
			return farS1.getTraceLength();
		maxDist = 0d;
		for (Location l1 : farS1.getFaultTrace())
			for (Location l2 : farS2.getFaultTrace())
				maxDist = Math.max(maxDist, LocationUtils.horzDistanceFast(l1, l2));
		return maxDist;
	}
	
	private static class EventScore {
		RSQSimEvent event;
		Map<Metric, Double> scores;
		
		public EventScore(RSQSimEvent event, Map<Metric, Double> scores) {
			this.event = event;
			this.scores = scores;
		}
	}
	
	private static class ScoreComparator implements Comparator<EventScore> {
		
		private Metric metric;

		public ScoreComparator(Metric metric) {
			this.metric = metric;
		}

		@Override
		public int compare(EventScore o1, EventScore o2) {
			Double v1 = o1.scores.get(metric);
			Double v2 = o2.scores.get(metric);
			return v2.compareTo(v1);
		}
		
	}
	
	private static double getSubSectDist(Map<IDPairing, Double> distsCache, FaultSectionPrefData s1, FaultSectionPrefData s2) {
		if (s1.getSectionId() == s2.getSectionId())
			return 0d;
		IDPairing pairing = new IDPairing(s1.getSectionId(), s2.getSectionId());
		Double cachedDist = distsCache.get(pairing);
		if (cachedDist != null)
			return cachedDist;
		double minDist = Double.POSITIVE_INFINITY;
		for (Location loc1 : s1.getFaultTrace()) {
			for (Location loc2 : s2.getFaultTrace()) {
				double dist = LocationUtils.horzDistance(loc1, loc2);
				if (dist < minDist)
					minDist = dist;
			}
		}
		distsCache.put(pairing, minDist);
		distsCache.put(pairing.getReversed(), minDist);
		return minDist;
	}
	
	private static double getElementDist2D(Map<IDPairing, Double> distsCache, SimulatorElement e1, SimulatorElement e2) {
		if (e1.getID() == e2.getID())
			return 0d;
		IDPairing pairing = new IDPairing(e1.getID(), e2.getID());
		Double cachedDist = distsCache.get(pairing);
		if (cachedDist != null)
			return cachedDist;
		double minDist = Double.POSITIVE_INFINITY;
		for (Location p1 : e1.getVertices())
			for (Location p2 : e2.getVertices())
				minDist = Math.min(minDist, LocationUtils.horzDistanceFast(p1, p2));
		distsCache.put(pairing, minDist);
		distsCache.put(pairing.getReversed(), minDist);
		return minDist;
	}
	
	private static boolean plotScoreHistogram(Metric metric, List<EventScore> eventScores, File resourcesDir, String prefix)
			throws IOException {
		double[] scores = new double[eventScores.size()];
		for (int i=0; i<eventScores.size(); i++)
			scores[i] = eventScores.get(i).scores.get(metric);
		double min = StatUtils.min(scores);
		double max = StatUtils.max(scores);
		if ((float)min == (float)max)
			return false;
		System.out.println("Plotting histogram for "+metric.name+" with min="+(float)min+" and max="+(float)max);
		
		HistogramFunction hist;
		if (min >= 0 && max <= 500) {
			double delta = max - min;
			if (delta < 1)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 0.1d);
			else if (delta < 2)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 0.2d);
			else if (delta < 5)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 0.5d);
			else if (delta < 15)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 1d);
			else if (delta < 25)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 2d);
			else if (delta < 50)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 5d);
			else if (delta < 105)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 10);
			else if (delta < 250)
				hist = HistogramFunction.getEncompassingHistogram(min, max, 20);
			else
				hist = HistogramFunction.getEncompassingHistogram(min, max, 50);
		} else {
			hist = new HistogramFunction(min, max, 11);
		}
		System.out.println("Hist bounds: ["+hist.getMinX()+" "+hist.getMaxX()+"], num="+hist.size()+", delta="+hist.getDelta());
		
		for (double score : scores)
			hist.add(score, 1d);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		String xAxisLabel = metric.name;
		if (metric.units != null)
			xAxisLabel += " ("+metric.units+")";
		PlotSpec spec = new PlotSpec(funcs, chars, metric.name, xAxisLabel, "Num Events");
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		gp.drawGraphPanel(spec, false, false, new Range(hist.getMinX()-0.5*hist.getDelta(), hist.getMaxX()+0.5*hist.getDelta()), null);
		gp.getChartPanel().setSize(600, 500);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		
		return true;
	}
	
	private static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File mainOutputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		
//		RSQSimCatalog catalog = Catalogs.JG_tuneBase1m.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance(baseDir);
		
		RSQSimUtils.populateFaultIDWithParentIDs(catalog.getElements(), catalog.getU3SubSects());
		
		double skipYears = 5000d;
		double minMag = 7d;
		
		Metric[] metrics = Metric.values();
		int numToPlot = 5;
		
		File catalogOutputDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File outputDir = new File(catalogOutputDir, "extreme_events");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<RSQSimEvent> events = catalog.loader().skipYears(skipYears).minMag(minMag).load();
		
		System.out.println("Found "+events.size()+" events with M>"+minMag);
		
		List<EventScore> eventScores = new ArrayList<>();
		
		Map<IDPairing, Double> sectDistsCache = catalog.getSubSectDistsCache();
		Map<IDPairing, Double> elemDistsCache = new HashMap<>();
		
		int subSectIDOffset = RSQSimUtils.getSubSectIndexOffset(catalog.getElements(), catalog.getU3SubSects());
		
		for (RSQSimEvent event : events) {
			Map<Metric, Double> scores = new HashMap<>();
			
			List<FaultSectionPrefData> subSects = catalog.getSubSectsForRupture(event);
			
			for (Metric metric : metrics)
				scores.put(metric, metric.calculate(catalog, event, subSects, subSectIDOffset, sectDistsCache, elemDistsCache));
			
			eventScores.add(new EventScore(event, scores));
		}
		
		List<String> lines = new ArrayList<>();
		
		// header
		lines.add("# "+catalog.getName()+" Extreme Events");
		lines.add("");
		lines.add("*Subsections participate in a rupture if at least "
				+(float)(catalog.getMinSubSectFractForInclusion()*100d)+" % of its area ruptures*");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (Metric metric : metrics) {
			lines.add("## "+metric.name);
			lines.add(topLink); lines.add("");
			lines.add("");
			lines.add(metric.description);
			lines.add("");
			if (!Double.isNaN(metric.threshold)) {
				int numAbove = 0;
				for (EventScore score : eventScores) {
					if (score.scores.get(metric) >= metric.threshold)
						numAbove++;
				}
				lines.add("");
				String line = numAbove+" events above threshold of "+optionalDigitDF.format(metric.threshold);
				if (metric.units != null)
					line += " ["+metric.units+"]";
				lines.add(line);
				if (numAbove == 0)
					continue;
			}
			String histPrefix = metric.prefix+"_hist";
			if (plotScoreHistogram(metric, eventScores, resourcesDir, histPrefix)) {
				lines.add("### "+metric.name+" Histogram");
//				lines.add(topLink); lines.add("");
				
				lines.add("![Histogram]("+resourcesDir.getName()+"/"+histPrefix+".png)");
			}
			lines.add("### "+metric.name+" Events");
			lines.add(topLink); lines.add("");
			TableBuilder table = MarkdownUtils.tableBuilder();
			table.addLine("Event ID", metric.name, "Plot");
			ScoreComparator comp = new ScoreComparator(metric);
			eventScores.sort(comp);
			for (int i=0; i<numToPlot && i<eventScores.size(); i++) {
				EventScore score = eventScores.get(i);
				double val = score.scores.get(metric);
				
				if (!Double.isNaN(metric.threshold) && val < metric.threshold)
					break;
				
				String plotPrefix = "event_"+score.event.getID();
				File plotFile = new File(resourcesDir, plotPrefix+".png");
				
				if (!plotFile.exists()) {
					RuptureSurface surf = catalog.getMappedSubSectRupture(score.event).getRuptureSurface();
//					RupturePlotGenerator.OTHER_SURF_COLOR = Color.RED;
					RupturePlotGenerator.OTHER_SURF_COLOR = new Color(139, 69, 19).darker();
					RupturePlotGenerator.OTHER_SURF_STROKE = PlotLineType.DOTTED;
					double[] elementSlips = score.event.getAllElementSlips();
					CPT slipCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, StatUtils.max(elementSlips));
//					for (CPTVal v : slipCPT) {
//						v.minColor = v.minColor.darker();
//						v.maxColor = v.maxColor.darker();
//					}
					RupturePlotGenerator.writeMapPlot(catalog.getElements(), score.event, null, resourcesDir, plotPrefix,
							null, null, surf, elementSlips, slipCPT, "Slip (m)");
					File pdfFile = new File(resourcesDir, plotPrefix+".pdf");
					if (pdfFile.exists())
						pdfFile.delete();
				}
				String valStr;
				if (val > 1e4 || val < 1e-3)
					valStr = (float)val+"";
				else
					valStr = optionalDigitDF.format(val);
				if (metric.units != null)
					valStr += " ("+metric.units+")";
				table.addLine("**"+score.event.getID()+"**", "**"+valStr+"**",
						"![Rupture Plot]("+resourcesDir.getName()+"/"+plotPrefix+".png)");
			}
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);

		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}

}
