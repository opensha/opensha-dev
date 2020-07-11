package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYDrawableAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public abstract class SourceSiteDetailPageGen {

	private SimulationRotDProvider<RSQSimEvent> simProv;
	private int[] parentSectIDs;
	private HashSet<Integer> parentIDsSet;
	private String sourceName;
	private RSQSimCatalog catalog;
	private List<RSQSimEvent> events;
	private List<Site> sites;
	
	private Location refFrame1;
	private Location refFrame2;
	private double refFrameDist;
	private double refFrameAz;
	
	private double minDeltaX = 10d;
	private double maxDeltaX = 120d;
	private double minDeltaY = 30d;
	private double maxDeltaY = 70d;

	public SourceSiteDetailPageGen(SimulationRotDProvider<RSQSimEvent> simProv, String sourceName,
			int[] parentSectIDs, RSQSimCatalog catalog, List<RSQSimEvent> events, List<Site> sites)
					throws IOException {
		this.simProv = simProv;
		this.sourceName = sourceName;
		this.parentSectIDs = parentSectIDs;
		parentIDsSet = new HashSet<>();
		for (int id : parentSectIDs)
			parentIDsSet.add(id);
		this.catalog = catalog;
		this.events = events;
		this.sites = sites;
		
		// figure out reference frame
		System.out.println("Calculating reference frame");
		List<FaultSection> subSects = new ArrayList<>();
		for (FaultSection sect : catalog.getSubSectMapper().getSubSections()) {
			if (parentIDsSet.contains(sect.getParentSectionId()))
				subSects.add(sect);
		}
		double maxDist = 0d;
		for (int i=0; i<subSects.size(); i++) {
			FaultTrace trace1 = subSects.get(i).getFaultTrace();
			Location p11 = trace1.first();
			Location p12 = trace1.last();
			double dist = LocationUtils.horzDistanceFast(p11, p12);
			if (dist > maxDist) {
				refFrame1 = p11;
				refFrame2 = p12;
				maxDist = dist;
			}
			for (int j=i; j<subSects.size(); j++) {
				FaultTrace trace2 = subSects.get(j).getFaultTrace();
				Location p21 = trace2.first();
				Location p22 = trace2.last();
				dist = LocationUtils.horzDistanceFast(p11, p21);
				if (dist > maxDist) {
					refFrame1 = p11;
					refFrame2 = p21;
					maxDist = dist;
				}
				dist = LocationUtils.horzDistanceFast(p11, p22);
				if (dist > maxDist) {
					refFrame1 = p11;
					refFrame2 = p22;
					maxDist = dist;
				}
				dist = LocationUtils.horzDistanceFast(p12, p21);
				if (dist > maxDist) {
					refFrame1 = p12;
					refFrame2 = p21;
					maxDist = dist;
				}
				dist = LocationUtils.horzDistanceFast(p12, p22);
				if (dist > maxDist) {
					refFrame1 = p12;
					refFrame2 = p22;
					maxDist = dist;
				}
			}
		}
		refFrameAz = LocationUtils.azimuth(refFrame1, refFrame2);
		if (refFrameAz < 0d || refFrameAz > 180d) {
			// flip them
			Location tmp = refFrame1;
			refFrame1 = refFrame2;
			refFrame2 = tmp;
			refFrameAz = LocationUtils.azimuth(refFrame1, refFrame2);
		}
		refFrameDist = maxDist;
		System.out.println("Reference frame: "+(float)maxDist+" km, "+
				(float)refFrameAz+" degrees");
		System.out.println("\t"+refFrame1);
		System.out.println("\t"+refFrame2);
	}
	
	public void generatePage(File outputDir, List<String> metadataLines, IMT[] imts) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<String> lines = new LinkedList<>();
		
		lines.add("# "+sourceName+" Distributions");
		
		if (metadataLines != null && !metadataLines.isEmpty()) {
			lines.addAll(metadataLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		List<Range> magRanges = new ArrayList<>();
//		magRanges.add(new Range(6, 6.5));
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
		
		// plot MNDs
		System.out.println("Plotting MNDs");
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		IncrementalMagFreqDist particMND = new IncrementalMagFreqDist(5.05, 41, 0.1d);
		IncrementalMagFreqDist nuclMND = new IncrementalMagFreqDist(5.05, 41, 0.1d);
		
		Map<RSQSimEvent, Location> eventHypos = new HashMap<>();
		
		for (RSQSimEvent event : events) {
			List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings(event);
			boolean match = false; // make sure that it fully maps, not just one element
			for (List<SubSectionMapping> parentMappings : mappings) {
				for (SubSectionMapping mapping : parentMappings) {
					if (parentIDsSet.contains(mapping.getSubSect().getParentSectionId())) {
						match = true;
						break;
					}
				}
			}
			if (!match)
				continue;
			int magIndex = particMND.getClosestXIndex(event.getMagnitude());
			particMND.add(magIndex, 1d);
			double earliestTime = Double.POSITIVE_INFINITY;
			SimulatorElement hypoElem = null;
			for (EventRecord rec : event) {
				List<SimulatorElement> patches = rec.getElements();
				double[] patchTimes = rec.getElementTimeFirstSlips();
				for (int i=0; i<patches.size(); i++) {
					if (patchTimes[i] < earliestTime) {
						earliestTime = patchTimes[i];
						hypoElem = patches.get(i);
					}
				}
			}
			eventHypos.put(event, hypoElem.getCenterLocation());
			if (parentIDsSet.contains(mapper.getMappedSection(hypoElem).getParentSectionId()))
				nuclMND.add(magIndex, 1d);
		}
		
		plotMNDs(resourcesDir, "mag_num_dist", particMND, nuclMND);
		
		lines.add("## Magnitude-Number Distribution");
		lines.add(topLink); lines.add("");
		lines.add("![MND](resources/mag_num_dist.png)");
		
		List<Range> xRanges = new ArrayList<>();
		List<Range> yRanges = new ArrayList<>();
		List<PlotSpec> magHypoSpecs = new ArrayList<>();
		List<List<RSQSimEvent>> magEventLists = new ArrayList<>();
		List<String> magEventPlotPrefixes = new ArrayList<>();
		List<Double> estHeightEaches = new ArrayList<>();
		List<HistogramFunction> hypoHists = new ArrayList<>();
		
		System.out.println("Plotting events for mag bins");
		
		double estBufferEach = 20d;
		double estStandardBuffer = 100d;
		int width = 1000;
		
		for (int m=0; m<magRanges.size(); m++) {
			List<RSQSimEvent> matches = new ArrayList<>();
			
			Range magRange = magRanges.get(m);

			HashSet<SimulatorElement> particElems = new HashSet<>();
			Map<SimulatorElement, Integer> particCountsMap = new HashMap<>();
			// count of times where this element is to the left of the hypocenter in the
			// current reference frame
			Map<SimulatorElement, Integer> particLeftCountsMap = new HashMap<>();
			
			for (RSQSimEvent event : events) {
				if (magRange.contains(event.getMagnitude())) {
					matches.add(event);
					Point2D hypoPt = refFramePt(eventHypos.get(event));
					ArrayList<SimulatorElement> elems = event.getAllElements();
					for (SimulatorElement elem : elems) {
						particElems.add(elem);
						Integer prevCount = particCountsMap.get(elem);
						if (prevCount == null)
							prevCount = 0;
						particCountsMap.put(elem, prevCount+1);
						Point2D myPt = refFramePt(elem.getCenterLocation());
						if (myPt.getX() < hypoPt.getX()) {
							// we're left of the hypocenter
							prevCount = particLeftCountsMap.get(elem);
							if (prevCount == null)
								prevCount = 0;
							particLeftCountsMap.put(elem, prevCount+1);
						}
					}
				}
			}
			magEventLists.add(matches);
			
			String magStr = magStrs.get(m);
			System.out.println(magStr+" has "+matches.size()+" events");
			
			if (matches.isEmpty()) {
				xRanges.add(null);
				yRanges.add(null);
				magHypoSpecs.add(null);
				magEventPlotPrefixes.add(null);
				estHeightEaches.add(null);
				hypoHists.add(null);
				continue;
			}
			
			Range xRange = new Range(-maxDeltaX, refFrameDist + maxDeltaX);
			Range yRange = new Range(-maxDeltaY, maxDeltaY);
			
			List<DefaultXY_DataSet> particElemXYs = getElemXYs(particElems, xRange, yRange);
			// now shrink ranges if we can
			MinMaxAveTracker xTrack = new MinMaxAveTracker();
			MinMaxAveTracker yTrack = new MinMaxAveTracker();
			for (DefaultXY_DataSet xy : particElemXYs) {
				for (Point2D pt : xy) {
					xTrack.addValue(pt.getX());
					yTrack.addValue(pt.getY());
				}
			}
			double minX = xRange.getLowerBound();
			if (xTrack.getMin() > minX)
				minX = Math.min(xTrack.getMin(), -minDeltaX);
			double maxX = xRange.getUpperBound();
			if (xTrack.getMax() < maxX)
				maxX = Math.max(xTrack.getMax(), refFrameDist + minDeltaX);
			xRange = new Range(minX, maxX);
			xRanges.add(xRange);
			double minY = yRange.getLowerBound();
			if (yTrack.getMin() > minY)
				minY = Math.min(yTrack.getMin(), -minDeltaY);
			double maxY = yRange.getUpperBound();
			if (yTrack.getMax() < maxY)
				maxY = Math.max(yTrack.getMax(), minDeltaY);
			yRange = new Range(minY, maxY);
			yRanges.add(yRange);
			
			List<DefaultXY_DataSet> allElems = getElemXYs(catalog.getElements(), xRange, yRange);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			for (int i=0; i<allElems.size(); i++) {
				DefaultXY_DataSet xy = allElems.get(i);
				if (i == 0)
					xy.setName("All Elems");
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(220, 220, 220)));
			}
			
			for (int i=0; i<particElemXYs.size(); i++) {
				DefaultXY_DataSet xy = particElemXYs.get(i);
				if (i == 0)
					xy.setName("Participating Elems");
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(150, 150, 150)));
			}
			
			List<SimulatorElement> particFaultElems = new ArrayList<>();
			for (SimulatorElement elem : particElems)
				if (parentIDsSet.contains(mapper.getMappedSection(elem).getParentSectionId()))
					particFaultElems.add(elem);
			List<DefaultXY_DataSet> particFaultXYs = getElemXYs(particFaultElems, xRange, yRange);
			
			for (int i=0; i<particFaultXYs.size(); i++) {
				DefaultXY_DataSet xy = particFaultXYs.get(i);
				if (i == 0)
					xy.setName("Source Elems");
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
			}
			

			DefaultXY_DataSet hypos = new DefaultXY_DataSet();
			for (RSQSimEvent event : matches) {
				Point2D pt = refFramePt(eventHypos.get(event));
				if (pt.getX() < xRange.getLowerBound())
					pt.setLocation(xRange.getLowerBound(), pt.getY());
				else if (pt.getX() > xRange.getUpperBound())
					pt.setLocation(xRange.getUpperBound(), pt.getY());
				hypos.set(pt);
			}
			
			hypos.setName("Hypocenters");
			funcs.add(hypos);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, Color.BLUE.darker()));
			
			System.out.println("Plotting...");
			
			PlotSpec spec = new PlotSpec(funcs, chars, magStr+" Event Locations",
					"DAS (km)", "FN Dist (km)");
			spec.setLegendVisible(true);
			
			magHypoSpecs.add(spec);
			
			// hypocenter count
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			HistogramFunction hypoFunc = new HistogramFunction(
					xRange.getLowerBound(), xRange.getUpperBound(), 50);
			for (Point2D pt : hypos)
				hypoFunc.add(hypoFunc.getClosestXIndex(pt.getX()), 1d);
			hypoHists.add(hypoFunc);
			
			Range hypoYRange = new Range(0, hypoFunc.getMaxY()*1.05);
			
			funcs.add(hypoFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE.darker()));
			
			PlotSpec hypoSpec = new PlotSpec(funcs, chars, magStr+" Event Locations",
					"DAS (km)", "Hypocenter Count");
			
			// plot participation rates
			int maxCount = 0;
			for (Integer count : particCountsMap.values())
				maxCount = Integer.max(maxCount, count);
			
			CPT countCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
			countCPT = countCPT.rescale(0, Math.ceil(Math.log10(maxCount)));
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			for (int i=0; i<allElems.size(); i++) {
				DefaultXY_DataSet xy = allElems.get(i);
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(220, 220, 220)));
			}
			
			for (SimulatorElement elem : ComparablePairing.getSortedData(particCountsMap)) {
				for (DefaultXY_DataSet xy : getElemXYs(Lists.newArrayList(elem), xRange, yRange)) {
					Color c = countCPT.getColor((float)Math.log10(particCountsMap.get(elem).floatValue()));
					funcs.add(xy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
				}
			}
			
			PlotSpec particSpec = new PlotSpec(funcs, chars, spec.getTitle(),
					spec.getXAxisLabel(), spec.getYAxisLabel());
			
			// plot directivities rates
			CPT directivityCPT = new CPT(-1, 1d,
					new Color(0, 0, 140), new Color(0, 60, 200 ), new Color(0, 120, 255),
					new Color(100, 200, 255), Color.WHITE, new Color(255, 200, 100),
					new Color(255, 120, 0), new Color(200, 60, 0), new Color(140, 0, 0));
			
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			for (int i=0; i<allElems.size(); i++) {
				DefaultXY_DataSet xy = allElems.get(i);
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(220, 220, 220)));
			}
			
			for (SimulatorElement elem : ComparablePairing.getSortedData(particCountsMap)) {
				for (DefaultXY_DataSet xy : getElemXYs(Lists.newArrayList(elem), xRange, yRange)) {
					int particCount = particCountsMap.get(elem);
					Integer leftOfHypoCount = particLeftCountsMap.get(elem);
					if (leftOfHypoCount == null)
						leftOfHypoCount = 0;
					double directivityRightFract = 1d - (double)leftOfHypoCount/(double)particCount;
					// -1 means all leftward directivity
					// 0 means balanced
					// 1 means all rightward directivity
					double directivityVal = directivityRightFract*2d - 1;
					Color c = directivityCPT.getColor((float)directivityVal);
					funcs.add(xy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
				}
			}
			
			PlotSpec directivitySpec = new PlotSpec(funcs, chars, spec.getTitle(),
					spec.getXAxisLabel(), spec.getYAxisLabel());
//			new XY
//			directivitySpec.addPlotAnnotation(annotation);
			
			List<PlotSpec> specs = Lists.newArrayList(particSpec, directivitySpec, hypoSpec, spec);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(24);
			gp.setLegendFontSize(24);
			gp.setBackgroundColor(Color.WHITE);
			
			PlotPreferences prefs = gp.getPlotPrefs();

			double aspectRatio = xRange.getLength()/yRange.getLength();
			double estHeightEach = (width-100d)/aspectRatio;
			estHeightEaches.add(estHeightEach);
			double height = specs.size()*(estHeightEach + estBufferEach) + estStandardBuffer;
			
			double scaleAnnX = xRange.getCentralValue();
			double scaleWidth = 500d;
			double scaleHeight = 75d;
			double scaleCoordsHeight = yRange.getLength()*scaleHeight/estHeightEach;
			double scaleAnnY = yRange.getLowerBound() + yRange.getLength()*0.95d - 0.5*scaleCoordsHeight;
			
//			System.out.println("\taspect ratio: "+aspectRatio);
//			System.out.println("\test height each: "+estHeightEach);
			
			PaintScaleLegend particBar = XYZGraphPanel.getLegendForCPT(countCPT, "Log₁₀ Participation Count",
					prefs.getAxisLabelFontSize(), prefs.getTickLabelFontSize(), 1d, RectangleEdge.TOP);
			particBar.setBackgroundPaint(new Color(255, 255, 255, 127));
			XYDrawableAnnotation particBarAnn = new XYDrawableAnnotation(
					scaleAnnX, scaleAnnY, scaleWidth, scaleHeight, particBar);
			particSpec.addPlotAnnotation(particBarAnn);
			
			PaintScaleLegend dirBar = XYZGraphPanel.getLegendForCPT(directivityCPT,
					"⇦Leftward   Ave-Prop-Dir  Rightward⇨",
					prefs.getAxisLabelFontSize(), prefs.getTickLabelFontSize(), 1d, RectangleEdge.TOP);
			dirBar.setBackgroundPaint(new Color(255, 255, 255, 127));
			XYDrawableAnnotation dirBarAnn = new XYDrawableAnnotation(
					scaleAnnX, scaleAnnY, scaleWidth, scaleHeight, dirBar);
			directivitySpec.addPlotAnnotation(dirBarAnn);
//			XYTextAnnotation leftDirAnn = new XYTextAnnotation("Leftward", x, y);
//			XYTextAnnotation rightDirAnn = new XYTextAnnotation("Rightward", x, y);
//			directivitySpec.addPlotAnnotation(annotation);

			gp.drawGraphPanel(specs, false, false, Lists.newArrayList(xRange),
					Lists.newArrayList(yRange, yRange, hypoYRange, yRange));
			
			String prefix = magPrefixes.get(m)+"_events";
			magEventPlotPrefixes.add(prefix);
			
			File file = new File(resourcesDir, prefix);
			
			gp.getChartPanel().setSize(width, (int)height);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
		}
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		for (int m=0; m<magStrs.size(); m++)
			if (magPrefixes.get(m) != null)
				table.addColumn("**"+magStrs.get(m)+"**");
		table.finalizeLine();
		table.initNewLine();
		for (int m=0; m<magStrs.size(); m++)
			if (magPrefixes.get(m) != null)
				table.addColumn("![Map plot](resources/"+magEventPlotPrefixes.get(m)+".png)");
		table.finalizeLine();

		lines.add("");
		lines.add("## Event Section & Hypocenter Map");
		lines.add(topLink); lines.add("");
		lines.addAll(table.wrap(2, 0).build());
		
		for (Site site : sites) {
			lines.add("");
			lines.add("## "+site.getName());
			lines.add(topLink); lines.add("");
			
			Point2D sitePt = refFramePt(site.getLocation());

			List<String> myMagTitles = new ArrayList<>();
			List<String> sitePrefixes = new ArrayList<>();
			
			magLoop:
			for (int m=0; m<magRanges.size(); m++) {
				String magStr = magStrs.get(m);
				
				Range xRange = xRanges.get(m);
				if (xRange == null)
					continue;
				Range yRange = yRanges.get(m);
				double aspectRatio = xRange.getLength()/yRange.getLength();
				HistogramFunction hypoHist = hypoHists.get(m);

				String title = site.getName()+", "+magStr;
				System.out.println("Doing "+title);
				
				PlotSpec mapSpec = magHypoSpecs.get(m);
				
				List<PlotSpec> specs = new ArrayList<>();
				List<Range> myYRanges = new ArrayList<>();
				
				for (IMT imt : imts) {
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					List<RSQSimEvent> magEventList = magEventLists.get(m);
					int spanAlpha = Integer.max(30, (int)(255d*20d/magEventList.size()));
					spanAlpha = Integer.min(spanAlpha, 127);
					
					DefaultXY_DataSet xy = new DefaultXY_DataSet();
					for (RSQSimEvent event : magEventList) {
						double val;
						try {
							val = Math.log10(simProv.getValue(site, event, imt, 0));
						} catch (Exception e) {
//							e.printStackTrace();
							continue; // possible for rupture to not be applicable to site
						}
						double hypoX = refFramePt(eventHypos.get(event)).getX();
						if (hypoX < xRange.getLowerBound())
							hypoX = xRange.getLowerBound();
						else if (hypoX > xRange.getUpperBound())
							hypoX = xRange.getUpperBound();
						xy.set(hypoX, val);

						DefaultXY_DataSet spanX = new DefaultXY_DataSet();
						double minX = Double.POSITIVE_INFINITY;
						double maxX = Double.NEGATIVE_INFINITY;
						for (SimulatorElement elem : event.getAllElements()) {
							for (Point2D pt : elemPointsCache.get(elem)) {
								minX = Math.min(minX, pt.getX());
								maxX = Math.max(maxX, pt.getX());
							}
						}
						minX = Math.max(minX, xRange.getLowerBound());
						maxX = Math.min(maxX, xRange.getUpperBound());
						spanX.set(minX, val);
						spanX.set(maxX, val);
						funcs.add(spanX);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f,
								new Color(0, 0, 0, spanAlpha)));
					}
					
					if (xy.size() == 0) {
						System.out.println("No site ruptures");
						continue magLoop;
					}
					
					funcs.add(xy);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, Color.BLUE.darker()));

					double minY = 0.5*Math.floor(2d*xy.getMinY());
					double maxY = 0.5*Math.ceil(2d*xy.getMaxY());
					myYRanges.add(new Range(minY, maxY));
					
					if (xRange.contains(sitePt.getX())) {
						DefaultXY_DataSet siteLine = new DefaultXY_DataSet();
						siteLine.set(sitePt.getX(), minY);
						siteLine.set(sitePt.getX(), maxY);
						
						funcs.add(0, siteLine);
						chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GREEN.darker()));
					}
					
					String yAxisLabel = "Log10 "+imt.getDisplayName();
					PlotSpec spec = new PlotSpec(funcs, chars, title, mapSpec.getXAxisLabel(), yAxisLabel);
					
					int histNum = (int)(hypoHist.size()/aspectRatio + 0.5);
					HistogramFunction imtHist = new HistogramFunction(minY, maxY, histNum+1);
					imtHist = new HistogramFunction(minY + 0.5*imtHist.getDelta(), histNum, imtHist.getDelta());
					for (Point2D pt : xy)
						imtHist.add(imtHist.getClosestXIndex(pt.getY()), 1);
					double scalarY = 0.2*yRange.getLength()/imtHist.getMaxY();
					List<XYAnnotation> anns = new ArrayList<>();
					anns.addAll(buildHistAnns(imtHist, true, xRange, yRange, scalarY, new Color(255, 0, 0, 80)));
					spec.setPlotAnnotations(anns);
					
					specs.add(spec);
				}
				List<XYAnnotation> anns = new ArrayList<>();
				if (xRange.contains(sitePt.getX()) && yRange.contains(sitePt.getY())) {
					// add site annotation
					XYTextAnnotation nameAnn = new XYTextAnnotation(
							" "+site.getName(), sitePt.getX(), sitePt.getY());
					nameAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
					nameAnn.setTextAnchor(TextAnchor.TOP_LEFT);
					anns.add(nameAnn);
					
					Shape shape = new Ellipse2D.Double(sitePt.getX(), sitePt.getY(), 3f, 3f);
					XYShapeAnnotation shapeAnn = new XYShapeAnnotation(shape,
							null, null, Color.GREEN.darker());
					anns.add(shapeAnn);
				}
				// add hist annotation
				double scalarX = 0.2*yRange.getLength()/hypoHist.getMaxY();
				anns.addAll(buildHistAnns(hypoHist, false, xRange, yRange, scalarX, new Color(0, 0, 255, 80)));
				mapSpec.setPlotAnnotations(anns);
				specs.add(mapSpec);
				myYRanges.add(yRange);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(24);
				gp.setPlotLabelFontSize(24);
				gp.setLegendFontSize(24);
				gp.setBackgroundColor(Color.WHITE);
				
				gp.drawGraphPanel(specs, false, false, Lists.newArrayList(xRange), myYRanges);
				
				String prefix = site.getName()+"_"+magPrefixes.get(m)+"_gms";
				
				File file = new File(resourcesDir, prefix);
				
				double estHeightEach = estHeightEaches.get(m);
				double height = specs.size()*(estHeightEach + estBufferEach) + estStandardBuffer;
				gp.getChartPanel().setSize(width, (int)height);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				
				sitePrefixes.add(prefix);
				myMagTitles.add(title);
			}
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (String title : myMagTitles)
				table.addColumn("**"+title+"**");
			table.finalizeLine();
			table.initNewLine();
			for (String prefix : sitePrefixes)
				table.addColumn("![Ground Motions](resources/"+prefix+".png)");
			table.finalizeLine();
			
			lines.addAll(table.wrap(2, 0).build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static void plotMNDs(File resourcesDir, String prefix,
			IncrementalMagFreqDist particMND, IncrementalMagFreqDist nuclMND) throws IOException {
		EvenlyDiscretizedFunc particCumFunc = particMND.getCumRateDistWithOffset();
		particCumFunc.setName(null);
		EvenlyDiscretizedFunc nuclCumFunc = nuclMND.getCumRateDistWithOffset();
		nuclCumFunc.setName(null);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		particMND.setName("Participation");
		funcs.add(particMND);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));

		nuclMND.setName("Nucleation");
		funcs.add(nuclMND);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, new Color(100, 100, 255)));
		
		funcs.add(particCumFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		funcs.add(nuclCumFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = 0d;
		for (int i=0; i<particMND.size(); i++) {
			if (particMND.getY(i) > 0d) {
				double x = particMND.getX(i);
				if (x < minMag)
					minMag = 0.5*Math.floor(2d*x + particMND.getDelta());
				maxMag = 0.5*Math.ceil(2d*x + particMND.getDelta());
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Magnitude Number Distribution",
				"Magnitude", "Count");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		double maxY = Math.pow(10, Math.ceil(Math.log10(particCumFunc.getMaxY())));
		
		Range xRange = new Range(minMag, maxMag);
		Range yRange = new Range(0.95d, maxY);
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		File file = new File(resourcesDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private Point2D refFramePt(Location loc) {
		double dist = LocationUtils.horzDistanceFast(refFrame1, loc);
		double distToLine = LocationUtils.distanceToLineFast(refFrame1, refFrame2, loc);
		double y = -distToLine;
		double x = Math.sqrt(dist*dist - distToLine*distToLine);
		double az = LocationUtils.azimuth(refFrame1, loc);
		if (RSQSimBBP_Config.absAngleDiff(az, refFrameAz) > 90d)
			x = -x;
		return new Point2D.Double(x, y);
	}
	
	private List<XYAnnotation> buildHistAnns(HistogramFunction hist, boolean vertical,
			Range xAxis, Range yAxis, double scalar, Color color) {
		List<XYAnnotation> anns = new ArrayList<>();
		
		double baselineValue = vertical ? xAxis.getLowerBound() : yAxis.getLowerBound();
		
		double delta = hist.getDelta();
		double halfDelta = 0.5*delta;
		
		for (int i=0; i<hist.size(); i++) {
			double middle = hist.getX(i);
			double start = middle - halfDelta;
			double end = middle + halfDelta;
			
			double rawHeight = hist.getY(i);
			if (rawHeight == 0d)
				continue;
			double top = baselineValue + scalar*rawHeight;
			
			double x0, y0, x1, y1;
			if (vertical) {
				x0 = baselineValue;
				x1 = top;
				y0 = start;
				y1 = end;
			} else {
				y0 = baselineValue;
				y1 = top;
				x0 = start;
				x1 = end;
			}
			
			anns.add(new XYBoxAnnotation(x0, y0, x1, y1, null, null, color));
		}
		
		return anns;
	}
	
	private Map<SimulatorElement, List<Point2D>> elemPointsCache = new HashMap<>();
	
	private List<DefaultXY_DataSet> getElemXYs(Collection<SimulatorElement> elems, Range xRange,
			Range yRange) {
		List<DefaultXY_DataSet> ret = new ArrayList<>();
		
		for (SimulatorElement elem : elems) {
			List<Point2D> elemPoints = elemPointsCache.get(elem);
			if (elemPoints == null) {
				elemPoints = new ArrayList<>();
				for (Location loc : elem.getVertices())
					elemPoints.add(refFramePt(loc));
				elemPointsCache.put(elem, elemPoints);
			}
			boolean inside = false;
			for (Point2D pt : elemPoints) {
				if (xRange.contains(pt.getX()) && yRange.contains(pt.getY())) {
					inside = true;
					break;
				}
			}
			if (inside) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				for (Point2D pt : elemPoints)
					xy.set(pt);
				xy.set(elemPoints.get(0));
				ret.add(xy);
			}
		}
		return ret;
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

}
