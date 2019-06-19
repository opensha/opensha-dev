package scratch.kevin.simulators.plots;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.DAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectDAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class ElasticReboundTriggeringPlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private double minMag;
	private double[] maxTimes;
	private double globalMaxTime;

	private HistogramFunction[] hists;
	private HistogramFunction[] reRuptureHists;
	
	private List<RuptureMappings> prevRups;
	
	private RuptureMappings[] exampleHypoReRup;
	private RuptureMappings[] exampleNewHypo;

	public ElasticReboundTriggeringPlot(RSQSimSubSectionMapper mapper, double minMag, double[] maxTimes) {
		this.mapper = mapper;
		this.minMag = minMag;
		this.maxTimes = maxTimes;

		hists = new HistogramFunction[maxTimes.length];
		reRuptureHists = new HistogramFunction[maxTimes.length];
		for (int i=0; i<hists.length; i++) {
			hists[i] = HistogramFunction.getEncompassingHistogram(0d, 3d, 0.1);
			reRuptureHists[i] = HistogramFunction.getEncompassingHistogram(0d, 3d, 0.1);
		}
		
		mapper.trackSlipOnSections();
		
		globalMaxTime = StatUtils.max(maxTimes);
		prevRups = new LinkedList<>();
	}
	
	private class RuptureMappings {
		private final List<SubSectionMapping> allMappings;
		private final Map<Integer, SubSectionMapping> idToMappings;
		private final HashSet<SimulatorElement> allElems;
		private final double timeYears;
		private final double length;
		private final double magnitude;
		private final SimulatorElement hypocenter;
		
		public RuptureMappings(List<SubSectionMapping> allMapings, HashSet<SimulatorElement> allElems,
				double magnitude, SimulatorElement hypocenter, double timeYears) {
			this.allMappings = allMapings;
			this.allElems = allElems;
			this.magnitude = magnitude;
			this.timeYears = timeYears;
			this.hypocenter = hypocenter;
			
			idToMappings = new HashMap<>();
			double length = 0d;
			for (SubSectionMapping mapping : allMapings) {
				idToMappings.put(mapping.getSubSect().getSectionId(), mapping);
				length += mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN);
			}
			this.length = length;
		}
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < minMag)
			return;
		
		double time = e.getTimeInYears();
		double cullTime = time - globalMaxTime;
		while (!prevRups.isEmpty() && prevRups.get(0).timeYears < cullTime)
			prevRups.remove(0);
		
		RSQSimEvent event = (RSQSimEvent)e;
		List<List<SubSectionMapping>> bundled = mapper.getFilteredSubSectionMappings(event);
		if (bundled.isEmpty())
			return;
		List<SubSectionMapping> allMappings = new ArrayList<>();
		for (List<SubSectionMapping> bundle : bundled)
			allMappings.addAll(bundle);
		SimulatorElement hypoEl = null;
		double minTime = Double.POSITIVE_INFINITY;
		HashSet<SimulatorElement> allElems = new HashSet<>();
		for (EventRecord rec : event) {
			RSQSimEventRecord rRec = (RSQSimEventRecord)rec;
			double[] times = rRec.getElementTimeFirstSlips();
			List<SimulatorElement> elems = rRec.getElements();
			allElems.addAll(elems);
			for (int i=0; i<times.length; i++) {
				if (times[i] < minTime) {
					minTime = times[i];
					hypoEl = elems.get(i);
				}
			}
		}
		FaultSectionPrefData hypoSect = mapper.getMappedSection(hypoEl);
		int hypoSectIndex = -1;
		boolean allVert = true;
		for (int i=0; i<allMappings.size(); i++) {
			FaultSectionPrefData subSect = allMappings.get(i).getSubSect();
			if (subSect == hypoSect)
				hypoSectIndex = i;
			allVert = allVert && subSect.getAveDip() == 90d;
		}
		double hypoSectDAS = Double.NaN;
		if (hypoSectIndex >= 0) {
			// hypocenter is on a mapped subsection
			SubSectionMapping mapping = allMappings.get(hypoSectIndex);
			DAS_Record das = mapper.getElemSubSectDAS(hypoEl);
			if (mapping.isReversed())
				das = das.getReversed(mapping.getSubSect().getFaultTrace().getTraceLength());
			hypoSectDAS = das.midDAS;
		}
		RuptureMappings rupMapping = new RuptureMappings(allMappings, allElems, mag, hypoEl, time);
		
		if (hypoSectIndex >= 0) {
			rupLoop:
				for (RuptureMappings prevMapping : prevRups) {
					// see if we encompass it
					boolean haveAll = true;
					for (Integer id : prevMapping.idToMappings.keySet())
						haveAll = haveAll && rupMapping.idToMappings.containsKey(id);
					if (haveAll && rupMapping.idToMappings.size() > prevMapping.idToMappings.size()) {
						// we encompass it
						double curDAS = 0d;
						double prevStartDAS = -1;
						double prevEndDAS = -1;
						double hypoDAS = -1;
						for (int i=0; i<allMappings.size(); i++) {
							SubSectionMapping mapping = allMappings.get(i);
							int sectID = mapping.getSubSect().getSectionId();
							if (i == hypoSectIndex)
								hypoDAS = curDAS + hypoSectDAS;
							double sectLen = mapping.getSubSect().getFaultTrace().getTraceLength();
							double myStartDAS;
							double myEndDAS;
							if (i == 0) {
								DAS_Record myDAS = mapping.getDASforSlip(SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN);
								if (mapping.isReversed())
									myDAS = myDAS.getReversed(sectLen);
								myStartDAS = myDAS.startDAS;
								myEndDAS = sectLen;
							} else  if (i == allMappings.size()-1) {
								DAS_Record myDAS = mapping.getDASforSlip(SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN);
								if (mapping.isReversed())
									myDAS = myDAS.getReversed(sectLen);
								myStartDAS = 0d;
								myEndDAS = myDAS.endDAS;
							} else {
								myStartDAS = 0d;
								myEndDAS = sectLen;
							}
							SubSectionMapping prevSectMapping = prevMapping.idToMappings.get(sectID);
							if (prevSectMapping != null) {
								DAS_Record myDAS = prevSectMapping.getDASforSlip(SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN);
								if (mapping.isReversed())
									myDAS = myDAS.getReversed(sectLen);
								if (prevStartDAS < 0) {
									prevStartDAS = curDAS + myDAS.startDAS;
									if (i == 0 && (float)prevStartDAS <= (float)myStartDAS)
										// doesn't actually encompass
										continue rupLoop;
								}
								prevEndDAS = curDAS + myDAS.endDAS;
								if (i == allMappings.size()-1 && prevEndDAS >= (curDAS + myEndDAS))
									// doesn't actually encompass
									continue rupLoop;
							}
							if (i == 0 && myStartDAS > 0) {
								// trim it
								myEndDAS -= myStartDAS;
								hypoDAS -= myStartDAS; // doesn't matter if NaN
								prevStartDAS -= myStartDAS; // doesn't matter if not yet set and negative
								prevEndDAS -= myStartDAS; // doesn't matter if not yet set and negative
								myStartDAS = 0;
							}
							curDAS += myEndDAS - myStartDAS;
						}
						double timeDiff = time - prevMapping.timeYears;
//						System.out.println("Found one "+timeDiff+" years later!");
						Preconditions.checkState(prevEndDAS > prevStartDAS,
								"prevStart=%s >= prevEnd=%s", prevStartDAS, prevEndDAS);
						if (prevStartDAS == 0d)
							continue;
						Preconditions.checkState(prevStartDAS > 0 && prevStartDAS < curDAS,
								"Bad prevStart=%s, curDAS=%s", prevStartDAS, curDAS);
						Preconditions.checkState(prevEndDAS > 0 && prevEndDAS < curDAS,
								"Bad prevEnd=%s, curDAS=%s", prevEndDAS, curDAS);
						Preconditions.checkState(hypoDAS >= 0 && hypoDAS <= curDAS,
								"Bad hypoDAS=%s, curDAS=%s", hypoDAS, curDAS);
						double hypoRelative;
						if (hypoDAS < prevStartDAS)
							hypoRelative = hypoDAS / prevStartDAS;
						else if (hypoDAS < prevEndDAS)
							hypoRelative = 1d + (hypoDAS-prevStartDAS)/(prevEndDAS-prevStartDAS);
						else
							hypoRelative = 2d + (hypoDAS - prevEndDAS)/(curDAS-prevEndDAS);
						double hypoRelative2 = (3d-hypoRelative);
						int hypoIndex1 = hists[0].getClosestXIndex(hypoRelative);
						int hypoIndex2 = hists[0].getClosestXIndex(hypoRelative2);
						// now see if our hypocenter element ruptured in that event
						boolean hypoPrevRuptured = prevMapping.allElems.contains(hypoEl);
						for (int i=0; i<maxTimes.length; i++) {
							if (timeDiff <= maxTimes[i]) {
								hists[i].add(hypoIndex1, 1d);
								hists[i].add(hypoIndex2, 1d);
								if (hypoPrevRuptured) {
									reRuptureHists[i].add(hypoIndex1, 1d);
									reRuptureHists[i].add(hypoIndex2, 1d);
								}
							}
						}
						
						double fractLenOverlap = prevMapping.length/rupMapping.length;
						if (allVert && fractLenOverlap > 0.2 && fractLenOverlap < 0.7
								&& rupMapping.allMappings.size() - prevMapping.allMappings.size() >= 1
								&& timeDiff <= StatUtils.min(maxTimes)) {
							// example candidate
							if (hypoPrevRuptured) {
								if (exampleHypoReRup == null)
									exampleHypoReRup = new RuptureMappings[] { prevMapping, rupMapping };
							} else {
								if (exampleNewHypo == null)
									exampleNewHypo = new RuptureMappings[] { prevMapping, rupMapping };
							}
						}
					}
				}
		}
		
		prevRups.add(rupMapping);
	}

	@Override
	public void finalizePlot() throws IOException {
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();
		
		for (int i=0; i<maxTimes.length; i++) {
			HistogramFunction hist = hists[i];
			HistogramFunction reRupHist = reRuptureHists[i];
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double sum = hist.calcSumOfY_Vals();
			if (sum > 0) {
				hist.scale(1d/sum);
				reRupHist.scale(1d/sum);
			}
			double histMax = hist.getMaxY();
			
			Range xRange = new Range(0d, 3d);
			Range yRange = new Range(0d, histMax > 0 ? 1.4*histMax : 1);
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			funcs.add(reRupHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
			
			double lineBottom = 0.01*yRange.getUpperBound();
			double insideTop = 0.89*yRange.getUpperBound();
			double outsideTop = 0.99*yRange.getUpperBound();
			
			XY_DataSet outsideLine = new DefaultXY_DataSet();
			double outsideLeft = 0.007*xRange.getUpperBound();
			double outsideRight = 0.993*xRange.getUpperBound();
			outsideLine.set(outsideLeft, lineBottom);
			outsideLine.set(outsideLeft, outsideTop);
			outsideLine.set(outsideRight, outsideTop);
			outsideLine.set(outsideRight, lineBottom);
			
			Color prevColor = Color.BLUE.darker();
			Color curColor = Color.RED.darker();
			
			funcs.add(outsideLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, curColor));
			
			XY_DataSet insideLine = new DefaultXY_DataSet();
			insideLine.set(1d, lineBottom);
			insideLine.set(1d, insideTop);
			insideLine.set(2d, insideTop);
			insideLine.set(2d, lineBottom);
			
			funcs.add(insideLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, prevColor));
			
			List<XYAnnotation> anns = new ArrayList<>();
			double outsideY = 0.98*outsideTop;
			double spacingY = outsideTop - insideTop;
			XYTextAnnotation curAnn = new XYTextAnnotation("Encompassing Rupture", 1.5, outsideY);
			curAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			curAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			curAnn.setPaint(curColor);
			anns.add(curAnn);
			double insideY = outsideY - spacingY;
			XYTextAnnotation prevAnn = new XYTextAnnotation("Previous Rupture", 1.5, insideY);
			prevAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			prevAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			prevAnn.setPaint(prevColor);
			anns.add(prevAnn);
			
			XYTextAnnotation occurAnn = new XYTextAnnotation((int)(sum/2d)+" Occurrences ",
					outsideRight, outsideY);
			occurAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 16));
			occurAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			occurAnn.setPaint(Color.BLACK);
			anns.add(occurAnn);
			
			double fractWithin = 0d;
			for (Point2D pt : hist)
				if (pt.getX() > 1 && pt.getX() < 2)
					fractWithin += pt.getY();
			XYTextAnnotation fractAnn = new XYTextAnnotation(optionalDigitDF.format(fractWithin*100d)+"% Within ",
					outsideRight, outsideY - 0.5*spacingY);
			fractAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 16));
			fractAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			fractAnn.setPaint(Color.GRAY);
			anns.add(fractAnn);
			
			double fractRerup = 0d;
			for (Point2D pt : reRupHist)
				if (pt.getX() > 1 && pt.getX() < 2)
					fractRerup += pt.getY();
			XYTextAnnotation fractRerupAnn = new XYTextAnnotation(optionalDigitDF.format(fractRerup*100d)+"% Reruptured ",
					outsideRight, outsideY - spacingY);
			fractRerupAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 16));
			fractRerupAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			anns.add(fractRerupAnn);
			
			String title = "Encompassing Rupture Triggering, "+optionalDigitDF.format(maxTimes[i])+" yr";
			if (maxTimes[i] != 1d)
				title += "s";
			if (minMag > 0)
				title += ", M≥"+optionalDigitDF.format(minMag);
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Normalized Distance Along Strike", "Hypocenter Density");
			spec.setPlotAnnotations(anns);
			
			String myPrefix = prefix+"_"+optionalDigitDF.format(maxTimes[i])+"yr";
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(800, 500);
			gp.getPlot().getDomainAxis().setTickLabelsVisible(false);
			gp.saveAsTXT(new File(outputDir, myPrefix+".txt").getAbsolutePath());
			gp.saveAsPNG(new File(outputDir, myPrefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, myPrefix+".pdf").getAbsolutePath());
		}

		if (exampleHypoReRup != null)
			plotExample(outputDir, prefix+"_example_re_rup", exampleHypoReRup[0], exampleHypoReRup[1]);
		if (exampleNewHypo != null)
			plotExample(outputDir, prefix+"_example_new_hypo", exampleNewHypo[0], exampleNewHypo[1]);
	}
	
	private void plotExample(File outputDir, String prefix, RuptureMappings prevRup, RuptureMappings encompRup) throws IOException {
		List<XYAnnotation> anns = new ArrayList<>();
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		HashSet<SimulatorElement> prevElems = new HashSet<>(prevRup.allElems);
		HashSet<SimulatorElement> encompElems = new HashSet<>(encompRup.allElems);
		
		double curDAS = 0d;
		
		Stroke thickElemStroke = new BasicStroke(1.5f);
		Stroke hypoElemStroke = new BasicStroke(10f);
		
		List<XYAnnotation> hypoAnns = new ArrayList<>();
		
		double maxDepth = 0d;
		
		for (SubSectionMapping mapping : encompRup.allMappings) {
			FaultSectionPrefData subSect = mapping.getSubSect();
			
			double sectLen = subSect.getTraceLength();
			
			for (SimulatorElement elem : mapper.getElementsForSection(subSect)) {
				boolean prev = prevElems.contains(elem);
				boolean encomp = encompElems.contains(elem);
				
				SubSectDAS_Record elemDAS = mapper.getElemSubSectDAS(elem);
				DefaultXY_DataSet elemXY = new DefaultXY_DataSet();
				Vertex[] verts = elem.getVertices();
				for (int i=0; i<verts.length; i++) {
					double vertDAS = elemDAS.vertDASs[i];
					if (mapping.isReversed())
						vertDAS = sectLen - vertDAS;
					double das = curDAS + vertDAS;
					elemXY.set(das, verts[i].getDepth());
				}
				
				if (prev || encomp) {
					boolean prevHypo = prevRup.hypocenter == elem;
					boolean encompHypo = encompRup.hypocenter == elem;
					
					double[] polyElems = new double[verts.length*2];
					int ind = 0;
					for (Point2D pt : elemXY) {
						polyElems[ind++] = pt.getX();
						polyElems[ind++] = pt.getY();
						maxDepth = Math.max(maxDepth, pt.getY());
					}
					
					Color color;
					if (prev && encomp)
						color = new Color(255, 0, 255);
					else if (prev)
						color = Color.BLUE;
					else
						color = Color.RED;
					Stroke stroke = prevHypo || encompHypo ? hypoElemStroke : thickElemStroke;
					XYPolygonAnnotation ann = new XYPolygonAnnotation(polyElems, stroke, Color.BLACK, color);
					anns.add(ann);
					if (stroke == hypoElemStroke)
						hypoAnns.add(ann);
				} else {
					funcs.add(elemXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
				}
			}
			
			curDAS += sectLen;
		}
		anns.addAll(hypoAnns);
		
		String title = "Encompassing Rupture Example: M"+optionalDigitDF.format(prevRup.magnitude)
			+" followed by M"+optionalDigitDF.format(encompRup.magnitude)
			+" ("+getTimeShortLabel(encompRup.timeYears-prevRup.timeYears)+" later)";
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Distance Along Strike (km)", "Depth (km)");
		spec.setPlotAnnotations(anns);
		
		Range xRange = new Range(0d, curDAS);
		Range yRange = new Range(0d, maxDepth+1);
		
		double xPerY = xRange.getLength()/yRange.getLength();
		int ySize = 800;
		int xSize = (int)Math.round(ySize*xPerY);
		
		HeadlessGraphPanel gp = buildGraphPanel();
		gp.setyAxisInverted(true);
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(xSize, ySize);
		gp.getPlot().setDomainGridlinesVisible(false);
		gp.getPlot().setRangeGridlinesVisible(false);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		double skipYears = 5000;
		
		double minMag = 7d;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		double[] maxTimes = { 1d, 10d, 100d };
		
		ElasticReboundTriggeringPlot plot = new ElasticReboundTriggeringPlot(catalog.getSubSectMapper(), minMag, maxTimes);
		plot.initialize(catalog.getName(), outputDir, "trigger_within_prev");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
