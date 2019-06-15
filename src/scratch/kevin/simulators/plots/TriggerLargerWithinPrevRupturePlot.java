package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.Font;
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
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
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
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.DAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectDAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class TriggerLargerWithinPrevRupturePlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private double minMag;
	private double[] maxTimes;
	private double globalMaxTime;
	
	private HistogramFunction[] hists;
	
	private List<RuptureMappings> prevRups;

	public TriggerLargerWithinPrevRupturePlot(RSQSimSubSectionMapper mapper, double minMag, double[] maxTimes) {
		this.mapper = mapper;
		this.minMag = minMag;
		this.maxTimes = maxTimes;
		
		hists = new HistogramFunction[maxTimes.length];
		for (int i=0; i<hists.length; i++)
			hists[i] = HistogramFunction.getEncompassingHistogram(0d, 3d, 0.1);
		
		mapper.trackSlipOnSections();
		
		globalMaxTime = StatUtils.max(maxTimes);
		prevRups = new LinkedList<>();
	}
	
	private class RuptureMappings {
		private final List<SubSectionMapping> allMapings;
		private final Map<Integer, SubSectionMapping> idToMappings;
		private final double timeYears;
		
		public RuptureMappings(List<SubSectionMapping> allMapings, double timeYears) {
			this.allMapings = allMapings;
			this.timeYears = timeYears;
			
			idToMappings = new HashMap<>();
			for (SubSectionMapping mapping : allMapings)
				idToMappings.put(mapping.getSubSect().getSectionId(), mapping);
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
		List<SubSectionMapping> allMappings = new ArrayList<>();
		for (List<SubSectionMapping> bundle : bundled)
			allMappings.addAll(bundle);
		SimulatorElement hypoEl = null;
		double minTime = Double.POSITIVE_INFINITY;
		for (EventRecord rec : event) {
			RSQSimEventRecord rRec = (RSQSimEventRecord)rec;
			double[] times = rRec.getElementTimeFirstSlips();
			List<SimulatorElement> elems = rRec.getElements();
			for (int i=0; i<times.length; i++) {
				if (times[i] < minTime) {
					minTime = times[i];
					hypoEl = elems.get(i);
				}
			}
		}
		FaultSectionPrefData hypoSect = mapper.getMappedSection(hypoEl);
		int hypoSectIndex = -1;
		for (int i=0; i<allMappings.size(); i++) {
			if (allMappings.get(i).getSubSect() == hypoSect) {
				hypoSectIndex = i;
				break;
			}
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
		RuptureMappings rupMapping = new RuptureMappings(allMappings, time);
		
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
							if (i == 0 && prevStartDAS <= myStartDAS)
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
//				System.out.println("Found one "+timeDiff+" years later!");
				Preconditions.checkState(prevEndDAS > prevStartDAS);
				Preconditions.checkState(prevStartDAS > 0 && prevStartDAS < curDAS);
				Preconditions.checkState(prevEndDAS > 0 && prevEndDAS < curDAS);
				Preconditions.checkState(hypoDAS >= 0 && hypoDAS <= curDAS);
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
				for (int i=0; i<maxTimes.length; i++) {
					if (timeDiff <= maxTimes[i]) {
						hists[i].add(hypoIndex1, 1d);
						hists[i].add(hypoIndex2, 1d);
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
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double sum = hist.calcSumOfY_Vals();
			if (sum > 0)
				hist.scale(1d/sum);
			double histMax = hist.getMaxY();
			
			Range xRange = new Range(0d, 3d);
			Range yRange = new Range(0d, histMax > 0 ? 1.4*histMax : 1);
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
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
			outsideLine.set(outsideLeft, lineBottom);
			
			Color prevColor = Color.BLUE.darker();
			Color curColor = Color.RED.darker();
			
			funcs.add(outsideLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, curColor));
			
			XY_DataSet insideLine = new DefaultXY_DataSet();
			insideLine.set(1d, lineBottom);
			insideLine.set(1d, insideTop);
			insideLine.set(2d, insideTop);
			insideLine.set(2d, lineBottom);
			insideLine.set(1d, lineBottom);
			
			funcs.add(insideLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, prevColor));
			
			List<XYAnnotation> anns = new ArrayList<>();
			double outsideY = 0.98*outsideTop;
			XYTextAnnotation curAnn = new XYTextAnnotation("Encompassing Rupture", 1.5, outsideY);
			curAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			curAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			curAnn.setPaint(curColor);
			anns.add(curAnn);
			double insideY = 0.98*insideTop;
			XYTextAnnotation prevAnn = new XYTextAnnotation("Previous Rupture", 1.5, 0.98*insideY);
			prevAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			prevAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			prevAnn.setPaint(prevColor);
			anns.add(prevAnn);
			double fractWithin = 0d;
			for (Point2D pt : hist)
				if (pt.getX() > 1 && pt.getX() < 2)
					fractWithin += pt.getY();
			XYTextAnnotation fractAnn = new XYTextAnnotation(optionalDigitDF.format(fractWithin*100d)+"% Within",
					1.5, insideY-2d*(outsideY - insideY));
			fractAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			fractAnn.setTextAnchor(TextAnchor.TOP_CENTER);
			anns.add(fractAnn);
			
			
			String title = "Encompassing Rupture Triggering ("+optionalDigitDF.format(maxTimes[i])+" yr";
			if (maxTimes[i] != 1d)
				title += "s";
			title += ")";
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
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		double skipYears = 5000;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		double[] maxTimes = { 1d, 10d, 100d };
		
		TriggerLargerWithinPrevRupturePlot plot = new TriggerLargerWithinPrevRupturePlot(catalog.getSubSectMapper(), 6.5, maxTimes);
		plot.initialize(catalog.getName(), outputDir, "trigger_within_prev");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
