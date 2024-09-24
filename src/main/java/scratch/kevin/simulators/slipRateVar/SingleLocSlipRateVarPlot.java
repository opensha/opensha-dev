package scratch.kevin.simulators.slipRateVar;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SingleLocSlipRateVarPlot {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5892.instance();
		List<? extends FaultSection> subSects = catalog.getSubSects();
		
		int[] parentIDs = {
				FaultSectionUtils.findParentSectionID(subSects, "Garlock", "center")
		};
		
		double windowSize = 10000;
		
		CPT colors = GMT_CPT_Files.CATEGORICAL_TAB10_NOGRAY.instance();
//		int numTopToColor = colors.size();
		int numTopToColor = 5;
		double topWindowSize = 3000;
		double topWindowCalcDiscr = 10;
		
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		SlipAlongSectAlgorithm slipAlg = SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN;
		mapper.trackSlipOnSections();
		
		for (int parentID : parentIDs) {
			List<FaultSection> matchingSubSects = new ArrayList<>();
			for (FaultSection sect : subSects)
				if (sect.getParentSectionId() == parentID)
					matchingSubSects.add(sect);
			Preconditions.checkState(!matchingSubSects.isEmpty());
			FaultSection middle = matchingSubSects.get(matchingSubSects.size()/2);
			
			DiscretizedFunc cumulativeSlipFunc = new ArbitrarilyDiscretizedFunc();
			double cumulativeSlip = 0d;
			
			List<RSQSimEvent> events = catalog.loader().skipYears(20000).forSections(false, middle.getSectionId()).minMag(5d).load();
			System.out.println("Loaded "+events.size()+" events for "+middle.getSectionName());
			
			double firstEventTime = 0d;
			double lastEventTime = 0d;
			int numMatches = 0;
			for (RSQSimEvent event : events) {
				List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings(event);
				
				for (List<SubSectionMapping> subMappings : mappings) {
					for (SubSectionMapping mapping : subMappings) {
						if (mapping.getSubSect().getSectionId() == middle.getSectionId()) {
							double slip = mapping.getAverageSlip(slipAlg);
							cumulativeSlip += slip;
							double years = event.getTimeInYears();
							cumulativeSlipFunc.set(years, cumulativeSlip);
							if (firstEventTime == 0d)
								firstEventTime = years;
							lastEventTime = years;
							numMatches++;
						}
					}
				}
			}
			
			double duration = lastEventTime - firstEventTime;
			// add half a seismic cycle
			double ri = duration/numMatches;
			duration += 0.5*ri;
			
			double slipRate = cumulativeSlip/duration; // m/yr
			
			System.out.println("Slip rate: "+(slipRate*1e3)+" mm/yr");
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			List<VariationContainer> funcMaxVariations = new ArrayList<>();
			
			double maxFinalSlip = 0d;
			
			int prevIndex = 0;
			double prevLastSlip = 0d;
			
			MinMaxAveTracker fullWindowDiffTrack = new MinMaxAveTracker();
			MinMaxAveTracker fullWindowEventCountTrack = new MinMaxAveTracker();
			MinMaxAveTracker subWindowDiffTrack = new MinMaxAveTracker();
			MinMaxAveTracker subWindowEventCountTrack = new MinMaxAveTracker();
			for (double startTime=firstEventTime; startTime+windowSize<lastEventTime; startTime+=windowSize) {
				double endTime = startTime + windowSize;
				DefaultXY_DataSet subFunc = new DefaultXY_DataSet();
				double slipAtStart = prevLastSlip;
				subFunc.set(0d, 0d);
				double prevRelSlip = 0d;
				int eventCount = 0;
				for (int i=prevIndex; i<cumulativeSlipFunc.size(); i++) {
					double time = cumulativeSlipFunc.getX(i);
					if (time < startTime)
						continue;
					if (time > endTime)
						break;
					eventCount++;
					prevIndex = i;
					double slip = cumulativeSlipFunc.getY(i);
					prevLastSlip = slip;
					double relSlip = slip - slipAtStart;
					double relTime = time - startTime;
					subFunc.set(relTime, prevRelSlip);
					prevRelSlip = relSlip;
					Preconditions.checkState(relTime >= 0d);
					subFunc.set(relTime, relSlip);
				}
				fullWindowEventCountTrack.addValue(eventCount);
				double finalSlip = subFunc.getMaxY();
				subFunc.set(windowSize, finalSlip);
				maxFinalSlip = Math.max(maxFinalSlip, finalSlip);
				funcs.add(subFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
				
				// see how far that got away from the mean
				double maxRateDiff = 0d;
				double slipRateForMaxDiff = 0d;
				for (double subStart=0d; (float)subStart<=(float)(windowSize-topWindowSize); subStart+=topWindowCalcDiscr) {
					double subEnd = subStart + topWindowSize;
					double startSlip = 0d;
					double endSlip = 0d;
					double prevSlip = 0d;
					int subEventCount = 0;
					for (Point2D pt : subFunc) {
						if ((float)pt.getX() <= (float)subStart)
							startSlip = pt.getY();
						if ((float)pt.getX() <= (float)subEnd)
							endSlip = pt.getY();
						else
							break;
						if (pt.getX() > subStart && pt.getX() <= subEnd) {
							if (pt.getY() > prevSlip)
								// this was an event
								subEventCount++;
						}
						prevSlip = pt.getY();
					}
					subWindowEventCountTrack.addValue(subEventCount);
					double totSlip = endSlip - startSlip;
					Preconditions.checkState(totSlip > 0d);
					double subSlipRate = totSlip / topWindowSize;
					double diff = subSlipRate - slipRate;
					if (Math.abs(diff) > Math.abs(maxRateDiff)) {
						maxRateDiff = diff;
						slipRateForMaxDiff = subSlipRate;
					}
				}
				double windowSlipRate = finalSlip/windowSize;
				double windowSlipDiff = windowSlipRate - slipRate;
				fullWindowDiffTrack.addValue(Math.abs(windowSlipDiff));
				subWindowDiffTrack.addValue(Math.abs(maxRateDiff));
				funcMaxVariations.add(new VariationContainer(subFunc, windowSlipRate, windowSlipDiff, slipRateForMaxDiff, maxRateDiff));
			}
			
			funcs.get(0).setName("Individual "+(int)windowSize+"yr Windows");
			
			DecimalFormat slipDF = new DecimalFormat("0.0");
			DecimalFormat pDF = new DecimalFormat("0%");
			
			System.out.println("Average full ("+(int)windowSize+"yr) slip rate difference: "+pDF.format(fullWindowDiffTrack.getAverage()/slipRate));
			System.out.println("Average short ("+(int)topWindowSize+"yr) slip rate difference: "+pDF.format(subWindowDiffTrack.getAverage()/slipRate));
			System.out.println("Full window event counts: "+fullWindowEventCountTrack);
			System.out.println("Short window event counts: "+subWindowEventCountTrack);
			
			ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
			oneToOne.set(0d, 0d);
			oneToOne.set(windowSize, slipRate*windowSize);
			
			oneToOne.setName("Overall Slip Rate: "+slipDF.format(slipRate*1e3)+" mm/yr");
			
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
			
			Collections.sort(funcMaxVariations);
			if (funcMaxVariations.size() > numTopToColor)
				funcMaxVariations = funcMaxVariations.subList(0, numTopToColor);
			
			for (int i=0; i<funcMaxVariations.size(); i++) {
				Color color = colors.get(i% colors.size()).minColor;
				VariationContainer var = funcMaxVariations.get(i);
				XY_DataSet clone = new DefaultXY_DataSet(var.func.xValues(), var.func.yValues());
				String shortPercentStr = pDF.format(var.shortSlipRateDiff/slipRate);
				if (var.shortSlipRateDiff >= 0d)
					shortPercentStr = "+"+shortPercentStr;
				String windowPercentStr = pDF.format(var.windowSlipRateDiff/slipRate);
				if (var.windowSlipRateDiff >= 0d)
					windowPercentStr = "+"+windowPercentStr;
				clone.setName((int)windowSize+"yr: "+slipDF.format(var.windowSlipRate*1e3)+" mm/yr ("+windowPercentStr
						+"); "+(int)topWindowSize+"yr: "+slipDF.format(var.shortSlipRate*1e3)+" mm/yr ("+shortPercentStr+")");
				funcs.add(clone);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, middle.getParentSectionName(), "Windowed Time (years)", "Cumulative Slip (m)");
			spec.setLegendInset(RectangleAnchor.TOP_LEFT);
			
			Range xRange = new Range(0d, windowSize);
			Range yRange = new Range(0d, Math.ceil(maxFinalSlip/5d)*5d);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			String prefix = "slips_"+middle.getParentSectionName().replaceAll("\\W+", "_");
			while (prefix.endsWith("_"))
				prefix = prefix.substring(0, prefix.length()-1);
			PlotUtils.writePlots(new File("/tmp"), prefix, gp, 1000, 850, true, false, false);
		}
	}
	
	private static class VariationContainer implements Comparable<VariationContainer> {
		public final XY_DataSet func;
		public final double windowSlipRate;
		public final double windowSlipRateDiff;
		public final double shortSlipRate;
		public final double shortSlipRateDiff;
		private VariationContainer(XY_DataSet func, double windowSlipRate, double windowSlipRateDiff,
				double shortSlipRate, double shortSlipRateDiff) {
			super();
			this.func = func;
			this.windowSlipRate = windowSlipRate;
			this.windowSlipRateDiff = windowSlipRateDiff;
			this.shortSlipRate = shortSlipRate;
			this.shortSlipRateDiff = shortSlipRateDiff;
		}
		@Override
		public int compareTo(VariationContainer o) {
			return -Double.compare(Math.abs(shortSlipRateDiff), Math.abs(o.shortSlipRateDiff));
		}
	}

}
