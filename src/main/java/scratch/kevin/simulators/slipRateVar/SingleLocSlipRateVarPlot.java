package scratch.kevin.simulators.slipRateVar;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
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
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			int prevIndex = 0;
			for (double startTime=firstEventTime; startTime+windowSize<lastEventTime; startTime+=windowSize) {
				double endTime = startTime + windowSize;
				DiscretizedFunc subFunc = new ArbitrarilyDiscretizedFunc();
				double slipAtStart = 0d;
				for (int i=prevIndex; i<cumulativeSlipFunc.size(); i++) {
					double time = cumulativeSlipFunc.getX(i);
					if (time < startTime)
						continue;
					if (time > endTime)
						break;
					prevIndex = i;
					if (slipAtStart==0 && startTime > firstEventTime)
						slipAtStart = cumulativeSlipFunc.getY(i-1);
					double slip = cumulativeSlipFunc.getY(i) - slipAtStart;
					double relTime = time - startTime;
					Preconditions.checkState(relTime >= 0d);
					subFunc.set(relTime, slip);
				}
				funcs.add(subFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
			}
			
			ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
			oneToOne.set(0d, 0d);
			oneToOne.set(windowSize, slipRate*windowSize);
			
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, middle.getParentSectionName(), "Windowed Time (years)", "Cumulative Slip (mm/yr)");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec);
			
			String prefix = "slips_"+middle.getParentSectionName().replaceAll("\\W+", "_");
			PlotUtils.writePlots(new File("/tmp"), prefix, gp, 1000, 850, true, false, false);
		}
	}

}
