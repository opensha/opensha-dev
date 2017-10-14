package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.FaultIDIden;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RuptureSearch {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2273.instance(baseDir);
		
		RSQSimUtils.populateFaultIDWithParentIDs(catalog.getElements(), catalog.getU3SubSects());
		
		int parentID = 301; // Mojave S
		double minMag = 7;
		
		int numToPlot = 5;
		File outputDir = new File(catalog.getCatalogDir(), "search_events");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		double aveElemArea = catalog.getAveArea();
		double aveSectArea = 0d;
		int numSects = 0;
		for (FaultSectionPrefData sect : catalog.getU3SubSects()) {
			if (sect.getParentSectionId() == parentID) {
				aveSectArea += sect.getTraceLength() * sect.getOrigDownDipWidth();
				numSects++;
			}
		}
		aveSectArea /= numSects;
		double elementOffPenalty = 20d * aveElemArea / aveSectArea;
		System.out.print("Stray penalty (each): "+elementOffPenalty);
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).matches(
				new FaultIDIden("FID", catalog.getElements(), parentID)).load();
		
		System.out.println("Found "+events.size()+" events with M>"+minMag+" on parent "+parentID);
		
		List<EventScore> scores = new ArrayList<>();
		
		for (RSQSimEvent event : events) {
			List<FaultSectionPrefData> sects = catalog.getSubSectsForRupture(
					event, RSQSimBBP_Config.MIN_SUB_SECT_FRACT);
			int numOn = 0;
			int numOff = 0;
			for (FaultSectionPrefData sect : sects) {
				if (sect.getParentSectionId() == parentID)
					numOn++;
				else
					numOff++;
			}
			int numStrays = 0;
			double strayElemPenalty = 0d;
			for (SimulatorElement elem : event.getAllElements()) {
				if (elem.getFaultID() != parentID) {
					strayElemPenalty += elementOffPenalty;
					numStrays++;
				}
			}
			scores.add(new EventScore(event, numOn, numOff, strayElemPenalty, numStrays));
		}
		
		Collections.sort(scores);
		
		for (int i=0; i<numToPlot && i<scores.size(); i++) {
			EventScore score = scores.get(i);
			RSQSimEvent event = score.event;
			System.out.print(i+". Event "+event.getID()+", M="+(float)event.getMagnitude());
			System.out.println("\tscore = "+score.score()+" = "+score.numOn+" on - "+score.numOff+" off - "
					+(float)score.strayElemPenalty+" stray penalty ("+score.numStrays+" strays)");
			
			RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
			
			String prefix = "search_"+parentID+"_num"+i+"_score"+(float)score.score()+"_event_"+event.getID();
			RupturePlotGenerator.writeSlipPlot(event, func, outputDir, prefix);
			RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, func, outputDir, prefix+"_map");
		}
	}
	
	private static class EventScore implements Comparable<EventScore> {
		RSQSimEvent event;
		int numOn;
		int numOff;
		double strayElemPenalty;
		int numStrays;
		
		public EventScore(RSQSimEvent event, int numOn, int numOff, double strayElemPenalty, int numStrays) {
			this.event = event;
			this.numOn = numOn;
			this.numOff = numOff;
			this.strayElemPenalty = strayElemPenalty;
			this.numStrays = numStrays;
		}
		
		public double score() {
			return numOn - numOff - strayElemPenalty;
		}

		@Override
		public int compareTo(EventScore o) {
			return Double.compare(o.score(), score());
		}
	}

}
