package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_A_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Somerville_2006_MagAreaRel;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.EventRecord;
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

	public static void main(String[] args) throws IOException, DocumentException {
		File mainOutputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_5044.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(baseDir);
		
		RSQSimUtils.populateFaultIDWithParentIDs(catalog.getElements(), catalog.getU3SubSects());
		
		int parentID = 301; // Mojave S
		String parentName = "SAF South Mojave";
		double minMag = 7.0;
		double maxMag = 7.5;
		double maxMagDiff = 0.3;
//		MagAreaRelationship magAreaRel = new Somerville_2006_MagAreaRel();
		MagAreaRelationship magAreaRel = new Ellsworth_A_WG02_MagAreaRel();
//		Region hypocenterReg = new Region(new Location(34.25, -117.5), 20d);
		Region hypocenterReg = null;
		
		File catalogOutputDir = new File(mainOutputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File outputDir = new File(catalogOutputDir, "search_events");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		// header
		List<String> lines = new ArrayList<>();
		lines.add("# "+catalog.getName()+" Event Search");
		lines.add("");
		lines.add("This page lists and plots candidate events which match the following criteria:");
		lines.add("");
		lines.add("* Almostly exclusively on "+parentName);
		lines.add("* Magnitude in the range ["+(float)minMag+","+(float)maxMag+"]");
		if (maxMagDiff > 0)
			lines.add("* Magnitude no more than "+(float)maxMagDiff
					+" different from that predicted by "+magAreaRel.getName());
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		int numToPlot = 10;
		double aveElemArea = catalog.getAveArea();
		double aveSectArea = 0d;
		int numSects = 0;
		for (FaultSection sect : catalog.getU3SubSects()) {
			if (sect.getParentSectionId() == parentID) {
				aveSectArea += sect.getTraceLength() * sect.getOrigDownDipWidth();
				numSects++;
			}
		}
		aveSectArea /= numSects;
		double elementOffPenalty = 50d * aveElemArea / aveSectArea;
		System.out.println("Stray penalty (each): "+elementOffPenalty);
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).maxMag(maxMag)
				.forParentSections(true, parentID).load();
		
		System.out.println("Found "+events.size()+" events with M>"+minMag+" on parent "+parentID);
		
		List<EventScore> scores = new ArrayList<>();
		
		for (RSQSimEvent event : events) {
			if (hypocenterReg != null) {
				Location hypoLoc = null;
				double minTime = Double.POSITIVE_INFINITY;
				for (EventRecord e : event) {
					List<SimulatorElement> elems = e.getElements();
					double[] times = e.getElementTimeFirstSlips();
					for (int i=0; i<elems.size(); i++) {
						if (times[i] < minTime) {
							minTime = times[i];
							hypoLoc = elems.get(i).getCenterLocation();
						}
					}
				}
				if (!hypocenterReg.contains(hypoLoc))
					continue;
			}
			if (maxMagDiff > 0) {
				double calcMag = magAreaRel.getMedianMag(event.getArea() * 1e-6);
				double magDiff = Math.abs(calcMag - event.getMagnitude());
				System.out.println("M"+(float)event.getMagnitude()+" vs M/A M"+(float)calcMag+", diff="+(float)magDiff);
				if (magDiff > maxMagDiff)
					continue;
			}
			List<FaultSection> sects = catalog.getSubSectsForRupture(event);
			int numOn = 0;
			int numOff = 0;
			for (FaultSection sect : sects) {
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
		
		boolean hasTrans = false;
		try {
			catalog.getTransitions();
			hasTrans = true;
			System.out.println("We have transitions");
		} catch (Exception e1) {}
		
		for (int i=0; i<numToPlot && i<scores.size(); i++) {
			EventScore score = scores.get(i);
			RSQSimEvent event = score.event;
			double calcMag = magAreaRel.getMedianMag(event.getArea() * 1e-6);
			System.out.print(i+". Event "+event.getID()+", M="+(float)event.getMagnitude()+", M/A calc M="+(float)calcMag);
			System.out.println("\tscore = "+score.score()+" = "+score.numOn+" on - "+score.numOff+" off - "
					+(float)score.strayElemPenalty+" stray penalty ("+score.numStrays+" strays)");
			double aveRupVel = calcAvgRupVel(event);
			System.out.println("\taverage velocity: "+(float)aveRupVel);
			
			RSQSimEventSlipTimeFunc func = null;
			if (hasTrans)
				func = catalog.getSlipTimeFunc(event);
			
			String prefix = "search_"+parentID+"_num"+i+"_score"+(float)score.score()+"_event_"+event.getID()
				+"_m"+(float)event.getMagnitude()+"_calc_mag_m"+(float)calcMag;
			RupturePlotGenerator.writeSlipPlot(event, func, resourcesDir, prefix);
			RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, func, resourcesDir, prefix+"_map");
			
			lines.add("## Event "+event.getID()+", M"+(float)event.getMagnitude());
			lines.add(topLink); lines.add("");
			lines.add("* Event Time: "+(float)event.getTimeInYears()+" (yrs)");
			lines.add("* Area: "+event.getArea()+" m^2");
			lines.add("* "+magAreaRel.getName()+" mag: "+(float)calcMag);
			lines.add("* Average rupture velocity: "+(float)aveRupVel);
			lines.add("* Match score: "+(float)score.score()+" = "+score.numOn+" on - "+score.numOff+" off - "
					+(float)score.strayElemPenalty+" stray penalty ("+score.numStrays+" strays)");
			lines.add("");
			lines.add("![map plot](resources/"+prefix+"_map.png)");
			lines.add("");
			lines.add("![slip plot](resources/"+prefix+".png)");
			
			if (hasTrans) {
				// add animation
				File rupAnim = new File(resourcesDir, prefix+"_anim.gif");
				RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 10);
				lines.add("");
				lines.add("![slip anim](resources/"+rupAnim.getName()+")");
			}
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
	
	private static double calcAvgRupVel(RSQSimEvent event) {
		double aveVel = 0d;
		
		double minTime = Double.POSITIVE_INFINITY;
		int hypoID = -1;
		Location hypo = null;
		for (EventRecord rec : event) {
			double[] times = rec.getElementTimeFirstSlips();
			List<SimulatorElement> elems = rec.getElements();
			Preconditions.checkNotNull(times, "Event doesn't have timing information");
			int[] ids = rec.getElementIDs();
			
			for (int i=0; i<ids.length; i++) {
				if (times[i] < minTime) {
					minTime = times[i];
					hypoID = ids[i];
					hypo = elems.get(i).getCenterLocation();
				}
			}
		}
		
		int num = 0;
		
		for (EventRecord rec : event) {
			double[] times = rec.getElementTimeFirstSlips();
			int[] ids = rec.getElementIDs();
			List<SimulatorElement> elems = rec.getElements();
			
			for (int i=0; i<ids.length; i++) {
				if (ids[i] != hypoID) {
					double dist = LocationUtils.linearDistanceFast(hypo, elems.get(i).getCenterLocation());
					double tDelta = times[i] - minTime;
					if (tDelta == 0)
						continue;
					double vel = dist/(tDelta);
					Preconditions.checkState(Double.isFinite(vel) && vel > 0,
							"Bad velocity! vel = %s / %s = %s", dist, tDelta, vel);
					aveVel += vel;
					num++;
				}
			}
		}
		
		if (num > 0)
			aveVel /= num;
		return aveVel;
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
