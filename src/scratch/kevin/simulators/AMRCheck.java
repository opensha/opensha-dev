package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.QuietPeriodIdenMatcher;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.simulators.erf.SubSectionBiulder;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;

public class AMRCheck {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		List<SimulatorElement> elements = tools.getElementsList();
		System.out.println("Done loading events.");
		
		MagRangeRuptureIdentifier sevenPlusIden = new MagRangeRuptureIdentifier(7, 10);
		
		// we don't want aftershock sequences - crude approximation, find out no other M7's within 5 years before
		// -5 here because we'er actually filtering out foreshocks
		QuietPeriodIdenMatcher quietIden = new QuietPeriodIdenMatcher(sevenPlusIden, -5, 0, sevenPlusIden);
		
		List<? extends SimulatorEvent> matches = quietIden.getMatches(events);
		
		System.out.println("# Matches: "+matches.size());
		
		int daysBefore = 100;
		double distThresh = 100d;
		
		SubSectionBiulder subSectBuild = new SubSectionBiulder(elements);
		Map<Integer, Integer> elemToSubSectsMap = subSectBuild.getElemIDToSubSectsMap();
		List<FaultSectionPrefData> subSects = subSectBuild.getSubSectsList();
		
		EvenlyDiscretizedFunc beforeFunc = new EvenlyDiscretizedFunc(-100, daysBefore, 1d);
		
		List<List<LocationList>> tracesForEvents = Lists.newArrayList();
		for (int i=0; i<events.size(); i++) {
			SimulatorEvent event = events.get(i);
			List<LocationList> tracesForEvent = Lists.newArrayList();
			for (EventRecord rec : event) {
				HashSet<Integer> ssIDs = new HashSet<Integer>();
				for (SimulatorElement e : rec.getElements()) {
					ssIDs.add(elemToSubSectsMap.get(e.getID()));
				}
				List<Integer> ssIDsList = Lists.newArrayList(ssIDs);
				Collections.sort(ssIDsList);
				LocationList ll = new LocationList();
				for (int ssID : ssIDsList)
					ll.add(subSects.get(ssID).getFaultTrace().get(0));
				ll.add(subSects.get(ssIDsList.get(ssIDsList.size()-1)).getFaultTrace().get(1));
				tracesForEvent.add(ll);
			}
			Preconditions.checkState(!tracesForEvent.isEmpty());
			tracesForEvents.add(tracesForEvent);
		}
		
//		List<List<List<Location>>> centerLocsCache = Lists.newArrayList();
//		for (int i=0; i<events.size(); i++) {
//			List<List<Location>> locsList = Lists.newArrayList();
//			EQSIM_Event event = events.get(i);
//			for (EventRecord rec : event) {
//				List<Location> locs = Lists.newArrayList();
//				for (RectangularElement e : rec.getRectangularElements()) {
////					Preconditions.checkState(e.getNumDownDip() >= 0);
////					if (e.getNumDownDip() == 0)
//						locs.add(e.getCenterLocation());
//				}
//				Preconditions.checkState(!locs.isEmpty());
//				locsList.add(locs);
//			}
//			Preconditions.checkState(!locsList.isEmpty());
//			centerLocsCache.add(locsList);
//		}
		
		Table<Location, Location, Double> distsTable = HashBasedTable.create();
		
		int startInd = 0;
		
		for (int j=0; j<matches.size(); j++) {
			SimulatorEvent e = matches.get(j);
			if (j % 1000 == 0)
				System.out.println("Processing match "+j);
			double eventTime = e.getTimeInYears();
			double startTime = eventTime - (double)daysBefore * BatchPlotGen.DAYS_PER_YEAR;
			
//			List<Region> regions = Lists.newArrayList();
//			for (List<Location> centerLocsList1 : centerLocsCache.get(j)) {
//				LocationList ll = new LocationList();
//				ll.addAll(centerLocsList1);
//				Region sub = new Region(ll, distThresh);
//				regions.add(sub);
//			}
			
			List<Region> regions = Lists.newArrayList();
			for (LocationList traceForEvent : tracesForEvents.get(j)) {
				Region sub = new Region(traceForEvent, distThresh);
				regions.add(sub);
			}
			
			for (int i=startInd; i<events.size(); i++) {
				SimulatorEvent o = events.get(i);
				double oTime = o.getTimeInYears();
				if (oTime < startTime) {
					startInd = i;
					continue;
				} else if (oTime >= eventTime) {
					break;
				}
				double mag = o.getMagnitude();
				double distForSkip;
				if (mag < 6)
					distForSkip = distThresh * 1.2;
				else if (mag < 7)
					distForSkip = distThresh * 2;
				else
					distForSkip = Double.POSITIVE_INFINITY;
				
				// see if it's within the distance cutoff
				boolean within = false;
//				for (List<Location> centerLocs2 : centerLocsCache.get(i)) {
//					for (Region r : regions) {
//						for (Location loc2 : centerLocs2) {
//							if (r.contains(loc2)) {
//								within = true;
//								break;
//							}
//						}
//					}
//				}
				// see if we can skip
				if (regions.get(0).distanceToLocation(tracesForEvents.get(i).get(0).get(0)) < distForSkip) {
					for (LocationList traceForEvent : tracesForEvents.get(i)) {
						for (Region r : regions) {
							for (Location loc : traceForEvent) {
								if (r.contains(loc)) {
									within = true;
									break;
								}
							}
						}
					}
				}
//				outsideLoop:
//				for (Location loc1 : centerLocs1) {
//					for (Location loc2 : centerLocs2) {
//						if (loc1 == loc2) {
//							within = true;
//							break outsideLoop;
//						}
//						Double dist;
////						Double dist = distsTable.get(loc1, loc2);
////						if (dist == null) {
////							dist = distsTable.get(loc2, loc1);
////							if (dist == null) {
//								dist = LocationUtils.horzDistanceFast(
//										loc1, loc2);
////								distsTable.put(loc1, loc2, dist);
////							}
////						}
//						if (dist <= distThresh) {
//							within = true;
//							break outsideLoop;
//						}
//					}
//				}
				if (within) {
					int dayBin = (int)((oTime - startTime) / BatchPlotGen.DAYS_PER_YEAR);
					beforeFunc.add(dayBin, 1d);
				}
			}
		}
		beforeFunc.scale(1d/(double)matches.size());
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(beforeFunc);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		
		GraphWindow gw = new GraphWindow(funcs, "# Events Per Day Before Mainshock", chars);
	}

}
