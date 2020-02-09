package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class NegativeSlipDebugPageGen {

	public static void main(String[] args) throws IOException {
		File gitDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_4660.instance();
		
		File catalogOutputDir = new File(gitDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File negSlipDir = new File(catalogOutputDir, "negative_slip_debug");
		Preconditions.checkState(negSlipDir.exists() || negSlipDir.mkdir());
		
		File resourcesDir = new File(negSlipDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		Map<SimulatorElement, Integer> negSlipCounts = new HashMap<>();
		Map<SimulatorElement, Integer> negEventDisplacementCounts = new HashMap<>();
		
		RSQSimStateTransitionFileReader transFile = catalog.getTransitions();
		Preconditions.checkState(transFile.isVariableSlipSpeed());
		
		List<RSQSimEvent> negDisplEvents = new ArrayList<>();
		
		double minTime = Double.POSITIVE_INFINITY;
		double maxTime = Double.NEGATIVE_INFINITY;
		
		int totalEvents = 0;
		for (RSQSimEvent event : catalog.loader().skipYears(2000).iterable()) {
			double time = event.getTime();
			minTime = Math.min(time, minTime);
			maxTime = Math.max(time, maxTime);
			totalEvents++;
			List<SimulatorElement> elems = event.getAllElements();
			double[] slips = event.getAllElementSlips();
			int numFound = 0;
			for (int i=0; i<slips.length; i++) {
				if (slips[i] < 0) {
					Integer prevCount = negEventDisplacementCounts.get(elems.get(i));
					if (prevCount == null)
						prevCount = 0;
					negEventDisplacementCounts.put(elems.get(i), prevCount+1);
					numFound++;
				}
			}
			if (numFound > 0)
				negDisplEvents.add(event);
		}
		
		System.out.println("Found "+negDisplEvents.size()+" events with a negative displacement (on "
				+negEventDisplacementCounts.size()+" elements)");
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Negative Slip Debug");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Summary Statistics");
		lines.add(topLink); lines.add("");
		
		int numElems = catalog.getElements().size();
		
		float percent = 100f*negDisplEvents.size()/(float)totalEvents;
		lines.add("* "+negDisplEvents.size()+"/"+totalEvents+" events have at least one total negative displacement ("+percent+" %)");
		percent = 100f*negEventDisplacementCounts.size()/(float)numElems;
		lines.add("* "+negEventDisplacementCounts.size()+" elements have at least one total negative displacement ("+percent+" %). "
				+ "Stats of those elements:");
		MinMaxAveTracker negDisplTrack = stats(negEventDisplacementCounts);
		lines.add("  * min="+(int)negDisplTrack.getMin()+", max="+(int)negDisplTrack.getMax()
			+", mean="+(float)negDisplTrack.getAverage());
		
		for (double startTime=minTime; startTime<maxTime; startTime += SimulatorUtils.SECONDS_PER_YEAR) {
			double endTime = startTime + SimulatorUtils.SECONDS_PER_YEAR;
			List<RSQSimStateTime> transList = transFile.getTransitions(startTime, endTime);
			for (RSQSimStateTime trans : transList) {
				
			}
		}
		
		String negDisplPrefix = "negative_displacements";
		List<SimulatorElement> scaledElems = new ArrayList<>();
		for (SimulatorElement elem : negEventDisplacementCounts.keySet())
			scaledElems.add(elem);
		double[] scalars = new double[scaledElems.size()];
		for (int i=0; i<scalars.length; i++)
			scalars[i] = negEventDisplacementCounts.get(scaledElems.get(i)).doubleValue();
		CPT scalarCPT = new CPT(1d, 10d, Color.BLACK, Color.BLACK);
		RupturePlotGenerator.writeMapPlot(catalog.getElements(), null, null, 
				resourcesDir, negDisplPrefix, null, null, null, scaledElems, scalars, scalarCPT, " ", null);
		
		lines.add("");
		
		lines.add("## Negative Displacement Map");
		lines.add(topLink); lines.add("");
		
		lines.add("![map](resources/"+negDisplPrefix+".png)");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, negSlipDir);
	}
	
	private static MinMaxAveTracker stats(Map<SimulatorElement, Integer> countsMap) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (SimulatorElement elem : countsMap.keySet())
			track.addValue(countsMap.get(elem));
		return track;
	}

}
