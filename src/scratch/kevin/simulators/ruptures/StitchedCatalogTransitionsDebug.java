package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;

public class StitchedCatalogTransitionsDebug {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/data/kevin/simulators/catalogs/rundir2585_1myr/stitch_debug");
		File stitchedDir = new File("/data/kevin/simulators/catalogs/rundir2585_1myr");
		File beforeDir = new File("/data/kevin/simulators/catalogs/restart2585");
		File afterDir = new File("/data/kevin/simulators/catalogs/rundir2585extend");
		
		double stitchTime = 26772011003374.0078125d;
//		double minMag = 6.5;
		double minMag = 0d;
		
		double bufferYears = 50d;
		double tStart = stitchTime - bufferYears*SimulatorUtils.SECONDS_PER_YEAR;
		double tEnd = stitchTime + bufferYears*SimulatorUtils.SECONDS_PER_YEAR;
		
		RSQSimCatalog stitchedCatalog = new RSQSimCatalog(stitchedDir, "Stitched", null, null, null, null, null);
		RSQSimCatalog beforeCatalog = new RSQSimCatalog(beforeDir, "Before", null, null, null, null, null);
		RSQSimCatalog afterCatalog = new RSQSimCatalog(afterDir, "After", null, null, null, null, null);
		
		System.out.println("Loading stitched catalog");
		List<RSQSimEvent> stitchedEvents = stitchedCatalog.loader().minMag(minMag).withinTimeRange(tStart, tEnd).load();
		System.out.println("Loaded "+stitchedEvents.size()+" events");
		System.out.println("Loading before catalog");
		List<RSQSimEvent> beforeEvents = beforeCatalog.loader().minMag(minMag).withinTimeRange(tStart, tEnd).load();
		System.out.println("Loaded "+beforeEvents.size()+" events");
		System.out.println("Loading after catalog");
		List<RSQSimEvent> afterEvents = afterCatalog.loader().minMag(minMag).withinTimeRange(tStart, tEnd).load();
		System.out.println("Loaded "+afterEvents.size()+" events");
		
		List<String> lines = new ArrayList<>();
		for (RSQSimEvent event : stitchedEvents) {
			double stitchDeltaYears = (event.getTime() - stitchTime) / SimulatorUtils.SECONDS_PER_YEAR;
			lines.add((float)stitchDeltaYears+" yrs");
			lines.add("\tStitch:\t"+getEventString(stitchedCatalog, event, -1));
			lines.add("\tBefore:\t"+getEventString(beforeCatalog, findEvent(event, beforeEvents), event.getTime()));
			lines.add("\tAfter:\t"+getEventString(afterCatalog, findEvent(event, afterEvents), event.getTime()));
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File outputFile;
		if (minMag > 0)
			outputFile = new File(outputDir, "event_match_m"+(float)+minMag+".txt");
		else
			outputFile = new File(outputDir, "event_match_all.txt");
		FileWriter fw = new FileWriter(outputFile);
		
		for (String line : lines) {
			System.out.println(line);
			fw.write(line+"\n");
		}
		
		fw.close();
		
		// now plot trans vs time
		List<RSQSimStateTime> stitchedTrans = getTransitions(stitchedCatalog, tStart, tEnd);
		List<RSQSimStateTime> beforeTrans = getTransitions(beforeCatalog, tStart, tEnd);
		List<RSQSimStateTime> afterTrans = getTransitions(afterCatalog, tStart, tEnd);
		
		for (boolean cumulative : new boolean[] {false, true}) {
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			DiscretizedFunc stitchedFunc, beforeFunc, afterFunc;
			String title, yAxisLabel, prefix;
			if (cumulative) {
				stitchedFunc = getCumulativeTransFunc(stitchedTrans, stitchTime);
				beforeFunc = getCumulativeTransFunc(beforeTrans, stitchTime);
				afterFunc = getCumulativeTransFunc(afterTrans, stitchTime);
				yAxisLabel = "Cumulative Transitions Count";
				title = "Cumulative Transitions Debug";
				prefix = "trans_count_cumulative";
			} else {
				stitchedFunc = getBinnedIncrTransFunc(stitchedTrans, stitchTime, tStart, tEnd, 500);
				beforeFunc = getBinnedIncrTransFunc(beforeTrans, stitchTime, tStart, tEnd, 500);
				afterFunc = getBinnedIncrTransFunc(afterTrans, stitchTime, tStart, tEnd, 500);
				yAxisLabel = "Binned Transitions Count";
				title = "Binned Transitions Debug";
				prefix = "trans_count_binned";
			}
			
			beforeFunc.setName("Before Catalog");
			funcs.add(beforeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
			
			afterFunc.setName("After Catalog");
			funcs.add(afterFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));

			stitchedFunc.setName("Stitched Catalog");
			funcs.add(stitchedFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK));
			
			double maxY = Math.max(stitchedFunc.getMaxY(), Math.max(beforeFunc.getMaxY(), afterFunc.getMaxY()));
			maxY *= 1.1;
			
			DefaultXY_DataSet stitchLine = new DefaultXY_DataSet();
			stitchLine.set(0d, 0d);
			stitchLine.set(0d, maxY);
			
			funcs.add(stitchLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, Color.RED));
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Δ Time (s)", yAxisLabel);
			spec.setLegendVisible(true);
			
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			
			Range xRange = new Range(tStart-stitchTime, tEnd-stitchTime);
			Range yRange = new Range(0, maxY);
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(600, 500);
			File pngFile = new File(outputDir, prefix+".png");
			gp.saveAsPNG(pngFile.getAbsolutePath());
			
			GraphWindow gw = new GraphWindow(new GraphWidget(spec, plotPrefs, false, false, xRange, yRange));
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
	}
	
	private static String getEventString(RSQSimCatalog catalog, RSQSimEvent event, double targetTime) throws IOException {
		if (event == null)
			return "(missing)";
		Map<Integer, List<RSQSimStateTime>> trans = catalog.getTransitions().getTransitions(event);
		int numTrans = 0;
		for (List<RSQSimStateTime> states : trans.values())
			numTrans += states.size();
		String str = "ID="+event.getID()+"\tM="+(float)event.getMagnitude()+"\tt="+event.getTime()+"\t#Trans="+numTrans;
		if (targetTime >= 0)
			str += "\tΔt="+(float)(event.getTime()-targetTime);
		return str;
	}
	
	private static RSQSimEvent findEvent(RSQSimEvent target, List<RSQSimEvent> catalog) {
		double mag = target.getMagnitude();
		Map<Integer, Double> patchSlipMap = new HashMap<>();
		int[] ids = target.getAllElementIDs();
		double[] slips = target.getAllElementSlips();
		for (int i=0; i<ids.length; i++) {
			Preconditions.checkState(!patchSlipMap.containsKey(ids[i]));
			patchSlipMap.put(ids[i], slips[i]);
		}
		for (RSQSimEvent e : catalog) {
			if ((float)mag == (float)e.getMagnitude()) {
				ids = target.getAllElementIDs();
				slips = target.getAllElementSlips();
				if (ids.length != patchSlipMap.size())
					continue;
				for (int i=0; i<ids.length; i++) {
					if (!patchSlipMap.containsKey(ids[i]))
						continue;
					if ((float)slips[i] != patchSlipMap.get(ids[i]).floatValue())
						continue;
				}
				return e;
			}
//			if (id == e.getID() && (float)mag == (float)e.getMagnitude())
//				return e;
//			if ((float)mag == (float)e.getMagnitude() && (float)time == (float)e.getTime())
//				return e;
		}
		return null;
	}
	
	private static List<RSQSimStateTime> getTransitions(RSQSimCatalog catalog, double tStart, double tEnd) throws IOException {
		List<RSQSimStateTime> trans = catalog.getTransitions().getTransitions(tStart, tEnd);
		if (trans.isEmpty())
			return trans;
		int numInitializers = 0;
		double tFirst = trans.get(0).getStartTime();
		for (int i=0; i<trans.size(); i++) {
			RSQSimStateTime st = trans.get(i);
			if (st.getStartTime() == tFirst && st.getState() == RSQSimState.LOCKED)
				numInitializers++;
			else
				break;
		}
		if (numInitializers > 100) {
			// this is the start of the file, initializing all to state 0. strip those out
			System.out.println("Stripping out "+numInitializers+" transition initializations");
			trans = trans.subList(numInitializers, trans.size());
		}
		return trans;
	}
	
	private static LightFixedXFunc getCumulativeTransFunc(List<RSQSimStateTime> trans, double stitchTime) {
		double[] times = new double[trans.size()];
		double[] counts = new double[trans.size()];
		
		boolean inOrder = true;
		for (int i=0; i<trans.size(); i++) {
			RSQSimStateTime st = trans.get(i);
			times[i] = st.getStartTime()-stitchTime;
			counts[i] = i+1;
			if (i > 0 && times[i]<times[i-1])
				inOrder = false;
		}
		if (inOrder)
			System.out.println("Transitions are in order!");
		else
			System.out.println("WARNING: transitions not in order!");
		
		return new LightFixedXFunc(times, counts);
	}

	private static EvenlyDiscretizedFunc getBinnedIncrTransFunc(List<RSQSimStateTime> trans,
			double stitchTime, double tStart, double tEnd, int numBins) {
		EvenlyDiscretizedFunc bins = new EvenlyDiscretizedFunc(tStart-stitchTime, tEnd-stitchTime, numBins);
		
		for (int i=0; i<trans.size(); i++) {
			RSQSimStateTime st = trans.get(i);
			int index = bins.getClosestXIndex(st.getStartTime()-stitchTime);
			bins.add(index, 1d);
		}
		
		return bins;
	}

}
