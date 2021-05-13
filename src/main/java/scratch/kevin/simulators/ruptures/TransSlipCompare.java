package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class TransSlipCompare {

	public static void main(String[] args) throws IOException {
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4950.instance();
//		File stderrFile = new File(catalog.getCatalogDir(), "write.706045.0.err");
//		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/singleSS"),
//				"Single SS", null, null);
//		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/singleSS_tri"),
//				"Single SS Triangles", null, null);
		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/bruce/rundirtest2"),
				"Bruce test", null, null);
		
		int maxEventID = -1;
//		int maxEventID = 41499;
		File stderrFile = null;
		
		Map<Integer, HashSet<Integer>> eventPinnedPatches = null;
		DefaultXY_DataSet pinnedAveSlipXY = null;
		MinMaxAveTracker pinnedAveSlipRatioTrack = null;
		DefaultXY_DataSet noPinAveSlipXY = null;
		MinMaxAveTracker noPinAveSlipRatioTrack = null;
		DefaultXY_DataSet futurePinAveSlipXY = null;
		MinMaxAveTracker futurePinAveSlipRatioTrack = null;
		HashSet<Integer> pinnedPatches = null;
		if (stderrFile != null) {
			eventPinnedPatches = loadPinnedEvents(stderrFile);
			pinnedAveSlipXY = new DefaultXY_DataSet();
			pinnedAveSlipRatioTrack = new MinMaxAveTracker();
			noPinAveSlipXY = new DefaultXY_DataSet();
			noPinAveSlipRatioTrack = new MinMaxAveTracker();
			
			pinnedPatches = new HashSet<>();
			for (HashSet<Integer> patches : eventPinnedPatches.values())
				pinnedPatches.addAll(patches);
			futurePinAveSlipXY = new DefaultXY_DataSet();
			futurePinAveSlipRatioTrack = new MinMaxAveTracker();
		}
		
		RSQSimStateTransitionFileReader trans = catalog.getTransitions();
		trans.setQuiet(true);
		
		DefaultXY_DataSet aveSlipXY = new DefaultXY_DataSet();
		MinMaxAveTracker aveSlipRatioTrack = new MinMaxAveTracker();
		
		int minRatioID = -1;
		int maxRatioID = -1;
		
		int numSlips = 0;
		double totalTimeSlipping = 0d;
		
//		List<RSQSimEvent> events = catalog.loader().minMag(6.5).skipYears(5000).hasTransitions().load();
//		for (RSQSimEvent event : events) {
//			if (event.getID() >= 163608) {
//				boolean found = false;
//				for (int id : event.getAllElementIDs())
//					found = found || id == 432203;
//				if (event.getID() == 163608)
//					System.out.println("Next event at t="+event.getTime()+", hasPatch ? "+found);
//				else if (found)
//					System.out.println("Finally participates in event "+event.getID()+" at "+event.getTime());
//				if (found)
//					break;
//			}
//		}
//		System.out.println("Loaded "+events.size()+" events");
//		for (RSQSimEvent e : events) {
		Loader loader;
		if (catalog.getFaultModel() != null)
			loader = catalog.loader().minMag(6.5).skipYears(5000);
		else
			loader = catalog.loader();
		if (maxEventID > 0)
			loader.maxEventID(maxEventID);
		for (RSQSimEvent e : loader.hasTransitions().iterable()) {
//			System.out.println("event "+e.getID()+" at year "+e.getTimeInYears());
			ArrayList<SimulatorElement> elems = e.getAllElements();
			double[] slips = e.getAllElementSlips();
			List<RSQSimStateTime> eventTrans = new ArrayList<>();
			double[] transSlips = calcTransSlips(catalog, e, trans, eventTrans);
			
			for (RSQSimStateTime eTrans : eventTrans) {
				if (eTrans.state == RSQSimState.EARTHQUAKE_SLIP) {
					numSlips++;
					totalTimeSlipping += eTrans.getDuration();
				}
			}
			
			double totArea = 0;
			double totSlipArea = 0;
			
			double totTransSlipArea = 0;
			
			for (int i=0; i<elems.size(); i++) {
				double area = elems.get(i).getArea(); // m^2
				totArea += area;
				
				totSlipArea += slips[i]*area;
				totTransSlipArea += transSlips[i]*area;
			}
			
			double aveSlip = totSlipArea/totArea;
			double transAveSlip = totTransSlipArea/totArea;
			
			aveSlipXY.set(aveSlip, transAveSlip);
			double slipRatio = transAveSlip/aveSlip;
			if (slipRatio > aveSlipRatioTrack.getMax())
				maxRatioID = e.getID();
			if (slipRatio < aveSlipRatioTrack.getMin())
				minRatioID = e.getID();
			aveSlipRatioTrack.addValue(slipRatio);
			
			if (eventPinnedPatches != null) {
				boolean pinned = eventPinnedPatches.containsKey(e.getID());
				
				if (pinned) {
					pinnedAveSlipXY.set(aveSlip, transAveSlip);
					pinnedAveSlipRatioTrack.addValue(slipRatio);
				} else {
					// see if any of the elements will be pinned
					boolean futurePin = false;
					for (int patchID : e.getAllElementIDs()) {
						if (pinnedPatches.contains(patchID)) {
							futurePin = true;
							break;
						}
					}
					if (futurePin) {
						futurePinAveSlipXY.set(aveSlip, transAveSlip);
						futurePinAveSlipRatioTrack.addValue(slipRatio);
					} else {
						noPinAveSlipXY.set(aveSlip, transAveSlip);
						noPinAveSlipRatioTrack.addValue(slipRatio);
					}
				}
			}
			
			if (aveSlipRatioTrack.getNum() % 1000 == 0) {
				System.out.println("Stats after "+aveSlipRatioTrack.getNum()+":");
				System.out.println("\tAve Slip Ratios: "+aveSlipRatioTrack);
			}
		}
		
		System.out.println("Final Stats after "+aveSlipRatioTrack.getNum()+":");
		System.out.println("\tAve Slip Ratios: "+aveSlipRatioTrack);
		System.out.println("Min ratio event ID: "+minRatioID);
		System.out.println("Max ratio event ID: "+maxRatioID);
		
		System.out.println("Total number of slips: "+numSlips);
		System.out.println("Total time slipping: "+totalTimeSlipping);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (eventPinnedPatches == null) {
			funcs.add(aveSlipXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		} else {
			funcs.add(noPinAveSlipXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.GRAY));
			
			funcs.add(futurePinAveSlipXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 2f, Color.GREEN));
			
			funcs.add(pinnedAveSlipXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 2f, Color.BLUE));
			
			System.out.println("No-pin stats:");
			System.out.println("\tAve Slip Ratios: "+noPinAveSlipRatioTrack);
			System.out.println("Future-pin stats:");
			System.out.println("\tAve Slip Ratios: "+futurePinAveSlipRatioTrack);
			System.out.println("Pinned stats:");
			System.out.println("\tAve Slip Ratios: "+pinnedAveSlipRatioTrack);
		}
		
		double maxSlip = Double.max(aveSlipXY.getMaxX(), aveSlipXY.getMaxY());
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Ave Slip Comparison", "From List Files", "From Transitions");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(0, Math.ceil(maxSlip));
		Range yRange = xRange;
		
		XY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(xRange.getLowerBound(), xRange.getLowerBound());
		oneToOne.set(xRange.getUpperBound(), xRange.getUpperBound());
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f,
				eventPinnedPatches == null ? Color.GRAY : Color.BLACK));
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File file = new File(catalog.getCatalogDir(), "trans_slip_compare");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	static double[] calcTransSlips(RSQSimCatalog catalog, RSQSimEvent event,
			RSQSimStateTransitionFileReader transReader) throws IOException {
		return calcTransSlips(catalog, event, transReader, null);
	}
	
	static double[] calcTransSlips(RSQSimCatalog catalog, RSQSimEvent event,
			RSQSimStateTransitionFileReader transReader, List<RSQSimStateTime> transToFillIn)
					throws IOException {
		ArrayList<SimulatorElement> elems = event.getAllElements();
		double[] transSlips = new double[elems.size()];
		
		Map<Integer, List<RSQSimStateTime>> trans = transReader.getTransitions(event, transToFillIn);
		
		for (int i=0; i<elems.size(); i++) {
			List<RSQSimStateTime> patchTrans = trans.get(elems.get(i).getID());
			for (int j=0; j<patchTrans.size(); j++) {
				RSQSimStateTime thisTrans = patchTrans.get(j);
				if (thisTrans.state == RSQSimState.EARTHQUAKE_SLIP) {
					double duration = thisTrans.getDuration();
					double vel = thisTrans.velocity;
					transSlips[i] += vel*duration;
				}
			}
		}
		return transSlips;
	}
	
	private static Map<Integer, HashSet<Integer>> loadPinnedEvents(File stderrFile) throws IOException {
		Map<Integer, HashSet<Integer>> ret = new HashMap<>();
		
		BufferedReader read = new BufferedReader(new FileReader(stderrFile));
		
		int eventID = -1;
		
		int patchCount = 0;
		
		String line;
		while ((line = read.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("event number") && line.contains("just started")) {
				String[] split = line.split(" ");
				eventID = Integer.parseInt(split[2])+1; // these are zero-based
			}
			// remove potential duplicate that doesn't have the patch ID
			line = line.replaceAll("Eliminating element Restressing", "");
			if (line.contains("Eliminating element")) {
				line = line.substring(line.indexOf("Eliminating element")+("Eliminating element".length())).trim();
				int patchID = Integer.parseInt(line.split(" ")[0]); // these are zero-based
				HashSet<Integer> patches = ret.get(eventID);
				if (patches == null) {
					patches = new HashSet<>();
					ret.put(eventID, patches);
				}
				patchCount++;
				patches.add(patchID);
			}
		}
		
		System.out.println(ret.size()+" events have pinned patches, "+patchCount+" patches in total");
		
		read.close();
		
		return ret;
	}

}
