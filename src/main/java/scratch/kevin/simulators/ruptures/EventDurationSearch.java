package scratch.kevin.simulators.ruptures;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class EventDurationSearch {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5685_SUB.instance();
		double minMag = 8.5;
		
		double minPrintSecs = 400;
		
		List<RSQSimEvent> events = catalog.loader().skipYears(2000).minMag(minMag).load();
		
		catalog.getTransitions().setQuiet(true);
		
		RSQSimEvent longest = null;
		double longestDuration = 0d;
		RSQSimEventSlipTimeFunc longestSlipTime = null;
		
		double[] thresholds = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
		int[] threshCounts = new int[thresholds.length];
		
		for (RSQSimEvent event : events) {
			RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(event);
			
			double firstSlip = slipTimeFunc.getStartTime();
			double lastSlip = slipTimeFunc.getEndTime();
			double duration = lastSlip - firstSlip;
			if (duration > minPrintSecs)
				System.out.println(eventDurStr(event, duration, slipTimeFunc));
			for (int i=0; i<thresholds.length; i++)
				if (duration > thresholds[i])
					threshCounts[i]++;
			if (duration > longestDuration) {
				longestDuration = duration;
				longest = event;
				longestSlipTime = slipTimeFunc;
			}
		}
		System.out.println("Longest event:\n"+eventDurStr(longest, longestDuration, longestSlipTime));
		for (int i=0; i<threshCounts.length; i++) {
			if (thresholds[i] > longestDuration)
				break;
			System.out.println("Events over "+(int)thresholds[i]+" s ("+twoDigits.format(thresholds[i]/60d)+" m):\t"+threshCounts[i]);
		}
	}
	
	private static final DecimalFormat twoDigits = new DecimalFormat("0.00");
	private static String eventDurStr(RSQSimEvent event, double duration, RSQSimEventSlipTimeFunc slipTimeFunc) {
		String str = "event "+event.getID()+", M"+twoDigits.format(event.getMagnitude())+", "
				+twoDigits.format(duration)+" s = "+twoDigits.format(duration/60d)+" m";
		
		Set<Integer> patches = slipTimeFunc.getPatchIDs();
		Map<Integer, Double> areas = new HashMap<>(patches.size());
		for (SimulatorElement elem : event.getAllElements())
			areas.put(elem.getID(), elem.getArea());
		double totMoment = 0d;
		double minTime = slipTimeFunc.getStartTime();
		double maxTime = slipTimeFunc.getEndTime();
		for (int patchID : patches)
			totMoment += FaultMomentCalc.getMoment(areas.get(patchID), slipTimeFunc.getCumulativeEventSlip(patchID, maxTime));
		
		String timeLine = "";
		String fractLine = "";
		for (double thresh=50d; thresh<duration; thresh += 50d) {
			double time = minTime + thresh;
			double timeMoment = 0;
			for (int patchID : patches)
				timeMoment += FaultMomentCalc.getMoment(areas.get(patchID), slipTimeFunc.getCumulativeEventSlip(patchID, time));
			double fract = timeMoment/totMoment;
			timeLine += "\t"+(int)(thresh+0.5)+" s";
			fractLine += "\t"+twoDigits.format(fract);
		}
		str += ";\tMoment Fractions:\n"+timeLine+"\n"+fractLine;
		return str;
	}

}
