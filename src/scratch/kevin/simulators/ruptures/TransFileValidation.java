package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class TransFileValidation {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		if (args.length == 1) {
			File dir = new File(args[0]);
			catalog = new RSQSimCatalog(dir, "Temp", "None", null, "Metadata", null, null);
		} else {
			catalog = Catalogs.BRUCE_2585.instance(new File("/home/kevin/Simulators/catalogs"));
		}
		
		double printThreshold = 2;
		
		int numZeroInTrans = 0;
		int numZeroInList = 0;
		SummaryStatistics stats = new SummaryStatistics();
		SummaryStatistics absDiffStats = new SummaryStatistics();
		int numAbove0p1 = 0;
		int numAbove1 = 0;
		long testCount = 0;
		
		List<SimulatorElement> allElems = catalog.getElements();
		Map<Integer, Double> patchAreas = new HashMap<>();
		for (SimulatorElement elem : allElems)
			patchAreas.put(elem.getID(), elem.getArea());
		
		SummaryStatistics momentStats = new SummaryStatistics();
		SummaryStatistics momentAbsDiffStats = new SummaryStatistics();
		
		Iterable<RSQSimEvent> eventsIt = RSQSimFileReader.getEventsIterable(catalog.getCatalogDir(), allElems);
		int count = 0;
		for (RSQSimEvent event : eventsIt) {
			if (count % 1000 == 0)
				System.out.println("Processing event: "+count);
			count++;
			RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(event);
			int[] elems = event.getAllElementIDs();
			double[] slips = event.getAllElementSlips();
			
			double listMoment = 0d;
			for (EventRecord rec : event)
				listMoment += rec.getMoment();
			
			double transMoment = 0d;
			
			double time = slipTimeFunc.getEndTime();
			for (int i=0; i<elems.length; i++) {
				transMoment += FaultMomentCalc.getMoment(patchAreas.get(elems[i]), slips[i]);
				testCount++;
				double transSlip = slipTimeFunc.getCumulativeEventSlip(elems[i], time);
				if (transSlip == 0 && slips[i] > 0) {
					numZeroInTrans++;
				} else if (transSlip > 0 && slips[i] == 0) {
					numZeroInList++;
				} else {
					double pDiff = DataUtils.getPercentDiff(transSlip, slips[i]);
					stats.addValue(pDiff);
					double diff = Math.abs(transSlip - slips[i]);
					absDiffStats.addValue(diff);
					if (pDiff >= 1d)
						numAbove1++;
					if (pDiff >= 0.1d)
						numAbove0p1++;
					if (pDiff > printThreshold) {
						System.out.println("Bad slip for event "+event.getID()+" (M="+(float)event.getMagnitude()+") at t="
								+event.getTime()+" ("+(float)event.getTimeInYears()+" yrs)");
						System.out.println("\tList File Slip:  "+(float)slips[i]);
						System.out.println("\tTrans File Slip: "+(float)transSlip);
						System.out.println("\t% Diff: "+(float)pDiff);
						System.out.println("\tAbs Diff: "+(float)diff);
					}
				}
			}
			double pDiff = DataUtils.getPercentDiff(transMoment, listMoment);
			momentStats.addValue(pDiff);
			momentAbsDiffStats.addValue(Math.abs(transMoment - listMoment));
			
			if (count % 10000 == 0) {
				System.out.println("=== SLIP STATS ===");
				System.out.println("Average % diff: "+(float)stats.getMean());
				System.out.println("Max % diff: "+(float)stats.getMax());
				System.out.println("Average diff: "+(float)absDiffStats.getMean());
				System.out.println("Max diff: "+(float)absDiffStats.getMax());
				System.out.println("=== MOMENT STATS ===");
				System.out.println("Average % diff: "+(float)momentStats.getMean());
				System.out.println("Max % diff: "+(float)momentStats.getMax());
				System.out.println("Average diff: "+(float)momentAbsDiffStats.getMean());
				System.out.println("Max diff: "+(float)momentAbsDiffStats.getMax());
				System.out.println("====================");
			}
		}
		double percentZeroTrans = 100d*(double)numZeroInTrans/(double)testCount;
		System.out.println(numZeroInTrans+"/"+testCount+" ("+(float)percentZeroTrans+" %) zero in trans, nonzero in list");
		double percentZeroList = 100d*(double)numZeroInList/(double)testCount;
		System.out.println(numZeroInList+"/"+testCount+" ("+(float)percentZeroList+" %) nonzero in trans, zero in list");
		System.out.println("Average % diff: "+(float)stats.getMean());
		System.out.println("Max % diff: "+(float)stats.getMax());
		System.out.println("% diff std. dev.: "+(float)stats.getStandardDeviation());
		System.out.println("Average diff: "+(float)absDiffStats.getMean());
		System.out.println("Max diff: "+(float)absDiffStats.getMax());
		double percentAbove1 = 100d*(double)numAbove1/(double)testCount;
		System.out.println(numAbove1+"/"+testCount+" ("+(float)percentAbove1+" %) above 1 %");
		double percentAbove0p1 = 100d*(double)numAbove0p1/(double)testCount;
		System.out.println(numAbove0p1+"/"+testCount+" ("+(float)percentAbove0p1+" %) above 0.1 %");
		System.out.println("");
		System.out.println("Trans moment ave % diff: "+(float)momentStats.getMean());
		System.out.println("Trans moment max % diff: "+(float)momentStats.getMax());
		System.out.println("Trans moment ave diff: "+(float)momentAbsDiffStats.getMean());
		System.out.println("Trans moment max diff: "+(float)momentAbsDiffStats.getMax());
	}

}
