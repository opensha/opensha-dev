package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class TransVsListSlipDebug {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4950.instance();
		int eventID = 368122;
		int patchID = -1;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
//		int eventID = 8324165;
//		int patchID = -1;
		
//		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/singleSS"),
//				"Single SS", null, null);
//		int eventID = -1;
//		int patchID = -1;
		
		RSQSimStateTransitionFileReader transReader = catalog.getTransitions();
		transReader.setQuiet(true);
		
		RSQSimEvent event;
		if (eventID < 0) {
			event = findLargestEventAveSlipRatio(catalog, transReader);
			eventID = event.getID();
			patchID = -1;
		} else {
			event = catalog.loader().byID(eventID);
		}
		
		if (patchID < 0) {
			patchID = findLargestPatchDelta(catalog, event, transReader);
		}
		
		double eventTime = event.getTime();
		
		double patchStartTime = Double.NaN;
//		double patchEndTime = Double.NaN;
		double patchListSlip = Double.NaN;
		for (EventRecord rec : event) {
			Preconditions.checkState(rec instanceof RSQSimEventRecord);
			RSQSimEventRecord rsRecord = (RSQSimEventRecord)rec;
			int[] recPatchIDs = rsRecord.getElementIDs();
			Preconditions.checkNotNull(recPatchIDs);
			double[] startSlipTimes = rsRecord.getElementTimeFirstSlips();
//			double[] nextSlipTimes = rsRecord.getNextSlipTimes();
//			Preconditions.checkNotNull(nextSlipTimes);
			double[] slips = rsRecord.getElementSlips();
			for (int i=0; i<recPatchIDs.length; i++) {
				if (recPatchIDs[i] == patchID) {
					patchStartTime = startSlipTimes[i];
//					patchEndTime = nextSlipTimes[i];
					patchListSlip = slips[i];
				}
			}
		}

		System.out.println("Event occurs at "+eventTime);
		List<RSQSimStateTime> eventTrans = new ArrayList<>();
		transReader.getTransitions(event, eventTrans);
		double eventEndTime = eventTrans.get(eventTrans.size()-1).getStartTime();
		System.out.println("Event ends at "+eventEndTime);
		System.out.println("Event duration: "+(eventEndTime-eventTime));
		System.out.println("Next event is at "+event.getNextEventTime());
		System.out.println("\twhich is "+(event.getNextEventTime()-eventEndTime)+" s after end");
		
		double startTime = eventTime - 3600; // 1h before
		double endTime = eventEndTime + 3600; // 1h after
		List<RSQSimStateTime> allTrans = transReader.getTransitions(startTime, endTime);

		System.out.println("List file patch times:");
		System.out.println("\tPatch start time: "+patchStartTime);
//		System.out.println("\tPatch next time: "+patchEndTime);
//		System.out.println("\tPatch duration: "+(patchEndTime - patchStartTime));
		
		double patchVel = transReader.isVariableSlipSpeed() ?
				Double.NaN : catalog.getSlipVelocities().get(patchID);
		
		boolean before = true;
		boolean during = false;
		
		System.out.println();
		System.out.println("Transitions Listing");
		System.out.println("*************** BEFORE ***************");
		
		double beforeSlip = 0;
		double duringSlip = 0;
		double afterSlip = 0;
		
		List<RSQSimStateTime> patchTrans = new ArrayList<>();
		for (RSQSimStateTime trans : allTrans)
			if (trans.getPatchID() == patchID)
				patchTrans.add(trans);
		
		for (int i=0; i<patchTrans.size(); i++) {
			RSQSimStateTime trans = patchTrans.get(i);
			
			double time = trans.getStartTime();
			double relTime = time - eventTime;
			
			if (before && time >= eventTime) {
				before = false;
				during = true;
				System.out.println("*************** END BEFORE ***************");
				System.out.println("Before slip: "+beforeSlip);
				System.out.println("*************** DURING ***************");
			}
			if (during && time > eventEndTime) {
				during = false;
				System.out.println("*************** END DURING ***************");
				System.out.println("During slip: "+duringSlip);
				System.out.println("*************** AFTER ***************");
			}
			
			double slip = Double.NaN;
			double myEndTime = trans.getEndTime();
			double duration = trans.getDuration();
			if (trans.getState() == RSQSimState.EARTHQUAKE_SLIP) {
				double vel = transReader.isVariableSlipSpeed() ? trans.getVelocity() : patchVel;
				slip = vel*duration;
				
				if (before)
					beforeSlip += slip;
				else if (during)
					duringSlip += slip;
				else
					afterSlip += slip;
			}
			String str = relTime+" => "+(myEndTime-eventTime)+" ("+duration+"): "+trans.getState();
			if (trans.getState() == RSQSimState.EARTHQUAKE_SLIP) {
				str += "\tslip="+slip;
				if (transReader.isVariableSlipSpeed())
					str += "\tvel="+trans.getVelocity();
			}
			System.out.println(str);
		}
		if (during) {
			System.out.println("*************** END DURING ***************");
			System.out.println("During slip: "+duringSlip);
			System.out.println("*************** AFTER ***************");
		}
		System.out.println("*************** END AFTER ***************");
		System.out.println("After slip: "+afterSlip);
		System.out.println("*************************************");
		
		System.out.println("List file slip: "+patchListSlip);
	}
	
	private static RSQSimEvent findLargestEventAveSlipRatio(RSQSimCatalog catalog,
			RSQSimStateTransitionFileReader transReader) throws IOException {
		double maxRatio = 1;
		RSQSimEvent maxEvent = null;
		
		Loader loader;
		if (catalog.getFaultModel() != null)
			loader = catalog.loader().minMag(6.5).skipYears(5000);
		else
			loader = catalog.loader();
		for (RSQSimEvent e : loader.iterable()) {
			ArrayList<SimulatorElement> elems = e.getAllElements();
			double[] slips = e.getAllElementSlips();
			double[] transSlips = TransSlipCompare.calcTransSlips(catalog, e, transReader);
			
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
			
			double ratio = aveSlip > transAveSlip ? aveSlip / transAveSlip : transAveSlip / aveSlip;
			if (ratio > maxRatio) {
				maxRatio = ratio;
				maxEvent = e;
			}
		}
		
		System.out.println("Max event: "+maxEvent.getID()+" with ratio="+(float)maxRatio);
		
		return maxEvent;
	}
	
	private static int findLargestPatchDelta(RSQSimCatalog catalog, RSQSimEvent event,
			RSQSimStateTransitionFileReader transReader) throws IOException {
		int[] patchIDs = event.getAllElementIDs();
		double[] slips = event.getAllElementSlips();
		double[] transSlips = TransSlipCompare.calcTransSlips(catalog, event, transReader);
		
		double max = 0d;
		int maxID = -1;
		double slip = Double.NaN;
		double transSlip = Double.NaN;
		for (int i=0; i<slips.length; i++) {
			double diff = Math.abs(slips[i]-transSlips[i]);
			if (diff > max) {
				max = diff;
				maxID = patchIDs[i];
				slip = slips[i];
				transSlip = transSlips[i];
			}
		}
		
		System.out.println("Max patch mismatch of "+(float)max+" on patch "+maxID
				+". slip="+slip+", transSlip="+transSlip);
		
		return maxID;
	}

}
