package scratch.kevin.simulators.ruptures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimState;
import org.opensha.sha.simulators.srf.RSQSimStateTime;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SlipTimeDebug {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5597_CRUSTAL.instance();
		int eventID = 99110;
		
		RSQSimEvent event = catalog.loader().byID(eventID);
		
		System.out.println("Event "+eventID+" is an M"+(float)event.getMagnitude()+" at "+(float)event.getTimeInYears()+" yrs");
		
		RSQSimEventSlipTimeFunc origSlipTimeFunc = catalog.getSlipTimeFunc(event);
		RSQSimEventSlipTimeFunc slipTimeFunc = origSlipTimeFunc.asRelativeTimeFunc();
		
		List<Integer> patchIDs = new ArrayList<>(slipTimeFunc.getPatchIDs());
		patchIDs.sort(new Comparator<Integer>() {

			@Override
			public int compare(Integer p1, Integer p2) {
				return Double.compare(slipTimeFunc.getTimeOfFirstSlip(p1), slipTimeFunc.getTimeOfFirstSlip(p2));
			}
		});
		
		int maxSlipPatch = -1;
		double overallMaxSlip = 0d;
		boolean firstDebug = true;
		
		for (int p=0; p<patchIDs.size(); p++) {
			boolean print = p < 50 || patchIDs.size() <= 200 || p >= patchIDs.size()-52;
			if (p == 50 && patchIDs.size() > 200) {
				int newP = patchIDs.size()-52;
				System.out.println("(...skipping "+(newP-(p+1))+" patches...)");
				p = newP;
			}
			int patchID = patchIDs.get(p);
			double timeFirst = slipTimeFunc.getTimeOfFirstSlip(patchID);
			double timeLast = slipTimeFunc.getTimeOfLastSlip(patchID);
			double maxSlip = slipTimeFunc.getCumulativeEventSlip(patchID, timeLast);
			if (maxSlip > overallMaxSlip) {
				overallMaxSlip = maxSlip;
				maxSlipPatch = patchID;
			}
			if (print) {
				double duration = timeLast - timeFirst;
				List<RSQSimStateTime> states = slipTimeFunc.getTransitions(patchID);
				int numSlipStates = 0;
				for (RSQSimStateTime state : states)
					if (state.state == RSQSimState.EARTHQUAKE_SLIP)
						numSlipStates++;
				System.out.println("Patch "+patchID+" at "+(float)+timeFirst+"s => "+(float)timeLast+"s");
				System.out.println("\t"+(float)maxSlip+" m in "+(float)duration+" s");
				System.out.println("\t"+states.size()+" states, "+numSlipStates+" eq slip");
				if (maxSlip == 0d && numSlipStates > 0 && firstDebug) {
					System.out.println("\tDebugging zero slip with "+numSlipStates+" slip states");
					firstDebug = false;
					for (RSQSimStateTime state : states) {
						if (state.state == RSQSimState.EARTHQUAKE_SLIP) {
							System.out.println("\t\t"+state);
						}
					}
					DiscretizedFunc slipTime = slipTimeFunc.getSlipFunc(patchID);
					System.out.println(slipTime);
				}
			}
		}
		if (maxSlipPatch == -1) {
			System.out.println("No slip!");
		} else {
			System.out.println("Max slip was on "+maxSlipPatch+": "+(float)overallMaxSlip
					+" m starting at "+(float)slipTimeFunc.getTimeOfFirstSlip(maxSlipPatch)+" s");
		}
	}

}
