package scratch.kevin.simulators.ruptures.multifaultSeparate;

import java.io.File;
import java.io.IOException;
import java.util.BitSet;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class MultifaultTransFileDebug {

	public static void main(String[] args) throws IOException {
		File dir = new File("/data/kevin/simulators/catalogs/bruce/rundir5413_multifault_separate");
		RSQSimCatalog catalog = new RSQSimCatalog(dir, "Test Catalog", FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		RSQSimStateTransitionFileReader transReader = catalog.getTransitions();
		
		Iterator<RSQSimStateTime> transIt = transReader.getTransitionsIterable(0d, Double.MAX_VALUE).iterator();
		RSQSimStateTime curTrans = null;
		
		BitSet patchesSet = new BitSet(catalog.getElements().size()+1);
		
		RSQSimStateTime prevLastTrans = null;
		for (RSQSimEvent event : catalog.loader().load()) {
			int id = event.getID();
			double eventTime = event.getTime();
			
			System.out.println("Event "+id+" at "+eventTime+" s ("+event.getTimeInYears()+" yrs)");
			
			patchesSet.clear();
			for (int elemID : event.getAllElementIDs())
				patchesSet.set(elemID);
			
			if (prevLastTrans != null)
				Preconditions.checkState(prevLastTrans.absoluteTime < eventTime,
						"Previous last transition is after next event? nextID=%s, nextTime=%s\n\tprevLastTrans: %s",
						id, eventTime, prevLastTrans);
			
			if (curTrans == null)
				curTrans = transIt.next();
			
			Preconditions.checkState(curTrans.eventID == id,
					"First transition for %s has unexpected eventID:\n\t%s", id, curTrans);
			double firstTransTime = curTrans.absoluteTime;
			Preconditions.checkState(firstTransTime >= eventTime,
					"First transition for %s is before the event? eventTime=%s, delta=%s\n\ttrans: %s",
					id, eventTime, curTrans.absoluteTime - eventTime, curTrans);
			double prevTransTime = firstTransTime;
			Preconditions.checkState(patchesSet.get(curTrans.patchID), "Unexpected patchID=%s for event %s, trans: %s",
					curTrans.patchID, id, curTrans);
			
			int transCount = 1;
			prevLastTrans = curTrans;
			while (transIt.hasNext()) {
				curTrans = transIt.next();
				if (curTrans.eventID != id)
					break;
				transCount++;
				
				Preconditions.checkState(curTrans.absoluteTime >= prevTransTime,
						"Transitions out of order for %s: prevTransTime=%s, delta=%s;\n\ttrans: ",
						id, prevTransTime, curTrans.absoluteTime - prevTransTime, curTrans);
				prevTransTime = curTrans.absoluteTime;
				
				Preconditions.checkState(patchesSet.get(curTrans.patchID), "Unexpected patchID=%s for event %s\n\ttrans: %s",
						curTrans.patchID, id, curTrans);
				
				prevLastTrans = curTrans;
			}
			double lastTransTime = prevTransTime;
			System.out.println("\t"+transCount+" transitions for "+id+", duration: "+(float)(lastTransTime - firstTransTime)+" s");
		}
	}

}
