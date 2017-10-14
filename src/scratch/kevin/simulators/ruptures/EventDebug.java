package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTime;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class EventDebug {
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2273.instance(baseDir);
		int eventID = 816469;
		
		RSQSimEvent event = catalog.loader().byID(eventID);
		
		System.out.println("Event "+eventID+", M="+event.getMagnitude()
			+", time="+event.getTime()+"s = "+event.getTimeInYears()+" yr");
		ArrayList<SimulatorElement> elements = event.getAllElements();
		System.out.println(elements.size()+" elements");
		List<String> sectNames = new ArrayList<>();
		for (SimulatorElement elem : elements)
			if (!sectNames.contains(elem.getSectionName()))
				sectNames.add(elem.getSectionName());
		System.out.println("Sections:");
		for (String sectName : sectNames)
			System.out.println("\t"+sectName);
		
		// load transitions
		System.out.println("Catalog first trans: "+catalog.getTransitions().getFirstTransitionTime());
		System.out.println("Catalog last trans: "+catalog.getTransitions().getLastTransitionTime());
		Map<Integer, List<RSQSimStateTime>> trans = catalog.getTransitions().getTransitions(event);
		System.out.println("Transitions: "+trans.size());
		for (int id : trans.keySet()) {
			if (trans.get(id).isEmpty())
				continue;
			System.out.println("\tElement "+id);
			for (RSQSimStateTime t : trans.get(id))
				System.out.println("\t\t"+t);
		}
		RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
		System.out.println("Func start: "+func.getStartTime());
		System.out.println("Func end: "+func.getStartTime());
	}

}
