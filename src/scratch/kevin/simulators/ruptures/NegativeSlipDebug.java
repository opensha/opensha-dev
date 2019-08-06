package scratch.kevin.simulators.ruptures;

import java.io.IOException;

import org.opensha.sha.simulators.RSQSimEvent;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class NegativeSlipDebug {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance();
		
		for (RSQSimEvent event : catalog.loader().iterable()) {
			int[] ids = event.getAllElementIDs();
			double[] slips = event.getAllElementSlips();
			for (int i=0; i<ids.length; i++) {
				if (slips[i] < 0) {
					System.out.println("NEGATIVE SLIP:\teventID="+event.getID()+"\tM="+(float)event.getMagnitude()
						+"\tpatchID="+ids[i]+"\tslip="+slips[i]);
				}
			}
		}
	}

}
