package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class CatalogLoadDemo {

	public static void main(String[] args) throws IOException {
		// input catalog dir. At a minimum, this directory must contain the following 6 files
		// 		<prefix>.dList
		// 		<prefix>.eList
		// 		<prefix>.tList
		// 		<prefix>.pList
		// 		<parameter file>.in
		// 		<fault file>.flt or <fault fault>.in
		// 
		// The prefix for the 4 list files can be empty (Bruce usually leaves it empty), so those files
		// could be called simply .pList, .eList, etc..., which can be annoying, especially if your
		// file browser or ls alias doesn't show hidden files. You'll see them with ls -la though.
		// 
		// There is no convention for the name of the parameter file, except that is usually ends in ".in,"
		// though there can be multiple ".in" files and you'll need to make sure that the correct one is there.
		// Bruce names his parameter files "multiparam.in" while Jacqui gives it a name similar to <prefix>.
		//
		// Similar for the fault input file. It's name will be listed in the 'faultFname' field in the parameter
		// file. Buce usually uses "zfault_Deepen.in" and Jacqui something ending in .flt
		File catalogDir = new File("/data/kevin/simulators/catalogs/rundir2585_1myr");
		
		RSQSimCatalog catalog = new RSQSimCatalog(catalogDir, "2585 Million Year  Extension",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		// so far all have been FM 3.1/GEOL, but this could change in the future. used for mappings
		
		// now load the events using catalog.loader()
		Loader loader = catalog.loader();
		
		// you'll almost always want to filter them as they are loaded to reduce memory requirements.
		// set this as high as you can to do what you need
		loader.minMag(6.5);
		
		// we usually skip the first 5k years to avoid spin up time
		loader.skipYears(5000);
		
		// now load them inloader
		List<RSQSimEvent> events = loader.load();
		System.out.println("Loaded "+events.size()+" events!");
		
		// you could have done this in a single line:
		// 
		// List<RSQSimEvent> events = catalog.loader().minMag(6.5).skipYears(5000).load()
		//
		// or, to process the events as they are loaded (without needing to be able to store the
		// entire catalog in memory at once:
		// 
		// for (RSQSimEvent event : loader.iterable()) { // do stuff }
		
		// RSQSimEvent extends SimulatorEvent, which I think you wrote, so it should be familiar
		
		// Here's how to get the UCERF3 mapped rupture for a given event:
		RSQSimSubSectEqkRupture rup = catalog.getMappedSubSectRupture(events.get(0), 0.2);
		// this extends EqkRupture with some additions, including getSubSections() to get the list
		// of mapped subsections. The '0.2' argument was the fraction of a subsection which is required
		// to rupture (by area) before calling it a match. If no such subsections match this criterion,
		// it currently falls back to matching all sections without the criterion, that could be changed though.
	}

}
