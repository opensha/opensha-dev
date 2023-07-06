package scratch.kevin.simulators.ruptures.subduction;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileWriter;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class BruceSubductionFilter {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5566.instance();
		File catalogDir = catalog.getCatalogDir();
		
		HashSet<Integer> subductionPatchIDs = new HashSet<>();
		HashSet<Integer> crustalPatchIDs = new HashSet<>();
		
		for (SimulatorElement elem : catalog.getElements()) {
			if (elem.getSectionName().trim().startsWith("nn1"))
				subductionPatchIDs.add(elem.getID());
			else
				crustalPatchIDs.add(elem.getID());
		}
		System.out.println("Found "+subductionPatchIDs.size()+" subduction elements");
		System.out.println("Found "+crustalPatchIDs.size()+" crustal elements");
		
		boolean includeOtherCorupturing = false;
		
//		RuptureIdentifier coruptureIden = null;
		RuptureIdentifier coruptureIden = new LogicalAndRupIden(
				new ElementIden("Subduction", new ArrayList<>(subductionPatchIDs)),
				new ElementIden("Crustal", new ArrayList<>(crustalPatchIDs)));
		
		List<RSQSimEvent> events = catalog.loader().magRange(6.5, 11d).load();
		System.out.println("Loaded "+events.size()+" events");
		
		if (coruptureIden != null) {
			int numMatches = 0;
			int numNuclCrustal = 0;
			int numNuclSubduction = 0;
			for (RSQSimEvent event : catalog.loader().minMag(6.5).matches(coruptureIden).iterable()) {
				SimulatorElement hypoEl = RSQSimUtils.getHypocenterElem(event);
				boolean subHypo = subductionPatchIDs.contains(hypoEl.getID());
				if (subHypo)
					numNuclSubduction++;
				else
					numNuclCrustal++;
				numMatches++;
			}
			DecimalFormat pDF = new DecimalFormat("0.00%");
			System.out.println(numNuclSubduction+"/"+numMatches+" ("
					+pDF.format((double)numNuclSubduction/(double)numMatches)+") nucleate on subduction");
			System.out.println(numNuclCrustal+"/"+numMatches+" ("
					+pDF.format((double)numNuclCrustal/(double)numMatches)+") nucleate on crustal");
		}
		System.exit(0);
		
		RSQSimStateTransitionFileReader trans = catalog.getTransitions();
		
		File paramFile = new File(catalogDir, "multiparam.in");
		File elemFile = new File(catalogDir, "zfault_Deepen.in");
		
		for (boolean subduction : new boolean[] {false,true}) {
			String prefix;
			HashSet<Integer> includePatches;
			if (subduction) {
				prefix = "subduction";
				includePatches = subductionPatchIDs;
			} else {
				prefix = "crustal";
				includePatches = crustalPatchIDs;
			}
			if (coruptureIden != null)
				prefix += "_corupture";
			System.out.println("Writing "+prefix);
			String dirName = catalogDir.getName()+"_"+prefix;
			File outputDir = new File(catalogDir.getParentFile(), dirName);
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			
			RSQSimFileWriter.patchFilterCatalog(events, catalog.getElements(), includePatches, includeOtherCorupturing,
					coruptureIden, outputDir, prefix, false, trans);
			
			Files.copy(paramFile, new File(outputDir, paramFile.getName()));
			Files.copy(elemFile, new File(outputDir, elemFile.getName()));
		}
	}

}
