package scratch.bill.rsqsim;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;

public class ASCII_CatalogWriter {

	public static void main(String[] args) throws IOException {
		// directory that contains the .eList, .tList, .pList, .dList, geometry, and input files
		// you can find many catalogs here: /home/scec-00/rsqsim/catalogs
		// descriptions of theses catalogs are on Kevin's GitHub: https://github.com/kevinmilner/rsqsim-analysis/tree/master/catalogs/
		File catalogDir = new File("/home/kevin/Simulators/catalogs/rundir2585_1myr");
		
		// UCERF3 fault and deformation models.
		// currently we only use FM3.1, Geologic
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		
		File outputDir = catalogDir;
		
		double minMag = 0d;
		double skipYears = 5000; // skip the first 5000 years
		
		RSQSimCatalog catalog = new RSQSimCatalog(catalogDir, catalogDir.getName(), fm, dm);
		
		// get iterable events
		// you can also load the full event list as a list with .load(), but it requires more memory
		Iterable<RSQSimEvent> catalogIterable = catalog.loader().minMag(minMag).skipYears(skipYears).iterable();
		
		File outputFile = new File(outputDir, "ascii_hypocenter_catalog.txt.gz");
		
//		OutputStream out = new GZIPOutputStream(new FileOutputStream(outputFile));
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile))));
		
		GriddedRegion gridReg = new GriddedRegion(new CaliforniaRegions.RELM_TESTING(), 0.1, null);
		GriddedGeoDataSet nuclXYZ = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet particXYZ = new GriddedGeoDataSet(gridReg, false);
		
		RSQSimSubSectionMapper subsectionMapper = catalog.getSubSectMapper();
		
		out.write("# time-seconds magnitude hypo-lat hypo-lon hypo-depth");
		for (RSQSimEvent event : catalogIterable) {
			Location hypo = RSQSimUtils.getHypocenter(event);
			out.write(event.getTime()+"\t"+(float)event.getMagnitude()+"\t"+(float)hypo.getLatitude()+"\t"+(float)hypo.getLongitude()
				+"\t"+(float)hypo.getDepth()+"\n");
			
			// track nucleation with hypocenter
			int hypoIndex = gridReg.indexForLocation(hypo);
			if (hypoIndex >= 0)
				// will be -1 if outside of the relm region, possible for mendocino
				nuclXYZ.set(hypoIndex, nuclXYZ.get(hypoIndex)+1);
			
			// track participation over all elements
			ArrayList<SimulatorElement> elems = event.getAllElements();
			HashSet<Integer> nodes = new HashSet<>(); // unique set of grid node indexes which contain participating elements
			for (SimulatorElement elem : elems) {
				Location center = elem.getCenterLocation();
				int centerIndex = gridReg.indexForLocation(center);
				if (centerIndex >= 0)
					// will be -1 if outside of the relm region, possible for mendocino
					nodes.add(centerIndex);
			}
			for (int nodeIndex : nodes)
				particXYZ.set(nodeIndex, particXYZ.get(nodeIndex)+1);
			
			// subsection mappings
			// outer lits are parent sections, inner lists subsections
			
//			List<List<SubSectionMapping>> allSubSectionMappings = subsectionMapper.getFilteredSubSectionMappings(event, 0.2);
//			allSubSectionMappings.get(0).get(0).getSubSect().getSectionId(); // these are 0-based indexes
//			allSubSectionMappings.get(0).get(0).getSubSect().getParentSectionId(); // these are unique IDs
		}
		
		out.close();
		
		double duration = catalog.getDurationYears();
		// scale XYZ's from count to annual rate
		nuclXYZ.scale(1d/duration);
		particXYZ.scale(1d/duration);
		
		// write XYZs
		GriddedGeoDataSet.writeXYZFile(nuclXYZ, new File(outputDir, "gridded_nucleation_rates.xyz"));
		GriddedGeoDataSet.writeXYZFile(particXYZ, new File(outputDir, "gridded_participation_rates.xyz"));
	}

}
