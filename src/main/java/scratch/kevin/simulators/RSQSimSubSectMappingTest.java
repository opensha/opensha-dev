package scratch.kevin.simulators;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;

public class RSQSimSubSectMappingTest {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
//		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_3164.instance(baseDir);
		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File catalogDir = catalog.getCatalogDir();
		File outputFile = new File(catalogDir, "sub_sect_mappings.txt");
		FileWriter fw = new FileWriter(outputFile);
		
		double spacing = 0.1d;
		
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		
		List<? extends FaultSection> subSects = catalog.getU3SubSects();
		List<SimulatorElement> elems = catalog.getElements();
		int offset = RSQSimUtils.getSubSectIndexOffset(elems, subSects);
		fw.write("Catalog: "+catalog.getName()+"\n");
		if (offset != 0)
			fw.write("Using sub-section index offset: "+offset+"\n");
		fw.write("All distances are horizontal distances in kilometers to the UCERF3 surface (rJB), which is discretized with "
				+(float)spacing+"km spacing\n");
		fw.write("\n");
		
		DecimalFormat df = new DecimalFormat("0.000");
		List<FaultSection> missingSubSects = new ArrayList<>();
		MinMaxAveTracker totTrack = new MinMaxAveTracker();
		for (int s=0; s<subSects.size(); s++) {
			MinMaxAveTracker sectTrack = new MinMaxAveTracker();
			FaultSection subSect = subSects.get(s);
			RuptureSurface surf = subSect.getFaultSurface(spacing, false, false);
			for (SimulatorElement elem : mapper.getElementsForSection(subSect)) {
				Location loc = elem.getCenterLocation();
				double dist = surf.getDistanceJB(new Location(loc.getLatitude(), loc.getLongitude()));
				sectTrack.addValue(dist);
				totTrack.addValue(dist);
			}
			System.out.println(s+". "+subSect.getName());
			fw.write(s+". "+subSect.getName()+"\n");
			if (sectTrack.getNum() == 0) {
				missingSubSects.add(subSect);
				fw.write("\tno mapped elements\n");
			} else {
				fw.write("\t"+sectTrack.getNum()+" mapped elements\n");
				fw.write("\tDistances:\tmin="+df.format(sectTrack.getMin())+"\tmax="+df.format(sectTrack.getMax())
					+"\tmean="+df.format(sectTrack.getAverage())+"\n");
			}
			System.out.println("\t"+sectTrack);
		}
		fw.write("Total Distances:\tmin="+df.format(totTrack.getMin())+"\tmax="+df.format(totTrack.getMax())
			+"\tmean="+df.format(totTrack.getAverage())+"\n");
		fw.write(missingSubSects.size()+" sub-sections without mappings"+(missingSubSects.isEmpty() ? "" : ":")+"\n");
		for (FaultSection missing : missingSubSects)
			fw.write("\t"+missing.getSectionId()+". "+missing.getName()+"\n");
		fw.close();
	}

}
