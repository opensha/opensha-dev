package scratch.kevin.simulators;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimUtils;

public class RSQSimSubSectMappingTest {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_3014.instance(baseDir);
//		RSQSimCatalog catalog = RSQSimCatalog.Catalogs.BRUCE_2585.instance(baseDir);
		
		File catalogDir = catalog.getCatalogDir();
		File outputFile = new File(catalogDir, "sub_sect_mappings.txt");
		FileWriter fw = new FileWriter(outputFile);
		
		double spacing = 0.1d;
		
		List<FaultSectionPrefData> subSects = catalog.getU3SubSects();
		List<SimulatorElement> elems = catalog.getElements();
		int offset = RSQSimUtils.getSubSectIndexOffset(elems, subSects);
		fw.write("Catalog: "+catalog.getName()+"\n");
		if (offset != 0)
			fw.write("Using sub-section index offset: "+offset+"\n");
		fw.write("All distances are horizontal distances in kilometers to the UCERF3 surface (rJB), which is discretized with "
				+(float)spacing+"km spacing\n");
		fw.write("\n");
		
		DecimalFormat df = new DecimalFormat("0.000");
		int numMissing = 0;
		MinMaxAveTracker totTrack = new MinMaxAveTracker();
		for (int s=0; s<subSects.size(); s++) {
			MinMaxAveTracker sectTrack = new MinMaxAveTracker();
			FaultSectionPrefData subSect = subSects.get(s);
			StirlingGriddedSurface surf = subSect.getStirlingGriddedSurface(spacing, false, false);
			for (SimulatorElement elem : elems) {
				int mappedIndex = elem.getSectionID()+offset;
				if (mappedIndex != s)
					continue;
				Location loc = elem.getCenterLocation();
				double dist = surf.getDistanceJB(new Location(loc.getLatitude(), loc.getLongitude()));
				sectTrack.addValue(dist);
				totTrack.addValue(dist);
			}
			System.out.println(s+". "+subSect.getName());
			fw.write(s+". "+subSect.getName()+"\n");
			if (sectTrack.getNum() == 0) {
				numMissing++;
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
		fw.write(numMissing+" sub-sections without mappings\n");
		fw.close();
	}

}
