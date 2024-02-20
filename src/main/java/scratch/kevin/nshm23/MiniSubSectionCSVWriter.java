package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class MiniSubSectionCSVWriter {

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v1p4;
		
		List<? extends FaultSection> fullSects = fm.getFaultSections();
		
		List<? extends FaultSection> subSects = fm.getDefaultDeformationModel().build(fm);
		Map<Integer, List<FaultSection>> subSectIDmap = subSects.stream().collect(
				Collectors.groupingBy(s -> s.getParentSectionId()));
		
		CSVFile<String> csv = new CSVFile<>(false);
		
		csv.addLine("Section ID", "Section Name", "Minisection Index", "Start Lat", "Start Lon",
				"End Lat", "End Lon", "Subsection Index(es)");
		
		MinMaxAveTracker sectMappingStats = new MinMaxAveTracker();
		
		for (FaultSection sect : fullSects) {
			FaultTrace trace = sect.getFaultTrace();
			
			List<FaultSection> mySubSects = subSectIDmap.get(sect.getSectionId());
			
			for (int i=1; i<trace.size(); i++) {
				Location start = trace.get(i-1);
				Location end = trace.get(i);
				
				List<String> line = new ArrayList<>();
				line.add(sect.getSectionId()+"");
				line.add(sect.getSectionName());
				line.add((i-1)+""); // minsection index, 0-based
				line.add(start.getLatitude()+"");
				line.add(start.getLongitude()+"");
				line.add(end.getLatitude()+"");
				line.add(end.getLongitude()+"");
				int numFound = 0;
				for (FaultSection subSect : mySubSects) {
					for (Location loc : subSect.getFaultTrace()) {
						if (loc.equals(start) || loc.equals(end) || LocationUtils.areSimilar(start, loc)
								|| LocationUtils.areSimilar(end, loc)) {
							numFound++;
							line.add(subSect.getSectionId()+"");
							break;
						}
					}
				}
				Preconditions.checkState(numFound > 0);
				sectMappingStats.addValue(numFound);
				csv.addLine(line);
			}
		}
		System.out.println("Sub-sect mapping stats: "+sectMappingStats);
		csv.writeToFile(new File("/tmp/"+fm.getFilePrefix()+"_mini_and_sub_sects.csv"));
	}

}
