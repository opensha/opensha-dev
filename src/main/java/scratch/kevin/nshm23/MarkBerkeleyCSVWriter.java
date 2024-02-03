package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class MarkBerkeleyCSVWriter {

	public static void main(String[] args) throws IOException {
		Location loc = new Location(37.8707, -122.2508);
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		File outputDir = new File("/tmp/berkeley_participating_rups");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minDist = Double.POSITIVE_INFINITY;
		FaultSection closestSect = null;
		double traceConsiderDist = 50d;
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			boolean consider = false;
			for (Location traceLoc : sect.getFaultTrace()) {
				if (LocationUtils.horzDistanceFast(traceLoc, loc) < traceConsiderDist) {
					consider = true;
					break;
				}
			}
			if (consider) {
				for (Location traceLoc : sect.getFaultSurface(1d).getEvenlyDiscritizedUpperEdge()) {
					double dist = LocationUtils.horzDistanceFast(traceLoc, loc);
					if (dist < minDist) {
						minDist = dist;
						closestSect = sect;
					}
				}
			}
		}
		Preconditions.checkNotNull(closestSect);
		System.out.println("Closest section: "+closestSect.getSectionId()+". "
				+closestSect.getSectionName()+" ("+(float)minDist+" km away)");
		
		CSVFile<String> csv = new CSVFile<>(false);
		
		csv.addLine("Rupture Index", "Magnitude", "Annual Rate", "Average Rake",
				"Average Dip", "Num Sections", "Section Index 1", "Section Index N");
		
		for (int rupIndex : rupSet.getRupturesForSection(closestSect.getSectionId())) {
			List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
			List<String> line = new ArrayList<>(6+sects.size());
			
			line.add(rupIndex+"");
			line.add((float)rupSet.getMagForRup(rupIndex)+"");
			line.add(sol.getRateForRup(rupIndex)+"");
			line.add((float)rupSet.getAveRakeForRup(rupIndex)+"");
			line.add((float)rupSet.getSurfaceForRupture(rupIndex, 1d).getAveDip()+"");
			line.add(sects.size()+"");
			for (int sect : sects)
				line.add(sect+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, "participating_ruptures.csv"));
		GeoJSONFaultReader.writeFaultSections(new File(outputDir, "all_sections.geojson"), rupSet.getFaultSectionDataList());
	}

}
