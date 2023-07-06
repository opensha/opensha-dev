package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.CubedGriddedRegion;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;

import com.google.common.base.Preconditions;

public class GridAndCubeLocWriter {
	
	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File(invDir,
				"results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		FaultCubeAssociations cubeAssoc = rupSet.requireModule(FaultCubeAssociations.class);
		
		File outputDir = new File(invDir, "associations");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>();
		archive.addModule(cubeAssoc);
		archive.write(new File(outputDir, "fault_cube_associations.zip"));
		
		CSVFile<String> nodeCSV = new CSVFile<>(true);
		nodeCSV.addLine("Grid Node Index", "Latitude", "Longitude");
		GriddedRegion reg = cubeAssoc.getRegion();
		for (int i=0; i<reg.getNodeCount(); i++) {
			Location loc = reg.getLocation(i);
			nodeCSV.addLine(i+"", (float)loc.lat+"", (float)loc.lon+"");
		}
		nodeCSV.writeToFile(new File(outputDir, "grid_node_locations.csv"));
		
		CSVFile<String> cubeCSV = new CSVFile<>(true);
		cubeCSV.addLine("Cube Index", "Latitude", "Longitude", "Depth (km)");
		CubedGriddedRegion cgr = cubeAssoc.getCubedGriddedRegion();
		for (int i=0; i<cgr.getNumCubes(); i++) {
			int[] sects = cubeAssoc.getSectsAtCube(i);
			if (sects != null && sects.length > 0) {
				Location loc = cgr.getCubeLocationForIndex(i);
				cubeCSV.addLine(i+"", (float)loc.lat+"", (float)loc.lon+"", (float)loc.depth+"");
			}
		}
		cubeCSV.writeToFile(new File(outputDir, "cube_locations.csv"));
		
		// also write subsections
		GeoJSONFaultReader.writeFaultSections(new File(outputDir, "sub_sections.geojson"), rupSet.getFaultSectionDataList());
	}

}
