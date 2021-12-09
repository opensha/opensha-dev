package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.faultSysSolution.modules.PolygonFaultGridAssociations;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;

public class RegionInsidePolysWriter {

	public static void main(String[] args) throws IOException {
		PolygonFaultGridAssociations associations = FaultPolyMgr.loadSerializedUCERF3(FaultModels.FM3_1);
		GriddedRegion reg = associations.getRegion();
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Latitude", "Longitude", "Fraction Associated w/ Faults");
		for (int i=0; i<reg.getNodeCount(); i++) {
			double fract = associations.getNodeFraction(i);
			Location loc = reg.getLocation(i);
			csv.addLine((float)loc.getLatitude()+"", (float)loc.getLongitude()+"", (float)fract+"");
		}
		
		csv.writeToFile(new File("/tmp/grid_reg_fault_associations.csv"));
	}

}
