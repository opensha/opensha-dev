package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class CatInsidePolys {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		ArrayList<FaultSectionPrefData> sects = FaultModels.FM3_1.fetchFaultSections();
		
//		ObsEqkRupList cat = UCERF3_CatalogParser.loadCatalog(new File("/home/kevin/OpenSHA/UCERF3/UCERF3_Catalog3_0.txt"));
		ObsEqkRupList cat = UCERF3_CatalogParser.loadCatalog(new File("/home/kevin/OpenSHA/UCERF3/Felzer_UCERF3_Catalog4_0.txt"));
		
		int tot = 0;
		int contains = 0;
		for (ObsEqkRupture rup : cat) {
			Location loc = rup.getHypocenterLocation();
			tot++;
			
			for (FaultSectionPrefData sect : sects) {
				Region poly = sect.getZonePolygon();
				if (poly != null && poly.contains(loc)) {
					contains++;
					break;
				}
			}
		}
		
		System.out.println(contains+"/"+tot+" = "+(100d * (double)contains / (double)tot));
	}

}
