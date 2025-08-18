package scratch.kevin.prvi25.figures;

import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;

public class BunceDistanceCalc {

	public static void main(String[] args) throws IOException {
		GriddedRegion gridReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_MapExtents(), 0.01, GriddedRegion.ANCHOR_0_0);
		
		GriddedGeoDataSet mask = MapSourceTypeDisagg.buildLandMask(gridReg);
		
		List<? extends FaultSection> sects = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p2.getFaultSections();
		
		double minDist = Double.POSITIVE_INFINITY;
		String closestName = null;
		for (FaultSection sect : sects) {
			String name = sect.getName();
			if (name.contains("Bunce")) {
				System.out.println("Calculating distances to "+name);
				double myClosest = Double.POSITIVE_INFINITY;
				RuptureSurface surf = sect.getFaultSurface(1d);
				for (int i=0; i<mask.size(); i++) {
					if (mask.get(i) > 0) {
						double dist = surf.getDistanceJB(mask.getLocation(i));
						if (dist < myClosest)
							myClosest = dist;
					}
				}
				System.out.println("\tClosest distance is:\t"+(float)myClosest+" km");
				if (myClosest < minDist) {
					minDist = myClosest;
					closestName = name;
				}
			} else {
				System.out.println("(skipping "+name+")");
			}
		}
		
		System.out.println("Overall closest distance is:\t"+(float)minDist+" km");
		System.out.println("Closest section is:\t"+closestName);
	}

}
