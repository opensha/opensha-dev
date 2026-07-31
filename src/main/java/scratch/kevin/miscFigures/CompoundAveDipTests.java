package scratch.kevin.miscFigures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupCartoonGenerator;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;

public class CompoundAveDipTests {

	public static void main(String[] args) throws IOException {
		List<GeoJSONFaultSection> sects = new ArrayList<>();
		HashMap<Integer, Boolean> flips = new HashMap<>();
		sects.add(new GeoJSONFaultSection.Builder(sects.size(), "Fault 1",
				FaultTrace.of(new Location(0d, 0.48d), new Location(0d, 0d)))
				.dip(60).lowerDepth(10).upperDepth(0).rake(0)
				.build());
		flips.put(sects.get(sects.size()-1).getSectionId(), true);
		sects.add(new GeoJSONFaultSection.Builder(sects.size(), "Fault 2",
				FaultTrace.of(new Location(0d, 0.52d), new Location(0d, 1d)))
				.dip(45).lowerDepth(10).upperDepth(0).rake(0)
				.build());
		
		List<RuptureSurface> surfs = new ArrayList<>();
		for (GeoJSONFaultSection sect : sects)
			surfs.add(sect.getFaultSurface(1d));
		
		CompoundSurface.Simple surf = new CompoundSurface.Simple(surfs, sects);
		
		ClusterRupture cRup = null;
		int subSectIndex = 0;
		for (GeoJSONFaultSection sect : sects) {
			double ddw = sect.getOrigDownDipWidth();
			List<? extends FaultSection> subSects = sect.getSubSectionsList(ddw/2d, subSectIndex);
			if (flips.containsKey(sect.getSectionId())) {
				subSects = new ArrayList<>(subSects);
				Collections.reverse(subSects);
			}
			subSectIndex += subSects.size();
			if (cRup == null) {
				cRup = new ClusterRupture(new FaultSubsectionCluster(subSects));
			} else {
				FaultSubsectionCluster prev = cRup.clusters[cRup.clusters.length-1];
				FaultSection from = null;
				FaultSection to = null;
				double minDist = Double.POSITIVE_INFINITY;
				for (FaultSection testFrom : prev.subSects) {
					for (FaultSection testTo : subSects) {
						double dist = Double.POSITIVE_INFINITY;
						for (Location l1 : testFrom.getFaultTrace())
							for (Location l2 : testTo.getFaultTrace())
								dist = Math.min(dist, LocationUtils.horzDistanceFast(l2, l1));
						if (dist < minDist) {
							minDist = dist;
							from = testFrom;
							to = testTo;
						}
					}
				}
				cRup = cRup.take(new Jump(from, prev, to, new FaultSubsectionCluster(subSects), 0));
			}
		}
		
		DecimalFormat dipDF = new DecimalFormat("0");
		String title = "Oriented dip: "+dipDF.format(surf.getAveOrientedDip())+"; Average dip: "+dipDF.format(surf.getAveDip());
		
		RupCartoonGenerator.plotRupture(new File("/tmp"), "compound_dip_example", cRup, title, false, false);
		
		
	}

}
