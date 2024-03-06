package scratch.kevin.nshm23;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;

public class FullParentSectMagCalc {

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		NSHM23_ScalingRelationships[] scales = {
				NSHM23_ScalingRelationships.AVERAGE,
				NSHM23_ScalingRelationships.LOGA_C4p3,
				NSHM23_ScalingRelationships.LOGA_C4p2,
				NSHM23_ScalingRelationships.LOGA_C4p1,
				NSHM23_ScalingRelationships.WIDTH_LIMITED,
//				NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN,
//				NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD
		};
		
		List<? extends FaultSection> sects = dm.build(fm);
		
		int[] parentIDs = {
//				FaultSectionUtils.findParentSectionID(sects, "Compton")
//				FaultSectionUtils.findParentSectionID(sects, "Raymond")
				FaultSectionUtils.findParentSectionID(sects, "Puente Hills", "Los Angeles"),
				FaultSectionUtils.findParentSectionID(sects, "Puente Hills", "Santa Fe Springs"),
				FaultSectionUtils.findParentSectionID(sects, "Puente Hills", "Coyote Hills"),
		};
		
		List<FaultSection> matchingSects = new ArrayList<>();
		
		for (FaultSection sect : sects)
			for (int parentID : parentIDs)
				if (sect.getParentSectionId() == parentID)
					matchingSects.add(sect);
		
		System.out.println("Found "+matchingSects.size()+" matching sects");
		for (FaultSection sect : matchingSects)
			System.out.println("\t"+sect.getSectionName());
		
		double totArea = 0;
		double totOrigArea = 0;
		double totLen = 0;
		List<Double> rakes = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		for (FaultSection sect : matchingSects) {
			totOrigArea += sect.getArea(false);	// sq-m
			totArea += sect.getArea(true);	// sq-m
			totLen += sect.getTraceLength()*1e3; // km -> m
			rakes.add(sect.getAveRake());
			weights.add(sect.getArea(true));
		}
		double origDDW = totOrigArea / totLen;
		double aveRake = FaultUtils.getInRakeRange(FaultUtils.getScaledAngleAverage(weights, rakes));
		
		for (NSHM23_ScalingRelationships scale : scales) {
			double mag = scale.getMag(totArea, totLen, totArea/totLen, origDDW, aveRake);
			System.out.println(scale.getShortName()+":\t"+(float)mag);
		}
	}

}
