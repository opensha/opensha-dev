package scratch.kevin.segmentation;

import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator.AggregationMethod;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchAlignment;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;
import org.opensha.sha.util.FocalMech;

public class RelativeCFFTests {

	public static void main(String[] args) {
		Location origin = new Location(0d, 0d);
		double len1 = 10d;
//		FocalMech mech1 = FocalMech.STRIKE_SLIP;
		FocalMech mech1 = FocalMech.REVERSE;
		double strike1 = 0d;
		
//		FocalMech mech2 = FocalMech.STRIKE_SLIP;
		FocalMech mech2 = FocalMech.REVERSE;
		double jumpAngle = 0d;
//		double strike2 = 0d;
		double strike2 = 0d;
		double distance = 1;
//		double len2 = len1;
		double len2 = 10d;
		
		double upper = 0d;
		double lower = 10d;
		
		double stiffGridSpacing = 2d;
		double coeffOfFriction = 0.5;
//		double coeffOfFriction = 0.0;
		double lameLambda = 3e4;
		double lameMu = 3e4;
		
		Location end1 = LocationUtils.location(origin, Math.toRadians(strike1), len1);
		
		GeoJSONFaultSection sect1 = new GeoJSONFaultSection.Builder(1, "Source Fault",
				FaultTrace.of(origin, end1))
				.upperDepth(upper).lowerDepth(lower)
				.dip(mech1.dip()).rake(mech1.rake()).build();
		
		Location start2 = LocationUtils.location(end1, Math.toRadians(jumpAngle), distance);
		Location end2 = LocationUtils.location(start2, Math.toRadians(strike2), len2);
		GeoJSONFaultSection sect2 = new GeoJSONFaultSection.Builder(2, "Destination Fault",
				FaultTrace.of(start2, end2))
				.upperDepth(upper).lowerDepth(lower)
				.dip(mech2.dip()).rake(mech2.rake()).build();
		
		GeoJSONFaultSection sect2_equivRef = new GeoJSONFaultSection.Builder(1002, "Reference Destination Fault",
				FaultTrace.of(end1, LocationUtils.location(end1, Math.toRadians(strike1), len2)))
				.upperDepth(upper).lowerDepth(lower)
				.dip(mech1.dip()).rake(mech1.rake()).build();
		
		List<GeoJSONFaultSection> sects = List.of(sect1, sect2, sect2_equivRef);
		
//		SubSectStiffnessCalculator calc = new SubSectStiffnessCalculator(null, distance, len2, upper, lower)
		SubSectStiffnessCalculator stiffnessCalc = new SubSectStiffnessCalculator(
				sects, stiffGridSpacing, lameLambda, lameMu, coeffOfFriction, PatchAlignment.FILL_OVERLAP, 1d);
		AggregatedStiffnessCalculator sumAgg = new AggregatedStiffnessCalculator(StiffnessType.CFF, stiffnessCalc, true,
				AggregationMethod.FLATTEN, AggregationMethod.SUM, AggregationMethod.SUM, AggregationMethod.SUM);
		
		Shaw07JumpDistProb[] segModels = {
				NSHM23_SegmentationModels.HIGH.getShawModel(),
				NSHM23_SegmentationModels.MID.getShawModel(),
				NSHM23_SegmentationModels.LOW.getShawModel()
		};
		Shaw07JumpDistProb vanillaSeg = new Shaw07JumpDistProb(1d, 3d);
		
		System.out.println("Distance: "+(float)distance+" km");
		
		double refStiffness = sumAgg.calc(sect1, sect2_equivRef);
		System.out.println("Ref stiffness: "+(float)refStiffness);
		double sect2Stiffness = sumAgg.calc(sect1, sect2);
		System.out.println("Actual stiffness: "+(float)sect2Stiffness);
		System.out.println("Actual relative: "+(float)(sect2Stiffness/refStiffness));
		
		System.out.print("NSHM23 segs:");
		for (Shaw07JumpDistProb seg : segModels)
			System.out.print(" "+(float)seg.calcJumpProbability(distance));
		System.out.println();
		System.out.println("Vanilla ("+vanillaSeg+"): "+(float)vanillaSeg.calcJumpProbability(distance));
	}

}
