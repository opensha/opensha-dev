package scratch.kevin.prvi25.figures;

import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class SouthLajasSlipIn03Calc {

	public static void main(String[] args) throws IOException {
		List<? extends FaultSection> sects = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1.getFaultSections();
		FaultSection sect = FaultSectionUtils.findSection(sects, "Lajas");

		double mag = 7d;
		double ri = 3500;
		
		double origDDW = sect.getReducedDownDipWidth()*1e3;
		double reducedDDW = sect.getReducedDownDipWidth()*1e3;
		double sumWeight = 0d;
		double sumWeightedSlipRate = 0d;
		for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values()) {
			double weight = scale.getNodeWeight(null);
			if (weight == 0d)
				continue;
			System.out.println(scale);
			double len = calcLenForMag(scale, origDDW, reducedDDW, mag);
			double calcMag = scale.getMag(len*reducedDDW, len, reducedDDW, origDDW, Double.NaN);
			System.out.println("\tCalculated len="+(float)(len*1e-3)+" km to get M="+(float)calcMag);
			double slip = scale.getAveSlip(len*reducedDDW, len, reducedDDW, origDDW, Double.NaN);
			System.out.println("\tCalculated slip="+(float)slip+" m");
			sumWeight += weight;
			double slipRate = 1e3*slip/ri;
			System.out.println("\tCalculated slip rate="+(float)slipRate+" mm/yr");
			sumWeightedSlipRate += slipRate*weight;
		}
		System.out.println("Average slip rate: "+(float)(sumWeightedSlipRate/sumWeight)+" mm/yr");
	}
	
	private static double calcLenForMag(NSHM23_ScalingRelationships scale, double origDDW, double reducedDDW, double mag) {
		double minLen = 0d;
		double maxLen = 1000*1e3; // in m
		int numLen = 100;
		double bestLen = Double.NaN;
		for (int i=0; i<5; i++) {
			EvenlyDiscretizedFunc lenDist = new EvenlyDiscretizedFunc(minLen, maxLen, numLen);
			int lastBelow = -1;
			int firstAbove = -1;
			for (int l=0; l<numLen; l++) {
				double len = lenDist.getX(l);
				double calcMag = scale.getMag(len*reducedDDW, len, reducedDDW, origDDW, Double.NaN);
				if (calcMag <= mag)
					lastBelow = l;
				else if (firstAbove < 0)
					firstAbove = l;
				lenDist.set(l, calcMag);
			}
			Preconditions.checkState(lastBelow >= 0);
			Preconditions.checkState(firstAbove > lastBelow);
			bestLen = lenDist.getFirstInterpolatedX(mag);
			minLen = lenDist.getX(lastBelow);
			maxLen = lenDist.getX(firstAbove);
		}
		return bestLen;
	}

}
