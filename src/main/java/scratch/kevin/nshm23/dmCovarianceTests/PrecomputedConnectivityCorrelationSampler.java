package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.List;

import org.opensha.commons.util.Interpolate;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class PrecomputedConnectivityCorrelationSampler extends SectionCovarianceSampler {

	private SectionDistanceAzimuthCalculator distCalc;
	private double maxDist;
	private double zeroDistCoeff;
	private double negativeCorrMaxDist;
	
	private double[] totRates;
	private double[][] coruptureRates;

	public PrecomputedConnectivityCorrelationSampler(List<? extends FaultSection> subSects, FaultSystemSolution sol,
			SectionDistanceAzimuthCalculator distCalc,
			double maxDist, double zeroDistCoeff, double negativeCorrMaxDist) {
		super(subSects);
		this.distCalc = distCalc;
		this.maxDist = maxDist;
		this.zeroDistCoeff = zeroDistCoeff;
		this.negativeCorrMaxDist = negativeCorrMaxDist;
		
		calcCoruptureRates(sol);
	}
	
	private void calcCoruptureRates(FaultSystemSolution sol) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		int numSects = rupSet.getNumSections();
		int numRups = rupSet.getNumRuptures();
		totRates = new double[numSects];
		coruptureRates = new double[numSects][numSects];
		
		for (int r=0; r<numRups; r++) {
			double rate = sol.getRateForRup(r);
			if (rate == 0d)
				continue;
			List<Integer> rupSects = rupSet.getSectionsIndicesForRup(r);
			int numRupSects = rupSects.size();
			
			for (int i=0; i<numRupSects; i++) {
				int id1 = rupSects.get(i);
				totRates[id1] += rate;
				coruptureRates[id1][id1] += rate;
				for (int j=i+1; j<numRupSects; j++) {
					int id2 = rupSects.get(j);
					coruptureRates[id1][id2] += rate;
					coruptureRates[id2][id1] += rate;
				}
			}
		}
	}

	@Override
	public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
		double dist = distCalc.getDistance(sect1, sect2);
		if (dist > maxDist)
			return 0;
		int id1 = sect1.getSectionId();
		int id2 = sect2.getSectionId();
		double connRate = coruptureRates[id1][id2];
		if (connRate == 0d) {
			// unconnected
			if (dist < negativeCorrMaxDist) {
				// anticorrelated
				return Interpolate.findY(0d, -1d, negativeCorrMaxDist, 0d, dist);
			} else {
				// uncorrelated
				return 0d;
			}
		} else {
			// correlated, multiply connected fraction with distance fraction
			double avgRate = 0.5*(totRates[id1]+totRates[id2]);
			double connFract = connRate/avgRate;
			Preconditions.checkState(connFract >= 0 && connFract <= 1d);
			return connFract * Interpolate.findY(0d, zeroDistCoeff, maxDist, 0d, dist);
		}
	}

	@Override
	protected String getSamplerPrefix() {
		String prefix = "conn_corr_dist"+(float)maxDist+"km_zeroCoeff"+(float)zeroDistCoeff;
		if (negativeCorrMaxDist > 0)
			prefix += "_negCorrDist"+(float)negativeCorrMaxDist+"km";
		return prefix;
	}

}
