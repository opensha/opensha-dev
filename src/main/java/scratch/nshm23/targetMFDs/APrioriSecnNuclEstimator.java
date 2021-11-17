package scratch.nshm23.targetMFDs;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class APrioriSecnNuclEstimator implements DataSectNucleationRateEstimator {
	
	private FaultSystemRupSet rupSet;
	private List<Integer> rupIndexes;
	private double[] rates;
	
	private HashSet<Integer> sectIDs;

	public APrioriSecnNuclEstimator(FaultSystemRupSet rupSet, List<Integer> rupIndexes, double totRate) {
		this(rupSet, rupIndexes, distributeEqually(rupIndexes.size(), totRate));
	}
	
	private static double[] distributeEqually(int numRates, double totRate) {
		double[] rates = new double[numRates];
		Arrays.fill(rates, totRate/(double)numRates);
		return rates;
	}
	
	public APrioriSecnNuclEstimator(FaultSystemRupSet rupSet, List<Integer> rupIndexes, double[] rates) {
		this.rupSet = rupSet;
		this.rupIndexes = rupIndexes;
		this.rates = rates;
		Preconditions.checkState(rates.length > 0);
		Preconditions.checkState(rates.length == rupIndexes.size());
		
		sectIDs = new HashSet<>();
		for (int rupIndex : rupIndexes)
			sectIDs.addAll(rupSet.getSectionsIndicesForRup(rupIndex));
	}

	@Override
	public boolean appliesTo(FaultSection sect) {
		return sectIDs.contains(sect.getSectionId());
	}

	@Override
	public double estimateNuclRate(FaultSection sect, IncrementalMagFreqDist curSectSupraSeisMFD) {
		Preconditions.checkState(appliesTo(sect));
		// assume that the nucleation rate is only from the a-priori ruptures, good enough for parkfield at least
		
		double sectArea = rupSet.getAreaForSection(sect.getSectionId());
		double nuclRate = 0d;
		for (int r=0; r<rates.length; r++) {
			int rupIndex = rupIndexes.get(r);
			if (rupSet.getSectionsIndicesForRup(rupIndex).contains(sect.getSectionId()))
				nuclRate += rates[r]*sectArea/rupSet.getAreaForRup(rupIndex);
		}
		return nuclRate;
	}
	
	
}