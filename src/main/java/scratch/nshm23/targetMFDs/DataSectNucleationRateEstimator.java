package scratch.nshm23.targetMFDs;

import java.util.Arrays;
import java.util.List;

import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public abstract class DataSectNucleationRateEstimator {
	public abstract boolean appliesTo(FaultSection sect);
	
	public abstract IncrementalMagFreqDist estimateNuclMFD(FaultSection sect, IncrementalMagFreqDist curSectSupraSeisMFD,
			List<Integer> availableRupIndexes, List<Double> availableRupMags, UncertainDataConstraint sectMomentRate,
			boolean sparseGR);
	
	protected IncrementalMagFreqDist buildGRFromRate(IncrementalMagFreqDist refFunc, List<Double> rupMags,
			double rate, double moRate, boolean sparseGR) {
		double minMag = Double.POSITIVE_INFINITY;
		for (double mag : rupMags)
			minMag = Math.min(minMag, mag);
		int minIndex = refFunc.getClosestXIndex(minMag);
		minMag = refFunc.getX(minIndex);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
				refFunc.getMinX(), refFunc.size(), refFunc.getDelta());
		gr.setAllButBvalue(minMag, refFunc.getMaxX(), moRate, rate);
		
		if (sparseGR)
			return SparseGutenbergRichterSolver.getEquivGR(refFunc, rupMags, moRate, gr.get_bValue());
		return gr;
	}
	
	protected IncrementalMagFreqDist buildGRFromBVal(IncrementalMagFreqDist refFunc, List<Double> rupMags,
			double bVal, double moRate, boolean sparseGR) {
		double minMag = Double.POSITIVE_INFINITY;
		for (double mag : rupMags)
			minMag = Math.min(minMag, mag);
		int minIndex = refFunc.getClosestXIndex(minMag);
		minMag = refFunc.getX(minIndex);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
				refFunc.getMinX(), refFunc.size(), refFunc.getDelta());
		gr.setAllButTotCumRate(minMag, gr.getMaxX(), moRate, bVal);
		
		if (sparseGR)
			return SparseGutenbergRichterSolver.getEquivGR(refFunc, rupMags, moRate, gr.get_bValue());
		return gr;
	}
	
	protected UncertainIncrMagFreqDist getBounded(UncertaintyBoundType boundType, IncrementalMagFreqDist... mfds) {
		IncrementalMagFreqDist middle = new IncrementalMagFreqDist(mfds[0].getMinX(), mfds[0].size(), mfds[0].getDelta());
		IncrementalMagFreqDist lower = new IncrementalMagFreqDist(middle.getMinX(), middle.size(), middle.getDelta());
		IncrementalMagFreqDist upper = new IncrementalMagFreqDist(middle.getMinX(), middle.size(), middle.getDelta());
		for (int i=0; i<middle.size(); i++) {
			// these bounds aren't actually upper/lower, one will be lower at low M and highter at high M
			double[] vals = new double[mfds.length];
			for (int j=0; j<vals.length; j++)
				vals[j] = mfds[j].getY(i);
			Arrays.sort(vals);
			middle.set(i, vals.length == 3 ? vals[1] : DataUtils.median_sorted(vals));
			lower.set(i, vals[0]);
			upper.set(i, vals[2]);
		}
		return new UncertainBoundedIncrMagFreqDist(middle, lower, upper, boundType);
	}

}
