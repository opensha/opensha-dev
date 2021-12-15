package scratch.kevin.nshm23;

import java.text.DecimalFormat;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

/**
 * This constraint can force systemwide jump distances of a particular functional form
 * 
 * Probably not the way to do, as Shaw & Dieterich doesn't actually give us a region-wide target model
 * 
 * @author kevin
 *
 */
public class RelativeRupJumpDistConstraint extends InversionConstraint {
	
	private static final DecimalFormat oDF = new DecimalFormat("0.#");
	
	private transient FaultSystemRupSet rupSet;
	private transient EvenlyDiscretizedFunc distBinning;
	private transient int[] ruptureBins;
	private transient double[] ruptureCondProbs;
	private transient double maxJumpDist;

	private double r0;
	private double distDelta;

	public RelativeRupJumpDistConstraint(FaultSystemRupSet rupSet, double r0, double distDelta, double weight, boolean inequality) {
		super("Relative Rupture Jump Dist, R₀="+oDF.format(r0), "RelRupJumpR0"+oDF.format(r0),
				weight, inequality, ConstraintWeightingType.UNNORMALIZED);
		this.rupSet = rupSet;
		Preconditions.checkState(r0 > 0);
		Preconditions.checkState(distDelta > 0);
		this.r0 = r0;
		this.distDelta = distDelta;
	}
	
	private synchronized void checkInitDistBinning() {
		if (distBinning == null) {
			ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
			
			double[] rupTotDists = new double[cRups.size()];
			
			for (int r=0; r<cRups.size(); r++) {
				ClusterRupture rup = cRups.get(r);
				
				double dist = 0d;
				for (Jump jump : rup.getJumpsIterable())
					dist += jump.distance;
				rupTotDists[r] = dist;
				maxJumpDist = Math.max(maxJumpDist, dist);
			}

			Preconditions.checkState(distDelta > 0);
			distBinning = HistogramFunction.getEncompassingHistogram(
					distDelta*1e-2, Math.max(maxJumpDist, 1.9d*distDelta), distDelta);

			ruptureBins = new int[rupTotDists.length];
			ruptureCondProbs = new double[rupTotDists.length];
			for (int r=0; r<rupTotDists.length; r++) {
				ruptureBins[r] = distBinning.getClosestXIndex(rupTotDists[r]);
				ruptureCondProbs[r] = Shaw07JumpDistProb.calcJumpProbability(rupTotDists[r], 1d, r0);
			}
		}
	}

	@Override
	public int getNumRows() {
		checkInitDistBinning();
		return distBinning.size();
	}

	@Override
	public long encode(DoubleMatrix2D A, double[] d, int startRow) {
		long count = 0;
		
		checkInitDistBinning();
		
		// estimate total supra-seis event rate
		IncrementalMagFreqDist targetMFD = null;
		if (rupSet.hasModule(InversionTargetMFDs.class))
			targetMFD = rupSet.requireModule(InversionTargetMFDs.class).getTotalOnFaultSupraSeisMFD();
		if (targetMFD == null)
			targetMFD = Inversions.inferTargetGRFromSlipRates(rupSet, 1d);
		double estTotRate = targetMFD.getTotalIncrRate();
		
		System.out.println("Estimated total supra-seis rate: "+(float)estTotRate);
		
		// scale weight by that estimated total event rate
		double weight = this.weight/estTotRate;
		
		double sumProb = StatUtils.sum(ruptureCondProbs);
		
		Preconditions.checkState(maxJumpDist > 0d);
		
//		double totalRange = Shaw07JumpDistProb.calcJumpProbability(0d, 1d, r0)
//				- Shaw07JumpDistProb.calcJumpProbability(maxJumpDist, 1d, r0);
		
		for (int i=0; i<distBinning.size(); i++) {
			int row = startRow + i;
			
			double binCenter = distBinning.getX(i);
			double binStart = binCenter - 0.5*distBinning.getDelta();
			double binEnd = binCenter + 0.5*distBinning.getDelta();
			
			double sumInBin = 0d;
			for (int r=0; r<ruptureBins.length; r++)
				if (ruptureBins[r] == i)
					sumInBin += ruptureCondProbs[r];
			
			double fractInBin = sumInBin/sumProb;
			
//			double p0 = Shaw07JumpDistProb.calcJumpProbability(binStart, 1d, r0);
//			double p1 = Shaw07JumpDistProb.calcJumpProbability(binEnd, 1d, r0);
//			
//			double fractInBin = (p0-p1)/totalRange;
			
			double scalarIn = weight*(1d-fractInBin);
			double scalarOut = -weight*fractInBin;
			
			int numInBin = 0;
			for (int r=0; r<ruptureBins.length; r++) {
				if (ruptureBins[r] == i) {
					setA(A, row, r, scalarIn);
					numInBin++;
				} else {
					setA(A, row, r, scalarOut);
				}
			}
			count += ruptureBins.length;
			System.out.println("Bin "+i+", rups="+numInBin+", ["+(float)binStart+","+(float)binEnd
					+"], sumInBin="+(float)sumInBin+", fract in bin: "+(float)fractInBin);
//			System.out.println("Bin "+i+", rups="+numInBin+", ["+(float)binStart+","+(float)binEnd
//					+"], prob range: ["+(float)p0+","+(float)p1+"], prob in bin: "+(float)fractInBin);
			
			d[row] = 0d;
		}
		
		return count;
	}

	@Override
	public void setRuptureSet(FaultSystemRupSet rupSet) {
		this.rupSet = rupSet;
	}

}
