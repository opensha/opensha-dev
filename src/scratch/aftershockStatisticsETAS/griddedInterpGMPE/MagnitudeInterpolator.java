package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.util.HashSet;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;

import com.google.common.base.Preconditions;

//public class MagnitudeInterpolator extends AbstractGMPEInterpolation.EvenlySpacedDouble {
//	
//	private final double tolerance;
//
//	public MagnitudeInterpolator(double minMag, double maxMag, int numBins) {
//		super("Magnitude", minMag, maxMag, numBins, false, false);
//		double delta = (maxMag - minMag) / (numBins - 1);
//		// set tolerance to only get ones that line up nicely
//		tolerance = 0.4*delta;
//	}
//
//	@Override
//	public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {}
//
//	@Override
//	public Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index) {
//		double mag = getValue(index);
//		HashSet<EqkRupture> viableRups = new HashSet<>();
//		double closestMag = Double.NaN;
//		double closestDiff = Double.POSITIVE_INFINITY;
//		for (EqkRupture rup : ruptures) {
//			double diff = Math.abs(rup.getMag() - mag);
//			if (diff < tolerance)
//				viableRups.add(rup);
//			if (diff < closestDiff) {
//				closestDiff = diff;
//				closestMag = rup.getMag();
//			}
//		}
//		Preconditions.checkState(!viableRups.isEmpty(), "No rups matching M=%s found. Closest was %s away at M=%s",
//				(float)mag, (float)closestDiff, (float)closestMag);
//		return viableRups;
//	}
//
//	@Override
//	public Double detectCurrentVal(ScalarIMR gmpe, Site site) {
//		return gmpe.getEqkRupture().getMag();
//	}
//
//}
