package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.util.Interpolate;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc.BinaryJumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.GRParticRateEstimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class BvalAndSegConnectivityCorrelationSampler extends SectionCovarianceSampler {
	
	private SectionDistanceAzimuthCalculator distCalc;
	private double maxDist;
	private double zeroDistCoeff;
	private double negativeCorrMaxDist;
	private boolean fromMFDs;
	
	private FaultSystemRupSet rupSet;
	private double bVal;
	private NSHM23_SegmentationModels segModel;
	
	private double[] totRates;
	private double[][] coruptureRates;
	
	public BvalAndSegConnectivityCorrelationSampler(List<? extends FaultSection> subSects, FaultSystemRupSet rupSet,
			SectionDistanceAzimuthCalculator distCalc,
			double maxDist, double zeroDistCoeff, double negativeCorrMaxDist,
			double bVal, NSHM23_SegmentationModels segModel, boolean fromMFDs) {
		super(subSects);
		this.rupSet = rupSet;
		this.distCalc = distCalc;
		this.maxDist = maxDist;
		this.zeroDistCoeff = zeroDistCoeff;
		this.negativeCorrMaxDist = negativeCorrMaxDist;
		this.bVal = bVal;
		this.segModel = segModel;
		this.fromMFDs = fromMFDs;
	}
	
	private void calcCoruptureRates(FaultSystemRupSet rupSet, double bVal, NSHM23_SegmentationModels segModel) {
		System.out.println("Estimating rupture rates for rupSet with "+rupSet.getNumRuptures()+" ruptures, "
				+rupSet.getNumSections()+" sections, b="+(float)bVal+", and segmentation branch "+segModel);
		int numSects = rupSet.getNumSections();
		int numRups = rupSet.getNumRuptures();
		double[] totRates = new double[numSects];
		double[][] coruptureRates = new double[numSects][numSects];
		
		if (fromMFDs) {
			// build corupture rates for each in the view of section, then average them
			SupraSeisBValInversionTargetMFDs.Builder builder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal);
			if (segModel != null) {
				JumpProbabilityCalc segModelCalc = segModel.getModel(rupSet, null);
//				if (segModelCalc instanceof BinaryJumpProbabilityCalc)
//					builder.forBinaryRupProbModel((BinaryJumpProbabilityCalc)segModel.getExclusionModel(rupSet, null));
//				else
					builder.adjustTargetsForData(new ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb(segModelCalc, 100, true));
			}
			builder.applyDefModelUncertainties(false);
			SupraSeisBValInversionTargetMFDs mfds = builder.build();
			
			System.out.println("Calculating estimated connection rates for "+numSects+" sects");

			for (int s=0; s<numSects; s++) {
				IncrementalMagFreqDist nuclGR = mfds.getOnFaultSupraSeisNucleationMFDs().get(s);

				List<Integer> rups = mfds.getRupturesForSect(s);
				List<Double> rupMags = new ArrayList<>(rups.size());
				int[] rupsPerBin = new int[nuclGR.size()];
				for (int r : rups) {
					double rupMag = rupSet.getMagForRup(r);
					rupMags.add(rupMag);
					rupsPerBin[nuclGR.getClosestXIndex(rupMag)]++;
				}

				if (rups.isEmpty()) {
					Preconditions.checkState(nuclGR.calcSumOfY_Vals() == 0d);
					continue;
				}

				double sectArea = rupSet.getAreaForSection(s);

				// spread to all ruptures evenly to get partic rate
				double calcRate = 0d;
				for (int r=0; r<rups.size(); r++) {
					int rupIndex = rups.get(r);
					int bin = nuclGR.getClosestXIndex(rupMags.get(r));
					/// this is a nucleation rate
					double nuclRate = nuclGR.getY(bin)/(double)rupsPerBin[bin];
					// turn back into participation rate
					double particRate = nuclRate*rupSet.getAreaForRup(rupIndex)/sectArea;
					// adjust for visibility
					calcRate += particRate;

					for (int oSect : rupSet.getSectionsIndicesForRup(rupIndex))
						coruptureRates[s][oSect] += particRate;
				}
				totRates[s] = calcRate;
			}
		} else {
			// estimate rupture rates
			GRParticRateEstimator est = new GRParticRateEstimator(
					rupSet, bVal, segModel == null ? null : segModel.getModel(rupSet, null));

			double[] rupRates = est.estimateRuptureRates();

			for (int r=0; r<numRups; r++) {
				double rate = rupRates[r];
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
		this.coruptureRates = coruptureRates;
		this.totRates = totRates;
	}

	@Override
	public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
		if (totRates == null) {
			synchronized (this) {
				if (totRates == null) {
					calcCoruptureRates(rupSet, bVal, segModel);
				}
			}
		}
		double dist = distCalc.getDistance(sect1, sect2);
		if (dist > maxDist)
			return 0;
		int id1 = sect1.getSectionId();
		int id2 = sect2.getSectionId();
		if (coruptureRates[id1][id2] == 0d) {
			// unconnected
			Preconditions.checkState(coruptureRates[id2][id1] == 0d);
			if (dist < negativeCorrMaxDist) {
				// anticorrelated
				return Interpolate.findY(0d, -1d, negativeCorrMaxDist, 0d, dist);
			} else {
				// uncorrelated
				return 0d;
			}
		} else {
			// correlated, multiply connected fraction with distance fraction
			Preconditions.checkState(coruptureRates[id2][id1] > 0d);
			// in the world of sect 1
			double connRate1 = coruptureRates[id1][id2]/totRates[id1];
			// in the world of sect 2
			double connRate2 = coruptureRates[id2][id1]/totRates[id2];
			double connFract = 0.5*(connRate1 + connRate2);
			Preconditions.checkState(connFract >= 0 && connFract <= 1d);
			return connFract * Interpolate.findY(0d, zeroDistCoeff, maxDist, 0d, dist);
		}
	}

	@Override
	protected String getSamplerPrefix() {
		String prefix = "conn_corr_b"+(float)bVal+"_"+segModel.getFilePrefix();
		if (fromMFDs)
			prefix += "_fromMFDs";
		prefix += "_dist"+(float)maxDist+"km_zeroCoeff"+(float)zeroDistCoeff;
		if (negativeCorrMaxDist > 0)
			prefix += "_negCorrDist"+(float)negativeCorrMaxDist+"km";
		return prefix;
	}
	

}
