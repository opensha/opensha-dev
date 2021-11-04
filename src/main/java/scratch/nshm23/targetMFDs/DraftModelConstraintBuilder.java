package scratch.nshm23.targetMFDs;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.LaplacianSmoothingInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.RupRateMinimizationConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSystemRupSetCalc;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;

public class DraftModelConstraintBuilder {
	
	private List<InversionConstraint> constraints;
	private FaultSystemRupSet rupSet;
	
	public DraftModelConstraintBuilder(FaultSystemRupSet rupSet) {
		this.rupSet = rupSet;
		this.constraints = new ArrayList<>();
	}
	
	public List<InversionConstraint> build() {
		return constraints;
	}
	
	public DraftModelConstraintBuilder defaultConstraints() {
		return defaultConstraints(0.8);
	}
	
	public DraftModelConstraintBuilder defaultConstraints(double supraBVal) {
		// TODO: section rate constraints
		return defaultDataConstraints(supraBVal).defaultMetaConstraints();
	}
	
	public DraftModelConstraintBuilder except(Class<? extends InversionConstraint> clazz) {
		for (int i=constraints.size(); --i>=0;)
			if (clazz.isAssignableFrom(constraints.get(i).getClass()))
				constraints.remove(i);
		return this;
	}
	
	public DraftModelConstraintBuilder defaultDataConstraints() {
		return defaultDataConstraints(0.8);
	}
	
	public DraftModelConstraintBuilder defaultDataConstraints(double supraBVal) {
		return slipRates().paleo().parkfield().supraBValMFDs(supraBVal).sectSupraRates(supraBVal);
	}
	
	public DraftModelConstraintBuilder slipRates() {
		constraints.add(new SlipRateInversionConstraint(1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, rupSet));
		return this;
	}
	
	public DraftModelConstraintBuilder paleo() {
		return paleoRates().paleoSlips();
	}
	
	public DraftModelConstraintBuilder paleoRates() {
		PaleoseismicConstraintData data = rupSet.requireModule(PaleoseismicConstraintData.class);
		constraints.add(new PaleoRateInversionConstraint(rupSet, 1d, data.getPaleoRateConstraints(), data.getPaleoProbModel()));
		return this;
	}
	
	public DraftModelConstraintBuilder paleoSlips() {
		PaleoseismicConstraintData data = rupSet.requireModule(PaleoseismicConstraintData.class);
		constraints.add(new PaleoSlipInversionConstraint(rupSet, 1d, data.getPaleoSlipConstraints(), data.getPaleoSlipProbModel(), true));
		return this;
	}
	
	private InversionTargetMFDsFromBValAndDefModel getTargetMFDs(double supraBVal) {
		InversionTargetMFDsFromBValAndDefModel target = rupSet.getModule(InversionTargetMFDsFromBValAndDefModel.class);
		if (target == null || target.getSupraSeisBValue() != supraBVal) {
			target = new InversionTargetMFDsFromBValAndDefModel(rupSet, supraBVal);
			rupSet.addModule(target);
		}
		return target;
	}
	
	public DraftModelConstraintBuilder supraBValMFDs(double supraBVal) {
		InversionTargetMFDsFromBValAndDefModel target = getTargetMFDs(supraBVal);
		
		List<? extends IncrementalMagFreqDist> origMFDs = target.getMFD_Constraints();
		List<UncertainIncrMagFreqDist> uncertainMFDs = new ArrayList<>();
		System.err.println("WARNING: temporary relative standard deviation of 0.1 set for all MFD bins"); // TODO
		for (IncrementalMagFreqDist mfd : origMFDs) {
			// TODO fake standard deviation
			UncertainIncrMagFreqDist uMFD = UncertainIncrMagFreqDist.constantRelStdDev(mfd, 0.1);
//			for (int i=0; i<uMFD.size(); i++)
//				System.out.println(uMFD.getX(i)+"\t"+uMFD.getY(i)+"\t"+uMFD.getStdDev(i));
			uncertainMFDs.add(uMFD);
		}
		constraints.add(new MFDInversionConstraint(rupSet, 1d, false,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, uncertainMFDs));
		return this;
	}
	
	public DraftModelConstraintBuilder sectSupraRates(double supraBVal) {
		InversionTargetMFDsFromBValAndDefModel target = getTargetMFDs(supraBVal);
		
		double[] targetRates = new double[rupSet.getNumSections()];
		double[] targetRateStdDevs = new double[rupSet.getNumSections()];
		
		System.err.println("WARNING: temporary relative standard deviation of 0.1 set for all section target rates"); // TODO
		
		List<IncrementalMagFreqDist> sectSupraMFDs = target.getSectSupraSeisNuclMFDs();
		for (int s=0; s<targetRates.length; s++) {
			targetRates[s] = sectSupraMFDs.get(s).calcSumOfY_Vals();
			targetRateStdDevs[s] = 0.1*targetRates[s];
//			System.out.println(rupSet.getFaultSectionData(s).getSectionName()+": total rate: "+targetRates[s]);
		}
		
		constraints.add(new SectionTotalRateConstraint(rupSet, 1d,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, targetRates, targetRateStdDevs, true));
		return this;
	}
	
	public DraftModelConstraintBuilder parkfield() {
		double parkfieldMeanRate = 1.0/25.0; // Bakun et al. (2005)
		System.err.println("WARNING: temporary relative standard deviation of 0.1 set for parkfield constraint"); // TODO
		double parkfieldStdDev = 0.1d*parkfieldMeanRate;
		
		// Find Parkfield M~6 ruptures
		List<Integer> parkfieldRups = UCERF3InversionInputGenerator.findParkfieldRups(rupSet);
		constraints.add(new ParkfieldInversionConstraint(1d, parkfieldMeanRate, parkfieldRups,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, parkfieldStdDev));
		return this;
	}
	
	public DraftModelConstraintBuilder defaultMetaConstraints() {
		return minimizeBelowSectMinMag();
	}
	
	public DraftModelConstraintBuilder minimizeBelowSectMinMag() {
		ModSectMinMags modMinMags = rupSet.requireModule(ModSectMinMags.class);
		Preconditions.checkNotNull(modMinMags, "Rupture set must supply ModSectMinMags if minimization constraint is enabled");
		List<Integer> belowMinIndexes = new ArrayList<>();
		for (int r=0; r<rupSet.getNumRuptures(); r++)
//			if (rupSet.isRuptureBelowSectMinMag(r))
			if (FaultSystemRupSetCalc.isRuptureBelowSectMinMag(rupSet, r, modMinMags))
				belowMinIndexes.add(r);
		constraints.add(new RupRateMinimizationConstraint(100000, belowMinIndexes));
		return this;
	}
	
	public DraftModelConstraintBuilder supraSmooth() {
		constraints.add(new LaplacianSmoothingInversionConstraint(rupSet, 1000));
		return this;
	}
	
	/**
	 * sets the weight of the most recently added constraint
	 * 
	 * @param weight
	 * @return
	 */
	public DraftModelConstraintBuilder weight(double weight) {
		Preconditions.checkState(!constraints.isEmpty());
		constraints.get(constraints.size()-1).setWeight(weight);
		return this;
	}

}
