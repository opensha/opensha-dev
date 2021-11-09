package scratch.nshm23.targetMFDs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
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
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBValuePlot;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSystemRupSetCalc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.UCERF2_ComparisonSolutionFetcher;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.UCERF2_A_FaultMapper;

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
	
	private List<Integer> getRupIndexesBelowMinMag() {
		ModSectMinMags modMinMags = rupSet.requireModule(ModSectMinMags.class);
		Preconditions.checkNotNull(modMinMags, "Rupture set must supply ModSectMinMags if minimization constraint is enabled");
		
		// we want to only grab ruptures with magnitudes below the MFD bin in which the section minimium magnitude resides
		EvenlyDiscretizedFunc refMagFunc = new EvenlyDiscretizedFunc(InversionTargetMFDsFromBValAndDefModel.MIN_MAG,
				InversionTargetMFDsFromBValAndDefModel.NUM_MAG, InversionTargetMFDsFromBValAndDefModel.DELTA_MAG);
		double halfDelta = refMagFunc.getDelta()*0.5;
		double[] mfdMappedSectMins = new double[rupSet.getNumSections()];
		for (int s=0; s<mfdMappedSectMins.length; s++)
			mfdMappedSectMins[s] = refMagFunc.getX(refMagFunc.getClosestXIndex(modMinMags.getMinMagForSection(s)))-halfDelta;
		
		List<Integer> belowMinIndexes = new ArrayList<>();
		float maxMin = (float)StatUtils.max(modMinMags.getMinMagForSections());
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double mag = rupSet.getMagForRup(r);
			if ((float)mag >= maxMin)
				continue;
			boolean below = false;
			for (int s : rupSet.getSectionsIndicesForRup(r)) {
				if ((float)mag < (float)mfdMappedSectMins[s]) {
					below = true;
					break;
				}
//				int m = refMagFunc.getClosestXIndex(mag);
//				float sectMin = (float)modMinMags.getMinMagForSection(s);
//				if (sectMin > (float)(refMagFunc.getX(m)-halfDelta)) {
//					below = true;
//					break;
//				}
			}
			if (below)
				belowMinIndexes.add(r);
		}
		System.out.println("Found "+belowMinIndexes.size()+" ruptures below sect min mags");
		return belowMinIndexes;
	}
	
	public IntegerPDF_FunctionSampler getSkipBelowMinSampler() {
		double[] vals = new double[rupSet.getNumRuptures()];
		Arrays.fill(vals, 1d);
		for (int index : getRupIndexesBelowMinMag())
			vals[index] = 0d;
		return new IntegerPDF_FunctionSampler(vals);
	}
	
	public DraftModelConstraintBuilder minimizeBelowSectMinMag() {
		List<Integer> belowMinIndexes = getRupIndexesBelowMinMag();
		constraints.add(new RupRateMinimizationConstraint(100000, belowMinIndexes));
		return this;
	}
	
	public DraftModelConstraintBuilder supraSmooth() {
		constraints.add(new LaplacianSmoothingInversionConstraint(rupSet, 1000));
		return this;
	}
	
	public DraftModelConstraintBuilder supraPaleoSmooth() {
		HashSet<Integer> paleoParentIDs = new HashSet<>();
		
		PaleoseismicConstraintData paleoData = rupSet.requireModule(PaleoseismicConstraintData.class);
		List<SectMappedUncertainDataConstraint> paleos = new ArrayList<>();
		if (paleoData.hasPaleoRateConstraints())
			paleos.addAll(paleoData.getPaleoRateConstraints());
		if (paleoData.hasPaleoSlipConstraints())
			paleos.addAll(paleoData.getPaleoSlipConstraints());
		Preconditions.checkState(!paleos.isEmpty());
		for (SectMappedUncertainDataConstraint paleo : paleos)
			paleoParentIDs.add(rupSet.getFaultSectionData(paleo.sectionIndex).getParentSectionId());
		constraints.add(new LaplacianSmoothingInversionConstraint(rupSet, 1000, paleoParentIDs));
		return this;
	}
	
	/**
	 * sets the weight of all constraints of the given type
	 * 
	 * @param weight
	 * @return
	 */
	public DraftModelConstraintBuilder weight(Class<? extends InversionConstraint> clazz, double weight) {
		for (int i=constraints.size(); --i>=0;) {
			InversionConstraint constraint = constraints.get(i);
			if (clazz.isAssignableFrom(constraint.getClass()))
				constraint.setWeight(weight);
		}
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
	
	public DraftModelConstraintBuilder testFlipBVals(FaultSystemSolution prevSol, double targetBVal) {
		Preconditions.checkState(rupSet.isEquivalentTo(prevSol.getRupSet()));
		InversionTargetMFDsFromBValAndDefModel targetMFDs = getTargetMFDs(targetBVal);
		
		List<IncrementalMagFreqDist> origSupraNuclMFDs = targetMFDs.getSectSupraSeisNuclMFDs();
		
		double[] solRupMoRates = SectBValuePlot.calcRupMomentRates(prevSol);

		double[] targetNuclRates = new double[rupSet.getNumSections()];
		double[] targetNuclRateStdDevs = new double[rupSet.getNumSections()];
		
		MinMaxAveTracker solBVals = new MinMaxAveTracker();
		MinMaxAveTracker flippedBVals = new MinMaxAveTracker();
		for (int s=0; s<targetNuclRates.length; s++) {
			IncrementalMagFreqDist origMFD = origSupraNuclMFDs.get(s);
			
			double minMag = Double.POSITIVE_INFINITY;
			double maxMag = 0d;
			for (int i=0; i<origMFD.size(); i++) {
				if (origMFD.getY(i) > 0) {
					double x = origMFD.getX(i);
					minMag = Math.min(minMag, x);
					maxMag = Math.max(maxMag, x);
				}
			}
			
			double solSupraNuclRate = 0d;
			double solSupraMoRate = 0d;
			double sectArea = rupSet.getAreaForSection(s);
			for (int r : rupSet.getRupturesForSection(s)) {
				double rate = prevSol.getRateForRup(r);
				if (rate > 0) {
					double fract = sectArea/rupSet.getAreaForRup(r);
					solSupraNuclRate += rate*fract;
					solSupraMoRate += solRupMoRates[r]*fract;
				}
			}
			
			int numOrig = 1+origMFD.getClosestXIndex(maxMag)-origMFD.getClosestXIndex(minMag);
			GutenbergRichterMagFreqDist grWithSolRate = new GutenbergRichterMagFreqDist(
					minMag, numOrig, origMFD.getDelta());
			grWithSolRate.setAllButBvalue(minMag, maxMag, solSupraMoRate, solSupraNuclRate);
			
			double solBVal = grWithSolRate.get_bValue();
			solBVals.addValue(solBVal);
			
			double flippedBVal = targetBVal + (targetBVal-solBVal);
			flippedBVals.addValue(flippedBVal);
			
			GutenbergRichterMagFreqDist grWithFlippedB = new GutenbergRichterMagFreqDist(
					minMag, numOrig, origMFD.getDelta());
			grWithFlippedB.setAllButTotCumRate(minMag, maxMag, origMFD.getTotalMomentRate(), flippedBVal);
			
			double flippedNuclRate = grWithFlippedB.getTotalIncrRate();
			
			System.out.println(s+". targetB="+(float)targetBVal+"\ttargetRate="+(float)origMFD.getTotalIncrRate()
					+"\ttargetMoRate="+(float)origMFD.getTotalMomentRate()+"\tsolMoRate="+(float)solSupraMoRate
					+"\tsolRate="+(float)solSupraNuclRate+"\tsolBVal="+(float)solBVal
					+"\tflippedBVal="+(float)flippedBVal+"\tflippedRate="+(float)flippedNuclRate);
			targetNuclRates[s] = flippedNuclRate;
			targetNuclRateStdDevs[s] = 0.1*targetNuclRates[s];
		}
		System.out.println("Solution b-values: "+solBVals);
		System.out.println("Flipped b-values: "+flippedBVals);
		constraints.add(new SectionTotalRateConstraint(rupSet, 1d,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, targetNuclRates, targetNuclRateStdDevs, true));
		return this;
	}
	
	public DraftModelConstraintBuilder testSameBVals(FaultSystemSolution prevSol) {
		Preconditions.checkState(rupSet.isEquivalentTo(prevSol.getRupSet()));
		double[] targetNuclRates = new double[rupSet.getNumSections()];
		double[] targetNuclRateStdDevs = new double[rupSet.getNumSections()];
		for (int s=0; s<targetNuclRates.length; s++) {
			double solSupraNuclRate = 0d;
			double sectArea = rupSet.getAreaForSection(s);
			for (int r : rupSet.getRupturesForSection(s))
				solSupraNuclRate += prevSol.getRateForRup(r)*sectArea/rupSet.getAreaForRup(r);
			
			targetNuclRates[s] = solSupraNuclRate;
			targetNuclRateStdDevs[s] = 0.1*targetNuclRates[s];
		}
		constraints.add(new SectionTotalRateConstraint(rupSet, 1d,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, targetNuclRates, targetNuclRateStdDevs, true));
		return this;
	}
	
	public IntegerPDF_FunctionSampler testGetSampleAllNew(FaultSystemSolution prevSol, boolean skipBelow) {
		Preconditions.checkState(rupSet.isEquivalentTo(prevSol.getRupSet()));
		double[] weights = new double[rupSet.getNumRuptures()];
		Arrays.fill(weights, 1d);
		if (skipBelow)
			for (int r : getRupIndexesBelowMinMag())
				weights[r] = 0;
		for (int r=0; r<weights.length; r++)
			if (prevSol.getRateForRup(r) > 0d)
				weights[r] = 0;
		return new IntegerPDF_FunctionSampler(weights);
	}
	
	public DraftModelConstraintBuilder u2NuclBVals(boolean aFaults) {
		U3LogicTreeBranch branch = rupSet.requireModule(U3LogicTreeBranch.class);
		AveSlipModule aveSlipModule = rupSet.requireModule(AveSlipModule.class);
		double fractGR = 0.33333;
				
		FaultSystemSolution UCERF2_FltSysSol = null;
		if (aFaults)
			UCERF2_FltSysSol = UCERF2_ComparisonSolutionFetcher.getUCERF2Solution(
					rupSet, branch.getValue(FaultModels.class), aveSlipModule);
		
		double[] targetNuclRates = new double[rupSet.getNumSections()];
		double[] targetNuclRateStdDevs = new double[rupSet.getNumSections()];
		
		// targets with b=1
		InversionTargetMFDsFromBValAndDefModel targetB1 = new InversionTargetMFDsFromBValAndDefModel(rupSet, 1d);

		MinMaxAveTracker implBValTrack = new MinMaxAveTracker();
		MinMaxAveTracker implBValAFaultTrack = new MinMaxAveTracker();
		
		for (int s=0;s <rupSet.getNumSections(); s++) {
			FaultSection data = rupSet.getFaultSectionData(s);
			IncrementalMagFreqDist sectNuclB1 = targetB1.getSectSupraSeisNuclMFDs().get(s);
			int minIndex = -1;
			int maxIndex = 0;
			for (int i=0; i<sectNuclB1.size(); i++) {
				if (sectNuclB1.getY(i) > 0d) {
					maxIndex = i;
					if (minIndex < 0)
						minIndex = i;
				}
			}
			double maxMag = sectNuclB1.getX(maxIndex);
			double totSectMoment = sectNuclB1.getTotalMomentRate();
			
			boolean aFault = false;
			
			if (aFaults && UCERF2_A_FaultMapper.wasUCERF2_TypeAFault(data.getParentSectionId())) {
				// type a fault, use UCERF2
				double solSupraNuclRate = 0d;
				double sectArea = rupSet.getAreaForSection(s);
				for (int r : UCERF2_FltSysSol.getRupSet().getRupturesForSection(s))
					solSupraNuclRate += UCERF2_FltSysSol.getRateForRup(r)*sectArea/UCERF2_FltSysSol.getRupSet().getAreaForRup(r);
				if (solSupraNuclRate == 0d) {
					System.err.println("WARNING: no a-fault nucl rate for "+s+". "+data.getName()+", reverting to b-fault formula");
				} else {
					Preconditions.checkState(solSupraNuclRate > 0d, "nucl rate is zero for a-fault: %s. %s", s, data.getName());
					targetNuclRates[s] = solSupraNuclRate;
					targetNuclRateStdDevs[s] = 0.1*targetNuclRates[s];
					aFault = true;
				}
			}
			
			if (targetNuclRates[s] == 0d) {
				// b fault, use functional form
				GutenbergRichterMagFreqDist charMFD = new GutenbergRichterMagFreqDist(
						sectNuclB1.getMinX(), sectNuclB1.size(), sectNuclB1.getDelta());
				
				double charMoment = (1d-fractGR)*totSectMoment;
				charMFD.setAllButTotCumRate(maxMag, maxMag, charMoment, 1d);
				IncrementalMagFreqDist grMFD = sectNuclB1.deepClone();
				grMFD.scaleToTotalMomentRate(totSectMoment*fractGR);
				Preconditions.checkState((float)(grMFD.getTotalMomentRate()+charMFD.getTotalMomentRate()) == (float)totSectMoment);
				SummedMagFreqDist u2Target = new SummedMagFreqDist(charMFD.getMinX(), charMFD.size(), charMFD.getDelta());
				u2Target.addIncrementalMagFreqDist(grMFD);
				u2Target.addIncrementalMagFreqDist(charMFD);
				
				targetNuclRates[s] = u2Target.getTotalIncrRate();
				Preconditions.checkState(targetNuclRates[s] > 0d, "nucl rate is zero for b-fault: %s. %s", s, data.getName());
				targetNuclRateStdDevs[s] = 0.1*targetNuclRates[s];
			}
			GutenbergRichterMagFreqDist grWithSolRate = new GutenbergRichterMagFreqDist(
					sectNuclB1.getMinX(), sectNuclB1.size(), sectNuclB1.getDelta());
			grWithSolRate.setAllButBvalue(sectNuclB1.getX(minIndex), maxMag, totSectMoment, targetNuclRates[s]);
			
			double implB = grWithSolRate.get_bValue();
			
			System.out.println(s+".\trate="+(float)targetNuclRates[s]+"\tb="+(float)implB+"\taFault="+aFault);
			
			if (aFault)
				implBValAFaultTrack.addValue(implB);
			else
				implBValTrack.addValue(implB);
		}
		if (implBValAFaultTrack.getNum() > 0) {
			System.out.println("UCERF2-style implied b-values:");
			System.out.println("\t"+implBValAFaultTrack.getNum()+" a-faults: "+implBValAFaultTrack);
			System.out.println("\t"+implBValTrack.getNum()+" b-faults: "+implBValTrack);
		} else {
			System.out.println("UCERF2-style implied b-values: "+implBValTrack);
		}
		
		constraints.add(new SectionTotalRateConstraint(rupSet, 1d,
				ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, targetNuclRates, targetNuclRateStdDevs, true));
		return this;
	}

}
