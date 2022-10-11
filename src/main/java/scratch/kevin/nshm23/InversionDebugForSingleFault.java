package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.IntegerSampler;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ColumnOrganizedAnnealingData;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.InversionState;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class InversionDebugForSingleFault {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/solution.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		int parentID = 85;
		AveSlipModule slips = rupSet.requireModule(AveSlipModule.class);
		
		LogicTreeBranch<?> branch = rupSet.requireModule(LogicTreeBranch.class);
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		BinaryRuptureProbabilityCalc exclusion = NSHM23_InvConfigFactory.getExclusionModel(
				rupSet, branch, cRups);
		
		double rateSum = 0d;
		for (int rupIndex : rupSet.getRupturesForParentSection(parentID)) {
			ClusterRupture rup = cRups.get(rupIndex);
			System.out.println(rupIndex+"\tmag="+(float)rupSet.getMagForRup(rupIndex)
				+"\tslip="+(float)slips.getAveSlip(rupIndex)+"\trate="+(float)sol.getRateForRup(rupIndex)
				+"\texclusion="+(exclusion == null ? "null" : exclusion.isRupAllowed(rup, false)));
			rateSum += sol.getRateForRup(rupIndex);
		}
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, 16);
		
		for (InversionConstraint constr : config.getConstraints()) {
			if (constr instanceof MFDInversionConstraint)
				constr.setWeight(1e-8);
//			if (constr instanceof TotalRateInversionConstraint)
//				constr.setWeight(1e-8);
		}
		
		for (InversionConstraint constr : config.getConstraints())
			System.out.println(constr.getName()+": "+constr.getWeight());
		
		InversionInputGenerator inputs = new InversionInputGenerator(rupSet, config);
		inputs.generateInputs(true);
		inputs.columnCompress();
		
		ColumnOrganizedAnnealingData equalityData = new ColumnOrganizedAnnealingData(inputs.getA(), inputs.getD());
		ColumnOrganizedAnnealingData inequalityData = null;
		if (inputs.getA_ineq() != null)
			inequalityData = new ColumnOrganizedAnnealingData(inputs.getA_ineq(), inputs.getD_ineq());
		SerialSimulatedAnnealing sa = new SerialSimulatedAnnealing(equalityData, inequalityData, sol.getRateForAllRups(), 0d);
		
		double[] energy = sa.calculateEnergy(sol.getRateForAllRups());
		System.out.println("Energy: "+energy[0]);
		
		ColumnOrganizedAnnealingData data = new ColumnOrganizedAnnealingData(inputs.getA(), inputs.getD());
		double[] misfits = new double[data.nRows];
		SerialSimulatedAnnealing.calculateMisfit(data, sol.getRateForAllRups(), misfits);
		ConstraintRange slipRange = null;
		ConstraintRange nuclRateRange = null;
		for (ConstraintRange range : inputs.getConstraintRowRanges()) {
			if (range.name.contains("Slip Rate")) {
				slipRange = range;
			} else if (range.name.contains("Subsection Total")) {
				nuclRateRange = range;
			}
		}
		printSectMisfits(rupSet, parentID, misfits, slipRange, nuclRateRange);
		InversionTargetMFDs targets = rupSet.requireModule(InversionTargetMFDs.class);
		UncertainIncrMagFreqDist targetMFD = (UncertainIncrMagFreqDist)targets.getTotalOnFaultMFD();
		printTotalMFDMisfits(sol, targetMFD);
//		
//			
//				System.out.println(range.name);
//				double[] trimmed = new double[range.endRow - range.startRow];
//				System.out.println("\trows: "+trimmed.length);
//				System.arraycopy(misfits, range.startRow, trimmed, 0, trimmed.length);
//				MisfitStats stats = new MisfitStats(trimmed, range);
//				System.out.println("\tMAD: "+stats.absMean);
//				
//			}
//		}
		
		// now re-invert only allowing it to use those particular ruptures
		IntegerSampler sampler = new IntegerSampler.ArbitraryIntegerSampler(rupSet.getRupturesForParentSection(parentID));
		
		sa.setRuptureSampler(sampler);
		InversionState state = sa.iterate(1000000l);
		double[] finalEnergy = sa.getBestEnergy();
		System.out.println("Perturbs kept: "+state.numPerturbsKept);
		System.out.println("Final energy: "+finalEnergy[0]);
		double[] finalSol = sa.getBestSolution();
		
		double finalRateSum = 0d;
		System.out.println("Final rates:");
		for (int rupIndex : rupSet.getRupturesForParentSection(parentID)) {
			System.out.println(rupIndex+"\tmag="+(float)rupSet.getMagForRup(rupIndex)
			+"\tslip="+(float)slips.getAveSlip(rupIndex)+"\trate="+(float)sol.getRateForRup(rupIndex)+" -> "+(float)finalSol[rupIndex]);
			finalRateSum += finalSol[rupIndex];
		}
		
		System.out.println("Total rate change: "+rateSum+" -> "+finalRateSum+", diff="+(finalRateSum-rateSum));
		double[] modMisfits = sa.getBestMisfit();
		printSectMisfits(rupSet, parentID, modMisfits, slipRange, nuclRateRange);
		printTotalMFDMisfits(new FaultSystemSolution(rupSet, finalSol), targetMFD);
	}
	
	private static void printSectMisfits(FaultSystemRupSet rupSet, int parentID, double[] misfits,
			ConstraintRange slipRange, ConstraintRange nuclRateRange) {
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == parentID) {
				System.out.println("\t"+sect.getSectionId()+". "+sect.getSectionName()
					+": slipMisfit="+(float)misfits[slipRange.startRow+sect.getSectionId()]
					+": nuclMisfit="+(float)misfits[nuclRateRange.startRow+sect.getSectionId()]);
			}
		}
	}
	
	private static void printTotalMFDMisfits(FaultSystemSolution sol, UncertainIncrMagFreqDist targetMFD) {
		IncrementalMagFreqDist solMFD = sol.calcTotalNucleationMFD(targetMFD.getMinX(), targetMFD.getMaxX(), targetMFD.getDelta());
		System.out.println("Target MFD & misfits");
		for (int i=0; i<targetMFD.size(); i++) {
			double solY = solMFD.getY(i);
			double targetY = targetMFD.getY(i);
			if (solY == 0 && targetY == 0)
				continue;
			double targetSD = targetMFD.getStdDev(i);
			double solZ = (solY - targetY)/targetSD;
			System.out.println((float)solMFD.getX(i)+"\ttarget="+(float)targetY+"\tsd="+(float)targetSD
					+"\tfracSD="+(float)(targetSD/targetY)+"\tsol="+(float)solY+"\tz="+solZ);
		}
	}

}
