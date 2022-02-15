package scratch.kevin.nshm23.segModelTests;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class NuclMFD_SegDebug {

	public static void main(String[] args) throws IOException {
		File inputFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_02_08-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-SubB1-2000ip/results/"
				+ "FM3_1_CoulombRupSet_U3_ABM_EllBsqrtLen_DsrUni_SupraB0.0_NuclMFD_SubB1_ShawR0_2/solution.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		factory.buildInversionConfig(rupSet, rupSet.requireModule(LogicTreeBranch.class), 16);
		
		SupraSeisBValInversionTargetMFDs targetMFDs = rupSet.requireModule(SupraSeisBValInversionTargetMFDs.class);
		
		List<UncertainIncrMagFreqDist> supraSeisTargets = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
		
		factory.setAdjustTargetsForSegmentation(false);
		factory.buildInversionConfig(rupSet, rupSet.requireModule(LogicTreeBranch.class), 16);
		List<UncertainIncrMagFreqDist> noAdjSupraSeisTargets = rupSet.requireModule(SupraSeisBValInversionTargetMFDs.class)
				.getOnFaultSupraSeisNucleationMFDs();
		
		double printThreshold = 1000d;
		
		DecimalFormat eDF = new DecimalFormat("0.0E0");
		DecimalFormat zDF = new DecimalFormat("0.0");
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			UncertainIncrMagFreqDist targetMFD = supraSeisTargets.get(s);
			UncertainIncrMagFreqDist noAdjTargetMFD = noAdjSupraSeisTargets.get(s);
			IncrementalMagFreqDist solMFD = sol.calcNucleationMFD_forSect(s, targetMFD.getMinX(), targetMFD.getMaxX(), targetMFD.size());
			
			boolean exceeded = false;
			for (int i=0; i<targetMFD.size(); i++) {
				double targetY = targetMFD.getY(i);
				double solY = solMFD.getY(i);
				double targetStdDev = targetMFD.getStdDev(i);
				
				double noAdjTargetY = noAdjTargetMFD.getY(i);
				double noAdjTargetStdDev = noAdjTargetMFD.getStdDev(i);
				
				if (targetY > 0 || solY > 0) {
					double zScore = (solY - targetY) / targetStdDev;
					if (Math.abs(zScore) > printThreshold || targetStdDev <= 0d) {
						if (!exceeded) {
							System.out.println("Exceeded for s="+s+", "+rupSet.getFaultSectionData(s).getSectionName());
							System.out.println("\tMag\tTarget\tSol\tStd.\tZ\tNoAdjY\tStd\tZ");
						}
						exceeded = true;

						double noAdjzScore = (solY - noAdjTargetY) / noAdjTargetStdDev;
						System.out.println("\tM="+(float)targetMFD.getX(i)+"\t"+eDF.format(targetY)
							+"\t"+eDF.format(solY)+"\t"+eDF.format(targetStdDev)+"\t"
							+(Double.isFinite(zScore) ? zDF.format(zScore)+"" : (zScore+"").replace("Infinity", "Inf"))
							+"\t"+eDF.format(noAdjTargetY)+"\t"+eDF.format(noAdjTargetStdDev)+"\t"
							+(Double.isFinite(noAdjzScore) ? zDF.format(noAdjzScore)+"" : (noAdjzScore+"").replace("Infinity", "Inf")));
					}
				}
			}
		}
		
		
	}

}
