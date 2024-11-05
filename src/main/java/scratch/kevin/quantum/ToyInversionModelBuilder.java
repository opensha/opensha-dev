package scratch.kevin.quantum;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.CSVWriter;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionSolver;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.RuptureSubSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.GRParticRateEstimator;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class ToyInversionModelBuilder {
	
	public static void main(String[] args) throws IOException {
		boolean rebuildRupSet = false;
		
		File outputDir = new File("/home/kevin/markdown/inversions/2024_10_15-quantum_test_problem");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File rsFile = new File(outputDir, "rup_set.zip");
		FaultSystemRupSet rupSet;
		if (rsFile.exists() && !rebuildRupSet) {
			rupSet = FaultSystemRupSet.load(rsFile);
		} else {
			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
			factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
			LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
			branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
			rupSet = factory.buildRuptureSet(branch, 32);
			System.out.println("Building connectivity clusters");
			ConnectivityClusters clusters = ConnectivityClusters.build(rupSet);
			
			int targetParentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Likely");
			FaultSection targetSect = null;
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
				if (sect.getParentSectionId() == targetParentID) {
					targetSect = sect;
					break;
				}
			}
			Preconditions.checkNotNull(targetSect);
			ConnectivityCluster cluster = null;
			for (ConnectivityCluster testCluster : clusters) {
				if (testCluster.containsSect(targetSect)) {
					cluster = testCluster;
					break;
				}
			}
			
			System.out.println("Cluster for "+targetSect.getSectionName()+" has "+cluster.getNumRuptures()
					+" ruptures on "+cluster.getNumSections()+" sub-sections");
			
			rupSet = rupSet.getForSectionSubSet(cluster.getSectIDs());
			rupSet.removeModuleInstances(NamedFaults.class);
			rupSet.write(rsFile);
		}
		
		
		
		// write out a-priori connectivity estimates
		GRParticRateEstimator estimator = new GRParticRateEstimator(rupSet, 0.5, NSHM23_SegmentationModels.MID.getModel(rupSet, null));
		double[][] estCoeffs = calcConnCoeffs(rupSet, estimator.estimateRuptureRates());
		writeCoeffCSV(estCoeffs, new File(outputDir, "conn_coeffs_a_priori.csv"));
		
		List<InversionConstraint> constraints = new ArrayList<>();
		constraints.add(new SlipRateInversionConstraint(1d, ConstraintWeightingType.NORMALIZED, rupSet));
		
		CompletionCriteria completion = new IterationCompletionCriteria(1000000l);
//		CompletionCriteria completion = new IterationCompletionCriteria(1000l);
		
		InversionConfiguration config = InversionConfiguration.builder(constraints, completion)
//				.avgThreads(4, new IterationCompletionCriteria(10000l))
//				.threads(16).subCompletion(new IterationCompletionCriteria(1000l))
				.build();
		InversionInputGenerator inputGen = new InversionInputGenerator(rupSet, config);
		inputGen.generateInputs(true);
		
		InversionSolver.Default solver = new InversionSolver.Default();
		FaultSystemSolution sol = solver.run(rupSet, config, inputGen, null);

		sol.write(new File(outputDir, "solution.zip"));
		
		inputGen.writeArchive(new File(outputDir, "constraints.zip"), sol.getRateForAllRups(), false);
		
		double[][] solCoeffs = calcConnCoeffs(rupSet, sol.getRateForAllRups());
		writeCoeffCSV(solCoeffs, new File(outputDir, "conn_coeffs_sol.csv"));
	}
	
	public static double[][] calcConnCoeffs(FaultSystemRupSet rupSet, double[] rates) {
		int numSects = rupSet.getNumSections();
		int numRups = rupSet.getNumRuptures();
		Preconditions.checkState(numRups == rates.length);
		double[][] coeffs = new double[numSects][numSects];
		BitSet[] sectRups = new BitSet[numSects];
		for (int s=0; s<numSects; s++) {
			sectRups[s] = new BitSet(numRups);
			for (int r : rupSet.getRupturesForSection(s))
				sectRups[s].set(r);
		}
		
		for (int s1=0; s1<numSects; s1++) {
			for (int s2=0; s2<numSects; s2++) {
				if (s1 == s2) {
					coeffs[s1][s2] = 1d;
				} else {
					double sectParticRate = 0d;
					double coruptureRate = 0d;
					for (int r=0; r<numRups; r++) {
						if (sectRups[s1].get(r)) {
							sectParticRate += rates[r];
							if (sectRups[s2].get(r))
								coruptureRate += rates[r];
						}
					}
					coeffs[s1][s2] = coruptureRate/sectParticRate;
				}
			}
		}
		
		return coeffs;
	}
	
	private static void writeCoeffCSV(double[][] coeffs, File outputFile) throws IOException {
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
		CSVWriter csv = new CSVWriter(out, true);
		
		csv.write(List.of("Primary Section Index", "Connected Section Index", "Connection Coefficient"));
		for (int s1=0; s1<coeffs.length; s1++)
			for (int s2=0; s2<coeffs.length; s2++)
				csv.write(List.of(s1+"", s2+"", (float)coeffs[s1][s2]+""));
		
		csv.flush();
		out.close();
	}

}
