package scratch.kevin.quantum;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.RuptureSubSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;

public class SubsetReportGen {

	public static void main(String[] args) throws IOException {
		File inputSol = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSol);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		ConnectivityClusters clusters = rupSet.getModule(ConnectivityClusters.class);
		if (clusters == null)
			clusters = ConnectivityClusters.build(rupSet);
		List<ConnectivityCluster> clustersList = clusters.get();
		Collections.sort(clustersList);
		ConnectivityCluster cluster = clusters.get(clusters.size()-2);
		System.out.println("Cluster has "+cluster.getNumSections()+" sects");
		FaultSystemRupSet rupSubSet = rupSet.getForSectionSubSet(cluster.getSectIDs());
		RuptureSubSetMappings mappings = rupSubSet.requireModule(RuptureSubSetMappings.class);
		
		double[] ratesSubSet = new double[rupSubSet.getNumRuptures()];
		for (int r=0; r<ratesSubSet.length; r++) {
			ratesSubSet[r] = sol.getRateForRup(mappings.getOrigRupID(r));
		}
		FaultSystemSolution solSubset = new FaultSystemSolution(rupSubSet, ratesSubSet);
		solSubset.write(new File("/home/kevin/markdown/inversions/2024_11_06-quantum_test_rup_set_problem-213_sects/orig_sol.zip"));
		
	}

}
