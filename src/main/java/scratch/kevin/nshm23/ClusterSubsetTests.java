package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRuptureBuilder;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;

public class ClusterSubsetTests {

	public static void main(String[] args) throws IOException {
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		int threads = 16;

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());

//		U3InversionConfigFactory factory = new U3InversionConfigFactory.OriginalCalcParams();
//		dirName += "-u3_orig_params";
//		
//		LogicTreeBranch<U3LogicTreeBranchNode<?>> branch = U3LogicTreeBranch.DEFAULT;

		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		dirName += "-nshm23";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.ClusterSpecific();
//		dirName += "-nshm23-cluster_specific";
		
//		LogicTreeBranch<LogicTreeNode> branch = NSHM18_LogicTreeBranch.DEFAULT;
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT;
		
		RupturePlausibilityModels plausibilityBranch = RupturePlausibilityModels.UCERF3_REDUCED;
		branch.setValue(plausibilityBranch);
		dirName += "-u3_reduced";
		
		dirName += "-keddie_debug";
		int targetSectID = 999;
		
		FaultSystemRupSet fullRupSet = factory.buildRuptureSet(branch, threads);
		
		List<ConnectivityCluster> clusters = ConnectivityCluster.build(fullRupSet);
		ConnectivityCluster matchingCluster = null;
		for (ConnectivityCluster cluster : clusters) {
			if (cluster.containsSect(targetSectID)) {
				matchingCluster = cluster;
				break;
			}
		}
		Preconditions.checkNotNull(matchingCluster);
		System.out.println("Found matching cluster: "+matchingCluster);
		
		// build the cluster as a subset
		FaultSystemRupSet rupSet = fullRupSet.getForSectionSubSet(matchingCluster.getSectIDs());
		
		FaultSystemSolution solution = Inversions.run(rupSet, factory, branch, threads, null);
		
		File outputDir = new File(parentDir, dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		solution.write(new File(outputDir, "subset_solution.zip"));
		
		// now build it individually
		RupSetConfig rsConfig = plausibilityBranch.getConfig(rupSet.getFaultSectionDataList(),
				branch.requireValue(RupSetScalingRelationship.class));
		
		rupSet = rsConfig.build(threads);
		
		PlausibilityConfiguration plausibility = rupSet.getModule(PlausibilityConfiguration.class);
		RupSetScalingRelationship scale = branch.requireValue(RupSetScalingRelationship.class);
		
		ClusterRuptures cRups = rupSet.getModule(ClusterRuptures.class);
		if (cRups == null) {
			rupSet = FaultSystemRupSet.builder(rupSet.getFaultSectionDataList(), rupSet.getSectionIndicesForAllRups())
					.forScalingRelationship(scale).build();
			if (plausibility != null)
				rupSet.addModule(plausibility);
			rupSet.addModule(ClusterRuptures.singleStranged(rupSet));
		} else {
			rupSet = ClusterRuptureBuilder.buildClusterRupSet(scale, rupSet.getFaultSectionDataList(), plausibility, cRups.getAll());
		}
		
		SlipAlongRuptureModels slipAlong = branch.requireValue(SlipAlongRuptureModels.class);
		rupSet.addModule(slipAlong.getModel());
		
		// add other modules
		factory.getSolutionLogicTreeProcessor().processRupSet(rupSet, branch);
		
		// remove paleo
		rupSet.removeModuleInstances(PaleoseismicConstraintData.class);
		
		solution = Inversions.run(rupSet, factory, branch, threads, null);
		
		solution.write(new File(outputDir, "rebuilt_solution.zip"));
	}

}
