package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoProbabilityModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.U3_UncertAddDeformationModels;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;

public class U3SerializationDebug {

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT.copy();
		
		// CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/results/FM3_1_CoulombRupSet_U3_ABM_Shaw09Mod_DsrUni_SupraB0.
		// 0_TotNuclRate_NoRed_UnderFitPaleo_Classic_ThreshAvgIterRelGR
		branch.setValue(U3_UncertAddDeformationModels.U3_ABM);
		branch.setValue(ScalingRelationships.SHAW_2009_MOD);
		branch.setValue(SlipAlongRuptureModels.UNIFORM);
		branch.setValue(SupraSeisBValues.B_0p0);
		branch.setValue(NSHM23_PaleoUncertainties.UNDER_FIT);
		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		
		InversionConfigurationFactory factory = new NSHM23_InvConfigFactory.FullSysInv();
		factory.setCacheDir(new File("/data/kevin/nshm23/rup_sets/cache"));
		
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 20);
		InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, 20);
		config = InversionConfiguration.builder(config)
				.avgThreads(4, TimeCompletionCriteria.getInSeconds(20))
				.completion(TimeCompletionCriteria.getInMinutes(1))
				.build();
		
		FaultSystemSolution sol = Inversions.run(rupSet, config);
		
		sol.write(new File("/tmp/test_sol.zip"));
	}

}
