package scratch.kevin.mfdInversion;

import java.io.File;
import java.io.IOException;

import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PosteriorSectionBValueDistributions;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;

import com.google.common.base.Preconditions;

public class PaleoBValueEstimatorRunner {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/paleo_b_value");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(NSHM23_DeformationModels.AVERAGE);
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		System.out.println(branch);
		
		NSHM23_InvConfigFactory.APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT = false;
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
		ContinuousDistribution priorDist = UniformContinuousDistribution.of(0d, 1d);
//		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 11);
		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 21);
		
		PaleoBValueEstimator estimator = new PaleoBValueEstimator(priorDist, bVals, factory);
		
		FaultSystemRupSet rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		PosteriorSectionBValueDistributions result = estimator.calculate(rs, branch, true);
		
		PaleoBValueEstimator.plotSectDistributions(outputDir, rs, result);
		estimator.plotAlongStrike(outputDir, "along_strike", rs, result, rs.requireModule(NamedFaults.class));
	}
	
}