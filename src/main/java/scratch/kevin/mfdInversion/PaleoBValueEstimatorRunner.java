package scratch.kevin.mfdInversion;

import java.io.File;
import java.io.IOException;

import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.PosteriorSectionBValueDistributions;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator.BackgroundDistribution;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator.PosteriorCombination;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator.WeightingType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class PaleoBValueEstimatorRunner {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/nshm23/paleo_b_value");
		Preconditions.checkState(mainDir.exists() || mainDir.mkdir());
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(NSHM23_DeformationModels.AVERAGE);
		branch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		branch.setValue(NSHM23_SegmentationModels.AVERAGE);
		System.out.println(branch);
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
		String distAdd = "";
		ContinuousDistribution priorDist = UniformContinuousDistribution.of(0d, 1d);
//		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 11);
		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 21);
//		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0.05d, 1d, 20);
		
//		String distAdd = "-widePrior";
//		ContinuousDistribution priorDist = UniformContinuousDistribution.of(0d, 2d);
////		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 1d, 11);
//		EvenlyDiscretizedFunc bVals = new EvenlyDiscretizedFunc(0d, 2d, 41);
		
		PaleoBValueEstimator estimator = new PaleoBValueEstimator(priorDist, bVals, factory);
		
//		estimator.setWeightingType(WeightingType.CONNECTED_FRACTION);
//		estimator.setPosteriorCombination(PosteriorCombination.WEIGHTED_SITE_LIKELIHOODS);
//		estimator.setPosteriorCombination(PosteriorCombination.WEIGHTED_SITE_DISTRIBUTIONS);
//		estimator.setBackgroundDistribution(BackgroundDistribution.NAMED_FAULT_AVERAGE);
//		estimator.setBackgroundDistribution(BackgroundDistribution.PALEO_AVERAGE);
		
//		estimator.setStudentsTLogLikelihood(4);
//		estimator.setLikelihoodInRecurrence(true);
//		String dirName ="students_t_4_ri";
		
//		estimator.setStudentsTLogLikelihood(5);
//		estimator.setLikelihoodInRecurrence(true);
//		String dirName = "students_t_5_ri";
		
		estimator.setStudentsTLogLikelihood(5);
		estimator.setLikelihoodInRecurrence(true);
		estimator.setPosteriorCombination(PosteriorCombination.WEIGHTED_SITE_DISTRIBUTIONS);
		String dirName = "students_t_5_ri_dist_avg";
		
//		estimator.setStudentsTLogLikelihood(4);
//		estimator.setLikelihoodInRecurrence(false);
//		File outputDir = new File(mainDir, "students_t_4_rate");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
//		estimator.setGaussianLogLikelihood();
//		estimator.setLikelihoodInRecurrence(false);
//		File outputDir = new File(mainDir, "gaussian_rate");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
//		estimator.setGaussianLogLikelihood();
//		estimator.setLikelihoodInRecurrence(true);
//		String dirName = "gaussian_ri";
		
//		estimator.setGaussianLogLikelihood();
//		estimator.setLikelihoodInRecurrence(true);
//		estimator.setPosteriorCombination(PosteriorCombination.WEIGHTED_SITE_DISTRIBUTIONS);
//		String dirName = "gaussian_ri_dist_avg";
		
//		estimator.setStudentsTLogLikelihood(5);
//		estimator.setLikelihoodInRecurrence(true);
//		branch.setValue(NSHM23_DeformationModels.EVANS);
//		branch.setValue(NSHM23_PaleoUncertainties.UNDER_FIT);
//		branch.setValue(NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD);
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
//		String dirName = "quick_test";
		
		File outputDir = new File(mainDir, dirName+distAdd);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		estimator.setTargetMFDCalc(NSHM23_InvConfigFactory.getMFDCalcFunc());
		
		FaultSystemRupSet rs = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
//		for (SectMappedUncertainDataConstraint data : rs.requireModule(PaleoseismicConstraintData.class).getPaleoRateConstraints()) {
//			FaultSection sect = rs.getFaultSectionData(data.sectionIndex);
//			System.out.println(data.name+"\t\t"+sect.getSectionName()+"\t\trake="+(int)sect.getAveRake());
//		}
//		System.exit(0);
		
		PosteriorSectionBValueDistributions result = estimator.calculate(rs, branch, true);
		
		NamedFaults named = rs.requireModule(NamedFaults.class);
		PaleoBValueEstimator.plotSectDistributions(outputDir, rs, result);
		PaleoBValueEstimator.plotFaultSiteAverages(outputDir, "combined", rs, result, named);
		estimator.plotAlongStrike(outputDir, "along_strike", rs, result, named);
		
		result.writeToFile(new File(outputDir, "posteriors.json"));
	}
	
}