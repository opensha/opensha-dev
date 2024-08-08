package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoDataLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_WasatchSegmentationData;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.U3_UncertAddDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.random.RandomBValSampler;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.random.RandomSegModelSampler;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.prior2018.NSHM18_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SparseGutenbergRichterSolver;
import org.opensha.sha.magdist.SparseGutenbergRichterSolver.SpreadingMethod;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.U3InversionConfigFactory;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;
import scratch.kevin.nshm23.dmCovarianceTests.DefModSamplingEnabledInvConfig;
import scratch.kevin.nshm23.dmCovarianceTests.RandomDefModSampleLevel;
import scratch.kevin.nshm23.dmCovarianceTests.SectionCovarianceSampler;

public class HardcodedInversionFactoryRunner {

	public static void main(String[] args) throws IOException {
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		int threads = 16;

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());

//		U3InversionConfigFactory factory = new U3InversionConfigFactory.OriginalCalcParams();
//		dirName += "-u3_orig_params";
//		
//		LogicTreeBranch<U3LogicTreeBranchNode<?>> branch = U3LogicTreeBranch.DEFAULT;

//		InversionConfigurationFactory factory = new U3InversionConfigFactory.ForceNewPaleo();
//		dirName += "-u3-new_paleo";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.FullSysInv();
//		dirName += "-nshm23-full_sys";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.NSHM18_UseU3Paleo();
//		dirName += "-nshm18-u3_paleo";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.FullSysInv() {
//
//			@Override
//			public InversionConfiguration buildInversionConfig(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch,
//					int threads) {
//				InversionConfiguration config = super.buildInversionConfig(rupSet, branch, threads);
//				// disable reweighting
//				config = InversionConfiguration.builder(config).reweight(null).build();
//				// override weights
//				// these are from the 2022_12_08 U3 ref branch runs, top is even fit paleo, bottom under fit
//				for (InversionConstraint constr : config.getConstraints()) {
//					if (constr instanceof SlipRateInversionConstraint)
////						constr.setWeight(0.9081427);
//						constr.setWeight(0.92429525);
//					else if (constr instanceof PaleoRateInversionConstraint)
////						constr.setWeight(8.173598);
//						constr.setWeight(0.05);
//					else if (constr instanceof ParkfieldInversionConstraint)
////						constr.setWeight(15.864809);
//						constr.setWeight(17.998356);
//					else if (constr instanceof MFDInversionConstraint)
////						constr.setWeight(7.9695096);
//						constr.setWeight(6.446304);
//					else if (constr instanceof SectionTotalRateConstraint)
////						constr.setWeight(0.57912546);
//						constr.setWeight(0.67924696);
//				}
//				
//				return config;
//			}
//			
//		};
//		dirName += "-nshm23-full_sys-prev_weights";
//		NSHM23_PaleoDataLoader.INCLUDE_U3_PALEO_SLIP = false;
//		dirName += "-no_paleo_slip";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
//		dirName += "-nshm23";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.PaleoSlipInequality();
//		dirName += "-nshm23-paleo_slip_ineq";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.ForceNewPaleo();
//		dirName += "-nshm23-new_paleo";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.MatchFullBA();
//		dirName += "-nshm23-nucl_match_ba";
//		NSHM23_InvConfigFactory factory = new DefModSamplingEnabledInvConfig.ConnDistB0p5MidSegCorrCapSigma();
//		dirName += "-nshm23-dm_sample_cap_sigma";
		
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		dirName += "-prvi25";
		PRVI25_InvConfigFactory.SUB_SECT_DDW_FRACT = 0.25; dirName += "-quarter_len_sub_sects";
		
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
		boolean writeRS = true;
		
//		LogicTreeBranch<U3LogicTreeBranchNode<?>> branch = U3LogicTreeBranch.DEFAULT;
//		LogicTreeBranch<LogicTreeNode> branch = NSHM18_LogicTreeBranch.DEFAULT; dirName += "-2018_inputs";
//		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT; dirName += "-u3";
//		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT;
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTreeBranch.DEFAULT_CRUSTAL_ON_FAULT;
		branch = branch.copy();
		
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		
		// seg/b sampling
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levelsOnFault;
//		levels = new ArrayList<>(levels);
//		boolean randB = true;
//		boolean randSeg = true;
//		int origSize = levels.size();
//		List<LogicTreeNode> values = new ArrayList<>();
//		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
//			values.add(node);
//		for (int i=levels.size(); --i>=0;) {
//			if (randB && SupraSeisBValues.class.isAssignableFrom(levels.get(i).getType())) {
//				levels.remove(i);
//				values.remove(i);
//			}
//			if (randSeg && SegmentationModelBranchNode.class.isAssignableFrom(levels.get(i).getType())) {
//				levels.remove(i);
//				values.remove(i);
//			}
//		}
//		Preconditions.checkState(levels.size() < origSize);
//		if (randB) {
//			dirName += "-randB";
//			RandomBValSampler.Level level = new RandomBValSampler.Level(1);
//			levels.add(level);
//			values.add(level.getNodes().get(0));
//		}
//		if (randSeg) {
//			dirName += "-randSeg";
//			RandomSegModelSampler.Level level = new RandomSegModelSampler.Level(1);
//			levels.add(level);
//			values.add(level.getNodes().get(0));
//		}
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels, values);
//		branch.setValue(NSHM23_DeformationModels.EVANS);
		
		// DM sampling
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levelsOnFault;
//		levels = new ArrayList<>(levels);
//		List<LogicTreeNode> values = new ArrayList<>();
//		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
//			values.add(node);
//		levels.add(new RandomDefModSampleLevel(1));
//		values.add(levels.get(levels.size()-1).getNodes().get(0));
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels, values);
//		branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		
		// all branch averaged
//		dirName += "-all_ba";
//		branch.setValue(NSHM23_DeformationModels.AVERAGE);
//		branch.setValue(NSHM23_ScalingRelationships.AVERAGE);
//		branch.setValue(SupraSeisBValues.AVERAGE);
//		branch.setValue(NSHM23_SegmentationModels.AVERAGE);
//		branch.setValue(NSHM23_PaleoUncertainties.AVERAGE);
		
//		NSHM23_InvConfigFactory.MFD_MIN_FRACT_UNCERT = 0.1;
//		dirName += "-mfd_min_uncert_0.1";
//		NSHM23_InvConfigFactory.MFD_MIN_FRACT_UNCERT = 0.2;
//		dirName += "-mfd_min_uncert_0.2";
		
//		SparseGutenbergRichterSolver.METHOD_DEFAULT = SpreadingMethod.NEAREST;
//		dirName += "-sparse_gr_nearest";
		
//		dirName += "-quick_phase_in";
		
//		NSHM23_InvConfigFactory.PARKFIELD_INITIAL = false;
//		dirName += "-no_parkfield_initial";
		
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_U3_HybridLogicTreeBranch.levels;
////		dirName += "-u3_ref";
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
//		for (LogicTreeNode node : NSHM23_U3_HybridLogicTreeBranch.DEFAULT)
//			branch.setValue(node);
////		branch.setValue(U3_UncertAddDeformationModels.U3_GEOL);
////		branch.setValue(ScalingRelationships.HANKS_BAKUN_08);
////		branch.setValue(SupraSeisBValues.B_0p0);
////		branch.setValue(NSHM23_PaleoUncertainties.OVER_FIT);
////		branch.setValue(NSHM23_SegmentationModels.MID);
//		branch.setValue(U3_UncertAddDeformationModels.U3_ZENG);
//		branch.setValue(ScalingRelationships.ELLB_SQRT_LENGTH);
//		branch.setValue(SupraSeisBValues.B_0p0);
//		branch.setValue(NSHM23_PaleoUncertainties.OVER_FIT);
//		branch.setValue(NSHM23_SegmentationModels.NONE);
		
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM18_LogicTreeBranch.levelsNewScale;
//		dirName += "-nshm18_dms";
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
//		for (LogicTreeNode node : NSHM18_LogicTreeBranch.DEFAULT_NEW_SCALE)
//			branch.setValue(node);
//		branch.setValue(NSHM18_FaultModels.NSHM18_WUS_PlusU3_FM_3p1);
//		branch.setValue(NSHM18_DeformationModels.BRANCH_AVERAGED);
		
//		levels = new ArrayList<>(levels);
//		int origSize = levels.size();
//		for (int i=levels.size(); --i>=0;)
//			if (levels.get(i).getType().isAssignableFrom(ScalingRelationships.class))
//				levels.remove(i);
//		Preconditions.checkState(levels.size() < origSize);
//		levels.add(NSHM23_LogicTreeBranch.SCALE);
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
//		for (LogicTreeNode node : NSHM23_U3_HybridLogicTreeBranch.DEFAULT)
//			if (!(node instanceof ScalingRelationships))
//				branch.setValue(node);
		
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levelsOnFault;
		
//		dirName += "-single_state";
//		levels = new ArrayList<>(levels);
//		levels.add(NSHM23_LogicTreeBranch.SINGLE_STATES);
		
//		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
//		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
//			branch.setValue(node);
		
//		branch.setValue(NSHM23_SingleStates.UT);
//		branch.setValue(NSHM23_SingleStates.NM);
//		branch.setValue(NSHM23_SingleStates.CA);
		
//		NSHM23_WasatchSegmentationData.APPLY_TO_ALL_JUMPS_FROM_LOC = false;
//		dirName += "-prevWasatchSeg";
		
//		branch.setValue(RupturePlausibilityModels.UCERF3);
//		branch.setValue(RupturePlausibilityModels.UCERF3_REDUCED);
		
//		branch.setValue(NSHM23_FaultModels.NSHM23_v2);
		
//		branch.setValue(RupturePlausibilityModels.AZIMUTHAL_REDUCED);
		
//		branch.setValue(NSHM23_DeformationModels.ZENG);
//		branch.setValue(NSHM23_DeformationModels.EVANS);
//		branch.setValue(NSHM23_DeformationModels.SHEN_BIRD);
//		branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
//		branch.setValue(NSHM23_DeformationModels.AVERAGE);
//		branch.setValue(NSHM23_DeformationModels.MEDIAN);
		
//		branch.setValue(ScalingRelationships.MEAN_UCERF3);
		
//		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		
//		branch.setValue(NSHM23_ScalingRelationships.LOGA_C4p1);
//		branch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		
//		branch.setValue(SupraSeisBValues.B_0p5);
		
//		branch.setValue(NSHM23_PaleoUncertainties.EVEN_FIT);
//		branch.setValue(NSHM23_PaleoUncertainties.UNDER_FIT);
		
////		branch.setValue(NSHM23_ScalingRelationships.LOGA_C4p1);
////		branch.setValue(NSHM23_ScalingRelationships.LOGA_C4p2);
////		branch.setValue(NSHM23_ScalingRelationships.WIDTH_LIMITED);
////		branch.setValue(NSHM23_ScalingRelationships.LOGA_C4p1_SQRT_LEN);
////		branch.setValue(NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN);
//		branch.setValue(NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD);
//		dirName += "-scale"+branch.getValue(NSHM23_ScalingRelationships.class).getShortName();
		
		boolean first = true;
		for (int i=0; i<branch.size(); i++) {
			LogicTreeNode node = branch.getValue(i);
			if (node != null) {
				if (!(node instanceof RandomlySampledNode) && node.getNodeWeight(branch) > 0d) {
					// only include its name if there are other alternatives (unless we have chosen a zero-weight option)
					boolean hasOthers = false;
					for (LogicTreeNode oNode : branch.getLevel(i).getNodes()) {
						if (oNode != node && oNode.getNodeWeight(branch) > 0d) {
							hasOthers = true;
							break;
						}
					}
					if (!hasOthers)
						continue;
				}
				if (first) {
					dirName += "-";
					first = false;
				} else {
					dirName += "_";
				}
				dirName += node.getFilePrefix();
			}
		}
		
		System.out.println("Logic tree branch: "+branch);
		
		File outputDir = new File(parentDir, dirName);
		System.out.println("Will save results in: "+outputDir.getAbsolutePath());
		
		FaultSystemSolution solution;
		if (writeRS) {
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, threads);
			
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			
			rupSet.write(new File(outputDir, "rupSet.zip"));
			
			solution = Inversions.run(rupSet, factory, branch, threads, null);
		} else {
			solution = Inversions.run(factory, branch, threads);
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		solution.write(new File(outputDir, "solution.zip"));
		
//		System.out.println("Currently loaded modules:");
//		IncrementalMagFreqDist target1 = null;
//		for (OpenSHA_Module module : solution.getRupSet().getModules(false)) {
//			System.out.println(module.getName()+" ("+module.getClass()+")");
//			if (module instanceof InversionTargetMFDs)
//				target1 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
//		}
//		for (OpenSHA_Module module : solution.getModules(false))
//			System.out.println(module.getName()+" ("+module.getClass()+")");
//		System.out.println("Loading all available....");
//		IncrementalMagFreqDist target2 = null;
//		for (OpenSHA_Module module : solution.getRupSet().getModules(true)) {
//			System.out.println(module.getName()+" ("+module.getClass()+")");
//			if (module instanceof InversionTargetMFDs)
//				target2 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
//		}
//		for (OpenSHA_Module module : solution.getModules(true))
//			System.out.println(module.getName()+" ("+module.getClass()+")");
//		
//		System.out.println("Processing rupture set manually...");
//		factory.getSolutionLogicTreeProcessor().processRupSet(solution.getRupSet(), branch);
//		IncrementalMagFreqDist target3 = null;
//		for (OpenSHA_Module module : solution.getRupSet().getModules(true)) {
//			System.out.println(module.getName()+" ("+module.getClass()+")");
//			if (module instanceof InversionTargetMFDs)
//				target3 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
//		}
//		
//		solution.getRupSet().removeModuleInstances(InversionTargetMFDs.class);
//		factory.buildInversionConfig(solution.getRupSet(), branch, threads);
//		IncrementalMagFreqDist target4 = null;
//		target4 = solution.getRupSet().requireModule(InversionTargetMFDs.class).getTotalOnFaultSupraSeisMFD();
//		
//		System.out.println("Original target, if loaded:\n"+target1);
//		System.out.println("Target after loading available:\n"+target2);
//		System.out.println("Target after processing:\n"+target3);
//		System.out.println("Target after full constraint build:\n"+target4);
	}

}
