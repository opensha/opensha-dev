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
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_PaleoUncertainties;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.U3InversionConfigFactory;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;

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
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		dirName += "-nshm23";
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.ForceNewPaleo();
//		dirName += "-nshm23-new_paleo";
		
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		
//		LogicTreeBranch<U3LogicTreeBranchNode<?>> branch = U3LogicTreeBranch.DEFAULT;
//		LogicTreeBranch<LogicTreeNode> branch = NSHM18_LogicTreeBranch.DEFAULT;
//		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT; dirName += "-u3";
//		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT;
		
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_U3_HybridLogicTreeBranch.levels;
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
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = NSHM23_LogicTreeBranch.levels;
		dirName += "-single_state";
		levels = new ArrayList<>(levels);
		levels.add(NSHM23_LogicTreeBranch.SINGLE_STATES);
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT)
			branch.setValue(node);
		
		branch.setValue(NSHM23_SingleStates.UT);
//		branch.setValue(NSHM23_SingleStates.NM);
		
//		branch.setValue(RupturePlausibilityModels.UCERF3_REDUCED);
		
//		branch.setValue(NSHM23_FaultModels.NSHM23_v2);
		
//		branch.setValue(RupturePlausibilityModels.AZIMUTHAL_REDUCED);
		
//		branch.setValue(NSHM23_DeformationModels.ZENG);
//		branch.setValue(NSHM23_DeformationModels.EVANS);
//		branch.setValue(NSHM23_DeformationModels.SHEN_BIRD);
		branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
		
//		branch.setValue(ScalingRelationships.MEAN_UCERF3);
		
		branch.setValue(NSHM23_SegmentationModels.CLASSIC);
		
		branch.setValue(NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD);
		
		branch.setValue(SupraSeisBValues.B_0p0);
		
//		branch.setValue(NSHM23_PaleoUncertainties.OVER_FIT);
		
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
				if (node.getNodeWeight(branch) > 0d) {
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
		
		File outputDir = new File(parentDir, dirName);
		System.out.println("Will save results in: "+outputDir.getAbsolutePath());
		
		FaultSystemSolution solution = Inversions.run(factory, branch, threads);
		
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
