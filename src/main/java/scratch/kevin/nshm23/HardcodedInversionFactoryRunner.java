package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class HardcodedInversionFactoryRunner {

	public static void main(String[] args) throws IOException {
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		int threads = 16;

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());

//		U3InversionConfigFactory factory = new U3InversionConfigFactory.OriginalCalcParams();
//		dirName += "-u3_orig_params";
//		
//		LogicTreeBranch<U3LogicTreeBranchNode<?>> branch = U3LogicTreeBranch.DEFAULT;

//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
//		dirName += "-nshm23";
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.ClusterSpecific();
		dirName += "-nshm23-cluster_specific";
		
//		LogicTreeBranch<LogicTreeNode> branch = NSHM18_LogicTreeBranch.DEFAULT;
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT;
		
		branch.setValue(RupturePlausibilityModels.UCERF3_REDUCED);
		dirName += "-u3_reduced";
		
		FaultSystemSolution solution = Inversions.run(factory, branch, threads);
		
		System.out.println("Currently loaded modules:");
		IncrementalMagFreqDist target1 = null;
		for (OpenSHA_Module module : solution.getRupSet().getModules(false)) {
			System.out.println(module.getName()+" ("+module.getClass()+")");
			if (module instanceof InversionTargetMFDs)
				target1 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
		}
		for (OpenSHA_Module module : solution.getModules(false))
			System.out.println(module.getName()+" ("+module.getClass()+")");
		System.out.println("Loading all available....");
		IncrementalMagFreqDist target2 = null;
		for (OpenSHA_Module module : solution.getRupSet().getModules(true)) {
			System.out.println(module.getName()+" ("+module.getClass()+")");
			if (module instanceof InversionTargetMFDs)
				target2 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
		}
		for (OpenSHA_Module module : solution.getModules(true))
			System.out.println(module.getName()+" ("+module.getClass()+")");
		
		System.out.println("Processing rupture set manually...");
		factory.getSolutionLogicTreeProcessor().processRupSet(solution.getRupSet(), branch);
		IncrementalMagFreqDist target3 = null;
		for (OpenSHA_Module module : solution.getRupSet().getModules(true)) {
			System.out.println(module.getName()+" ("+module.getClass()+")");
			if (module instanceof InversionTargetMFDs)
				target3 = ((InversionTargetMFDs)module).getTotalOnFaultSupraSeisMFD();
		}
		
		solution.getRupSet().removeModuleInstances(InversionTargetMFDs.class);
		factory.buildInversionConfig(solution.getRupSet(), branch, threads);
		IncrementalMagFreqDist target4 = null;
		target4 = solution.getRupSet().requireModule(InversionTargetMFDs.class).getTotalOnFaultSupraSeisMFD();
		
		System.out.println("Original target, if loaded:\n"+target1);
		System.out.println("Target after loading available:\n"+target2);
		System.out.println("Target after processing:\n"+target3);
		System.out.println("Target after full constraint build:\n"+target4);
		
		File outputDir = new File(parentDir, dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		solution.write(new File(outputDir, "solution.zip"));
	}

}
