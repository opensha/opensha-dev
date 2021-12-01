package scratch.kevin.nshm23;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;

public class U3InversionConfigFactory implements InversionConfigurationFactory {

	@Override
	public FaultSystemRupSet buildRuptureSet(LogicTreeBranch<?> branch, int threads) {
		FaultSystemRupSet rupSet = new RuptureSets.U3RupSetConfig(branch.requireValue(FaultModels.class),
				branch.requireValue(ScalingRelationships.class)).build(threads);
		
		U3LogicTreeBranch u3Branch;
		if (branch instanceof U3LogicTreeBranch) {
			u3Branch = (U3LogicTreeBranch)branch;
		} else {
			List<U3LogicTreeBranchNode<?>> nodes = new ArrayList<>();
			for (LogicTreeNode node : branch) {
				Preconditions.checkState(node instanceof U3LogicTreeBranchNode<?>);
				nodes.add((U3LogicTreeBranchNode<?>) node);
			}
			u3Branch = U3LogicTreeBranch.fromValues(nodes);
		}
		
		return FaultSystemRupSet.buildFromExisting(rupSet).forU3Branch(u3Branch).build();
	}

	@Override
	public InversionConfiguration buildInversionConfig(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch,
			CommandLine cmd, int threads) {
		List<InversionConstraint> constraints;
		try {
			constraints = InversionsCLI.getU3Constraints(rupSet);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		int avgThreads = threads / 4;
		
		return InversionConfiguration.builder(constraints, TimeCompletionCriteria.getInHours(5l))
				.threads(threads)
				.avgThreads(avgThreads, TimeCompletionCriteria.getInMinutes(5l))
				.perturbation(GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE)
				.nonNegativity(NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN)
				.forCommandLine(cmd).build();
	}

	@Override
	public SolutionLogicTree initSolutionLogicTree(LogicTree<?> logicTree) {
		return new SolutionLogicTree.UCERF3(logicTree);
	}

}
