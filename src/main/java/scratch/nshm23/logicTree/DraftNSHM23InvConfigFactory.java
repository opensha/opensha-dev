package scratch.nshm23.logicTree;

import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;
import scratch.kevin.nshm23.InversionConfigurationFactory;
import scratch.nshm23.targetMFDs.DraftModelConstraintBuilder;
import scratch.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.SubSeisMoRateReduction;

public class DraftNSHM23InvConfigFactory implements InversionConfigurationFactory {

	@Override
	public FaultSystemRupSet buildRuptureSet(LogicTreeBranch<?> branch) {
		FaultSystemRupSet rupSet = new RuptureSets.U3RupSetConfig(branch.requireValue(FaultModels.class),
				branch.requireValue(ScalingRelationships.class)).build(FaultSysTools.defaultNumThreads());
		
		// create equivalent U3 branch
		U3LogicTreeBranch u3Branch = U3LogicTreeBranch.DEFAULT.copy();
		for (LogicTreeNode node : branch)
			if (node instanceof U3LogicTreeBranchNode<?>)
				u3Branch.setValue((U3LogicTreeBranchNode<?>)node);
		
		System.out.println("Equivalent U3 branch: "+u3Branch);
		rupSet = FaultSystemRupSet.buildFromExisting(rupSet).forU3Branch(u3Branch).build();
		rupSet.addModule(branch);
		return rupSet;
	}

	@Override
	public InversionConfiguration buildInversionConfig(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch,
			CommandLine cmd) {
		double bVal = branch.requireValue(SupraSeisBValue.class).bValue;
		DraftModelConstraintBuilder constrBuilder = new DraftModelConstraintBuilder(rupSet, bVal,
				true, false, true);
		
		constrBuilder.subSeisMoRateReduction(SubSeisMoRateReduction.SYSTEM_AVG_IMPLIED_FROM_SUPRA_B);
		
		SubSectConstraintModel constrModel = branch.requireValue(SubSectConstraintModel.class);
		
		double slipWeight = 1d;
		double paleoWeight = 5;
		double parkWeight = 100;
		double mfdWeight = 10;
		double nuclWeight = constrModel == SubSectConstraintModel.TOT_NUCL_RATE ? 0.5 : 0d;
		double nuclMFDWeight = constrModel == SubSectConstraintModel.NUCL_MFD ? 0.1 : 0d;
		double paleoSmoothWeight = paleoWeight > 0 ? 10000 : 0;
		
		constrBuilder.magDepRelStdDev(M->0.1*Math.pow(10, bVal*0.5*(M-6)));
		
		if (slipWeight > 0d)
			constrBuilder.slipRates().weight(slipWeight);
		
		if (paleoWeight > 0d) {
			constrBuilder.paleoRates().weight(paleoWeight);
			constrBuilder.paleoSlips().weight(paleoWeight);
		}
		
		if (parkWeight > 0d)
			constrBuilder.parkfield().weight(parkWeight);
		
		if (mfdWeight > 0d)
			constrBuilder.supraBValMFDs().weight(mfdWeight);
		
		if (nuclWeight > 0d)
			constrBuilder.sectSupraRates().weight(nuclWeight);
		
		if (nuclMFDWeight > 0d)
			constrBuilder.sectSupraNuclMFDs().weight(nuclMFDWeight);
		
		if (paleoSmoothWeight > 0d)
			constrBuilder.supraPaleoSmooth().weight(paleoSmoothWeight);
		
		IntegerPDF_FunctionSampler sampler = constrBuilder.getSkipBelowMinSampler();
		
		List<InversionConstraint> constraints = constrBuilder.build();
		
		int totThreads = Runtime.getRuntime().availableProcessors();
		int avgThreads = totThreads / 4;
		
		CompletionCriteria completion;
		if (constrModel == SubSectConstraintModel.NUCL_MFD)
			completion = TimeCompletionCriteria.getInHours(5l);
		else
			completion = TimeCompletionCriteria.getInHours(2l);
		
		InversionConfiguration.Builder builder = InversionConfiguration.builder(constraints, completion)
				.threads(totThreads)
				.avgThreads(avgThreads, TimeCompletionCriteria.getInMinutes(5l))
				.perturbation(GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE)
				.nonNegativity(NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN)
				.forCommandLine(cmd).sampler(sampler);
		
		return builder.build();
	}

}
