package scratch.kevin.segmentation.normalFault;

import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Transition;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationRuptureExemptions.Policy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels.MaxLengthSegmentationModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupsThroughCreepingSectBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;

import scratch.kevin.segmentation.normalFault.NormalFaultGeometricSegmentationProb.CoefficientSet;

/**
 * Dev-only experimental replacement for the NSHM23 segmentation branching level. Both choices
 * use a 10-degree bend threshold and unrestricted perpendicular deviation. Normal
 * uses OpenSHA's rake classification; geologic rakes are supplied by the GEOLOGIC
 * deformation model. Other/mixed mechanisms use MID jump rules without hardcoded
 * Wasatch segmentation. Retains MID's global rupture-length limit.
 *
 * <p>Replace the original segmentation level with an enum level for this class;
 * do not add it alongside the original segmentation level. Probabilities are
 * averaged across crossing ruptures in both the constraint and target-MFD paths.
 * The two preliminary fits have equal exploratory weights, not calibrated weights. Internal bends exempt
 * tied minimum-size ruptures on a side only when no eligible non-crossing alternative exists entirely within
 * the same parent fault; larger crossings remain constrained. This enforcement policy is shared by the
 * inversion constraint and MFD adjustment.
 */
@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public enum NSHM23_NormalFaultSegmentationModels
implements SegmentationModelBranchNode<Transition>, RupsThroughCreepingSectBranchNode {
	MANUAL_60K("Normal Fault Geometry, Manual 1:60k", "NormalGeomManual60k", CoefficientSet.MANUAL_60K),
	AUTOMATED_5KM("Normal Fault Geometry, Automated 5 km", "NormalGeomAuto5km", CoefficientSet.AUTOMATED_5KM);

	private final String name;
	private final String shortName;
	private final CoefficientSet coefficients;

	private NSHM23_NormalFaultSegmentationModels(String name, String shortName, CoefficientSet coefficients) {
		this.name = name;
		this.shortName = shortName;
		this.coefficients = coefficients;
	}

	public CoefficientSet getCoefficientSet() { return coefficients; }
	public double getMinimumBendAngleDegrees() { return 10d; }
	public NSHM23_SegmentationModels getBackgroundChoice() { return NSHM23_SegmentationModels.MID; }

	@Override
	public Policy getRuptureExemptionPolicy() { return Policy.MINIMUM_INTERNAL_IF_NO_ALTERNATIVE; }

	@Override
	public NormalFaultGeometricSegmentationModel getModel(
			FaultSystemRupSet rupSet, LogicTreeBranch<?> branch) {
		NSHM23_SegmentationModels background = getBackgroundChoice();
		JumpProbabilityCalc jumps = NSHM23_SegmentationModels.buildModel(rupSet, background.getShawR0(),
				background.getShawShift(), 1d, background.getCreepingSectPassthrough(), Double.NaN, false);
		return new NormalFaultGeometricSegmentationModel(coefficients, getMinimumBendAngleDegrees(), jumps);
	}

	@Override
	public BinaryRuptureProbabilityCalc getExclusionModel(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch) {
		BinaryRuptureProbabilityCalc transitions = SegmentationModelBranchNode.super.getExclusionModel(rupSet, branch);
		double maxLength = getBackgroundChoice().getMaxRuptureLength();
		if (!Double.isFinite(maxLength) || maxLength <= 0d) return transitions;
		BinaryRuptureProbabilityCalc length = new MaxLengthSegmentationModel(maxLength);
		return transitions == null ? length : new RuptureProbabilityCalc.LogicalAnd(transitions, length);
	}

	@Override public boolean isExcludeRupturesThroughCreepingSect() { return false; }
	@Override public String getName() { return name; }
	@Override public String getShortName() { return shortName; }
	@Override public String getFilePrefix() { return shortName; }
	@Override public double getNodeWeight(LogicTreeBranch<?> branch) { return 0.5; }
}
