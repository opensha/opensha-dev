package scratch.kevin.segmentation;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupsThroughCreepingSectBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;

public enum TestRelativeSegModels implements SegmentationModelBranchNode, RupsThroughCreepingSectBranchNode {
	RELATIVE_CFF_MIDDLE("Relative (Middle)", NSHM23_SegmentationModels.MID);
	
	private String name;
	private NSHM23_SegmentationModels segModel;

	private TestRelativeSegModels(String name, NSHM23_SegmentationModels segModel) {
		this.name = name;
		this.segModel = segModel;
		
	}

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return 1;
	}

	@Override
	public String getFilePrefix() {
		return name();
	}

	@Override
	public String getShortName() {
		return name;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public boolean isExcludeRupturesThroughCreepingSect() {
		return segModel.isExcludeRupturesThroughCreepingSect();
	}

	@Override
	public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch) {
		JumpProbabilityCalc upstreamModel = NSHM23_SegmentationModels.buildModel(rupSet, Double.NaN, Double.NaN,
				segModel.getWasatchPassthrough(), segModel.getCreepingSectPassthrough(), Double.NaN, false);
		return null;
	}

}
