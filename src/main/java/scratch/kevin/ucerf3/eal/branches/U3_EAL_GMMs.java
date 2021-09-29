package scratch.kevin.ucerf3.eal.branches;

import org.opensha.sha.imr.AttenRelRef;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;

public enum U3_EAL_GMMs implements U3LogicTreeBranchNode<U3_EAL_GMMs> {
	ASK_2014(AttenRelRef.ASK_2014, 0.22),
	BSSA_2014(AttenRelRef.BSSA_2014, 0.22),
	CB_2014(AttenRelRef.CB_2014, 0.22),
	CY_2014(AttenRelRef.CY_2014, 0.22),
	IDRISS_2014(AttenRelRef.IDRISS_2014, 0.12);
	
	private AttenRelRef ref;
	private double weight;
	private U3_EAL_GMMs(AttenRelRef ref, double weight) {
		this.ref = ref;
		this.weight = weight;
	}

	@Override
	public String getShortName() {
		return name();
	}

	@Override
	public String getName() {
		return ref.getName();
	}

	@Override
	public double getRelativeWeight(InversionModels im) {
		return weight;
	}

	@Override
	public String encodeChoiceString() {
		return getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "Ground Motion Model";
	}

	@Override
	public String getShortBranchLevelName() {
		return "GMM";
	}

}
