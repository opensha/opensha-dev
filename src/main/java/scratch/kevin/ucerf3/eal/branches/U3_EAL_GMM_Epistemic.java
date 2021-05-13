package scratch.kevin.ucerf3.eal.branches;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

public enum U3_EAL_GMM_Epistemic implements LogicTreeBranchNode<U3_EAL_GMM_Epistemic> {
	UPPER("Upper", "UPPER", 0.185),
	NONE("None", "NONE", 0.630),
	LOWER("Lower", "LOWER", 0.185);
	
	private String name;
	private String shortName;
	private double weight;
	
	private U3_EAL_GMM_Epistemic(String name, String shortName, double weight) {
		this.name = name;
		this.shortName = shortName;
		this.weight = weight;
	}

	@Override
	public String getShortName() {
		return shortName;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public double getRelativeWeight(InversionModels im) {
		return weight;
	}

	@Override
	public String encodeChoiceString() {
		return "Epi"+getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "GMM Additional Epistemic Uncertainty";
	}

}
