package scratch.kevin.ucerf3.eal.branches;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

public enum U3_EAL_Vs30Model implements LogicTreeBranchNode<U3_EAL_Vs30Model> {
	WILLS_2015("Wills (2015)", "Wills2015", 0.5),
	WALD_ALLEN("Wald & Allen (2007,2008)", "WaldAllen", 0.5);

	private String name;
	private String shortName;
	private double weight;
	
	private U3_EAL_Vs30Model(String name, String shortName, double weight) {
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
		return getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "Vs30 Model";
	}

	@Override
	public String getShortBranchLevelName() {
		return "Vs30";
	}

}
