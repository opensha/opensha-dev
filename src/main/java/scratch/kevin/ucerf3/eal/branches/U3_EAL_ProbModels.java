package scratch.kevin.ucerf3.eal.branches;

import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;

public enum U3_EAL_ProbModels implements LogicTreeBranchNode<U3_EAL_ProbModels> {
	POISSON("Poisson", "POISSON", null, FaultSystemSolutionERF.PREF_BLEND_POISSON_WEIGHT),
	COV_LOW("Low COV Values", "LOW_VALUES", null, FaultSystemSolutionERF.PREF_BLEND_COV_LOW_WEIGHT),
	COV_MID("Mid COV Values", "MID_VALUES", null, FaultSystemSolutionERF.PREF_BLEND_COV_MID_WEIGHT),
	COV_HIGH("High COV Values", "HIGH_VALUES", null, FaultSystemSolutionERF.PREF_BLEND_COV_HIGH_WEIGHT);
	
	private String name;
	private String shortName;
	private MagDependentAperiodicityOptions cov;
	private double weight;
	
	private U3_EAL_ProbModels(String name, String shortName, MagDependentAperiodicityOptions cov, double weight) {
		this.name = name;
		this.shortName = shortName;
		this.cov = cov;
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
		return "Prob"+getShortName();
	}

	@Override
	public String getBranchLevelName() {
		return "ERF Probability Model";
	}

	@Override
	public String getShortBranchLevelName() {
		return "Prob";
	}

}
