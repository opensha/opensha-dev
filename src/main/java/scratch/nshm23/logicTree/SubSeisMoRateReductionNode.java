package scratch.nshm23.logicTree;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;

import scratch.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.SubSeisMoRateReduction;

public enum SubSeisMoRateReductionNode implements LogicTreeNode {
	FAULT_SPECIFIC("Fault-Specific", "FaultSpec", .5d, SubSeisMoRateReduction.FAULT_SPECIFIC_IMPLIED_FROM_SUPRA_B),
	SYSTEM_AVG("System-Average", "SysAvg", .5d, SubSeisMoRateReduction.SYSTEM_AVG_IMPLIED_FROM_SUPRA_B);
	
	private String name;
	private String shortName;
	private double weight;
	private SubSeisMoRateReduction choice;

	private SubSeisMoRateReductionNode(String name, String shortName, double weight, SubSeisMoRateReduction choice) {
		this.name = name;
		this.shortName = shortName;
		this.weight = weight;
		this.choice = choice;
	}
	
	public SubSeisMoRateReduction getChoice() {
		return choice;
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
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return weight;
	}

	@Override
	public String getFilePrefix() {
		return shortName;
	}

}
