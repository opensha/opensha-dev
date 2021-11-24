package scratch.nshm23.logicTree;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;

public enum SupraSeisBValue implements LogicTreeNode {
	
	B_0p0(0d, 0.05),
	B_0p3(0.3, 0.1),
	B_0p5(0.5, 0.15),
	B_0p7(0.7, 0.2),
	B_0p8(0.8, 0.25),
	B_0p9(0.9, 0.15),
	B_1p0(1d, 0.1);
	
	public final double bValue;
	public final double weight;

	private SupraSeisBValue(double bValue, double weight) {
		this.bValue = bValue;
		this.weight = weight;
	}

	@Override
	public String getShortName() {
		return getName();
	}

	@Override
	public String getName() {
		return "b="+(float)bValue;
	}

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return weight;
	}

	@Override
	public String getFilePrefix() {
		return "supra_b="+(float)bValue;
	}
	
	public static void main(String[] args) {
		double totWeight = 0d;
		for (SupraSeisBValue b : values())
			totWeight += b.weight;
		System.out.println("tot weight: "+(float)totWeight);
	}

}
