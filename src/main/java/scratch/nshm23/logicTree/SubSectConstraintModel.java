package scratch.nshm23.logicTree;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;

public enum SubSectConstraintModel implements LogicTreeNode {
	TOT_NUCL_RATE("Total Nucleation Rate", "TotNuclRate", .5d),
	NUCL_MFD("Nucleation MFD", "NuclMFD", .5d),
	NONE("None", "None", 0d);
	
	private String name;
	private String shortName;
	private double weight;

	private SubSectConstraintModel(String name, String shortName, double weight) {
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
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return weight;
	}

	@Override
	public String getFilePrefix() {
		return shortName;
	}

}
