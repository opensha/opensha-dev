package scratch.nshm23.logicTree;

import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;

@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public enum SegmentationModel implements LogicTreeNode {
	SHAW_R0_1("Shaw & Dieterich (2007) R₀=1", "ShawR₀=1", 0.1d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 1);
		}
	},
	SHAW_R0_2("Shaw & Dieterich (2007) R₀=2", "ShawR₀=2", 0.2d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 2);
		}
	},
	SHAW_R0_3("Shaw & Dieterich (2007) R₀=3", "ShawR₀=3", 0.4d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 3);
		}
	},
	SHAW_R0_4("Shaw & Dieterich (2007) R₀=4", "ShawR₀=4", 0.2d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 4);
		}
	},
	SHAW_R0_5("Shaw & Dieterich (2007) R₀=5", "ShawR₀=5", 0.0d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 5);
		}
	},
	SHAW_R0_6("Shaw & Dieterich (2007) R₀=6", "ShawR₀=6", 0.0d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return new Shaw07JumpDistProb(1, 6);
		}
	},
	NONE("None", "None", .1d) {
		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet) {
			return null;
		}
	};
	
	private String name;
	private String shortName;
	private double weight;

	private SegmentationModel(String name, String shortName, double weight) {
		this.name = name;
		this.shortName = shortName;
		this.weight = weight;
	}
	
	public abstract JumpProbabilityCalc getModel(FaultSystemRupSet rupSet);

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
		return shortName.replace("R₀=", "R0_");
	}

}
