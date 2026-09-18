package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Objects;

import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlyGeneratedNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;

@Affects(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public class RandomDefModSampleNode extends RandomlyGeneratedNode {

	private RandomDefModSampleNode() {
		super();
	}

	public RandomDefModSampleNode(String name, String shortName, String prefix, double weight, long seed) {
		super(name, shortName, prefix, weight, seed);
	}
}
