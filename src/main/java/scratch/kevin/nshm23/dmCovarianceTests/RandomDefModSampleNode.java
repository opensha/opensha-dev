package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Objects;

import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;

@Affects(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public class RandomDefModSampleNode implements RandomlySampledNode {

	private String name;
	private String shortName;
	private String prefix;
	private double weight;
	private long seed;

	@SuppressWarnings("unused") // for deserialization
	private RandomDefModSampleNode() {}
	
	RandomDefModSampleNode(int index, long seed, double weight) {
		init("Deformation Model Sample "+index, "DMSample"+index, "DMSample"+index, weight, seed);
	}

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return weight;
	}

	@Override
	public String getFilePrefix() {
		return prefix;
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
	public long getSeed() {
		return seed;
	}

	@Override
	public void init(String name, String shortName, String prefix, double weight, long seed) {
		this.name = name;
		this.shortName = shortName;
		this.prefix = prefix;
		this.weight = weight;
		this.seed = seed;
	}

	@Override
	public int hashCode() {
		return Objects.hash(name, prefix, seed, shortName, weight);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		RandomDefModSampleNode other = (RandomDefModSampleNode) obj;
		return Objects.equals(name, other.name) && Objects.equals(prefix, other.prefix) && seed == other.seed
				&& Objects.equals(shortName, other.shortName)
				&& Double.doubleToLongBits(weight) == Double.doubleToLongBits(other.weight);
	}

}
