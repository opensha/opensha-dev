package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Random;

import org.opensha.commons.logicTree.LogicTreeLevel.RandomlyGeneratedLevel;

public class RandomDefModSampleLevel extends RandomlyGeneratedLevel<RandomDefModSampleNode> {
	
	public RandomDefModSampleLevel() {
		
	}
	
	public RandomDefModSampleLevel(int numSamples) {
		this(numSamples, new Random().nextLong());
	}
	
	public RandomDefModSampleLevel(int numSamples, long seed) {
		super("Random Deformation Model Sample", "DMSample",
				"Deformation Model Sample ", "DMSample", "DMSample");
		build(seed, numSamples);
	}

	@Override
	public RandomDefModSampleNode buildNodeInstance(int index, long seed, double weight) {
		return new RandomDefModSampleNode(getNodeName(index), getNodeShortName(index), getNodeFilePrefix(index), weight, seed);
	}

	@Override
	public Class<? extends RandomDefModSampleNode> getType() {
		return RandomDefModSampleNode.class;
	}

}
