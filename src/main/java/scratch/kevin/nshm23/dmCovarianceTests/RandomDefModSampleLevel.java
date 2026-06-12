package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Random;

import org.opensha.commons.logicTree.LogicTreeLevel.RandomlyGeneratedLevel;

public class RandomDefModSampleLevel extends RandomlyGeneratedLevel<RandomDefModSampleNode> {
	
	public RandomDefModSampleLevel(String name, String shortName) {
		super(name, shortName);
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
	public Class<? extends RandomDefModSampleNode> getType() {
		return RandomDefModSampleNode.class;
	}

	@Override
	public RandomDefModSampleNode build(Long seed, double weight, String name, String shortName, String filePrefix) {
		return new RandomDefModSampleNode(name, shortName, filePrefix, weight, seed);
	}

}
