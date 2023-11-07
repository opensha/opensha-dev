package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Random;

import org.opensha.commons.logicTree.LogicTreeLevel.RandomlySampledLevel;

public class RandomDefModSampleLevel extends RandomlySampledLevel<RandomDefModSampleNode> {
	
	public RandomDefModSampleLevel() {
		
	}
	
	public RandomDefModSampleLevel(int numSamples) {
		this(numSamples, new Random());
	}
	
	public RandomDefModSampleLevel(int numSamples, Random r) {
		buildNodes(r, numSamples);
	}

	@Override
	public String getShortName() {
		return "DMSample";
	}

	@Override
	public String getName() {
		return "Random Deformation Model Sample";
	}

	@Override
	public RandomDefModSampleNode buildNodeInstance(int index, long seed, double weight) {
		return new RandomDefModSampleNode(index, seed, weight);
	}

	@Override
	public Class<? extends RandomDefModSampleNode> getType() {
		return RandomDefModSampleNode.class;
	}

}
