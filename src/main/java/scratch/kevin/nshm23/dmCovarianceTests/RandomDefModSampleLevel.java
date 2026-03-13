package scratch.kevin.nshm23.dmCovarianceTests;

import java.util.Random;

import org.opensha.commons.logicTree.LogicTreeLevel.RandomlyGeneratedLevel;

public class RandomDefModSampleLevel extends RandomlyGeneratedLevel<RandomDefModSampleNode> {
	
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
		return new RandomDefModSampleNode(getNodeName(index), getNodeShortName(index), getNodeFilePrefix(index), weight, seed);
	}

	@Override
	public Class<? extends RandomDefModSampleNode> getType() {
		return RandomDefModSampleNode.class;
	}

	@Override
	protected String getNodeNamePrefix() {
		return "Deformation Model Sample ";
	}

	@Override
	protected String getNodeShortNamePrefix() {
		return "DMSample";
	}

	@Override
	protected String getNodeFilePrefix() {
		return "DMSample";
	}

}
