package scratch.kevin.nshm26;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.ContinuousDistribution.Sampler;
import org.apache.commons.statistics.distribution.CorrTruncatedNormalDistribution;
import org.apache.commons.statistics.distribution.TruncatedNormalDistribution;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlyGeneratedNode;
import org.opensha.commons.logicTree.LogicTreeNode.ValuedLogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;

import com.google.common.base.Preconditions;

public class UpdatedRandTreeSerialzationTests {

	public static void main(String[] args) throws IOException {
		List<LogicTree<LogicTreeNode>> inTrees = new ArrayList<>();
		
		inTrees.add(LogicTree.read(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2025_01_17-prvi25_crustal_branches-dmSample10x/logic_tree.json")));

		TestValuedLevel testValLevel = new TestValuedLevel();
		testValLevel.build(12345l, 1000);
		inTrees.add(LogicTree.buildExhaustive(List.of(testValLevel), true));
		
		ContinuousDistribution dist = CorrTruncatedNormalDistribution.of(1d, 0.1, 0.7, 1.3);
		Sampler sampler = dist.createSampler(RandomSource.XO_RO_SHI_RO_128_PP.create(123456l));
		MinMaxAveTracker sampleTrack = new MinMaxAveTracker();
		for (int i=0; i<100000; i++) {
			sampleTrack.addValue(sampler.sample());
		}
		System.out.println(sampleTrack);
		System.out.println("Dist mean: "+dist.getMean());
		System.out.println("Dist bounds: "+dist.getSupportLowerBound()+", "+dist.getSupportUpperBound());
		TestDistLevel testDistLevel = new TestDistLevel(dist);
		testDistLevel.build(12345l, 1000);
		inTrees.add(LogicTree.buildExhaustive(List.of(testDistLevel), true));
		
		for (int t=0; t<inTrees.size(); t++) {
			LogicTree<LogicTreeNode> tree = inTrees.get(t);
			System.out.println("Processing tree "+t+" with "+tree.size()+" branches");
			
			File treeOutFile = new File("/tmp/tree_out_"+t+".json");
			System.out.println("\tWriting tree to: "+treeOutFile.getAbsolutePath());
			tree.write(treeOutFile);
			
			File branchOutFile = new File("/tmp/branch_out_"+t+".json");
			System.out.println("\tWriting branch 0 to: "+branchOutFile.getAbsolutePath());
			tree.getBranch(0).writeToFile(branchOutFile);
			
			System.out.println("\tReading tree back in");
			LogicTree<LogicTreeNode> tree2 = LogicTree.read(treeOutFile);
			
			System.out.println("\tLevels:");
			for (LogicTreeLevel<?> level : tree2.getLevels())
				System.out.println("\t\t"+level.getName()+" ("+level.getClass()+")");

			System.out.println("\tReading branch 0 back in");
			LogicTreeBranch<LogicTreeNode> loadedBranch0 = LogicTreeBranch.read(branchOutFile);
			for (int i=-1; i<tree2.size(); i++) {
				LogicTreeBranch<LogicTreeNode> branch1, branch2;
				if (i < 0) {
					branch1 = tree.getBranch(0);
					branch2 = loadedBranch0;
				} else {
					branch1 = tree.getBranch(i);
					branch2 = tree2.getBranch(i);
				}
				List<RandomlyGeneratedNode> randGenNodes1 = branch1.getValues(RandomlyGeneratedNode.class);
				List<RandomlyGeneratedNode> randGenNodes2 = branch2.getValues(RandomlyGeneratedNode.class);
				if (randGenNodes1 == null) {
					Preconditions.checkState(randGenNodes2 == null);
				} else {
					Preconditions.checkState(randGenNodes1.size() == randGenNodes2.size());
					for (int n=0; n<randGenNodes1.size(); n++) {
						RandomlyGeneratedNode gen1 = randGenNodes1.get(n);
						RandomlyGeneratedNode gen2 = randGenNodes2.get(n);
						Preconditions.checkState(gen1.getSeed() == gen2.getSeed());
						Preconditions.checkState(gen1.equals(gen2));
					}
				}
				List<ValuedLogicTreeNode<?>> valNodes1 = branch1.getValues(ValuedLogicTreeNode.class);
				List<ValuedLogicTreeNode<?>> valNodes2 = branch2.getValues(ValuedLogicTreeNode.class);
				if (valNodes1 == null) {
					Preconditions.checkState(valNodes2 == null);
				} else {
					Preconditions.checkState(valNodes1.size() == valNodes2.size(),
							"Branch %s: Valued node list size mismatch: %s vs %s", i, valNodes1.size(), valNodes2.size());
					for (int n=0; n<valNodes1.size(); n++) {
						ValuedLogicTreeNode<?> gen1 = valNodes1.get(n);
						ValuedLogicTreeNode<?> gen2 = valNodes2.get(n);
						Preconditions.checkState(gen1.getValue().equals(gen2.getValue()));
						Preconditions.checkState(gen1.equals(gen2));
					}
				}
			}
		}
		
	}
	
	private static class TestValuedLevel extends LogicTreeLevel.RandomlySampledLevel<Double> {

		public TestValuedLevel() {
			super("Test Valued Level", "TVL", "Random Node ", "RN", "RN");
		}

		@Override
		protected void doBuild(long seed, int numNodes, double weightEach) {
			Random rand = new Random(seed);
			super.build(()->rand.nextDouble(), numNodes, weightEach);
		}

		@Override
		public Class<? extends Double> getValueType() {
			return Double.class;
		}
		
	}
	
	private static class TestDistLevel extends LogicTreeLevel.ContinuousDistributionSampledLevel {
		
		private TestDistLevel() {
			this(null);
		};

		public TestDistLevel(ContinuousDistribution dist) {
			super("Test Distribution Level", "TestDistLevel", dist,
					"Distribution Sample ", "DistSample", "DistSample");
		}
		
	}

}

