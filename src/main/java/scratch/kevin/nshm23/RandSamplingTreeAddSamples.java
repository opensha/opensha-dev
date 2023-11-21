package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.RandomlySampledLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;

public class RandSamplingTreeAddSamples {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions");
		
		File origFile = new File(mainDir, "2023_11_14-nshm23_branches-dm_sampling-randB-randSeg-mini_one_fifth-"
				+ "NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/logic_tree.json");
		LogicTree<?> origTree = LogicTree.read(origFile);
		File destFile = new File(mainDir, "2023_11_20-nshm23_branches-dm_sampling-randB-randSeg-"
				+ "NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/logic_tree.json");
		int sampleMult = 5;
		
		int origSize = origTree.size();
		int newSize = origTree.size()*sampleMult;
		int numAdded = origTree.size()*(sampleMult-1);
		
		System.out.println("Will add "+numAdded+" branches: "+origSize+" -> "+newSize);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> modLevels = new ArrayList<>();
		
		Random rand = new Random(987654321l);
		
		int numLevels = origTree.getLevels().size();
		List<LinkedList<RandomlySampledNode>> randNodeStacks = new ArrayList<>();
		for (int l=0; l<numLevels; l++) {
			LogicTreeLevel<?> sourceLevel = origTree.getLevels().get(l);
			if (sourceLevel instanceof RandomlySampledLevel<?>) {
				System.out.println("Building node list for "+sourceLevel.getName());
				HashSet<Long> prevSeeds = new HashSet<>(newSize);
				List<Long> seeds = new ArrayList<>();
				// add the original nodes
				for (LogicTreeNode sourceNode : sourceLevel.getNodes()) {
					Preconditions.checkState(sourceNode instanceof RandomlySampledNode);
					RandomlySampledNode randNode = (RandomlySampledNode)sourceNode;
					Preconditions.checkState(!prevSeeds.contains(randNode.getSeed()));
					prevSeeds.add(randNode.getSeed());
					seeds.add(randNode.getSeed());
				}
				// add additional nodes
				for (int i=0; i<numAdded; i++) {
					long seed = rand.nextLong();
					Preconditions.checkState(!prevSeeds.contains(seed));
					prevSeeds.add(seed);
					seeds.add(seed);
				}
				try {
					Constructor<? extends RandomlySampledLevel> constructor = ((RandomlySampledLevel<?>)sourceLevel).getClass().getConstructor();
					constructor.setAccessible(true);
					RandomlySampledLevel<?> modLevel = constructor.newInstance();
					modLevel.buildNodes(seeds, 1d);
					modLevels.add(modLevel);
					randNodeStacks.add(new LinkedList<>(modLevel.getNodes()));
				} catch (InstantiationException | IllegalAccessException | IllegalArgumentException
						| InvocationTargetException | NoSuchMethodException | SecurityException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			} else {
				modLevels.add(sourceLevel);
				randNodeStacks.add(null);
			}
		}
		
		List<LogicTreeBranch<LogicTreeNode>> modBranches = new ArrayList<>(newSize);
		for (int i=0; i<newSize; i++)
			modBranches.add(null);
		double sumWeight = 0d;
		for (int round=0; round<sampleMult; round++) {
			for (int origIndex=0; origIndex<origSize; origIndex++) {
				LogicTreeBranch<?> origBranch = origTree.getBranch(origIndex);
				LogicTreeBranch<LogicTreeNode> modBranch = new LogicTreeBranch<>(modLevels);
				for (int l=0; l<numLevels; l++) {
					LogicTreeNode node = origBranch.getValue(l);
					if (node instanceof RandomlySampledNode) {
						RandomlySampledNode stackNode = randNodeStacks.get(l).removeFirst();
						if (round == 0) {
							Preconditions.checkState(stackNode.getSeed() == ((RandomlySampledNode)node).getSeed());
							Preconditions.checkState(stackNode.getFilePrefix().equals(node.getFilePrefix()));
						}
						modBranch.setValue(l, stackNode);
					} else {
						modBranch.setValue(l, node);
					}
				}
//				int index = round*origSize + origIndex;
				int index = origIndex*sampleMult + round;
				System.out.println("Branch "+index+": "+modBranch);
				modBranches.set(index, modBranch);
				sumWeight += modBranch.getBranchWeight();
			}
		}
		double origWeight = 0d;
		for (LogicTreeBranch<?> branch : origTree)
			origWeight += branch.getBranchWeight();
		System.out.println("Done building tree. origWeight="+(float)origWeight+", modWeight="+(float)sumWeight);
		
		LogicTree<?> modTree = LogicTree.fromExisting(modLevels, modBranches);
		modTree.write(destFile);
	}

}
