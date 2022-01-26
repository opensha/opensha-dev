package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;

import com.google.common.base.Preconditions;

public class TestReweightFromMisfits {

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");

//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-max_dist-FM3_1-CoulombRupSet-TotNuclRate");
//		File mainDir = new File(invDir, "2021_12_17-nshm23_draft_branches-no_seg-FM3_1-CoulombRupSet");
//		File mainDir = new File(invDir, "2021_12_17-u3_branches-coulomb-FM3_1-5h");
//		File mainDir = new File(invDir, "2022_01_07-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-CoulombRupSet-SubB1-175_samples");
//		File mainDir = new File(invDir, "2022_01_10-nshm23_draft_branches-no_seg-reweighted_even_fit-conserve-FM3_1-CoulombRupSet-SubB1-105_samples");
//		File mainDir = new File(invDir, "2022_01_11-nshm23_draft_branches-no_seg-reweighted_even_fit-conserve-aggressive-FM3_1-CoulombRupSet-SubB1-105_samples");
//		File mainDir = new File(invDir, "2022_01_11-nshm23_draft_branches-no_seg-reweighted_even_fit-conserve-aggressiver-FM3_1-CoulombRupSet-SubB1-105_samples");
		File mainDir = new File(invDir, "2022_01_18-nshm23_draft_branches-no_seg-reweighted_even_fit-FM3_1-U3RupSet-SubB1-5000ip");
		
		File resultsFile = new File(mainDir, "results.zip");
		
		File outputDir = new File(mainDir, "logic_tree_misfits");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
		Map<LogicTreeLevel<?>, HashSet<LogicTreeNode>> levelNodes = new HashMap<>();
		List<? extends LogicTreeLevel<?>> levels = tree.getLevels();
		for (LogicTreeLevel<?> level : levels)
			levelNodes.put(level, new HashSet<>());
		
		for (LogicTreeBranch<?> branch : tree) {
			for (int i=0; i<levels.size(); i++) {
				LogicTreeNode node = branch.getValue(i);
				levelNodes.get(levels.get(i)).add(node);
			}
		}
		
		List<LogicTreeBranch<?>> branches = new ArrayList<>();
		List<InversionMisfitStats> branchStats = new ArrayList<>();
		List<Double> branchAvgMisfits = new ArrayList<>();
		
		Quantity targetQuantity = Quantity.MAD;
		boolean onlyUncertWeighted = true;
		
		ZipFile	zip = new ZipFile(resultsFile);
		
		for (LogicTreeBranch<?> branch : tree) {
			branches.add(branch);
			
			String entryName = "solution_logic_tree/";
			for (int i=0; i<branch.size(); i++) {
				LogicTreeLevel<?> level = branch.getLevel(i);
				if (level.affects(InversionMisfitStats.MISFIT_STATS_FILE_NAME, true))
					entryName += branch.getValue(i).getFilePrefix()+"/";
			}
			entryName += InversionMisfitStats.MISFIT_STATS_FILE_NAME;
			System.out.println("Loading "+entryName);
			ZipEntry entry = zip.getEntry(entryName);
			Preconditions.checkNotNull(entry, "Entry not found: %s", entryName);
			
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			InversionMisfitStats stats = new InversionMisfitStats(null);
			stats.initFromCSV(csv);
			branchStats.add(stats);
			
			double avgVal = 0;
			int num = 0;
			for (MisfitStats misfits : stats.getStats()) {
				if (onlyUncertWeighted && misfits.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
					continue;
				
				avgVal += misfits.get(targetQuantity);
				num++;
			}
			avgVal /= (double)num;
			
			branchAvgMisfits.add(avgVal);
		}
		
		zip.close();
		
		double overallAverage = branchAvgMisfits.stream().mapToDouble(d->d).average().getAsDouble();
		System.out.println("Overall average misfit: "+overallAverage);
		
		List<Double> weights = new ArrayList<>();
		double sumWeight = 0d;
		for (int i=0; i<branchAvgMisfits.size(); i++) {
			double weight = overallAverage/branchAvgMisfits.get(i);
			weights.add(weight);
			sumWeight += weight;
		}
		// normalize
		for (int i=0; i<weights.size(); i++)
			weights.set(i, weights.get(i)/sumWeight);
		sumWeight = 1d;
		
		System.out.println(targetQuantity+" targeted branch weights");
		
		for (LogicTreeLevel<?> level : levels) {
			List<LogicTreeNode> nodes = new ArrayList<>(levelNodes.get(level));
			if (nodes.size() < 2)
				continue;
			nodes.sort(new Comparator<LogicTreeNode>() {

				@Override
				public int compare(LogicTreeNode o1, LogicTreeNode o2) {
					return o1.getName().compareTo(o2.getName());
				}
			});
			System.out.println(level.getName());
			
			for (LogicTreeNode node : nodes) {
				double nodeSumWeight = 0d;
				for (int i=0; i<branches.size(); i++) {
					if (branches.get(i).hasValue(node)) {
						nodeSumWeight += weights.get(i);
					}
				}
				double calcWeight = nodeSumWeight/sumWeight;
				double origWeight = node.getNodeWeight(null);
				System.out.println("\t"+node.getShortName()+":\torig="+(float)origWeight+"\tcalc="+(float)calcWeight);
			}
		}
	}

}
