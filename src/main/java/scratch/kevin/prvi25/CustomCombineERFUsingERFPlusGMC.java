package scratch.kevin.prvi25;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.FileBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.FileBackedNode;
import org.opensha.commons.logicTree.treeCombiner.AbstractLogicTreeCombiner;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.SiteHazardCurveCombinationProcessor;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.imr.logicTree.ScalarIMR_ParamsLogicTreeNode;
import org.opensha.sha.imr.logicTree.ScalarIMRsLogicTreeNode;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;

public class CustomCombineERFUsingERFPlusGMC extends AbstractLogicTreeCombiner {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/project2/scec_608/kmilner/fss_inversions");
		File crustalDir = new File(baseDir, "2025_09_12-prvi25_crustal_branches-dmSample10x");
		File subductionDir = new File(baseDir, "2025_09_12-prvi25_subduction_branches");
		File gmcDir = new File(baseDir, "2025_09_12-prvi25_crustal_subduction_combined_branches-gmTreeCalcs-vs760");
		File outputDir = new File(baseDir, "2025_09_12-prvi25_crustal_subduction_combined_branches-erfSampleWithCombTree");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File gmcTreeFile = new File(gmcDir, "logic_tree_full_gridded_sampled.json");
		LogicTree<?> gmcTree = LogicTree.read(gmcTreeFile);
		
		LogicTree<?> crustalTree = LogicTree.read(new File(crustalDir, "logic_tree_full_gridded.json"));
		LogicTree<?> subductionTree = LogicTree.read(new File(subductionDir, "logic_tree_full_gridded.json"));
		
		List<? extends LogicTreeLevel<?>> gmcLevels = gmcTree.getLevels();
		List<Integer> crustalLevelIndexes = new ArrayList<>();
		List<Integer> subLevelIndexes = new ArrayList<>();
		List<Integer> commonLevelIndexes = new ArrayList<>();
		boolean isSubLevel = false;
		for (int l=0; l<gmcLevels.size(); l++) {
			LogicTreeLevel<?> level = gmcLevels.get(l);
			if (isLevelGMM(level))
				continue;
			if (!isSubLevel && level.getName().toLowerCase().contains("interface"))
				isSubLevel = true;
			if (isSubLevel) {
				subLevelIndexes.add(l);
			} else {
				crustalLevelIndexes.add(l);
				if (level.getName().contains("Epoch"))
					commonLevelIndexes.add(l);
			}
		}
		for (int l : commonLevelIndexes) {
			String name = gmcLevels.get(l).getName();
			int matchingL = -1;
			for (int i=0; i<subductionTree.getLevels().size(); i++) {
				LogicTreeLevel<?> level = subductionTree.getLevels().get(i);
				String subName = level.getName();
				if (name.endsWith(level.getName())) {
					System.out.println("Matched common level "+name+" to "+subName);
					matchingL = i;
				}
			}
			Preconditions.checkState(matchingL >= 0);
			subLevelIndexes.add(matchingL, l);
		}
		System.out.println("Crustal levels:");
		for (int i=0; i<crustalLevelIndexes.size(); i++) {
			int l = crustalLevelIndexes.get(i);
			String gmcName = gmcLevels.get(l).getName();
			String crustalName = crustalTree.getLevels().get(i).getName();
			System.out.println("\t"+l+". "+gmcName+" -> "+crustalName);
			Preconditions.checkState(gmcName.endsWith(crustalName));
		}
		System.out.println("Subduction levels:");
		for (int i=0; i<subLevelIndexes.size(); i++) {
			int l = subLevelIndexes.get(i);
			String gmcName = gmcLevels.get(l).getName();
			String subName = subductionTree.getLevels().get(i).getName();
			System.out.println("\t"+l+". "+gmcName+" -> "+subName);
			Preconditions.checkState(gmcName.endsWith(subName));
		}
		Preconditions.checkState(crustalLevelIndexes.size() == crustalTree.getLevels().size());
		Preconditions.checkState(subLevelIndexes.size() == subductionTree.getLevels().size());
		
		List<LogicTreeLevel<? extends LogicTreeNode>> commonLevels = new ArrayList<>();
		for (int l : commonLevelIndexes) {
			int crustalIndex = crustalLevelIndexes.indexOf(l);
			LogicTreeLevel<?> level = crustalTree.getLevels().get(crustalIndex);
			System.out.println("Common level: "+level.getName());
			commonLevels.add(level);
		}
		
		Map<String, Integer> crustalIndexesMap = new HashMap<>();
		for (int i=0; i<crustalTree.size(); i++)
			crustalIndexesMap.put(crustalTree.getBranch(i).buildFileName(), i);
		String crustalPrefix0 = crustalIndexesMap.keySet().iterator().next();
		Map<String, Integer> subIndexesMap = new HashMap<>();
		for (int i=0; i<subductionTree.size(); i++)
			subIndexesMap.put(subductionTree.getBranch(i).buildFileName(), i);
		String subPrefix0 = subIndexesMap.keySet().iterator().next();
		
		List<Integer> crustalIndexes = new ArrayList<>(gmcTree.size());
		List<Integer> subductionIndexes = new ArrayList<>(gmcTree.size());
		List<Double> weights = new ArrayList<>(gmcTree.size());
		Joiner join = Joiner.on("_");
		for (LogicTreeBranch<?> gmcBranch : gmcTree) {
			List<String> crustalPrefixes = new ArrayList<>(crustalLevelIndexes.size());
			for (int l : crustalLevelIndexes) {
				String prefix = gmcBranch.getValue(l).getFilePrefix();
				if (prefix.startsWith("Crust"))
					prefix = prefix.substring(5);
				crustalPrefixes.add(prefix);
			}
			List<String> subPrefixes = new ArrayList<>(crustalLevelIndexes.size());
			for (int l : subLevelIndexes) {
				String prefix = gmcBranch.getValue(l).getFilePrefix();
				if (prefix.startsWith("Sub"))
					prefix = prefix.substring(3);
				subPrefixes.add(prefix);
			}
			String crustalPrefix = join.join(crustalPrefixes);
			Preconditions.checkState(crustalIndexesMap.containsKey(crustalPrefix),
					"No matching crustal prefix found:\n\tPrefix: %s\n\tBranch: %s\n\tPrefix example: %s",
					crustalPrefix, gmcBranch, crustalPrefix0);
			int crustalIndex = crustalIndexesMap.get(crustalPrefix);
			String subPrefix = join.join(subPrefixes);
			Preconditions.checkState(subIndexesMap.containsKey(subPrefix),
					"No matching sub prefix found:\n\tPrefix: %s\n\tBranch: %s\n\tPrefix example: %s",
					subPrefix, gmcBranch, subPrefix0);
			int subIndex = subIndexesMap.get(subPrefix);
			int combIndex = weights.size();
			crustalIndexes.add(crustalIndex);
			subductionIndexes.add(subIndex);
			weights.add(gmcTree.getBranchWeight(gmcBranch));
			if (weights.size() <= 10) {
				System.out.println("Comb branch "+combIndex+": "+gmcBranch);
				System.out.println("\tCrustal branch "+crustalIndex+": "+crustalTree.getBranch(crustalIndex));
				System.out.println("\tSub branch "+subIndex+": "+subductionTree.getBranch(subIndex));
			}
		}
		System.out.println("Matched all "+weights.size()+" indexes");
		
		CustomCombineERFUsingERFPlusGMC combiner = new CustomCombineERFUsingERFPlusGMC(
				crustalTree, subductionTree, commonLevels, null);
		
		combiner.pairwiseSampleExplicit(crustalIndexes, subductionIndexes, weights);
		
		File crustalSiteFile = new File(crustalDir, "results_hazard_sites_full_gridded.zip");
		File subductionSiteFile = new File(subductionDir, "results_hazard_sites_full_gridded.zip");
		System.out.println("Will combine site hazard curves");
		File sitesOutputFile = new File(outputDir, "results_hazard_sites_full_gridded.zip");
		combiner.addProcessor(new SiteHazardCurveCombinationProcessor(crustalSiteFile, subductionSiteFile, sitesOutputFile));
		
		LogicTree<LogicTreeNode> combTree = combiner.getCombTree();
		File treeOutputFile = new File(outputDir, "logic_tree_full_gridded.json");
		combTree.write(treeOutputFile);
		
		try {
			combiner.processCombinations();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public CustomCombineERFUsingERFPlusGMC(LogicTree<?> outerLT, LogicTree<?> innerLT,
			List<LogicTreeLevel<? extends LogicTreeNode>> commonLevels, List<LogicTreeLevel<?>> averageAcrossLevels) {
		super(outerLT, innerLT, commonLevels, averageAcrossLevels);
	}
	
	private static boolean isLevelGMM(LogicTreeLevel<?> level) {
		return ScalarIMR_ParamsLogicTreeNode.class.isAssignableFrom(level.getType())
				|| ScalarIMRsLogicTreeNode.class.isAssignableFrom(level.getType())
				|| level.getName().contains("GMM Epistemic")
				|| level.getName().contains("GMM Sigma");
	}
	
	private static GriddedRegion loadGridReg(File regFile) throws IOException {
		Preconditions.checkState(regFile.exists(), "Supplied region file doesn't exist: %s", regFile.getAbsolutePath());
		if (regFile.getName().toLowerCase().endsWith(".zip")) {
			// it's a zip file, assume it's a prior hazard calc
			ZipFile zip = new ZipFile(regFile);
			ZipEntry regEntry = zip.getEntry(MPJ_LogicTreeHazardCalc.GRID_REGION_ENTRY_NAME);
			System.out.println("Reading gridded region from zip file: "+regEntry.getName());
			BufferedReader bRead = new BufferedReader(new InputStreamReader(zip.getInputStream(regEntry)));
			GriddedRegion region = GriddedRegion.fromFeature(Feature.read(bRead));
			zip.close();
			return region;
		} else {
			Feature feature = Feature.read(regFile);
			return GriddedRegion.fromFeature(feature);
		}
	}
	
	private static void remapTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps, String nameAdd, String shortNameAdd) {
		for (LogicTreeLevel<?> level : tree.getLevels()) {
			String name = level.getName();
			String lowerName = name.toLowerCase();
			if (level.getName().equals(PRVI25_LogicTree.SEIS_EPOCH.getName()))
				// common to both, don't remap
				continue;
			if (lowerName.contains("crustal") || lowerName.contains("subduction")
					|| lowerName.contains("interface") || lowerName.contains("slab")
					|| lowerName.contains("muertos") || lowerName.contains("caribbean")) {
				// keep it as is
				levelRemaps.put(level, level);
				for (LogicTreeNode node : level.getNodes())
					nodeRemaps.put(node, node);
			} else {
				// remap to a file backed
				List<FileBackedNode> modNodes = new ArrayList<>();
				for (LogicTreeNode node : level.getNodes()) {
					FileBackedNode modNode = new FileBackedNode(nameAdd+" "+node.getName(), node.getShortName(),
							node.getNodeWeight(null), shortNameAdd+node.getFilePrefix());
					modNodes.add(modNode);
					nodeRemaps.put(node, modNode);
				}
				FileBackedLevel modLevel = new FileBackedLevel(nameAdd+" "+name, shortNameAdd+level.getShortName(), modNodes);
				modLevel.setAffected(level.getAffected(), level.getNotAffected(), false);
				levelRemaps.put(level, modLevel);
			}
		}
	}

	@Override
	protected void remapOuterTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps) {
		remapTree(tree, levelRemaps, nodeRemaps, "Crustal", "Crust");
	}

	@Override
	protected void remapInnerTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps) {
		remapTree(tree, levelRemaps, nodeRemaps, "Subduction", "Sub");
	}

}
