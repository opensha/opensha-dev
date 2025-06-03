package scratch.kevin.prvi25;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.FileBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.FileBackedNode;
import org.opensha.commons.logicTree.treeCombiner.AbstractLogicTreeCombiner;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.HazardMapCombinationProcessor;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.SiteHazardCurveCombinationProcessor;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.SolutionLogicTreeCombinationProcessor;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.imr.logicTree.ScalarIMR_ParamsLogicTreeNode;
import org.opensha.sha.imr.logicTree.ScalarIMRsLogicTreeNode;

import com.google.common.base.Preconditions;

public class CrustalSubductionLogicTreeCombine extends AbstractLogicTreeCombiner {
	
	private static final boolean REVERSE = false;
	private static final IncludeBackgroundOption GRID_SEIS_DEFAULT = IncludeBackgroundOption.EXCLUDE;
	private static final String RESULTS_FILE_NAME_DEFAULT = "results.zip";
	private static final String HAZARD_DIR_NAME_DEFAULT = "results";
	private static final String GRID_REG_FILE_NAME_DEFAULT = "results_hazard.zip"; // get region from hazard zip

	public static void main(String[] args) throws IOException {
		CommandLine cmd = FaultSysTools.parseOptions(createOptions(), args, CrustalSubductionLogicTreeCombine.class);
		args = cmd.getArgs();
		Preconditions.checkArgument(args.length == 3, "USAGE: <crustal-dir> <subduction-dir> <output-dir>");
		ModuleContainer.VERBOSE_DEFAULT = false;
		File crustalDir = new File(args[0]);
		Preconditions.checkState(crustalDir.exists(),
				"Crustal directory doesn't exist: %s", crustalDir.getAbsolutePath());
		File subductionDir = new File(args[1]);
		Preconditions.checkState(subductionDir.exists(),
				"Subduction directory doesn't exist: %s", subductionDir.getAbsolutePath());
		File outputDir = new File(args[2]);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Output directory doesn't exist and couldn't be created: %s", outputDir.getAbsolutePath());
		
		String resultsFileName = RESULTS_FILE_NAME_DEFAULT;
		if (cmd.hasOption("results-file-name"))
			resultsFileName = cmd.getOptionValue("results-file-name");
		String hazardDirName = HAZARD_DIR_NAME_DEFAULT;
		if (cmd.hasOption("hazard-dir-name"))
			hazardDirName = cmd.getOptionValue("hazard-dir-name");
		String gridRegFileName = GRID_REG_FILE_NAME_DEFAULT;
		if (cmd.hasOption("grid-reg-file-name"))
			gridRegFileName = cmd.getOptionValue("grid-reg-file-name");
		IncludeBackgroundOption bgOp = GRID_SEIS_DEFAULT;
		if (cmd.hasOption("gridded-seis"))
			bgOp = IncludeBackgroundOption.valueOf(cmd.getOptionValue("gridded-seis"));
		
		LogicTree<?> crustalLT;
		SolutionLogicTree crustalSLT;
		LogicTree<?> subductionLT;
		SolutionLogicTree subductionSLT;
		if (cmd.hasOption("logic-tree-file-name")) {
			String logicTreeFileName = cmd.getOptionValue("logic-tree-file-name");
			File crustalLogicTreefile = new File(crustalDir, logicTreeFileName);
			crustalLT = LogicTree.read(crustalLogicTreefile);
			File subductionLogicTreeFile = new File(subductionDir, logicTreeFileName);
			subductionLT = LogicTree.read(subductionLogicTreeFile);
			File crustalSLTFile = new File(crustalDir, resultsFileName);
			if (crustalSLTFile.exists())
				crustalSLT = SolutionLogicTree.load(crustalSLTFile, crustalLT);
			else
				crustalSLT = null;
			File subductionSLTFile = new File(subductionDir, resultsFileName);
			if (subductionSLTFile.exists())
				subductionSLT = SolutionLogicTree.load(subductionSLTFile, subductionLT);
			else
				subductionSLT = null;
		} else {
			crustalSLT = SolutionLogicTree.load(new File(crustalDir, resultsFileName));
			crustalLT = crustalSLT.getLogicTree();
			subductionSLT = SolutionLogicTree.load(new File(subductionDir, resultsFileName));
			subductionLT = subductionSLT.getLogicTree();
		}
		
		if (crustalSLT != null)
			crustalSLT.setVerbose(false);
		if (subductionSLT != null)
			subductionSLT.setVerbose(false);
		
		GriddedRegion gridReg = loadGridReg(new File(crustalDir, gridRegFileName));
		Preconditions.checkState(gridReg.equalsRegion(loadGridReg(new File(subductionDir, gridRegFileName))),
				"Crustal and subduction gridded regions differ");
		
		String outputHazardFileName;
		if (cmd.hasOption("output-hazard-file-name"))
			outputHazardFileName = cmd.getOptionValue("output-hazard-file-name");
		else if (bgOp == IncludeBackgroundOption.EXCLUDE)
			outputHazardFileName = "results_hazard.zip";
		else if (bgOp == IncludeBackgroundOption.INCLUDE)
			outputHazardFileName = "results_hazard_full_gridded.zip";
		else if (bgOp == IncludeBackgroundOption.ONLY)
			outputHazardFileName = "results_hazard_gridded_only.zip";
		else
			throw new IllegalStateException();
		
		boolean averageERF = cmd.hasOption("average-erf");
		boolean averageGMM = cmd.hasOption("average-gmm");
		boolean averageEither = averageERF || averageGMM;
		
		List<LogicTreeLevel<? extends LogicTreeNode>> averageAcrossLevels = null;
		
		if (averageERF) {
			Preconditions.checkState(!averageGMM);
			averageAcrossLevels = new ArrayList<>();
			for (LogicTreeLevel<?> level : crustalLT.getLevels())
				if (!isLevelGMM(level))
					averageAcrossLevels.add(level);
			for (LogicTreeLevel<?> level : subductionLT.getLevels())
				if (!isLevelGMM(level))
					averageAcrossLevels.add(level);
			System.out.println("Will average across "+averageAcrossLevels.size()+" ERF levels");
		} else if (averageGMM) {
			averageAcrossLevels = new ArrayList<>();
			for (LogicTreeLevel<?> level : crustalLT.getLevels())
				if (isLevelGMM(level))
					averageAcrossLevels.add(level);
			for (LogicTreeLevel<?> level : subductionLT.getLevels())
				if (isLevelGMM(level))
					averageAcrossLevels.add(level);
			System.out.println("Will average across "+averageAcrossLevels.size()+" GMM levels");
		}
		
//		CrustalSubductionLogicTreeCombine combiner = new CrustalSubductionLogicTreeCombine(
//				crustalLT, crustalSLT, new File(crustalDir, hazardDirName), bgOp,
//				subductionLT, subductionSLT, new File(subductionDir, hazardDirName), bgOp,
//				resultsOutputFile,
//				new File(outputDir, outputHazardFileName), gridReg);
		CrustalSubductionLogicTreeCombine combiner;
		if (REVERSE)
			combiner = new CrustalSubductionLogicTreeCombine(subductionLT, crustalLT, averageAcrossLevels);
		else
			combiner = new CrustalSubductionLogicTreeCombine(crustalLT, subductionLT, averageAcrossLevels);
		
		if (!cmd.hasOption("disable-write-results") && !averageEither) {
			File resultsOutputFile = new File(outputDir, resultsFileName);
			Preconditions.checkNotNull(crustalSLT);
			Preconditions.checkNotNull(subductionSLT);
			boolean serializeGridded = false; // TODO option?
			if (REVERSE)
				combiner.addProcessor(new SolutionLogicTreeCombinationProcessor(subductionSLT, crustalSLT,
						resultsOutputFile, true, true, serializeGridded));
			else
				combiner.addProcessor(new SolutionLogicTreeCombinationProcessor(crustalSLT, subductionSLT,
						resultsOutputFile, true, true, serializeGridded));
		}
		
		File crustalHazardDir = new File(crustalDir, hazardDirName);
		File subductionHazardDir = new File(subductionDir, hazardDirName);
		if (!cmd.hasOption("no-maps") && crustalHazardDir.exists() && subductionHazardDir.exists()) {
			System.out.println("Will combine hazard maps");
			HazardMapCombinationProcessor processor;
			if (REVERSE)
				processor = new HazardMapCombinationProcessor(subductionHazardDir, bgOp,
						crustalHazardDir, bgOp, new File(outputDir, outputHazardFileName), gridReg);
			else
				processor = new HazardMapCombinationProcessor(crustalHazardDir, bgOp,
						subductionHazardDir, bgOp, new File(outputDir, outputHazardFileName), gridReg);
			
			if (cmd.hasOption("disable-preload"))
				processor.setPreloadInnerCurves(false);
			
			combiner.addProcessor(processor);
		}
		
		String siteHazardFileName;
		if (cmd.hasOption("site-hazard-file-name"))
			siteHazardFileName = cmd.getOptionValue("site-hazard-file-name");
		else if (bgOp == IncludeBackgroundOption.EXCLUDE)
			siteHazardFileName = "results_hazard_sites.zip";
		else if (bgOp == IncludeBackgroundOption.INCLUDE)
			siteHazardFileName = "results_hazard_sites_full_gridded.zip";
		else if (bgOp == IncludeBackgroundOption.ONLY)
			siteHazardFileName = "results_hazard_sites_only.zip";
		else
			throw new IllegalStateException();
		
		File crustalSiteFile = new File(crustalDir, siteHazardFileName);
		File subductionSiteFile = new File(subductionDir, siteHazardFileName);
		if (!cmd.hasOption("no-curves") && crustalHazardDir.exists() && subductionHazardDir.exists() && !averageEither) {
			System.out.println("Will combine site hazard curves");
			File sitesOutputFile = new File(outputDir, siteHazardFileName);
			if (REVERSE)
				combiner.addProcessor(new SiteHazardCurveCombinationProcessor(subductionSiteFile, crustalSiteFile, sitesOutputFile));
			else
				combiner.addProcessor(new SiteHazardCurveCombinationProcessor(crustalSiteFile, subductionSiteFile, sitesOutputFile));
		}
		
		long treeSize = (long)crustalLT.size()*(long)subductionLT.size();
		if (cmd.hasOption("pairwise-samples")) {
			// more memory efficient pairwise-sampling version
			int samples = Integer.parseInt(cmd.getOptionValue("pairwise-samples"));
			long seed = cmd.hasOption("rand-seed") ? Long.parseLong(cmd.getOptionValue("rand-seed")) : treeSize*(long)samples;
			combiner.pairwiseSampleTree(samples, seed);
		}
		
		if (cmd.hasOption("samples")) {
			int samples = Integer.parseInt(cmd.getOptionValue("samples"));
			if (samples < treeSize) {
				long seed = cmd.hasOption("rand-seed") ? Long.parseLong(cmd.getOptionValue("rand-seed")) : treeSize*(long)samples;
				combiner.sampleTree(samples, seed);
			}
		}
		
		LogicTree<LogicTreeNode> combTree = combiner.getCombTree();
		if (!averageERF && !averageGMM) {
			File treeOutputFile;
			if (cmd.hasOption("logic-tree-file-name"))
				treeOutputFile = new File(outputDir, cmd.getOptionValue("logic-tree-file-name"));
			else
				treeOutputFile = new File(outputDir, "logic_tree_"+bgOp.name()+".json");
			combTree.write(treeOutputFile);
		}
		
		try {
			combiner.processCombinations();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static Options createOptions() {
		Options ops = new Options();

		ops.addOption(null, "pairwise-samples", true, "Number of pairwise-samples to draw. This is more memory efficient and doesn't require "
				+ "instantiating the entire tree in memory first.");
		ops.addOption(null, "samples", true, "Number of samples to draw");
		ops.addOption(null, "rand-seed", true, "Random seed for use when sampling");
		ops.addOption(null, "gridded-seis", true, "Gridded seismicity option. One of "
				+FaultSysTools.enumOptions(IncludeBackgroundOption.class)+". Default: "+GRID_SEIS_DEFAULT.name());
		ops.addOption(null, "results-file-name", true, "Results file name. Default: "+RESULTS_FILE_NAME_DEFAULT);
		ops.addOption(null, "logic-tree-file-name", true, "Logic tree file name. Default uses that from results file");
		ops.addOption(null, "hazard-dir-name", true, "Hazard dir name. Default: "+HAZARD_DIR_NAME_DEFAULT);
		ops.addOption(null, "grid-reg-file-name", true, "Gridded region file name. Default: "+GRID_REG_FILE_NAME_DEFAULT);
		ops.addOption(null, "disable-write-results", false, "Flag to disable writing the combined solution logic tree fule");
		ops.addOption(null, "output-hazard-file-name", true, "Name for output hazard file; default depends on background seismicity.");
		ops.addOption(null, "site-hazard-file-name", true, "Name of site hazard file; default depends on background seismicity.");
		ops.addOption(null, "disable-preload", false, "Flag to disable pre-loading of slab hazard curves (to quickly get to the build loop)");
		ops.addOption(null, "no-maps", false, "Flag to disable combining hazard maps");
		ops.addOption(null, "no-curves", false, "Flag to disable combining site hazard curves");
		ops.addOption(null, "average-erf", false, "Flag to average all ERF branches (retain GMM variability only)");
		ops.addOption(null, "average-gmm", false, "Flag to average all GMM branches (retain ERF variability only)");
		
		return ops;
	}

	public CrustalSubductionLogicTreeCombine(LogicTree<?> outerLT, LogicTree<?> innerLT, List<LogicTreeLevel<?>> averageAcrossLevels) {
		super(outerLT, innerLT, null, averageAcrossLevels);
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
		if (REVERSE)
			remapTree(tree, levelRemaps, nodeRemaps, "Subduction", "Sub");
		else
			remapTree(tree, levelRemaps, nodeRemaps, "Crustal", "Crust");
	}

	@Override
	protected void remapInnerTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps) {
		if (REVERSE)
			remapTree(tree, levelRemaps, nodeRemaps, "Crustal", "Crust");
		else
			remapTree(tree, levelRemaps, nodeRemaps, "Subduction", "Sub");
	}

}
