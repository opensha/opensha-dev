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
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.hazard.AbstractLogicTreeHazardCombiner;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.imr.logicTree.ScalarIMR_ParamsLogicTreeNode;
import org.opensha.sha.imr.logicTree.ScalarIMRsLogicTreeNode;

import com.google.common.base.Preconditions;

public class SubductionInterfaceSlabHazardLogicTreeCombine extends AbstractLogicTreeHazardCombiner {

	public static void main(String[] args) throws IOException {
		CommandLine cmd = FaultSysTools.parseOptions(createOptions(), args, SubductionInterfaceSlabHazardLogicTreeCombine.class);
		args = cmd.getArgs();
		Preconditions.checkArgument(args.length == 1, "USAGE: <hazard-dir>");
		ModuleContainer.VERBOSE_DEFAULT = false;
		File hazardDir = new File(args[0]);
		Preconditions.checkState(hazardDir.exists(), "Directory doesn't exist: %s", hazardDir.getAbsolutePath());
		
		LogicTree<?> interfaceTree = LogicTree.read(new File(hazardDir, "logic_tree_full_gridded_interface_only.json"));
		LogicTree<?> slabTree = LogicTree.read(new File(hazardDir, "logic_tree_gridded_slab_only.json"));
		
		File interfaceHazardMapDir = new File(hazardDir, "results_full_gridded_interface");
		File slabHazardMapDir = new File(hazardDir, "results_gridded_only_slab");
		
		File outputDir = new File(hazardDir, "results_full_gridded");
		GriddedRegion gridReg = loadGridReg(new File(hazardDir, "gridded_region.geojson"));
		
		SubductionInterfaceSlabHazardLogicTreeCombine combiner = new SubductionInterfaceSlabHazardLogicTreeCombine(
				interfaceTree, slabTree);
		
		if (cmd.hasOption("disable-preload"))
			combiner.setPreloadInnerCurves(false);
		
		if (!cmd.hasOption("no-maps") && interfaceHazardMapDir.exists()) {
			System.out.println("Will combine hazard maps");
			combiner.setCombineHazardMaps(interfaceHazardMapDir, IncludeBackgroundOption.INCLUDE, slabHazardMapDir,
					IncludeBackgroundOption.ONLY, outputDir, gridReg);
		}
		
		File interfaceHazardSitesFile = new File(hazardDir, "results_full_gridded_sites_interface.zip");
		File slabHazardSitesFile = new File(hazardDir, "results_gridded_only_sites_slab.zip");
		
		if (!cmd.hasOption("no-curves") && interfaceHazardSitesFile.exists()) {
			System.out.println("Will combine hazard sites");
			File sitesOutputFile = new File(hazardDir, "results_full_gridded_sites");
			combiner.setCombineHazardCurves(interfaceHazardSitesFile, slabHazardSitesFile, sitesOutputFile);
		}
		
		System.out.println();
		List<LogicTreeLevel<?>> commonLevels = getCommonLevels(slabTree);
		System.out.println("Interface tree has "+interfaceTree.size()+" branches; Levels:");
		for (LogicTreeLevel<?> level : interfaceTree.getLevels())
			System.out.println("\t"+level.getName()+" (common="+commonLevels.contains(level)+")");
		System.out.println("Slab tree has "+slabTree.size()+" branches; Levels:");
		for (LogicTreeLevel<?> level : slabTree.getLevels())
			System.out.println("\t"+level.getName()+" (common="+commonLevels.contains(level)+")");
		System.out.println("Raw tree product is "+(interfaceTree.size()*slabTree.size()));
		System.out.println();
		
		if (cmd.hasOption("samples")) {
			int samples = Integer.parseInt(cmd.getOptionValue("samples"));
			if (cmd.hasOption("rand-seed"))
				combiner.sampleTree(samples, Long.parseLong(cmd.getOptionValue("rand-seed")));
			else
				combiner.sampleTree(samples);
		}
		
		combiner.getCombTree().write(new File(hazardDir, "logic_tree_full_gridded.json"));
		
		try {
			combiner.build();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public SubductionInterfaceSlabHazardLogicTreeCombine(LogicTree<?> interfaceTree, LogicTree<?> slabTree) {
		super(interfaceTree, slabTree, getCommonLevels(slabTree), null);
//		super(interfaceTree, interfaceHazardDir, IncludeBackgroundOption.INCLUDE,
//				slabTree, slabHazardDir, IncludeBackgroundOption.ONLY,
//				outputHazardFile, gridReg);
	}
	
	public static Options createOptions() {
		Options ops = new Options();
		
		ops.addOption(null, "samples", true, "Number of samples to draw");
		ops.addOption(null, "rand-seed", true, "Random seed for use when sampling");
		ops.addOption(null, "disable-preload", false, "Flag to disable pre-loading of slab hazard curves (to quickly get to the build loop)");
		ops.addOption(null, "no-maps", false, "Flag to disable combining hazard maps");
		ops.addOption(null, "no-curves", false, "Flag to disable combining site hazard curves");
		
		return ops;
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
			if (ScalarIMR_ParamsLogicTreeNode.class.isAssignableFrom(level.getType())
					&& !name.toLowerCase().contains("interface") && !name.toLowerCase().contains("slab")) {
				// need to remap
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
			} else {
				// don't need to remap
				levelRemaps.put(level, level);
				for (LogicTreeNode node : level.getNodes())
					nodeRemaps.put(node, node);
			}
		}
	}

	@Override
	protected void remapOuterTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps) {
		remapTree(tree, levelRemaps, nodeRemaps, "Interface", "");
	}

	@Override
	protected void remapInnerTree(LogicTree<?> tree, Map<LogicTreeLevel<?>, LogicTreeLevel<?>> levelRemaps,
			Map<LogicTreeNode, LogicTreeNode> nodeRemaps) {
		remapTree(tree, levelRemaps, nodeRemaps, "Intraslab", "");
	}

	@Override
	protected boolean doesOuterSupplySols() {
		return false;
	}

	@Override
	protected boolean doesInnerSupplySols() {
		return false;
	}

	@Override
	protected boolean isSerializeGridded() {
		return false;
	}
	
	private static List<LogicTreeLevel<?>> getCommonLevels(LogicTree<?> innerTree) {
		List<LogicTreeLevel<?>> commonLevels = new ArrayList<>();
		for (LogicTreeLevel<?> level : innerTree.getLevels()) {
			// everything is common except GMM branches
			if (!ScalarIMRsLogicTreeNode.class.isAssignableFrom(level.getType())
					&& !ScalarIMR_ParamsLogicTreeNode.class.isAssignableFrom(level.getType()))
				commonLevels.add(level);
		}
		return commonLevels;
	}

}
