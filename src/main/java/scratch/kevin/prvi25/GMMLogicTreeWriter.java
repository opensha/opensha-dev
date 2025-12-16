package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.FileBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.FileBackedNode;
import org.opensha.commons.logicTree.treeCombiner.AbstractLogicTreeCombiner;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SiteLogicTreeHazardCurveCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class GMMLogicTreeWriter {

	public static void main(String[] args) throws IOException {
		String jobSuffix = "";
		String outputSuffix = "";
		String dirSuffix = "";
		String forceInputFileName = null;
		String logicTreeOutputName = null;
		boolean combineOnly = false;
		
//		GriddedRegion gridReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_ModelBroad(), 0.1, GriddedRegion.ANCHOR_0_0);
//		GriddedRegion gridReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_MapExtents(), 0.025, GriddedRegion.ANCHOR_0_0);
		// use for rough ERF+GMC
//		GriddedRegion gridReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_Tight(), 0.05, GriddedRegion.ANCHOR_0_0);
		// use for GMC and final ERF+GMC
		GriddedRegion gridReg = new GriddedRegion(PRVI25_RegionLoader.loadPRVI_Tight(), 0.025, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		Double vs30 = null;
		
//		vs30 = 760d; dirSuffix = "-vs760";
		vs30 = 260d; dirSuffix = "-vs260";
		Double sigmaTrunc = 3d;
		boolean supersample = true;
		int erfSamples = -1;
		int gmmSamplesPerERF = -1;
		double[] periods = { 0d, 0.2d, 1d, 5d };
		
		/*
		 * Active crustal
		 *
		 * do supra-seis and then include
		 * 
		 * fetch logic trees first
		 */
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsCrustalGMM;
//		File sourceDir = CRUSTAL_DIR;
//		File outputDir = new File(sourceDir.getParentFile(), sourceDir.getName()+"-gmTreeCalcs"+dirSuffix);
//		// supra-seis only
////		File sourceTreeFile = new File(sourceDir, "logic_tree.json");
////		int mins = 1440*5;
////		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
//		// including gridded
//		int mins = 1440*5;
//		File sourceTreeFile = new File(sourceDir, "logic_tree_full_gridded.json");
//		erfSamples = 20000; gmmSamplesPerERF = 1; jobSuffix = "_sampled"; logicTreeOutputName = "logic_tree_full_gridded_sampled.json";
////		File sourceTreeFile = new File(sourceDir, "logic_tree_full_gridded_sampled.json"); jobSuffix = "_sampled";
//		IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
		
		/*
		 * Interface separate slab and interface
		 * 
		 * do supra-seis, then each gridded-only, then interface combine
		 * 
		 * then need to separately combine the slab logic tree
		 * 
		 * fetch logic trees first
		 */
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsInterfaceGMM;
//		File sourceDir = SUBDUCTION_DIR;
//		File outputDir = new File(sourceDir.getParentFile(), sourceDir.getName()+"-gmTreeCalcs"+dirSuffix);
//		// supra-seis only
//		File sourceTreeFile = new File(sourceDir, "logic_tree.json");
//		int mins = 1440*4;
//		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
//		// interface gridded only
////		int mins = 1440*4;
//////		File sourceTreeFile = new File(sourceDir, "logic_tree_gridded_only.json");
////		File sourceTreeFile = new File(sourceDir, "logic_tree_full_gridded_for_only_calc.json");
////		logicTreeOutputName = "logic_tree_gridded_interface_only.json";
////		IncludeBackgroundOption bgOp = IncludeBackgroundOption.ONLY;
////		// this was for if gridded only depended on FM but it also depends on scale
//////		forceInputFileName = "results_gridded_branches_interface_only.zip";
////		// use this one because it has scaling relationship specific gridded models
////		forceInputFileName = "results_full_gridded_interface_only.zip";
////		jobSuffix = "_interface";
////		outputSuffix = jobSuffix;
//		// interface both (combine only)
////		combineOnly = true;
////		int mins = 1440*4;
////		forceInputFileName = "results_full_gridded_interface_only.zip";
////		File sourceTreeFile = new File(sourceDir, "logic_tree_full_gridded.json");
////		logicTreeOutputName = "logic_tree_full_gridded_interface_only.json";
////		IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
////		jobSuffix = "_interface";
////		outputSuffix = jobSuffix;
		
		/*
		 * Slab
		 * 
		 * the logic_tree_gridded_only.json file you need is generated by subduction_slt_split.slurm
		 */
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsSlabGMM;
//		File sourceDir = SUBDUCTION_DIR;
//		File outputDir = new File(sourceDir.getParentFile(), sourceDir.getName()+"-gmTreeCalcs"+dirSuffix);
//		// always slab gridded only
//		int mins = 1440*4;
//		File sourceTreeFile = new File(sourceDir, "logic_tree_gridded_only.json");
//		logicTreeOutputName = "logic_tree_gridded_slab_only.json";
//		IncludeBackgroundOption bgOp = IncludeBackgroundOption.ONLY;
//		forceInputFileName = "results_gridded_branches_slab_only.zip";
//		jobSuffix = "_slab";
//		outputSuffix = jobSuffix;
		
		/*
		 * Branch averaged (GMC-only)
		 */
		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsCombinedGMM;
		File sourceDir = COMBINED_DIR;
//		File sourceDir = new File(INV_DIR, "2025_01_02-prvi25_crustal_subduction_combined_branches");
		File outputDir = new File(sourceDir.getParentFile(), sourceDir.getName()+"-ba_only-gmTreeCalcs"+dirSuffix);
		// write out a SLT that only contains that node
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File sourceTreeFile = new File(outputDir, "fake_erf_logic_tree.json");
		FileBackedLevel fakeLevel = new FileBackedLevel("ERF Model", "ERF",
				List.of(new FileBackedNode("Branch Averaged ERF", "BranchAveragedERF", 1d, "BA_ERF")));
		LogicTree<?> tempTree = LogicTree.buildExhaustive(List.of(fakeLevel), true);
		Preconditions.checkState(tempTree.size() == 1);
		File sourceFile = new File(outputDir, "fake_erf_slt.zip");
		SolutionLogicTree.FileBuilder builder = new SolutionLogicTree.FileBuilder(sourceFile);
		builder.setSerializeGridded(true);
		builder.solution(FaultSystemSolution.load(new File(sourceDir, COMBINED_SOL.getName())), tempTree.getBranch(0));
		builder.close();
		forceInputFileName = sourceFile.getName();
		tempTree.write(sourceTreeFile);
		logicTreeOutputName = "logic_tree.json";
		sourceDir = outputDir;
		int mins = 1440;
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
		
		
		// FOR ALL
		System.out.println("Output dir: "+outputDir.getName());
		
		LogicTree<?> erfTree = LogicTree.read(sourceTreeFile);
		System.out.println("Read "+erfTree.size()+" ERF branches");
		
		LogicTree<?> gmmTree = LogicTree.buildExhaustive(gmmLevels, true);
		System.out.println("Built "+gmmTree.size()+" GMM branches");
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File gridRegFile = new File(outputDir, "gridded_region.geojson");
		Feature.write(gridReg.toFeature(), gridRegFile);
		
		LogicTree<?> logicTree;
		if (gmmSamplesPerERF > 0) {
			System.out.println("Pairwise sampling with pairwise="+gmmSamplesPerERF+" branches");
			if (erfSamples <= 0)
				erfSamples = erfTree.size();
			logicTree = AbstractLogicTreeCombiner.pairwiseSampleLogicTrees(erfTree, gmmTree, erfSamples, gmmSamplesPerERF);
		} else {
			List<LogicTreeLevel<? extends LogicTreeNode>> combLevels = new ArrayList<>();
			combLevels.addAll(erfTree.getLevels());
			combLevels.addAll(gmmLevels);
			
			List<LogicTreeBranch<LogicTreeNode>> combBranches = new ArrayList<>(erfTree.size()*gmmTree.size());
			
			for (LogicTreeBranch<?> branch : erfTree) {
				for (LogicTreeBranch<?> gmmBranch : gmmTree) {
					LogicTreeBranch<LogicTreeNode> combBranch = new LogicTreeBranch<>(combLevels);
					for (LogicTreeNode node : branch)
						combBranch.setValue(node);
					for (LogicTreeNode node : gmmBranch)
						combBranch.setValue(node);
					combBranches.add(combBranch);
				}
			}
			
			logicTree = LogicTree.fromExisting(combLevels, combBranches);
		}
		
		File localLogicTree = new File(outputDir, logicTreeOutputName == null ? sourceTreeFile.getName() : logicTreeOutputName);
		logicTree.write(localLogicTree);
		System.out.println("Wrote "+logicTree.size()+" branches");
		
		String dirName = outputDir.getName();
		File localDir = outputDir;
		
		File remoteMainDir = new File("/project2/scec_608/kmilner/fss_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
//		JavaShellScriptWriter mpjWrite = new NoMPJSingleNodeShellScriptWriter(USC_CARC_ScriptWriter.JAVA_BIN,
//				remoteTotalMemGB*1024, null); nodes = 1; remoteInversionsPerBundle = 2;
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		mpjWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		mpjWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		mpjWrite.setEnvVar("ORIG_DIR", mainDirPath+"/"+sourceDir.getName());
		String origDirPath = "$ORIG_DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		if (mpjWrite instanceof NoMPJSingleNodeShellScriptWriter)
			classpath.add(new File("/project2/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		
		mpjWrite.setClasspath(classpath);
		if (mpjWrite instanceof MPJExpressShellScriptWriter)
			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		else if (mpjWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		String mapScriptName = "batch_hazard_"+bgOp.name()+jobSuffix+".slurm";
		String siteScriptName = "batch_hazard_sites_"+bgOp.name()+jobSuffix+".slurm"; 
		
		// now write hazard map script
		String inputFileName;
		if (forceInputFileName != null)
			inputFileName = origDirPath+"/"+forceInputFileName; // use supplied
		else if (bgOp == IncludeBackgroundOption.EXCLUDE)
			inputFileName = origDirPath+"/results.zip"; // use the zip file
		else
			inputFileName = origDirPath+"/results"; // use the results dir
		String argz = "--input-file "+inputFileName;
		
		String outputDirName, outputFilePrefix;
		if (bgOp == IncludeBackgroundOption.INCLUDE) {
			argz += " --combine-with-dir "+dirPath+"/results";
			argz += " --combine-with-dir "+dirPath+"/results_gridded_only"+outputSuffix;
			outputDirName = "results_full_gridded";
			outputFilePrefix = outputDirName;
		} else if (bgOp == IncludeBackgroundOption.ONLY) {
			outputDirName = "results_gridded_only";
			outputFilePrefix = outputDirName;
		} else {
			outputDirName = "results";
			outputFilePrefix = "results_hazard";
		}
		argz += " --output-dir "+dirPath+"/"+outputDirName+outputSuffix;
		argz += " --output-file "+dirPath+"/"+outputFilePrefix+outputSuffix+".zip";
		argz += " --gridded-seis "+bgOp.name();
		String logicTreePath = dirPath+"/"+localLogicTree.getName();
		argz += " --logic-tree "+logicTreePath;
		if (!inputFileName.contains("interface_only") && (bgOp == IncludeBackgroundOption.ONLY || bgOp == IncludeBackgroundOption.INCLUDE))
			argz += " --quick-grid-calc";
		if (combineOnly)
			argz += " --combine-only";
		argz += " --region "+dirPath+"/"+gridRegFile.getName();
		if (vs30 != null)
			argz += " --vs30 "+vs30.floatValue();
		if (supersample)
			argz += " --supersample-quick";
		if (sigmaTrunc != null)
			argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
		if (periods != null) {
			argz += " --periods ";
			for (int p=0; p<periods.length; p++) {
				if (p > 0)
					argz += ",";
				argz += (float)periods[p];
			}
		}
		argz += " "+MPJTaskCalculator.argumentBuilder().maxDispatch(100).threads(remoteTotalThreads).build();
		List<String> script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
		pbsWrite.writeScript(new File(localDir, mapScriptName), script, mins, nodes, remoteTotalThreads, queue);
		
		// now write hazard curve script
		CSVFile<String> csv = CSVFile.readStream(PRVI25_CrustalFaultModels.class.getResourceAsStream("/data/erf/prvi25/sites/prvi_sites.csv"), true);
		List<Site> sites = new ArrayList<>();
		for (int row=1; row<csv.getNumRows(); row++) {
			String name = csv.get(row, 0);
			Location loc = new Location(csv.getDouble(row, 2), csv.getDouble(row, 1));
			sites.add(new Site(loc, name));
		}
		csv = new CSVFile<>(true);
		csv.addLine("Name", "Latitude", "Longitude");
		for (Site site : sites)
			csv.addLine(site.getName(), site.getLocation().lat+"", site.getLocation().lon+"");
		File localSitesFile = new File(localDir, "hazard_sites.csv");
		csv.writeToFile(localSitesFile);
		
		argz = "--input-file "+inputFileName;
		argz += " --logic-tree "+logicTreePath;
		argz += " --output-dir "+dirPath+"/"+outputFilePrefix+"_sites"+outputSuffix;
		argz += " --sites-file "+dirPath+"/"+localSitesFile.getName();
		argz += " --gridded-seis "+bgOp.name();
		if (vs30 != null)
			argz += " --vs30 "+vs30.floatValue();
		if (supersample)
			argz += " --supersample-quick";
		if (sigmaTrunc != null)
			argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
		if (periods != null) {
			argz += " --periods ";
			for (int p=0; p<periods.length; p++) {
				if (p > 0)
					argz += ",";
				argz += (float)periods[p];
			}
		}
		argz += " "+MPJTaskCalculator.argumentBuilder().minDispatch(1).maxDispatch(10).threads(remoteTotalThreads).build();
		script = mpjWrite.buildScript(MPJ_SiteLogicTreeHazardCurveCalc.class.getName(), argz);
		pbsWrite.writeScript(new File(localDir, siteScriptName), script, mins, nodes, remoteTotalThreads, queue);
		
//		combTree = LogicTree.read(new File(outputDir, "logic_tree.json"));
//		
//		LogicTreeBranch<?> prevBranch = null;
//		for (int i=0; i<combTree.size(); i++) {
//			LogicTreeBranch<?> branch = combTree.getBranch(i);
//			if (prevBranch != null) {
//				boolean onlyGmmDifferent = true;
//				for (int l=0; l<branch.size(); l++) {
//					if (!branch.getValue(l).equals(prevBranch.getValue(l))) {
//						if (!isGMMLevel(branch.getLevel(l))) {
//							onlyGmmDifferent = false;
//							break;
//						}
//					}
//				}
//				if (onlyGmmDifferent)
//					System.out.println("Only GMM branches differ for "+i+" (relative to "+(i-1)+")");
//			}
//			prevBranch = branch;
//		}
	}
	
//	private static boolean isGMMLevel(LogicTreeLevel<?> level) {
//		return ScalarIMRsLogicTreeNode.class.isAssignableFrom(level.getType()) ||
//				ScalarIMR_ParamsLogicTreeNode.class.isAssignableFrom(level.getType());
//	}

}
