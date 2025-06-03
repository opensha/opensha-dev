package scratch.kevin.prvi25;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.CoulombRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupSetTectonicRegimes;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.SolutionLogicTreeCombinationProcessor;
import org.opensha.sha.earthquake.faultSysSolution.treeCombiners.SolutionLogicTreeCombinationProcessor.CombinedRupSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolModuleStripper;
import org.opensha.sha.earthquake.faultSysSolution.util.TrueMeanSolutionCreator;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_MaxMagOffFault;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class BowinFaultAddTest {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions");

		File subBADir = new File(invsDir, "2025_01_17-prvi25_subduction_branches");

//		Double customSlipRate = null;
//		File bowinTestWithDir = new File(invsDir, "2025_05_08-prvi25-bowin_test_with");
		Double customSlipRate = 5d;
		File bowinTestWithDir = new File(invsDir, "2025_05_08-prvi25-bowin_test_with-5mm");
		File bowinTestWithoutDir = new File(invsDir, "2025_05_08-prvi25-bowin_test_without");
		
		Feature feature = FeatureCollection.read(new File("/home/kevin/Downloads/westBowin_v1.geojson")).features.get(0);
		feature.properties.set(GeoJSONFaultSection.SLIP_RATE,
				customSlipRate == null ? feature.properties.get("PrefRate") : customSlipRate);
		GeoJSONFaultSection bowinSect = GeoJSONFaultSection.fromFeature(feature);
		bowinSect.setSlipRateStdDev(PRVI25_CrustalDeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND
				*bowinSect.getOrigAveSlipRate());
		System.out.println("Parsed GeoJSON:\n"+bowinSect.toFeature().toJSON());
		Preconditions.checkState(bowinSect.getOrigAveSlipRate() > 0d);
		Preconditions.checkState(bowinSect.getOrigDownDipWidth() > 0d);
		
		double minDist = Double.POSITIVE_INFINITY;
		for (XY_DataSet outline : PoliticalBoundariesData.loadUSState("Puerto Rico")) {
			for (Point2D pt : outline) {
				Location outlineLoc = new Location(pt.getY(), pt.getX());
				for (Location traceLoc : bowinSect.getFaultTrace()) {
					minDist = Math.min(LocationUtils.horzDistanceFast(traceLoc, outlineLoc), minDist);
				}
			}
		}
		System.out.println("Bowin is "+(float)minDist+" km away from Puerto Rico");
		
		System.exit(0);

		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();

		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));

		int threads = FaultSysTools.defaultNumThreads();
		
		boolean rebuild = false;

		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		levels.addAll(PRVI25_LogicTreeBranch.levelsOnFault);
		levels.addAll(PRVI25_LogicTreeBranch.levelsCrustalOffFault);
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		
		// start with on-fault defaults
		for (LogicTreeNode node : PRVI25_LogicTreeBranch.DEFAULT_CRUSTAL_ON_FAULT)
			branch.setValue(node);
		
		PRVI25_CrustalFaultModels crustalFM = branch.requireValue(PRVI25_CrustalFaultModels.class); 

		// set on-fault average values
		branch.setValue(NSHM23_SegmentationModels.AVERAGE);
		branch.setValue(SupraSeisBValues.AVERAGE);
		branch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		branch.setValue(PRVI25_CrustalDeformationModels.GEOLOGIC_DIST_AVG);
		
		// set off-fault averages
		branch.setValue(PRVI25_CrustalSeismicityRate.AVERAGE);
		branch.setValue(PRVI25_DeclusteringAlgorithms.AVERAGE);
		branch.setValue(PRVI25_SeisSmoothingAlgorithms.AVERAGE);
		branch.setValue(NSHM23_MaxMagOffFault.MAG_7p6);
		
		Map<PRVI25_SubductionFaultModels, FaultSystemSolution> subductionBASols = new HashMap<>();
		subductionBASols.put(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE,
				FaultSystemSolution.load(new File(subBADir, "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip")));
		subductionBASols.put(PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL,
				FaultSystemSolution.load(new File(subBADir, "results_PRVI_SUB_FM_SMALL_branch_averaged_gridded.zip")));
		
		List<LogicTreeLevel<? extends LogicTreeNode>> combLevels = new ArrayList<>();
		List<LogicTreeBranch<LogicTreeNode>> combBranches = new ArrayList<>();
		
		combLevels.add(PRVI25_LogicTreeBranch.CRUSTAL_FM);
		combLevels.add(PRVI25_LogicTreeBranch.SUB_FM);
		
		for (PRVI25_SubductionFaultModels subductionFM : subductionBASols.keySet())
			combBranches.add(new LogicTreeBranch<>(combLevels, List.of(crustalFM, subductionFM)));
		
		LogicTree<?> tree = LogicTree.fromExisting(combLevels, combBranches);
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
		int mins = 600;
//		int nodes = 18;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter parallelMPJWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
//		File remoteMainDir = new File("/caldera/hovenweep/projects/usgs/hazards/ehp/kmilner/nshm23/batch_inversions");
//		int remoteTotalThreads = 128;
//		int remoteTotalMemGB = 448;
//		String queue = null;
//		int nodes = 4;
//		int mins = 180;
////		int nodes = 18;
////		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
////				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
//		JavaShellScriptWriter parallelMPJWrite = new FastMPJShellScriptWriter(
//				HovenweepScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, HovenweepScriptWriter.FMPJ_HOME);
//		JavaShellScriptWriter singleMPJWrite = new NoMPJSingleNodeShellScriptWriter(HovenweepScriptWriter.JAVA_BIN,
//				remoteTotalMemGB*1024, null);
//		BatchScriptWriter pbsWrite = new HovenweepScriptWriter();
		
		if (parallelMPJWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)parallelMPJWrite).setUseLaunchWrapper(true);
		
		parallelMPJWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		parallelMPJWrite.setClasspath(classpath);

		for (boolean with : new boolean[] {false, true}) {
			File dir = with ? bowinTestWithDir : bowinTestWithoutDir;
			Preconditions.checkState(dir.exists() || dir.mkdir());
			
			File rupSetFile = new File(dir, "rup_set.zip");
			File solutionFile = new File(dir, "solution.zip");
			File trueMeanFile = new File(dir, "solution_comb_with_subduction.zip");

			FaultSystemRupSet rupSet;
			
			if (!rebuild && rupSetFile.exists()) {
				rupSet = FaultSystemRupSet.load(rupSetFile);
			} else {
				if (solutionFile.exists())
					solutionFile.delete();
				if (with) {
					// need to build it custom
					List<FaultSection> subSects = new ArrayList<>(branch.requireValue(PRVI25_CrustalDeformationModels.class).build(branch));
					
					List<GeoJSONFaultSection> bowinSubSects = bowinSect.getSubSectionsList(0.5*bowinSect.getOrigDownDipWidth(), subSects.size(), 2);
					subSects.addAll(bowinSubSects);
					
					RupSetConfig config = RupturePlausibilityModels.COULOMB.getConfig(subSects, branch.requireValue(NSHM23_ScalingRelationships.class));
					((CoulombRupSetConfig)config).setConnectProxyFaults(PRVI25_InvConfigFactory.ALLOW_CONNECTED_PROXY_FAULTS);
					rupSet = config.build(threads);
					factory.getSolutionLogicTreeProcessor().processRupSet(rupSet, branch);
				} else {
					rupSet = factory.buildRuptureSet(branch, threads);
				}

				rupSet.write(rupSetFile);

				System.out.println("Built rupture set "+(with?"with":"without")+" Bowin fault: "+rupSet.getNumRuptures()+" ruptures");
			}
			
			FaultSystemSolution solution;
			GridSourceList crustalGridded;
			if (!rebuild && solutionFile.exists()) {
				solution = FaultSystemSolution.load(solutionFile);
				crustalGridded = solution.requireModule(GridSourceList.class);
			} else {
				if (trueMeanFile.exists())
					trueMeanFile.delete();
				solution = Inversions.run(rupSet, factory, branch, threads, null);
				
				// now add gridded seismicity
				crustalGridded = PRVI25_GridSourceBuilder.buildCrustalGridSourceProv(solution, branch);
				solution.setGridSourceProvider(crustalGridded);
				
				solution.write(solutionFile);
			}
			
			Preconditions.checkState(solution.getRupSet().hasModule(RupSetTectonicRegimes.class), "Crustal solution doesn't have TRTs");
			
			if (rebuild || !trueMeanFile.exists()) {
				// simplify crustal to distribute proxies
				solution = SolModuleStripper.stripModules(solution, 5d, true, false);
				
				// now combine
				TrueMeanSolutionCreator creator = new TrueMeanSolutionCreator(tree);
				creator.setDoGridProv(true);
				for (LogicTreeBranch<?> combBranch : tree) {
					FaultSystemSolution subductionSol = subductionBASols.get(combBranch.requireValue(PRVI25_SubductionFaultModels.class));
					Preconditions.checkState(subductionSol.getRupSet().hasModule(RupSetTectonicRegimes.class), "Subduction solution doesn't have TRTs");
					
					FaultSystemSolution combined = SolutionLogicTreeCombinationProcessor.combineSols(solution, subductionSol, true);
					Preconditions.checkState(combined.getRupSet().hasModule(RupSetTectonicRegimes.class), "Combined solution doesn't have TRTs");
					GridSourceList subductionGridded = subductionSol.requireModule(GridSourceList.class);
					CombinedRupSetMappings mappings = combined.getRupSet().requireModule(CombinedRupSetMappings.class);
					crustalGridded = GridSourceList.remapAssociations(crustalGridded, mappings.getInnerSectMappings());
					subductionGridded = GridSourceList.remapAssociations(subductionGridded, mappings.getOuterSectMappings());
					combined.setGridSourceProvider(GridSourceList.combine(subductionGridded, crustalGridded));
					
					creator.addSolution(combined, combBranch);
				}
				
				FaultSystemSolution trueMean = creator.build();
				Preconditions.checkState(trueMean.getRupSet().hasModule(RupSetTectonicRegimes.class), "True mean solution doesn't have TRTs");
				trueMean.write(trueMeanFile);
			}
			
			parallelMPJWrite.setEnvVar("DIR", mainDirPath+"/"+dir.getName());
			
			Region region = PRVI25_RegionLoader.loadPRVI_IntermediateModelMapExtents();
			double gridSpacing = 0.05;
			
			AttenRelRef[] gmms = new AttenRelRef[] { AttenRelRef.USGS_PRVI_ACTIVE, AttenRelRef.USGS_PRVI_INTERFACE, AttenRelRef.USGS_PRVI_SLAB };
			double[] periods = new double[] { 0d, 0.2d, 1d, 5d };
			boolean supersample = true;
			Double sigmaTrunc = 3d;
			IncludeBackgroundOption bgOp = IncludeBackgroundOption.INCLUDE;
			
			GriddedRegion gridReg = new org.opensha.commons.geo.GriddedRegion(region, gridSpacing, GriddedRegion.ANCHOR_0_0);
			
			// write the region
			File localReg = new File(dir, "gridded_region.json");
			Feature.write(gridReg.toFeature(), localReg);
			
			String resultsPath = dirPath+"/results";
			String regPath = dirPath+"/"+localReg.getName();
			
			String argz = "--input-file "+dirPath+"/"+trueMeanFile.getName();
			argz += " --output-dir "+resultsPath;
			argz += " --output-file "+resultsPath+"_hazard_"+bgOp.name()+".zip";
			argz += " --region "+regPath;
			
			argz += " --gridded-seis "+bgOp.name();
			if (gmms != null)
				for (AttenRelRef gmm : gmms)
					argz += " --gmpe "+gmm.name();
			if (periods != null) {
				argz += " --periods ";
				for (int p=0; p<periods.length; p++) {
					if (p > 0)
						argz += ",";
					argz += (float)periods[p];
				}
			}
			if (supersample)
				argz += " --supersample";
			if (sigmaTrunc != null)
				argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
			argz += " "+MPJTaskCalculator.argumentBuilder().minDispatch(remoteTotalThreads)
					.maxDispatch(500).threads(remoteTotalThreads).build();
			
			File jobFile = new File(dir, "batch_hazard_"+bgOp.name()+".slurm");
			
			List<String> script = parallelMPJWrite.buildScript(MPJ_SingleSolHazardCalc.class.getName(), argz);
			
			System.out.println("Writing "+jobFile.getAbsolutePath());
			
			pbsWrite.writeScript(jobFile, script, mins, nodes, remoteTotalThreads, queue);
		}
	}

}
