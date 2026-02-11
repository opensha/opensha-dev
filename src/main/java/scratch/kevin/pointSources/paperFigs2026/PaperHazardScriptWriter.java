package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class PaperHazardScriptWriter {
	
	public static void main(String[] args) throws IOException {
		IncludeBackgroundOption[] bgOps = IncludeBackgroundOption.values();
		
		boolean linkFromBase = true;
		Double vs30 = null;
		Double sigmaTrunc = 3d;
		
		double[] periods = { 0d, 1d };
		AttenRelRef[] gmms = new AttenRelRef[] {AttenRelRef.USGS_NSHM23_ACTIVE};
		
		File refDir = ORIG_SOL_DIR;
		String refDirName = refDir.getName();
		String solFileName = ORIG_SOL_FILE.getName();
		String modSolFileName = "mod_input_solution.zip";
		
		Models modelToWriteInputs = null;
		String refInputsName = HAZARD_MODEL_PREFIX+"-ref_exclude";
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		FaultSystemSolution inputSol = FaultSystemSolution.load(ORIG_SOL_FILE);
		
		inputSol.removeModuleInstances(BranchRegionalMFDs.class);
		inputSol.removeModuleInstances(BranchSectBVals.class);
		inputSol.removeModuleInstances(BranchSectNuclMFDs.class);
		inputSol.removeModuleInstances(BranchSectParticMFDs.class);
		inputSol.removeModuleInstances(BranchParentSectParticMFDs.class);
		
		Integer forceMaxDispatch = 100;
		
		System.out.println("Region has "+FULL_GRID_REG.getNodeCount()+" nodes");
		System.out.println("Zoom region has "+ZOOM_GRID_REG.getNodeCount()+" nodes");
		
		File remoteMainDir = new File("/project2/scec_608/kmilner/fss_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
		int mins = 5*1440;
//		int nodes = 18;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter parallelMPJWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
		JavaShellScriptWriter singleMPJWrite = new NoMPJSingleNodeShellScriptWriter(USC_CARC_ScriptWriter.JAVA_BIN,
				remoteTotalMemGB*1024, null);
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
		singleMPJWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		parallelMPJWrite.setClasspath(classpath);
		
		List<File> singleClasspath = new ArrayList<>(classpath);
		singleClasspath.add(new File("/project2/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		singleMPJWrite.setClasspath(singleClasspath);
		
		GridSourceList origGridded = inputSol.requireModule(GridSourceList.class);
		
		for (Models model : Models.values()) {
			System.out.println("Writing "+model.name+" to "+model.getMapDir().getName());
			String gridArgz = model.getGridArgs();
			Function<GridSourceList, GridSourceList> gridFunc = model.getGridModFunction();
			System.out.println("\t"+gridArgz);
			if (gridFunc != null) {
				System.out.println("\tHas custom gridded ruptures properties");
				GridSourceList modGridded = gridFunc.apply(origGridded);
				
				inputSol.setGridSourceProvider(modGridded);
			}
			
//			boolean[] zooms = model.finiteNum < 50 ? new boolean[] {false,true} : new boolean[] {false};
			boolean[] zooms = {false, true};
			
			for (boolean zoom : zooms) {
				File localDir = zoom ? model.getZoomMapDir() : model.getMapDir();
				GriddedRegion gridReg = zoom ? ZOOM_GRID_REG : FULL_GRID_REG;
				
				Preconditions.checkState(localDir.exists() || localDir.mkdir());
				
				parallelMPJWrite.setEnvVar("DIR", mainDirPath+"/"+localDir.getName());
				singleMPJWrite.setEnvVar("DIR", mainDirPath+"/"+localDir.getName());
				
				// write the region
				File localReg = new File(localDir, "gridded_region.json");
				Feature.write(gridReg.toFeature(), localReg);
				
				String resultsPath = dirPath+"/results";
				String regPath = dirPath+"/"+localReg.getName();
				
				String solFilePath;
				
				List<String> setupLines = new ArrayList<>();
				if (gridFunc == null) {
					solFilePath = "$DIR/"+solFileName;
					setupLines.add("if [[ ! -e "+solFilePath+" ]];then");
					setupLines.add("  ln -s $MAIN_DIR/"+refDirName+"/"+solFileName+" "+solFilePath);
					setupLines.add("fi");
				} else {
					solFilePath = "$DIR/"+modSolFileName;
					
					File modSolFile = new File(localDir, modSolFileName);
					System.out.println("\t\tWriting modified solution to "+modSolFile.getAbsolutePath());
					
					// don't copy unkonwn files
					inputSol.write(ArchiveOutput.getDefaultOutput(modSolFile, inputSol.getArchive().getInput()), false);
				}
				setupLines.add("if [[ ! -e results ]];then");
				setupLines.add("  cp -r $MAIN_DIR/"+refInputsName+"/results .");
				setupLines.add("fi");
				parallelMPJWrite.setCustomSetupLines(setupLines);
				singleMPJWrite.setCustomSetupLines(setupLines);
				
				
				for (IncludeBackgroundOption bgOp : bgOps) {
					if (bgOp == IncludeBackgroundOption.EXCLUDE && model != modelToWriteInputs)
						continue;
					int myNodes;
					JavaShellScriptWriter mpjWrite;
					String dispatchArgs;
					if (bgOp == IncludeBackgroundOption.INCLUDE && bgOps.length == 3) {
						// do single node
						myNodes = 1;
						mpjWrite = singleMPJWrite;
						dispatchArgs = MPJTaskCalculator.argumentBuilder().exactDispatch(gridReg.getNodeCount()).threads(remoteTotalThreads).build();
					} else {
						// do parallel
						myNodes = nodes;
						mpjWrite = parallelMPJWrite;
						
						int maxDispatch;
						if (forceMaxDispatch != null)
							maxDispatch = forceMaxDispatch;
						else if (gridReg.getNodeCount() > 50000)
							maxDispatch = Integer.max(remoteTotalThreads*20, 1000);
						else if (gridReg.getNodeCount() > 10000)
							maxDispatch = Integer.max(remoteTotalThreads*10, 500);
						else if (gridReg.getNodeCount() > 5000)
							maxDispatch = remoteTotalThreads*5;
						else
							maxDispatch = remoteTotalThreads*3;
						dispatchArgs = MPJTaskCalculator.argumentBuilder().minDispatch(remoteTotalThreads)
								.maxDispatch(maxDispatch).threads(remoteTotalThreads).build();
					}
					
					String argz = "--input-file "+solFilePath;
					argz += " --output-dir "+resultsPath;
					argz += " --output-file "+resultsPath+"_hazard_"+bgOp.name()+".zip";
					argz += " --region "+regPath;
					
					argz += " --gridded-seis "+bgOp.name();
					argz += " "+gridArgz;
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
					if (vs30 != null)
						argz += " --vs30 "+vs30.floatValue();
					if (sigmaTrunc != null)
						argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
					argz += " "+dispatchArgs;
					
					File jobFile = new File(localDir, "batch_hazard_"+bgOp.name()+".slurm");
					
					List<String> script = mpjWrite.buildScript(MPJ_SingleSolHazardCalc.class.getName(), argz);
					
//					System.out.println("\t\tWriting "+jobFile.getAbsolutePath());
					
					pbsWrite.writeScript(jobFile, script, mins, myNodes, remoteTotalThreads, queue);
				}
			}
		}
		
		System.out.println("DONE");
	}

}
