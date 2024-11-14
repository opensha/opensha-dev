package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.HovenweepScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class BranchAveragedHazardScriptWriter {

	public static void main(String[] args) throws IOException {
		
		IncludeBackgroundOption[] bgOps = IncludeBackgroundOption.values();
		
		boolean linkFromBase = true;
		Double vs30 = null;
		double gridSpacing = 0.1;
		
		double[] periods = { 0d, 0.2d, 1d, 5d };
		AttenRelRef[] gmms = null;
		
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		/*
		 * NSHM23
		 */
//		String baseDirName = "2024_02_02-nshm23_branches-WUS_FM_v3";
//		String baseDirName = "2023_11_20-nshm23_branches-dm_sampling-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2023_11_17-nshm23_branches-dm_sampling-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2023_11_16-nshm23_branches-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2024_05_07-nshm23_branches-WUS_FM_v3-AvgSupraB-AvgSeg";
		
//		String suffix = "true_mean";
//		String solFileName = "true_mean_solution.zip";
		
//		String suffix = "ba_only";
//		String solFileName = "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
////		String solFileName = "results_WUS_FM_v3_branch_averaged_gridded.zip";
		
//		String suffix = "ba_only-mod_gridded";
//		String solFileName = "results_WUS_FM_v3_branch_averaged_mod_gridded.zip";
		
		/*
		 * PRVI
		 */
////		String baseDirName = "2024_07_31-prvi25_subduction_branches";
////		String baseDirName = "2024_09_04-prvi25_crustal_subduction_combined_branches";
////		String baseDirName = "2024_09_04-prvi25_crustal_branches-dmSample5x";
//		String baseDirName = "2024_09_04-prvi25_subduction_branches";
//		region = PRVI25_RegionLoader.loadPRVI_ModelBroad();
//		gridSpacing = 0.01;
//		
//		gmms = new AttenRelRef[] { AttenRelRef.USGS_PRVI_ACTIVE, AttenRelRef.USGS_PRVI_INTERFACE, AttenRelRef.USGS_PRVI_SLAB };
//		
////		String suffix = "ba_only-LARGE";
////		String solFileName = "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip";
//		
////		String suffix = "ba_only-LARGE-true_pt_src";
////		String solFileName = "results_PRVI_SUB_FM_LARGE_branch_averaged_gridded.zip";
//		
////		String suffix = "ba_only";
////		String solFileName = "combined_branch_averaged_solution.zip";
//		
////		String suffix = "ba_only";
////		String solFileName = "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip";
//		
////		String suffix = "ba_only-SLAB_only";
////		String solFileName = "results_PRVI_SLAB_ONLY_branch_averaged_gridded.zip";
////		bgOps = new IncludeBackgroundOption[] { IncludeBackgroundOption.ONLY };
//		
////		String suffix = "ba_only-INTERFACE_only";
////		String solFileName = "results_PRVI_INTERFACE_ONLY_branch_averaged_gridded.zip";
//		
//		String suffix = "ba_only-both_fms";
//		String solFileName = "results_PRVI_SUB_FMs_combined_branch_averaged_gridded.zip";
		
		/*
		 * RSQSim
		 */
		String suffix = null;
		
		String baseDirName = "2024_11_12-rsqsim-wus-5895";
		String solFileName = "fss_m6_skip10000_sectArea0.5_minSubSects2.zip";
		
//		String baseDirName = "2024_11_12-rsqsim-wus-5892";
//		String solFileName = "fss_m6_skip20000_sectArea0.5_minSubSects2.zip";
		linkFromBase = false;
		
//		periods = new double[] {0d, 1d, 2d, 3d, 5d};
//		vs30 = 500d; suffix = "vs30_500";
		periods = new double[] {2d, 3d, 5d};
		gridSpacing = 0.5; suffix = "bbp_"+(float)gridSpacing;
		
		boolean noMFDs = false;
		
		GriddedRegion gridReg = new GriddedRegion(
				region, gridSpacing, GriddedRegion.ANCHOR_0_0);
//		GriddedRegion gridReg = new GriddedRegion(
//				PRVI25_RegionLoader.loadPRVI_ModelBroad(), gridSpacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		String dirName = baseDirName;
		if (suffix != null && !suffix.isBlank())
			dirName += "-"+suffix;
		if (noMFDs)
			dirName += "-no_mfds";
		
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
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
		
		parallelMPJWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		singleMPJWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		parallelMPJWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		singleMPJWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		parallelMPJWrite.setClasspath(classpath);
		
		List<File> singleClasspath = new ArrayList<>(classpath);
		singleClasspath.add(new File("/project/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		singleMPJWrite.setClasspath(singleClasspath);

		// write the region
		File localReg = new File(localDir, "gridded_region.json");
		Feature.write(gridReg.toFeature(), localReg);
		
		String resultsPath = dirPath+"/results";
		String regPath = dirPath+"/"+localReg.getName();
		
		String solFilePath = "$DIR/"+solFileName;
		
		List<String> setupLines = new ArrayList<>();
		if (linkFromBase) {
			setupLines.add("if [[ ! -e "+solFilePath+" ]];then");
			setupLines.add("  ln -s $MAIN_DIR/"+baseDirName+"/"+solFileName+" "+solFilePath);
			setupLines.add("fi");
		}
		parallelMPJWrite.setCustomSetupLines(setupLines);
		singleMPJWrite.setCustomSetupLines(setupLines);
		
		
		for (IncludeBackgroundOption bgOp : bgOps) {
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
				if (gridReg.getNodeCount() > 50000)
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
			
			String argz = "--input-file "+dirPath+"/"+solFileName;
			argz += " --output-dir "+resultsPath;
			argz += " --output-file "+resultsPath+"_hazard_"+bgOp.name()+".zip";
			argz += " --region "+regPath;
			if (noMFDs)
				argz += " --no-mfds";
			if (vs30 != null)
				argz += " --vs30 "+vs30.floatValue();
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
			argz += " "+dispatchArgs;
			
			File jobFile = new File(localDir, "batch_hazard_"+bgOp.name()+".slurm");
			
			List<String> script = mpjWrite.buildScript(MPJ_SingleSolHazardCalc.class.getName(), argz);
			
			System.out.println("Writing "+jobFile.getAbsolutePath());
			
			pbsWrite.writeScript(jobFile, script, mins, myNodes, remoteTotalThreads, queue);
		}
	}

}
