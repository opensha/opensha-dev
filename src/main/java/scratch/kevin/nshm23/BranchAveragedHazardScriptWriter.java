package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.HovenweepScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class BranchAveragedHazardScriptWriter {

	public static void main(String[] args) throws IOException {
		
		IncludeBackgroundOption[] bgOps = IncludeBackgroundOption.values();
		
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
		
		boolean linkFromBase = true;
		Double vs30 = null;
		Double sigmaTrunc = null;
		double gridSpacing = 0.1;
		boolean supersample = false;
		boolean supersampleFinite = false;
		
		double[] periods = { 0d, 0.2d, 1d, 5d };
		AttenRelRef[] gmms = null;
		
		PointSourceDistanceCorrections distCorr = null;
		BackgroundRupType bgRupType = null;
		Integer bgFiniteNum = null;
		String extraGridArgs = null;
		Double pointFiniteMinMag = null;
		
		Integer forceMaxDispatch = null;
		
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		String forceOutputName = null;
		
		/*
		 * NSHM23
		 */
		String baseDirName = "2024_02_02-nshm23_branches-WUS_FM_v3";
//		String baseDirName = "2023_11_20-nshm23_branches-dm_sampling-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2023_11_17-nshm23_branches-dm_sampling-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2023_11_16-nshm23_branches-randB-randSeg-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR";
//		String baseDirName = "2024_05_07-nshm23_branches-WUS_FM_v3-AvgSupraB-AvgSeg";
		
//		String suffix = "true_mean";
//		String solFileName = "true_mean_solution.zip";
		
//		String suffix = "ba_only";
////		String solFileName = "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
//		String solFileName = "results_WUS_FM_v3_branch_averaged_gridded.zip";
		
//		String suffix = "ba_only-mod_gridded";
//		String solFileName = "results_WUS_FM_v3_branch_averaged_mod_gridded.zip";
		
		forceOutputName = "2025_04_01-nshm23_pt_src_tests";
		String solFileName = "results_WUS_FM_v3_branch_averaged_gridded.zip";
		gmms = new AttenRelRef[] {AttenRelRef.USGS_NSHM23_ACTIVE};
//		gmms = new AttenRelRef[] {AttenRelRef.NGAWest_2014_AVG_NOIDRISS};
		region = new Region(new Location(36, -120), new Location(39, -117)); gridSpacing = 0.02; forceOutputName += "-zoom";
//		region = new Region(new Location(37, -119), new Location(37.6, -118.4)); gridSpacing = 0.01; forceOutputName += "-tiny_zoom";
		sigmaTrunc = 3d;
		supersample = true;
		supersampleFinite = true; // only applies if supersample == true
		pointFiniteMinMag = 5d;
		forceMaxDispatch = 100;
		if (supersample)
			forceOutputName += "-supersample";
		
//		String suffix = "pure_pt_src";
//		distCorr = PointSourceDistanceCorrections.NONE;
//		bgRupType = BackgroundRupType.POINT;
		
//		String suffix = "pt_src_corr-2013";
//		distCorr = PointSourceDistanceCorrections.NSHM_2013;
//		bgRupType = BackgroundRupType.POINT;
//		pointFiniteMinMag = 6d;
		
//		String suffix = "finite_single";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 1;
		
//		String suffix = "finite_crosshair";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 2;
		
//		String suffix = "finite_quad_crosshair";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 4;
//		mins = 24*60;
		
//		String suffix = "finite_quad_crosshair-sample_along";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 4;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
//		mins = 24*60;
		
//		String suffix = "finite_oct_crosshair";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 8;
//		mins = 24*60;
		
//		String suffix = "finite_dodec_crosshair";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 12;
//		mins = 2*24*60;
		
		String suffix = "finite_dodec_crosshair-sample_along";
		bgRupType = BackgroundRupType.FINITE;
		bgFiniteNum = 12;
		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
		mins = 24*60;
		
//		String suffix = "finite_dodec_crosshair_to_m6";
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 12;
//		extraGridArgs = "--point-finite-min-mag 6";
		
//		String suffix = "pt_src_corr-5pt";
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_RJB_DIST;
//		bgRupType = BackgroundRupType.POINT;
		
//		String suffix = "pt_src_corr-5pt-approx_ss";
//		distCorr = PointSourceDistanceCorrections.SUPERSAMPLING_0p1_FIVE_POINT_RJB_DIST;
//		bgRupType = BackgroundRupType.POINT;
//		Preconditions.checkState(!supersample);
		
//		String suffix = "pt_src_corr-5pt-sample_along";
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_RJB_DIST_ALONG;
//		bgRupType = BackgroundRupType.POINT;
		
//		String suffix = "pt_src_corr-5pt-approx_ss-sample_along";
//		distCorr = PointSourceDistanceCorrections.SUPERSAMPLING_0p1_FIVE_POINT_RJB_DIST_ALONG;
//		bgRupType = BackgroundRupType.POINT;
//		Preconditions.checkState(!supersample);
		
		/*
		 * PRVI
		 */
////		region = PRVI25_RegionLoader.loadPRVI_ModelBroad();
//		region = PRVI25_RegionLoader.loadPRVI_MapExtents();
//		gridSpacing = 0.01;
//		
//		gmms = new AttenRelRef[] { AttenRelRef.USGS_PRVI_ACTIVE, AttenRelRef.USGS_PRVI_INTERFACE, AttenRelRef.USGS_PRVI_SLAB };
//		periods = new double[] { 0d, 0.2d, 1d, 5d };
//		supersample = true;
//		sigmaTrunc = 3d;
//		
////		String baseDirName = "2025_05_21-prvi25_crustal_subduction_combined_branches";
////		String suffix = "ba_only";
////		String solFileName = "combined_branch_averaged_solution.zip";
//		
//		String baseDirName = "2025_05_21-prvi25_crustal_subduction_combined_branches";
//		String suffix = "ba_only-slab_mc_7p4";
//		String solFileName = "combined_branch_averaged_solution.zip";
//		gridSpacing = 0.02;
//		
////		region = PRVI25_RegionLoader.loadPRVI_IntermediateModelMapExtents();
////		gridSpacing = 0.02;
////		periods = new double[] { 0d, 1d };
////		String baseDirName = "2025_05_21-prvi25_crustal_subduction_combined_branches";
////		String suffix = "ba_only-wider_region";
////		String solFileName = "combined_branch_averaged_solution.zip";
//		
////		String baseDirName = "2025_05_21-prvi25_crustal_branches-dmSample10x";
////		String suffix = "ba_only";
////		String solFileName = "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip";
//		
////		String baseDirName = "2025_05_21-prvi25_subduction_branches";
////		// slab (gridded only)
//////		String suffix = "ba_only-SLAB_only";
//////		String solFileName = "results_PRVI_SLAB_ONLY_branch_averaged_gridded.zip";
//////		bgOps = new IncludeBackgroundOption[] { IncludeBackgroundOption.ONLY };
////		// interface (will do fault + gridded)
//////		String suffix = "ba_only-INTERFACE_only";
//////		String solFileName = "results_PRVI_INTERFACE_ONLY_branch_averaged_gridded.zip";
////		// both
////		String suffix = "ba_only-both_fms";
////		String solFileName = "results_PRVI_SUB_FMs_combined_branch_averaged_gridded.zip";
//		
//		vs30 = 760d; suffix += "-vs760";
////		vs30 = 260d; suffix += "-vs260";
		
		/*
		 * RSQSim
		 */
//		String suffix = null;
//		
//		String baseDirName = "2024_11_12-rsqsim-wus-5895";
//		String solFileName = "fss_m6_skip10000_sectArea0.5_minSubSects2.zip";
//		
////		String baseDirName = "2024_11_12-rsqsim-wus-5892";
////		String solFileName = "fss_m6_skip20000_sectArea0.5_minSubSects2.zip";
//		linkFromBase = false;
//		
////		periods = new double[] {0d, 1d, 2d, 3d, 5d};
////		vs30 = 500d; suffix = "vs30_500";
//		periods = new double[] {2d, 3d, 5d};
//		gridSpacing = 0.5; suffix = "bbp_"+(float)gridSpacing;
		
		boolean noMFDs = false;
		
		GriddedRegion gridReg = new GriddedRegion(
				region, gridSpacing, GriddedRegion.ANCHOR_0_0);
//		GriddedRegion gridReg = new GriddedRegion(
//				PRVI25_RegionLoader.loadPRVI_ModelBroad(), gridSpacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		String dirName = forceOutputName == null ? baseDirName : forceOutputName;
		if (suffix != null && !suffix.isBlank())
			dirName += "-"+suffix;
		if (noMFDs)
			dirName += "-no_mfds";
		
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		if (parallelMPJWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)parallelMPJWrite).setUseLaunchWrapper(true);
		
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
			
			String argz = "--input-file "+dirPath+"/"+solFileName;
			argz += " --output-dir "+resultsPath;
			argz += " --output-file "+resultsPath+"_hazard_"+bgOp.name()+".zip";
			argz += " --region "+regPath;
			if (noMFDs)
				argz += " --no-mfds";
			
			argz += " --gridded-seis "+bgOp.name();
			if (bgOp != IncludeBackgroundOption.EXCLUDE) {
				if (distCorr != null)
					argz += " --dist-corr "+distCorr.name();
				if (bgRupType != null)
					argz += " --point-source-type "+bgRupType.name();
				if (bgFiniteNum != null)
					argz += " --point-finite-num-rand-surfaces "+bgFiniteNum;
				if (pointFiniteMinMag != null)
					argz += " --point-finite-min-mag "+pointFiniteMinMag.floatValue();
				if (extraGridArgs != null && !extraGridArgs.isBlank())
					argz += " "+extraGridArgs;
			}
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
			if (supersample) {
				if (supersampleFinite)
					argz += " --supersample-finite";
				else
					argz += " --supersample";
			}
			if (sigmaTrunc != null)
				argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
			argz += " "+dispatchArgs;
			
			File jobFile = new File(localDir, "batch_hazard_"+bgOp.name()+".slurm");
			
			List<String> script = mpjWrite.buildScript(MPJ_SingleSolHazardCalc.class.getName(), argz);
			
			System.out.println("Writing "+jobFile.getAbsolutePath());
			
			pbsWrite.writeScript(jobFile, script, mins, myNodes, remoteTotalThreads, queue);
		}
	}

}
