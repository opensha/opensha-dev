package scratch.kevin.pointSources;

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
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.kevin.pointSources.paperFigs2026.CalcPaths;

public class PointSourceBranchAveragedHazardScriptWriter {

	public static void main(String[] args) throws IOException {
		
		IncludeBackgroundOption[] bgOps = IncludeBackgroundOption.values();
		
		boolean linkFromBase = true;
		Double vs30 = null;
		Double sigmaTrunc = null;
		boolean supersample = false;
		boolean supersampleQuick = false;
		boolean supersampleFinite = false;
		boolean quickGridded = false;
		
		double[] periods = { 0d, 1d };
		AttenRelRef[] gmms = null;
		
		PointSourceDistanceCorrections distCorr = null;
		BackgroundRupType bgRupType = null;
		Integer bgFiniteNum = null;
		String extraGridArgs = null;
		Double pointFiniteMinMag = null;
		
		Integer forceMaxDispatch = null;
		
		String forceOutputName = null;
		String baseDirName = "2024_02_02-nshm23_branches-WUS_FM_v3";
		
		String suffix = "ba_only";
		forceOutputName = "2025_10_07-nshm23_pt_src_tests";
		String solFileName = "results_WUS_FM_v3_branch_averaged_gridded.zip";
		gmms = new AttenRelRef[] {AttenRelRef.USGS_NSHM23_ACTIVE};
//		gmms = new AttenRelRef[] {AttenRelRef.NGAWest_2014_AVG_NOIDRISS};
		
//		GriddedRegion gridReg = CalcPaths.FULL_GRID_REG;
		GriddedRegion gridReg = CalcPaths.ZOOM_GRID_REG; forceOutputName += "-zoom";
//		GriddedRegion gridReg = new GriddedRegion(new Region(new Location(37, -119), new Location(37.6, -118.4)),
//				0.01, GriddedRegion.ANCHOR_0_0); forceOutputName += "-tiny_zoom";
		sigmaTrunc = 3d;
		supersample = true;
		supersampleFinite = true; // only applies if supersample == true
		
//		forceOutputName += "-as_published";
//		pointFiniteMinMag = 6d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.NSHM_2013;
		
//		forceOutputName += "-finite_20x";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 20;
		
//		forceOutputName += "-finite_1x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 1;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
		
//		forceOutputName += "-finite_1x_along_alt_rand";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 1;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip "
//				+ "-D"+PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME+"="+System.currentTimeMillis();
		
//		forceOutputName += "-finite_2x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 2;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
		
//		forceOutputName += "-finite_2x_along_alt_rand";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 2;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip "
//				+ "-D"+PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME+"="+System.currentTimeMillis();
		
//		forceOutputName += "-finite_20x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 20;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
		
//		forceOutputName += "-finite_40x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 40;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip";
		
//		forceOutputName += "-finite_50x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 50;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip "
//				+ "-D"+PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME+"="+System.currentTimeMillis();
		
//		forceOutputName += "-finite_100x_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.FINITE;
//		bgFiniteNum = 100;
//		extraGridArgs = "--point-finite-sample-along-strike --point-finite-sample-down-dip "
//				+ "-D"+PointSurfaceBuilder.GLOBAL_SEED_PROP_NAME+"="+System.currentTimeMillis();
		
//		forceOutputName += "-average_spinning";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING;
		
//		forceOutputName += "-average_spinning_m6";
//		pointFiniteMinMag = 6d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING;
		
//		forceOutputName += "-average_spinning_along";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.AVERAGE_SPINNING_ALONG;
		
//		forceOutputName += "-five_pt_spinning";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST;
		
//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-evenly_spaced_dist";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
		
//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-m4";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		extraGridArgs = "--point-min-mag 4";
		
//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-m3.5";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		extraGridArgs = "--point-min-mag 3.5";
		
//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-m3";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		extraGridArgs = "--point-min-mag 3";
		
		forceOutputName += "-five_pt_spinning_along";
		forceOutputName += "-alt_grid_depths";
		pointFiniteMinMag = 5d;
		bgRupType = BackgroundRupType.POINT;
		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
		linkFromBase = false;
		solFileName = "results_WUS_FM_v3_branch_averaged_gridded_mod_grid_depths.zip";
		
//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-alt_wc_lengths";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		linkFromBase = false;
//		solFileName = "results_WUS_FM_v3_branch_averaged_gridded_mod_wc_lengths.zip";

//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-alt_leonard_lengths";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		linkFromBase = false;
//		solFileName = "results_WUS_FM_v3_branch_averaged_gridded_mod_leonard_lengths.zip";

//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-proposed_grid_mods";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		linkFromBase = false;
//		solFileName = "results_WUS_FM_v3_branch_averaged_gridded_mod_proposed.zip";

//		forceOutputName += "-five_pt_spinning_along";
//		forceOutputName += "-proposed_grid_mods-m3.5";
//		pointFiniteMinMag = 5d;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
//		linkFromBase = false;
//		solFileName = "results_WUS_FM_v3_branch_averaged_gridded_mod_proposed.zip";
//		extraGridArgs = "--point-min-mag 3.5";
		
//		forceOutputName += "-five_pt_spinning_along";
//		pointFiniteMinMag = 5d;
//		supersample = false;
//		bgRupType = BackgroundRupType.POINT;
//		distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG;
		
//		pointFiniteMinMag = 5d;
		forceMaxDispatch = 100;
		if (supersample)
			forceOutputName += "-supersample";
		
		boolean noMFDs = false;
		
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		String dirName = forceOutputName == null ? baseDirName : forceOutputName;
		if (suffix != null && !suffix.isBlank())
			dirName += "-"+suffix;
		if (noMFDs)
			dirName += "-no_mfds";
		
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
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
		parallelMPJWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		singleMPJWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		parallelMPJWrite.setClasspath(classpath);
		
		List<File> singleClasspath = new ArrayList<>(classpath);
		singleClasspath.add(new File("/project2/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
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
			if (supersampleQuick || supersample) {
				if (supersampleQuick)
					argz += " --supersample-quick";
				if (supersampleFinite)
					argz += " --supersample-finite";
				else if (!supersampleQuick)
					argz += " --supersample";
			}
			if (sigmaTrunc != null)
				argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
			if (quickGridded && bgOp != IncludeBackgroundOption.EXCLUDE)
				argz += " --quick-grid-calc";
			argz += " "+dispatchArgs;
			
			File jobFile = new File(localDir, "batch_hazard_"+bgOp.name()+".slurm");
			
			List<String> script = mpjWrite.buildScript(MPJ_SingleSolHazardCalc.class.getName(), argz);
			
			System.out.println("Writing "+jobFile.getAbsolutePath());
			
			pbsWrite.writeScript(jobFile, script, mins, myNodes, remoteTotalThreads, queue);
		}
	}

}
