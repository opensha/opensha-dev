package scratch.kevin.nshm23.wrapper;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_WrapperHazardCalcScriptWriter {

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		int nodes = 36;
		IncludeBackgroundOption griddedOp = IncludeBackgroundOption.ONLY;
		boolean subduction = false;
		
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;
//		AttenRelRef gmpeRef = AttenRelRef.NGAWest_2014_AVG;
		
//		String erfPrefix = "nshm18";
//		String tagName = "nshm-conus-5.3.0"; // NSHM18
		
		String erfPrefix = "nshm23-wrapped";
		String tagName = "nshm-conus-6.a.6"; // NSHM23 draft
		
		double gridSpacing = 0.2d; int mins = 2000;
//		double gridSpacing = 0.1d; int mins = 2000;
		
//		String regName = "wus";
//		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		String regName = "conus";
		Region region = NSHM23_RegionLoader.loadFullConterminousUS();
		
		String extGridProvPath = null;
//		String extGridProvPath = "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
//		erfPrefix += "-grid_src_from_23"; griddedOp = IncludeBackgroundOption.INCLUDE;
//		String extGridProvPath = "2022_12_07-nshm23_branches-no_paleo_slip-mod_dm_weights-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
//		erfPrefix += "-grid_src_from_modWeightDM_23"; gridded = true;
//		String extGridProvPath = "2022_12_07-nshm23_branches-no_paleo_slip-NSHM23_v2-CoulombRupSet-MEDIAN-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
//				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
//		erfPrefix += "-grid_src_from_medDM_23"; gridded = true;
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
//		JavaShellScriptWriter mpjWrite = new NoMPJSingleNodeShellScriptWriter(USC_CARC_ScriptWriter.JAVA_BIN,
//				remoteTotalMemGB*1024, null); nodes = 1; remoteInversionsPerBundle = 2;
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
		
		dirName += "-"+erfPrefix+"-"+regName+"-hazard-"+gmpeRef.getShortName().toLowerCase()+"-"+(float)gridSpacing+"deg";
		if (!subduction)
			dirName += "-noSub";
		if (griddedOp == IncludeBackgroundOption.EXCLUDE)
			dirName += "-faultOnly";
		else if (griddedOp == IncludeBackgroundOption.ONLY)
			dirName += "-griddedOnly";
		
		System.out.println(dirName);
		
		GriddedRegion gridReg = new GriddedRegion(region, gridSpacing, GriddedRegion.ANCHOR_0_0);
		
		File localDir = new File(localMainDir, dirName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		mpjWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		mpjWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		mpjWrite.setEnvVar("MODEL_PATH", mainDirPath+"/nshmp-haz-models/"+tagName);
		
		String remoteERFsPath = "$MODEL_PATH";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		if (mpjWrite instanceof NoMPJSingleNodeShellScriptWriter)
			classpath.add(new File("/project/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		
		File localReg = new File(localDir, "gridded_region.json");
		Feature.write(gridReg.toFeature(), localReg);
		
		mpjWrite.setClasspath(classpath);
		if (mpjWrite instanceof MPJExpressShellScriptWriter)
			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		else if (mpjWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		String resultsPath = dirPath+"/results";
		String regPath = dirPath+"/"+localReg.getName();
		
		String argz = "--input-dir "+remoteERFsPath+" --region "+regPath;
		argz += " --output-dir "+resultsPath;
		if (griddedOp == IncludeBackgroundOption.INCLUDE || griddedOp == IncludeBackgroundOption.ONLY
				|| extGridProvPath != null) {
//			argz += " --max-distance 200";
			if (extGridProvPath != null)
				argz += " --external-grid-prov $MAIN_DIR/"+extGridProvPath;
		}
		argz += " --gridded-seis "+griddedOp.name();
		if (!subduction)
			argz += " --no-subduction";
		argz += " --gmpe "+gmpeRef.name();
		argz += " "+MPJTaskCalculator.argumentBuilder().minDispatch(remoteTotalThreads).build();
		List<String> script = mpjWrite.buildScript(MPJ_WrapperHazardCalc.class.getName(), argz);
		
		// make sure to not exceed 1 week
		mins = Integer.min(mins, 60*24*7 - 1);
		pbsWrite.writeScript(new File(localDir, "batch_wrapper_calc.slurm"), script, mins, nodes, remoteTotalThreads, queue);
	}

}
