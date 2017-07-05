package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.kevin.ucerf3.MPJ_UCERF3_ShakeMapPrecalcScriptGen;

public class MPJ_ETAS_HazardMapCalcScriptGen {
	
	public static final DateFormat df = new SimpleDateFormat("yyyy_MM_dd");

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard");
		
		String longTermDurations = null;
//		String longTermDurations = "THREE";
		
		// 5/17 USGS Exercise
//		String etasSimName = "2017_05_17-USGS_Exercise_Catalog-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont";
		String etasSimName = "2017_05_17-USGS_Exercise_Catalog-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont-sect-reset-1pm";
		String etasFileName = "results_m5.bin";
		String etasShortName = "2017_05-usgs_exercise-1pm";
		TestScenario scenario = null;
		// Haywired Fault
//		String etasSimName = "2016_06_15-haywired_m7-10yr-full_td-no_ert-combined";
//		String etasFileName = "results_descendents_m5.bin";
//		String etasShortName = "haywired_m7_combined_descendents";
//		TestScenario scenario = TestScenario.HAYWIRED_M7;
		// Haywired Gridded
//		String etasSimName = "2017_01_02-haywired_m7-10yr-gridded-only-200kcombined";
//		String etasFileName = "results_descendents_m5_preserve.bin";
//		String etasShortName = "haywired_m7_gridded_descendents";
//		TestScenario scenario = TestScenario.HAYWIRED_M7;
		// Northridge Fault
//		String etasSimName = "2017_02_01-northridge-m6.7-10yr-full_td-no_ert-combined";
//		String etasFileName = "results_descendents_m5.bin";
//		String etasShortName = "northridge_combined_descendents";
//		TestScenario scenario = TestScenario.NORTHRIDGE;
		// Northride Gridded
//		String etasSimName = "2017_02_01-northridge-m6.7-10yr-gridded-only-combined200k";
//		String etasFileName = "results_descendents_m5_preserve.bin";
//		String etasShortName = "northridge_gridded_descendents";
//		TestScenario scenario = TestScenario.NORTHRIDGE;
		// Mojave Fault
//		String etasSimName = "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined";
//		String etasFileName = "results_descendents_m5_preserve.bin";
//		String etasShortName = "mojave_m7_combined_descendents";
//		TestScenario scenario = TestScenario.MOJAVE_M7;
		// Mojave Gridded
//		String etasSimName = "2016_12_03-mojave_m7-10yr-gridded-only";
//		String etasFileName = "results_descendents_m5_preserve.bin";
//		String etasShortName = "mojave_m7_gridded_descendents";
//		TestScenario scenario = TestScenario.MOJAVE_M7;
		// 2016 Bombay Swarm Fault
//		String etasSimName = "2016_10_27-2016_bombay_swarm-10yr-full_td-no_ert-combined";
//		String etasFileName = "results_m5_preserve.bin";
//		String etasShortName = "2016_bombay_swarm_combined_descendents";
//		TestScenario scenario = null;
		
		// --------------------
		// for fault precalc
//		String shakemapRunName = "2017_02_23-NGAWest_2014_NoIdr-spacing0.05-site-effects-with-basin";
//		String gmpeFileName = "NGAWest_2014_NoIdr.xml";
//		boolean siteEffects = true;
//		boolean basinDepth = true;
//		AttenuationRelationship gmpe = null;
		// for faults on the fly
		String shakemapRunName = null;
		AttenRelRef gmpeRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		AttenuationRelationship gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		boolean siteEffects = true;
		boolean basinDepth = true;
		// --------------------
		String[] imts = { "pgv", "pga" };
		double[] periods = { Double.NaN, Double.NaN };
		double spacing = 0.02;
		String shakemapShortName = "NGA2-"+(float)spacing;
		if (siteEffects) {
			shakemapShortName += "-site-effects";
			if (basinDepth)
				shakemapShortName += "-with-basin";
			else
				shakemapShortName += "-no-basin";
		} else {
			shakemapShortName += "-no-site-effects";
		}
		
		double griddedSpacing = 0.01;
		
		String dateStr = df.format(new Date());
//		String dateStr = "2017_03_23";
		String jobName = dateStr+"-"+etasShortName+"-"+shakemapShortName;
		
//		boolean stampede = true;
//		boolean knl = false;
//		int nodes = 34;
//		int hours = 24;
//		int threads = 16;
////		boolean knl = true;
////		int nodes = 10;
////		int hours = 10;
////		int threads = 272;
//		String queue = null;
//		boolean scecLarge = false;
		
		boolean stampede = false;
		boolean knl = false;
		boolean scecLarge = false;
		int nodes;
		if (scecLarge)
			nodes = 4;
		else
			nodes = 34;
		int hours = 24;
		int threads = 20;
		String queue = "scec";
		
		int minDispatch = threads;
		if (minDispatch < MPJTaskCalculator.MIN_DISPATCH_DEFAULT)
			minDispatch = MPJTaskCalculator.MIN_DISPATCH_DEFAULT;
		int maxDispatch = MPJTaskCalculator.MAX_DISPATCH_DEFAULT;
		if (spacing <= 0.05)
			maxDispatch = 500;
		
		File localDir = new File(localMainDir, jobName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		System.out.println("Job name: "+jobName);
		System.out.println("Dir: "+localDir.getAbsolutePath());
		
		JavaShellScriptWriter mpjWrite;
		BatchScriptWriter pbsWrite;
		
		int memGigs;
		int ppn;
		File remoteMainDir, remoteShakemapDir, remoteETASDir, remoteFSSFile;
		if (stampede) {
			if (knl) {
				memGigs = 90;
				ppn = 68;
			} else {
				memGigs = 26;
				ppn = 16;
			}
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
					null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			pbsWrite = new StampedeScriptWriter(knl);
			
			remoteMainDir = new File("/work/00950/kevinm/ucerf3/etas_hazard");
			
			remoteShakemapDir = null;
			
			remoteETASDir = new File("/work/00950/kevinm/ucerf3/etas_sim");
			remoteFSSFile = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		} else {
			if (queue == null) {
				memGigs = 9;
				ppn = 8;
			} else {
				if (scecLarge)
					memGigs = 200;
				else
					memGigs = 60;
				ppn = 20;
			}
			boolean fmpj = nodes < 25;
			fmpj = false;
			if (fmpj) {
				mpjWrite = new FastMPJShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.FMPJ_HOME);
				((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			} else {
				mpjWrite = new MPJExpressShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.MPJ_HOME);
			}
			pbsWrite = new USC_HPCC_ScriptWriter();
			
			remoteMainDir = new File("/home/scec-02/kmilner/ucerf3/etas_hazard");
			
			remoteShakemapDir = new File("/home/scec-02/kmilner/ucerf3/shakemap_precalc");
			
//			remoteETASDir = new File("/home/scec-00/kmilner/ucerf3_etas_results_stampede/");
			remoteETASDir = new File("/home/scec-02/kmilner/ucerf3/etas_sim/");
			remoteFSSFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		}
		File remoteEtasCatalogFile = new File(new File(remoteETASDir, etasSimName), etasFileName);
		File remoteJobDir = new File(remoteMainDir, jobName);
		
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteMainDir, "commons-cli-1.2.jar"));
		classpath.add(new File(remoteJobDir, "OpenSHA_complete.jar"));
		mpjWrite.setClasspath(classpath);
		
		File sitesFile;
		File remoteShakemapRunDir = null;
		if (shakemapRunName == null) {
			// calculating on the fly
			ArrayList<SiteData<?>> siteData = null;
			if (siteEffects)
				siteData = MPJ_UCERF3_ShakeMapPrecalcScriptGen.getSiteDataProviders(basinDepth);
			File localSitesFile = new File(localDir, "sites.xml");
			GriddedRegion reg = new CaliforniaRegions.RELM_TESTING_GRIDDED(spacing);
			MPJ_UCERF3_ShakeMapPrecalcScriptGen.writeSitesXML(reg, siteData, gmpe, localSitesFile);
			sitesFile = new File(remoteJobDir, localSitesFile.getName());
		} else {
			remoteShakemapRunDir = new File(remoteShakemapDir, shakemapRunName);
			sitesFile = new File(remoteShakemapRunDir, "sites.xml");
		}
		
		for (int i=0; i<imts.length; i++) {
			String imt = imts[i];
			double period = periods[i];
			String imtName = imt;
			if (!Double.isNaN(period) && period > 0)
				imtName += "_"+(float)period+"s";
			
			String argz = "--min-dispatch "+minDispatch+" --max-dispatch "+maxDispatch;
			if (threads > 0)
				argz += " --threads "+threads;
			argz += " --catalogs "+remoteEtasCatalogFile.getAbsolutePath();
			if (remoteShakemapRunDir == null) {
				// calc on the fly
				argz += " --solution-file "+remoteFSSFile.getAbsolutePath();
			} else {
				// use precalc
				File remoteShakemapFile = new File(remoteShakemapRunDir, "results_"+imtName+".bin");
				argz += " --fault-data-file "+remoteShakemapFile.getAbsolutePath();
			}
			argz += " --spacing "+(float)spacing;
			argz += " --gmpe "+gmpeRef.name();
			argz += " --sites-file "+sitesFile.getAbsolutePath();
			argz += " --imt "+imt;
			if (!Double.isNaN(period) && period > 0)
				argz += " --period"+(float)period;
			argz += " --gridded-spacing "+(float)griddedSpacing;
			argz += " --output-dir "+remoteJobDir.getAbsolutePath();
			if (etasSimName.contains("gridded-only"))
				argz += " --no-fault";
			if (longTermDurations != null && !longTermDurations.isEmpty()) {
				argz += " --calc-long-term --long-term-durations "+longTermDurations;
				if (scenario != null)
					argz += " --elastic-rebound "+scenario.name();
			}
			
			List<String> script = mpjWrite.buildScript(MPJ_ETAS_HazardMapCalc.class.getName(), argz);
			
			int mins = hours*60;
			script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
			String scriptName;
			if (stampede) {
				if (knl)
					scriptName = imt+"_stampede_knl.pbs";
				else
					scriptName = imt+"_stampede.pbs";
			} else {
				scriptName = imt+".pbs";
			}
			pbsWrite.writeScript(new File(localDir, scriptName), script);
		}
	}

}
