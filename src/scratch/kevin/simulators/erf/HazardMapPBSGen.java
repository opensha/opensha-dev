package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVMHBasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.AsciiFileCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.components.CurveResultsArchiver;
import org.opensha.sha.calc.hazardMap.mpj.MPJHazardCurveDriver;
import org.opensha.sha.calc.hazus.parallel.HazusJobWriter;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.util.TectonicRegionType;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter.RunSites;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class HazardMapPBSGen {
	
	public static DateFormat df = new SimpleDateFormat("yyyy_MM_dd");

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		RunSites runSite = RunSites.HPCC;
		BatchScriptWriter writer = new USC_HPCC_ScriptWriter();
		int nodes = 30;
		int mins = 60 * 5;
		int heapMB = 9000;
//		String runName = "rsqsim-long-la-box-0.05-1sec-cb2008-30yrs";
		String runName = "rsqsim-long-la-box-0.05-1sec-cb2008-30yrs-quietmojave156";
//		String runName = "ucerf2-compare-time-1sec-indep";
		String pbsName = runName;
		runName = df.format(new Date())+"-"+runName;
		
		int durationYears = 30;
		
//		String solFileName = "simulators_long_sol.zip";
		String solFileName = "simulators_long_sol_mojave_trigger_quiet_156_wind_30_yr.zip";
		
		boolean ucerf2_compare = false;
		String ucerf_prob_model = UCERF2.PROB_MODEL_POISSON;
//		String ucerf_prob_model = MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND;
		
		File localDir = new File("/home/kevin/Simulators/maps");
		File remoteDir = new File("/home/scec-02/kmilner/hazMaps");
		
		File localRunDir = new File(localDir, runName);
		if (!localRunDir.exists())
			localRunDir.mkdir();
		File remoteRunDir = new File(remoteDir, runName);
		
		GriddedRegion region = new CaliforniaRegions.LA_BOX_GRIDDED(0.01);
		
		System.out.println("Creating input file for "+region.getNodeCount()+" sites");
		
		ScalarIMR imr = new CB_2008_AttenRel(null);
		
		imr.setParamDefaults();
		imr.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
		imr.getParameter(SigmaTruncLevelParam.NAME).setValue(2d);
//		imr.setIntensityMeasure(PGA_Param.NAME);
		imr.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), 1d);
		List<ScalarIMR> imrs = Lists.newArrayList(imr);
		
		ArrayList<SiteData<?>> provs = Lists.newArrayList();
		provs.add(new WillsMap2006());
		provs.add(new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		provs.add(new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		List<Site> sites = HazusJobWriter.loadSites(provs, region.getNodeList(), imrs, false, false, null);
		
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps = Lists.newArrayList();
		Map<TectonicRegionType, ScalarIMR> imrMap = Maps.newHashMap();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
		imrMaps.add(imrMap);
		
		IMT_Info imtInfo = new IMT_Info();
		ArbitrarilyDiscretizedFunc xValues = imtInfo.getDefaultHazardCurve(imr.getIntensityMeasure());
		
		double maxSourceDistance = 200;
		
		CalculationSettings calcSettings = new CalculationSettings(xValues, maxSourceDistance);
		File curveDir = new File(remoteRunDir, "curves");
		CurveResultsArchiver archiver = new AsciiFileCurveArchiver(curveDir.getAbsolutePath(), true, false);
		
		File solFile = new File(remoteRunDir, solFileName);
		
		AbstractERF erf;
		if (ucerf2_compare) {
			erf = new MeanUCERF2();
			erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
			erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, ucerf_prob_model);
			erf.getTimeSpan().setDuration(durationYears, TimeSpan.YEARS);
			
			if (ucerf_prob_model.equals(MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND))
				erf.getTimeSpan().setStartTime(2013);
			
		} else {
			erf = new FaultSystemSolutionERF();
			erf.getParameter(FaultSystemSolutionERF.FILE_PARAM_NAME).setValue(solFile);
			erf.getParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME).setValue(false);
			erf.getTimeSpan().setDuration(durationYears, TimeSpan.YEARS);
		}
		
		CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
		
		File localInputsFile = new File(localRunDir, "inputs.xml");
		File remoteInputsFile = new File(remoteRunDir, "inputs.xml");
		XMLUtils.writeObjectToXMLAsRoot(inputs, localInputsFile);
		
		String cliArgs = remoteInputsFile.getAbsolutePath();
		
		FastMPJShellScriptWriter mpj = new FastMPJShellScriptWriter(runSite.getJAVA_BIN(), heapMB,
				LogicTreePBSWriter.getClasspath(remoteRunDir, remoteRunDir), runSite.getMPJ_HOME());
		
		List<String> script = mpj.buildScript(MPJHazardCurveDriver.class.getName(), cliArgs);
		
		script = writer.buildScript(script, mins, nodes, 0, null);
		
		File pbsFile = new File(localRunDir, pbsName+".pbs");
		JavaShellScriptWriter.writeScript(pbsFile, script);
	}

}
