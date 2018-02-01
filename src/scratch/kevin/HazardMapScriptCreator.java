package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.mpj.MPJHazardCurveDriver;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class HazardMapScriptCreator {

	public static void main(String[] args) throws IOException {
		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(PGA_Param.NAME);
//		gmpe.setIntensityMeasure(SA_Param.NAME);
//		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 1d);
		
		File fssFile = new File("/auto/scec-02/kmilner/ucerf3/erf_cache/cached_dep100.0_depMean_rakeMean.zip");
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF();
		erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, fssFile.getAbsolutePath());
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(1d);
		
//		Region region = new CaliforniaRegions.LA_BOX();
		Region region = new Region(new Location(34.5, -117.2), new Location(33.6, -119));
		double spacing = 0.005;
		
		List<SiteData<?>> dataProvs = new ArrayList<>();
//		dataProvs.add(new WillsMap2015());
//		dataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
//		dataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		
		double maxSourceDist = 200d;
		
		int mins = 60*10;
		int nodes = 32;
		int ppn = 20;
		String queue = "scec";
		int maxMemMB = 60000;
		
		DateFormat df = new SimpleDateFormat("yyyy_MM_dd");
//		String jobName = df.format(new Date())+"-ucerf3-ba-la-zoomed-site-effects-pga";
		String jobName = df.format(new Date())+"-ucerf3-ba-la-zoomed-pga-rock";
		System.out.println(jobName);
		
		File localDir = new File("/home/kevin/OpenSHA/maps");
		File remoteDir = new File("/auto/scec-02/kmilner/hazMaps");
		
		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		File jarFile = new File(remoteDir, "git/opensha-dev/build/libs/opensha-dev-all.jar");
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(jarFile);
		
		MPJExpressShellScriptWriter mpj = new MPJExpressShellScriptWriter(javaBin, maxMemMB, classpath,
				USC_HPCC_ScriptWriter.MPJ_HOME);
		
		File localJobDir = new File(localDir, jobName);
		Preconditions.checkState(localJobDir.exists() || localJobDir.mkdir());
		File remoteJobDir = new File(remoteDir, jobName);
		
		List<Site> sites = new ArrayList<>();
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		List<SiteDataValueList<?>> allVals = new ArrayList<>();
		for (SiteData<?> prov : dataProvs) {
			System.out.println("Fetching site data for: "+prov.getName());
			allVals.add(prov.getAnnotatedValues(gridReg.getNodeList()));
		}
		SiteTranslator trans = new SiteTranslator();
		System.out.println("Creating "+gridReg.getNodeCount()+" sites/setting site params");
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			List<SiteDataValue<?>> siteDatas = new ArrayList<>();
			for (SiteDataValueList<?> provVals : allVals)
				siteDatas.add(provVals.getValue(i));
			Site site = new Site(gridReg.getLocation(i));
			for (Parameter<?> param : gmpe.getSiteParams()) {
				param = (Parameter<?>) param.clone();
				trans.setParameterValue(param, siteDatas);
				site.addParameter(param);
			}
			sites.add(site);
		}
		
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps = new ArrayList<>();
		Map<TectonicRegionType, ScalarIMR> imrMap = new HashMap<>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, gmpe);
		imrMaps.add(imrMap);
		
		DiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		CalculationSettings calcSettings = new CalculationSettings(xValues, maxSourceDist);
		
		Map<String, DiscretizedFunc> xValsMap = new HashMap<>();
		xValsMap.put("curves", xValues);
		BinaryCurveArchiver archiver = new BinaryCurveArchiver(new File(remoteJobDir, "curves"), sites.size(), xValsMap);
		
		CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
		
		File localInputsFile = new File(localJobDir, "inputs.xml");
		File remoteInputsFile = new File(remoteJobDir, "inputs.xml");
		XMLUtils.writeObjectToXMLAsRoot(inputs, localInputsFile);
		
		int maxDispatch = 1000;
		if (sites.size() < 50000)
			maxDispatch = 500;
		if (sites.size() < 10000)
			maxDispatch = 100;
		
//		String cliArgs = "--max-dispatch 1000 --mult-erfs "+inputsFile.getAbsolutePath();
		String cliArgs = "--max-dispatch "+maxDispatch+" "+remoteInputsFile.getAbsolutePath();
		
		List<String> script = mpj.buildScript(MPJHazardCurveDriver.class.getName(), cliArgs);
		USC_HPCC_ScriptWriter writer = new USC_HPCC_ScriptWriter();
		
		script = writer.buildScript(script, mins, nodes, ppn, queue);
		
		File pbsFile = new File(localJobDir, jobName+".pbs");
		System.out.println("Writing "+pbsFile.getAbsolutePath());
		JavaShellScriptWriter.writeScript(pbsFile, script);
	}

}
