package scratch.kevin.cybershake;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVMHBasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_CCAi6BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.AsciiFileCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.components.CurveResultsArchiver;
import org.opensha.sha.calc.hazardMap.mpj.MPJHazardCurveDriver;
import org.opensha.sha.calc.hazus.parallel.HazusJobWriter;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.gui.controls.CyberShakePlotFromDBControlPanel;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class CyberShakeBaseMapGen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 9 || args.length > 10) {
			System.out.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(CyberShakeBaseMapGen.class)
					+" <IMRs> <SA period> <spacing> <CVM4/CVMH/CVMHnGTL/BBP/CVM4i26/CCAi6/CCA1D/null> <constrainBasinMin>"
					+" <jobName> <minutes> <nodes> <queue> [<LA/CCA/CA>]");
			System.exit(2);
		}
		
		// TODO args
		String imrNames = args[0];
		double period = Double.parseDouble(args[1]);
		double spacing = Double.parseDouble(args[2]);
		String cvmName = args[3];
		boolean constrainBasinMin = Boolean.parseBoolean(args[4]);
		String jobName = args[5];
		int mins = Integer.parseInt(args[6]);
		int nodes = Integer.parseInt(args[7]);
		int ppn = 8;
		String queue = args[8];
		if (queue.equals("scec"))
			ppn = 20;
		else if (queue.toLowerCase().equals("null"))
			queue = null;
		GriddedRegion region = new CaliforniaRegions.CYBERSHAKE_MAP_GRIDDED(spacing);
		if (args.length > 9) {
			String regName = args[9].toLowerCase();
			if (regName.equals("la"))
				region = new CaliforniaRegions.CYBERSHAKE_MAP_GRIDDED(spacing);
			else if (regName.equals("ccs"))
				region = new CaliforniaRegions.CYBERSHAKE_CCA_MAP_GRIDDED(spacing);
			else if (regName.equals("ca"))
				region = new CaliforniaRegions.RELM_TESTING_GRIDDED(spacing);
			else
				throw new IllegalArgumentException("Unknown region: "+args[9]);
		}
		
		File hazMapsDir = new File("/home/scec-02/kmilner/hazMaps");
		
		File jobDir = new File(hazMapsDir, jobName);
		if (!jobDir.exists())
			jobDir.mkdir();
		
		ArrayList<ScalarIMR> imrs = Lists.newArrayList();
		for (String imrName : Splitter.on(",").split(imrNames)) {
			ScalarIMR imr = AttenRelRef.valueOf(imrName).instance(null);
			imr.setParamDefaults();
			imr.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
			imr.getParameter(SigmaTruncLevelParam.NAME).setValue(3d);
			if (period > 0) {
				imr.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
			} else {
				imr.setIntensityMeasure(PGA_Param.NAME);
			}
			imrs.add(imr);
		}
		
		List<SiteData<?>> provs = Lists.newArrayList();
//		provs.add(new WillsMap2006());
		provs.add(new WillsMap2015());
		boolean nullBasin = false;
		if (cvmName.toLowerCase().equals("cvm4")) {
			provs.add(new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			provs.add(new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		} else if (cvmName.toLowerCase().equals("cvm4i26")) {
			provs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			provs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		} else if (cvmName.toLowerCase().equals("ccai6")) {
			provs.add(new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			provs.add(new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		} else if (cvmName.toLowerCase().equals("cca1d")) {
			provs = HazardCurvePlotter.getCCA_1D_Providers();
		} else if (cvmName.toLowerCase().equals("cvmh")) {
			provs.add(new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			provs.add(new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		} else if (cvmName.toLowerCase().equals("cvmhngtl")) {
			CVMHBasinDepth z25 = new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
			CVMHBasinDepth z10 = new CVMHBasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
			z25.getAdjustableParameterList().setValue(CVMHBasinDepth.GTL_PARAM_NAME, false);
			z10.getAdjustableParameterList().setValue(CVMHBasinDepth.GTL_PARAM_NAME, false);
			provs.add(z25);
			provs.add(z10);
		} else if (cvmName.toLowerCase().equals("bbp")) {
			// this will also clear Wills 2006 from the list
			provs = HazardCurvePlotter.getBBP_1D_Providers();
		} else if (cvmName.toLowerCase().equals("null")){
			nullBasin = true;
		} else {
			System.err.println("Unknown basin model: "+cvmName);
		}
		
		ArrayList<Site> sites = HazusJobWriter.loadSites(provs, region.getNodeList(), imrs, nullBasin, constrainBasinMin, null);
		
		ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		
		ArbitrarilyDiscretizedFunc xValues = CyberShakePlotFromDBControlPanel.createUSGS_PGA_Function();
		double maxSourceDistance = 200;
		
		Map<String, DiscretizedFunc> xValsMap = Maps.newHashMap();
		xValsMap.put("imrs1", xValues);
		CalculationSettings calcSettings = new CalculationSettings(xValues, maxSourceDistance);
		
		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		File svnDir = new File(hazMapsDir, "svn");
		File distDir = new File(svnDir, "dist");
		File libDir = new File(svnDir, "lib");
		File jarFile = new File(distDir, "OpenSHA_complete.jar");
		
		ArrayList<File> classpath = new ArrayList<File>();
		classpath.add(jarFile);
		classpath.add(new File(libDir, "commons-cli-1.2.jar"));
		
		MPJExpressShellScriptWriter mpj = new MPJExpressShellScriptWriter(javaBin, 60000, classpath,
				USC_HPCC_ScriptWriter.MPJ_HOME);
		
		for (ScalarIMR imr : imrs) {
			List<Map<TectonicRegionType, ScalarIMR>> imrMaps = Lists.newArrayList();
			
			HashMap<TectonicRegionType, ScalarIMR> map = Maps.newHashMap();
			map.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
			imrMaps.add(map);
			
			File imrDir = new File(jobDir, imr.getShortName());
			if (!imrDir.exists())
				imrDir.mkdir();
			
			File curveDir = new File(imrDir, "curves");
//			CurveResultsArchiver archiver = new AsciiFileCurveArchiver(
//				curveDir.getAbsolutePath()+File.separator, true, false);
			CurveResultsArchiver archiver = new BinaryCurveArchiver(curveDir, sites.size(), xValsMap);
			
			CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
			
			File inputsFile = new File(imrDir, "inputs.xml");
			XMLUtils.writeObjectToXMLAsRoot(inputs, inputsFile);
			
			String cliArgs = "--max-dispatch 1000 --mult-erfs "+inputsFile.getAbsolutePath();
			
			List<String> script = mpj.buildScript(MPJHazardCurveDriver.class.getName(), cliArgs);
			USC_HPCC_ScriptWriter writer = new USC_HPCC_ScriptWriter();
			
			script = writer.buildScript(script, mins, nodes, ppn, queue);
			
			File pbsFile = new File(imrDir, imr.getShortName().toLowerCase()+".pbs");
			JavaShellScriptWriter.writeScript(pbsFile, script);
		}
	}

}
