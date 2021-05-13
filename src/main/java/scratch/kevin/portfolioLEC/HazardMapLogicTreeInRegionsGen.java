package scratch.kevin.portfolioLEC;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.RegionUtils;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.HazardCurveDriver;
import org.opensha.sha.calc.hazardMap.components.AsciiFileCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.components.CurveResultsArchiver;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2_TimeDependentEpistemicList;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.util.SiteTranslator;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class HazardMapLogicTreeInRegionsGen {
	
	static List<GriddedRegion> getRegions() throws MalformedURLException, DocumentException {
		return getRegions(new File(
				"/home/kevin/OpenSHA/portfolio_lec/SPA Risk LLC (10 Jul 2012) Geovera Add1 study areas.kml"));
	}
	
	static List<GriddedRegion> getRegions(File regionFile) throws MalformedURLException, DocumentException {
		List<GriddedRegion> regions = Lists.newArrayList();
		
		Document doc = XMLUtils.loadDocument(regionFile);
		
		Element root = doc.getRootElement().element("Document");
		Element folderEl = root.element("Folder");
		Iterator<Element> placemarkIt = folderEl.elementIterator("Placemark");
		while (placemarkIt.hasNext()) {
			Element placeEl = placemarkIt.next();
			String name = placeEl.elementText("name");
			double gridSpacing;
			if (name.contains("1-km"))
				gridSpacing = 0.01;
			else
				gridSpacing = 0.05;
			Element coordinatesEl = placeEl.element("Polygon").element("outerBoundaryIs")
					.element("LinearRing").element("coordinates");
			String coordsStr = coordinatesEl.getTextTrim();
			coordsStr = coordsStr.replaceAll("\n", " ");
			while (coordsStr.contains("  "))
				coordsStr.replaceAll("  ", " ");
			if (coordsStr.isEmpty()) {
				System.out.println("WARNING: empty polygon: "+name);
				continue;
			}
			String[] locStrs = coordsStr.split(" ");
			LocationList locs = new LocationList();
			for (String locStr : locStrs) {
				String[] strs = locStr.split(",");
				double lon = Double.parseDouble(strs[0]);
				double lat = Double.parseDouble(strs[1]);
				Double.parseDouble(strs[2]); // depth not needed but parsed as validation
				locs.add(new Location(lat, lon));
			}
			GriddedRegion region = new GriddedRegion(locs, BorderType.GREAT_CIRCLE, gridSpacing, null);
			System.out.println("Loaded region '"+name+"' with spacing "+gridSpacing+". locs: "+region.getNumLocations());
			regions.add(region);
		}
		return regions;
	}

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws DocumentException, IOException {
		File writeDir =  new File("/home/kevin/OpenSHA/portfolio_lec/region_map_runs"); // TODO
		File outputDir =  new File("/auto/scec-02/kmilner/hazMaps/2012_07_16-porter-map-runs"); // TODO
		int mins = 500;
		int nodes = 1;
		String queue = "nbns";
		int ppn = 8;
		
		List<GriddedRegion> regions = getRegions();
		int tot = 0;
		for (GriddedRegion reg : regions)
			tot += reg.getNumLocations();
		System.out.println("TOTAL: "+tot);
//		LocationList locs = new LocationList();
//		for (GriddedRegion reg : regions)
//			locs.addAll(reg.getNodeList());
//		RegionUtils.locListToKML(locs, "keith_locs.kml", Color.RED);
		
		EpistemicListERF listERF = new UCERF2_TimeDependentEpistemicList();
		listERF.getTimeSpan().setDuration(50);
		int startYear = 2013;
		listERF.updateForecast();
		
		// TODO order?
		AttenRelRef[] imrRefs = { AttenRelRef.CB_2008, AttenRelRef.BA_2008, AttenRelRef.CY_2008, AttenRelRef.AS_2008 };
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps = Lists.newArrayList();
		ParameterList siteParams = new ParameterList();
		for (AttenRelRef imrRef : imrRefs) {
			for (int i=0; i<4; i++) {
				// once for each IMT
				HashMap<TectonicRegionType, ScalarIMR> imrMap = Maps.newHashMap();
				ScalarIMR imr = imrRef.instance(null);
				imr.setParamDefaults();
				for (Parameter<?> param : imr.getSiteParams())
					if (!siteParams.containsParameter(param.getName()))
						siteParams.addParameter(param);
				imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
				imrMaps.add(imrMap);
			}
		}
		List<Parameter<Double>> imts = Lists.newArrayList();
		
		for (int i=0; i<imrMaps.size(); i += 4) {
			ScalarIMR imr = imrMaps.get(i).values().iterator().next();
//			System.out.println(imr.getName());
			
			// PGA
			imr.setIntensityMeasure(PGA_Param.NAME);
			imts.add(imr.getIntensityMeasure());
			
			// SA03
			imr.setIntensityMeasure(SA_Param.NAME);
			SA_Param saParam = (SA_Param)imr.getIntensityMeasure();
			SA_Param.setPeriodInSA_Param(saParam, 0.3);
			imts.add((Parameter<Double>) saParam.clone());
			
			// SA10
			imr.setIntensityMeasure(SA_Param.NAME);
			saParam = (SA_Param) imr.getIntensityMeasure();
			SA_Param.setPeriodInSA_Param(saParam, 1.0);
			imts.add((Parameter<Double>) saParam.clone());
			
			// PGV
			imr.setIntensityMeasure(PGV_Param.NAME);
			imts.add(imr.getIntensityMeasure());
		}
		
		System.out.println("IMRs: "+imrMaps.size()+"\tIMTs: "+imts.size());
		
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		
		List<Site> sites = Lists.newArrayList();
		for (GriddedRegion r : regions) {
			for (Location l : r.getNodeList()) {
				Site s = new Site(l);
				for (Parameter<?> param : siteParams)
					s.addParameter((Parameter<?>) param.clone());
				sites.add(s);
			}
		}
		
		ArrayList<SiteDataValueList<?>> datasList = provs.getAllAvailableData(sites);
		SiteTranslator trans = new SiteTranslator();
		
		for (int i=0; i<sites.size(); i++) {
			List<SiteDataValue<?>> datas = Lists.newArrayList();
			
			for (int j=0; j<datasList.size(); j++)
				datas.add(datasList.get(j).getValue(i));
			
			for (Parameter<?> param : sites.get(i)) {
				boolean success = trans.setParameterValue(param, datas);
				if (!success && (param.getName().equals(DepthTo2pt5kmPerSecParam.NAME)
						|| param.getName().equals(DepthTo2pt5kmPerSecParam.NAME)))
					param.setValue(null); // set to "null" if it is a miss for basin depth
			}
		}
		
		IMT_Info imtInfo = new IMT_Info();
		HashMap<String, DiscretizedFunc> imtXValMap = new HashMap<String, DiscretizedFunc>();
		imtXValMap.put(PGA_Param.NAME, imtInfo.getDefaultHazardCurve(PGA_Param.NAME));
		imtXValMap.put(PGV_Param.NAME, imtInfo.getDefaultHazardCurve(PGV_Param.NAME));
		imtXValMap.put(SA_Param.NAME, imtInfo.getDefaultHazardCurve(SA_Param.NAME));
		CalculationSettings calcSettings = new CalculationSettings(imtXValMap, 200d);
		
		int digits = ((listERF.getNumERFs()-1)+"").length();
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(new File("/home/scec-02/kmilner/hazMaps/svn/dist/OpenSHA_complete.jar"));
		classpath.add(new File("/home/scec-02/kmilner/hazMaps/svn/lib/commons-cli-1.2.jar"));
		
		for (int i=0; i<listERF.getNumERFs(); i++) {
			ERF erf = listERF.getERF(i);
			if (!erf.getTimeSpan().getStartTimePrecision().equals(TimeSpan.NONE)) {
				// can only set time span for BPT, not empirical model
				try {
					erf.getTimeSpan().setStartTime(startYear);
				} catch (RuntimeException e) {
					for (Parameter<?> param : erf.getAdjustableParameterList())
						System.out.println(param.getName()+":\t"+param.getValue());
					throw e;
				}
			}
			
			String numStr = i+"";
			while (numStr.length() < digits)
				numStr = "0"+numStr;
			
			File curveDir = new File(outputDir, "curves_"+numStr);
			CurveResultsArchiver archiver = new AsciiFileCurveArchiver(curveDir.getAbsolutePath(), true, false);
			
			CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, imts, sites, calcSettings, archiver);
			
			Document doc = XMLUtils.createDocumentWithRoot();
			inputs.toXMLMetadata(doc.getRootElement());
			
			String prefix = "inputs_"+numStr;
			String fname = prefix+".xml";
			File remoteFName = new File(outputDir, fname);
			
			XMLUtils.writeDocumentToFile(new File(writeDir, fname), doc);
			// TODO write pbs
			JavaShellScriptWriter javaWriter = new JavaShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, 7000, classpath);
			List<String> script = javaWriter.buildScript(HazardCurveDriver.class.getName(),
					"--threaded "+remoteFName.getAbsolutePath());
			
			File pbsFile = new File(writeDir, prefix+".pbs");
			
			USC_HPCC_ScriptWriter pbsWriter = new USC_HPCC_ScriptWriter();
			pbsWriter.writeScript(pbsFile, script, mins, nodes, ppn, queue);
		}
	}

}
