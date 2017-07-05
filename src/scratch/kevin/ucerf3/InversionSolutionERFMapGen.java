package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.CVM4BasinDepth;
import org.opensha.commons.data.siteData.impl.USGSBayAreaBasinDepth;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.AsciiFileCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.SiteTranslator;
import org.opensha.sha.util.TectonicRegionType;

import scratch.UCERF3.oldStuff.OldInversionSolutionERF;

public class InversionSolutionERFMapGen {
	
	public static void main(String[] args) throws IOException {
//		File solutionFile = new File("/home/scec-02/kmilner/ucerf3/" +
//				"inversion_solutions/dsa_4threads_50nodes_FAST_SA_dSub200_sub100_run3.xml");
//		File outputDir = new File("/home/scec-02/kmilner/hazMaps/ucerf3_inv_state_run_1/curves");
		File solutionFile = new File("/home/kevin/2011_09_08-morgan-CS_fixed.zip");
		File outputDir = new File("/tmp/hazMapTest");
		
		double spacing = 0.2;
		
		OldInversionSolutionERF erf = new OldInversionSolutionERF();
		erf.setSolutionFile(solutionFile);
		
		ScalarIMR imr = new CB_2008_AttenRel(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), 1.0);
		
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps =
			new ArrayList<Map<TectonicRegionType,ScalarIMR>>();
		HashMap<TectonicRegionType, ScalarIMR> imrMap = new HashMap<TectonicRegionType, ScalarIMR>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
		imrMaps.add(imrMap);
		
		GriddedRegion region = new CaliforniaRegions.RELM_TESTING_GRIDDED(spacing);
		ArrayList<SiteData<?>> provList = new ArrayList<SiteData<?>>();
		provList.add(new WillsMap2006());
		provList.add(new WaldAllenGlobalVs30());
		provList.add(new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		provList.add(new USGSBayAreaBasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		
		System.out.println("Region has "+region.getNodeCount()+" sites");
		
		
		ArrayList<Site> sites = new ArrayList<Site>();
		for (Location loc : region) {
			Site site = new Site(loc);
			
			site.addParameter((Parameter)imr.getParameter(Vs30_Param.NAME).clone());
			site.addParameter((Parameter)imr.getParameter(DepthTo2pt5kmPerSecParam.NAME).clone());
			
			sites.add(site);
		}
		OrderedSiteDataProviderList provs = new OrderedSiteDataProviderList(provList);
		System.out.println("Fetching data!");
		ArrayList<SiteDataValueList<?>> vals  = provs.getAllAvailableData(sites);
		System.out.println("done.");
		
		SiteTranslator trans = new SiteTranslator();
		for (int i=0; i<sites.size(); i++) {
			ArrayList<SiteDataValue<?>> siteVals = new ArrayList<SiteDataValue<?>>();
			for (SiteDataValueList<?> valList : vals) {
				siteVals.add(valList.getValue(i));
			}
			Iterator<Parameter<?>> it = sites.get(i).getParametersIterator();
			while (it.hasNext()) {
				trans.setParameterValue(it.next(), siteVals);
			}
		}
		
		CalculationSettings calcSettings = new CalculationSettings(IMT_Info.getUSGS_SA_Function(), 200);
		calcSettings.setSerializeERF(false);
		AsciiFileCurveArchiver archiver = new AsciiFileCurveArchiver(outputDir.getAbsolutePath(), true, false);
		
		CalculationInputsXMLFile xml = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
		System.out.println("Writing XML!");
		XMLUtils.writeObjectToXMLAsRoot(xml, new File("/tmp/hazard_map.xml"));
		System.out.println("DONE!");
	}

}
