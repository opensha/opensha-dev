package scratch.kevin.cybershake;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_CCAi6BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.data.xyz.AbstractGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.cybershake.ModProbConfig;
import org.opensha.sha.cybershake.bombay.ModProbConfigFactory;
import org.opensha.sha.cybershake.bombay.ScenarioBasedModProbConfig;
import org.opensha.sha.cybershake.maps.CyberShake_GMT_MapGenerator;
import org.opensha.sha.cybershake.maps.GMT_InterpolationSettings;
import org.opensha.sha.cybershake.maps.HardCodedInterpDiffMapCreator;
import org.opensha.sha.cybershake.maps.InterpDiffMap;
import org.opensha.sha.cybershake.maps.ProbGainCalc;
import org.opensha.sha.cybershake.maps.InterpDiffMap.InterpDiffMapType;
import org.opensha.sha.cybershake.maps.servlet.CS_InterpDiffMapServletAccessor;

import com.google.common.collect.Lists;

public class BatchBaseMapPlot {

	public static void main(String[] args) throws Exception {
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_08-cvm4-cs-nga2-1sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_08-cvm4-cs-nga-1sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_08-cvm4-cs-nga2-pga");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_08-cvm4-cs-nga-pga");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_11-cvm4-cs-nga2-5sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2013_11_11-cvm4-cs-nga-5sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_03-bbp-cs-nga-3sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_03-cvm4i26-cs-nga-3sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_05-cvmhnogtl-cs-nga-3sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_17-cvm4i26-cs-nga2-3sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga2-1sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga2-5sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga2-pga");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga-1sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga-5sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2014_03_18-cvm4i26-cs-nga-pga");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2015_05_27-cvm4i26-cs-nga-2sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2015_06_12-cvm4i26-cs-nga-10sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2016_08_31-ccai6-cs-nga2avg-3sec");
//		File dir = new File("/home/kevin/CyberShake/baseMaps/2017_04_10-ccai6-cs-nga2-2sec");
		File dir = new File("/home/kevin/CyberShake/baseMaps/2017_04_11-cca-nobasin-cs-nga2-2sec");
		
		boolean ratios = false;
		
		List<SiteData<Double>> siteDatas = Lists.newArrayList();
		siteDatas.add(new WillsMap2015());
//		siteDatas.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
//		siteDatas.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		siteDatas.add(new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		siteDatas.add(new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
		
//		Region region = new CaliforniaRegions.CYBERSHAKE_MAP_REGION();
		Region region = new CaliforniaRegions.CYBERSHAKE_CCA_MAP_REGION();
		
//		String imtFileLabel = "1sec";
//		String label = "1sec SA";
//		Double customMax = 3d; // for 1 sec
		
		String imtFileLabel = "2sec";
		String label = "2sec SA";
		Double customMax = 1.0; // for 2 sec
		
//		String imtFileLabel = "pga";
//		String label = "PGA";
//		Double customMax = 3d; // for PGA
		
//		String imtFileLabel = "3sec";
//		String label = "3sec SA";
//		Double customMax = 1.4; // for 3 sec
		
//		String imtFileLabel = "5sec";
//		String label = "5sec SA";
//		Double customMax = 1d; // for 5 sec
		
//		String imtFileLabel = "10sec";
//		String label = "10sec SA";
//		Double customMax = 0.6d; // for 10 sec
		
		boolean isProbAt_IML = false;
		double val = 0.0004;
		String probLabel = "2% in 50 yrs";
		String probFileLabel = "2p_in_50";
//		double val = 0.0002;
//		String probLabel = "1% in 50 yrs";
//		String probFileLabel = "1p_in_50";
//		double val = 0.0001;
//		String probLabel = "1% in 100 yrs";
//		String probFileLabel = "1p_in_100";
		
		label += ", "+probLabel;
		
		Double customMin = 0d;
		
		List<GeoDataSet> maps = Lists.newArrayList();
		
		List<String> names = Lists.newArrayList();
		
		for (File subDir : dir.listFiles()) {
			// this will be the IMR name
			if (!subDir.isDirectory())
				continue;
			
			String name = subDir.getName();
			
			File binFile = new File(new File(subDir, "curves"), "imrs1.bin");
			if (!binFile.exists())
				continue;
			
			File outputFile = new File(dir, name.toLowerCase()+"_"+imtFileLabel+"_"+probFileLabel+".png");
			boolean skip = outputFile.exists();
			if (!ratios && skip)
				continue;
			
			System.out.println("Loading "+name);
			
			BinaryHazardCurveReader reader = new BinaryHazardCurveReader(binFile.getAbsolutePath());
			Map<Location, ArbitrarilyDiscretizedFunc> curves = reader.getCurveMap();
			
			GeoDataSet baseMap = HazardDataSetLoader.extractPointFromCurves(curves, isProbAt_IML, val);
			maps.add(baseMap);
			names.add(name);
			if (skip)
				continue;
//			int nanCount = 0;
//			for (int i=0; i<baseMap.size(); i++) {
//				if (Double.isNaN(baseMap.get(i))) {
//					System.out.println("NaN!");
//					System.out.print("X:");
//					for (Point2D pt : curves.get(baseMap.getLocation(i)))
//						System.out.print((float)pt.getX());
//					System.out.println();
//					System.out.print("Y:");
//					for (Point2D pt : curves.get(baseMap.getLocation(i)))
//						System.out.print((float)pt.getY());
//					System.out.println();
//					System.out.println(baseMap.getLocation(i));
//					System.exit(0);
//					nanCount++;
//				}
//			}
//			System.out.println("NaN count: "+nanCount);
			
			String metadata = "isProbAt_IML: " + isProbAt_IML + "\n" +
					"val: " + val + "\n";
			
			System.out.println("Plotting "+name);
			plot(outputFile, baseMap, region, customMin, customMax, name+" "+label, metadata, false);
		}
		
		// now site data
		for (SiteData<Double> data : siteDatas)
			checkMakeSiteDataPlot(data, region, dir);
		
		if (ratios && maps.size() > 1) {
			for (int i=0; i<maps.size(); i++) {
				GeoDataSet map1 = maps.get(i);
				String name1 = names.get(i);
				for (int j=i+1; j<maps.size(); j++) {
					GeoDataSet map2 = maps.get(j);
					String name2 = names.get(j);
					
					String fName = "ratio_"+name1+"_"+name2+"_"+probFileLabel+".png";
					File outputFile = new File(dir, fName);
					if (outputFile.exists())
						continue;
					
					System.out.println("Calculating ratio: "+name1+" vs "+name2);
					
					GeoDataSet ratio = GeoDataSetMath.divide(map1, map2);
					int numIdentical = 0;
					for (int k=0; k<ratio.size(); k++)
						if (ratio.get(k) == 1)
							numIdentical++;
					if (numIdentical == ratio.size())
						System.out.println("They're identical!");
					else
						System.out.println(numIdentical+"/"+ratio.size()+" numerically identical");
					
					System.out.println("Plotting ratio");
					
					plot(outputFile, ratio, region, null, null, label+" Ratio", "asdf", true);
				}
			}
		}
	}
	
	private static CPT getCPT(boolean ratio) throws IOException {
		if (ratio)
			return CyberShake_GMT_MapGenerator.getRatioCPT();
		else
//			return CPT.loadFromStream(HardCodedInterpDiffMapCreator.class.getResourceAsStream(
//				"/resources/cpt/MaxSpectrum2.cpt"));
			return CyberShake_GMT_MapGenerator.getHazardCPT();
	}
	
	private static void plot(File outputFile, GeoDataSet baseMap, Region region,
			Double customMin, Double customMax, String label, String metadata, boolean ratio)
					throws IOException, ClassNotFoundException {
		plot(outputFile, baseMap, region, customMin, customMax, label, metadata, getCPT(ratio), !ratio);
	}
	
	private static void plot(File outputFile, GeoDataSet baseMap, Region region,
			Double customMin, Double customMax, String label, String metadata, CPT cpt, boolean rescaleCPT)
					throws IOException, ClassNotFoundException {
		
		double baseMapRes = 0.005;
		System.out.println("Loading basemap...");
		
		System.out.println("Fetching curves...");
		AbstractGeoDataSet scatterData = null;
		
		System.out.println("Creating map instance...");
		GMT_InterpolationSettings interpSettings = GMT_InterpolationSettings.getDefaultSettings();
		
		InterpDiffMapType[] mapTypes = {InterpDiffMapType.BASEMAP};
		
		InterpDiffMap map = new InterpDiffMap(region, baseMap, baseMapRes, cpt, scatterData, interpSettings, mapTypes);
		map.setCustomLabel(label);
		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setLogPlot(false);
		map.setDpi(300);
		map.setXyzFileName("base_map.xyz");
		map.setCustomScaleMin(customMin);
		map.setCustomScaleMax(customMax);
		map.setRescaleCPT(rescaleCPT);
		
		System.out.println("Making map...");
		String url = CS_InterpDiffMapServletAccessor.makeMap(null, map, metadata);
		
		FileUtils.downloadURL(url+"/basemap.150.png", outputFile);
	}
	
	private static void checkMakeSiteDataPlot(SiteData<Double> prov, Region region, File dir) throws IOException, ClassNotFoundException {
		String shortType;
		Double customMin, customMax;
		if (prov.getDataType().equals(SiteData.TYPE_VS30)) {
			customMax = 1000d;
			customMin = 180d;
			shortType = "Vs30";
		} else if (prov.getDataType().equals(SiteData.TYPE_DEPTH_TO_1_0)) {
			customMax = 2d;
			customMin = 0d;
			shortType = "Z1.0";
		} else if (prov.getDataType().equals(SiteData.TYPE_DEPTH_TO_2_5)) {
			customMax = 10d;
			customMin = 0d;
			shortType = "Z2.5";
		} else {
			throw new IllegalStateException();
		}
		
		File outputFile = new File(dir, prov.getShortName()+"_"+shortType+".png");
		if (outputFile.exists())
			return;
		
		GriddedRegion gridReg = new GriddedRegion(region, 0.005, null);
		List<Double> vals = prov.getValues(gridReg.getNodeList());
		GriddedGeoDataSet data = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<vals.size(); i++)
			data.set(i, vals.get(i));
		
		CPT cpt = CPT.loadFromStream(HardCodedInterpDiffMapCreator.class.getResourceAsStream(
				"/resources/cpt/MaxSpectrum2.cpt"));
		
		plot(outputFile, data, region, customMin, customMax, prov.getShortName()+" "+shortType, prov.getMetadata(), cpt, true);
	}

}
