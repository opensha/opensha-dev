package scratch.kevin;

import java.awt.Color;
import java.io.IOException;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.GMT_Map.HighwayFile;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.cybershake.maps.GMT_InterpolationSettings;
import org.opensha.sha.cybershake.maps.InterpDiffMap;
import org.opensha.sha.cybershake.maps.InterpDiffMap.InterpDiffMapType;
import org.opensha.sha.cybershake.maps.servlet.CS_InterpDiffMapServletAccessor;

public class PatrickZ2500Plotter {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws RuntimeException 
	 * @throws GMT_MapException 
	 * @throws ClassNotFoundException 
	 */
	public static void main(String[] args) throws IOException, GMT_MapException, RuntimeException, ClassNotFoundException {
		GMT_MapGenerator gen = new GMT_MapGenerator();
		gen.getAdjustableParamsList().getParameter(GMT_MapGenerator.LOG_PLOT_NAME).setValue(false);
		
		Region region = new CaliforniaRegions.CYBERSHAKE_MAP_REGION();
		region =
			new Region(new Location(region.getMaxLat(), region.getMinLon()),
					new Location(region.getMinLat(), region.getMaxLon()));
		
		GeoDataSet data1 = ArbDiscrGeoDataSet.loadXYZFile("/home/kevin/OpenSHA/ucvm/ucvm_cvmh_z2500.data", false);
		GeoDataSet data2 = ArbDiscrGeoDataSet.loadXYZFile("/home/kevin/OpenSHA/ucvm/ucvm_cvms_z2500.data", false);
		
//		double min = data1.getMinZ();
//		double max = data1.getMaxZ();
//		double min2 = data2.getMinZ();
//		double max2 = data2.getMaxZ();
//		if (min2 < min)
//			min = min2;
//		if (max2 > max)
//			max = max2;
		double min = 0;
		double max = 6000;
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		
//		gen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CUSTOM_SCALE_LABEL_PARAM_CHECK_NAME).setValue(true);
//		gen.getAdjustableParamsList().getParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME).setValue(min);
//		gen.getAdjustableParamsList().getParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME).setValue(max);
//		gen.getAdjustableParamsList().getParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME).setValue(0.01);
		
		GMT_InterpolationSettings interpSettings = GMT_InterpolationSettings.getDefaultSettings();
		InterpDiffMapType[] mapTypes = {InterpDiffMapType.BASEMAP};
		
//		GMT_Map map1 = gen.getGMTMapSpecification(data1);
//		GMT_Map map2 = gen.getGMTMapSpecification(data2);
//		GMT_Map map1 = new GMT_Map(region, data1, 0.01, cpt);
		InterpDiffMap map1 = new InterpDiffMap(region, data1, 0.01, cpt, null, interpSettings, mapTypes);
		map1.setCustomLabel("Z2500 (m), CVM-H11.2");
		setupMap(map1, min, max);
//		GMT_Map map2 = new GMT_Map(region, data2, 0.01, cpt);
		InterpDiffMap map2 = new InterpDiffMap(region, data2, 0.01, cpt, null, interpSettings, mapTypes);
		map2.setCustomLabel("Z2500 (m), CVM-S4");
		setupMap(map2, min, max);
		
//		System.out.println(gen.makeMapUsingServlet(map1, "", null).replace("jpg", "png"));
//		System.out.println(gen.makeMapUsingServlet(map2, "", null).replace("jpg", "png"));
		
		System.out.println(CS_InterpDiffMapServletAccessor.makeMap(null, map1, ""));
		System.out.println(CS_InterpDiffMapServletAccessor.makeMap(null, map2, ""));
	}
	
	private static void setupMap(GMT_Map map, double min, double max) {
		map.setCustomScaleMin(min);
		map.setCustomScaleMax(max);
		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
//		map.setDpi(300);
		map.setHighwayFile(HighwayFile.ALL);
		map.setCoast(new CoastAttributes(Color.WHITE, 2d));
//		map.set
	}

}
