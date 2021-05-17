package scratch.kevin.miscFigures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CS_Study18_8_BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_CCAi6BasinDepth;
import org.opensha.commons.data.siteData.impl.USGSBayAreaBasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class ChristineSiteDataMapPlot {

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException, GMT_MapException {
		Region reg = new CaliforniaRegions.RELM_TESTING();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.02, null);
		
		WillsMap2015 wills = new WillsMap2015();
		
		GriddedGeoDataSet willsXYZ = build(gridReg, wills);
		
		USGSBayAreaBasinDepth usgs_1p0 = new USGSBayAreaBasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
		USGSBayAreaBasinDepth usgs_2p5 = new USGSBayAreaBasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
		
//		CVM4i26BasinDepth s426_1p0 = new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
//		CVM4i26BasinDepth s426_2p5 = new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
//		
//		CVM_CCAi6BasinDepth cca_1p0 = new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
//		CVM_CCAi6BasinDepth cca_2p5 = new CVM_CCAi6BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
//		
//		GriddedGeoDataSet s426_cca_usgs_1p0 = build(gridReg, s426_1p0, cca_1p0, usgs_1p0);
//		GriddedGeoDataSet s426_cca_usgs_2p5 = build(gridReg, s426_2p5, cca_2p5, usgs_2p5);
//		
//		GriddedGeoDataSet s426_usgs_cca_1p0 = build(gridReg, s426_1p0, usgs_1p0, cca_1p0);
//		GriddedGeoDataSet s426_usgs_cca_2p5 = build(gridReg, s426_2p5, usgs_2p5, cca_2p5);
		
		CS_Study18_8_BasinDepth cs18_8_z1p0 = new CS_Study18_8_BasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
		CS_Study18_8_BasinDepth cs18_8_z2p5 = new CS_Study18_8_BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
		
		GriddedGeoDataSet s426_cca_usgs_1p0 = build(gridReg, cs18_8_z1p0, usgs_1p0);
		GriddedGeoDataSet s426_cca_usgs_2p5 = build(gridReg, cs18_8_z2p5, usgs_2p5);
		
		GriddedGeoDataSet s426_usgs_cca_1p0 = build(gridReg, usgs_1p0, cs18_8_z1p0);
		GriddedGeoDataSet s426_usgs_cca_2p5 = build(gridReg, usgs_2p5, cs18_8_z2p5);
		
		File outputDir = new File("/tmp/christine_site_data");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		String willsPrefix = "wills_2015_vs30";
		plot(outputDir, willsPrefix, "Wills (2015) Vs30 (m/s)", willsXYZ);
		writeCSV(outputDir, willsPrefix, "Vs30 (m/s)", willsXYZ);
		
		String z1_order1_Prefix = "s426_cca_usgs_z1p0";
		plot(outputDir, z1_order1_Prefix, "Stiched Z1.0, S4.26,CCA,USGS (km)", s426_cca_usgs_1p0);
		writeCSV(outputDir, z1_order1_Prefix, "Z1.0 (km)", s426_cca_usgs_1p0);
		
		String z2p5_order1_Prefix = "s426_cca_usgs_z2p5";
		plot(outputDir, z2p5_order1_Prefix, "Stiched Z2.5, S4.26,CCA,USGS (km)", s426_cca_usgs_2p5);
		writeCSV(outputDir, z2p5_order1_Prefix, "Z2.5 (km)", s426_cca_usgs_2p5);
		
		String z1_order2_Prefix = "s426_usgs_cca_z1p0";
		plot(outputDir, z1_order2_Prefix, "Stiched Z1.0, S4.26,USGS,CGS (km)", s426_usgs_cca_1p0);
		writeCSV(outputDir, z1_order2_Prefix, "Z1.0 (km)", s426_usgs_cca_1p0);
		
		String z2p5_order2_Prefix = "s426_usgs_cca_z2p5";
		plot(outputDir, z2p5_order2_Prefix, "Stiched Z2.5, S4.26,USGS,CGS (km)", s426_usgs_cca_2p5);
		writeCSV(outputDir, z2p5_order2_Prefix, "Z2.5 (km)", s426_usgs_cca_2p5);
	}
	
	@SuppressWarnings("unchecked")
	private static GriddedGeoDataSet build(GriddedRegion gridReg, SiteData<Double>... provs)
			throws IOException {
		List<List<Double>> valsLists = new ArrayList<>();
		
		for (SiteData<Double> prov : provs)
			valsLists.add(prov.getValues(gridReg.getNodeList()));
		
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<xyz.size(); i++) {
			double val = Double.NaN;
			for (int j=0; j<provs.length; j++) {
				val = valsLists.get(j).get(i);
				if (provs[j].isValueValid(val))
					break;
			}
			xyz.set(i, val);
		}
		
		return xyz;
	}
	
	private static final boolean TOPOGRAPHY = false;
	
	private static void plot(File outputDir, String prefix, String label, GriddedGeoDataSet xyz)
			throws GMT_MapException, IOException {
		double minFinite = Double.POSITIVE_INFINITY;
		double maxFinite = 0d;
		for (int i=0; i<xyz.size(); i++) {
			double val = xyz.get(i);
			if (Double.isFinite(val)) {
				minFinite = Math.min(minFinite, val);
				maxFinite = Math.max(maxFinite, val);
			}
		}
		System.out.println(label+" Range: "+minFinite+" "+maxFinite);
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(minFinite, maxFinite);
		GMT_Map map = new GMT_Map(xyz.getRegion(), xyz, xyz.getRegion().getSpacing(), cpt);
		map.setCustomLabel(label);
		if (TOPOGRAPHY) {
			map.setTopoResolution(TopographicSlopeFile.US_SIX);
		} else {
			map.setTopoResolution(null);
			map.setTopoResolution(null);
			map.setUseGMTSmoothing(false);
		}
		map.setLogPlot(false);
		map.setDpi(300);
		map.setXyzFileName("site_data.xyz");
		map.setCustomScaleMin(minFinite);
		map.setCustomScaleMax(maxFinite);
		map.setBlackBackground(false);
		map.setRescaleCPT(false);
		map.setJPGFileName(null);
		map.setPDFFileName(null);
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
	}
	
	private static void writeCSV(File outputDir, String prefix, String value, GriddedGeoDataSet xyz)
			throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Longitude", "Latitude", value);
		for (int i=0; i<xyz.size(); i++) {
			Location loc = xyz.getLocation(i);
			csv.addLine((float)loc.getLongitude()+"", (float)loc.getLatitude()+"", (float)xyz.get(i)+"");
		}
		csv.writeToFile(new File(outputDir, prefix+".csv"));
	}

}
