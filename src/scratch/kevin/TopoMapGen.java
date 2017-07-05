package scratch.kevin;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.siteData.impl.SRTM30PlusTopography;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol.Symbol;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class TopoMapGen {

	public static void main(String[] args) throws IOException, GMT_MapException {
		Region region = new Region(new Location(35, -119), new Location(33, -115));
		double spacing = 0.01;
		boolean forcePositive = true;
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/tom_rancho_mirage");
		String prefix = "topo";
		
		List<Location> annotations = Lists.newArrayList();
		annotations.add(new Location(33.739683, -116.412925));
		
		List<FaultSectionPrefData> faults = FaultModels.FM3_1.fetchFaultSections();
		float faultThickness = 1f;
		
		GriddedRegion gridReg = new GriddedRegion(region, spacing, null);
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		SRTM30PlusTopography fetch = new SRTM30PlusTopography();
		System.out.println("Fetching Data");
		ArrayList<Double> data = fetch.getValues(xyz.getLocationList());
		double minVal = forcePositive ? 0 : Double.NEGATIVE_INFINITY;
		for (int i=0; i<xyz.size(); i++)
			xyz.set(i, Math.max(data.get(i), minVal));
		
		CPT cpt = GMT_CPT_Files.GMT_GLOBE.instance();
		if (forcePositive) {
			for (int i=cpt.size(); --i>=0;)
				if (cpt.get(i).start < 0f)
					cpt.remove(i);
		}
		
		GMT_Map map = new GMT_Map(region, xyz, spacing, cpt);
		
		map.setLogPlot(false);
//		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setBlackBackground(false);
		map.setCustomScaleMin((double)cpt.getMinValue());
		map.setCustomScaleMax((double)cpt.getMaxValue());
		map.setCustomLabel("Elevation (m)");
		map.setRescaleCPT(false);
		map.setJPGFileName(null);
		
		if (faults != null) {
			for (FaultSectionPrefData sect : faults) {
				for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), Color.BLACK, faultThickness))
					map.addPolys(poly);
			}
		}
		
		if (annotations != null) {
			for (Location loc : annotations) {
				java.awt.geom.Point2D.Double pt = new Point2D.Double(loc.getLongitude(), loc.getLatitude());
				map.addSymbol(new PSXYSymbol(pt, Symbol.INVERTED_TRIANGLE, 0.1f, 0f, null, Color.BLACK));
			}
		}
		
		System.out.println("Making Map");
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
		System.out.println("DONE");
	}

}
