package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.data.siteData.impl.SRTM30PlusTopography;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.PSText;
import org.opensha.commons.mapping.gmt.elements.PSText.Justify;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol.Symbol;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public class SiteMapGenerator {

//	private static Region region = new Region(new Location(35.1, -114.5), new Location(32, -120));
//	private static Region region = new Region(new Location(35, -115), new Location(33, -120));
//	private static Region region = new Region(new Location(35, -116), new Location(33.25, -119.5));
	private static Region region = new Region(new Location(34.8, -116.7), new Location(33.6, -119.5));
	
	private static GriddedGeoDataSet fetchTopo(TopographicSlopeFile topoRes) throws IOException {
		SRTM30PlusTopography topo = new SRTM30PlusTopography();
		double discr = Math.max((double)topoRes.resolution()/3600d, topo.getResolution());
		System.out.println("Topo discretization: "+discr);
		GriddedRegion gridReg = new GriddedRegion(region, discr, null);
		System.out.println("Fetching "+gridReg.getNodeCount()+" topo values");
		ArrayList<Double> vals = new SRTM30PlusTopography().getValues(gridReg.getNodeList());
		GriddedGeoDataSet topoXYZ = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<topoXYZ.size(); i++)
			topoXYZ.set(i, vals.get(i));
		return topoXYZ;
	}
	
	public static void main(String[] args) throws IOException, GMT_MapException {
		File outputDir = new File("/tmp");
		String prefix = "site_map";
		FaultBasedMapGen.LOCAL_MAPGEN = true;
//	private static void plotMap(File outputDir, String prefix, double arraySpacing, FaultModels fm,
//			MeanUCERF3 meanU3, ScalarType scalarType, boolean plotLinearTrace) throws IOException, GMT_MapException {
//		TopographicSlopeFile topoRes = TopographicSlopeFile.CA_THREE;
//		GriddedGeoDataSet topoXYZ = fetchTopo(topoRes);
		GriddedGeoDataSet topoXYZ = null;
		CPT topoCPT = new CPT(-100d, 3400, Color.WHITE, Color.DARK_GRAY);
		GMT_Map map = new GMT_Map(region, topoXYZ, 3d/3600d, topoCPT);
		
		FaultModels fm = FaultModels.FM2_1;
		
		map.setBlackBackground(false);
		map.setRescaleCPT(false);
		map.setCustomScaleMin((double)topoCPT.getMinValue());
		map.setCustomScaleMax((double)topoCPT.getMaxValue());
//		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(160, 200, 200)));
		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(112, 131, 147)));
		map.setCustomLabel(null);
		map.setUseGMTSmoothing(true);
//		map.setTopoResolution(topoRes);
		map.setCustomIntenPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grad");
		map.setCustomGRDPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grd");
//		map.setLabelSize(MAP_LABEL_SIZE);
//		map.setLabelTickSize(MAP_LABEL_TICK_SIZE);
//		map.setlabel
		map.setDpi(150);
//		map.setImageWidth(imageWidth);
		map.setDrawScaleKM(false);
		map.setGMT_Param("MAP_FRAME_TYPE", "FANCY");
		map.setGMT_Param("MAP_FRAME_WIDTH", "0.04i");
		map.setGMT_Param("MAP_TICK_LENGTH_PRIMARY", "0.04i");
		
		for (FaultSectionPrefData sect : fm.fetchFaultSections())
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), new Color(100, 100, 100), 0.3))
				map.addPolys(poly);
		
		Table<String, Justify, PSXYSymbol> cities = HashBasedTable.create();
		// San Diego, Los Angeles, Santa Barbara, Victorville, Palm Springs
//		cities.put("San Diego", new Location(32.70, -117.15));
		double symbolWidth = 0.07;
		int fontSize = 16;
		Color fontColor = Color.BLACK;
//		Color fontColor = new Color(60, 15, 15);
//		cities.put("Los Angeles", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-118.25, 34.05),
//				Symbol.CIRCLE, symbolWidth, 0d, null, Color.BLACK));
//		cities.put("Palm Springs", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-116.545593, 33.830368),
//				Symbol.CIRCLE, symbolWidth, 0d, null, Color.BLACK));
//		cities.put("STNI", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-118.17881, 33.93088),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("PDU", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-117.63808, 34.1207),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("Bombay Beach M6", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-115.728, 33.3172),
//				Symbol.STAR, symbolWidth*1.5, 0d, null, Color.RED.darker()));
//		cities.put("Mojave M6", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-117.80177, 34.42295),
//				Symbol.STAR, symbolWidth*1.5, 0d, null, Color.RED.darker()));
		

//		cities.put("STNI", Justify.LEFT_TOP, new PSXYSymbol(new Point2D.Double(-118.17881, 33.93088),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("WNGC", Justify.LEFT_TOP, new PSXYSymbol(new Point2D.Double(-118.0653, 34.041824),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("USC", Justify.RIGHT_BOTTOM, new PSXYSymbol(new Point2D.Double(-118.286, 34.0192),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("SBSM", Justify.LEFT_TOP, new PSXYSymbol(new Point2D.Double(-117.29201, 34.064987),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
//		cities.put("SMCA", Justify.RIGHT_TOP, new PSXYSymbol(new Point2D.Double(-118.48939, 34.00909),
//				Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
		for (BBP_Site site : RSQSimBBP_Config.getCyberShakeVs500LASites()) {
			Location loc = site.getLoc();
			Justify just = Justify.LEFT_TOP;
			if (site.getName().equals("SMCA"))
				just = Justify.RIGHT_TOP;
			cities.put(site.getName(), just, new PSXYSymbol(
					new Point2D.Double(loc.getLongitude(), loc.getLatitude()),
					Symbol.INVERTED_TRIANGLE, symbolWidth*1.5, 0d, null, Color.BLUE.darker()));
		}
		fontSize = 14;
		
		for (Cell<String, Justify, PSXYSymbol> cell : cities.cellSet()) {
			String city = cell.getRowKey();
			PSXYSymbol symbol = cell.getValue();
			Point2D pt = symbol.getPoint();
			Justify justify = cell.getColumnKey();
			double addX = 0;
			if (justify == Justify.LEFT || justify == Justify.LEFT_TOP || justify == Justify.LEFT_BOTTOM)
				addX = 0.04;
			if (justify == Justify.RIGHT || justify == Justify.RIGHT_TOP || justify == Justify.RIGHT_BOTTOM)
				addX = -0.04;
			Point2D textPT = new Point2D.Double(pt.getX()+addX, pt.getY());
			
			map.addSymbol(symbol);
			map.addText(new PSText(textPT, fontColor, fontSize, city, justify));
		}
		
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
	}

}
