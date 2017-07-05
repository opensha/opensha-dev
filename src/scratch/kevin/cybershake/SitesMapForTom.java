package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.impl.SRTM30PlusTopography;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol.Symbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbolSet;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.SiteInfo2DB;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class SitesMapForTom {

	public static void main(String[] args) throws IOException, GMT_MapException {
		Region laReg = new CaliforniaRegions.CYBERSHAKE_MAP_REGION();
		Region ccaReg = new CaliforniaRegions.CYBERSHAKE_CCA_MAP_REGION();
		
		DBAccess db = new DBAccess("moment.usc.edu", "CyberShake");
		SiteInfo2DB sites2db = new SiteInfo2DB(db);
		List<CybershakeSite> allSites = sites2db.getAllSitesFromDB();
		RunIDFetcher laRunIDFetch = new RunIDFetcher(db, 61, 21);
		HashSet<String> ccaSiteNames = new HashSet<String>();
		CSVFile<String> csv = CSVFile.readFile(new File("/home/kevin/CyberShake/study_16_9_sites.csv"), true);
		for (int row=0; row<csv.getNumRows(); row++)
			ccaSiteNames.add(csv.get(row, 2));
		
		HashSet<CybershakeSite> laSites = new HashSet<CybershakeSite>();
		HashSet<CybershakeSite> ccaSites = new HashSet<CybershakeSite>();
		HashSet<CybershakeSite> bothSites = new HashSet<CybershakeSite>();
		for (CybershakeSite site : allSites) {
			Location loc = site.createLocation();
			boolean la = laReg.contains(loc) && laRunIDFetch.getSite(site.short_name) != null;
			boolean cca = ccaReg.contains(loc) && ccaSiteNames.contains(site.short_name);
			
			if (la)
				laSites.add(site);
			if (cca)
				ccaSites.add(site);
			if (la && cca)
				bothSites.add(site);
//			if (site.type_id == CybershakeSite.TYPE_GRID_05_KM) {
//				// only continue if run exists
//				if (laRunIDFetch.getSite(site.short_name) == null)
//					continue;
//			} else if (site.type_id == CybershakeSite.TYPE_TEST_SITE) {
//				continue;
//			}
//			
//			if (laReg.contains(site.createLocation()))
//				laSites.add(site);
//			if (ccaReg.contains(site.createLocation())) {
//				ccaSites.add(site);
//				if (laRunIDFetch.getSite(site.short_name) != null)
//					bothSites.add(site);
//			}
		}
		int numLA = laSites.size();
		int numCCA = ccaSites.size();
		int numBoth = bothSites.size();
		
		System.out.println("LA: "+numLA);
		System.out.println("CCA: "+numCCA);
		System.out.println("Both: "+numBoth);
		
		db.destroy();
		
//		double maxLat = Double.max(laReg.getMaxLat(), ccaReg.getMaxLat());
//		double minLat = Double.min(laReg.getMinLat(), ccaReg.getMinLat());
//		double maxLon = Double.max(laReg.getMaxLon(), ccaReg.getMaxLon());
//		double minLon = Double.min(laReg.getMinLon(), ccaReg.getMinLon());
		double maxLat = 38;
		double minLat = 32.5;
		double maxLon = -116;
		double minLon = -123;
		
		Region plotRegion = new Region(new Location(maxLat, maxLon), new Location(minLat, minLon));
		
		double spacing = 0.01;
		GriddedRegion gridReg = new GriddedRegion(new Location(maxLat+spacing, maxLon+spacing),
				new Location(minLat-spacing, minLon-spacing), spacing, null);
		
		// fetch topography
		System.out.println("Fetching topography");
		SRTM30PlusTopography topo = new SRTM30PlusTopography();
		LocationList nodes = gridReg.getNodeList();
		ArrayList<Double> vals = topo.getValues(nodes);
		
		GriddedGeoDataSet topoXYZ = new GriddedGeoDataSet(gridReg, false);
		
		for (int i=0; i<nodes.size(); i++)
			topoXYZ.set(i, vals.get(i));
		
//		CPT cpt = GMT_CPT_Files.GMT_RELIEF.instance();
		CPT cpt = CPT.loadFromFile(new File("/home/kevin/CyberShake/dem.cpt"));
		for (CPTVal val : cpt) {
			val.minColor = saturate(val.minColor);
			val.maxColor = saturate(val.maxColor);
		}
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setAboveMaxColor(cpt.getMaxColor());
		cpt = cpt.rescale(0d, 3000d);
		
		GMT_Map map = new GMT_Map(plotRegion, topoXYZ, spacing, cpt);
		map.setRescaleCPT(false);
		
		System.out.println("Building symbols");
		PSXYSymbolSet symbols = new PSXYSymbolSet();
		Color laColor = Color.BLUE.darker();
		Color ccaColor = Color.BLACK;
//		Color ccaColor = Color.WHITE;
		Color bothColor = Color.RED.darker();
		double laVal = 0d;
		double ccaVal = 1d;
		double bothVal = 2d;
		CPT symbolCPT = new CPT(laVal, bothVal, laColor, ccaColor, bothColor);
		symbols.setCpt(symbolCPT);
		Symbol symbol = Symbol.INVERTED_TRIANGLE;
		float width = 0.07f;
		for (CybershakeSite site : laSites) {
			if (bothSites.contains(site))
				continue;
			Point2D pt = new Point2D.Double(site.lon, site.lat);
			symbols.addSymbol(new PSXYSymbol(pt, symbol, width), laVal);
		}
		for (CybershakeSite site : ccaSites) {
			if (bothSites.contains(site))
				continue;
			Point2D pt = new Point2D.Double(site.lon, site.lat);
			symbols.addSymbol(new PSXYSymbol(pt, symbol, width), ccaVal);
		}
		for (CybershakeSite site : bothSites) {
			Point2D pt = new Point2D.Double(site.lon, site.lat);
			symbols.addSymbol(new PSXYSymbol(pt, symbol, width), bothVal);
		}
		
		map.setSymbolSet(symbols);
		
		double lineWidth = 1.2;
		System.out.println("Adding outlines");
		addOutline(laReg, laColor, lineWidth, map);
		addOutline(ccaReg, ccaColor, lineWidth, map);
		
		map.setBlackBackground(false);
		
		System.out.println("Plotting map");
		FaultBasedMapGen.plotMap(new File("/tmp"), "cs_sites", true, map);
	}
	
	private static Color saturate(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		
		int saturationSteps = 1;
		
		for (int i=0; i<saturationSteps; i++) {
			r = (int)(0.5d*(r + 255d)+0.5);
			g = (int)(0.5d*(g + 255d)+0.5);
			b = (int)(0.5d*(b + 255d)+0.5);
		}
		
		return new Color(r, g, b);
	}
	
	private static void addOutline(Region reg, Color c, double width, GMT_Map map) {
		LocationList border = reg.getBorder();
		LocationList borderCompleted = new LocationList();
		borderCompleted.addAll(border);
		borderCompleted.add(border.get(0));
		PSXYPolygon poly = new PSXYPolygon(borderCompleted);
		poly.setFillColor(null);
		poly.setPenColor(c);
		poly.setPenWidth(width);
		map.addPolys(poly);
	}

}