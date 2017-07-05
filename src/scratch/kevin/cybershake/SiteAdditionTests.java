package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.SRTM30PlusTopoSlope;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol.Symbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbolSet;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.maps.CyberShake_GMT_MapGenerator;

import scratch.UCERF3.analysis.FaultBasedMapGen;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class SiteAdditionTests {

	public static void main(String[] args) throws IOException, GMT_MapException {
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
		
		List<CybershakeSite> allSites = sites2db.getAllSitesFromDB();
		Map<String, CybershakeSite> sitesByName = Maps.newHashMap();
		for (CybershakeSite site : allSites)
			sitesByName.put(site.short_name, site);
		
		Region csReg = new CaliforniaRegions.CYBERSHAKE_MAP_REGION();
		
		List<CybershakeSite> sitesInRegion = Lists.newArrayList();
		
		for (CybershakeSite site : allSites)
			if (csReg.contains(site.createLocation()))
				sitesInRegion.add(site);
		
		System.out.println(sitesInRegion.size()+" total CS sites in CS Region");
		
		List<CybershakeSite> currentSitesInRegion = Lists.newArrayList();
		HazardCurveFetcher fetch = new HazardCurveFetcher(db, 35, 21);
		for (CybershakeSite site : fetch.getCurveSites()) {
			if (site.type_id != CybershakeSite.TYPE_TEST_SITE && csReg.contains(site.createLocation()))
				currentSitesInRegion.add(site);
		}
		System.out.println(currentSitesInRegion.size()+" current sites in CS Region");
		
		int numGrid20Sites = 0;
		int numGrid10Sites = 0;
		int numGrid05Sites = 0;
		int numBasinProfileSites = 0;
		
		List<CybershakeSite> grid05Sites = Lists.newArrayList();
		
		for (CybershakeSite site : sitesInRegion) {
			if (site.type_id == CybershakeSite.TYPE_GRID_20_KM)
				numGrid20Sites++;
			else if (site.type_id == CybershakeSite.TYPE_GRID_10_KM)
				numGrid10Sites++;
			else if (site.type_id == CybershakeSite.TYPE_GRID_05_KM) {
				numGrid05Sites++;
				grid05Sites.add(site);
			} else if (site.name.startsWith("p") && Character.isDigit(site.name.charAt(1)))
				numBasinProfileSites++;
		}
		
		LocationList grid05SiteLocs = new LocationList();
		for (CybershakeSite site : grid05Sites)
			grid05SiteLocs.add(site.createLocation());
		WillsMap2006 wills2006 = new WillsMap2006();
		CVM4i26BasinDepth cvm = new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
		SRTM30PlusTopoSlope topoSlope = new SRTM30PlusTopoSlope();
		
		System.out.println(numGrid20Sites+" 20km sites in Region");
		System.out.println(numGrid10Sites+" 10km sites in Region");
		System.out.println(numGrid05Sites+" 05km sites in Region");
		System.out.println(numBasinProfileSites+" basin profile sites in Region");
		
		List<Double> grid05Vs30 = wills2006.getValues(grid05SiteLocs);
		List<Double> grid05BasinDepth = cvm.getValues(grid05SiteLocs);
		List<Double> grid05TopoSlopes = cvm.getValues(grid05SiteLocs);
		
		double vsCutoff = 300;
		double depthCutoff = 2.5;
		
		int numGrid05VS = 0;
		int numGrid05Depth = 0;
		int numGrid05Both = 0;
		int numTopoDepthCombo = 0;
		
		for (int i=0; i<grid05SiteLocs.size(); i++) {
			if (grid05Vs30.get(i) <= vsCutoff)
				numGrid05VS++;
			if (grid05BasinDepth.get(i) >= depthCutoff)
				numGrid05Depth++;
			if (grid05Vs30.get(i) <= vsCutoff && grid05BasinDepth.get(i) >= depthCutoff)
				numGrid05Both++;
		}
		
		System.out.println(numGrid05VS+" 05km sites with Vs30 <= "+(int)vsCutoff);
		System.out.println(numGrid05Depth+" 05km sites with Z2.5 >= "+(float)depthCutoff);
		System.out.println(numGrid05Both+" 05km sites with both cutoffs");
		
//		LocationList laBorder = new LocationList();
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		laBorder.add(new Location(34.0313, -118.5322));
//		Region laBsinReg = new Region(laBorder, null);
		
		List<String> additions = Lists.newArrayList();
		
		// BASINS
		additions.add("s346");
		additions.add("s366");
		additions.add("s388");
		additions.add("s453");
		additions.add("s493");
		additions.add("s344");
		additions.add("s365");
		additions.add("s328");
		additions.add("s348");
		additions.add("s292");
		additions.add("s071");
		additions.add("s048");
		additions.add("s001");
		additions.add("s003");
		additions.add("s043");
		additions.add("s647");
		additions.add("s668");
		additions.add("s689");
		additions.add("s710");
		additions.add("s731");
		additions.add("s624");
		additions.add("s666");
		additions.add("s541");
		additions.add("s765");
		additions.add("s660");
		additions.add("s794");
		additions.add("s795");
		additions.add("s451");
		additions.add("s491");
		additions.add("s531");
		additions.add("s410");
		
		// SAF
		additions.add("s586");
		additions.add("s588");
		additions.add("s545");
		additions.add("s547");
		additions.add("s505");
		additions.add("s507");
		additions.add("s465");
		additions.add("s467");
		additions.add("s424");
		additions.add("s378");
		additions.add("s380");
		additions.add("s339");
		additions.add("s302");
		additions.add("s266");
		additions.add("s228");
		additions.add("s187");
		additions.add("s145");
		additions.add("s081");
		additions.add("s035");
		
		// check for duplicates
		HashSet<String> runningSet = new HashSet<String>();
		for (String addition : additions) {
			Preconditions.checkState(!runningSet.contains(addition), "duplicate addition: "+addition);
			runningSet.add(addition);
		}
		
		System.out.println(additions.size()+" hand picked additions");
		
		File baseMapFile = new File("/home/kevin/CyberShake/cache/ar_curves_5_21_false_4.0E-4_21.txt");
		ArbDiscrGeoDataSet baseMap = ArbDiscrGeoDataSet.loadXYZFile(baseMapFile.getAbsolutePath(), true);
		GMT_Map map = new GMT_Map(csReg, baseMap, 0.005,
				CyberShake_GMT_MapGenerator.getHazardCPT().rescale(0d, 1d));
		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setDpi(300);
		map.setCustomLabel("GMPE Basemap, 3sec SA");
		
		// add current sites as white
		CPT symSetCPT = new CPT();
		symSetCPT.add(new CPTVal(0f, Color.WHITE, 1f, Color.WHITE));
		symSetCPT.add(new CPTVal(1f, Color.BLACK, 2f, Color.BLACK));
		ArrayList<PSXYSymbol> symbols = Lists.newArrayList();
		ArrayList<Double> symVals = Lists.newArrayList();
		float symSize = 0.1f;
		for (CybershakeSite site : currentSitesInRegion) {
			Point2D pt = new Point2D.Double(site.lon, site.lat);
			symbols.add(new PSXYSymbol(pt, Symbol.INVERTED_TRIANGLE, symSize, 0f, Color.WHITE, Color.WHITE));
			symVals.add(0.5);
		}
		// add new prospective sites as cyan
		for (String addition : additions) {
			CybershakeSite site = sitesByName.get(addition);
			Preconditions.checkNotNull(site, "no site found named "+addition);
			Point2D pt = new Point2D.Double(site.lon, site.lat);
			symbols.add(new PSXYSymbol(pt, Symbol.INVERTED_TRIANGLE, symSize, 0f, Color.CYAN, Color.CYAN));
			symVals.add(1.5);
		}
		PSXYSymbolSet symSet = new PSXYSymbolSet(symSetCPT, symbols, symVals);
		map.setSymbolSet(symSet);
		
		System.out.println("Making map...");
		FaultBasedMapGen.plotMap(new File("/tmp"), "cs_sites", false, map);
		System.out.println("DONE.");
		
		System.out.println(additions.size()+" additional sites:");
		for (String addition : additions) {
			System.out.println(addition);
		}
		
		db.destroy();
	}

}
