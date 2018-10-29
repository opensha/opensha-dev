package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.DepthIden;
import org.opensha.sha.simulators.iden.FocalMechIden;
import org.opensha.sha.simulators.iden.LinearRuptureIden;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class BBP_CatalogPartBValidationConfig {
	
	public enum Scenario {
		M6p6_VERT_SS_SURFACE("M6.6, vertical strike slip with surface rupture", "m6p6_vert_ss_surface",
				new String[] { "M=[6.55,6.65]", "Ztor=[0,1]", "Rake=[-180,-170] or [-10,10] or [170,180]",
						"Dip=90", "Linear rupture (max 0.5km deviation from ideal)"}) {
			@Override
			public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException {
				Loader loader = catalog.loader().skipYears(skipYears);
				loader.minMag(6.55).maxMag(6.65);
				loader.matches(new DepthIden(Range.closed(0d, 1d), null));
				loader.matches(FocalMechIden.builder().strikeSlip(10).forDip(90).build());
				loader.matches(new LinearRuptureIden(0.5d));
				return loader.load();
			}
		},
		M6p6_REVERSE("M6.6, Reverse, Dip=45, Ztor=3", "m6p6_reverse",
				new String[] { "M=[6.55,6.65]", "Ztor=[2,4]", "Rake=[80,100]", "Dip=[40,50]"}) {
			@Override
			public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException {
				Loader loader = catalog.loader().skipYears(skipYears);
				loader.minMag(6.55).maxMag(6.65);
				loader.matches(new DepthIden(Range.closed(1d, 5d), null));
				loader.matches(FocalMechIden.builder().forRake(75, 105).forDip(35, 55).build());
				return loader.load();
			}
		};
		
		private String name;
		private String prefix;
		private String[] matchCriteria;

		private Scenario(String name, String prefix, String[] matchCriteria) {
			this.name = name;
			this.prefix = prefix;
			this.matchCriteria = matchCriteria;
		}
		
		public abstract List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException;

		public String getName() {
			return name;
		}

		public String getPrefix() {
			return prefix;
		}

		public String[] getMatchCriteria() {
			return matchCriteria;
		}
	}
	
	public static double[] DISTANCES = { 20d, 50d };
	
	public static Location[] selectSitesSites(int num, double distance, boolean randomAz, RSQSimCatalog catalog, RSQSimEvent event) {
		// start with GMPE surface in order to determine footwall
		RuptureSurface rupSurf = catalog.getGMPE_Rupture(
				event, RSQSimBBP_Config.MIN_SUB_SECT_FRACT).getRuptureSurface();
		Location firstLoc = rupSurf.getFirstLocOnUpperEdge();
		Location lastLoc = rupSurf.getLastLocOnUpperEdge();
		// footwall will be on the left of the line from first to last
		double strike = LocationUtils.vector(firstLoc, lastLoc).getAzimuthRad();
//		System.out.println("First Loc: "+firstLoc);
//		System.out.println("Last Loc: "+lastLoc);
//		System.out.println("Strike: "+(float)strike+" ("+(float)Math.toDegrees(strike)+")");
		
		List<SimulatorElement> elems = event.getAllElements();
		
		double minDepth = Double.POSITIVE_INFINITY;
		
		double aveLat = 0d;
		double aveLon = 0d;
		double aveDep = 0d;
		
		for (SimulatorElement elem : elems) {
			Location loc = elem.getCenterLocation();
			aveLat += loc.getLatitude();
			aveLon += loc.getLongitude();
			aveDep += loc.getDepth();
			for (Location l : elem.getVertices())
				minDepth = Math.min(minDepth, l.getDepth());
		}
		
		aveLat /= elems.size();
		aveLon /= elems.size();
		aveDep /= elems.size();
		Preconditions.checkState(aveDep >= minDepth);
		
		Location centerLoc = new Location(aveLat, aveLon, aveDep);
		
		// now we want to move up the fault
		
		// find closest trace loc
		Location closestTraceLoc = null;
		double traceDist = Double.POSITIVE_INFINITY;
		for (Location traceLoc : rupSurf.getEvenlyDiscritizedUpperEdge()) {
			double dist = LocationUtils.linearDistanceFast(centerLoc, traceLoc);
			if (dist < traceDist) {
				traceDist = dist;
				closestTraceLoc = traceLoc;
			}
		}
		if (traceDist > 0) {
			LocationVector vectorToTrace = LocationUtils.vector(centerLoc, closestTraceLoc);
			double depthRatio = (aveDep - minDepth)/(aveDep - closestTraceLoc.getDepth());
			vectorToTrace.setVertDistance(vectorToTrace.getVertDistance()*depthRatio);
			vectorToTrace.setHorzDistance(vectorToTrace.getHorzDistance()*depthRatio);
			
			centerLoc = LocationUtils.location(centerLoc, vectorToTrace);
		}
		centerLoc = new Location(centerLoc.getLatitude(), centerLoc.getLongitude());
		
		Location[] ret = new Location[num];
		
		for (int i=0; i<num; i++) {
			// to the left will be in the range [azimuth+PI,azimuth+2*PI]
			double azimuth;
			if (randomAz)
				azimuth = strike + Math.PI + Math.random()*Math.PI;
			else
				azimuth = strike + Math.PI + ((double)i/(double)(num-1))*Math.PI;
			
			// go distance in random location from center
			Location startLoc = LocationUtils.location(centerLoc, azimuth, distance);
			
			// that loc will likely be closer than intended though, so adjust to actual distances
			Location closestLoc = null;
			double minDist = Double.POSITIVE_INFINITY;
			
			for (SimulatorElement elem : elems) {
				for (Location loc : elem.getVertices()) {
					double dist = LocationUtils.linearDistanceFast(loc, startLoc);
					if (dist < minDist) {
						minDist = dist;
						closestLoc = loc;
					}
				}
			}
			
			LocationVector vector = LocationUtils.vector(closestLoc, startLoc);
			double vertDist = vector.getVertDistance();
			// change horzonatal distance such that 3-D dist matches target
			double horzDist = Math.sqrt(distance*distance - vertDist*vertDist);
			vector.setHorzDistance(horzDist);
			
			Location loc = LocationUtils.location(closestLoc, vector);
			if (loc.getDepth() != 0d) {
				Preconditions.checkState(Math.abs(loc.getDepth()) < 0.01,
						"Bad site depth!\n\tSite loc: %s\n\tClosest loc:%s\n\tVector:%s", loc, closestLoc, vector);
				loc = new Location(loc.getLatitude(), loc.getLongitude());
			}
			double calcDist = LocationUtils.linearDistanceFast(loc, closestLoc);
			Preconditions.checkState(Math.abs(calcDist - distance) < 0.01,
					"Bad site distance: %s != %s\n\tSite loc: %s\n\tClosest loc:%s\n\tVector:%s",
					(float)calcDist, (float)distance, loc, closestLoc, vector);
			ret[i] = loc;
		}
		
		return ret;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		int skipYears = 2000;
		
		int numToPlot = 20;
		double locRectWidth = 0.01;
		Color[] distColors = {Color.BLUE.darker(), Color.GREEN.darker()};
		int numSites = 100;
		boolean randomAz = false;
		
		File outputDir = new File(catalog.getCatalogDir(), "bbp_part_b");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (Scenario scenario : Scenario.values()) {
			System.out.println("Scenario: "+scenario.getName());
			
			List<RSQSimEvent> events = scenario.getMatches(catalog, skipYears);
			System.out.println("Found "+events.size()+" matches!");
			
			File scenarioDir = new File(outputDir, scenario.getPrefix());
			Preconditions.checkState(scenarioDir.exists() || scenarioDir.mkdir());
			
			for (int i=0; i<events.size() && i<numToPlot; i++) {
				String idStr = i+"";
				while (idStr.length() < ((numToPlot-1)+"").length())
					idStr = "0"+idStr;
				
				RSQSimEvent event = events.get(i);
				
				List<XYAnnotation> anns = new ArrayList<>();
				for (int d=0; d<DISTANCES.length; d++) {
					double distance = DISTANCES[d];
					Color c = distColors[d];
					
					for (Location loc : selectSitesSites(numSites, distance, randomAz, catalog, event)) {
						double[] poly = new double[10];
						double lat = loc.getLatitude();
						double lon = loc.getLongitude();
						double ux = lon+0.5*locRectWidth;
						double lx = lon-0.5*locRectWidth;
						double uy = lat+0.5*locRectWidth;
						double ly = lat-0.5*locRectWidth;
						poly[0] = ux;
						poly[1] = uy;
						poly[2] = lx;
						poly[3] = uy;
						poly[4] = lx;
						poly[5] = ly;
						poly[6] = ux;
						poly[7] = ly;
						poly[8] = ux;
						poly[9] = uy;
						XYPolygonAnnotation ann = new XYPolygonAnnotation(poly, null, null, c);
						anns.add(ann);
					}
				}
				
				String prefix = "match_"+idStr+"_event_"+event.getID()+"_m"+(float)event.getMagnitude();
				RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, null, scenarioDir, prefix,
						null, null, null, null, null, null, anns);
			}
		}
	}

}
