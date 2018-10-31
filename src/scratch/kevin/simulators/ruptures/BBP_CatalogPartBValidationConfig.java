package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
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
	
	// these and values in the get*NGA2 functions before from 
	// https://github.com/SCECcode/bbp/blob/dev/bbp/utils/batch/gmpe_boxplot_gen.py
	private static double[] bbp_periods = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25,
	                                       0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10 };
	private static double bbp_max_acceptance_period = 3d;
	
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

			@Override
			double[] getMeanNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.123231246, 0.125469445, 0.136076373,
                            0.175016514, 0.222432391, 0.248653998,
                            0.271641827, 0.25852406, 0.232697091,
                            0.207732879, 0.169808845, 0.140483291,
                            0.092461139, 0.06701949, 0.04059643,
                            0.02825896, 0.017440236, 0.012132,
                            0.00921392, 0.005018336, 0.003001676 };
				if (distance == 50d)
					return new double[] { 0.044451045, 0.045187755, 0.049722252,
                            0.063144032, 0.079112127, 0.088121596,
                            0.096865495, 0.093803595, 0.086597108,
                            0.079054594, 0.06593417, 0.055200516,
                            0.037401207, 0.027011192, 0.016746435,
                            0.011705015, 0.007208154, 0.00512826,
                            0.003914806, 0.002287936, 0.001397211 };
				throw new IllegalStateException("Unsupported distance: "+distance);
			}

			@Override
			double[] getUpperNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.179247367, 0.18480203, 0.215915347, 0.292453367,
		                     0.356440194, 0.375056915, 0.391896317, 0.372971379,
		                     0.335710939, 0.299695194, 0.244982378, 0.202674547,
		                     0.133393226, 0.096688686, 0.058568268, 0.040769061,
		                     0.025160942, 0.01757914, 0.014688678, 0.008873136,
		                     0.005018356 };
				if (distance == 50d)
					return new double[] { 0.06661183, 0.068101077, 0.078444572, 0.103073492,
		                     0.122849466, 0.129733734, 0.13974737, 0.135329981,
		                     0.124933219, 0.11405167, 0.0951229, 0.079637511,
		                     0.053958536, 0.038968912, 0.024159999, 0.016886766,
		                     0.010399168, 0.007978404, 0.006734335, 0.004488779,
		                     0.002598649 };
				throw new IllegalStateException("Unsupported distance: "+distance);
			}

			@Override
			double[] getLowerNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.085417391, 0.086968792, 0.094320954, 0.117282459,
		                     0.152663489, 0.172353817, 0.188287767, 0.179195223,
		                     0.161293332, 0.143989459, 0.117702522, 0.097375597,
		                     0.064089178, 0.046454371, 0.028139301, 0.019587618,
		                     0.011989636, 0.007364112, 0.004795379, 0.002143994,
		                     0.001177935 };
				if (distance == 50d)
					return new double[] { 0.030572971, 0.031321765, 0.034464839, 0.043768108,
		                     0.054836348, 0.061081236, 0.067142045, 0.065019697,
		                     0.060024542, 0.054796469, 0.045702084, 0.038262082,
		                     0.025924541, 0.018722731, 0.011607744, 0.008113298,
		                     0.004697369, 0.002904508, 0.001895507, 0.000848924,
		                     0.000466798 };
				throw new IllegalStateException("Unsupported distance: "+distance);
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

			@Override
			double[] getMeanNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.132602154, 0.135177304, 0.149903826,
                            0.191287243, 0.244288563, 0.275286188,
                            0.302126379, 0.28824219, 0.260950763,
                            0.234128463, 0.190300147, 0.156459605,
                            0.102110994, 0.073483962, 0.0439006,
                            0.029943494, 0.017940568, 0.01204706,
                            0.008983379, 0.004841892, 0.002874356 };
				if (distance == 50d)
					return new double[] { 0.048775377, 0.049552538, 0.054407432,
                            0.06857154, 0.086312568, 0.096842118,
                            0.106774396, 0.103615204, 0.096257186,
                            0.088352351, 0.073281556, 0.060968573,
                            0.040944855, 0.029355217, 0.017953993,
                            0.012300317, 0.007364016, 0.005059772,
                            0.003793167, 0.002196334, 0.001331835 };
				throw new IllegalStateException("Unsupported distance: "+distance);
			}

			@Override
			double[] getUpperNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.194889368, 0.200576566, 0.234480754, 0.312298198,
		                     0.382137301, 0.406009288, 0.435876229, 0.415845578,
		                     0.376472372, 0.337775973, 0.274545078, 0.225723496,
		                     0.147315025, 0.106014948, 0.063335178, 0.04319933,
		                     0.025882768, 0.01757914, 0.014688678, 0.008873136,
		                     0.005018356 };
				if (distance == 50d)
					return new double[] { 0.072425493, 0.073914951, 0.085211981, 0.110196595,
		                     0.131902644, 0.140616093, 0.154042892, 0.149485142,
		                     0.138869765, 0.127465499, 0.105722937, 0.087959058,
		                     0.059070939, 0.042350627, 0.025902137, 0.017745606,
		                     0.010624029, 0.007978404, 0.006734335, 0.004488779,
		                     0.002598649 };
				throw new IllegalStateException("Unsupported distance: "+distance);
			}

			@Override
			double[] getLowerNGA2(double distance) {
				if (distance == 20d)
					return new double[] { 0.091912809, 0.093697767, 0.103905415, 0.132085813,
		                     0.169327928, 0.190813845, 0.209418048, 0.199794261,
		                     0.180877286, 0.162285484, 0.13190601, 0.108449534,
		                     0.070777948, 0.050935201, 0.030429577, 0.020755248,
		                     0.012435454, 0.007913061, 0.005022391, 0.002186587,
		                     0.001190133 };
				if (distance == 50d)
					return new double[] { 0.033119304, 0.034025707, 0.037341911, 0.047530169,
		                     0.059827313, 0.067125841, 0.074010372, 0.071772064,
		                     0.066391977, 0.061241183, 0.050794904, 0.042260195,
		                     0.028380811, 0.020347486, 0.01244476, 0.00852593,
		                     0.005104347, 0.003121034, 0.00198524, 0.000865789,
		                     0.000471632 };
				throw new IllegalStateException("Unsupported distance: "+distance);
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
		
		abstract double[] getMeanNGA2(double distance);
		abstract double[] getUpperNGA2(double distance);
		abstract double[] getLowerNGA2(double distance);
		
		public UncertainArbDiscDataset getAcceptanceCriteria(double distance) {
			double[] avgVals = getMeanNGA2(distance);
			double[] lowerVals = getLowerNGA2(distance);
			double[] upperVals = getUpperNGA2(distance);
			Preconditions.checkState(bbp_periods.length == avgVals.length);
			Preconditions.checkState(bbp_periods.length >= lowerVals.length);
			Preconditions.checkState(bbp_periods.length >= upperVals.length);
			
			DiscretizedFunc avgFunc = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
			for (int p=0; p<avgVals.length && bbp_periods[p] <= bbp_max_acceptance_period; p++) {
				avgFunc.set(bbp_periods[p], avgVals[p]);
				lowerFunc.set(bbp_periods[p], lowerVals[p]);
				upperFunc.set(bbp_periods[p], upperVals[p]);
			}
			UncertainArbDiscDataset func = new UncertainArbDiscDataset(avgFunc, lowerFunc, upperFunc);
			func.setName("NGA-W2 Acceptance Criteria");
			return func;
		}
		
		public DiscretizedFunc getMeanPrediction(double distance) {
			double[] avgVals = getMeanNGA2(distance);
			Preconditions.checkState(bbp_periods.length == avgVals.length);
			DiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
			func.setName("NGA-W2 Mean Prediction");
			for (int p=0; p<avgVals.length; p++)
				func.set(bbp_periods[p], avgVals[p]);
			return func;
		}

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
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
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
