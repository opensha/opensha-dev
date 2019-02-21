package scratch.kevin.simulators.ruptures;

import static org.opensha.commons.geo.GeoTools.TO_RAD;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.BSSA_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CB_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CY_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.imr.attenRelImpl.ngaw2.IMT;
import org.opensha.sha.imr.attenRelImpl.ngaw2.Idriss_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_GMM;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.DepthIden;
import org.opensha.sha.simulators.iden.FocalMechIden;
import org.opensha.sha.simulators.iden.LinearRuptureIden;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Range;
import com.google.common.collect.Table;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class BBP_PartBValidationConfig {
	
	// these and values in the get*NGA2 functions before from 
	// https://github.com/SCECcode/bbp/blob/dev/bbp/utils/batch/gmpe_boxplot_gen.py
//	private static double[] BBP_PERIODS = { 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25,
//	                                       0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10 };
	public static double BBP_MAX_ACCEPTANCE_PERIOD = 3d;
	public static double BBP_MIN_ACCEPTANCE_PERIOD = 1d;
	
	/**
	 * If true then the sites should be rJB away from the rupture, otherwise rRup
	 */
	public static final boolean DIST_JB = false;
	
	public enum Scenario {
		M6p6_VERT_SS_SURFACE("M6.6, Vertical Strike-Slip with Surface Rupture", "M6.6 SS", "m6p6_vert_ss_surface",
				new String[] { "M=[6.55,6.65]", "Ztor=[0,1]", "Rake=[-180,-170] or [-10,10] or [170,180]",
						"Dip=90", "Linear rupture (max 0.5km deviation from ideal)"},
				6.6, FaultStyle.STRIKE_SLIP, 90, 0, 5.62, false, true) {
			@Override
			public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException {
				Loader loader = catalog.loader().skipYears(skipYears);
				loader.minMag(6.55).maxMag(6.65);
				loader.matches(new DepthIden(Range.closed(0d, 1d), null));
				loader.matches(FocalMechIden.builder().strikeSlip(10).forDip(90).build());
				loader.matches(new LinearRuptureIden(0.5d));
				try {
					loader.hasTransitions();
				} catch (Exception e) {
					System.out.println("Warning, couldn't force events with transitions. Missing trans file? "+e.getMessage());
				}
				return loader.load();
			}

			@Override
			public Map<String, String> getPublishedComparisonURLs(double distance) {
				Map<String, String> figs = new HashMap<>();
				if (distance == 20d) {
					figs.put("UCSB", getFigureURL(16));
					figs.put("ExSIM", getFigureURL(17));
					figs.put("G&P", getFigureURL(18));
					figs.put("SDSU", getFigureURL(19));
				} else if (distance == 50d) {
					figs.put("UCSB", getFigureURL(20));
					figs.put("ExSIM", getFigureURL(21));
					figs.put("G&P", getFigureURL(22));
					figs.put("SDSU", getFigureURL(23));
				} else {
					return null;
				}
				return figs;
			}
		},
		M6p6_REVERSE("M6.6, Reverse, Dip=45, Ztor=3", "M6.6 Reverse", "m6p6_reverse",
				new String[] { "M=[6.55,6.65]", "Ztor=[1,5]", "Rake=[75,105]", "Dip=[35,55]"},
				6.6, FaultStyle.REVERSE, 45, 3, 15, false, true) {
			@Override
			public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException {
				Loader loader = catalog.loader().skipYears(skipYears);
				loader.minMag(6.55).maxMag(6.65);
				loader.matches(new DepthIden(Range.closed(1d, 5d), null));
				loader.matches(FocalMechIden.builder().forRake(75, 105).forDip(35, 55).build());
				try {
					loader.hasTransitions();
				} catch (Exception e) {
					System.out.println("Warning, couldn't force events with transitions. Missing trans file? "+e.getMessage());
				}
				return loader.load();
			}

			@Override
			public Map<String, String> getPublishedComparisonURLs(double distance) {
				Map<String, String> figs = new HashMap<>();
				if (distance == 20d) {
					figs.put("UCSB", getFigureURL(24));
					figs.put("ExSIM", getFigureURL(25));
					figs.put("G&P", getFigureURL(26));
					figs.put("SDSU", getFigureURL(27));
				} else if (distance == 50d) {
					figs.put("UCSB", getFigureURL(28));
					figs.put("ExSIM", getFigureURL(29));
					figs.put("G&P", getFigureURL(30));
					figs.put("SDSU", getFigureURL(31));
				} else {
					return null;
				}
				return figs;
			}
		},
		M7p2_VERT_SS_SURFACE("M7.2, Vertical Strike-Slip with Surface Rupture", "M7.2 SS", "m7p2_vert_ss_surface",
				new String[] { "M=[7.15,7.25]", "Ztor=[0,1]", "Rake=[-180,-170] or [-10,10] or [170,180]",
						"Dip=90", "Linear rupture (max 0.5km deviation from ideal)"},
				7.2, FaultStyle.STRIKE_SLIP, 90, 0, 12, false, false) {
			@Override
			public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException {
				Loader loader = catalog.loader().skipYears(skipYears);
				loader.minMag(7.15).maxMag(7.25);
				loader.matches(new DepthIden(Range.closed(0d, 1d), null));
				loader.matches(FocalMechIden.builder().strikeSlip(10).forDip(90).build());
				loader.matches(new LinearRuptureIden(0.5d));
				try {
					loader.hasTransitions();
				} catch (Exception e) {
					System.out.println("Warning, couldn't force events with transitions. Missing trans file? "+e.getMessage());
				}
				return loader.load();
			}

			@Override
			public Map<String, String> getPublishedComparisonURLs(double distance) {
				return null;
			}
		};
		
		private String name;
		private String shortName;
		private String prefix;
		private String[] matchCriteria;

		private Table<Double, Double, UncertainArbDiscDataset> rawCriteriaCache;
		private Table<Double, Double, UncertainArbDiscDataset> trimmedCriteriaCache;
		
		private double mag;
		private FaultStyle style;
		private double dip;
		private double zTor;
		private double width;
		private boolean rXpositive;
		
		private boolean officialCriteria;

		private Scenario(String name, String shortName, String prefix, String[] matchCriteria, double mag, FaultStyle style,
				double dip, double zTor, double width, boolean rXpositive, boolean officialCriteria) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.matchCriteria = matchCriteria;
			this.mag = mag;
			this.style = style;
			this.dip = dip;
			this.zTor = zTor;
			this.width = width;
			this.rXpositive = rXpositive;
			this.officialCriteria = officialCriteria;
		}
		
		public double getMagnitude() {
			return mag;
		}
		
		public boolean isOfficialCriteria() {
			return officialCriteria;
		}
		
		public abstract List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears) throws IOException;
		
		public abstract Map<String, String> getPublishedComparisonURLs(double distance);
		
		private static String getFigureURL(int figNum) {
			return "http://www.seismosoc.org/Publications/SRL/SRL_86/srl_86-1_dreger_et_al-esupp/SRL_2014118_esupp_Figure_S"+figNum+".png";
		}
		
		private synchronized UncertainArbDiscDataset calcLoacRawCriterion(double vs30, double distance) {
			if (rawCriteriaCache == null)
				rawCriteriaCache = HashBasedTable.create();
			UncertainArbDiscDataset criterion = rawCriteriaCache.get(vs30, distance);
			if (criterion == null) {
				double rRup, rJB;
				if (DIST_JB) {
					rJB = distance;
					rRup = Math.sqrt(zTor*zTor + rJB*rJB);
				} else {
					rRup = distance;
					rJB = Math.sqrt(rRup*rRup - zTor*zTor);
				}
				double rX = rXpositive ? rJB : -rJB;
				criterion = calcNGA2_Criterion(mag, rRup, rJB, rX, style, dip, zTor, width, vs30);
				criterion.setName("NGA-W2 Mean Prediction");
				rawCriteriaCache.put(vs30, distance, criterion);
			}
			return criterion;
		}
		
		public synchronized UncertainArbDiscDataset getAcceptanceCriteria(double vs30, double distance) {
			UncertainArbDiscDataset criterion = trimmedCriteriaCache.get(vs30, distance);
			if (criterion == null) {
				UncertainArbDiscDataset rawCriterion = calcLoacRawCriterion(vs30, distance);
				
				DiscretizedFunc avgFunc = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
				
				for (int p=0; p<rawCriterion.size(); p++) {
					double period = rawCriterion.getX(p);
					if ((float)p > (float)BBP_MAX_ACCEPTANCE_PERIOD)
						break;
					avgFunc.set(period, rawCriterion.getY(p));
					lowerFunc.set(period, rawCriterion.getLowerY(p));
					upperFunc.set(period, rawCriterion.getUpperY(p));
				}
				criterion = new UncertainArbDiscDataset(avgFunc, lowerFunc, upperFunc);
				criterion.setName("NGA-W2 Acceptance Criteria");
				trimmedCriteriaCache.put(vs30, distance, criterion);
			}
			
			return criterion;
		}
		
		public DiscretizedFunc getMeanPrediction(double vs30, double distance) {
			return calcLoacRawCriterion(vs30, distance);
		}

		public String getName() {
			return name;
		}

		public String getShortName() {
			return shortName;
		}

		public String getPrefix() {
			return prefix;
		}

		public String[] getMatchCriteria() {
			return matchCriteria;
		}
	}
	
	public static double[] OFFICIAL_DISTANCES = { 20d, 50d };
	public static Scenario[] OFFICIAL_SCENARIOS = {Scenario.M6p6_VERT_SS_SURFACE, Scenario.M6p6_REVERSE};
	
	public static List<RSQSimEvent> getBestMatches(double targetMag, List<RSQSimEvent> matches, int maxNum) {
		if (matches.size() <= maxNum)
			return matches;
		matches = new ArrayList<>(matches);
		Collections.sort(matches, new MagDiffEventComparator(targetMag));
		matches = matches.subList(0, maxNum);
		Collections.sort(matches);
		return matches;
	}
	
	private static class MagDiffEventComparator implements Comparator<RSQSimEvent> {
		
		private double targetMag;

		public MagDiffEventComparator(double targetMag) {
			this.targetMag = targetMag;
		}

		@Override
		public int compare(RSQSimEvent e0, RSQSimEvent e1) {
			double diff1 = Math.abs(e0.getMagnitude() - targetMag);
			double diff2 = Math.abs(e1.getMagnitude() - targetMag);
			return Double.compare(diff1, diff2);
		}
		
	}
	
	private static UncertainArbDiscDataset calcNGA2_Criterion(double mag, double rRup, double rJB, double rX, FaultStyle style, double dip,
			double zTor, double width, double vs30) {
		NGAW2_GMM[] gmms = { new ASK_2014(), new BSSA_2014(), new CB_2014(), new CY_2014(), new Idriss_2014() };
		
//		double zHyp = zTor + Math.sin(dip * TO_RAD) * width / 2.0;
		double zHyp = 9.498;
//		System.out.println("zHyp: "+zHyp);
		
		HashSet<IMT> imtSet = null;
		
		for (NGAW2_GMM gmm : gmms) {
			gmm.set_dip(dip);
			gmm.set_fault(style);
			gmm.set_Mw(mag);
			gmm.set_rJB(rJB);
			gmm.set_rRup(rRup);
			gmm.set_rX(rX);
			gmm.set_vs30(vs30);
			gmm.set_vsInf(false);
			gmm.set_width(width);
			gmm.set_z1p0(Double.NaN);
			gmm.set_z2p5(Double.NaN);
			gmm.set_zHyp(zHyp);
			gmm.set_zTop(zTor);
			
			Collection<IMT> myIMTs = gmm.getSupportedIMTs();
			if (imtSet == null)
				imtSet = new HashSet<>(myIMTs);
			else
				imtSet.retainAll(myIMTs);
		}
		
		DiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc();
		
		for (IMT imt : imtSet) {
			Double period = imt.getPeriod();
			if (period == null)
				continue;
			double[] means = new double[gmms.length];
			for (int i=0; i<gmms.length; i++) {
				gmms[i].set_IMT(imt);
				ScalarGroundMotion gm = gmms[i].calc();
				means[i] = Math.exp(gm.mean());
//				if (period == 3d)
//					System.out.println(ClassUtils.getClassNameWithoutPackage(gmms[i].getClass())+": "+means[i]);
				Preconditions.checkState(Double.isFinite(means[i]), "Bad value for imt %s, gmm %s: %s", imt, gmms[i], means[i]);
			}
//			if (period == 3d)
//				System.exit(0);
			
			double mean = StatUtils.mean(means);
			
			double min1 = StatUtils.min(means)*0.85;
			double max1 = StatUtils.max(means)*1.15;
			double min2 = mean*Math.log(2);
			double max2 = mean/Math.log(2);
			
			double min = Math.min(min1, min2);
			double max = Math.max(max1, max2);
			
			minFunc.set(period, min);
			meanFunc.set(period, mean);
			maxFunc.set(period, max);
		}
		return new UncertainArbDiscDataset(meanFunc, minFunc, maxFunc);
	}
	
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
					double dist;
					if (DIST_JB)
						dist = LocationUtils.horzDistanceFast(loc, startLoc);
					else
						dist = LocationUtils.linearDistanceFast(loc, startLoc);
					if (dist < minDist) {
						minDist = dist;
						closestLoc = loc;
					}
				}
			}
			
			LocationVector vector = LocationUtils.vector(closestLoc, startLoc);
			if (DIST_JB) {
				vector.setHorzDistance(distance);
			} else {
				double vertDist = vector.getVertDistance();
				// change horzonatal distance such that 3-D dist matches target
				double horzDist = Math.sqrt(distance*distance - vertDist*vertDist);
				vector.setHorzDistance(horzDist);
			}
			
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
		System.out.println(calcNGA2_Criterion(6.6, 20, 20, -20, FaultStyle.STRIKE_SLIP, 90, 0, 5.62, 500));
//		System.out.println(getNGA2(6.6, 20, FaultStyle.STRIKE_SLIP, 90, 0, 863));
		System.exit(0);
		
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
		int skipYears = 2000;
		
		int numToPlot = 20;
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
				
				String prefix = "match_"+idStr+"_event_"+event.getID()+"_m"+(float)event.getMagnitude();
				
				plotEventAndSites(catalog, event, distColors, numSites, randomAz, scenarioDir, prefix);
			}
		}
	}

	static void plotEventAndSites(RSQSimCatalog catalog, RSQSimEvent event, Color[] distColors, int numSites,
			boolean randomAz, File outputDir, String prefix) throws IOException {
		double locRectWidth = 0.01;
		
		List<XYAnnotation> anns = new ArrayList<>();
		for (int d=0; d<OFFICIAL_DISTANCES.length; d++) {
			double distance = OFFICIAL_DISTANCES[d];
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
		
		RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, null, outputDir, prefix,
				null, null, null, null, null, null, anns);
	}

}