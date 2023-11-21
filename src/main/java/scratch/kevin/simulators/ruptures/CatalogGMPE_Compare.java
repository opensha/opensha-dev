package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.IntensityMeasureRelationship;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.FocalMechIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simCompare.MultiRupGMPE_ComparePageGen;
import scratch.kevin.simCompare.RuptureComparison;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import org.opensha.commons.util.MarkdownUtils;

class CatalogGMPE_Compare extends MultiRupGMPE_ComparePageGen<RSQSimEvent> {
	
	private RSQSimCatalog catalog;
	private BBP_CatalogSimZipLoader bbpZipFile;
	private List<BBP_Site> sites;
	private List<Site> highlightSites;
	private double minMag;

	private BiMap<BBP_Site, Site> sitesBBPtoGMPE;
	private BiMap<Site, BBP_Site> sitesGMPEtoBBP;
	private List<RSQSimEvent> events;
	private double catDurationYears;
	private Map<Integer, RSQSimEvent> eventsMap;
	private Map<Integer, EqkRupture> gmpeRupsMap;
	
	private boolean rupGen = false;
	
	private File gmpeCacheDir;

	public CatalogGMPE_Compare(RSQSimCatalog catalog, ZipFile bbpZipFile, List<BBP_Site> sites, double minMag, int skipYears,
			File gmpeCacheDir, RuptureIdentifier loadCriteria, VelocityModel vm, double maxDist) throws IOException {
		this.catalog = catalog;
		this.sites = sites;
		this.minMag = minMag;
		this.gmpeCacheDir = gmpeCacheDir;
		
		sitesBBPtoGMPE = HashBiMap.create();
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site site : sites) {
			Site gmpeSite = site.buildGMPE_Site(vm);
			gmpeSite.setName(RSQSimBBP_Config.siteCleanName(site));
			gmpeSites.add(gmpeSite);
			sitesBBPtoGMPE.put(site, gmpeSite);
		}
		sitesGMPEtoBBP = sitesBBPtoGMPE.inverse();
		
		Loader loader = catalog.loader();
		if (minMag > 0d)
			loader.minMag(minMag);
		if (skipYears > 0)
			loader.skipYears(skipYears);
		if (loadCriteria != null)
			loader.matches(loadCriteria);
		if (sites.size() < 20) {
			ArrayList<RegionIden> siteRegIdens = new ArrayList<>();
			for (BBP_Site site : sites)
				siteRegIdens.add(new RegionIden(new Region(site.getLoc(), maxDist)));
			loader.matches(new LogicalOrRupIden(siteRegIdens));
		}
		System.out.println("Loading events...");
		events = loader.minMag(minMag).load();
		System.out.println("Loaded "+events.size()+" events");
		eventsMap = new HashMap<>();
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		
		this.bbpZipFile = new BBP_CatalogSimZipLoader(bbpZipFile, sites, sitesBBPtoGMPE, eventsMap);
		
		// need to explicitly calculate duration from events used, as there could be more without transitions
		double startTime = Double.POSITIVE_INFINITY;
		double endTime = Double.NEGATIVE_INFINITY;
		for (BBP_Site site : sites) {
			for (Integer eventID : this.bbpZipFile.getEventIDs(site)) {
				double t = eventsMap.get(eventID).getTimeInYears();
				startTime = Math.min(startTime, t);
				endTime = Math.max(endTime, t);
			}
		}
		catDurationYears = endTime - startTime;
		
		gmpeRupsMap = new HashMap<>();
		
		double maxCatalogMag = 0d;
		for (RSQSimEvent e : events)
			maxCatalogMag = Math.max(maxCatalogMag, e.getMagnitude());
		
		init(this.bbpZipFile, gmpeSites, DIST_JB, RSQSimBBP_Config.MAX_DIST, minMag, maxCatalogMag);
	}
	
	public BBP_CatalogSimZipLoader getSimProv() {
		return bbpZipFile;
	}
	
	public List<RSQSimEvent> getEvents() {
		return events;
	}
	
	public Collection<Site> getGMPESites() {
		return sitesGMPEtoBBP.keySet();
	}
	
	public void setDoRupGen(boolean rupGen, File gmpeCacheDir) {
		this.rupGen = rupGen;
		this.gmpeCacheDir = gmpeCacheDir;
		gmpeRupsMap.clear();
	}
	
	public void setHighlightSites(String... highlightNames) {
		if (highlightNames == null || highlightNames.length == 0) {
			// no highlights. null would plot all sites, empty list is no highlights
			this.highlightSites = new ArrayList<>();
			return;
		}
		List<Site> highlightSites = new ArrayList<>();
		for (String name : highlightNames) {
			for (BBP_Site site : sites) {
				if (name.startsWith(site.getName())) {
					highlightSites.add(sitesBBPtoGMPE.get(site));
					break;
				}
			}
		}
		setHighlightSites(highlightSites);
	}
	
	public void setHighlightSites(List<Site> highlightSites) {
		this.highlightSites = highlightSites;
	}
	
	private EqkRupture getGMPE_Rup(RSQSimEvent event) {
		synchronized (gmpeRupsMap) {
			EqkRupture rup = gmpeRupsMap.get(event.getID());
			if (rup == null) {
				if (rupGen) {
					BBP_PlanarSurface plane = null;
					if (RSQSimBBP_Config.U3_SURFACES)
						plane = RSQSimBBP_Config.planarEquivalentU3Surface(catalog, event);
					else
						plane = RSQSimBBP_Config.estimateBBP_PlanarSurface(event);
					BBP_SourceFile source = RSQSimBBP_Config.buildBBP_Source(event, plane, RSQSimBBP_Config.DEFAULT_SEED);
					rup = new EqkRupture(event.getMagnitude(), source.getFocalMechanism().getRake(),
							plane.getQuadSurface(), source.getHypoLoc());
					rup.setRuptureSurface(plane.getQuadSurface());
				} else {
					rup = catalog.getEqkRupture(event);
				}
				gmpeRupsMap.put(event.getID(), rup);
			}
			return rup;
		}
	}
	
	private int loadCache(File file, Site site, Map<Integer, EventComparison> eventComps) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		
		int count = 0;
		
		for (int row=1; row<csv.getNumRows(); row++) {
			Integer eventID = Integer.parseInt(csv.get(row, 0));
			String imtStr = csv.get(row, 1);
			IMT imt = IMT.forString(imtStr);
			Double logMean = Double.parseDouble(csv.get(row, 2));
			Double stdDev = Double.parseDouble(csv.get(row, 3));
			Double rRup = Double.parseDouble(csv.get(row, 4));
			Double rJB = Double.parseDouble(csv.get(row, 5));
			
			EventComparison comp = eventComps.get(eventID);
			if (comp != null) {
				comp.addResult(site, imt, logMean, stdDev);
				comp.setDistances(site, rRup, rJB);
				count++;
			}
		}
		
		return count;
	}
	
	private void writeCache(File file, Site site, Map<Integer, EventComparison> comps) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Event ID", "Period (s)", "Ln(Mean)", "Std Dev", "rRup (km)", "rJB (km)");
		
		for (Integer eventID : sorted(comps.keySet())) {
			EventComparison comp = comps.get(eventID);
			Set<IMT> siteIMTs = comp.getIMTs(site);
			if (siteIMTs != null)
				for (IMT imt : sorted(siteIMTs))
					csv.addLine(eventID+"", imt.name(), comp.getLogMean(site, imt)+"", comp.getStdDev(site, imt)+"",
							comp.getDistanceRup(site)+"", comp.getDistanceJB(site)+"");
		}
		
		csv.writeToFile(file);
	}
	
	private <E extends Comparable<E>> Iterable<E> sorted(Set<E> vals) {
		List<E> list = new ArrayList<>(vals);
		Collections.sort(list);
		return list;
	}
	
	public synchronized void calcGMPE(Map<Integer, EventComparison> eventComps, Site site, AttenRelRef gmpeRef, IMT... imts) {
		List<Future<?>> futures = new ArrayList<>();
		
		File gmpeCacheFile = null;
		if (gmpeCacheDir != null) {
			String siteName = sitesGMPEtoBBP.get(site).getName();
			gmpeCacheFile = new File(gmpeCacheDir, siteName+"_"+gmpeRef.getShortName()+".csv");
			System.out.println(gmpeCacheFile.getAbsolutePath()+" exits ? "+gmpeCacheFile.exists());
			if (gmpeCacheFile.exists()) {
				try {
					System.out.println("Loading cache from "+gmpeCacheFile.getAbsolutePath());
					int added = loadCache(gmpeCacheFile, site, eventComps);
					System.out.println("Loaded "+added+" cached values");
				} catch (IOException e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
		}
		
		List<IMT> imtsList = Lists.newArrayList(imts);
		
		for (EventComparison comp : eventComps.values()) {
			if (!comp.isSiteApplicable(site))
				continue;
			List<IMT> myIMTs;
			Set<IMT> computedIMTs = null;
			if (comp.hasSite(site))
				computedIMTs = comp.getIMTs(site);
			if (computedIMTs == null || computedIMTs.isEmpty()) {
				// nothing cached
				myIMTs = imtsList;
			} else {
				myIMTs = new ArrayList<>();
				for (IMT imt : imts)
					if (!computedIMTs.contains(imt))
						myIMTs.add(imt);
			}
			if (!myIMTs.isEmpty())
				futures.add(exec.submit(new EventComparisonCalc(comp, gmpeRef, site, myIMTs)));
		}
		
		System.out.println("Calculating for "+futures.size()+" events, site "+site.getName());
		
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (Exception e) {
				if (e instanceof ExecutionException) {
					e.getCause().printStackTrace();
					System.err.flush();
					throw ExceptionUtils.asRuntimeException(e.getCause());
				}
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		if (gmpeCacheFile != null && !futures.isEmpty()) {
			System.out.println("Caching "+futures.size()+" new vals to "+gmpeCacheFile.getAbsolutePath());
			try {
				writeCache(gmpeCacheFile, site, eventComps);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	private class EventComparisonCalc implements Runnable {
		
		private EventComparison eventComp;
		private AttenRelRef gmpeRef;
		private Site site;
		private List<IMT> imts;
		
		public EventComparisonCalc(EventComparison eventComp, AttenRelRef gmpeRef, Site site, List<IMT> imts) {
			this.eventComp = eventComp;
			this.gmpeRef = gmpeRef;
			this.site = site;
			this.imts = imts;
		}

		@Override
		public void run() {
			ScalarIMR gmpe = checkOutGMPE(gmpeRef);
			
			eventComp.calculate(gmpe, site, imts.toArray(new IMT[0]));
			
			checkInGMPE(gmpeRef, gmpe);
		}
		
	}
	
	protected class EventComparison extends RuptureComparison.Cached<RSQSimEvent> {
		
		private Collection<Site> applicableSites;
		
		public EventComparison(RSQSimEvent event) {
			super(event);
			
			applicableSites = new HashSet<>();
		}
		
		public void addApplicableSite(Site site) {
			applicableSites.add(site);
		}
		
		public boolean isSiteApplicable(Site site) {
			return applicableSites.contains(site);
		}

		@Override
		public Collection<Site> getApplicableSites() {
			return applicableSites;
		}

		@Override
		public EqkRupture getGMPERupture() {
			return getGMPE_Rup(getRupture());
		}

		@Override
		public double getMagnitude() {
			return getRupture().getMagnitude();
		}

		@Override
		public double getAnnualRate() {
			return 1d/catDurationYears;
		}

		@Override
		public double getRuptureTimeYears() {
			return getRupture().getTimeInYears();
		}
	}
	
	public List<EventComparison> loadCalcComps(AttenRelRef gmpeRef, IMT[] imts) {
		Map<Integer, EventComparison> compsMap = new HashMap<>();
		
		for (BBP_Site site : sites) {
			Site gmpeSite = sitesBBPtoGMPE.get(site);
			for (Integer eventID : bbpZipFile.getEventIDs(site)) {
				EventComparison comp = compsMap.get(eventID);
				if (comp == null) {
					comp = new EventComparison(eventsMap.get(eventID));
					compsMap.put(eventID, comp);
				}
				if (eventID == 655368)
					System.out.println("Adding site for "+eventID+" "+site.getName());
				comp.addApplicableSite(gmpeSite);
			}
			calcGMPE(compsMap, sitesBBPtoGMPE.get(site), gmpeRef, imts);
		}
		return new ArrayList<>(compsMap.values());
	}
	
	public void generateGMPE_Page(File outputDir, AttenRelRef gmpeRef, IMT[] imts, List<EventComparison> comps)
			throws IOException {
		LinkedList<String> lines = new LinkedList<>();
		
		String distDescription = getDistDescription();
		
		// header
		if (rupGen) {
			lines.add("# "+catalog.getName()+" BBP Rupture Generator/"+gmpeRef.getShortName()+" GMPE Comparisons");
			lines.add("");
			lines.add("**NOTE: These tests use the BBP Rupture Generator, not RSQSim slip/time histories.**");
			lines.add("");
			lines.add("*All calculations here (both BBP & GMPE) use simple planar surface representations, which can be "
					+ "a poor approximation for many ruptures, but the same surface is used for both calculations.*");
		} else {
			lines.add("# "+catalog.getName()+" BBP/"+gmpeRef.getShortName()+" GMPE Comparisons");
		}
		lines.add("");
		lines.add("**GMPE: "+gmpeRef.getName()+"**");
		lines.add("");
		lines.add("Ruptures are binned by their moment magnitude (**Mw**) and the "+distDescription);
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		
		super.generateGMPE_Page(outputDir, lines, gmpeRef, imts, comps, highlightSites);
	}
	
	@Override
	protected double calcRupAzimuthDiff(RSQSimEvent event, int simIndex, Site site) {
		return calcRupAzimuthDiff(event, getGMPE_Rup(event), site.getLocation());
	}
	
	private static double calcRupAzimuthDiff(RSQSimEvent event, EqkRupture rup, Location loc) {
		Location loc1 = rup.getRuptureSurface().getFirstLocOnUpperEdge();
		Location loc2 = rup.getRuptureSurface().getLastLocOnUpperEdge();
		
		double aveLat = 0d;
		double aveLon = 0d;
		ArrayList<SimulatorElement> elems = event.getAllElements();
		for (SimulatorElement elem : elems) {
			Location l = elem.getCenterLocation();
			aveLat += l.getLatitude();
			aveLon += l.getLongitude();
		}
		aveLat /= elems.size();
		aveLon /= elems.size();
		Location centroid = new Location(aveLat, aveLon);
		
		Location hypo = rup.getHypocenterLocation();
		return calcRupAzimuthDiff(loc1, loc2, centroid, hypo, loc);
	}
	
	private void testPlotRupAzDiffs() throws IOException {
		for (RSQSimEvent event : events) {
			if (event.getMagnitude() > 7 && Math.random() < 0.001)
				testPlotRupAzDiff(event);
		}
		
		while (true);
	}
	
	private void testPlotRupAzDiff(RSQSimEvent event) throws IOException {
		EqkRupture rup = getGMPE_Rup(event);
		
		Location loc1 = rup.getRuptureSurface().getFirstLocOnUpperEdge();
		Location loc2 = rup.getRuptureSurface().getLastLocOnUpperEdge();
		LocationList line = new LocationList();
		line.add(loc1);
		line.add(loc2);
		
		Region reg = new Region(line, 100);
		
		GriddedRegion gridReg = new GriddedRegion(reg, 0.02, null);
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
		
		for (int i=0; i<xyz.size(); i++)
			xyz.set(i, calcRupAzimuthDiff(event, rup, xyz.getLocation(i)));
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		PlotCurveCharacterstics surfChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY);
		
		XY_DataSet hypoXY = new DefaultXY_DataSet();
		Location hypo = rup.getHypocenterLocation();
		hypoXY.set(hypo.getLongitude(), hypo.getLatitude());
		
		RupturePlotGenerator.addElementOutline(funcs, chars, event.getAllElements(), surfChar, reg);
		
		funcs.add(hypoXY);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_INV_TRIANGLE, 7f, Color.WHITE));
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, 180d);
		
		XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, "Event "+event.getID(), "Latitude", "Longitude", "Azimuth");
		spec.setXYElems(funcs);
		spec.setXYChars(chars);
		Range xRange = new Range(xyz.getMinLon(), xyz.getMaxLon());
		Range yRange = new Range(xyz.getMinLat(), xyz.getMaxLat());
		
		GraphWindow gw = new GraphWindow(spec, xRange, yRange);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	public void generateRotDRatioPage(File outputDir, double[] aggregatedPeriods, double[] scatterPeriods, AttenRelRef gmpeRef,
			List<EventComparison> gmpeComps) throws IOException {
		Preconditions.checkState(bbpZipFile.hasRotD100());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		String distDescription = getDistDescription();
		
		// header
		lines.add("# "+catalog.getName()+" BBP RotD100/RotD50 Ratios");
		lines.add("");
		lines.add("Distance dependent plots use the "+distDescription+" distance metric");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		
		IMT[] aggregatedIMTs = new IMT[aggregatedPeriods.length];
		for (int i=0; i<aggregatedIMTs.length; i++)
			aggregatedIMTs[i] = IMT.forPeriod(aggregatedPeriods[i]);
		
		if (gmpeComps == null) {
			gmpeComps = loadCalcComps(gmpeRef, aggregatedIMTs);
		} else {
			// might not have all periods required
			for (BBP_Site site : sites) {
				Map<Integer, EventComparison> compMap = new HashMap<>();
				for (EventComparison comp : gmpeComps)
					compMap.put(comp.getRupture().getID(), comp);
				calcGMPE(compMap, sitesBBPtoGMPE.get(site), gmpeRef, aggregatedIMTs);
			}
		}
		
		super.generateRotDRatioPage(outputDir, lines, aggregatedPeriods, scatterPeriods, gmpeRef, gmpeComps);
	}
	
	static final boolean DIST_JB = true;
	
	public static void main(String[] args) throws ZipException, IOException {
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_5413.instance();
		RSQSimCatalog catalog = Catalogs.BRUCE_5652.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_5566_CRUSTAL.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_5585_SUB.instance();
		
		boolean doGMPE = true;
		boolean doRotD = false;
		boolean doNonErgodicMaps = false;
		
		Region siteMapRegion = new CaliforniaRegions.CYBERSHAKE_MAP_REGION();
		Region sourceMapRegion = bufferRegion(siteMapRegion, 200d);
		
		boolean doGridded = true;
		
		double timeScale = 1d;
		boolean scaleVelocities = false;

//		VelocityModel forceVM = VelocityModel.LA_BASIN_500;
//		VelocityModel forceVM = VelocityModel.LA_BASIN_863;
		VelocityModel forceVM = null;
		
		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
//		IMT[] imts = { IMT.SA3P0 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.ASK_2014 };
		IMT[] imts = { IMT.PGV, IMT.SA2P0, IMT.SA3P0, IMT.SA5P0, IMT.SA10P0 };
		AttenRelRef rotDGMPE = AttenRelRef.ASK_2014;
		
//		AttenRelRef[] gmpeRefs = { AttenRelRef.AFSHARI_STEWART_2016 };
//		IMT[] imts = { IMT.DUR_5_75, IMT.DUR_5_95, IMT.DUR_20_80 };
//		AttenRelRef rotDGMPE = null;
		
//		AttenRelRef[] gmpeRefs = { AttenRelRef.ZHAO_2006 };
//		IMT[] imts = { IMT.SA2P0, IMT.SA3P0, IMT.SA5P0 };
//		AttenRelRef rotDGMPE = null;
		
		String[] highlightNames;
		if (doGridded)
			highlightNames = new String[0];
		else
			highlightNames = new String[] { "USC", "SBSM" };
		
		RuptureIdentifier loadIden = null;
		String loadIdenPrefix = null;
		
//		RuptureIdentifier loadIden = FocalMechIden.builder().strikeSlip(10d).forDip(90).build();
//		String loadIdenPrefix = "mech_vert_ss";
//		RuptureIdentifier loadIden = FocalMechIden.builder().forRake(75, 105).forDip(35, 55).build();
//		String loadIdenPrefix = "mech_reverse";
//		RuptureIdentifier loadIden = FocalMechIden.builder().forRake(-105, -75).forDip(35, 55).build();
//		String loadIdenPrefix = "mech_normal";

		boolean replotScatters = true;
		boolean replotZScores = true;
		boolean replotCurves = true;
		boolean replotResiduals = true;
		
//		boolean replotScatters = true;
//		boolean replotZScores = true;
//		boolean replotCurves = true;
//		boolean replotResiduals = true;
		
		if (timeScale != 1d)
			doRotD = false;
		
		boolean skipRGdirs = true;
		boolean rgOnlyIfPossible = true;
		
		double[] rotDPeriods = { 1, 2, 5, 7.5, 10 };
		
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		VelocityModel bbpVM = null;
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-all-m") && !name.contains("-all-maxDist")) {
				if (skipRGdirs && name.contains("-rg"))
					continue;
				if (doGridded && !name.contains("-gridded"))
					continue;
				if (!doGridded && name.contains("-gridded"))
					continue;
				if (timeScale == 1d && name.contains("-timeScale"))
					continue;
				if (timeScale != 1d) {
					if (!name.contains("-timeScale"+(float)timeScale))
						continue;
					if (scaleVelocities && !name.contains("-velScale"))
						continue;
					if (!scaleVelocities && name.contains("-velScale"))
						continue;
				}
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					System.out.println("Found candidate zip file: "+zipFile.getAbsolutePath());
					VelocityModel myVM = null;
					List<BBP_Site> sites = BBP_Site.readFile(dir);
					for (VelocityModel v : VelocityModel.values()) {
						if (name.contains(v.name())) {
							myVM = v;
							System.out.println("\tdectected VM: "+myVM);
							break;
						}
						if ((float)sites.get(0).getVs30() == v.getVs30()) {
							myVM = v;
							System.out.println("\tassuming VM from Vs30: "+myVM);
						}
					}
					if (forceVM != null && myVM != forceVM) {
						System.out.println("Skipping dir, wrong VM");
						continue;
					}
					
					bbpDir = dir;
					bbpZipFile = zipFile;
					bbpVM = myVM;
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("Velocity Model: "+bbpVM);
		
		File gmpeCacheDir = new File(bbpDir, "gmpe_cache");
		Preconditions.checkState(gmpeCacheDir.exists() || gmpeCacheDir.mkdir());
		
		String bbpDirName = bbpDir.getName();
		double minMag;
		if (bbpDirName.contains("-all-m")) {
			String magStr = bbpDirName.substring(bbpDirName.indexOf("-all-m")+"-all-m".length());
			if (magStr.contains("-"))
				magStr = magStr.substring(0, magStr.indexOf("-"));
			minMag = Double.parseDouble(magStr);
			System.out.println("Detected minMag="+minMag);
		} else {
			throw new IllegalStateException("Couldn't detect minMag from "+bbpDirName);
		}
		int skipYears = 0;
		if (bbpDirName.contains("-skipYears")) {
			String yearStr = bbpDirName.substring(bbpDirName.indexOf("-skipYears")+"-skipYears".length());
			if (yearStr.contains("-"))
				yearStr = yearStr.substring(0, yearStr.indexOf("-"));
			skipYears = Integer.parseInt(yearStr);
			System.out.println("Detected skipYears="+skipYears);
		}
		
		double maxDist = MPJ_BBP_CatalogSim.CUTOFF_DIST_DEFAULT;
		if (bbpDirName.contains("-maxDist")) {
			String distStr = bbpDirName.substring(bbpDirName.indexOf("-maxDist")+"-maxDist".length());
			if (distStr.contains("-"))
				distStr = distStr.substring(0, distStr.indexOf("-"));
			maxDist = Double.parseDouble(distStr);
			System.out.println("Detected maxDist="+maxDist);
		}
		
		boolean hasRG = bbpDirName.contains("-rg");
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		System.out.println("Zip file: "+bbpZipFile.getAbsolutePath());
		ZipFile zipFile = new ZipFile(bbpZipFile);
		
		CatalogGMPE_Compare comp = new CatalogGMPE_Compare(catalog, zipFile, sites, minMag, skipYears,
				gmpeCacheDir, loadIden, bbpVM, maxDist);
		comp.setReplotCurves(replotCurves);
		comp.setReplotResiduals(replotResiduals);
		comp.setReplotScatters(replotScatters);
		comp.setReplotZScores(replotZScores);
		System.out.println("Has RotD100? "+comp.bbpZipFile.hasRotD100());
		doRotD = doRotD && comp.bbpZipFile.hasRotD100();
		
		for (AttenRelRef ref : gmpeRefs)
			comp.checkTransSiteParams(ref);
//		comp.testPlotRupAzDiffs();
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File vmDir = new File(catalogOutputDir, "bbp_"+bbpVM.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		List<? extends FaultSection> subSects = null;
		Map<RSQSimEvent, List<FaultSection>> rupSectMappings = null;
		Map<RSQSimEvent, FaultSection> rupSectNucleations = null;
		
		try {
			List<EventComparison> rotD_GMPEcomps = null;
			if (doGMPE || doNonErgodicMaps) {
				if (!hasRG || !rgOnlyIfPossible) {
					for (AttenRelRef gmpeRef : gmpeRefs) {
						if (highlightNames != null) {
							if (gmpeRef == rotDGMPE)
								comp.setHighlightSites(highlightNames);
							else
								comp.setHighlightSites();
						}
						List<EventComparison> comps = comp.loadCalcComps(gmpeRef, imts);
						
						if (doGMPE) {
							String dirname = "gmpe_bbp_comparisons_"+gmpeRef.getShortName();
							if (doGridded)
								dirname += "_GriddedSites";
							if (timeScale != 1d) {
								dirname += "_timeScale"+(float)timeScale;
								if (scaleVelocities)
									dirname += "_velScale";
							}
							if (loadIden != null && loadIdenPrefix.length() > 0)
								dirname += "_"+loadIdenPrefix;
							File catalogGMPEDir = new File(vmDir, dirname);
							Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
							comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, imts, comps);
							if (gmpeRef == rotDGMPE)
								rotD_GMPEcomps = comps;
						}
						if (doNonErgodicMaps) {
							String dirname = "gmpe_bbp_non_ergodic_maps_"+gmpeRef.getShortName();
							File catalogGMPEDir = new File(vmDir, dirname);
							Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
							
							if (subSects == null) {
								subSects = catalog.getSubSects();
								rupSectMappings = new HashMap<>();
								rupSectNucleations = new HashMap<>();
								for (RSQSimEvent event : comp.events) {
									RSQSimSubSectEqkRupture mappedRup = catalog.getMappedSubSectRupture(event);
									rupSectMappings.put(event, new ArrayList<>(mappedRup.getSubSections()));
									rupSectNucleations.put(event, mappedRup.getNucleationSection());
								}
							}
							comp.generateNonErgodicMapPage(catalogGMPEDir, null, subSects, siteMapRegion,
									sourceMapRegion, gmpeRef, comps, rupSectMappings, rupSectNucleations, imts, comp.highlightSites);
						}
					}
				}
				if (hasRG) {
					// TODO
//					System.out.println("Has Rupture Generator!");
//					gmpeCacheDir = new File(bbpDir, "gmpe_cache_bbp_rg");
//					Preconditions.checkState(gmpeCacheDir.exists() || gmpeCacheDir.mkdir());
//					comp.setDoRupGen(true, gmpeCacheDir);
//					for (AttenRelRef gmpeRef : gmpeRefs) {
//						File catalogGMPEDir = new File(catalogOutputDir, "gmpe_bbp_rg_comparisons_"+gmpeRef.getShortName());
//						Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
//						comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, periods);
//					}
				}
			}
			if (doRotD) {
				File catalogRotDDir = new File(vmDir, "catalog_rotd_ratio_comparisons");
				Preconditions.checkState(catalogRotDDir.exists() || catalogRotDDir.mkdir());
				List<Double> allPeriods = new ArrayList<>();
				for (IMT imt : imts)
					if (imt.getParamName().equals(SA_Param.NAME))
						allPeriods.add(imt.getPeriod());
				comp.generateRotDRatioPage(catalogRotDDir, rotDPeriods, Doubles.toArray(allPeriods), rotDGMPE, rotD_GMPEcomps);
			}
			
			catalog.writeMarkdownSummary(catalogOutputDir, true, false);
			RSQSimCatalog.writeCatalogsIndex(outputDir);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		comp.shutdown();
	}

}
