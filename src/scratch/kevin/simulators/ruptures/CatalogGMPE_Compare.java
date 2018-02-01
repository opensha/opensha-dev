package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
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
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simCompare.GroundMotionScatterPlot;
import scratch.kevin.simCompare.MultiRupGMPE_ComparePageGen;
import scratch.kevin.simCompare.RuptureComparison;
import scratch.kevin.simCompare.RuptureComparisonFilter;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import scratch.kevin.simulators.hazard.HazardMapComparePlotter;
import scratch.kevin.util.MarkdownUtils;
import scratch.kevin.util.MarkdownUtils.TableBuilder;

class CatalogGMPE_Compare extends MultiRupGMPE_ComparePageGen<RSQSimEvent> {
	
	private RSQSimCatalog catalog;
	private BBP_CatalogSimZipLoader bbpZipFile;
	private List<BBP_Site> sites;
	private List<Site> highlightSites;
	private double minMag;
	private double minFractForInclusion;

	private BiMap<BBP_Site, Site> sitesBBPtoGMPE;
	private BiMap<Site, BBP_Site> sitesGMPEtoBBP;
	private List<RSQSimEvent> events;
	private double catDurationYears;
	private Map<Integer, RSQSimEvent> eventsMap;
	private Map<Integer, EqkRupture> gmpeRupsMap;
	
	private boolean rupGen = false;
	
	private File gmpeCacheDir;

	public CatalogGMPE_Compare(RSQSimCatalog catalog, ZipFile bbpZipFile, List<BBP_Site> sites, double minMag, int skipYears,
			VelocityModel vm, double minFractForInclusion, File gmpeCacheDir) throws IOException {
		this.catalog = catalog;
		this.sites = sites;
		this.minMag = minMag;
		this.minFractForInclusion = minFractForInclusion;
		this.gmpeCacheDir = gmpeCacheDir;
		
		ArrayList<RegionIden> siteRegIdens = new ArrayList<>();
		sitesBBPtoGMPE = HashBiMap.create();
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site site : sites) {
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), MPJ_BBP_CatalogSim.CUTOFF_DIST)));
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
		loader.matches(new LogicalOrRupIden(siteRegIdens));
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
	
	public void setDoRupGen(boolean rupGen, File gmpeCacheDir) {
		this.rupGen = rupGen;
		this.gmpeCacheDir = gmpeCacheDir;
		gmpeRupsMap.clear();
	}
	
	public void setHighlightSites(String... highlightNames) {
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
		synchronized(event) {
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
					rup = catalog.getGMPE_Rupture(event, minFractForInclusion);
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
			Double period = Double.parseDouble(csv.get(row, 1));
			Double logMean = Double.parseDouble(csv.get(row, 2));
			Double stdDev = Double.parseDouble(csv.get(row, 3));
			Double rRup = Double.parseDouble(csv.get(row, 4));
			Double rJB = Double.parseDouble(csv.get(row, 5));
			
			EventComparison comp = eventComps.get(eventID);
			if (comp != null) {
				comp.addResult(site, period, logMean, stdDev);
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
			for (Double period : sorted(comp.getPeriods(site))) {
				csv.addLine(eventID+"", period+"", comp.getLogMean(site, period)+"", comp.getStdDev(site, period)+"",
						comp.getDistanceRup(site)+"", comp.getDistanceJB(site)+"");
			}
		}
		
		csv.writeToFile(file);
	}
	
	private <E extends Comparable<E>> Iterable<E> sorted(Set<E> vals) {
		List<E> list = new ArrayList<>(vals);
		Collections.sort(list);
		return list;
	}
	
	public synchronized void calcGMPE(Map<Integer, EventComparison> eventComps, Site site, AttenRelRef gmpeRef, double... periods) {
		List<Future<?>> futures = new ArrayList<>();
		
		File gmpeCacheFile = null;
		if (gmpeCacheDir != null) {
			String siteName = sitesGMPEtoBBP.get(site).getName();
			gmpeCacheFile = new File(gmpeCacheDir, siteName+"_"+gmpeRef.getShortName()+".csv");
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
		
		List<Double> periodsList = Doubles.asList(periods);
		
		for (EventComparison comp : eventComps.values()) {
			if (!comp.isSiteApplicable(site))
				continue;
			List<Double> myPeriods;
			Set<Double> computedPeriods = null;
			if (comp.hasSite(site))
				computedPeriods = comp.getPeriods(site);
			if (computedPeriods == null || computedPeriods.isEmpty()) {
				// nothing cached
				myPeriods = periodsList;
			} else {
				myPeriods = new ArrayList<>();
				for (double period : periods)
					if (!computedPeriods.contains(period))
						myPeriods.add(period);
			}
			if (!myPeriods.isEmpty())
				futures.add(exec.submit(new EventComparisonCalc(comp, gmpeRef, site, myPeriods)));
		}
		
		System.out.println("Calculating for "+futures.size()+" events");
		
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (Exception e) {
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
		private List<Double> periods;
		
		public EventComparisonCalc(EventComparison eventComp, AttenRelRef gmpeRef, Site site, List<Double> periods) {
			this.eventComp = eventComp;
			this.gmpeRef = gmpeRef;
			this.site = site;
			this.periods = periods;
		}

		@Override
		public void run() {
			ScalarIMR gmpe = checkOutGMPE(gmpeRef);
			
			eventComp.calculate(gmpe, site, Doubles.toArray(periods));
			
			checkInGMPE(gmpeRef, gmpe);
		}
		
	}
	
	private class EventComparison extends RuptureComparison.Cached<RSQSimEvent> {
		
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
	}
	
	public List<EventComparison> loadCalcComps(AttenRelRef gmpeRef, double[] periods) {
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
			calcGMPE(compsMap, sitesBBPtoGMPE.get(site), gmpeRef, periods);
		}
		return new ArrayList<>(compsMap.values());
	}
	
	public void generateGMPE_Page(File outputDir, AttenRelRef gmpeRef, double[] periods, List<EventComparison> comps)
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
		
		super.generateGMPE_Page(outputDir, lines, gmpeRef, periods, comps, highlightSites);
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
		
		XYZPlotWindow gw = new XYZPlotWindow(spec, xRange, yRange);
		gw.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
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
		
		if (gmpeComps == null) {
			gmpeComps = loadCalcComps(gmpeRef, aggregatedPeriods);
		} else {
			// might not have all periods required
			for (BBP_Site site : sites) {
				Map<Integer, EventComparison> compMap = new HashMap<>();
				for (EventComparison comp : gmpeComps)
					compMap.put(comp.getRupture().getID(), comp);
				calcGMPE(compMap, sitesBBPtoGMPE.get(site), gmpeRef, aggregatedPeriods);
			}
		}
		
		super.generateRotDRatioPage(outputDir, lines, aggregatedPeriods, scatterPeriods, gmpeRef, gmpeComps);
	}
	
	static final boolean DIST_JB = true;
	
	public static void main(String[] args) throws ZipException, IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
//		RSQSimCatalog catalog = Catalogs.JG_modLoad_testB.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_2495.instance(baseDir);
		
		boolean doGMPE = true;
		boolean doRotD = true;
		
		String[] highlightNames = { "USC", "SBSM" };
//		String[] highlightNames = null;
		
		VelocityModel vm = VelocityModel.LA_BASIN;
		double minFractForInclusion = 0.2;
		boolean skipRGdirs = true;
		boolean rgOnlyIfPossible = true;
		
//		double[] periods = { 1, 2, 3, 5, 10 };
		double[] periods = { 1, 2, 5 };
		double[] rotDPeriods = { 1, 2, 5, 7.5, 10 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014,
//				AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.BSSA_2014 };
		AttenRelRef rotDGMPE = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-all")) {
				if (skipRGdirs && name.contains("-rg"))
					continue;
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					bbpDir = dir;
					bbpZipFile = zipFile;
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		
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
		
		boolean hasRG = bbpDirName.contains("-rg");
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		System.out.println("Zip file: "+bbpZipFile.getAbsolutePath());
		ZipFile zipFile = new ZipFile(bbpZipFile);
		
		CatalogGMPE_Compare comp = new CatalogGMPE_Compare(catalog, zipFile, sites, minMag, skipYears,
				vm, minFractForInclusion, gmpeCacheDir);
		System.out.println("Has RotD100? "+comp.bbpZipFile.hasRotD100());
		doRotD = doRotD && comp.bbpZipFile.hasRotD100();
		if (highlightNames != null)
			comp.setHighlightSites(highlightNames);
//		comp.testPlotRupAzDiffs();
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		try {
			List<EventComparison> rotD_GMPEcomps = null;
			if (doGMPE) {
				if (!hasRG || !rgOnlyIfPossible) {
					for (AttenRelRef gmpeRef : gmpeRefs) {
						List<EventComparison> comps = comp.loadCalcComps(gmpeRef, periods);
						File catalogGMPEDir = new File(catalogOutputDir, "gmpe_bbp_comparisons_"+gmpeRef.getShortName());
						Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
						comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, periods, comps);
						if (gmpeRef == rotDGMPE)
							rotD_GMPEcomps = comps;
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
				File catalogRotDDir = new File(catalogOutputDir, "catalog_rotd_ratio_comparisons");
				Preconditions.checkState(catalogRotDDir.exists() || catalogRotDDir.mkdir());
				comp.generateRotDRatioPage(catalogRotDDir, rotDPeriods, periods, rotDGMPE, rotD_GMPEcomps);
			}
			
			catalog.writeMarkdownSummary(catalogOutputDir, true, false);
			RSQSimCatalog.writeCatalogsIndex(outputDir);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		comp.shutdown();
	}

}
