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
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import scratch.kevin.simulators.hazard.HazardMapComparePlotter;

class CatalogGMPE_Compare {
	
	private RSQSimCatalog catalog;
	private BBP_CatalogSimZipLoader bbpZipFile;
	private List<BBP_Site> sites;
	private List<BBP_Site> highlightSites;
	private double minMag;
	private double minFractForInclusion;
	
	private Map<BBP_Site, Site> gmpeSites;
	private List<RSQSimEvent> events;
	private double catDurationYears;
	private Map<Integer, RSQSimEvent> eventsMap;
	private Map<Integer, EqkRupture> gmpeRupsMap;
	private Table<RSQSimEvent, BBP_Site, Double> rupSiteAzMap;
	
	private Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesInstancesCache;
	private ExecutorService exec;
	
	private boolean rupGen = false;
	
	private File gmpeCacheDir;
	
	private RSQSimBBP_HazardCurveCalc bbpCurveCalc = null;
	
	private List<Range> magRanges;
	private List<String> magLabels;
	private List<String> magFileLabels;
	private List<Range> distRanges;
	private List<String> distLabels;
	private List<String> distFileLabels;

	public CatalogGMPE_Compare(RSQSimCatalog catalog, ZipFile bbpZipFile, List<BBP_Site> sites, double minMag, int skipYears,
			VelocityModel vm, double minFractForInclusion, File gmpeCacheDir) throws IOException {
		this.catalog = catalog;
		this.bbpZipFile = new BBP_CatalogSimZipLoader(bbpZipFile, sites);
		this.sites = sites;
		this.minMag = minMag;
		this.minFractForInclusion = minFractForInclusion;
		this.gmpeCacheDir = gmpeCacheDir;
		
		ArrayList<RegionIden> siteRegIdens = new ArrayList<>();
		gmpeSites = new HashMap<>();
		for (BBP_Site site : sites) {
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), MPJ_BBP_CatalogSim.CUTOFF_DIST)));
			gmpeSites.put(site, site.buildGMPE_Site(vm));
		}
		
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
		
		gmpesInstancesCache = new HashMap<>();
		gmpeRupsMap = new HashMap<>();
		rupSiteAzMap = HashBasedTable.create();
		
		exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		
		initBins();
	}
	
	private void initBins() {

		
		double maxCatalogMag = 0d;
		for (RSQSimEvent e : events)
			maxCatalogMag = Math.max(maxCatalogMag, e.getMagnitude());
		
		magRanges = new ArrayList<>();
		if (minMag < 6d) {
			magRanges.add(new Range(minMag, 6d));
			magRanges.add(new Range(6d, 6.5d));
			magRanges.add(new Range(6.5d, 7d));
		} else if (minMag < 6.5d) {
			magRanges.add(new Range(minMag, 6.5d));
			magRanges.add(new Range(6.5d, 7d));
		} else if (minMag < 7d) {
			magRanges.add(new Range(minMag, 7d));
		}
		magRanges.add(new Range(7d, 7.5d));
		magRanges.add(new Range(7.5, 8d));
		if (maxCatalogMag > 8d) {
			if (maxCatalogMag < 8.5)
				magRanges.add(new Range(8d, 8.5));
			else
				magRanges.add(new Range(8d, 9));
		}
		magLabels = new ArrayList<>();
		for (Range magRange : magRanges)
			magLabels.add(optionalDigitDF.format(magRange.getLowerBound())+" < Mw < "
					+optionalDigitDF.format(magRange.getUpperBound()));
		magFileLabels = new ArrayList<>();
		for (Range magRange : magRanges)
			magFileLabels.add("mag_"+optionalDigitDF.format(magRange.getLowerBound())+"_"
					+optionalDigitDF.format(magRange.getUpperBound()));
		
		distRanges = new ArrayList<>();
		// 0-10, 10-20, 20-40, 40-80, 80-160
		distRanges.add(new Range(0d, 10d));
		distRanges.add(new Range(10d, 20d));
		distRanges.add(new Range(20d, 40d));
		distRanges.add(new Range(40d, 80d));
		distRanges.add(new Range(80d, 160d));
		distRanges.add(new Range(160d, MPJ_BBP_CatalogSim.CUTOFF_DIST));
		distLabels = new ArrayList<>();
		String distShortName = getDistShortName();
		for (Range distRange : distRanges)
			distLabels.add(optionalDigitDF.format(distRange.getLowerBound())+" km < "+distShortName+" < "
					+optionalDigitDF.format(distRange.getUpperBound())+" km");
		distFileLabels = new ArrayList<>();
		for (Range distRange : distRanges)
			distFileLabels.add("dist_"+optionalDigitDF.format(distRange.getLowerBound())+"_"
					+optionalDigitDF.format(distRange.getUpperBound()));
	}
	
	public void setDoRupGen(boolean rupGen, File gmpeCacheDir) {
		this.rupGen = rupGen;
		this.gmpeCacheDir = gmpeCacheDir;
		gmpeRupsMap.clear();
	}
	
	public void shutdown() {
		exec.shutdown();
	}
	
	public void setHighlightSites(String... highlightNames) {
		List<BBP_Site> highlightSites = new ArrayList<>();
		for (String name : highlightNames) {
			for (BBP_Site site : sites) {
				if (name.startsWith(site.getName())) {
					highlightSites.add(site);
					break;
				}
			}
		}
		setHighlightSites(highlightSites);
	}
	
	public void setHighlightSites(List<BBP_Site> highlightSites) {
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
	
	private ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
		synchronized (gmpesInstancesCache) {
			LinkedList<ScalarIMR> gmpes = gmpesInstancesCache.get(gmpeRef);
			if (gmpes == null) {
				gmpes = new LinkedList<>();
				gmpesInstancesCache.put(gmpeRef, gmpes);
			}
			if (!gmpes.isEmpty())
				return gmpes.pop();
		}
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		return gmpe;
	}
	
	private void checkInGMPE(AttenRelRef gmpeRef, ScalarIMR gmpe) {
		synchronized (gmpesInstancesCache) {
			gmpesInstancesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	private int loadCache(File file, BBP_Site site, AttenRelRef gmpeRef, Map<Integer, EventComparison> eventComps) throws IOException {
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
	
	private void writeCache(File file, BBP_Site site, Map<Integer, EventComparison> comps) throws IOException {
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
	
	public synchronized void calcGMPE(Map<Integer, EventComparison> eventComps, BBP_Site site, AttenRelRef gmpeRef, double... periods) {
		List<Future<?>> futures = new ArrayList<>();
		
		File gmpeCacheFile = null;
		if (gmpeCacheDir != null) {
			gmpeCacheFile = new File(gmpeCacheDir, site.getName()+"_"+gmpeRef.getShortName()+".csv");
			if (gmpeCacheFile.exists()) {
				try {
					System.out.println("Loading cache from "+gmpeCacheFile.getAbsolutePath());
					int added = loadCache(gmpeCacheFile, site, gmpeRef, eventComps);
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
		private BBP_Site site;
		private List<Double> periods;
		
		public EventComparisonCalc(EventComparison eventComp, AttenRelRef gmpeRef, BBP_Site site, List<Double> periods) {
			this.eventComp = eventComp;
			this.gmpeRef = gmpeRef;
			this.site = site;
			this.periods = periods;
		}

		@Override
		public void run() {
			ScalarIMR gmpe = checkOutGMPE(gmpeRef);
			
			SA_Param imt = (SA_Param)gmpe.getIntensityMeasure();
			SA_Param.setPeriodInSA_Param(imt, periods.get(0));
			EqkRupture gmpeRup = getGMPE_Rup(eventComp.getEvent());
			gmpe.setAll(gmpeRup, gmpeSites.get(site), imt);
			
			RuptureSurface surf = gmpeRup.getRuptureSurface();
			double distanceRup = surf.getDistanceRup(site.getLoc());
			double distanceJB = surf.getDistanceJB(site.getLoc());
			
			eventComp.setDistances(site, distanceRup, distanceJB);
			
			for (double period : periods) {
				SA_Param.setPeriodInSA_Param(imt, period);
				double logMean = gmpe.getMean();
				double stdDev = gmpe.getStdDev();
				eventComp.addResult(site, period, logMean, stdDev);
			}
			
			checkInGMPE(gmpeRef, gmpe);
		}
		
	}
	
	private class EventComparison {
		private final RSQSimEvent event;
		
		private List<BBP_Site> applicableSites;
		
		private Table<BBP_Site, Double, Double> logMeans;
		private Table<BBP_Site, Double, Double> stdDevs;
		private Map<BBP_Site, Double> distanceRups;
		private Map<BBP_Site, Double> distanceJBs;
		
		public EventComparison(RSQSimEvent event) {
			this.event = event;
			
			applicableSites = new ArrayList<>();
			
			logMeans = HashBasedTable.create();
			stdDevs = HashBasedTable.create();
			distanceRups = new HashMap<>();
			distanceJBs = new HashMap<>();
		}
		
		public void addApplicableSite(BBP_Site site) {
			applicableSites.add(site);
		}
		
		public boolean isSiteApplicable(BBP_Site site) {
			return applicableSites.contains(site);
		}
		
		public void addResult(BBP_Site site, double period, double logMean, double stdDev) {
			logMeans.put(site, period, logMean);
			stdDevs.put(site, period, stdDev);
		}
		
		public void setDistances(BBP_Site site, double distanceRup, double distanceJB) {
			distanceJBs.put(site, distanceJB);
			distanceRups.put(site, distanceRup);
		}
		
		public RSQSimEvent getEvent() {
			return event;
		}
		
		public double getLogMean(BBP_Site site, double period) {
			return logMeans.get(site, period);
		}
		
		public double getStdDev(BBP_Site site, double period) {
			return stdDevs.get(site, period);
		}
		
		public double getDistanceRup(BBP_Site site) {
			return distanceRups.get(site);
		}
		
		public double getDistanceJB(BBP_Site site) {
			return distanceJBs.get(site);
		}
		
		public boolean hasSite(BBP_Site site) {
			return logMeans.containsRow(site);
		}
		
		public boolean isComputed(BBP_Site site, double period) {
			return logMeans.contains(site, period);
		}
		
		public Set<Double> getPeriods(BBP_Site site) {
			return logMeans.row(site).keySet();
		}
	}
	
	private abstract class EventComparisonFilter {
		public abstract boolean matches(EventComparison comp, BBP_Site site);
	}
	
	private class DistJBFilter extends EventComparisonFilter {
		
		private double min, max;

		public DistJBFilter(double min, double max) {
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean matches(EventComparison comp, BBP_Site site) {
			if (!comp.hasSite(site))
				return false;
			double dist = comp.getDistanceJB(site);
			return dist >= min && dist <= max;
		}
		
	}
	
	private class DistRupFilter extends EventComparisonFilter {
		
		private double min, max;

		public DistRupFilter(double min, double max) {
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean matches(EventComparison comp, BBP_Site site) {
			if (!comp.hasSite(site))
				return false;
			double dist = comp.getDistanceRup(site);
			return dist >= min && dist <= max;
		}
		
	}
	
	public Map<Integer, EventComparison> buildEventComps() {
		Map<Integer, EventComparison> ret = new HashMap<>();
		
		for (BBP_Site site : sites) {
			for (Integer eventID : bbpZipFile.getEventIDs(site)) {
				EventComparison comp = ret.get(eventID);
				if (comp == null) {
					comp = new EventComparison(eventsMap.get(eventID));
					ret.put(eventID, comp);
				}
				comp.addApplicableSite(site);
			}
		}
		
		return ret;
	}
	
	public Map<Integer, EventComparison> forSite(Map<Integer, EventComparison> comps, BBP_Site site) {
		Map<Integer, EventComparison> ret = new HashMap<>();
		
		for (Integer eventID : comps.keySet()) {
			EventComparison comp = comps.get(eventID);
			if (comp.isSiteApplicable(site))
				ret.put(eventID, comp);
		}
		
		return ret;
	}
	
	public Map<Integer, EventComparison> forMag(Map<Integer, EventComparison> comps, double minMag, double maxMag) {
		Map<Integer, EventComparison> ret = new HashMap<>();
		
		for (Integer eventID : comps.keySet()) {
			EventComparison comp = comps.get(eventID);
			double mag = comp.getEvent().getMagnitude();
			if (mag >= minMag && mag <= maxMag)
				ret.put(eventID, comp);
		}
		
		return ret;
	}
	
	public boolean plotScatter(Collection<EventComparison> eventComps, BBP_Site site, double period, AttenRelRef gmpe,
			EventComparisonFilter filter, List<String> binDescriptions, File outputDir, String prefix)
					throws IOException {
		List<BBP_Site> sites;
		if (site == null) {
			sites = this.sites;
		} else {
			sites = new ArrayList<>();
			sites.add(site);
		}
		return plotScatter(eventComps, sites, period, gmpe, filter, binDescriptions, outputDir, prefix);
	}
	
	public boolean plotScatter(Collection<EventComparison> eventComps, List<BBP_Site> sites, double period, AttenRelRef gmpe,
			EventComparisonFilter filter, List<String> binDescriptions, File outputDir, String prefix)
					throws IOException {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		for (BBP_Site site : sites) {
			for (EventComparison comp : eventComps) {
				if (!comp.isComputed(site, period) || (filter != null && !filter.matches(comp, site)))
					continue;
				DiscretizedFunc[] spectras;
				if (rupGen)
					spectras = bbpZipFile.readRupGenRotD50(site, comp.getEvent().getID());
				else
					spectras = new DiscretizedFunc[] { bbpZipFile.readRotD50(site, comp.getEvent().getID()) };
				for (DiscretizedFunc spectra : spectras) {
					double rsVal = spectra.getY(period);
					double gmpeVal = Math.exp(comp.getLogMean(site, period));
					xy.set(gmpeVal, rsVal);
				}
			}
		}
		if (xy.size() == 0)
			return false;
		
		String title = gmpe.getShortName()+" Comparison Scatter";
		String xAxisLabel = gmpe.getShortName()+" "+(float)period+" s SA (g)";
		String yAxisLabel;
		if (rupGen)
			yAxisLabel = "BBP Rupture Generator "+(float)period+" s SA (g)";
		else
			yAxisLabel = "RSQSim "+(float)period+" s SA (g)";
		
		GroundMotionScatterPlot.PLOT_WIDTH = 600;
		GroundMotionScatterPlot.WRITE_PDF = false;
		GroundMotionScatterPlot.plot(xy, xAxisLabel, yAxisLabel, binDescriptions, title, outputDir, prefix);
		return true;
	}
	
	public boolean plotStandardNormal(Collection<EventComparison> eventComps, BBP_Site site, double[] periods, AttenRelRef gmpe,
			EventComparisonFilter filter, List<String> binDescriptions, File outputDir, String prefix)
			throws IOException {
		List<BBP_Site> sites;
		if (site == null) {
			sites = this.sites;
		} else {
			sites = new ArrayList<>();
			sites.add(site);
		}
		return plotStandardNormal(eventComps, sites, periods, gmpe, filter, binDescriptions, outputDir, prefix);
	}
	
	public boolean plotStandardNormal(Collection<EventComparison> eventComps, List<BBP_Site> sites, double[] periods,
			AttenRelRef gmpe, EventComparisonFilter filter, List<String> binDescriptions, File outputDir, String prefix)
					throws IOException {
		
		List<PlotSpec> specs = new ArrayList<>();
		double maxY = 0.7d;
		double numStdDev = 3.75;
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(new Range(-numStdDev, numStdDev));
		
		List<Double> means = new ArrayList<>();
		List<Double> stdDevs = new ArrayList<>();
		int numMatches = 0;
		for (BBP_Site site : sites) {
			for (EventComparison comp : eventComps) {
				if (!comp.isComputed(site, periods[0]) || (filter != null && !filter.matches(comp, site)))
					continue;
				numMatches++;
			}
		}
		
		int numBins;
		if (numMatches < 100)
			numBins = 10;
		else if (numMatches < 500)
			numBins = 40;
		else
			numBins = 100;
		
		Color stdDevColor = new Color(0, 150, 0);
		
		for (double period : periods) {
			HistogramFunction hist = new HistogramFunction(-numStdDev, numStdDev, numBins);
			
			double mean = 0d;
			List<Double> allVals = new ArrayList<>();
			int count = 0;
			for (BBP_Site site : sites) {
				for (EventComparison comp : eventComps) {
					if (!comp.isComputed(site, period) || (filter != null && !filter.matches(comp, site)))
						continue;
					DiscretizedFunc[] spectras;
					if (rupGen)
						spectras = bbpZipFile.readRupGenRotD50(site, comp.getEvent().getID());
					else
						spectras = new DiscretizedFunc[] { bbpZipFile.readRotD50(site, comp.getEvent().getID()) };
					for (DiscretizedFunc spectra : spectras) {
						// in log space
						double rsVal = Math.log(spectra.getY(period));
						double gmpeVal = comp.getLogMean(site, period);
						
						double val = (rsVal - gmpeVal)/comp.getStdDev(site, period);
						hist.add(hist.getClosestXIndex(val), 1d);
						mean += val;
						allVals.add(val);
						count++;
					}
				}
			}
			if (count == 0)
				return false;
			mean /= count;
			means.add(mean);
			double stdDev = Math.sqrt(StatUtils.variance(Doubles.toArray(allVals), mean));
			stdDevs.add(stdDev);
			
			double area = 0d;
			for (Point2D pt : hist)
				area += hist.getDelta()*pt.getY();
			hist.scale(1d/area);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			EvenlyDiscretizedFunc stdNormal = new EvenlyDiscretizedFunc(hist.getMinX(), hist.getMaxX(), 1000);
			double scalar = 1d/Math.sqrt(2d*Math.PI);
			for (int i=0; i<stdNormal.size(); i++) {
				double x = stdNormal.getX(i);
				double y = scalar*Math.exp(-0.5*x*x);
				stdNormal.set(i, y);
			}
			
			funcs.add(hist);
			if (rupGen)
				hist.setName("BBP Rupture Generator");
			else
				hist.setName("RSQSim");
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			funcs.add(stdNormal);
			stdNormal.setName("Standard Normal");
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
//			maxY = Math.max(maxY, Math.max(stdNormal.getMaxY(), hist.getMaxY()));
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			meanLine.set(mean, 0);
			meanLine.set(mean, maxY-0.1);
			meanLine.setName("Mean");
			funcs.add(meanLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
			
			for (double sigma=Math.ceil(-numStdDev); sigma<=numStdDev; sigma++) {
				DefaultXY_DataSet sigmaLine = new DefaultXY_DataSet();
				sigmaLine.set(sigma, 0);
				sigmaLine.set(sigma, maxY-0.075);
				funcs.add(sigmaLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, stdDevColor));
			}
			
			String title = gmpe.getShortName()+" Log-Normal Comparision";
			String xAxisLabel = "z-score (Standard Deviations)";
			String yAxisLabel = "Density";
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			spec.setLegendVisible(period == periods[periods.length-1]);
			
			specs.add(spec);
		}
		
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<periods.length; i++) {
			List<String> labels = new ArrayList<>(binDescriptions);
			labels.add(0, optionalDigitDF.format(periods[i])+"s SA");
			
			double yEach = maxY/8d;
			double x = -numStdDev + 0.2;
			double y = maxY - yEach*1.2;
			
			Font bigFont = new Font(Font.SANS_SERIF, Font.BOLD, 24);
			Font smallFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
			
			List<XYAnnotation> anns = new ArrayList<>();
			XYTextAnnotation meanAnn = new XYTextAnnotation("Mean = "+optionalDigitDF.format(means.get(i)), -x, y);
			meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			meanAnn.setFont(bigFont);
			anns.add(meanAnn);
			XYTextAnnotation stdDevAnn = new XYTextAnnotation("σ = "+optionalDigitDF.format(stdDevs.get(i)), -x, y-yEach);
			stdDevAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			stdDevAnn.setFont(bigFont);
			anns.add(stdDevAnn);
			
			for (int j=0; j<labels.size(); j++) {
				String label = labels.get(j);
				XYTextAnnotation ann = new XYTextAnnotation(label, x, y);
				y -= yEach;
				ann.setTextAnchor(TextAnchor.TOP_LEFT);
				if (j == 0) {
					ann.setFont(bigFont);
					yEach *= (double)smallFont.getSize()/(double)bigFont.getSize();
				} else {
					ann.setFont(smallFont);
				}
				anns.add(ann);
			}
			
			for (double sigma=Math.ceil(-numStdDev); sigma<=numStdDev; sigma++) {
				int s = (int)Math.round(Math.abs(sigma));
				String label;
				if (sigma < -0.1)
					label = "-"+s+" σ";
				else if (s == 0)
					label = s+" σ";
				else
					label = "+"+s+" σ";
				XYTextAnnotation ann = new XYTextAnnotation(label, sigma, maxY);
				ann.setTextAnchor(TextAnchor.TOP_CENTER);
				ann.setFont(bigFont);
				ann.setPaint(stdDevColor);
				anns.add(ann);
			}
			
			specs.get(i).setPlotAnnotations(anns);
			yRanges.add(new Range(0, maxY));
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 200 + 300*specs.size());
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		return true;
	}
	
	public List<File> plotHazardCurves(List<BBP_Site> sites, double period, AttenRelRef gmpeRef, File outputDir)
			throws IOException {
		List<Future<File>> futures = new ArrayList<>();
		
		for (BBP_Site site : sites)
			futures.add(exec.submit(new CurveCalcCallable(site, period, gmpeRef, outputDir)));
		
		List<File> files = new ArrayList<>();
		
		for (Future<File> future : futures) {
			try {
				files.add(future.get());
			} catch (Exception e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		return files;
	}
	
	private class CurveCalcCallable implements Callable<File> {
		private BBP_Site site;
		private double period;
		private AttenRelRef gmpeRef;
		private File outputDir;

		public CurveCalcCallable(BBP_Site site, double period, AttenRelRef gmpeRef, File outputDir) {
			this.site = site;
			this.period = period;
			this.gmpeRef = gmpeRef;
			this.outputDir = outputDir;
		}

		@Override
		public File call() throws Exception {
			return plotHazardCurve(site, period, gmpeRef, outputDir);
		}
		
	}
	
	private static int[] hazard_curve_rps = { 1000, 2500, 10000 };
	private static double[] gmpe_truncs = { 0d, 3d, 2d, 1d };
	private static PlotLineType[] gmpe_trunc_line_types = { PlotLineType.SOLID, PlotLineType.DASHED,
			PlotLineType.DOTTED, PlotLineType.DOTTED_AND_DASHED };
	
	@SuppressWarnings("unchecked")
	public File plotHazardCurve(BBP_Site site, double period, AttenRelRef gmpeRef, File outputDir)
			throws IOException {
		String xAxisLabel = (float)period+"s SA (g)";
		String yAxisLabel = "Annual Probability";
		double curveDuration = 1d;
		
		if (bbpCurveCalc == null)
			bbpCurveCalc = new RSQSimBBP_HazardCurveCalc(bbpZipFile, catDurationYears);
		
		System.out.println("Calculating BBP curve for "+site.getName()+", "+xAxisLabel);
		DiscretizedFunc bbpCurve = bbpCurveCalc.calc(site, period, curveDuration, rupGen);
		
		double rateEach = 1d/catDurationYears;
		double probEach = curveDuration/catDurationYears;
		
		// now calculate GMPE
		System.out.println("Calculating "+gmpeRef.getShortName()+" curve for "+site.getName()+", "+xAxisLabel);
		ScalarIMR gmpe = checkOutGMPE(gmpeRef);
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		DiscretizedFunc[] gmpeCurves = new DiscretizedFunc[gmpe_truncs.length];
		for (int t=0; t<gmpeCurves.length; t++) {
			gmpeCurves[t] = bbpCurve.deepClone();
			// init to 1, non-exceedance curves
			for (int i=0; i<gmpeCurves[t].size(); i++)
				gmpeCurves[t].set(i, 1d);
		}
		DiscretizedFunc gmpeMeanCurve = bbpCurve.deepClone();
		for (int i=0; i<gmpeMeanCurve.size(); i++)
			gmpeMeanCurve.set(i, 0);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : gmpeCurves[0])
			logXVals.set(Math.log(pt.getX()), 1d);
		logXVals = new LightFixedXFunc(logXVals);
		Site gmpeSite = gmpeSites.get(site);
		for (int eventID : bbpZipFile.getEventIDs(site)) {
			gmpe.setEqkRupture(getGMPE_Rup(eventsMap.get(eventID)));
			gmpe.setSite(gmpeSite);
			for (int i=0; i<gmpe_truncs.length; i++) {
				double trunc = gmpe_truncs[i];
				if (trunc > 0) {
					gmpe.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
					gmpe.getParameter(SigmaTruncLevelParam.NAME).setValue(trunc);
				} else {
					gmpe.getParameter(SigmaTruncTypeParam.NAME).setValue(SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_NONE);
				}
				gmpe.getExceedProbabilities(logXVals);
				for(int k=0; k<logXVals.size(); k++)
					gmpeCurves[i].set(k, gmpeCurves[i].getY(k)*Math.pow(1d-probEach, logXVals.getY(k)));
			}
			double mean = gmpe.getMean();
			Preconditions.checkState(Doubles.isFinite(mean));
			for(int k=0; k<logXVals.size(); k++)
				if (mean >= logXVals.getX(k))
					gmpeMeanCurve.set(k, gmpeMeanCurve.getY(k) + rateEach);
		}
		// convert to exceedance probabilities
		for (DiscretizedFunc gmpeCurve : gmpeCurves)
			for (int i=0; i<gmpeCurve.size(); i++)
				gmpeCurve.set(i, 1d-gmpeCurve.getY(i));
		// now mean curve -> probabilities
		for (int i=0; i<gmpeMeanCurve.size(); i++) {
			double rate = gmpeMeanCurve.getY(i);
			double prob = 1d - Math.exp(-rate*curveDuration);
			gmpeMeanCurve.set(i, prob);
		}
		checkInGMPE(gmpeRef, gmpe);
		
		// now plot
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = new Range(1e-3, 1e1);
		Range yRange = new Range(1e-8, 1e0);
		
		for (int i=0; i<gmpeCurves.length; i++) {
			DiscretizedFunc gmpeCurve = gmpeCurves[i];
			String name;
			float thickness;
			if (gmpe_truncs[i] == 0) {
				name = gmpeRef.getShortName();
				thickness = 3f;
			} else {
				name = optionalDigitDF.format(gmpe_truncs[i])+" σ";
				thickness = 2f;
			}
			gmpeCurve.setName(name);
			funcs.add(gmpeCurve);
			chars.add(new PlotCurveCharacterstics(gmpe_trunc_line_types[i % gmpe_trunc_line_types.length], thickness, Color.BLUE));
		}
		
		gmpeMeanCurve.setName("Mean Only");
		funcs.add(gmpeMeanCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GREEN.darker()));
		
		bbpCurve.setName("RSQSim/BBP");
		funcs.add(bbpCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		if (hazard_curve_rps != null && hazard_curve_rps.length > 0) {
			CPT rpCPT = HazardMapComparePlotter.getRPlogCPT(hazard_curve_rps);
			for (int rp : hazard_curve_rps) {
				Color color = rpCPT.getColor((float)Math.log10(rp));
				double probLevel = 1d/(double)rp;
				DiscretizedFunc probLine = new ArbitrarilyDiscretizedFunc();
				probLine.set(xRange.getLowerBound(), probLevel);
				probLine.set(xRange.getUpperBound(), probLevel);
				probLine.setName(rp+"yr");
				funcs.add(probLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color));
			}
		}
		
		String siteName = RSQSimBBP_Config.siteCleanName(site);

		PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Hazard Curves", xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		String prefix = siteName.replaceAll(" ", "_")+"_curves_"+(float)period+"s_"+gmpeRef.getShortName();
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(outputDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		System.out.println("DONE "+gmpeRef.getShortName()+", "+site.getName()+", "+xAxisLabel);
		return pngFile;
	}
	
	private synchronized double getRupAzimuthDiff(RSQSimEvent event, BBP_Site site) {
		Double az = rupSiteAzMap.get(event, site);
		if (az == null) {
			az = calcRupAzimuthDiff(event, getGMPE_Rup(event), site.getLoc());
			rupSiteAzMap.put(event, site, az);
		}
		return az;
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
		double hypoDist1 = LocationUtils.horzDistanceFast(hypo, loc1);
		double hypoDist2 = LocationUtils.horzDistanceFast(hypo, loc2);
		if (hypoDist2 < hypoDist1) {
			// flip so that hypocenter is closest to loc1
			Location t = loc1;
			loc1 = loc2;
			loc2 = t;
		}
		
//		double locDist1 = LocationUtils.horzDistanceFast(loc, loc1);
//		double locDist2 = LocationUtils.horzDistanceFast(loc, loc2);
//		
//		double refAz;
//		Location[] refLine;
//		if (locDist1 < locDist2) {
//			refAz = LocationUtils.azimuth(loc1, centroid);
//			refLine = new Location[] {loc1, centroid};
//		} else {
//			refAz = LocationUtils.azimuth(centroid, loc2);
//			refLine = new Location[] {centroid, loc2};
//		}
//		
//		double locAz = LocationUtils.azimuth(refLine[0], loc1);
//		
//		return locAz - refAz;
		
		double refAz = LocationUtils.azimuth(loc1, loc2);
		double locAz = LocationUtils.azimuth(centroid, loc);
		double diff1 = Math.abs(locAz - refAz);
		double diff2 = Math.abs((360+locAz) - refAz);
		double diff3 = Math.abs(locAz - (360+refAz));
		double diff = Math.min(diff1, Math.min(diff2, diff3));
		Preconditions.checkState(diff <= 180);
		return diff;
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
	
	public boolean plotAggregateRotDRatio(Range magRange, File outputDir, String prefix) throws IOException {
		List<DiscretizedFunc[]> rotDs = new ArrayList<>();
		
		for (BBP_Site site : sites) {
			for (Integer eventID : bbpZipFile.getEventIDs(site)) {
				RSQSimEvent event = eventsMap.get(eventID);
				if (magRange != null && !magRange.contains(event.getMagnitude()))
					continue;
				
//				double az = getRupAzimuthDiff(event, site);
//				if (azRange != null && !azRange.contains(az))
//					continue;
				
				rotDs.add(bbpZipFile.readRotD(site, eventID));
			}
		}
		if (rotDs.isEmpty())
			return false;
		
		String title = "RotD100/50 Ratio";
		if (magRange == null)
			title += ", All Mags";
		else
			title += ", "+optionalDigitDF.format(magRange.getLowerBound())
				+"≤M≤"+optionalDigitDF.format(magRange.getUpperBound());
		
		SpectraPlotter.plotRotDRatio(rotDs, "RSQSim", title, outputDir, prefix);
		
		return true;
	}
	
	public boolean plotRotDRatioDependence(Range magRange, Range distRange, double[] periods, boolean scatter, File outputDir, String prefix,
			boolean doDist, boolean doAzimuth, boolean doRotD50) throws IOException {
		Preconditions.checkState(doDist || doAzimuth || doRotD50);
		Preconditions.checkState(!(doDist && doAzimuth) && !(doAzimuth && doRotD50) && !(doDist && doRotD50));
		
		List<DiscretizedFunc[]> rotDs = new ArrayList<>();
		List<Double> scalars = new ArrayList<>();
		
		for (BBP_Site site : sites) {
			for (Integer eventID : bbpZipFile.getEventIDs(site)) {
				RSQSimEvent event = eventsMap.get(eventID);
				if (magRange != null && !magRange.contains(event.getMagnitude()))
					continue;
				
				double dist = Double.NaN;
				if (distRange != null || doDist) {
					RuptureSurface surf = getGMPE_Rup(event).getRuptureSurface();
					if (DIST_JB)
						dist = surf.getDistanceJB(site.getLoc());
					else
						dist = surf.getDistanceRup(site.getLoc());
					if (distRange != null && !distRange.contains(dist))
						continue;
				}
				
				if (doDist) {
					scalars.add(dist);
				} else if (doAzimuth) {
					scalars.add(getRupAzimuthDiff(event, site));
				} else if (doRotD50) {
					// filled in later for each period
					scalars = null;
				}
				
				rotDs.add(bbpZipFile.readRotD(site, eventID));
			}
		}
		if (rotDs.isEmpty())
			return false;
		
		String title = "RotD100/50";
		
		String scalarLabel;
		if (doDist) {
			scalarLabel = getDistShortName();
			title += " Distance Dependence";
		} else if (doAzimuth) {
			scalarLabel = "Azimuth";
			title += " Azimuth Dependence";
		} else if (doRotD50) {
			scalarLabel = "RotD50";
			title += " Amplitude Dependence";
		} else {
			throw new IllegalStateException();
		}
		
		if (magRange == null)
			title += ", All Mags";
		else
			title += ", "+optionalDigitDF.format(magRange.getLowerBound())
				+"≤M≤"+optionalDigitDF.format(magRange.getUpperBound());
		
		if (scatter) {
			for (double period : periods) {
				String myPrefix = prefix+"_"+optionalDigitDF.format(period)+"s";
				if (doRotD50) {
					Preconditions.checkState(scalars == null);
					SpectraPlotter.plotRotDRatioVsRd50Scatter(rotDs, period, "RSQSim", title, outputDir, myPrefix, 10);
				} else {
					SpectraPlotter.plotRotDRatioScatter(rotDs, scalars, scalarLabel, period,
							"RSQSim", title, outputDir, myPrefix, false, 10);
				}
			}
		} else {
			SpectraPlotter.plotRotDRatioDependence(rotDs, scalars, scalarLabel, 10, periods,
					"RSQSim", title, outputDir, prefix, doRotD50);
		}
		
		return true;
	}
	
	private String getDistDescription() {
		if (DIST_JB)
			return "Joyner-Boore distance (**rJB**), the shortest horizontal distance from a site to the "
					+ "surface projection of the rupture surface";
		return "rupture surface distance (**rRup**), the shortest 3-d distance from a site to the "
					+ "rupture surface";
	}
	
	private String getDistShortName() {
		if (DIST_JB)
			return "rJB";
		return "rRup";
	}
	
	private static int max_table_fig_columns = 3;
	
	public void generateGMPE_Page(File outputDir, AttenRelRef gmpeRef, double[] periods)
			throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
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
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		// precalc for all sites
		System.out.println("Precalculating/caching for all sites");
		Map<Integer, EventComparison> eventComps = buildEventComps();
		for (BBP_Site site : sites) {
			System.out.println("Calculating for "+site.getName()+", "+gmpeRef.getShortName());
			calcGMPE(eventComps, site, gmpeRef, periods);
		}
		
		List<BBP_Site> sites = new ArrayList<>();
		if (highlightSites == null)
			// make plots for all sites
			sites.addAll(this.sites);
		else
			// make plots for only the specified, use all for aggregations
			sites.addAll(highlightSites);
		// add null for all sites aggregated
		sites.add(0, null);
		
		for (BBP_Site site : sites) {
			Map<Integer, EventComparison> siteEventComps;
			String siteName;
			if (site == null) {
				System.out.println("All sites aggregated");
				lines.add("## All Sites Aggregated");
				lines.add(topLink); lines.add("");
				lines.add("**"+this.sites.size()+" sites**");
				lines.add("");
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine("Name", "Location", "# Events");
				for (BBP_Site s : this.sites) {
					table.initNewLine();
					table.addColumn(RSQSimBBP_Config.siteCleanName(s));
					Location loc = s.getLoc();
					table.addColumn("*"+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
					table.addColumn(bbpZipFile.getEventIDs(s).size()+"");
					table.finalizeLine();
				}
				lines.addAll(table.build());
				siteEventComps = eventComps;
				siteName = "All Sites";
				
				lines.add(siteEventComps.size()+" ruptures within "+(float)MPJ_BBP_CatalogSim.CUTOFF_DIST+" km of *any* site");
			} else {
				System.out.println("Site: "+site.getName());
				siteName = RSQSimBBP_Config.siteCleanName(site);
				lines.add("## Site "+site.getName());
				lines.add(topLink); lines.add("");
				Location loc = site.getLoc();
				lines.add("*Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
				siteEventComps = forSite(eventComps, site);
				
				lines.add(siteEventComps.size()+" ruptures within "+(float)MPJ_BBP_CatalogSim.CUTOFF_DIST+" km");
			}
			
			for (int m=0; m<magRanges.size(); m++) {
				Range magRange = magRanges.get(m);
				lines.add("### "+siteName+", "+magLabels.get(m));
				
				Map<Integer, EventComparison> magEventComps = forMag(siteEventComps, magRange.getLowerBound(), magRange.getUpperBound());
				lines.add(magEventComps.size()+" Ruptures");
				
				List<EventComparisonFilter> filters = new ArrayList<>();
				for (Range distRange : distRanges) {
					if (DIST_JB)
						filters.add(new DistJBFilter(distRange.getLowerBound(), distRange.getUpperBound()));
					else
						filters.add(new DistRupFilter(distRange.getLowerBound(), distRange.getUpperBound()));
				}
				
				lines.add("#### "+siteName+", "+magLabels.get(m)+", Scatter Plots");
				lines.add(topLink); lines.add("");
				lines.add("**Legend**");
				lines.add("* Red +: GMPE Mean/RSQSim single rupture comparison");
				lines.add("* Yellow Region: Factor of 2 above & below");
				lines.add("* Green Line: Linear Regression");
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.initNewLine().addColumn("**Distance Bin**");
				for (double period : periods)
					table.addColumn("**"+optionalDigitDF.format(period)+" s**");
				table.finalizeLine();
				
				for (int d=0; d<distRanges.size(); d++) {
					EventComparisonFilter filter = filters.get(d);
//					if (eventsForMagDist.isEmpty()) {
//						System.out.println("No events for "+site.getName()+", "+distLabels.get(d)+", "+magLabels.get(m));
//					}
					
					table.initNewLine().addColumn("**"+distLabels.get(d)+"**");
					
					for (double period : periods) {
						String prefix = siteName.replaceAll(" ", "_")+"_"+magFileLabels.get(m)+"_"+distFileLabels.get(d)
						+"_"+optionalDigitDF.format(period)+"s_"+gmpeRef.getShortName()+"_scatter";
				
						System.out.println("Plotting Scatter: "+prefix);
						
						List<String> binDescriptions = Lists.newArrayList(distLabels.get(d),
								magLabels.get(m), optionalDigitDF.format(period)+"s SA, "+gmpeRef.getShortName());
						boolean success = plotScatter(magEventComps.values(), site, period, gmpeRef, filter,
								binDescriptions, resourcesDir, prefix);
						if (success) {
							File scatterPlot = new File(resourcesDir, prefix+".png");
							Preconditions.checkState(scatterPlot.exists());
							table.addColumn("![Scatter Plot]("+resourcesDir.getName()
								+"/"+scatterPlot.getName()+")");
						} else {
							table.addColumn("N/A");
						}
					}
					
					table.finalizeLine();
				}
				lines.add("");
				lines.addAll(table.wrap(max_table_fig_columns, 1).build());
				
				lines.add("#### "+siteName+", "+magLabels.get(m)+", Standard Normal Plots");
				lines.add(topLink); lines.add("");
				lines.add("These plots compare RSQSim to the full GMPE log-normal distributions. "
						+ "Each rupture's GMPE distribution is converted to a standard log-normal "
						+ "distribution, and the z-score is computed for each rupture:");
				lines.add("");
				lines.add("**z-score**: (ln(*RSQSim*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
				lines.add("");
				lines.add("**Legend**");
				lines.add("* Black Line: Standard Normal distribution (in natural log space)");
				lines.add("* Gray Histogram: z-score for each rupture");
				lines.add("* Blue Dashed Line: RSQSim Mean");
				
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				
				for (int d=0; d<distRanges.size(); d++)
					table.addColumn("**"+distLabels.get(d)+"**");
				table.finalizeLine();
				
				table.initNewLine();
				for (int d=0; d<distRanges.size(); d++) {
					EventComparisonFilter filter = filters.get(d);
					String prefix = siteName.replaceAll(" ", "_")+"_"+magFileLabels.get(m)+"_"+distFileLabels.get(d)
						+"_"+gmpeRef.getShortName()+"_std_norm";
					
					System.out.println("Plotting Standard Normal: "+prefix);
					
					List<String> binDescriptions = Lists.newArrayList(distLabels.get(d), magLabels.get(m));
					boolean success = plotStandardNormal(magEventComps.values(), site, periods, gmpeRef, filter,
							binDescriptions, resourcesDir, prefix);
					if (success) {
						File plotFile = new File(resourcesDir, prefix+".png");
						Preconditions.checkState(plotFile.exists());
						table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
							+"/"+plotFile.getName()+")");
					} else {
						table.addColumn("N/A");
					}
				}
				table.finalizeLine();
				lines.add("");
				lines.addAll(table.wrap(max_table_fig_columns, 0).build());
			}
			String prefix = siteName.replaceAll(" ", "_")+"_all_mags_all_dists_"+gmpeRef.getShortName()+"_std_norm";
			boolean success = plotStandardNormal(siteEventComps.values(), site, periods, gmpeRef, null,
					new ArrayList<>(), resourcesDir, prefix);
			if (success) {
				lines.add("### "+siteName+", All Ruptures, Standard Normal Plots");
				lines.add(topLink); lines.add("");
				lines.add("");
				lines.add("z-score standard normal plots across all magnitudes/distances");
				lines.add("");
				lines.add("**z-score**: (ln(*RSQSim*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
				lines.add("");
				lines.add("**Legend**");
				lines.add("* Black Line: Standard Normal distribution (in natural log space)");
				lines.add("* Gray Histogram: z-score for each rupture");
				lines.add("* Blue Dashed Line: RSQSim Mean");
				
				lines.add("");
				File plotFile = new File(resourcesDir, prefix+".png");
				Preconditions.checkState(plotFile.exists());
				lines.add("![Standard Normal Plot]("+resourcesDir.getName()
					+"/"+plotFile.getName()+")");
			}
		}
		
		// now hazard curves
		List<List<File>> curveFiles = new ArrayList<>();
		List<BBP_Site> curveSites = this.sites;
		for (double period : periods)
			curveFiles.add(plotHazardCurves(curveSites, period, gmpeRef, resourcesDir));
		
		lines.add("## Hazard Curves");
		lines.add(topLink); lines.add("");
		lines.add("**Legend**:");
		lines.add("* Black Solid Line: RSQSim/BBP");
		for (int i=0; i<gmpe_truncs.length; i++) {
			String truncAdd = "";
			if (gmpe_truncs[i] > 0)
				truncAdd = " "+optionalDigitDF.format(gmpe_truncs[i])+"-sigma truncation";
			PlotLineType type = gmpe_trunc_line_types[i % gmpe_trunc_line_types.length];
			String lineType = type.name().replaceAll("_", " ");
			lineType = lineType.substring(0, 1).toUpperCase()+lineType.substring(1).toLowerCase();
			lines.add("* Blue "+lineType+" Line: "+gmpeRef.getShortName()+truncAdd);
		}
		lines.add("* Green Dashed Line: "+gmpeRef.getShortName()+" mean values only");
		lines.add("* Gray Dashed Lines: "+Joiner.on(" yr, ").join(Ints.asList(hazard_curve_rps))+" yr return periods");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine().addColumn("Site");
		for (double period : periods)
			table.addColumn(optionalDigitDF.format(period)+"s");
		table.finalizeLine();
		for (int s=0; s<curveSites.size(); s++) {
			BBP_Site site = curveSites.get(s);
			table.initNewLine().addColumn("**"+RSQSimBBP_Config.siteCleanName(site)+"**");
			for (int p=0; p<periods.length; p++) {
				File plotFile = curveFiles.get(p).get(s);
				table.addColumn("![Hazard Curve]("+resourcesDir.getName()
							+"/"+plotFile.getName()+")");
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	public void generateRotDRatioPage(File outputDir, double[] periods) throws IOException {
		Preconditions.checkState(bbpZipFile.hasRotD100());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		String distShortName = getDistShortName();
		String distDescription = getDistDescription();
		
		// header
		lines.add("# "+catalog.getName()+" BBP RotD100/RotD50 Ratios");
		lines.add("");
		lines.add("Distance dependent plots use the "+distDescription+" distance metric");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Site List");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Name", "Location", "# Events");
		for (BBP_Site s : this.sites) {
			table.initNewLine();
			table.addColumn(RSQSimBBP_Config.siteCleanName(s));
			Location loc = s.getLoc();
			table.addColumn("*"+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
			table.addColumn(bbpZipFile.getEventIDs(s).size()+"");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		
		Range totalDistRange = new Range(0d, RSQSimBBP_Config.MAX_DIST);
		
		for (int m=-1; m<magRanges.size(); m++) {
			Range magRange;
			String magLabel;
			String magFileLabel;
			if (m < 0) {
				magRange = null;
				magFileLabel = "all_mags";
				magLabel = "All Mags";
			} else {
				magRange = magRanges.get(m);
				if (magRange.getUpperBound() > 8)
					// not useful here
					continue;
				magFileLabel = magFileLabels.get(m);
				magLabel = magLabels.get(m);
			}

			lines.add("## "+magLabel);
			lines.add(topLink); lines.add("");
			
			// first total
			System.out.println("Plotting "+magLabel+" total RotD ratio");
			String prefix = "rot_d_ratio_"+magFileLabel+"_aggregated";
			Preconditions.checkState(plotAggregateRotDRatio(magRange, resourcesDir, prefix));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", "+distShortName+" Dependence");
			lines.add(topLink); lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio distance dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_"+distShortName+"_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, periods, false, resourcesDir, prefix, true, false, false));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", Azimuth Dependence");
			lines.add(topLink); lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio azimuth dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_az_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, periods, false, resourcesDir, prefix, false, true, false));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", RotD50 Dependence");
			lines.add(topLink); lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio amplitude dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_amplitude_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, periods, false, resourcesDir, prefix, false, false, true));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("#### "+magLabel+", RotD50 Dependence Scatters");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("Distance");
			for (double period : periods)
				table.addColumn("**"+optionalDigitDF.format(period)+"s**");
			table.finalizeLine();
			for (int d=0; d<distRanges.size(); d++) {
				Range distRange = distRanges.get(d);
				String distLabel = distLabels.get(d);
				String distFileLabel = distFileLabels.get(d);
				
				table.initNewLine().addColumn("**"+distLabel+"**");
				
				prefix = "rot_d_ratio_"+magFileLabel+"_"+distFileLabel+"_scatter";
				Preconditions.checkState(plotRotDRatioDependence(magRange, distRange, periods, true, resourcesDir, prefix,
						false, false, true));
				for (double period : periods) {
					File plotFile = new File(resourcesDir, prefix+"_"+optionalDigitDF.format(period)+"s.png");
					Preconditions.checkState(plotFile.exists());
					table.addColumn("![RotD Ratio]("+resourcesDir.getName()+"/"+plotFile.getName()+")");
				}
				
				table.finalizeLine();
			}
			lines.addAll(table.build());
			lines.add("");
//			table.addLine("Azimuth", "RotD100/RotD50 Ratio");
//			for (int a=0; a<azRanges.size(); a++) {
//				Range azRange = azRanges.get(a);
//				prefix = "rot_d_ratio_"+magFileLabel+"_"+azFileLabels.get(a);
//				table.initNewLine();
//				if (azRange == null)
//					table.addColumn("All");
//				else
//					table.addColumn(optionalDigitDF.format(azRange.getLowerBound())
//							+"°≤Az≤"+optionalDigitDF.format(azRange.getUpperBound())+"°");
//				if (plotRotDRatio(magRange, azRange, resourcesDir, prefix))
//					table.addColumn("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
//				else
//					table.addColumn("N/A");
//				table.finalizeLine();
//			}
//			lines.addAll(table.build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	static final boolean DIST_JB = true;
	
	public static void main(String[] args) throws ZipException, IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
		RSQSimCatalog catalog = Catalogs.JG_2194_K2.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2349.instance(baseDir);
		
		boolean doGMPE = false;
		boolean doRotD = true;
		
		String[] highlightNames = { "USC", "SBSM" };
//		String[] highlightNames = null;
		
		VelocityModel vm = VelocityModel.LA_BASIN;
		double minFractForInclusion = 0.2;
		boolean skipRGdirs = true;
		boolean rgOnlyIfPossible = true;
		
//		double[] periods = { 1, 2, 3, 5, 10 };
		double[] periods = { 1, 2, 5 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014,
//				AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.BSSA_2014 };
		
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
		if (highlightNames != null)
			comp.setHighlightSites(highlightNames);
//		comp.testPlotRupAzDiffs();
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		try {
			if (doGMPE) {
				if (!hasRG || !rgOnlyIfPossible) {
					for (AttenRelRef gmpeRef : gmpeRefs) {
						File catalogGMPEDir = new File(catalogOutputDir, "gmpe_bbp_comparisons_"+gmpeRef.getShortName());
						Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
						comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, periods);
					}
				}
				if (hasRG) {
					System.out.println("Has Rupture Generator!");
					gmpeCacheDir = new File(bbpDir, "gmpe_cache_bbp_rg");
					Preconditions.checkState(gmpeCacheDir.exists() || gmpeCacheDir.mkdir());
					comp.setDoRupGen(true, gmpeCacheDir);
					for (AttenRelRef gmpeRef : gmpeRefs) {
						File catalogGMPEDir = new File(catalogOutputDir, "gmpe_bbp_rg_comparisons_"+gmpeRef.getShortName());
						Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
						comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, periods);
					}
				}
			}
			if (doRotD && comp.bbpZipFile.hasRotD100()) {
				File catalogRotDDir = new File(catalogOutputDir, "catalog_rotd_ratio_comparisons");
				Preconditions.checkState(catalogRotDDir.exists() || catalogRotDDir.mkdir());
				comp.generateRotDRatioPage(catalogRotDDir, periods);
			}
			
			catalog.writeMarkdownSummary(catalogOutputDir, true, false);
			RSQSimCatalog.writeCatalogsIndex(outputDir);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		comp.shutdown();
	}

}
