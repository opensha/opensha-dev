package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingDeque;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.SwingUtilities;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol.Symbol;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.MMI_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.etas.ETAS_HazardMapCalc.MapType;

public class ETAS_HazardMapCalc {
	
	private static boolean force_serial = false;
	private static int debugCurvePlotModulus = 0;
	private static int debugStopIndex = 0;
	
	private boolean calcInLogSpace = true;
	private double distCutoff = 200;
	
	private List<List<ETAS_EqkRupture>> catalogs;
	private HashSet<Integer> faultIndexesTriggered;
	protected GriddedRegion region;
	private DiscretizedFunc xVals;
	private DiscretizedFunc calcXVals;
	
	private int printModulus = 100;
	
	// for precalc faults
	private DataInputStream in;
	int faultSiteIndex = 0;
	// for on the fly faults
	private FaultSystemSolution sol;
	private FaultSystemSolutionERF faultERF;
	private ProbEqkSource[] sourcesForFSSRuptures;
	
	private ETAS_CatalogGridSourceProvider gridSources;
	private AttenRelRef gmpeRef;
	protected String imtName;
	private double period;
	private Deque<ScalarIMR> gmpeDeque;
	private BlockingDeque<FaultSystemSolutionERF> erfDeque;
	
	private List<Site> sites;
	
	protected Duration[] durations;
	private long minOT;
	private int startYear;
	private long[] startOTs;
	private long[] endOTs;
	
	private Duration[] longTermCalcDurations;
	
	protected Table<Duration, MapType, DiscretizedFunc[]> curves;
	private int curvesCalculated;
	
	enum MapType {
		FAULT_ONLY("faults"),
		GRIDDED_ONLY("gridded"),
		COMBINED("combined"),
		U3TD("u3-td"),
		U3TI("u3-ti");
		
		final String fileName;
		private MapType(String fileName) {
			this.fileName = fileName;
		}
		
		public boolean isETAS() {
			return this == MapType.FAULT_ONLY || this == MapType.GRIDDED_ONLY || this == MapType.COMBINED;
		}
	}
	
	enum DurationConstants {
		FULL("Full Catalog", "full", 0d),
		DAY("1 Day", "day", 1d/365.25),
		THREE("Days 1-3", "days_0_3", 3d/365.25),
		WEEK("1 Week", "week", 7d/365.25),
		MONTH("1 Month", "month", 1d/12d),
		YEAR("1 Year", "year", 1d);
		
		final Duration duration;
		
		private DurationConstants(String plotName, String fileName, double years) {
			this(plotName, fileName, 0, years);
		}
		
		private DurationConstants(String plotName, String fileName, double startYears, double endYears) {
			this.duration = new Duration(plotName, fileName, startYears, endYears);
		}
	}
	
	static Duration[] getDurationDefaults() {
		List<Duration> defaults = Lists.newArrayList();
		
		// constants
		DurationConstants[] constants = DurationConstants.values();
		for (DurationConstants constant : constants)
			defaults.add(constant.duration);
		
		// 3 day offsets
		for (int startDay=0; startDay<10; startDay++)
			defaults.add(durationForDayRange(startDay, startDay+3));
		
		return defaults.toArray(new Duration[0]);
	}
	
	static Duration durationForDayRange(int startDay, int endDay) {
		Preconditions.checkState(endDay > startDay);
		String plotName;
		if (endDay == (startDay+1))
			plotName = "Day "+(startDay+1);
		else
			plotName = "Days "+(startDay+1)+"-"+endDay;
		return new Duration(plotName, "days_"+startDay+"_"+endDay,
				(double)startDay/365.25, (double)endDay/365.25);
	}
	
	static Duration durationForFileName(String fileName) {
		// first check constants
		for (DurationConstants c : DurationConstants.values())
			if (fileName.contains(c.duration.fileName+".") || fileName.contains(c.duration.fileName+"_"))
				return c.duration;
		// detect day range in format day_0_3 (for days 0-3)
		Preconditions.checkState(fileName.contains("days_"), "Can only detect constants or day ranges");
		fileName = fileName.substring(fileName.indexOf("days_")+5);
		Preconditions.checkState(fileName.contains("_"));
		int startDay = Integer.parseInt(fileName.substring(0, fileName.indexOf("_")));
		fileName = fileName.substring(fileName.indexOf("_")+1);
		String endDayStr = "";
		for (int i=0; i<fileName.length(); i++) {
			if (Character.isDigit(fileName.charAt(i)))
				endDayStr += fileName.charAt(i);
			else
				break;
		}
		int endDay = Integer.parseInt(endDayStr);
		return durationForDayRange(startDay, endDay);
	}
	
	static class Duration {
		final String plotName;
		final String fileName;
		final double startYears;
		final double endYears;
		
		public Duration(String plotName, String fileName, double years) {
			this(plotName, fileName, 0, years);
		}
		
		public Duration(String plotName, String fileName, double startYears, double endYears) {
			this.plotName = plotName;
			this.fileName = fileName;
			this.startYears = startYears;
			this.endYears = endYears;
		}
		
		public long getStartMillis(long catalogOT) {
			return catalogOT + (long)(startYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		}
		
		public long getEndMillis(long catalogOT) {
			if (endYears == 0)
				return Long.MAX_VALUE;
			return catalogOT + (long)(endYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			long temp;
			temp = Double.doubleToLongBits(endYears);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			result = prime * result + ((fileName == null) ? 0 : fileName.hashCode());
			result = prime * result + ((plotName == null) ? 0 : plotName.hashCode());
			temp = Double.doubleToLongBits(startYears);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Duration other = (Duration) obj;
			if (Double.doubleToLongBits(endYears) != Double.doubleToLongBits(other.endYears))
				return false;
			if (fileName == null) {
				if (other.fileName != null)
					return false;
			} else if (!fileName.equals(other.fileName))
				return false;
			if (plotName == null) {
				if (other.plotName != null)
					return false;
			} else if (!plotName.equals(other.plotName))
				return false;
			if (Double.doubleToLongBits(startYears) != Double.doubleToLongBits(other.startYears))
				return false;
			return true;
		}

		@Override
		public String toString() {
			return plotName;
		}
	}
	
	private boolean calcFaults;
	private boolean calcGridded;
	private boolean calcLongTerm;
	
	public ETAS_HazardMapCalc(List<List<ETAS_EqkRupture>> catalogs, GriddedRegion region, DiscretizedFunc xVals,
			File precalcFile, FaultSystemSolution sol, ETAS_CatalogGridSourceProvider gridSources,
			AttenRelRef gmpeRef, String imtName, double period, List<Site> sites, Duration[] durations) throws IOException {
		this.catalogs = catalogs;
		this.region = region;
		initXVals(xVals);
		this.gridSources = gridSources;
		this.gmpeRef = gmpeRef;
		this.imtName = imtName;
		this.period = period;
		this.sites = sites;
		this.sol = sol;
		
		if (durations == null || durations.length == 0)
			durations = new Duration[] {DurationConstants.FULL.duration};
		this.durations = durations;
		// claculate the min OT
		minOT = Long.MAX_VALUE;
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			if (catalog.isEmpty())
				continue;
			long ot = catalog.get(0).getOriginTime();
			if (ot < minOT)
				minOT = ot;
		}
		startYear = ETAS_MultiSimAnalysisTools.calcYearForOT(minOT);
		startOTs = new long[durations.length];
		endOTs = new long[durations.length];
		for (int i=0; i<durations.length; i++) {
			startOTs[i] = durations[i].getStartMillis(minOT);
			endOTs[i] = durations[i].getEndMillis(minOT);
		}
		
		calcGridded = (gridSources != null && gmpeRef != null && imtName != null && sites != null);
		Preconditions.checkState(sites == null || sites.size() == region.getNodeCount());
		
		// this is used to conserve memory and only load results for ruptures actually used
		faultIndexesTriggered = new HashSet<Integer>();
		for (List<ETAS_EqkRupture> catalog : catalogs)
			for (ETAS_EqkRupture rup : catalog)
				if (rup.getFSSIndex() > 0)
					faultIndexesTriggered.add(rup.getFSSIndex());
		
		if (precalcFile != null) {
			// load in precalculated fault data
			in = new DataInputStream(new BufferedInputStream(new FileInputStream(precalcFile)));
			int numSites = in.readInt();
			Preconditions.checkState(numSites == region.getNodeCount(), "Binary file has %s grid nodes, region has %s",
					numSites, region.getNodeCount());
			calcFaults = true;
		} else if (sol != null) {
			// calculate faults on the fly
			Preconditions.checkNotNull(sites);
			faultERF = new FaultSystemSolutionERF(sol);
			faultERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			faultERF.updateForecast();
			// organize by FSS index
			sourcesForFSSRuptures = new ProbEqkSource[sol.getRupSet().getNumRuptures()];
			for (int sourceID=0; sourceID<faultERF.getNumFaultSystemSources(); sourceID++) {
				int fssIndex = faultERF.getFltSysRupIndexForSource(sourceID);
				sourcesForFSSRuptures[fssIndex] = faultERF.getSource(sourceID);
			}
			calcFaults = true;
		}
	}
	
	public ETAS_HazardMapCalc(GriddedRegion region, Table<Duration, MapType, File> curveFiles) throws Exception {
		this(region, curveFiles, null, null, null, Double.NaN, null);
	}
	
	public ETAS_HazardMapCalc(GriddedRegion region, Table<Duration, MapType, File> curveFiles,
			DiscretizedFunc xVals, AttenRelRef gmpeRef, String imtName, double period, List<Site> sites)	throws Exception {
		this.region = region;
		this.gmpeRef = gmpeRef;
		this.imtName = imtName;
		this.period = period;
		this.sites = sites;
		
		curves = HashBasedTable.create();
		List<Duration> durations = Lists.newArrayList();
		List<Duration> longTermDurations = Lists.newArrayList();
		for (Cell<Duration, MapType, File> cell : curveFiles.cellSet()) {
			Duration duration = cell.getRowKey();
			MapType type = cell.getColumnKey();
			System.out.println("Loading "+region.getNodeCount()+" curves for "+duration.plotName+", "+type);
			curves.put(duration, type, loadCurves(cell.getValue()));
			if (type.isETAS())
				durations.add(duration);
			else
				longTermDurations.add(duration);
		}
		this.durations = durations.toArray(new Duration[0]);
		setCalcLongTerm(!longTermDurations.isEmpty());
		if (isCalcLongTerm())
			setLongTermCalcDurations(longTermDurations.toArray(new Duration[0]));
		
		Preconditions.checkArgument(!curves.isEmpty(), "Must supply at least one curve file");
		xVals = curves.values().iterator().next()[0];
		initXVals(xVals);
	}
	
	private void initXVals(DiscretizedFunc xVals) {
		this.xVals = xVals;
		if (calcInLogSpace) {
			calcXVals = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<xVals.size(); i++)
				calcXVals.set(Math.log(xVals.getX(i)), 1d);
		} else {
			calcXVals = xVals;
		}
	}
	
	public void setSol(FaultSystemSolution sol) {
		this.sol = sol;
	}
	
	private ETAS_HazardMapCalc() {}
	
	private DiscretizedFunc[] loadCurves(File curvesFile) throws Exception {
		BinaryHazardCurveReader reader = new BinaryHazardCurveReader(curvesFile.getAbsolutePath());
		
		DiscretizedFunc[] curves = new DiscretizedFunc[region.getNodeCount()];
		
		for (int i=0; i<curves.length; i++) {
			curves[i] = reader.nextLightCurve();
			Location loc = reader.currentLocation();
			Location gridLoc = region.getNodeList().get(i);
			Preconditions.checkState(loc.equals(gridLoc), "Unexpected grid location!\n\tFile: %s\n\tGrid: %s", loc, gridLoc);
			Preconditions.checkNotNull(curves[i]);
		}
		Preconditions.checkState(reader.nextCurve() == null, "More curves than expected!");
		
		return curves;
	}
	
	private ExecutorService createExecutor() {
		int threads = Runtime.getRuntime().availableProcessors();
		if (threads > region.getNodeCount())
			threads = region.getNodeCount();
		return createExecutor(threads);
	}
	
	ExecutorService createExecutor(int threads) {
//		ExecutorService executor = Executors.newFixedThreadPool(threads);
		// max tasks in the pool at any given time, prevents pre loading too much data and using all memory
		// while waiting for hazard calculations to finish. When the queue is full, it will be run in this
		// thread, effectively blocking
		int maxTasks = threads * 10;
		ExecutorService executor = new ThreadPoolExecutor(threads, threads,
                0L, TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<Runnable>(maxTasks), new ThreadPoolExecutor.CallerRunsPolicy());
		return executor;
	}
	
	public boolean isCalcFaults() {
		return calcFaults;
	}

	public void setCalcFaults(boolean calcFaults) {
		this.calcFaults = calcFaults;
	}

	public boolean isCalcGridded() {
		return calcGridded;
	}

	public void setCalcGridded(boolean calcGridded) {
		this.calcGridded = calcGridded;
	}

	public boolean isCalcLongTerm() {
		return calcLongTerm;
	}

	public void setCalcLongTerm(boolean calcLongTerm) {
		this.calcLongTerm = calcLongTerm;
	}
	
	public void setScenarioForElasticRebound(TestScenario scenario) {
		if (scenario.getFSS_Index() < 0)
			return;
		// figure out erf start time
		TimeSpan timeSpan = new TimeSpan(TimeSpan.YEARS, TimeSpan.YEARS);
		timeSpan.setStartTime(startYear);
		long erfStartMillis = timeSpan.getStartTimeInMillis();
		long eventTimeMillis = erfStartMillis - 1000l; // 1 second before
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (FaultSectionPrefData sect : rupSet.getFaultSectionDataForRupture(scenario.getFSS_Index()))
			// use startYear not minOT in case the calc started in the middle of the year
			sect.setDateOfLastEvent(eventTimeMillis);
	}
	
	synchronized Duration[] getLongTermCalcDurations() {
		if (longTermCalcDurations == null) {
			// all of them, except offsets
			List<Duration> durList = Lists.newArrayList();
			for (Duration dur : durations)
				if (dur.startYears == 0)
					durList.add(dur);
			longTermCalcDurations = durList.toArray(new Duration[0]);
		}
		return longTermCalcDurations;
	}
	
	public void setLongTermCalcDurations(Duration[] longTermCalcDurations) {
		this.longTermCalcDurations = longTermCalcDurations;
	}
	
	private Duration[] durationsForType(MapType type) {
		if (type.isETAS())
			return durations;
		Preconditions.checkState(isCalcLongTerm());
		return getLongTermCalcDurations();
	}

	public void calculate() throws IOException {
		if (in != null)
			Preconditions.checkState(faultSiteIndex == 0, "Can only process file once");
		
		curves = HashBasedTable.create();
		
		if (calcFaults) {
			for (Duration duration : durations)
				curves.put(duration, MapType.FAULT_ONLY, new DiscretizedFunc[region.getNodeCount()]);
		}
		if (calcGridded) {
			for (Duration duration : durations)
				curves.put(duration, MapType.GRIDDED_ONLY, new DiscretizedFunc[region.getNodeCount()]);
		}
		if (calcFaults && calcGridded) {
			for (Duration duration : durations)
				curves.put(duration, MapType.COMBINED, new DiscretizedFunc[region.getNodeCount()]);
		}
		if (calcLongTerm) {
			for (Duration duration : getLongTermCalcDurations()) { 
				curves.put(duration, MapType.U3TD, new DiscretizedFunc[region.getNodeCount()]);
				curves.put(duration, MapType.U3TI, new DiscretizedFunc[region.getNodeCount()]);
			}
		}
		
		ExecutorService executor = createExecutor();
		
		List<Future<Integer>> hazardFutures = Lists.newArrayList();
		
		Stopwatch watch = Stopwatch.createStarted();
		System.out.println("Calculating");
		
		curvesCalculated = 0;
		
		for (int index=0; index<region.getNodeCount(); index++) {
			if (index % printModulus == 0)
				System.out.println("Processing site "+index+"/"+region.getNodeCount());
			
			Map<Integer, double[]> precomputedFaultVals = null;
			if (calcFaults && in != null) {
				Preconditions.checkState(faultSiteIndex == index);
				precomputedFaultVals = loadNextSite();
			}
			
			Future<Integer> future = null;
			HazardCalcRunnable runnable = new HazardCalcRunnable(index, precomputedFaultVals);
			if (force_serial) {
				runnable.run();
			} else {
				future = executor.submit(runnable, index);
				hazardFutures.add(future);
			}
			
			if (debugStopIndex > 0 && index == debugStopIndex)
				break;
			
			if (debugCurvePlotModulus > 0 && index % debugCurvePlotModulus == 0) {
				if (future == null)
					plotCurve(index);
				else
					plotCurve(future); // asynchronous
			}
		}
		
		// wait until we're done
		System.out.println("Waiting for hazard calculations to finish");
		for (Future<?> future : hazardFutures) {
			try {
				future.get();
			} catch (InterruptedException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			} catch (ExecutionException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		System.out.println("Done");
		long secs = watch.elapsed(TimeUnit.SECONDS);
		System.out.println("Calculation took "+secs+" secs");
		watch.stop();
		double curvesPerSecond = (double)curvesCalculated/(double)secs;
		System.out.println((float)curvesPerSecond+" curves/sec");
		
		executor.shutdown();
	}
	
	private class HazardCalcRunnable implements Runnable {
		
		private int index;
		private Map<Integer, double[]> precomputedFaultVals;
		
		public HazardCalcRunnable(int index, Map<Integer, double[]> precomputedFaultVals) {
			this.index = index;
			this.precomputedFaultVals = precomputedFaultVals;
		}

		@Override
		public void run() {
			Table<Duration, MapType, DiscretizedFunc> result = calculateCurves(sites.get(index), precomputedFaultVals);
			for (Cell<Duration, MapType, DiscretizedFunc> cell : result.cellSet()) {
				Duration duration = cell.getRowKey();
				MapType type = cell.getColumnKey();
				DiscretizedFunc curve = result.get(duration, type);
				Preconditions.checkState(curve != null, "Curve not calculated for %s, %s! Size=%s", duration, type, curves.size());
//				System.out.println("Inserting curve for "+duration+", "+type);
				curves.get(duration, type)[index] = cell.getValue();
			}
			
			curvesCalculated++;
			if (curvesCalculated % printModulus == 0)
				System.out.println("Calculated "+curvesCalculated+"/"+region.getNodeCount()+" sites");
		}
	}
	
	/**
	 * Calculates a catalog based hazard curves. First computes a conditional exceedence curve for each individual
	 * catalog, then computes a hazard curve across all catalogs by summing the exceedence probabilities scaled
	 * by 1/numCatalogs.
	 * 
	 * @param site
	 * @param precomputedFaultVals
	 * @return HazardCalcResult instance
	 */
	Table<Duration, MapType, DiscretizedFunc> calculateCurves(Site site, Map<Integer, double[]> precomputedFaultVals) {
		Table<Duration, MapType, DiscretizedFunc> curves = getInitializedCurvesMap(xVals, 0d); // linear space
		
		if (calcFaults || calcGridded) {
			// prepare inputs
			ScalarIMR gmpe = null;
			Map<Integer, DiscretizedFunc> faultNonExceeds = null;
			if (calcFaults) {
				if (precomputedFaultVals == null) {
					// calculate them now
					gmpe = checkOutGMPE();
					faultNonExceeds = calcFaultIMs(site, gmpe);
				} else {
					// use precomputed
					faultNonExceeds = calcFaultExceeds(precomputedFaultVals);
				}
				complimentCurve(faultNonExceeds); // they were actually exceeds
			}
			
			Table<Integer, Integer, DiscretizedFunc> griddedNonExceeds = null;
			if (calcGridded) {
				// will calculate on the fly as needed
				griddedNonExceeds = HashBasedTable.create();
				Preconditions.checkState(gridSources.isConditional());
				if (gmpe == null)
					gmpe = checkOutGMPE();
				gmpe.setSite(site);
			}
			
			Preconditions.checkState(calcFaults || calcGridded);
			
			// now actually calculate
			double rateEach = 1d/catalogs.size();
			
			HashSet<Integer> ignoreGriddedNodes = new HashSet<Integer>();
			
			for (List<ETAS_EqkRupture> catalog : catalogs) {
				boolean hasM5 = false;
				for (ETAS_EqkRupture rup : catalog) {
					if (rup.getMag() > 5) {
						hasM5 = true;
						break;
					}
				}
				if (!hasM5)
					continue;
				Table<Duration, MapType, DiscretizedFunc> catCurves = getInitializedCurvesMap(calcXVals, 1d); // log space if applicable
				
				for (ETAS_EqkRupture rup : catalog) {
					DiscretizedFunc condNonExceed; // conditional non-exceedance probabilities
					MapType targetType; // hazard curve to apply this to
					if (rup.getFSSIndex() >= 0) {
						// fault based
						if (!calcFaults)
							continue;
						condNonExceed = faultNonExceeds.get(rup.getFSSIndex());
						if (condNonExceed == null)
							// not within cutoff dist
							continue;
						targetType = MapType.FAULT_ONLY;
					} else {
						// gridded
						if (!calcGridded)
							continue;
						int nodeIndex = gridSources.getNodeIndex(rup);
						int mfdIndex = gridSources.getMagIndex(rup);
						if (nodeIndex < 0 || mfdIndex < 0 || ignoreGriddedNodes.contains(nodeIndex))
							continue;
						double dist = LocationUtils.horzDistanceFast(site.getLocation(), rup.getHypocenterLocation());
						if (dist > distCutoff) {
							ignoreGriddedNodes.add(nodeIndex);
							continue;
						}
						condNonExceed = griddedNonExceeds.get(nodeIndex, mfdIndex);
						targetType = MapType.GRIDDED_ONLY;
						if (condNonExceed == null) {
							// calculate it
							// multiple ruptures with different focal mechanisms
							Iterable<ProbEqkRupture> rups = gridSources.getConditionalRuptures(rup);
							if (rups == null)
								continue;
							
							condNonExceed = calcXVals.deepClone();
							initializeCurve(condNonExceed, 1d);
							double sumRate = 0d;
							for (ProbEqkRupture subRup : rups) {
								double subMag = subRup.getMag();
								Preconditions.checkState(subMag >= rup.getMag()-0.06 && subMag <= rup.getMag()+0.06,
										"Unexpected mag in sub-rupture. Expected %s, got %s", rup.getMag(), subMag);
								gmpe.setEqkRupture(subRup);
								double rupProb = subRup.getProbability();
								double rupRate = -Math.log(1 - rupProb);
								sumRate += rupRate;
//								System.out.println(subRup.getAveRake()+": "+rupProb);
								
								for (int i=0; i<condNonExceed.size(); i++) {
									// TODO doing this right?
									double exceedProb = gmpe.getExceedProbability(condNonExceed.getX(i));
									// scale by the rate of this rupture
									condNonExceed.set(i, condNonExceed.getY(i)*(1-rupRate*exceedProb));
									// this way if treating it as poisson, but since it's an actual occurance, I don't
									// think that we should
//									condNonExceed.set(i, condNonExceed.getY(i)*Math.pow(1-rupProb, exceedProb));
								}
							}
//							System.out.println("SUM: "+sumProb);
//							System.exit(0);
							Preconditions.checkState((float)sumRate == 1f, "Rupture rates don't sum to 1! %s", sumRate);
							
							griddedNonExceeds.put(nodeIndex, mfdIndex, condNonExceed);
						}
					}
					
					long ot = rup.getOriginTime();
					
					// now add the rupture to the appropriate curves
					for (int i=0; i<durations.length; i++) {
						// duration of 0 means all
						if (ot < startOTs[i] || ot >= endOTs[i])
							// rup occurs outside of window, skip
							continue;
						DiscretizedFunc targetCurve = catCurves.get(durations[i], targetType);
						for (int k=0; k<targetCurve.size(); k++) {
							// multiply this into the total non-exceedance probability
							// (get the product of all non-eceedance probabilities)
							targetCurve.set(k, targetCurve.getY(k) * condNonExceed.getY(k));
						}
					}
				}
				
				if (catCurves.containsColumn(MapType.COMBINED)) {
					// build combined catalog curves
					for (Duration duration : catCurves.rowKeySet()) {
						DiscretizedFunc catFaultCurve = catCurves.get(duration, MapType.FAULT_ONLY);
						DiscretizedFunc catGriddedCurve = catCurves.get(duration, MapType.GRIDDED_ONLY);
						DiscretizedFunc catCombinedCurve = catCurves.get(duration, MapType.COMBINED);
						for (int k=0; k<catCombinedCurve.size(); k++)
							catCombinedCurve.set(k, catFaultCurve.getY(k) * catGriddedCurve.getY(k));
					}
				}
				
				for (Cell<Duration, MapType, DiscretizedFunc> cell : catCurves.cellSet()) {
					DiscretizedFunc catCurve = cell.getValue();
					// convert from total non-exceed prob to total exceed prob
					complimentCurve(catCurve);
					
					// add into total curves
					DiscretizedFunc totalCurve = curves.get(cell.getRowKey(), cell.getColumnKey());
					for (int k=0; k<xVals.size(); k++)
						totalCurve.set(k, totalCurve.getY(k) + rateEach*catCurve.getY(k));
				}
			}
			
			if (gmpe != null)
				checkInGMPE(gmpe);
		}
		
		if (calcLongTerm)
			calcLongTerm(site, getLongTermCalcDurations(), curves);
		
		return curves;
	}
	
	private void calcLongTerm(Site site, Duration[] calcDurations, Table<Duration, MapType, DiscretizedFunc> curves) {
		FaultSystemSolutionERF erf = checkOutERF();
		ScalarIMR gmpe = checkOutGMPE();
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		// first calculate a 1 year time independent curve which will be scaled to any other durations
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		DiscretizedFunc tiOneYearRateCurve = calc.getAnnualizedRates(calcLongTerm(erf, gmpe, site, calc), 1d);
		
		for (Duration duration : calcDurations) {
			Preconditions.checkState(duration.startYears == 0d,
					"Long term not enabled for offsets (not needed, will be practically identical)");
//			System.out.println("Calculating long term, duration: "+duration);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
			double myDur = duration.endYears;
			if (duration.endYears == 0)
				myDur = 10d;
			erf.getTimeSpan().setDuration(myDur);
			erf.getTimeSpan().setStartTime(startYear);
			erf.setParameter(HistoricOpenIntervalParam.NAME, (double)(startYear-1875));
			erf.updateForecast();
			
			DiscretizedFunc tdCurve = calcLongTerm(erf, gmpe, site, calc);
			curves.put(duration, MapType.U3TD, tdCurve);
			
			// scale TI curve to current duration
			DiscretizedFunc tiCurve = xVals.deepClone();
			for (int i=0; i<tiCurve.size(); i++) {
				double rate = tiOneYearRateCurve.getY(i); // already in rate space
				double prob = rateToProb(rate, myDur);
				tiCurve.set(i, prob);
			}
			curves.put(duration, MapType.U3TI, tiCurve);
		}
		
		checkInERF(erf);
		checkInGMPE(gmpe);
	}
	
	private static double rateToProb(double rate, double durationYears) {
		return 1 - Math.exp(-rate*durationYears);
	}
	
	private DiscretizedFunc calcLongTerm(FaultSystemSolutionERF erf, ScalarIMR gmpe, Site site, HazardCurveCalculator calc) {
		DiscretizedFunc curve = calcXVals.deepClone();
		calc.getHazardCurve(curve, site, gmpe, erf);
		DiscretizedFunc linearCurve = xVals.deepClone();
		for (int i=0; i<linearCurve.size(); i++)
			linearCurve.set(i, curve.getY(i));
		return linearCurve;
	}
	
	private Map<Integer, DiscretizedFunc> calcFaultIMs(Site site, ScalarIMR gmpe) {
		// used if no precomputed data file
		Map<Integer, DiscretizedFunc> rupVals = Maps.newHashMap();
		for (Integer fssIndex : faultIndexesTriggered) {
			ProbEqkSource source = sourcesForFSSRuptures[fssIndex];
			if (source == null)
				continue;
			Preconditions.checkState(source.getNumRuptures() == 1, "Must be a single rupture source");
			ProbEqkRupture rup = source.getRupture(0);
			
			double minDist = source.getMinDistance(site);
			if (minDist > distCutoff)
				continue;
			
			gmpe.setSite(site);
			gmpe.setEqkRupture(rup);
			
			rupVals.put(fssIndex, gmpe.getExceedProbabilities(calcXVals.deepClone()));
		}
		return rupVals;
	}
	
	private Map<Integer, DiscretizedFunc> calcFaultExceeds(Map<Integer, double[]> precomputedFaultVals) {
		Map<Integer, DiscretizedFunc> rupCondExceeds = Maps.newHashMap();
		for (int rupIndex : precomputedFaultVals.keySet()) {
			DiscretizedFunc condExceed = calcXVals.deepClone(); // log space if applicable
			double[] vals = precomputedFaultVals.get(rupIndex);
			double mean = vals[0];
			double stdDev = vals[1];
			
			for (int i=0; i<condExceed.size(); i++) {
				double exceedProb = AttenuationRelationship.getExceedProbability(
						mean, stdDev, condExceed.getX(i), null, null);
				condExceed.set(i, exceedProb);
			}
			rupCondExceeds.put(rupIndex, condExceed);
		}
		return rupCondExceeds;
	}
	
	private void complimentCurve(Map<Integer, DiscretizedFunc> curves) {
		for (DiscretizedFunc curve : curves.values())
			complimentCurve(curve);
	}
	
	private void complimentCurve(DiscretizedFunc curve) {
		for (int i=0; i<curve.size(); i++)
			curve.set(i, 1d-curve.getY(i));
	}
	
	private Map<Integer, double[]> loadNextSite() throws IOException {
		Map<Integer, double[]> rupVals = loadSiteFromInputStream(in, faultSiteIndex);
		
		faultSiteIndex++;
		return rupVals;
	}
	
	Map<Integer, double[]> loadSiteFromInputStream(DataInputStream in, int expectedIndex) throws IOException {
		int index = in.readInt();
		Preconditions.checkState(index == expectedIndex, "Bad site index. Expected %s, encountered %s", index, faultSiteIndex);
		double lat = in.readDouble();
		double lon = in.readDouble();
		Location myLoc = new Location(lat, lon);
		Location gridLoc = region.getLocation(index);
		Preconditions.checkState(gridLoc.equals(myLoc),
				"Grid locations don't match.\n\tFrom region: %s\n\tFrom file: %s", gridLoc, myLoc);
		int numRups = in.readInt();
		
		Map<Integer, double[]> rupVals = Maps.newHashMap();
		
		int fssIndex;
		double mean, stdDev;
		for (int i=0; i<numRups; i++) {
			fssIndex = in.readInt();
			mean = in.readDouble();
			stdDev = in.readDouble();
			if (!faultIndexesTriggered.contains(fssIndex))
				continue;
			
			Preconditions.checkState(!rupVals.containsKey(fssIndex));
			rupVals.put(fssIndex, new double[] { mean, stdDev });
		}
		
		return rupVals;
	}
	
	private void plotCurve(final Future<Integer> future) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {
				try {
					int index = future.get();
					plotCurve(index);
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	private void plotCurve(int index) {
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		Map<MapType, DiscretizedFunc[]> fullDurationCurves = curves.row(DurationConstants.FULL.duration);
		
		if (fullDurationCurves.containsKey(MapType.FAULT_ONLY)) {
			DiscretizedFunc curve = fullDurationCurves.get(MapType.FAULT_ONLY)[index];
			curve.setName("Fault Based");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		}
		if (fullDurationCurves.containsKey(MapType.GRIDDED_ONLY)) {
			DiscretizedFunc curve = fullDurationCurves.get(MapType.GRIDDED_ONLY)[index];
			curve.setName("Gridded");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		}
		if (fullDurationCurves.containsKey(MapType.COMBINED)) {
			DiscretizedFunc curve = fullDurationCurves.get(MapType.COMBINED)[index];
			curve.setName("Combined");
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Hazard curves for site "+index, imtName, "Probability of Exceedance");
		spec.setLegendVisible(true);
		GraphWindow gw = new GraphWindow(spec);
		gw.setXLog(true);
		gw.setYLog(true);
	}
	
	public void plotComparisons(MapType type, Duration duration, Location loc, double fullSimDuration, String title,
			File outputDir, String prefix) throws IOException {
		int index = region.indexForLocation(loc);
		Preconditions.checkState(index >= 0, "Couldn't detect index for location: "+loc);
		double dist = LocationUtils.horzDistanceFast(loc, region.getLocation(index));
		System.out.println("Mapped curve location error: "+(float)dist+" km");
		
		plotComparisons(type, duration, index, fullSimDuration, title, outputDir, prefix);
	}
	
//	public void plotComparisons(MapType type, Duration duration, int index, FaultSystemSolutionERF erf,
//			int startYear, double fullSimDuration, String title, File outputDir, String prefix) throws IOException {
	
	public void plotComparisons(MapType type, Duration duration, int index, double fullSimDuration, String title,
			File outputDir, String prefix) throws IOException {
		DiscretizedFunc etasCurve = curves.get(duration, type)[index];
		Preconditions.checkNotNull(etasCurve);
		etasCurve.setName("UCERF3-ETAS");
		Duration longTermDuration = getLongTermCompatibleDuration(duration);
		DiscretizedFunc tiCurve = curves.get(longTermDuration, MapType.U3TI)[index];
		Preconditions.checkNotNull(etasCurve);
		tiCurve.setName("UCERF3-TI");
		DiscretizedFunc tdCurve = curves.get(longTermDuration, MapType.U3TD)[index];
		Preconditions.checkNotNull(etasCurve);
		tdCurve.setName("UCERF3-TD");
		DiscretizedFunc etasCombCurve = getETASplusU3curve(etasCurve, tdCurve);
		etasCombCurve.setName("UCERF3-ETAS+TD");
		
//		if (outputDir.getAbsolutePath().contains("-gridded-only"))
//			type = MapType.COMBINED; // actually combined for comparison purposes
//		
//		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
//		switch (type) {
//		case COMBINED:
//			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
//			erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
//			break;
//		case FAULT_ONLY:
//			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
//			break;
//		case GRIDDED_ONLY:
//			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
//			erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
//			break;
//
//		default:
//			break;
//		}
		
		if (prefix == null || prefix.isEmpty())
			prefix = "";
		else
			prefix += "_";
		prefix += type.fileName+"_"+duration.fileName;
		
		double durationYears;
		String durationLabel;
		if (duration == DurationConstants.FULL.duration) {
			durationYears = fullSimDuration;
			durationLabel = (int)durationYears+" year";
		} else {
			durationYears = duration.endYears - duration.startYears;
			durationLabel = duration.plotName;
		}
		
//		ScalarIMR gmpe = checkOutGMPE();
//		Site site = sites.get(index);
//		HazardCurveCalculator calc = new HazardCurveCalculator();
//		
//		// calculate UCERF3-TI
//		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
//		erf.getTimeSpan().setDuration(durationYears);
//		erf.updateForecast();
//		
//		DiscretizedFunc tiCurve = calc.getHazardCurve(calcXVals.deepClone(), site, gmpe, erf);
//		tiCurve = getReplaceXVals(tiCurve, xVals);
//		tiCurve.setName("UCERF3-TI");
//		
//		// calculate UCERF3-TD
//		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
//		erf.setParameter(MagDependentAperiodicityParam.NAME, MagDependentAperiodicityOptions.MID_VALUES);
//		BPTAveragingTypeOptions aveType = BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE;
//		erf.setParameter(BPTAveragingTypeParam.NAME, aveType);
//		erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
//		erf.getTimeSpan().setDuration(durationYears);
//		erf.getTimeSpan().setStartTime(startYear);
//		erf.setParameter(HistoricOpenIntervalParam.NAME, (double)(startYear-1875));
//		erf.updateForecast();
//
//		DiscretizedFunc tdCurve = calc.getHazardCurve(calcXVals.deepClone(), site, gmpe, erf);
//		tdCurve = getReplaceXVals(tdCurve, xVals);
//		tdCurve.setName("UCERF3-TD");
//		
//		checkInGMPE(gmpe);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(tiCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		funcs.add(tdCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
//		funcs.add(etasCurve);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		funcs.add(etasCombCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		String xAxisLabel = imtName;
		if (xAxisLabel.equals(PGA_Param.NAME))
			xAxisLabel += " (g)";
		if (xAxisLabel.equals(PGV_Param.NAME))
			xAxisLabel += " (cm/s)";
		
		PlotSpec spec = new PlotSpec(funcs, chars, title+", "+durationLabel, xAxisLabel, "Probability of Exceedance");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setUserBounds(1e-2, 1e1, 1e-10, 1d);
		
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(22);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		
		gp.drawGraphPanel(spec, true, true);
//		gp.getChartPanel().setSize(1000, 800);
		gp.getChartPanel().setSize(600, 480);
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		gp.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
	}
	
	private static DiscretizedFunc getReplaceXVals(DiscretizedFunc curve, DiscretizedFunc xVals) {
		Preconditions.checkState(curve.size() == xVals.size());
		
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<curve.size(); i++)
			ret.set(xVals.getX(i), curve.getY(i));
		
		return ret;
	}
	
	public static List<Site> fetchSites(GriddedRegion region, ArrayList<SiteData<?>> provs, ScalarIMR gmpe) throws IOException {
		ArrayList<SiteDataValueList<?>> siteData = null;
		SiteTranslator siteTrans = null;
		if (provs != null) {
			System.out.print("Fetching site data...");
			siteData = new OrderedSiteDataProviderList(provs).getAllAvailableData(region.getNodeList());
			System.out.println("DONE.");
			siteTrans = new SiteTranslator();
		}
		
		List<Site> sites = Lists.newArrayList();
		for (int i=0; i<region.getNodeCount(); i++) {
			Site site = new Site(region.getLocation(i));
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>) param.clone());
			if (siteData != null) {
				List<SiteDataValue<?>> mySiteData = Lists.newArrayList();
				for (SiteDataValueList<?> vals : siteData)
					mySiteData.add(vals.getValue(i));
				for (Parameter<?> param : site)
					siteTrans.setParameterValue(param, mySiteData);
			}
			sites.add(site);
		}
		
		return sites;
	}
	
	private synchronized ScalarIMR checkOutGMPE() {
		if (gmpeDeque == null)
			gmpeDeque = new ArrayDeque<ScalarIMR>();
		if (gmpeDeque.isEmpty()) {
			// build a new one
			ScalarIMR gmpe = gmpeRef.instance(null);
			gmpe.setParamDefaults();
			gmpe.setIntensityMeasure(imtName);
			if (imtName.equals(SA_Param.NAME))
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
			return gmpe;
		}
		return gmpeDeque.pop();
	}
	
	private synchronized void checkInGMPE(ScalarIMR gmpe) {
		gmpeDeque.push(gmpe);
	}
	
	private static double mem_gb_per_erf = 1.5d; // actually closer to 0.5, but leave a healthy buffer for other memory usage
	private int maxNumERFs = -1;
	private int erfsBuilt = 0;
	
	private FaultSystemSolutionERF checkOutERF() {
		// need one ERF per thread because of background seismicity and different durations
		synchronized (this) {
			if (erfDeque == null) {
				erfDeque = new LinkedBlockingDeque<FaultSystemSolutionERF>();
				// if we have extreme parallelism, we need to limit the amount of ERFs that can ever exist at once
				Runtime runtime = Runtime.getRuntime();
				double maxMemGB = runtime.maxMemory()/1073741824d;
				maxNumERFs = (int)(maxMemGB/mem_gb_per_erf);
				System.out.println("Max number of ERFs to build: "+maxNumERFs);
				if (maxNumERFs < 4)
					maxNumERFs = 4;
			}
		}
		boolean canBuild = erfsBuilt < maxNumERFs;
		if (erfDeque.isEmpty() && canBuild) {
			// build a new one
			erfsBuilt++;
			Preconditions.checkState(sol != null, "Must supply solution");
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
			erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
			erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
			return erf;
		}
		try {
//			System.out.println("Taking ERF");
			FaultSystemSolutionERF erf = erfDeque.take();
//			System.out.println("Got ERF, returning");
			return erf;
		} catch (InterruptedException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	private synchronized void checkInERF(FaultSystemSolutionERF erf) {
//		System.out.println("Putting ERF");
		erfDeque.push(erf);
	}
	
	public GriddedGeoDataSet calcMap(MapType type1, MapType type2, Duration duration, boolean isProbAt_IML, double level) {
		GriddedGeoDataSet map = new GriddedGeoDataSet(region, false);
		
		Preconditions.checkState(curves.contains(duration, type1), "No curves for %s, %s", duration, type1);
		Duration ratioDuration = null;
		if (type2 != null)
			Preconditions.checkState(curves.contains(getLongTermCompatibleDuration(duration), type2), "No curves for %s, %s", ratioDuration, type2);
		
		for (int i=0; i<map.size(); i++) {
			DiscretizedFunc curve;
			if (type2 != null)
				curve = getETASplusU3curve(type1, type2, duration, i);
			else
				curve = curves.get(duration, type1)[i];
			
			double val = HazardDataSetLoader.getCurveVal(curve, isProbAt_IML, level);
			if (isProbAt_IML)
				Preconditions.checkState(val >= 0d && val <= 1d, "Bad probability (%s) for level (%s)", val, level);
			if (Double.isInfinite(val))
				val = Double.NaN;
			
			map.set(i, val);
		}
		
		return map;
	}
	
	private Duration getLongTermCompatibleDuration(Duration duration) {
		if (duration.startYears == 0)
			return duration;
		// find the no offset version
		float delta = (float)(duration.endYears - duration.startYears);
		for (Duration oDur : getLongTermCalcDurations()) {
			float oDelta = (float)(oDur.endYears - oDur.startYears);
			if (oDelta == delta)
				return oDur;
		}
		return null;
	}
	
	private Table<String, Duration, DiscretizedFunc> combCurvesCache;
	
	private synchronized DiscretizedFunc getETASplusU3curve(MapType type1, MapType type2, Duration duration, int index) {
		if (combCurvesCache == null)
			combCurvesCache = HashBasedTable.create();
		
		String keyStr = type1.name()+"_"+type2.name()+"_"+index;
		if (combCurvesCache.contains(keyStr, duration))
			return combCurvesCache.get(keyStr, duration);
		
		Preconditions.checkState(type1.isETAS(), "numerator must be ETAS");
		Preconditions.checkState(!type2.isETAS(), "denominator must NOT be ETAS");
		
		DiscretizedFunc etasCurve = curves.get(duration, type1)[index];
		DiscretizedFunc u3Curve = curves.get(getLongTermCompatibleDuration(duration), type2)[index];
		DiscretizedFunc combCurve = getETASplusU3curve(etasCurve, u3Curve);
		
		combCurvesCache.put(keyStr, duration, combCurve);
		return combCurve;
	}

	public static DiscretizedFunc getETASplusU3curve(DiscretizedFunc etasCurve, DiscretizedFunc u3Curve) {
		Preconditions.checkState(etasCurve.size() == u3Curve.size());
		
		DiscretizedFunc combCurve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<etasCurve.size(); i++) {
			Preconditions.checkState(etasCurve.getX(i) == u3Curve.getX(i));
			double etasProb = etasCurve.getY(i);
			double u3Prob = u3Curve.getY(i);
			double combProb = 1d - (1d - etasProb)*(1d - u3Prob);
			combCurve.set(etasCurve.getX(i), combProb);
		}
		return combCurve;
	}
	
	private static HashSet<Duration> saveXYZdurations = new HashSet<Duration>();
	static {
		saveXYZdurations.add(DurationConstants.THREE.duration);
	}
	
//	private static boolean shouldSaveXYZ(Duration duration, boolean isProbAt_IML, double level) {
//		if (!saveXYZmaps.contains(duration, level))
//			return false;
//		return saveXYZmaps.get(duration, level) == isProbAt_IML;
//	}
	
	public void plotMap(MapType type, Duration duration, boolean isProbAt_IML, double level, String label,
			File outputDir, String prefix, Region zoomRegion, List<Location> annotations, boolean faults,
			TestScenario scenario) throws IOException, GMT_MapException {
		plotMap(type, null, duration, isProbAt_IML, level, label, outputDir, prefix, zoomRegion, annotations, faults, scenario);
	}
	
	public void plotMap(MapType type1, MapType type2, Duration duration, boolean isProbAt_IML, double level, String label,
			File outputDir, String prefix, Region zoomRegion, List<Location> annotations, boolean faults,
			TestScenario scenario) throws IOException, GMT_MapException {
		GriddedGeoDataSet data = calcMap(type1, type2, duration, isProbAt_IML, level);
		String printDescription = type1.name();
		String typeFileName = type1.fileName;
		if (type2 != null) {
			printDescription = "Gain "+type1.name()+"/"+type2.name();
			typeFileName = type1.fileName+"_"+type2.fileName;
			// data data is combined, now get just the denomenator and compute gain
			GriddedGeoDataSet data2 = calcMap(type2, null, getLongTermCompatibleDuration(duration), isProbAt_IML, level);
			for (int i=0; i<data.size(); i++)
				data.set(i, data.get(i)/data2.get(i));
		}
		if (isProbAt_IML)
			System.out.println("Generating map for poe="+level+", "+printDescription+", "+duration.plotName);
		else
			System.out.println("Generating map for p="+level+", "+printDescription+", "+duration.plotName);
		System.out.println("Map range: "+data.getMinZ()+" "+data.getMaxZ());
		
		String fullPrefix = prefix+"_"+typeFileName+"_"+duration.fileName;
		
		if (type2 == null && saveXYZdurations.contains(duration))
			GriddedGeoDataSet.writeXYZFile(data, new File(outputDir, fullPrefix+".xyz"));
		
		CPT cpt;
		if (type2 != null) {
			data.log10();
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, 4);
//			cpt = GMT_CPT_Files.UCERF3_ETAS_GAIN.instance(); // from 0 to 1, we want to rescale the portion from 0.5 to 1
//			for (int i=cpt.size(); --i>=0;)
//				if (cpt.get(i).start <= 0.5f)
//					cpt.remove(i);
//			cpt = cpt.rescale(0d, 4d);
//			cpt.setBelowMinColor(cpt.getMinColor());
			
//			cpt = GMT_CPT_Files.GMT_POLAR.instance().rescale(0, 2d);
			label = "Log10 Gain ETAS/"+type2.name()+" "+label;
		} else {
			if (isProbAt_IML) {
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-6, 0);
				data.log10();
				label = "Log10("+label+")";
			} else if (prefix.toLowerCase().contains("_mmi")) {
				cpt = GMT_CPT_Files.SHAKEMAP.instance();
			} else {
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
				if (Double.isInfinite(data.getMaxZ()))
					cpt = cpt.rescale(0d, 1d); // no data
				else
					cpt = cpt.rescale(0d, data.getMaxZ());
			}
		}
		
		label = duration.plotName+", "+label;
		
		doPlotMap(isProbAt_IML, label, outputDir, zoomRegion, annotations, faults, scenario, data, fullPrefix,
				cpt);
	}
	
	public GriddedGeoDataSet calcPOEShakeMap(TestScenario scenario, Duration duration, MapType type) {
		GriddedGeoDataSet shakemap = calcShakeMap(scenario);
		GriddedGeoDataSet map = new GriddedGeoDataSet(region, false);
		
		for (int i=0; i<map.size(); i++) {
			DiscretizedFunc curve = curves.get(duration, type)[i];
			double shakemapVal = shakemap.get(i);
			
			double prob = HazardDataSetLoader.getCurveVal(curve, true, shakemapVal);
			Preconditions.checkState(prob >= 0d && prob <= 1d, "Bad probability (%s) for level (%s)", prob, shakemapVal);
			if (Double.isInfinite(prob))
				prob = Double.NaN;
			
			map.set(i, prob);
		}
		
		return map;
	}
	
	public void plotPOEShakeMap(TestScenario scenario, Duration duration, MapType type, File outputDir, Region zoomRegion,
			List<Location> annotations, boolean faults) throws IOException, GMT_MapException {
		GriddedGeoDataSet data = calcPOEShakeMap(scenario, duration, type);
		
		// write out CSV for sites within 200km
		CSVFile<String> csv = new CSVFile<String>(true);
		csv.addLine("Index", "Latitude", "Longitude", "Distance (km)", "ShakeMap Level", "POE ShakeMap Level");
		GriddedGeoDataSet shakemap = calcShakeMap(scenario);
		for (int i=0; i<data.size(); i++) {
			Location loc = data.getLocation(i);
			double shakemapVal = shakemap.get(i);
			if (Double.isNaN(shakemapVal)) {
				if (imtName.equals(MMI_Param.NAME))
					shakemapVal = 1d;
				else
					shakemapVal = 0d;
			}
			csv.addLine(i+"", loc.getLatitude()+"", loc.getLongitude()+"", siteDistances.get(i)+"",
					shakemapVal+"", data.get(i)+"");
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-4, 0);
		data.log10();
		boolean customTickInterval = true;
		
		String imtLabel = imtName;
		if (imtName.equals(PGA_Param.NAME))
			imtLabel += " (g)";
		else if (imtName.equals(SA_Param.NAME))
			imtLabel = (float)period+"s "+imtLabel+" (g)";
		else if (imtName.equals(PGV_Param.NAME))
			imtLabel += " (cm/s)";
		
		String label = "POE ShakeMap Level, "+duration.plotName+", "+imtLabel;
		
		String prefix = "poe_scenario_shakemap_"+imtName.toLowerCase();
		if (imtName.equalsIgnoreCase(SA_Param.NAME))
			prefix += "_"+(float)period;
		prefix += "_"+duration.fileName;
		
		csv.writeToFile(new File(outputDir, prefix+".csv"));
		
		doPlotMap(customTickInterval, label, outputDir, zoomRegion, annotations, faults, scenario, data, prefix, cpt);
	}
	
	private GriddedGeoDataSet shakemap = null;
	protected List<Double> siteDistances = null;
	
	public synchronized GriddedGeoDataSet calcShakeMap(TestScenario scenario) {
		if (shakemap != null)
			return shakemap;
		shakemap = new GriddedGeoDataSet(region, false);
		
		ScalarIMR gmpe = checkOutGMPE();
		
		EqkRupture rup;
		if (scenario.getFSS_Index() >= 0) {
			RuptureSurface surf = sol.getRupSet().getSurfaceForRupupture(scenario.getFSS_Index(), 1d, false);
			double mag = scenario.getMagnitude();
			if (Double.isNaN(mag))
				mag = sol.getRupSet().getMagForRup(scenario.getFSS_Index());
			rup = new EqkRupture(mag, sol.getRupSet().getAveRakeForRup(scenario.getFSS_Index()), surf, null);
		} else {
			rup = new EqkRupture(scenario.getMagnitude(), 0d, new PointSurface(scenario.getLocation()), scenario.getLocation());
		}
		
		gmpe.setEqkRupture(rup);
		siteDistances = Lists.newArrayList();
		
		for (int i=0; i<sites.size(); i++) {
			gmpe.setSite(sites.get(i));
			Location loc = sites.get(i).getLocation();
			siteDistances.add(rup.getRuptureSurface().getDistanceRup(loc));
			
			double val = Math.exp(gmpe.getMean());
			
			shakemap.set(i, val);
		}
		
		checkInGMPE(gmpe);
		
		return shakemap;
	}
	
	public void plotShakeMap(TestScenario scenario, File outputDir) throws IOException, GMT_MapException {
		GriddedGeoDataSet data = calcShakeMap(scenario);
		
		CPT cpt;
		boolean customTickInterval = false;
		if (imtName.equals(MMI_Param.NAME)) {
			cpt = GMT_CPT_Files.SHAKEMAP.instance();
			customTickInterval = true;
		} else {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, data.getMaxZ());
		}
		
		String imtLabel = imtName;
		if (imtName.equals(PGA_Param.NAME))
			imtLabel += " (g)";
		else if (imtName.equals(SA_Param.NAME))
			imtLabel = (float)period+"s "+imtLabel+" (g)";
		else if (imtName.equals(PGV_Param.NAME))
			imtLabel += " (cm/s)";
		
		String label = scenario.name()+" ShakeMap, "+imtLabel;
		
		String prefix = "scenario_shakemap_"+imtName.toLowerCase();
		if (imtName.equalsIgnoreCase(SA_Param.NAME))
			prefix += "_"+(float)period;
		
		doPlotMap(customTickInterval, label, outputDir, null, null, true, scenario, data, prefix, cpt);
	}

	private void doPlotMap(boolean customTickInterval, String label, File outputDir, Region zoomRegion,
			List<Location> annotations, boolean faults, TestScenario scenario, GriddedGeoDataSet data,
			String fullPrefix, CPT cpt) throws GMT_MapException, IOException {
		cpt.setNanColor(Color.WHITE);
		
		Region region;
		if (zoomRegion != null)
			region = zoomRegion;
		else
			region = this.region;
		
		GMT_Map map = new GMT_Map(region, data, this.region.getSpacing(), cpt);
		
		map.setLogPlot(false);
//		map.setTopoResolution(TopographicSlopeFile.CA_THREE);
		map.setTopoResolution(null);
		map.setUseGMTSmoothing(false);
		map.setBlackBackground(false);
		map.setCustomScaleMin((double)cpt.getMinValue());
		map.setCustomScaleMax((double)cpt.getMaxValue());
		map.setCustomLabel(label);
		map.setRescaleCPT(false);
		map.setJPGFileName(null);
		if (customTickInterval)
			map.setCPTCustomInterval(1d);
		
		float thickness;
		if (zoomRegion == null)
			thickness = 0.5f;
		else
			thickness = 1f;
		
		if (faults) {
			Preconditions.checkNotNull(sol, "Must have FSS if you want faults");
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			for (int s=0; s<rupSet.getNumSections(); s++) {
				FaultSectionPrefData sect = rupSet.getFaultSectionData(s);
				for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), Color.BLACK, thickness))
					map.addPolys(poly);
			}
		}
		if (scenario != null && scenario.getFSS_Index() >= 0) {
			Preconditions.checkNotNull(sol, "Must have FSS if you want faults");
			FaultSystemRupSet rupSet = sol.getRupSet();
			for (FaultSectionPrefData s : rupSet.getFaultSectionDataForRupture(scenario.getFSS_Index())) {
				for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(s.getFaultTrace(), Color.WHITE, 2f*thickness))
					map.addPolys(poly);
			}
		}
		
		if (annotations != null) {
			for (Location loc : annotations) {
				java.awt.geom.Point2D.Double pt = new Point2D.Double(loc.getLongitude(), loc.getLatitude());
				map.addSymbol(new PSXYSymbol(pt, Symbol.INVERTED_TRIANGLE, 0.1f, 0f, null, Color.BLACK));
			}
		}
		
		FaultBasedMapGen.plotMap(outputDir, fullPrefix, false, map);
	}
	
	private Table<Duration, MapType, DiscretizedFunc> getInitializedCurvesMap(DiscretizedFunc xVals, double initialVal) {
		Table<Duration, MapType, DiscretizedFunc> curves = HashBasedTable.create();
		
		if (calcFaults)
			for (Duration duration : durations)
				curves.put(duration, MapType.FAULT_ONLY, getInitializeClone(xVals, initialVal));
		
		if (calcGridded)
			for (Duration duration : durations)
				curves.put(duration, MapType.GRIDDED_ONLY, getInitializeClone(xVals, initialVal));
		
		if (calcFaults && calcGridded)
			for (Duration duration: durations)
				curves.put(duration, MapType.COMBINED, getInitializeClone(xVals, initialVal));
		
		return curves;
	}
	
	private static DiscretizedFunc getInitializeClone(DiscretizedFunc curve, double val) {
		curve = curve.deepClone();
		initializeCurve(curve, val);
		return curve;
	}
	
	private static void initializeCurve(DiscretizedFunc curve, double val) {
		for (int i=0; i<curve.size(); i++)
			curve.set(i, val);
	}
	
	private static Table<Duration, MapType, File> detectCurveFiles(File dir, String prefix) {
		Table<Duration, MapType, File> curveFiles = HashBasedTable.create();
		
		Preconditions.checkState(dir.exists() && dir.isDirectory());
		
		for (File file : dir.listFiles()) {
			String name = file.getName();
			
			if (!name.startsWith(prefix) || !name.endsWith(".bin"))
				continue;
			
			name = name.substring(prefix.length());
			
//			System.out.println("Processing: "+name);
			
			// detect duration
			Duration duration = durationForFileName(name);
//			System.out.println("Detected duration: "+duration);
			if (duration == null)
				// assume full
				duration = DurationConstants.FULL.duration;
			
			// detect type
			MapType type = null;
			for (MapType t : MapType.values()) {
				if (name.contains(t.fileName)) {
					type = t;
					break;
				}
			}
			Preconditions.checkNotNull(type, "Couldn't detect type from file: %s", file.getName());
			
			curveFiles.put(duration, type, file);
		}
		
		Preconditions.checkState(!curveFiles.isEmpty(), "No matching curve files with prefix: %s", prefix);
		
		return curveFiles;
	}
	
	static List<Site> loadSitesFile(AttenRelRef gmpeRef, File sitesFile) {
		Document siteDoc;
		try {
			siteDoc = XMLUtils.loadDocument(sitesFile);
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		Element sitesRoot = siteDoc.getRootElement();
		ArrayList<Parameter<?>> paramsToAdd = new ArrayList<Parameter<?>>();
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		for (Parameter<?> param : gmpe.getSiteParams())
			paramsToAdd.add(param);
		return Site.loadSitesFromXML(sitesRoot.element(Site.XML_METADATA_LIST_NAME), paramsToAdd);
	}
	
	private static class ETAS_MMI_HazardMapCalc extends ETAS_HazardMapCalc {
		
		private ETAS_HazardMapCalc pgaCalc, pgvCalc;
		
		public ETAS_MMI_HazardMapCalc(ETAS_HazardMapCalc pgaCalc, ETAS_HazardMapCalc pgvCalc) {
			this.pgaCalc = pgaCalc;
			this.pgvCalc = pgvCalc;
			this.imtName = MMI_Param.NAME;
			
			Preconditions.checkState(pgaCalc.region.equals(pgvCalc.region));
			this.region = pgaCalc.region;
			
			// fill in curves map with empty curve arrays, just so that we know what types/durations are available
			curves = HashBasedTable.create();
			for (Cell<Duration, MapType, DiscretizedFunc[]> cell : pgaCalc.curves.cellSet()) {
				Duration duration = cell.getRowKey();
				MapType type = cell.getColumnKey();
				if (pgvCalc.curves.contains(duration, type))
					curves.put(duration, type, new DiscretizedFunc[0]);
			}
			
			this.durations = pgaCalc.durations;
			if (pgaCalc.isCalcLongTerm()) {
				setCalcLongTerm(true);
				setLongTermCalcDurations(pgaCalc.getLongTermCalcDurations());
			}
		}

		@Override
		public GriddedGeoDataSet calcMap(MapType type1, MapType type2, Duration duration, boolean isProbAt_IML, double level) {
			if (isProbAt_IML) {
				double pga = Wald_MMI_Calc.getPGA(level);
				double pgv = Wald_MMI_Calc.getPGV(level);
				double weightPGV = Wald_MMI_Calc.getWeightVMMI(level);
				
				GriddedGeoDataSet pgaMap = pgaCalc.calcMap(type1, type2, duration, true, pga);
				GriddedGeoDataSet pgvMap = pgvCalc.calcMap(type1, type2, duration, true, pgv);
				
				GriddedGeoDataSet mmiMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
				for (int index=0; index<mmiMap.size(); index++) {
					double pgaProb = pgaMap.get(index);
					double pgvProb = pgvMap.get(index);
					double mmiProb = calcMMI_POE(weightPGV, pgaProb, pgvProb);
					mmiMap.set(index, mmiProb);
				}
				return mmiMap;
			} else {
				GriddedGeoDataSet pgaMap = pgaCalc.calcMap(type1, type2, duration, false, level);
				GriddedGeoDataSet pgvMap = pgvCalc.calcMap(type1, type2, duration, false, level);
				
				GriddedGeoDataSet mmiMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
				
				for (int index=0; index<mmiMap.size(); index++) {
					double pga = pgaMap.get(index);
					if (!Doubles.isFinite(pga))
						pga = 0;
					double pgv = pgvMap.get(index);
					if (!Doubles.isFinite(pgv))
						pgv = 0;
					double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
					Preconditions.checkState(Doubles.isFinite(mmi), "Bad MMI=%s for PGA=%s, PGV=%s", mmi, pga, pgv);
					if (mmi == 1d)
						// will speed things up
						mmi = Double.NaN;
					mmiMap.set(index, mmi);
				}
				return mmiMap;
			}
		}

		private double calcMMI_POE(double weightPGV, double pgaProb, double pgvProb) {
			if (!Doubles.isFinite(pgaProb))
				pgaProb = 0;
			Preconditions.checkState((float)pgaProb <= 1f, "Bad PGA prob: %s", pgaProb);
			
			if (!Doubles.isFinite(pgvProb))
				pgvProb = 0;
			Preconditions.checkState((float)pgvProb <= 1f, "Bad PGV prob: %s", pgvProb);
			
			double mmiProb = weightPGV*pgvProb + (1d-weightPGV)*pgaProb;
			Preconditions.checkState(mmiProb >= 0d && mmiProb <= 1d, "Bad MMI probability (%s) for pgaProb= %s and pgvProb=%s",
					mmiProb, pgaProb, pgvProb);
			if (mmiProb == 0d)
				// will speed things up
				mmiProb = Double.NaN;
			return mmiProb;
		}
		
		@Override
		public GriddedGeoDataSet calcPOEShakeMap(TestScenario scenario, Duration duration, MapType type) {
			GriddedGeoDataSet shakemap = calcShakeMap(scenario);
			GriddedGeoDataSet map = new GriddedGeoDataSet(region, false);
			
			for (int i=0; i<map.size(); i++) {
				DiscretizedFunc pgaCurve = pgaCalc.curves.get(duration, type)[i];
				DiscretizedFunc pgvCurve = pgvCalc.curves.get(duration, type)[i];
				double level = shakemap.get(i);
				if (Double.isNaN(level) || level < 1d)
					level = 1d;
				
				double pga = Wald_MMI_Calc.getPGA(level);
				double pgv = Wald_MMI_Calc.getPGV(level);
				double weightPGV = Wald_MMI_Calc.getWeightVMMI(level);
				
				double pgaProb = HazardDataSetLoader.getCurveVal(pgaCurve, true, pga);
				double pgvProb = HazardDataSetLoader.getCurveVal(pgvCurve, true, pgv);
				
				double prob = calcMMI_POE(weightPGV, pgaProb, pgvProb);
				if (Double.isInfinite(prob))
					prob = Double.NaN;
				
				map.set(i, prob);
			}
			
			return map;
		}
		
		private GriddedGeoDataSet mmiShakeMap = null;

		@Override
		public synchronized GriddedGeoDataSet calcShakeMap(TestScenario scenario) {
			if (mmiShakeMap != null)
				return mmiShakeMap;
			GriddedGeoDataSet pgaMap = pgaCalc.calcShakeMap(scenario);
			GriddedGeoDataSet pgvMap = pgvCalc.calcShakeMap(scenario);
			mmiShakeMap = new GriddedGeoDataSet(pgaMap.getRegion(), pgaMap.isLatitudeX());
			siteDistances = pgaCalc.siteDistances;
			
			for (int i=0; i<pgaMap.size(); i++) {
				double mmi = Wald_MMI_Calc.getMMI(pgaMap.get(i), pgvMap.get(i));
				if (mmi <= 1d)
					mmi = Double.NaN;
				mmiShakeMap.set(i, mmi);
			}
			
			return mmiShakeMap;
		}
		
	}
	
	private static class FakeGriddedRegion extends GriddedRegion {
		
		private LocationList locs;
		
		public FakeGriddedRegion(Region region, List<Location> locs) {
			super(region, 1d, null);
			this.locs = new LocationList();
			this.locs.addAll(locs);
		}

		@Override
		public int getNodeCount() {
			return locs.size();
		}

		@Override
		public LocationList getNodeList() {
			return locs;
		}

		@Override
		public Location locationForIndex(int index) {
			return locs.get(index);
		}

		@Override
		public Location getLocation(int index) {
			return locs.get(index);
		}

		@Override
		public int indexForLocation(Location loc) {
			return locs.indexOf(loc);
		}
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws Exception {
		force_serial = false;
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		boolean calcFault = true;
		boolean calcGridded = true;
		boolean calcLongTerm = true;
		boolean mapParallel = true;
		boolean plotCurves = false;
		boolean plotMaps = false;
		boolean plotLongTerm = false;
		boolean onlyGainLongTerm = false;
		boolean plotShakeMap = true;
		boolean plotShakeMapPOE = true;
		HashSet<MapType> mapTypePlotSubset = new HashSet<MapType>(Lists.newArrayList(MapType.COMBINED));
//		HashSet<MapType> mapTypePlotSubset = new HashSet<MapType>(Lists.newArrayList(MapType.COMBINED, MapType.U3TD, MapType.U3TI));
//		HashSet<MapType> mapTypePlotSubset = null;
		HashSet<Duration> durationPlotSubset = null;
//		HashSet<Duration> durationPlotSubset = new HashSet<Duration>(Lists.newArrayList(DurationConstants.THREE.duration));
//		calc.debugCurvePlotModulus = 10;
//		calc.debugStopIndex = 500;
		
//		File faultBasedPrecalc = null;
////		double spacing = 0.02d;
////		File precalcDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
////				+ "2017_02_28-mojave_m7_fulltd_descendents-NGA2-0.02-site-effects-with-basin");
//		double spacing = 0.01d;
//		File precalcDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
//				+ "2017_03_01-mojave_m7_fulltd_descendents-NGA2-0.01-site-effects-with-basin");
		
//		File faultBasedPrecalc = null;
//		double spacing = 1.0;
//		File precalcDir = new File("/tmp/etas_hazard_test");
		
		File faultBasedPrecalc = null;
		double spacing = 0.02;
		File precalcDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
				// orig set
//				+ "2017_03_03-haywired_m7_fulltd_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_20-haywired_m7_combined_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_03-haywired_m7_gridded_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_20-northridge_combined_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_20-northridge_gridded_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_21-mojave_m7_combined_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_22-mojave_m7_gridded_descendents-NGA2-0.02-site-effects-with-basin");
				// final Powell set
//				+ "2017_03_23-haywired_m7_combined_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_23-haywired_m7_gridded_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_23-mojave_m7_combined_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_23-mojave_m7_gridded_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_23-northridge_combined_descendents-NGA2-0.02-site-effects-with-basin");
				+ "2017_03_23-northridge_gridded_descendents-NGA2-0.02-site-effects-with-basin");
//				+ "2017_03_23-2016_bombay_swarm_combined_descendents-NGA2-0.02-site-effects-with-basin");
				// other
//				+ "2017_05_17-2017_05-usgs_exercise-1pm-NGA2-0.02-site-effects-with-basin");
		
//		File faultBasedPrecalc = null;
//		double spacing = 0.5;
//		File precalcDir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
//				+ "2017_03_23-haywired_m7_combined_descendents-NGA2-0.02-site-effects-with-basin/local_test");
		
//		File faultBasedPrecalc = null;
//		double spacing = 0.5;
//		File precalcDir = null; 
		
		TestScenario scenarioForRebound = null;
//		String etasDirName = "2016_06_15-haywired_m7-10yr-full_td-no_ert-combined"; scenarioForRebound = TestScenario.HAYWIRED_M7;
//		String etasDirName = "2017_01_02-haywired_m7-10yr-gridded-only-200kcombined"; scenarioForRebound = TestScenario.HAYWIRED_M7;
//		String etasDirName = "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined"; scenarioForRebound = TestScenario.MOJAVE_M7;
//		String etasDirName = "2016_12_03-mojave_m7-10yr-gridded-only"; scenarioForRebound = TestScenario.MOJAVE_M7;
//		String etasDirName = "2017_02_01-northridge-m6.7-10yr-full_td-no_ert-combined"; scenarioForRebound = TestScenario.NORTHRIDGE;
		String etasDirName = "2017_02_01-northridge-m6.7-10yr-gridded-only-combined200k"; scenarioForRebound = TestScenario.NORTHRIDGE;
//		String etasDirName = "2016_10_27-2016_bombay_swarm-10yr-full_td-no_ert-combined";
//		String etasDirName = "2017_05_17-USGS_Exercise_Catalog-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-noSpont-sect-reset-1pm";
//		String etasFileName = "results_descendents_m5_preserve.bin";
//		String etasFileName = "results_descendents_m5.bin";
//		String etasFileName = "results_m5_preserve.bin";
		String etasFileName = "results_m5.bin";
		double etasSimDuration = 10d;
		
		Region plotRegion = null;
		List<Location> annotations = null;
//		Region plotRegion = new Region(new Location(35, -119), new Location(33, -115));
//		List<Location> annotations = Lists.newArrayList();
//		annotations.add(new Location(33.739683, -116.412925));
		boolean faults = true;
		boolean highlightScenario = true && scenarioForRebound != null;
		
		boolean mapProbs = false;
		boolean mapPOE = true;
		
//		String imtName = PGA_Param.NAME;
//		double period = Double.NaN;
//		String imtLabel = "PGA";
//		String imtFileLabel = "pga";
//		String imtName = PGV_Param.NAME;
//		double period = Double.NaN;
//		String imtLabel = "PGV";
//		String imtFileLabel = "pgv";
//		String imtName = SA_Param.NAME;
//		double period = 1d;
//		String imtLabel = "1s Sa";
//		String imtFileLabel = "sa_1s";
		String imtName = MMI_Param.NAME;
		double period = Double.NaN;
		String imtLabel = "MMI";
		String imtFileLabel = "mmi";
		GriddedRegion region = new CaliforniaRegions.RELM_TESTING_GRIDDED(spacing);
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
		AttenRelRef gmpeRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
		ArrayList<SiteData<?>> provs = null;
		double griddedResolution = 0.01;
		boolean griddedConditional = true;
		
		List<Location> curveLocs = null;
		List<String> curveSiteNames = null;
		if (imtName.equals(MMI_Param.NAME)) {
			plotCurves = false;
			Preconditions.checkState(precalcDir != null);
		}
		
		if (plotCurves) {
			curveLocs = Lists.newArrayList();
			curveSiteNames = Lists.newArrayList();
			
			if (scenarioForRebound == TestScenario.HAYWIRED_M7) {
				curveLocs.add(NEHRP_TestCity.OAKLAND.location());
				curveSiteNames.add("Oakland");
				curveLocs.add(NEHRP_TestCity.SAN_FRANCISCO.location());
				curveSiteNames.add("San Francisco");
			} else {
				curveLocs.add(NEHRP_TestCity.LOS_ANGELES.location());
				curveSiteNames.add("Los Angeles");
				curveLocs.add(NEHRP_TestCity.SAN_BERNARDINO.location());
				curveSiteNames.add("San Bernardino");
				if (scenarioForRebound != TestScenario.NORTHRIDGE) {
					curveLocs.add(new Location(33.8303, -116.5453));
					curveSiteNames.add("Palm Springs");
				}
			} 
			
		}
		
		File etasCatalogs = new File(new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/", etasDirName), etasFileName);
		
		if (plotCurves && !plotMaps && !curveLocs.isEmpty() && precalcDir == null) {
			// only calculate for sites of interest
			System.out.println("Creating fake gridded region for "+curveLocs.size()+" locations");
			region = new FakeGriddedRegion(region, curveLocs);
		}
		
		File outputDir = new File(etasCatalogs.getParentFile(), "hazard_maps");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		if (faultBasedPrecalc == null && precalcDir == null) {
			if (region instanceof FakeGriddedRegion) {
				
			} else {
				outputDir = new File(outputDir, "testing");
			}
		} else if (faultBasedPrecalc == null) {
			outputDir = new File(outputDir, precalcDir.getName());
		} else {
			outputDir = new File(outputDir, faultBasedPrecalc.getParentFile().getName());
		}
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		
		if (!calcFault)
			faultBasedPrecalc = null;
		
		FaultSystemSolution sol = null;
		ETAS_HazardMapCalc calc;
		if (precalcDir != null) {
			if (mapTypePlotSubset != null && precalcDir.getName().contains("gridded")
					&& !mapTypePlotSubset.contains(MapType.GRIDDED_ONLY))
				mapTypePlotSubset.add(MapType.GRIDDED_ONLY);
			
			List<Site> sites = null;
			if (plotCurves || plotShakeMap) {
				// load in the sites, we'll need them
				File sitesFile = new File(precalcDir, "sites.xml");
				Preconditions.checkState(sitesFile.exists(), "Need sites, but no sites.xml in %s", precalcDir.getAbsolutePath());
				sites = loadSitesFile(gmpeRef, sitesFile);
			}
			
			if (imtName.equals(MMI_Param.NAME)) {
				ETAS_HazardMapCalc pgaCalc = new ETAS_HazardMapCalc(region, detectCurveFiles(precalcDir, "results_pga"),
						xVals, gmpeRef, PGA_Param.NAME, period, sites);
				ETAS_HazardMapCalc pgvCalc = new ETAS_HazardMapCalc(region, detectCurveFiles(precalcDir, "results_pgv"),
						xVals, gmpeRef, PGV_Param.NAME, period, sites);
				calc = new ETAS_MMI_HazardMapCalc(pgaCalc, pgvCalc);
				if (plotShakeMap) {
					sol = FaultSystemIO.loadSol(solFile);
					pgaCalc.setSol(sol);
					pgvCalc.setSol(sol);
					calc.setSol(sol);
				}
			} else {
				Table<Duration, MapType, File> curveFiles = detectCurveFiles(precalcDir, "results_"+imtFileLabel);
				calc = new ETAS_HazardMapCalc(region, curveFiles, xVals, gmpeRef, imtName, period, sites);
			}
		} else {
			calcGridded = calcGridded && gmpeRef != null;
			List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(etasCatalogs, 5d);
			
			if (plotCurves)
				calcLongTerm = true;
			
			if (calcFault && faultBasedPrecalc == null || calcLongTerm)
				sol = FaultSystemIO.loadSol(solFile);
			ETAS_CatalogGridSourceProvider gridSources = null;
			if (calcGridded)
				gridSources = new ETAS_CatalogGridSourceProvider(catalogs, griddedResolution, griddedConditional);
			List<Site> sites = null;
			if (calcGridded || sol != null) {
				// need sites
				AttenuationRelationship gmpe = gmpeRef.instance(null);
				gmpe.setParamDefaults();
				sites = fetchSites(region, provs, gmpe);
			}
			Duration[] durations;
			if (durationPlotSubset != null) {
				List<Duration> durLists = Lists.newArrayList(durationPlotSubset);
//				for (DurationConstants durConst : durationPlotSubset)
//					durLists.add(durConst.duration);
				durations = Lists.newArrayList(durLists).toArray(new Duration[0]);
			} else {
				durations = getDurationDefaults();
			}
			calc = new ETAS_HazardMapCalc(catalogs, region, xVals, faultBasedPrecalc, sol, gridSources,
					gmpeRef, imtName, period, sites, durations);
			if (!calcFault)
				calc.setCalcFaults(false);
			if (!calcGridded)
				calc.setCalcGridded(false);
			calc.setCalcLongTerm(calcLongTerm);
			if (scenarioForRebound != null)
				calc.setScenarioForElasticRebound(scenarioForRebound);
			
			Preconditions.checkState(calcFault || calcGridded);
			calc.calculate();
		}
		
//		for (int i=0; i<region.getNodeCount(); i+= 500)
//			calc.plotCurve(i);
		
		Set<MapType> types = calc.curves.columnKeySet();
		
		if (plotCurves) {
			File curveDir = new File(outputDir, "curves");
			Preconditions.checkState(curveDir.exists() || curveDir.mkdir());
			
			MapType type;
			if (types.contains(MapType.COMBINED))
				type = MapType.COMBINED;
			else
				type = types.iterator().next();
			for (int i=0; i<curveLocs.size(); i++) {
				Location loc = curveLocs.get(i);
				String name = curveSiteNames.get(i);
				System.out.println("Plotting curves for "+name);
				String prefix = name.toLowerCase().replaceAll(" ", "_")+"_"+imtFileLabel;
				File siteDir = new File(curveDir, prefix);
				Preconditions.checkState(siteDir.exists() || siteDir.mkdir());
				
				for (Duration duration : calc.durationsForType(type))
					calc.plotComparisons(type, duration, loc, etasSimDuration, name, siteDir, prefix);
			}
		}
		
		if (plotShakeMap) {
			calc.plotShakeMap(scenarioForRebound, outputDir);
		}
		
		if (plotShakeMapPOE) {
			MapType type;
			if (calc.curves.contains(calc.durations[0], MapType.COMBINED)) {
				type = MapType.COMBINED;
			} else {
				type = MapType.GRIDDED_ONLY;
				Preconditions.checkState(calc.curves.contains(calc.durations[0], type));
			}
			System.out.println("ShakeMap POE type: "+type);
			
			ExecutorService exec = null;
			List<Future<?>> futures = null;
			if (mapParallel) {
				if (FaultBasedMapGen.LOCAL_MAPGEN)
					exec = Executors.newFixedThreadPool(6);
				else
					exec = Executors.newFixedThreadPool(10);
				futures = Lists.newArrayList();
			}
			
			for (Duration duration : calc.durations) {
				File dir = getMapDir(outputDir, type, duration, plotRegion != null);
				
				ShakeMapPOEPlotRunnable runnable = new ShakeMapPOEPlotRunnable(type, duration, dir, calc, plotRegion,
						annotations, faults, scenarioForRebound);
				
				if (exec == null) {
					runnable.run();
				} else {
					futures.add(exec.submit(runnable));
					// sleep for just a bit to avoid dir name collisions
					Thread.sleep(1000);
				}
			}
			
			if (exec != null) {
				try {
					int num = futures.size();
					int index = 0;
					for (Future<?> future : futures) {
						future.get();
							System.out.println("Finished plot "+index+++"/"+num);
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				exec.shutdown();
			}
		}
		
		if (plotMaps) {
			Preconditions.checkState(!types.isEmpty());
			
			if ((faults || highlightScenario) && sol == null) {
				sol = FaultSystemIO.loadSol(solFile);
				calc.setSol(sol);
			}
			
			TestScenario scenario = null;
			if (highlightScenario)
				scenario = ETAS_MultiSimAnalysisTools.detectScenario(etasCatalogs.getParentFile());
			
			ExecutorService exec = null;
			List<Future<?>> futures = null;
			if (mapParallel) {
				if (FaultBasedMapGen.LOCAL_MAPGEN)
					exec = Executors.newFixedThreadPool(6);
				else
					exec = Executors.newFixedThreadPool(10);
				futures = Lists.newArrayList();
			}
			
			System.out.println("Plotting maps for "+Joiner.on(",").join(types));
			
			List<Double> levels = Lists.newArrayList();
			List<Boolean> isProbAtIMLs = Lists.newArrayList();
			if (mapPOE && imtName.equals(MMI_Param.NAME)) {
				levels.add(6d);
				isProbAtIMLs.add(true);
				levels.add(8d);
				isProbAtIMLs.add(true);
			}
			if (mapProbs) {
				double[] probVals = { 0.5, 0.25, 0.1d, 0.01d, 0.001 };
				for (double p : probVals) {
					levels.add(p);
					isProbAtIMLs.add(false);
				}
			}
			for (int i=0; i<levels.size(); i++) {
				double level = levels.get(i);
				boolean isProbAtIML = isProbAtIMLs.get(i);
				
				String label, prefix;
				if (isProbAtIML) {
					String levelStr;
					if (level == (int)level)
						levelStr = (int)level+"";
					else
						levelStr = (float)level+"";
					if (imtName.equals(MMI_Param.NAME))
						label = "POE "+imtLabel+" "+levelStr;
					else if (imtName.equals(PGV_Param.NAME))
						label = "POE "+levelStr+" cm/s "+imtLabel;
					else
						label = "POE "+levelStr+" G "+imtLabel;
					prefix = "map_"+imtFileLabel+"_poe_"+levelStr;
				} else {
					String probString;
					double p100 = level*100;
					if (p100 == Math.floor(p100))
						probString = (int)p100+"";
					else
						probString = (float)p100+"";
					label = imtLabel+" @@ "+probString+"% POE";
					prefix = "map_"+imtFileLabel+"_p"+(float)level;
				}
				
				List<MapType> type1s = Lists.newArrayList();
				List<MapType> type2s = Lists.newArrayList();
				for (MapType type : types) {
					if (!type.isETAS()) {
						// this is UCERF3 - plot gain instead
						MapType etas;
						if (types.contains(MapType.COMBINED)) {
							etas = MapType.COMBINED;
						} else {
							// this is gridded only
							Preconditions.checkState(types.contains(MapType.GRIDDED_ONLY));
							Preconditions.checkState(!types.contains(MapType.FAULT_ONLY));
							etas = MapType.GRIDDED_ONLY;
						}
						// add ratio map
						type1s.add(etas);
						type2s.add(type);
						// add long term only map
						if (plotLongTerm) {
							type1s.add(type);
							type2s.add(null);
						}
					} else {
						if (!onlyGainLongTerm) {
							type1s.add(type);
							type2s.add(null);
						}
					}
				}
				
				for (int t=0; t<type1s.size(); t++) {
					MapType type1 = type1s.get(t);
					MapType type2 = type2s.get(t);
					if (mapTypePlotSubset != null && !mapTypePlotSubset.contains(type1))
						continue;
					
					for (Duration duration : calc.durationsForType(type1)) {
						if (durationPlotSubset != null && !durationPlotSubset.contains(duration))
							continue;
						if (type2 != null && !calc.curves.contains(calc.getLongTermCompatibleDuration(duration), type2))
							// ratio doesn't exist for this duration
							continue;
						File durationDir = getMapDir(outputDir, type1, duration, plotRegion != null);
//						System.out.println("label: "+label+", "+type+", "+duration);
						String myPrefix = prefix;
						if (type2 != null)
							myPrefix = "gain_"+prefix;
						MapPlotRunnable runnable = new MapPlotRunnable(isProbAtIML, level, type1, type2, duration,
								label, durationDir, myPrefix, calc, plotRegion, annotations, faults, scenario);
						if (exec == null) {
							runnable.run();
						} else {
							futures.add(exec.submit(runnable));
							// sleep for just a bit to avoid dir name collisions
							Thread.sleep(1000);
						}
					}
				}
			}
			
			if (exec != null) {
				try {
					int num = futures.size();
					int index = 0;
					for (Future<?> future : futures) {
						future.get();
							System.out.println("Finished plot "+index+++"/"+num);
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				exec.shutdown();
			}
		}
	}
	
	private static File getMapDir(File outputDir, MapType type, Duration duration, boolean zoomed) {
		File typeDir = new File(outputDir, "maps_"+type.fileName);
		if (zoomed)
			typeDir = new File(typeDir.getAbsolutePath()+"_zoomed");
		Preconditions.checkState(typeDir.exists() || typeDir.mkdir());
		File durationDir = new File(typeDir, duration.fileName);
		Preconditions.checkState(durationDir.exists() || durationDir.mkdir());
		return durationDir;
	}
	
	private static class MapPlotRunnable implements Runnable {
		
		private boolean isProbAtIML;
		private double level;
		private MapType type1;
		private MapType type2;
		private Duration duration;
		private String label;
		private File outputDir;
		private String prefix;
		
		private ETAS_HazardMapCalc calc;
		private Region zoomRegion;
		private List<Location> annotations;
		private boolean faults;
		private TestScenario scenario;
		
		public MapPlotRunnable(boolean isProbAtIML, double level, MapType type1, MapType type2, Duration duration,
				String label, File outputDir, String prefix, ETAS_HazardMapCalc calc, Region zoomRegion,
				List<Location> annotations, boolean faults, TestScenario scenario) {
			super();
			this.isProbAtIML = isProbAtIML;
			this.level = level;
			this.type1 = type1;
			this.type2 = type2;
			this.duration = duration;
			this.label = label;
			this.outputDir = outputDir;
			this.prefix = prefix;
			this.calc = calc;
			this.zoomRegion = zoomRegion;
			this.annotations = annotations;
			this.faults = faults;
			this.scenario = scenario;
		}
		
		@Override
		public void run() {
			try {
				calc.plotMap(type1, type2, duration, isProbAtIML, level, label, outputDir, prefix, zoomRegion,
						annotations, faults, scenario);
				System.out.println("Done with "+prefix);
			} catch (Exception e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	private static class ShakeMapPOEPlotRunnable implements Runnable {
		
		private MapType type;
		private Duration duration;
		private File outputDir;
		
		private ETAS_HazardMapCalc calc;
		private Region zoomRegion;
		private List<Location> annotations;
		private boolean faults;
		private TestScenario scenario;
		
		public ShakeMapPOEPlotRunnable(MapType type, Duration duration, File outputDir, ETAS_HazardMapCalc calc,
				Region zoomRegion, List<Location> annotations, boolean faults, TestScenario scenario) {
			super();
			this.type = type;
			this.duration = duration;
			this.outputDir = outputDir;
			this.calc = calc;
			this.zoomRegion = zoomRegion;
			this.annotations = annotations;
			this.faults = faults;
			this.scenario = scenario;
		}


		@Override
		public void run() {
			try {
				calc.plotPOEShakeMap(scenario, duration, type, outputDir, zoomRegion, annotations, faults);
				System.out.println("Done with POE ShakeMap, "+duration);
			} catch (Exception e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}

}
