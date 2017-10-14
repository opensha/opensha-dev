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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class RSQSimCatalogGMPE_Compare {
	
	private RSQSimCatalog catalog;
	private BBP_CatalogSimZipLoader bbpZipFile;
	private List<BBP_Site> sites;
	private double minMag;
	private double minFractForInclusion;
	
	private Map<BBP_Site, Site> gmpeSites;
	private List<RSQSimEvent> events;
	private Map<Integer, RSQSimEvent> eventsMap;
	private Map<Integer, EqkRupture> gmpeRupsMap;
	private Map<Integer, DiscretizedFunc> rupSpectraMap;
	
	private Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesCache;
	private ExecutorService exec;
	
	private File gmpeCacheDir;

	public RSQSimCatalogGMPE_Compare(RSQSimCatalog catalog, ZipFile bbpZipFile, List<BBP_Site> sites, double minMag, int skipYears,
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
		
		gmpesCache = new HashMap<>();
		gmpeRupsMap = new HashMap<>();
		rupSpectraMap = new HashMap<>();
		
		exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	}
	
	public void shutdown() {
		exec.shutdown();
	}
	
	private EqkRupture getGMPE_Rup(RSQSimEvent event) {
		synchronized(event) {
			EqkRupture rup = gmpeRupsMap.get(event.getID());
			if (rup == null) {
				rup = catalog.getGMPE_Rupture(event, minFractForInclusion);
				gmpeRupsMap.put(event.getID(), rup);
			}
			return rup;
		}
	}
	
	private ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
		synchronized (gmpesCache) {
			LinkedList<ScalarIMR> gmpes = gmpesCache.get(gmpeRef);
			if (gmpes == null) {
				gmpes = new LinkedList<>();
				gmpesCache.put(gmpeRef, gmpes);
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
		synchronized (gmpesCache) {
			gmpesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	public DiscretizedFunc getEventSpectra(RSQSimEvent event, BBP_Site site) throws IOException {
		Integer eventID = event.getID();
		if (rupSpectraMap.containsKey(eventID))
			return rupSpectraMap.get(eventID);
		synchronized (rupSpectraMap) {
			// repeat here as could have been populated while waiting for the lock
			if (rupSpectraMap.containsKey(eventID))
				return rupSpectraMap.get(eventID);
			DiscretizedFunc spectra = bbpZipFile.readRotD50(site, event.getID());
			rupSpectraMap.put(event.getID(), spectra);
			return spectra;
		}
	}
	
	private Table<Integer, Double, EventComparison> loadCache(File file, BBP_Site site, AttenRelRef gmpeRef) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		
		Table<Integer, Double, EventComparison> ret = HashBasedTable.create();
		
		for (int row=1; row<csv.getNumRows(); row++) {
			Integer eventID = Integer.parseInt(csv.get(row, 0));
			Double period = Double.parseDouble(csv.get(row, 1));
			Double logMean = Double.parseDouble(csv.get(row, 2));
			Double stdDev = Double.parseDouble(csv.get(row, 3));
			Double rRup = Double.parseDouble(csv.get(row, 4));
			Double rJB = Double.parseDouble(csv.get(row, 5));
			
			ret.put(eventID, period, new EventComparison(eventsMap.get(eventID), gmpeRef, site, period, logMean, stdDev, rRup, rJB));
		}
		
		return ret;
	}
	
	private void writeCache(File file, Table<Integer, Double, EventComparison> table) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Event ID", "Period (s)", "Ln(Mean)", "Std Dev", "rRup (km)", "rJB (km)");
		
		for (Integer eventID : sorted(table.rowKeySet())) {
			Map<Double, EventComparison> eventVal = table.row(eventID);
			for (Double period : sorted(eventVal.keySet())) {
				EventComparison comp = table.get(eventID, period);
				csv.addLine(eventID+"", period+"", comp.logMean+"", comp.stdDev+"", comp.distanceRup+"", comp.distanceJB+"");
			}
		}
		
		csv.writeToFile(file);
	}
	
	private <E extends Comparable<E>> Iterable<E> sorted(Set<E> vals) {
		List<E> list = new ArrayList<>(vals);
		Collections.sort(list);
		return list;
	}
	
	public Table<Integer, Double, EventComparison> calcGMPE(BBP_Site site, AttenRelRef gmpeRef, double... periods) {
		List<Future<Map<Double, EventComparison>>> futures = new ArrayList<>();
		
		Table<Integer, Double, EventComparison> table = HashBasedTable.create();
		
		File gmpeCacheFile = null;
		if (gmpeCacheDir != null) {
			gmpeCacheFile = new File(gmpeCacheDir, site.getName()+"_"+gmpeRef.getShortName()+".csv");
			if (gmpeCacheFile.exists()) {
				try {
					System.out.println("Loading cache from "+gmpeCacheFile.getAbsolutePath());
					table.putAll(loadCache(gmpeCacheFile, site, gmpeRef));
					System.out.println("Loaded "+table.size()+" cached values");
				} catch (IOException e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
		}
		
		List<Double> periodsList = Doubles.asList(periods);
		
		for (int eventID : bbpZipFile.getEventIDs(site)) {
			List<Double> myPeriods;
			if (table.isEmpty()) {
				// nothing cached
				myPeriods = periodsList;
			} else {
				myPeriods = new ArrayList<>();
				for (double period : periods)
					if (!table.contains(eventID, period))
						myPeriods.add(period);
			}
			if (!myPeriods.isEmpty())
				futures.add(exec.submit(new EventComparisonCalc(eventsMap.get(eventID), gmpeRef, site, myPeriods)));
		}
		
		System.out.println("Calculating for "+futures.size()+" events");
		
		for (Future<Map<Double, EventComparison>> future : futures) {
			try {
				Map<Double, EventComparison> calcMap = future.get();
				for (Double period : calcMap.keySet()) {
					EventComparison calc = calcMap.get(period);
					table.put(calc.event.getID(), period, calc);
				}
			} catch (Exception e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		if (gmpeCacheFile != null && !futures.isEmpty()) {
			System.out.println("Caching "+futures.size()+" new vals to "+gmpeCacheFile.getAbsolutePath());
			try {
				writeCache(gmpeCacheFile, table);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		return table;
	}
	
	private class EventComparisonCalc implements Callable<Map<Double, EventComparison>> {
		
		private RSQSimEvent event;
		private AttenRelRef gmpeRef;
		private BBP_Site site;
		private List<Double> periods;
		
		public EventComparisonCalc(RSQSimEvent event, AttenRelRef gmpeRef, BBP_Site site, List<Double> periods) {
			this.event = event;
			this.gmpeRef = gmpeRef;
			this.site = site;
			this.periods = periods;
		}

		@Override
		public Map<Double, EventComparison> call() throws Exception {
			ScalarIMR gmpe = checkOutGMPE(gmpeRef);
			
			SA_Param imt = (SA_Param)gmpe.getIntensityMeasure();
			SA_Param.setPeriodInSA_Param(imt, periods.get(0));
			EqkRupture gmpeRup = getGMPE_Rup(event);
			gmpe.setAll(gmpeRup, gmpeSites.get(site), imt);
			
			Map<Double, EventComparison> comps = new HashMap<>();
			
			RuptureSurface surf = gmpeRup.getRuptureSurface();
			double distanceRup = surf.getDistanceRup(site.getLoc());
			double distanceJB = surf.getDistanceJB(site.getLoc());
			
			for (double period : periods) {
				SA_Param.setPeriodInSA_Param(imt, period);
				double logMean = gmpe.getMean();
				double stdDev = gmpe.getStdDev();
				comps.put(period, new EventComparison(event, gmpeRef, site, period, logMean, stdDev, distanceRup, distanceJB));
			}
			
			checkInGMPE(gmpeRef, gmpe);
			return comps;
		}
		
	}
	
	private class EventComparison {
		private final RSQSimEvent event;
		private final AttenRelRef gmpeRef;
		private final BBP_Site site;
		private final double period;
		
		private final double logMean;
		private final double stdDev;
		private final double distanceRup;
		private final double distanceJB;
		
		public EventComparison(RSQSimEvent event, AttenRelRef gmpeRef, BBP_Site site, double period,
				double logMean, double stdDev, double distanceRup, double distanceJB) {
			this.event = event;
			this.gmpeRef = gmpeRef;
			this.site = site;
			this.period = period;
			this.logMean = logMean;
			this.stdDev = stdDev;
			this.distanceRup = distanceRup;
			this.distanceJB = distanceJB;
		}
	}
	
	public List<RSQSimEvent> getEventsForSite(BBP_Site site) {
		return eventsForIndexes(bbpZipFile.getEventIDs(site));
	}
	
	public List<RSQSimEvent> eventsForIndexes(Collection<Integer> indexes) {
		List<RSQSimEvent> events = new ArrayList<>(indexes.size());
		for (Integer index : indexes)
			events.add(eventsMap.get(index));
		return events;
	}
	
	public List<RSQSimEvent> getWithinDistRup(List<RSQSimEvent> events, Table<Integer, Double, EventComparison> gmpeVals,
			double minDist, double maxDist) {
		List<RSQSimEvent> matches = new ArrayList<>();
		for (RSQSimEvent event : events) {
			Map<Double, EventComparison> row = gmpeVals.row(event.getID());
			double dist = row.values().iterator().next().distanceRup;
			if (dist >= minDist && dist <= maxDist)
				matches.add(event);
		}
		return matches;
	}
	
	public List<RSQSimEvent> getWithinDistJB(List<RSQSimEvent> events, Table<Integer, Double, EventComparison> gmpeVals,
			double minDist, double maxDist) {
		List<RSQSimEvent> matches = new ArrayList<>();
		for (RSQSimEvent event : events) {
			Map<Double, EventComparison> row = gmpeVals.row(event.getID());
			double dist = row.values().iterator().next().distanceJB;
			if (dist >= minDist && dist <= maxDist)
				matches.add(event);
		}
		return matches;
	}
	
	public void plotScatter(List<RSQSimEvent> events, BBP_Site site, double period, AttenRelRef gmpe,
			Map<Integer, EventComparison> gmpeVals, List<String> binDescriptions, File outputDir, String prefix)
					throws IOException {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		for (RSQSimEvent event : events) {
			DiscretizedFunc spectra = getEventSpectra(event, site);
			double rsVal = spectra.getY(period);
			EventComparison gmpeComp = gmpeVals.get(event.getID());
			double gmpeVal = Math.exp(gmpeComp.logMean);
			xy.set(gmpeVal, rsVal);
		}
		
		String title = gmpe.getShortName()+" Comparison Scatter";
		String xAxisLabel = gmpe.getShortName()+" "+(float)period+" s SA (g)";
		String yAxisLabel = "RSQSim "+(float)period+" s SA (g)";
		
		GroundMotionScatterPlot.PLOT_WIDTH = 600;
		GroundMotionScatterPlot.WRITE_PDF = false;
		GroundMotionScatterPlot.plot(xy, xAxisLabel, yAxisLabel, binDescriptions, title, outputDir, prefix);
	}
	
	public void plotStandardNormal(List<RSQSimEvent> events, BBP_Site site, double[] periods, AttenRelRef gmpe,
			Table<Integer, Double, EventComparison> gmpeVals, List<String> binDescriptions, File outputDir,
			String prefix) throws IOException {
		
		List<PlotSpec> specs = new ArrayList<>();
		double maxY = 0.7d;
		double numStdDev = 3.75;
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(new Range(-numStdDev, numStdDev));
		
		List<Double> means = new ArrayList<>();
		
		int numBins;
		if (events.size() < 100)
			numBins = 10;
		else if (events.size() < 500)
			numBins = 40;
		else
			numBins = 100;
		
		Color stdDevColor = new Color(0, 150, 0);
		
		for (double period : periods) {
			HistogramFunction hist = new HistogramFunction(-numStdDev, numStdDev, numBins);
			
			double mean = 0d;
			for (RSQSimEvent event : events) {
				DiscretizedFunc spectra = getEventSpectra(event, site);
				// in log space
				double rsVal = Math.log(spectra.getY(period));
				EventComparison gmpeComp = gmpeVals.get(event.getID(), period);
				double gmpeVal = gmpeComp.logMean;
				
				double val = (rsVal - gmpeVal)/gmpeComp.stdDev;
				hist.add(hist.getClosestXIndex(val), 1d);
				mean += val;
			}
			mean /= events.size();
			means.add(mean);
			
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
	}
	
	private static int max_table_fig_columns = 3;
	
	public void generateGMPE_Page(File outputDir, AttenRelRef gmpeRef, double[] periods,
			boolean distJB) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		String distShortName;
		String distDescription;
		if (distJB) {
			distShortName = "rJB";
			distDescription = "Joyner-Boore distance (**rJB**), the shortest horizontal distance from a site to the "
					+ "surface projection of the rupture surface";
		} else {
			distShortName = "rRup";
			distDescription = "rupture surface distance (**rRup**), the shortest 3-d distance from a site to the "
					+ "rupture surface";
		}
		
		// header
		lines.add("# "+catalog.getName()+" BBP/"+gmpeRef.getShortName()+" GMPE Comparisons");
		lines.add("");
		lines.add("**GMPE: "+gmpeRef.getName()+"**");
		lines.add("");
		lines.add("Ruptures are binned by their moment magnitude (**Mw**) and the "+distDescription);
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		double maxCatalogMag = 0d;
		for (RSQSimEvent e : events)
			maxCatalogMag = Math.max(maxCatalogMag, e.getMagnitude());
		
		List<Range> magRanges = new ArrayList<>();
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
		List<String> magLabels = new ArrayList<>();
		for (Range magRange : magRanges)
			magLabels.add(optionalDigitDF.format(magRange.getLowerBound())+" < Mw < "
					+optionalDigitDF.format(magRange.getUpperBound()));
		List<String> magFileLabels = new ArrayList<>();
		for (Range magRange : magRanges)
			magFileLabels.add("mag_"+optionalDigitDF.format(magRange.getLowerBound())+"_"
					+optionalDigitDF.format(magRange.getUpperBound()));
		
		List<Range> distRanges = new ArrayList<>();
		// 0-10, 10-20, 20-40, 40-80, 80-160
		distRanges.add(new Range(0d, 10d));
		distRanges.add(new Range(10d, 20d));
		distRanges.add(new Range(20d, 40d));
		distRanges.add(new Range(40d, 80d));
		distRanges.add(new Range(80d, 160d));
		distRanges.add(new Range(160d, MPJ_BBP_CatalogSim.CUTOFF_DIST));
		List<String> distLabels = new ArrayList<>();
		for (Range distRange : distRanges)
			distLabels.add(optionalDigitDF.format(distRange.getLowerBound())+" km < "+distShortName+" < "
					+optionalDigitDF.format(distRange.getUpperBound())+" km");
		List<String> distFileLabels = new ArrayList<>();
		for (Range distRange : distRanges)
			distFileLabels.add("dist_"+optionalDigitDF.format(distRange.getLowerBound())+"_"
					+optionalDigitDF.format(distRange.getUpperBound()));
		
		for (BBP_Site site : sites) {
			System.out.println("Site: "+site.getName());
			lines.add("## Site "+site.getName());
			lines.add(topLink); lines.add("");
			Location loc = site.getLoc();
			lines.add("*Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
			
			List<RSQSimEvent> eventsForSite = getEventsForSite(site);
			lines.add(eventsForSite.size()+" ruptures within "+(float)MPJ_BBP_CatalogSim.CUTOFF_DIST+" km");
			
			System.out.println("Calculating for "+site.getName()+", "+gmpeRef.getShortName());
			Table<Integer, Double, EventComparison> gmpeVals = calcGMPE(site, gmpeRef, periods);
			
			for (int m=0; m<magRanges.size(); m++) {
				Range magRange = magRanges.get(m);
				lines.add("### "+site.getName()+", "+magLabels.get(m));
				
				List<RSQSimEvent> eventsForMag = new ArrayList<>(eventsForSite);
				eventsForMag.removeIf((RSQSimEvent e) -> e.getMagnitude() < magRange.getLowerBound()
						|| e.getMagnitude() > magRange.getUpperBound());
				lines.add(eventsForMag.size()+" Ruptures");
				
				List<List<RSQSimEvent>> eventsForMagDists = new ArrayList<>();
				for (Range distRange : distRanges) {
					if (distJB)
						eventsForMagDists.add(getWithinDistJB(eventsForMag, gmpeVals,
								distRange.getLowerBound(), distRange.getUpperBound()));
					else
						eventsForMagDists.add(getWithinDistRup(eventsForMag, gmpeVals,
								distRange.getLowerBound(), distRange.getUpperBound()));
				}
				
				lines.add("#### "+site.getName()+", "+magLabels.get(m)+", Scatter Plots");
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
					List<RSQSimEvent> eventsForMagDist = eventsForMagDists.get(d);
					if (eventsForMagDist.isEmpty()) {
						System.out.println("No events for "+site.getName()+", "+distLabels.get(d)+", "+magLabels.get(m));
					}
					
					table.initNewLine().addColumn("**"+distLabels.get(d)+"**");
					
					for (double period : periods) {
						String prefix = site.getName()+"_"+magFileLabels.get(m)+"_"+distFileLabels.get(d)
						+"_"+optionalDigitDF.format(period)+"s_"+gmpeRef.getShortName()+"_scatter";
				
						System.out.println("Plotting Scatter: "+prefix);
						
						if (eventsForMagDist.isEmpty()) {
							table.addColumn("N/A");
						} else {
							List<String> binDescriptions = Lists.newArrayList(distLabels.get(d),
									magLabels.get(m), optionalDigitDF.format(period)+"s SA, "+gmpeRef.getShortName());
							plotScatter(eventsForMagDist, site, period, gmpeRef, gmpeVals.column(period),
									binDescriptions, resourcesDir, prefix);
							File scatterPlot = new File(resourcesDir, prefix+".png");
							Preconditions.checkState(scatterPlot.exists());
							table.addColumn("![Scatter Plot]("+resourcesDir.getName()
							+"/"+scatterPlot.getName()+")");
						}
					}
					
					table.finalizeLine();
				}
				lines.add("");
				lines.addAll(table.wrap(max_table_fig_columns, 1).build());
				
				lines.add("#### "+site.getName()+", "+magLabels.get(m)+", Standard Normal Plots");
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
				
				for (int d=0; d<distRanges.size(); d++) {
					List<RSQSimEvent> eventsForMagDist = eventsForMagDists.get(d);
					if (eventsForMagDist.isEmpty()) {
						System.out.println("No events for "+site.getName()+", "+distLabels.get(d)+", "+magLabels.get(m));
						continue;
					}
					table.addColumn("**"+distLabels.get(d)+"**");
				}
				table.finalizeLine();
				
				table.initNewLine();
				for (int d=0; d<distRanges.size(); d++) {
					List<RSQSimEvent> eventsForMagDist = eventsForMagDists.get(d);
					if (eventsForMagDist.isEmpty())
						continue;
					
					String prefix = site.getName()+"_"+magFileLabels.get(m)+"_"+distFileLabels.get(d)
						+"_"+gmpeRef.getShortName()+"_std_norm";
					
					System.out.println("Plotting Standard Normal: "+prefix);
					
					List<String> binDescriptions = Lists.newArrayList(distLabels.get(d), magLabels.get(m));
					plotStandardNormal(eventsForMagDist, site, periods, gmpeRef, gmpeVals,
							binDescriptions, resourcesDir, prefix);
					File plotFile = new File(resourcesDir, prefix+".png");
					Preconditions.checkState(plotFile.exists());
					table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
						+"/"+plotFile.getName()+")");
				}
				table.finalizeLine();
				lines.add("");
				lines.addAll(table.wrap(max_table_fig_columns, 0).build());
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	public static void main(String[] args) throws ZipException, IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_2273.instance(baseDir);
		
		VelocityModel vm = VelocityModel.LA_BASIN;
		double minFractForInclusion = 0.2;
		
//		double[] periods = { 1, 2, 3, 5, 10 };
		double[] periods = { 1, 2, 5 };
		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
		boolean distJB = true;
		
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-all") && new File(dir, "results.zip").exists())
				bbpDir = dir;
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
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		ZipFile zipFile = new ZipFile(new File(bbpDir, "results.zip"));
		
		RSQSimCatalogGMPE_Compare comp = new RSQSimCatalogGMPE_Compare(catalog, zipFile, sites, minMag, skipYears,
				vm, minFractForInclusion, gmpeCacheDir);
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		try {
			for (AttenRelRef gmpeRef : gmpeRefs) {
				File catalogGMPEDir = new File(catalogOutputDir, "gmpe_bbp_comparisons_"+gmpeRef.getShortName());
				Preconditions.checkState(catalogGMPEDir.exists() || catalogGMPEDir.mkdir());
				comp.generateGMPE_Page(catalogGMPEDir, gmpeRef, periods, distJB);
			}
			
			catalog.writeMarkdownSummary(catalogOutputDir);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		comp.shutdown();
	}

}
