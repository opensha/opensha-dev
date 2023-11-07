package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.PrimitiveArrayXY_Dataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.impl.WarningDoubleParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.MultiIMR_Averaged_AttenRel;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;

import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.simCompare.RuptureComparisonFilter.HasSiteFilter;
import scratch.kevin.simCompare.RuptureComparisonFilter.MatchesSiteFilter;
import scratch.kevin.simCompare.ZScoreHistPlot.ZScoreResult;

public abstract class MultiRupGMPE_ComparePageGen<E> {
	
	private SimulationRotDProvider<E> simProv;
	private String simName;
	private List<? extends Site> sites;
	private boolean distJB;
	private double cutoffDist;
	private double minMag;
	private double maxMag;
	
	private List<List<Site>> siteBundles = new ArrayList<>();
	private List<String> siteBundleNames = new ArrayList<>();
	
	private List<Range> magRanges;
	private List<RuptureComparisonFilter<E>> magFilters;
	private List<String> magLabels;
	private List<String> magFileLabels;
	
	private List<Range> distRanges;
	private List<RuptureComparisonFilter<E>> distFilters;
	private List<String> distLabels;
	private List<String> distFileLabels;
	
	private boolean doRakeBinned = true;
	private List<RuptureComparisonFilter<E>> rakeFilters;
	private List<String> rakeLabels;
	private List<String> rakeFileLabels;
	
	private Table<E, Site, Map<Integer, Double>> rupSiteAzMap;
	
	protected ExecutorService exec;
	
	private Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesInstancesCache;
	
	private boolean replotScatters = true;
	private boolean replotZScores = true;
	private boolean replotCurves = true;
	private boolean replotResiduals = true;
	
	private int sourceRupContributionNum;
	private Table<String, E, Double> sourceRupContributionFracts;
	
	protected void init(SimulationRotDProvider<E> simProv, List<? extends Site> sites, boolean distJB, double cutoffDist,
			double minMag, double maxMag) {
		this.simProv = simProv;
		this.simName = simProv.getName();
		this.sites = sites;
		this.distJB = distJB;
		this.cutoffDist = cutoffDist;
		this.minMag = minMag;
		this.maxMag = maxMag;
		
		initBins();
		
		gmpesInstancesCache = new HashMap<>();
		rupSiteAzMap = HashBasedTable.create();
		
		int threads = Runtime.getRuntime().availableProcessors();
		if (threads > 8)
			threads = Integer.min(32, threads-2);
		
		exec = Executors.newFixedThreadPool(threads);
	}
	
	protected void addSiteBundle(List<Site> sites, String name) {
		siteBundles.add(sites);
		siteBundleNames.add(name);
	}
	
	protected List<? extends Site> getSites() {
		return sites;
	}
	
	public void setReplotScatters(boolean replotScatters) {
		this.replotScatters = replotScatters;
	}

	public void setReplotZScores(boolean replotZScores) {
		this.replotZScores = replotZScores;
	}

	public void setReplotCurves(boolean replotCurves) {
		this.replotCurves = replotCurves;
	}

	public void setReplotResiduals(boolean replotResiduals) {
		this.replotResiduals = replotResiduals;
	}
	
	public void setSourceRupContributionFractions(Table<String, E, Double> sourceRupContribFracts,
			int numSourcesToPlot) {
		this.sourceRupContributionFracts = sourceRupContribFracts;
		this.sourceRupContributionNum = numSourcesToPlot;
	}

	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	protected String getDistDescription() {
		if (distJB)
			return "Joyner-Boore distance (**rJB**), the shortest horizontal distance from a site to the "
					+ "surface projection of the rupture surface";
		return "rupture surface distance (**rRup**), the shortest 3-d distance from a site to the "
					+ "rupture surface";
	}
	
	protected String getDistShortName() {
		if (distJB)
			return "rJB";
		return "rRup";
	}
	
	private void initBins() {
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
		if (maxMag > 8d) {
			if (maxMag < 8.5)
				magRanges.add(new Range(8d, 8.5));
			else
				magRanges.add(new Range(8d, 9));
		}
		magFilters = new ArrayList<>();
		for (Range magRange : magRanges)
			magFilters.add(new RuptureComparisonFilter.MagFilter<>(magRange.getLowerBound(), magRange.getUpperBound()));
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
		distRanges.add(new Range(10d, 40d));
		distRanges.add(new Range(40d, 80d));
		distRanges.add(new Range(80d, 160d));
		distRanges.add(new Range(160d, cutoffDist));
		distRanges.add(null);
		distFilters = new ArrayList<>();
		for (Range distRange : distRanges) {
			if (distRange == null)
				distFilters.add(new RuptureComparisonFilter.AcceptAllFilter<>());
			else if (distJB)
				distFilters.add(new RuptureComparisonFilter.DistJBFilter<>(distRange.getLowerBound(),
						distRange.getUpperBound()));
			else
				distFilters.add(new RuptureComparisonFilter.DistRupFilter<>(distRange.getLowerBound(),
						distRange.getUpperBound()));
		}
		distLabels = new ArrayList<>();
		distFileLabels = new ArrayList<>();
		String distShortName = getDistShortName();
		for (Range distRange : distRanges) {
			if (distRange == null) {
				distLabels.add("All Distances");
				distFileLabels.add("dist_all");
			} else {
				distLabels.add(optionalDigitDF.format(distRange.getLowerBound())+" km < "+distShortName+" < "
						+optionalDigitDF.format(distRange.getUpperBound())+" km");
				distFileLabels.add("dist_"+optionalDigitDF.format(distRange.getLowerBound())+"_"
						+optionalDigitDF.format(distRange.getUpperBound()));
			}
		}
		
		List<RakeFilter> indvRakeFilters = new ArrayList<>(3);
		rakeLabels = new ArrayList<>(4);
		rakeFileLabels = new ArrayList<>(4);
		
		indvRakeFilters.add(new RakeFilter(new Range(-180.01, -170), new Range(-10, 10), new Range(170, 180.01)));
		rakeLabels.add("Strike-Slip");
		rakeFileLabels.add("rake_ss");
		
		indvRakeFilters.add(new RakeFilter(new Range(80, 100)));
		rakeLabels.add("Reverse");
		rakeFileLabels.add("rake_reverse");
		
		indvRakeFilters.add(new RakeFilter(new Range(-100, -80)));
		rakeLabels.add("Normal");
		rakeFileLabels.add("rake_normal");
		
		rakeFilters = new ArrayList<>(4);
		rakeFilters.addAll(indvRakeFilters);
		rakeFilters.add(new ObliqueFilter(indvRakeFilters));
		rakeLabels.add("Oblique");
		rakeFileLabels.add("rake_oblique");
	}
	
	public boolean isDoRakeBinned() {
		return doRakeBinned;
	}

	public void setDoRakeBinned(boolean doRakeBinned) {
		this.doRakeBinned = doRakeBinned;
	}

	private class RakeFilter extends RuptureComparisonFilter<E> {
		
		private Range[] ranges;

		private RakeFilter(Range... ranges) {
			this.ranges = ranges;
			
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			double rake = simProv.getRake(comp.getRupture());
			
			for (Range range : ranges)
				if (range.contains(rake))
					return true;
			return false;
		}
		
	}
	
	private class ObliqueFilter extends RuptureComparisonFilter<E> {
		
		private List<Range> ranges;

		private ObliqueFilter(List<RakeFilter> rakeFilters) {
			ranges = new ArrayList<>();
			for (RakeFilter filter : rakeFilters)
				for (Range range : filter.ranges)
					ranges.add(range);
		}

		@Override
		public boolean matches(RuptureComparison<E> comp, Site site) {
			double rake = simProv.getRake(comp.getRupture());
			
			if (!Double.isFinite(rake))
				return false;
			
			// see if any of our other ranges contain it
			for (Range range : ranges)
				if (range.contains(rake))
					return false;
			return true;
		}
		
	}
	
	private static final int MAX_SCATTER_NUM_POINTS = 100000;
	
	public boolean plotScatter(Collection<? extends RuptureComparison<E>> eventComps, Collection<Site> sites, IMT imt,
			AttenRelRef gmpe, RuptureComparisonFilter<E> filter, List<String> binDescriptions, File outputDir, String prefix)
					throws IOException {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		int numComputed = 0;
		int numFiltered = 0;
		
		for (Site site : sites) {
			for (RuptureComparison<E> comp : eventComps) {
				if (!comp.isComputed(site, imt))
					continue;
				numComputed++;
//					System.out.println("Not computed? "+site.getName()+" "+((RSQSimEvent)comp.getRupture()).getID());
//				System.out.println("Scatter for "+comp.getMagnitude()+" "+site.getName()+" "+comp.getDistanceJB(site));
				if ((filter != null && !filter.matches(comp, site))) {
					numFiltered++;
					continue;
				}
				double gmpeVal = Math.exp(comp.getLogMean(site, imt));
				for (double simVal : simProv.getValues(site, comp.getRupture(), imt))
					xy.set(gmpeVal, simVal);
			}
		}
		System.out.println("XY size: "+xy.size());
		System.out.println(numFiltered+"/"+numComputed+" filtered out");
		if (xy.size() == 0)
			return false;
		if (xy.size() > MAX_SCATTER_NUM_POINTS) {
			Random r = new Random(xy.size());
			DefaultXY_DataSet oXY = new DefaultXY_DataSet();
			double p = (double)MAX_SCATTER_NUM_POINTS/xy.size();
			for (Point2D pt : xy)
				if (r.nextDouble() < p)
					oXY.set(pt);
			System.out.println("Filtered further down to "+oXY.size()+" points for plotting");
			xy = oXY;
		}
		
		String title = gmpe.getShortName()+" Comparison Scatter";
		String periodLabel = imt.getShortName()+" ("+imt.getUnits()+")";
		String xAxisLabel = gmpe.getShortName()+" "+periodLabel;
		String yAxisLabel = simName+" "+periodLabel;
		
		GroundMotionScatterPlot.PLOT_WIDTH = 600;
		GroundMotionScatterPlot.WRITE_PDF = false;
		GroundMotionScatterPlot.plot(xy, xAxisLabel, yAxisLabel, binDescriptions, title, outputDir, prefix);
		return true;
	}
	
	private static int max_table_fig_columns = 3;
	
	static int[] hazard_curve_rps = { 1000, 2500, 10000 };
	
	private static PlotLineType[] gmpe_alt_line_types = { PlotLineType.DASHED,
			PlotLineType.DOTTED, PlotLineType.DOTTED_AND_DASHED };
	
	private static double[] gmpe_truncs = { 3d, 2d, 1d };
//	private static double[] gmpe_fixed_sigmas = { 0d, 0.3 };
	private static double[] gmpe_fixed_sigmas = { 0d };
	
	private static PlotCurveCharacterstics simCurveChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
	private static PlotCurveCharacterstics simUncertainChar =
			new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120, 120, 120, 120));
	private static PlotCurveCharacterstics[] simCompCurveChars = {
			new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 2f, Color.DARK_GRAY)
	};
	
	public List<File> plotHazardCurves(List<SimulationHazardPlotter<E>> curvePlotters, IMT imt, File outputDir) throws IOException {
		List<Future<File>> futures = new ArrayList<>();
		
		for (SimulationHazardPlotter<E> curvePlotter : curvePlotters)
			futures.add(exec.submit(new CurveCalcCallable(curvePlotter, imt, outputDir)));
		
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
		private SimulationHazardPlotter<E> curvePlotter;
		private IMT imt;
		private File outputDir;

		public CurveCalcCallable(SimulationHazardPlotter<E> curvePlotter, IMT imt, File outputDir) {
			this.curvePlotter = curvePlotter;
			this.imt = imt;
			this.outputDir = outputDir;
		}

		@Override
		public File call() throws Exception {
			String siteName = curvePlotter.getSite().getName();
			String prefix = siteName.replaceAll(" ", "_")+"_curves_"+imt.getPrefix();
			prefix += "_"+curvePlotter.getGMPE().getShortName();
			File pngFile = new File(outputDir, prefix+".png");
			
			if (pngFile.exists() && !replotCurves)
				return pngFile;
			
			System.out.println("Calculating simulation curve for "+curvePlotter.getSite().getName()+", "+imt);
			return curvePlotter.plotHazardCurves(outputDir, prefix, imt);
		}
		
	}
	
	protected ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
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
		return gmpe;
	}
	
	protected void checkInGMPE(AttenRelRef gmpeRef, ScalarIMR gmpe) {
		synchronized (gmpesInstancesCache) {
			gmpesInstancesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	@SuppressWarnings("unchecked")
	public File plotHazardCurve(SimulationHazardPlotter<E> simCurvePlot, Site site, double period, AttenRelRef gmpeRef, File outputDir)
			throws IOException {
		String siteName = site.getName();
		String prefix = siteName.replaceAll(" ", "_")+"_curves_"+(float)period+"s_"+gmpeRef.getShortName();
		File pngFile = new File(outputDir, prefix+".png");
		
		if (pngFile.exists() && !replotCurves)
			return pngFile;
		
		System.out.println("Calculating simulation curve for "+site.getName()+", t="+(float)period);
		return simCurvePlot.plotHazardCurves(outputDir, prefix, period);
	}
	
	private static DecimalFormat twoSigFig = new DecimalFormat("0.00");
	
	private static List<String> getZscoreTableLine(String siteName, ZScoreResult[] scores, IMT[] imts) {
		List<String> line = new ArrayList<>();
		
		line.add(siteName);
		if (scores == null) {
			for (int i=0; i<imts.length; i++) {
				line.add("__N/A__");
				line.add("__N/A__");
			}
		} else {
			for (ZScoreResult score : scores) {
				line.add(twoSigFig.format(score.mean));
				line.add(twoSigFig.format(score.stdDevFract));
			}
		}
		
		return line;
	}
	
	protected void checkTransSiteParams(AttenRelRef gmpeRef) {
		SiteTranslator trans = null;
		ScalarIMR gmmInstance = gmpeRef.get();
		Site site0 = sites.get(0);
		if (site0.containsParameter(Vs30_Param.NAME)) {
			for (Parameter<?> gmmSiteParam : gmmInstance.getSiteParams()) {
				if (!site0.containsParameter(gmmSiteParam)) {
					System.out.println("Trying to translate site params for "+gmmSiteParam.getName()+", required by "+gmmInstance.getName());
					if (trans == null)
						trans = new SiteTranslator();
					
					for (Site site : sites) {
						String measurement = null;
						if (site.containsParameter(Vs30_TypeParam.NAME)) {
							if (site.getParameter(Vs30_TypeParam.NAME).getValue().equals(Vs30_TypeParam.VS30_TYPE_MEASURED))
								measurement = SiteData.TYPE_FLAG_MEASURED;
							else
								measurement = SiteData.TYPE_FLAG_INFERRED;
						}
						SiteDataValue<Double> val = new SiteDataValue<Double>(SiteData.TYPE_VS30, measurement,
								site.getParameter(Double.class, Vs30_Param.NAME).getValue());
						trans.setAllSiteParams(gmmInstance, val);
						for (Parameter<?> param : gmmInstance.getSiteParams()) {
							if (!site.containsParameter(param)) {
								Parameter<?> paramClone = (Parameter<?>)param.clone();
								System.out.println("\tTranslated Vs30="+val.getValue().floatValue()+" to '"
										+param.getValue()+"' for site "+site.getName());
								site.addParameter(paramClone);
							}
						}
					}
				}
			}
		}
	}
	
	public void generateGMPE_Page(File outputDir, List<String> headerLines, AttenRelRef gmpeRef, IMT[] imts,
			List<? extends RuptureComparison<E>> comps, List<Site> highlightSites) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.getLast().isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		// see if we need to translate any site types
		checkTransSiteParams(gmpeRef);
		
		List<Site> sites = new ArrayList<>();
		if (highlightSites == null)
			// make plots for all sites
			sites.addAll(this.sites);
		else
			// make plots for only the specified, use all for aggregations
			sites.addAll(highlightSites);
		// add null for all sites aggregated
		sites.add(0, null);
		
		lines.add("## Site Scatters/z-score Histograms");
		lines.add(topLink); lines.add("");
		
		lines.add("### z-score Summary Table");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Site");
		for (IMT imt : imts) {
			table.addColumn(imt.getShortName()+" Mean");
			table.addColumn(imt.getShortName()+" &sigma;-fract");
		}
		table.finalizeLine();
		
		ZScoreResult[] scores = ZScoreHistPlot.calcZScores(simProv, comps, this.sites, imts, null);
		table.addLine(getZscoreTableLine("All Sites", scores, imts));
		
		List<List<RuptureComparison<E>>> bundleSiteComps = null;
		if (this.sites.size() > 1) {
			HasSiteFilter<E> siteFilter = new RuptureComparisonFilter.HasSiteFilter<E>();
			if (!siteBundles.isEmpty()) {
				bundleSiteComps = new ArrayList<>();
				for (int s=0; s<siteBundles.size(); s++) {
					String name = siteBundleNames.get(s);
					
					List<Site> scatterSites = siteBundles.get(s);
					List<RuptureComparison<E>> siteComps = new ArrayList<>();
					
					for (RuptureComparison<E> comp : comps) {
						for (Site site : scatterSites) {
							if (siteFilter.matches(comp, site)) {
								siteComps.add(comp);
								break;
							}
						}
					}
					bundleSiteComps.add(siteComps);
					
					scores = ZScoreHistPlot.calcZScores(simProv, siteComps, scatterSites, imts, null);
					table.addLine(getZscoreTableLine(name, scores, imts));
				}
			}
			for (Site site : this.sites) {
				if (site == null)
					continue;
				List<Site> mySites = new ArrayList<>();
				mySites.add(site);
				
				List<? extends RuptureComparison<E>> siteComps = siteFilter.getMatches(comps, site);
				scores = ZScoreHistPlot.calcZScores(simProv, siteComps, mySites, imts, null);
				table.addLine(getZscoreTableLine(site.getName(), scores, imts));
			}
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		if (!siteBundles.isEmpty() && sites.size() > 1) {
			System.out.println("All sites aggregated");
			lines.add("### Bundled z-scores");
			lines.add(topLink); lines.add("");
			
			for (int s=0; s<siteBundles.size(); s++) {
				String name = siteBundleNames.get(s);
				String prefix = "site_bundle_"+name.replaceAll("\\W+", "_")+"_std_norm";
				
				System.out.println("Processing site bundle: "+name);
				
				List<Site> scatterSites = siteBundles.get(s);
				List<RuptureComparison<E>> siteComps = bundleSiteComps.get(s);
				
				File plotFile = new File(resourcesDir, prefix+".png");
				boolean success;
				if (plotFile.exists() && !replotZScores)
					success = true;
				else
					success = ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, imts,
						gmpeRef, null, new ArrayList<>(), resourcesDir, prefix);
				if (success) {
					lines.add("#### "+name+" z-scores");
					lines.add(topLink); lines.add("");
					lines.add("");
					lines.add("Site bundle with "+scatterSites.size()+" sites.");
					lines.add("");
					lines.add("z-score standard normal plots across all magnitudes/distances");
					lines.add("");
					lines.add("**z-score**: (ln(*"+simName+"*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
					lines.add("");
					lines.add("**Legend**");
					lines.add("* Black Line: Standard Normal distribution (in natural log space)");
					lines.add("* Gray Histogram: z-score for each rupture");
					lines.add("* Blue Dashed Line: "+simName+" Mean");
					
					lines.add("");
					Preconditions.checkState(plotFile.exists());
					
					String srcPrefix = prefix+"_source_contrib";
					File srcPlotFile = new File(resourcesDir, srcPrefix+".png");
					boolean srcSuccess = sourceRupContributionFracts != null &&
							((srcPlotFile.exists() && !replotZScores) ||
							ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, imts,
							gmpeRef, null, new ArrayList<>(), resourcesDir, srcPrefix, sourceRupContributionFracts,
							sourceRupContributionNum, false));
					
					if (srcSuccess) {
						// plot with source contributions
						table = MarkdownUtils.tableBuilder();
						table.initNewLine();
						table.addColumn("Total");
						table.addColumn("Source Contributions");
						table.finalizeLine();
						table.initNewLine();
						table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
								+"/"+plotFile.getName()+")");
						
						table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
								+"/"+srcPlotFile.getName()+")");
						table.finalizeLine();
						lines.addAll(table.build());

						for (IMT imt : imts) {
							String periodPrefix = prefix+"_"+imt.getPrefix();
							File periodPlotFile = new File(resourcesDir, periodPrefix+".png");
							if (!periodPlotFile.exists() || replotZScores)
								ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
										gmpeRef, null, new ArrayList<>(), resourcesDir, periodPrefix, null, 0, true);
							String srcPeriodPrefix = periodPrefix+"_source_contrib";
							File srcPeriodPlotFile = new File(resourcesDir, srcPeriodPrefix+".png");
							if (!srcPeriodPlotFile.exists() || replotZScores)
								ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
										gmpeRef, null, new ArrayList<>(), resourcesDir, srcPeriodPrefix,
										sourceRupContributionFracts, sourceRupContributionNum, true);
						}
					} else {
						
						lines.add("![Standard Normal Plot]("+resourcesDir.getName()
							+"/"+plotFile.getName()+")");

						for (IMT imt : imts) {
							String periodPrefix = prefix+"_"+imt.getPrefix();
							File periodPlotFile = new File(resourcesDir, periodPrefix+".png");
							if (!periodPlotFile.exists() || replotZScores)
								ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
										gmpeRef, null, new ArrayList<>(), resourcesDir, periodPrefix, null, 0, true);
						}
					}
					
				}
			}
		}
		
		for (Site site : sites) {
			List<? extends RuptureComparison<E>> siteComps;
			List<Site> scatterSites;
			String siteName;
			if (site == null) {
				System.out.println("All sites aggregated");
				lines.add("### All Sites Aggregated");
				lines.add(topLink); lines.add("");
				lines.add("**"+this.sites.size()+" sites**");
				lines.add("");
				table = MarkdownUtils.tableBuilder();
				table.addLine("Name", "Location", "# Ruptures", "Vs30 (m/s)", "Z1.0 (km)", "Z2.5 (km)");
				for (Site s : this.sites) {
					table.initNewLine();
					table.addColumn(s.getName());
					Location loc = s.getLocation();
					table.addColumn("*"+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
					table.addColumn(simProv.getRupturesForSite(s).size()+" ("+simProv.getNumSimlationsForSite(s)+" sims)");
					table.addColumn(optionalDigitDF.format(s.getParameter(Double.class, Vs30_Param.NAME).getValue()));
					Double z1 = s.getParameter(Double.class, DepthTo1pt0kmPerSecParam.NAME).getValue();
					if (z1 == null || Double.isNaN(z1))
						table.addColumn("N/A");
					else
						table.addColumn(optionalDigitDF.format(z1/1000d));
					Double z25 = s.getParameter(Double.class, DepthTo2pt5kmPerSecParam.NAME).getValue();
					if (z25 == null || Double.isNaN(z25))
						table.addColumn("N/A");
					else
						table.addColumn(optionalDigitDF.format(z25));
					table.finalizeLine();
				}
				lines.addAll(table.build());
				siteComps = comps;
				siteName = "All Sites";
				scatterSites = new ArrayList<>(this.sites);
				lines.add("");
				lines.add(siteComps.size()+" ruptures within "+optionalDigitDF.format(cutoffDist)+" km of *any* site");
			} else {
				System.out.println("Site: "+site.getName());
				siteName = site.getName();
				lines.add("### Site "+site.getName());
				lines.add(topLink); lines.add("");
				Location loc = site.getLocation();
				lines.add("*Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
				siteComps = new RuptureComparisonFilter.HasSiteFilter<E>().getMatches(comps, site);
				scatterSites = new ArrayList<>();
				scatterSites.add(site);
				
				lines.add(siteComps.size()+" ruptures within "+(float)cutoffDist+" km");
			}
			
			boolean[] magRakeBools = doRakeBinned ? new boolean[] { false, true } : new boolean[] { false };
			
			for (boolean isRake : magRakeBools) {
				List<RuptureComparisonFilter<E>> filters;
				List<String> labels;
				List<String> fileLabels;
				
				if (isRake) {
					filters = rakeFilters;
					labels = rakeLabels;
					fileLabels = rakeFileLabels;
				} else {
					filters = magFilters;
					labels = magLabels;
					fileLabels = magFileLabels;
				}
				
				for (int f=0; f<filters.size(); f++) {
					RuptureComparisonFilter<E> filter = filters.get(f);

					List<? extends RuptureComparison<E>> filteredEventComps = filter.getMatches(comps, null);
					if (filteredEventComps.isEmpty())
						continue;
					
					lines.add("#### "+siteName+", "+labels.get(f));
					lines.add("");
					lines.add(filteredEventComps.size()+" Ruptures");
					
					lines.add("##### "+siteName+", "+labels.get(f)+", Scatter Plots");
					lines.add(topLink); lines.add("");
					lines.add("**Legend**");
					lines.add("* Red +: GMPE Mean/"+simName+" single rupture comparison");
					lines.add("* Yellow Region: Factor of 2 above & below");
					lines.add("* Green Line: Linear Regression");
					
					table = MarkdownUtils.tableBuilder();
					table.initNewLine().addColumn("**Distance Bin**");
					for (IMT imt : imts)
						table.addColumn("**"+imt.getDisplayName()+"**");
					table.finalizeLine();
					
					for (int d=0; d<distFilters.size(); d++) {
						RuptureComparisonFilter<E> distFilter = distFilters.get(d);
//						if (eventsForMagDist.isEmpty()) {
//							System.out.println("No events for "+site.getName()+", "+distLabels.get(d)+", "+magLabels.get(m));
//						}
						
						table.initNewLine().addColumn("**"+distLabels.get(d)+"**");
						
						for (IMT imt : imts) {
							String prefix = siteName.replaceAll(" ", "_")+"_"+fileLabels.get(f)+"_"+distFileLabels.get(d);
							prefix += "_"+imt.getPrefix()+"_"+gmpeRef.getShortName()+"_scatter";
					
							System.out.println("Plotting Scatter: "+prefix);
							
							String label = imt.getShortName();
							
							List<String> binDescriptions = Lists.newArrayList(distLabels.get(d),
									labels.get(f), label+", "+gmpeRef.getShortName());
							File scatterPlot = new File(resourcesDir, prefix+".png");
							boolean success;
							if (scatterPlot.exists() && !replotScatters)
								success = true;
							else
								success = plotScatter(filteredEventComps, scatterSites, imt, gmpeRef, distFilter,
									binDescriptions, resourcesDir, prefix);
							if (success) {
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
					
					lines.add("##### "+siteName+", "+labels.get(f)+", z-Score Histograms");
					lines.add(topLink); lines.add("");
					lines.add("These plots compare "+simName+" to the full GMPE log-normal distributions. "
							+ "Each rupture's GMPE distribution is converted to a standard log-normal "
							+ "distribution, and the z-score is computed for each rupture:");
					lines.add("");
					lines.add("**z-score**: (ln(*"+simName+"*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
					lines.add("");
					lines.add("**Legend**");
					lines.add("* Black Line: Standard Normal distribution (in natural log space)");
					lines.add("* Gray Histogram: z-score for each rupture");
					lines.add("* Blue Dashed Line: "+simName+" Mean");
					
					table = MarkdownUtils.tableBuilder();
					table.initNewLine();
					
					for (int d=0; d<distFilters.size(); d++)
						table.addColumn("**"+distLabels.get(d)+"**");
					table.finalizeLine();
					
					table.initNewLine();
					for (int d=0; d<distFilters.size(); d++) {
						RuptureComparisonFilter<E> distFilter = distFilters.get(d);
						String prefix = siteName.replaceAll(" ", "_")+"_"+fileLabels.get(f)+"_"+distFileLabels.get(d)
							+"_"+gmpeRef.getShortName()+"_std_norm";
						
						System.out.println("Plotting Standard Normal: "+prefix);
						
						List<String> binDescriptions = Lists.newArrayList(distLabels.get(d), labels.get(f));
						File plotFile = new File(resourcesDir, prefix+".png");
						boolean success;
						if (plotFile.exists() && !replotZScores)
							success = true;
						else
							success = ZScoreHistPlot.plotStandardNormal(simProv, filteredEventComps, scatterSites, imts,
								gmpeRef, distFilter, binDescriptions, resourcesDir, prefix);
						if (success) {
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
			}
			
			String prefix = siteName.replaceAll(" ", "_")+"_all_mags_all_dists_"+gmpeRef.getShortName()+"_std_norm";
			File plotFile = new File(resourcesDir, prefix+".png");
			boolean success;
			if (plotFile.exists() && !replotZScores)
				success = true;
			else
				success = ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, imts,
					gmpeRef, null, new ArrayList<>(), resourcesDir, prefix);
			if (success) {
				lines.add("#### "+siteName+", All Ruptures, Z-Score Histograms");
				lines.add(topLink); lines.add("");
				lines.add("");
				lines.add("z-score standard normal plots across all magnitudes/distances");
				lines.add("");
				lines.add("**z-score**: (ln(*"+simName+"*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
				lines.add("");
				lines.add("**Legend**");
				lines.add("* Black Line: Standard Normal distribution (in natural log space)");
				lines.add("* Gray Histogram: z-score for each rupture");
				lines.add("* Blue Dashed Line: "+simName+" Mean");
				
				lines.add("");
				Preconditions.checkState(plotFile.exists());
				
				String srcPrefix = prefix+"_source_contrib";
				File srcPlotFile = new File(resourcesDir, srcPrefix+".png");
				boolean srcSuccess = sourceRupContributionFracts != null &&
						((srcPlotFile.exists() && !replotZScores) ||
						ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, imts,
						gmpeRef, null, new ArrayList<>(), resourcesDir, srcPrefix, sourceRupContributionFracts,
						sourceRupContributionNum, false));
				
				if (srcSuccess) {
					// plot with source contributions
					table = MarkdownUtils.tableBuilder();
					table.initNewLine();
					table.addColumn("Total");
					table.addColumn("Source Contributions");
					table.finalizeLine();
					table.initNewLine();
					table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
							+"/"+plotFile.getName()+")");
					
					table.addColumn("![Standard Normal Plot]("+resourcesDir.getName()
							+"/"+srcPlotFile.getName()+")");
					table.finalizeLine();
					lines.addAll(table.build());

					for (IMT imt : imts) {
						String periodPrefix = prefix+"_"+imt.getPrefix();
						File periodPlotFile = new File(resourcesDir, periodPrefix+".png");
						if (!periodPlotFile.exists() || replotZScores)
							ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
									gmpeRef, null, new ArrayList<>(), resourcesDir, periodPrefix, null, 0, true);
						String srcPeriodPrefix = periodPrefix+"_source_contrib";
						File srcPeriodPlotFile = new File(resourcesDir, srcPeriodPrefix+".png");
						if (!srcPeriodPlotFile.exists() || replotZScores)
							ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
									gmpeRef, null, new ArrayList<>(), resourcesDir, srcPeriodPrefix,
									sourceRupContributionFracts, sourceRupContributionNum, true);
					}
				} else {
					
					lines.add("![Standard Normal Plot]("+resourcesDir.getName()
					+"/"+plotFile.getName()+")");

					for (IMT imt : imts) {
						String periodPrefix = prefix+"_"+imt.getPrefix();
						File periodPlotFile = new File(resourcesDir, periodPrefix+".png");
						if (!periodPlotFile.exists() || replotZScores)
							ZScoreHistPlot.plotStandardNormal(simProv, siteComps, scatterSites, new IMT[] {imt},
									gmpeRef, null, new ArrayList<>(), resourcesDir, periodPrefix, null, 0, true);
					}
				}
				
			}
		}
		
		// now hazard curves
		List<List<File>> curveFiles = new ArrayList<>();
		List<Site> curveSites;
		if (this.sites.size() > 20 && highlightSites != null)
			curveSites = highlightSites;
		else
			curveSites = new ArrayList<>(this.sites);
		System.out.println("Have "+curveSites.size()+" curve sites (from "+sites.size()+" total sites)");
		
		if (curveSites != null && !curveSites.isEmpty()) {
			for (int s=curveSites.size(); --s>=0;) {
				if (simProv.getRupturesForSite(curveSites.get(s)).isEmpty()) {
					System.out.println("Skipping curves for site "+curveSites.get(s).getName()+", no ruptures");
					curveSites.remove(s);
				}
			}
			SimulationHazardCurveCalc<E> simCurveCalc = new SimulationHazardCurveCalc<>(simProv);
			List<SimulationHazardPlotter<E>> curvePlotters = new ArrayList<>();
			for (Site site : curveSites) {
				SimulationHazardPlotter<E> curvePlotter = new SimulationHazardPlotter<>(simCurveCalc, comps, site, 1d, gmpeRef);
				curvePlotter.setGMPE_FixedSigmas(gmpe_fixed_sigmas);
				curvePlotter.setGMPE_TruncationLevels(gmpe_truncs);
				curvePlotters.add(curvePlotter);
			}
			for (IMT imt : imts)
				curveFiles.add(plotHazardCurves(curvePlotters, imt, resourcesDir));
			
			lines.add("## Hazard Curves");
			lines.add(topLink); lines.add("");
//			lines.addAll(getCurveLegend(simName, gmpeRef.getShortName(), gmpe_truncs, gmpe_fixed_sigmas, 0));
			lines.addAll(curvePlotters.get(0).getCurveLegend(false, true, true, false, 0));
			lines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("Site");
			for (IMT imt : imts)
				table.addColumn(imt.getDisplayName());
			table.finalizeLine();
			for (int s=0; s<curveSites.size(); s++) {
				Site site = curveSites.get(s);
				table.initNewLine().addColumn("**"+site.getName()+"**");
				for (int p=0; p<imts.length; p++) {
					File plotFile = curveFiles.get(p).get(s);
					table.addColumn("![Hazard Curve]("+resourcesDir.getName()
								+"/"+plotFile.getName()+")");
				}
				table.finalizeLine();
			}
			lines.addAll(table.build());
			
			System.out.println("Done with curves");
		}
		
		// now GMPE residuals
		lines.add("## GMPE Residuals");
		lines.add(topLink); lines.add("");
		String residualLabel = "Residual: ln("+simProv.getName()+") - ln("+gmpeRef.getShortName()+")";
		lines.add("Residuals of simulation data ("+simProv.getName()+") in log space relative to GMPE log-mean");
		lines.add("");
		lines.add("**Legend**");
		lines.add("* Black Thick Line: Linear Least-Squares Fit to Residuals");
		lines.add("* Black Circles: Binned Linear Least-Squares Fit to Residuals");
		lines.add("  * Black Thin Dashes: binned mean ± data sigma");
		lines.add("  * Blue Thin Dotted: binned mean ± GMPE sigma");
		lines.add("");
		List<ResidualType> residualTypes = new ArrayList<>();
		residualTypes.add(ResidualType.MAG);
		ScalarIMR tempGMPE = checkOutGMPE(gmpeRef);
		if (tempGMPE.getPropagationEffectParams().containsParameter(DistanceJBParameter.NAME))
			residualTypes.add(ResidualType.DIST_JB);
		if (tempGMPE.getPropagationEffectParams().containsParameter(DistanceRupParameter.NAME))
			residualTypes.add(ResidualType.DIST_RUP);
		checkInGMPE(gmpeRef, tempGMPE);
		sites = new ArrayList<>(this.sites); // we previously added a null element for all sites above
		if (tempGMPE.getSiteParams().containsParameter(Vs30_Param.NAME) && doesParameterVary(Vs30_Param.NAME, sites))
			residualTypes.add(ResidualType.VS30);
		if (tempGMPE.getSiteParams().containsParameter(DepthTo1pt0kmPerSecParam.NAME) && doesParameterVary(DepthTo1pt0kmPerSecParam.NAME, sites))
			residualTypes.add(ResidualType.Z10);
		if (tempGMPE.getSiteParams().containsParameter(DepthTo2pt5kmPerSecParam.NAME) && doesParameterVary(DepthTo2pt5kmPerSecParam.NAME, sites))
			residualTypes.add(ResidualType.Z25);
		if (RuptureComparison.getRuptureTimeRange(comps) != null)
			residualTypes.add(ResidualType.OCCUR_TIME);
		checkInGMPE(gmpeRef, tempGMPE);
		
		lines.add("GMPE Residuals use the following values, averaged among all ruptures, for all paremeters which are not varied. "
				+ "All other parameters set to GMPE defaults");
		lines.add("");
		table = MarkdownUtils.tableBuilder().addLine("Name", "Average Value");
		Map<ResidualType, Double> residualDefaults = new HashMap<>();
		Map<ResidualType, Integer> regressValCounts = new HashMap<>();
		System.out.println("Determining GMPE default parameter values for residual analysis (averaged across all ruptures)");
		for (ResidualType type : ResidualType.values()) {
			Double meanVal = 0d;
			double sumWeights = 0d;
			int numVals = 0;
			for (RuptureComparison<E> comp : comps) {
				for (Site site : comp.getApplicableSites()) {
					Double x = type.getValue(comp, site);
					if (x == null) {
						Preconditions.checkState(meanVal == 0d || meanVal == null);
					} else {
						if (!Double.isNaN(type.minX) && x < type.minX)
							continue;
						meanVal += x*comp.getAnnualRate();
						sumWeights += comp.getAnnualRate();
						numVals += simProv.getNumSimulations(site, comp.getRupture());
					}
				}
			}
			regressValCounts.put(type, numVals);
			if (meanVal == null) {
				table.addLine(type.name, "(null)");
			} else {
				meanVal /= sumWeights;
				table.addLine(type.name, optionalDigitDF.format(meanVal));
			}
			residualDefaults.put(type, meanVal);
		}
		lines.addAll(table.build());
		lines.add("");
		
		// can be faster to loop over site in outer loop, rupture in inner loop
		Map<Site, List<RuptureComparison<E>>> siteCompsMap = new HashMap<>();
		for (RuptureComparison<E> comp : comps) {
			for (Site site : comp.getApplicableSites()) {
				List<RuptureComparison<E>> mySiteComps = siteCompsMap.get(site);
				if (mySiteComps == null) {
					mySiteComps = new ArrayList<>();
					siteCompsMap.put(site, mySiteComps);
				}
				mySiteComps.add(comp);
			}
		}
		
		lines.add("### Period-Dependent Residual Components");
		lines.add(topLink); lines.add("");
		lines.add("**Note: These are not corrected for covariance. "
				+ "Currently only useful for comparing relative phi and tau, not absolute values**");
		lines.add("");
//		lines.add("Phi is calculated as square root of the mean of the standard deviation");
		System.out.println("Calculating residual components");
		if (replotResiduals || !new File(resourcesDir, "period_residual_components.png").exists())
			ResidualScatterPlot.plotPeriodDependentSigma(resourcesDir, "period_residual_components", siteCompsMap, simProv, false, imts);
		lines.add("![Residual Components]("+resourcesDir.getName()+"/period_residual_components.png)");
		lines.add("");
		
		lines.add("### Detrended Period-Dependent Residual Components");
		lines.add(topLink); lines.add("");
		lines.add("**Note: These are not yet corrected for covariance. "
				+ "Currently only useful for comparing relative phi and tau, not absolute values**");
		lines.add("");
		lines.add("Residuals shown here are first detrended according to the following magnitude & log-distance dependent average residuals");
		lines.add("");
		
		table = MarkdownUtils.tableBuilder().initNewLine();
		for (IMT imt : imts)
			table.addColumn("**"+imt.getDisplayName()+"**");
		table.finalizeLine().initNewLine();
		for (IMT imt : imts)
			table.addColumn("![Detrend XYZ]("+resourcesDir.getName()+"/detrend_residuals_"+imt.getPrefix()+".png)");
		table.finalizeLine();
		table.initNewLine();
		for (IMT imt : imts)
			table.addColumn("![Detrend Std Dev XYZ]("+resourcesDir.getName()+"/detrend_std_dev_"+imt.getPrefix()+".png)");
		table.finalizeLine();
		table.initNewLine();
		for (IMT imt : imts)
			table.addColumn("[Download CSV]("+resourcesDir.getName()+"/detrend_residuals_"+imt.getPrefix()+".csv)");
		table.finalizeLine();
		lines.addAll(table.build());
		lines.add("");
		if (doRakeBinned) {
			boolean any = false;
			table = MarkdownUtils.tableBuilder().initNewLine();
			for (IMT imt : imts)
				table.addColumn("**"+imt.getDisplayName()+"**");
			table.finalizeLine();
			for (int r=0; r<rakeFilters.size(); r++) {
				RuptureComparisonFilter<E> rakeFilter = rakeFilters.get(r);
				Map<Site, List<RuptureComparison<E>>> siteRakeCompsMap = new HashMap<>();
				for (Site site : siteCompsMap.keySet()) {
					List<RuptureComparison<E>> siteRakeComps = rakeFilter.getMatches(siteCompsMap.get(site), site);
					if (siteRakeComps != null && !siteRakeComps.isEmpty())
						siteRakeCompsMap.put(site, siteRakeComps);
				}
				if (!siteRakeCompsMap.isEmpty()) {
					System.out.println("Calculating detrended residual components for "+rakeLabels.get(r));
					String prefix = "detrend_"+rakeFileLabels.get(r);
					boolean needToPlot = replotResiduals;
					any = true;
					table.initNewLine().addColumn("__"+rakeLabels.get(r)+"__");
					for (int i=1; i<imts.length; i++)
						table.addColumn("");
					table.finalizeLine().initNewLine();
					for (IMT imt : imts) {
						String fName = prefix+"_residuals_"+imt.getPrefix()+".png";
						needToPlot |= !(new File(resourcesDir, fName).exists());
						table.addColumn("![Detrend XYZ]("+resourcesDir.getName()+"/"+fName+")");
					}
					table.finalizeLine();
					table.initNewLine();
					for (IMT imt : imts) {
						String fName = prefix+"_std_dev_"+imt.getPrefix()+".png";
						needToPlot |= !(new File(resourcesDir, fName).exists());
						table.addColumn("![Detrend Std Dev XYZ]("+resourcesDir.getName()+"/"+fName+")");
					}
					table.finalizeLine();
					table.initNewLine();
					for (IMT imt : imts) {
						String fName = prefix+"_residuals_"+imt.getPrefix()+".csv";
						needToPlot |= !(new File(resourcesDir, fName).exists());
						table.addColumn("[Download CSV]("+resourcesDir.getName()+"/"+fName+")");
					}
					table.finalizeLine();
					
					if (needToPlot) {
						List<CSVFile<String>> csvs = new ArrayList<>();
						EvenlyDiscrXYZ_DataSet[][] detrends = ResidualScatterPlot.calcDetrendResiduals(siteRakeCompsMap, simProv, csvs, imts);
						ResidualScatterPlot.plotRedidualScatters(resourcesDir, prefix, imts, detrends[0], detrends[1]);
						Preconditions.checkState(csvs.size() == imts.length);
						for (int p=0; p<imts.length; p++)
							csvs.get(p).writeToFile(new File(resourcesDir, prefix+"_residuals_"+imts[p].getPrefix()+".csv"));
					}
				}
			}
			if (any) {
				lines.add("#### Rake-Specific Residual Scatter Plots");
				lines.add(topLink); lines.add("");
				lines.addAll(table.build());
				lines.add("");
			}
		}
		System.out.println("Calculating detrended residual components");
		if (replotResiduals || !new File(resourcesDir, "period_residual_detrend_components.png").exists()
				|| !new File(resourcesDir, "detrend_std_dev_"+imts[0].getPrefix()+".png").exists())
			ResidualScatterPlot.plotPeriodDependentSigma(resourcesDir, "period_residual_detrend_components", siteCompsMap, simProv, true, imts);
		lines.add("![Residual Components]("+resourcesDir.getName()+"/period_residual_detrend_components.png)");
		lines.add("");
		
		// for consistent iteration order
		List<Site> compSitesList = new ArrayList<>(siteCompsMap.keySet());
		for (ResidualType type : residualTypes) {
			System.out.println("Calculating GMPE residuals for: "+type.name);
			lines.add("### GMPE "+type.name+" Residuals");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (IMT imt : imts)
				table.addColumn("**"+imt.getDisplayName()+"**");
			table.finalizeLine();
			
			String[] residualPrefixes = new String[imts.length];
			boolean hasAllPlots = true;
			for (int i=0; i<imts.length; i++) {
				residualPrefixes[i] = "gmpe_residuals_"+type.name()+"_"+imts[i].getPrefix();
				hasAllPlots = hasAllPlots && new File(resourcesDir, residualPrefixes[i]+"_scatter.png").exists()
						&& new File(resourcesDir, residualPrefixes[i]+"_hist2d.png").exists();
			}
			
			ResidualScatterPlot[] residualPlots;
			if (!hasAllPlots || replotResiduals) {
				residualPlots = new ResidualScatterPlot[imts.length];
				PrimitiveArrayXY_Dataset prevScatter = null;
				for (int i=0; i<imts.length; i++) {
					IMT imt = imts[i];
					PrimitiveArrayXY_Dataset scatter;
					if (prevScatter == null)
						scatter = new PrimitiveArrayXY_Dataset(regressValCounts.get(type));
					else
						// avoids storing x array in memory twice, as it is identical
						scatter = new PrimitiveArrayXY_Dataset(prevScatter.getXArray());
					for (Site site : compSitesList) {
						for (RuptureComparison<E> comp : siteCompsMap.get(site)) {
							double x = type.getValue(comp, site);
							if (!Double.isNaN(type.minX) && x < type.minX)
								continue;
							if (type.log)
								Preconditions.checkState(x >= 0d);
							
							double gmpeVal = comp.getLogMean(site, imt);
							for (double simVal : simProv.getValues(site, comp.getRupture(), imt)) {
								// in log space
								simVal = Math.log(simVal);

								double val = (simVal - gmpeVal);
								scatter.set(x, val);
							}
						}
					}
					prevScatter = scatter;
					String xAxisLabel = type.name;
					if (type.units != null)
						xAxisLabel += " ("+type.units+")";
					ScalarIMR gmpe = checkOutGMPE(gmpeRef);
					gmpe.setParamDefaults();
					for (ResidualType t2 : residualTypes) {
						if (t2 == type || t2.parameterName == null)
							continue;
						double meanVal = residualDefaults.get(t2);
						if (gmpe instanceof MultiIMR_Averaged_AttenRel) {
							((MultiIMR_Averaged_AttenRel)gmpe).setParameterInIMRs(t2.parameterName, meanVal);
						} else {
							Parameter<Double> gmpeParameter = gmpe.getParameter(t2.parameterName);
							if (gmpeParameter instanceof WarningDoubleParameter)
								((WarningDoubleParameter)gmpeParameter).setValueIgnoreWarning(meanVal);
							else
								gmpeParameter.setValue(meanVal);
						}
					}
					imt.setIMT(gmpe);
					residualPlots[i] = new ResidualScatterPlot(scatter, xAxisLabel, type.log, type.deltaX, residualLabel,
							imt.getDisplayName()+" "+type.name+" Residuals", gmpe, type.parameterName);
					if (gmpe != null)
						checkInGMPE(gmpeRef, gmpe);
					residualPlots[i].setWritePDF(false);
//					residualPlots[i].setPlotLinearFit(false);
				}
			} else {
				residualPlots = null;
			}
			
			table.initNewLine();
			for (int i=0; i<imts.length; i++) {
				String prefix = residualPrefixes[i]+"_scatter";
				
				if (residualPlots != null) {
					System.out.println("Plotting GMPE residual: "+prefix);
					residualPlots[i].plotScatter(resourcesDir, prefix);
				}
				
				File plotFile = new File(resourcesDir, prefix+".png");
				Preconditions.checkState(plotFile.exists());
				table.addColumn("![Scatter]("+resourcesDir.getName()
					+"/"+plotFile.getName()+")");
			}
			table.finalizeLine();
			table.initNewLine();
			for (int i=0; i<imts.length; i++) {
				String prefix = residualPrefixes[i]+"_hist2d";
				
				if (residualPlots != null) {
					System.out.println("Plotting GMPE 2-D residual: "+prefix);
					residualPlots[i].plot2DHist(resourcesDir, prefix);
				}
				
				File plotFile = new File(resourcesDir, prefix+".png");
				Preconditions.checkState(plotFile.exists());
				table.addColumn("![2-D Hist]("+resourcesDir.getName()
					+"/"+plotFile.getName()+")");
			}
			table.finalizeLine();
			lines.add("");
			lines.addAll(table.wrap(max_table_fig_columns, 0).build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private boolean doesParameterVary(String paramName, Collection<Site> sites) {
		boolean allNull = true;
		Double prevVal = null;
		boolean first = true;
		for (Site site : sites) {
			Double val = site.getParameter(Double.class, paramName).getValue();
			allNull = allNull && val == null;
			if (!first && !allNull) {
				if (val == null)
					// previously not null, but null now
					return true;
				if (!val.equals(prevVal))
					return true;
			}
			first = false;
			prevVal = val;
		}
		return false;
	}
	
	private enum ResidualType {
		MAG("Magnitude", null, MagParam.NAME, 0.1) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return comp.getMagnitude();
			}
		},
		DIST_RUP("rRup", "km", DistanceRupParameter.NAME, true, 10d, Double.NaN) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return comp.getDistanceRup(site);
			}
		},
		DIST_JB("rJB", "km", DistanceJBParameter.NAME, true, 10d, Double.NaN) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return comp.getDistanceJB(site);
			}
		},
		VS30("Vs30", "m/s", Vs30_Param.NAME, true, 200d, Double.NaN) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return site.getParameter(Double.class, Vs30_Param.NAME).getValue();
			}
		},
		Z10("Z10", "m", DepthTo1pt0kmPerSecParam.NAME) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return site.getParameter(Double.class, DepthTo1pt0kmPerSecParam.NAME).getValue();
			}
		},
		Z25("Z25", "km", DepthTo2pt5kmPerSecParam.NAME) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return site.getParameter(Double.class, DepthTo2pt5kmPerSecParam.NAME).getValue();
			}
		},
		OCCUR_TIME("Occurrence Time", "yrs", null) {
			@Override
			public Double getValue(RuptureComparison<?> comp, Site site) {
				return comp.getRuptureTimeYears();
			}
		};

		private String name;
		private String units;
		private String parameterName;
		private boolean log;
		private double minX;
		private double deltaX;
		
		private ResidualType(String name, String units, String parameterName) {
			this(name, units, parameterName, false, Double.NEGATIVE_INFINITY, Double.NaN);
		}

		private ResidualType(String name, String units, String parameterName, double deltaX) {
			this(name, units, parameterName, false, Double.NEGATIVE_INFINITY, deltaX);
		}

		private ResidualType(String name, String units, String parameterName, boolean log, double minX, double deltaX) {
			this.name = name;
			this.units = units;
			this.parameterName = parameterName;
			this.log = log;
			this.minX = minX;
			this.deltaX = deltaX;
		}
		
		public abstract Double getValue(RuptureComparison<?> comp, Site site);
	}
	
	protected abstract double calcRupAzimuthDiff(E event, int simIndex, Site site);
	
	protected static double calcRupAzimuthDiff(Location rupStart, Location rupEnd, Location centroid, Location hypo,
			Location siteLoc) {
		double hypoDist1 = LocationUtils.horzDistanceFast(hypo, rupStart);
		double hypoDist2 = LocationUtils.horzDistanceFast(hypo, rupEnd);
		if (hypoDist2 < hypoDist1) {
			// flip so that hypocenter is closest to loc1
			Location t = rupStart;
			rupStart = rupEnd;
			rupEnd = t;
		}
		
		double refAz = LocationUtils.azimuth(rupStart, rupEnd);
		double locAz = LocationUtils.azimuth(centroid, siteLoc);
		double diff1 = Math.abs(locAz - refAz);
		double diff2 = Math.abs((360+locAz) - refAz);
		double diff3 = Math.abs(locAz - (360+refAz));
		double diff = Math.min(diff1, Math.min(diff2, diff3));
		Preconditions.checkState(diff <= 180);
		return diff;
	}
	
	private synchronized double getRupAzimuthDiff(E event, int simIndex, Site site) {
		Map<Integer, Double> azMap = rupSiteAzMap.get(event, site);
		if (azMap == null) {
			azMap = new HashMap<>();
			rupSiteAzMap.put(event, site, azMap);
		}
		Double az = azMap.get(simIndex);
		if (az == null) {
			az = calcRupAzimuthDiff(event, simIndex, site);
			azMap.put(simIndex, az);
		}
		return az;
	}
	
	public boolean plotAggregateRotDRatio(Range magRange, File outputDir, String prefix) throws IOException {
		List<DiscretizedFunc[]> ratios = new ArrayList<>();
		
		for (Site site : sites) {
			for (E rupture : simProv.getRupturesForSite(site)) {
				if (magRange != null && !magRange.contains(simProv.getMagnitude(rupture)))
					continue;
				
				for (DiscretizedFunc ratio : simProv.getRotDRatios(site, rupture)) {
					ratios.add(new DiscretizedFunc[] {ratio});
				}
			}
		}
		if (ratios.isEmpty())
			return false;
		
		String title = "RotD100/50 Ratio";
		if (magRange == null)
			title += ", All Mags";
		else
			title += ", "+optionalDigitDF.format(magRange.getLowerBound())
				+"≤M≤"+optionalDigitDF.format(magRange.getUpperBound());
		
		SpectraPlotter.plotRotDRatio(ratios, simName, title, outputDir, prefix);
		
		return true;
	}
	
	private enum RatioDependence {
		DIST("Distance Dependence", "Distance", false, false),
		AZIMUTH("Azimuth Dependence", "Hypocentral Azimuth", false, false),
		ROTD50("RotD50 Dependence", "RotD50", true, true),
		GMPE_Z("Z-Value Dependence", "GMPE Z", true, false);
		
		private String plotTitle;
		private String axisLabel;
		private boolean periodDependent;
		private boolean logX;

		private RatioDependence(String plotTitle, String axisLabel, boolean periodDependent, boolean logX) {
			this.plotTitle = plotTitle;
			this.axisLabel = axisLabel;
			this.periodDependent = periodDependent;
			this.logX = logX;
		}
		
		public boolean isPeriodDependent() {
			return periodDependent;
		}
		
		public String getPlotTitle() {
			return plotTitle;
		}
		
		public String getAxisLabel() {
			return axisLabel;
		}
		
		public boolean isLogX() {
			return logX;
		}
	}
	
	public boolean plotRotDRatioDependence(Range magRange, Range distRange, double[] periods, boolean scatter, int numPoints,
			File outputDir, String prefix, RatioDependence quantity, List<? extends RuptureComparison<E>> gmpeComps)
					throws IOException {
		
		List<DiscretizedFunc[]> rotDs = new ArrayList<>();
		List<List<Double>> scalars = new ArrayList<>();
		if (quantity.isPeriodDependent()) {
			// individual ones
			for (int i=0; i<periods.length; i++)
				scalars.add(new ArrayList<>());
		} else {
			// share it
			List<Double> vals = new ArrayList<>();
			for (int i=0; i<periods.length; i++)
				scalars.add(vals);
		}
		
		for (Site site : sites) {
			for (RuptureComparison<E> comp : new RuptureComparisonFilter.HasSiteFilter<E>().getMatches(gmpeComps, site)) {
				if (magRange != null && !magRange.contains(comp.getMagnitude()))
					continue;
				
				double dist = Double.NaN;
				if (distRange != null || quantity == RatioDependence.DIST) {
					if (distJB)
						dist = comp.getDistanceJB(site);
					else
						dist = comp.getDistanceRup(site);
					if (distRange != null && !distRange.contains(dist))
						continue;
				}
				
				for (int index=0; index<simProv.getNumSimulations(site, comp.getRupture()); index++) {
					DiscretizedFunc ratio = simProv.getRotDRatio(site, comp.getRupture(), index);
					DiscretizedFunc rd50Func = null;
					if (quantity == RatioDependence.ROTD50 || quantity == RatioDependence.GMPE_Z)
						rd50Func = simProv.getRotD50(site, comp.getRupture(), index);
					
					switch (quantity) {
					case DIST:
						scalars.get(0).add(dist);
						break;
					case AZIMUTH:
						scalars.get(0).add(getRupAzimuthDiff(comp.getRupture(), index, site));
						break;
					case ROTD50:
						for (int p=0; p<periods.length; p++)
							scalars.get(p).add(rd50Func.getInterpolatedY(periods[p]));
						break;
					case GMPE_Z:
						for (int p=0; p<periods.length; p++) {
							double simVal = Math.log(rd50Func.getInterpolatedY(periods[p]));
							IMT imt = IMT.forPeriod(periods[p]);
							double gmpeVal = comp.getLogMean(site, imt);
							double gmpeSigma = comp.getStdDev(site, imt);
							double z = (simVal - gmpeVal)/gmpeSigma;
							scalars.get(p).add(z);
						}
						break;

					default:
						throw new IllegalStateException("unknown quanitity: "+quantity);
					}
					
					rotDs.add(new DiscretizedFunc[] {ratio});
				}
			}
		}
		if (rotDs.isEmpty())
			return false;
		
		String title = "RotD100/50 "+quantity.getPlotTitle();
		
		String scalarLabel = quantity.getAxisLabel();
		if (quantity == RatioDependence.DIST)
			scalarLabel = getDistShortName()+" (km)";
		
		if (magRange == null)
			title += ", All Mags";
		else
			title += ", "+optionalDigitDF.format(magRange.getLowerBound())
				+"≤M≤"+optionalDigitDF.format(magRange.getUpperBound());
		
		if (scatter) {
			for (int i=0; i<periods.length; i++) {
				double period = periods[i];
				String myPrefix = prefix+"_"+optionalDigitDF.format(period)+"s";
				SpectraPlotter.plotRotDRatioScatter(rotDs, scalars.get(i), scalarLabel, period,
						simName, title, outputDir, myPrefix, quantity.isLogX(), numPoints);
			}
		} else {
			SpectraPlotter.plotRotDRatioPeriodDependence(rotDs, scalars, scalarLabel, numPoints, periods,
					simName, title, outputDir, prefix, quantity.isLogX());
		}
		
		return true;
	}
	
	public void generateRotDRatioPage(File outputDir, List<String> headerLines, double[] aggregatedPeriods, double[] scatterPeriods,
			AttenRelRef gmpeRef, List<? extends RuptureComparison<E>> gmpeComps) throws IOException {
		Preconditions.checkState(simProv.hasRotD100());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		String distShortName = getDistShortName();
		
		// header
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.getLast().isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Site List");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Name", "Location", "# Events");
		for (Site s : this.sites) {
			table.initNewLine();
			table.addColumn(s.getName());
			Location loc = s.getLocation();
			table.addColumn("*"+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
			table.addColumn(simProv.getRupturesForSite(s).size()+"");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		
		Range totalDistRange = new Range(0d, cutoffDist);
		
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
					magRange, totalDistRange, aggregatedPeriods, false, 20, resourcesDir, prefix, RatioDependence.DIST, gmpeComps));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", Azimuth Dependence");
			lines.add(topLink); lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio azimuth dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_az_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, aggregatedPeriods, false, 40, resourcesDir, prefix, RatioDependence.AZIMUTH, gmpeComps));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", RotD50 Dependence");
			lines.add(topLink); lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio amplitude dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_amplitude_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, aggregatedPeriods, false, 20, resourcesDir, prefix, RatioDependence.ROTD50, gmpeComps));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("### "+magLabel+", GMPE Z Dependence");
			lines.add(topLink); lines.add("");
			lines.add("**z-score**: (ln(*"+simName+"*) - ln(*GMPE-mean*)) / *GMPE-sigma*");
			lines.add("");
			lines.add("**GMPE**: "+gmpeRef.getName());
			lines.add("");
			System.out.println("Plotting "+magLabel+" RotD ratio GMPE Z dependence");
			prefix = "rot_d_ratio_"+magFileLabel+"_z_dependence";
			Preconditions.checkState(plotRotDRatioDependence(
					magRange, totalDistRange, aggregatedPeriods, false, 20, resourcesDir, prefix, RatioDependence.GMPE_Z, gmpeComps));
			lines.add("![RotD Ratio]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			lines.add("#### "+magLabel+", GMPE Z Dependence Scatters");
			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("Distance");
			for (double period : scatterPeriods)
				table.addColumn("**"+optionalDigitDF.format(period)+"s**");
			table.finalizeLine();
			for (int d=0; d<distRanges.size(); d++) {
				Range distRange = distRanges.get(d);
				String distLabel = distLabels.get(d);
				String distFileLabel = distFileLabels.get(d);
				
				table.initNewLine().addColumn("**"+distLabel+"**");
				
				prefix = "rot_d_ratio_"+magFileLabel+"_"+distFileLabel+"_scatter";
				Preconditions.checkState(plotRotDRatioDependence(magRange, distRange, scatterPeriods, true, 20, resourcesDir, prefix,
						RatioDependence.GMPE_Z, gmpeComps));
				for (double period : scatterPeriods) {
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
	
	protected static Region bufferRegion(Region reg, double bufferKM) {
		Location topLeft = new Location(reg.getMaxLat(), reg.getMinLon());
		Location botRight = new Location(reg.getMinLat(), reg.getMaxLon());
		
		return new Region(LocationUtils.location(topLeft, 1.75*Math.PI, bufferKM),
				LocationUtils.location(botRight, 0.75*Math.PI, bufferKM));
	}
	
	private static enum HazardRatioType {
		PDIFF,
		RATIO,
		LOG_RATIO
	}
	
	public void generateNonErgodicMapPage(File outputDir, List<String> headerLines, List<? extends FaultSection> subSects,
			Region siteRegion, Region sourceRegion, AttenRelRef gmpeRef, List<? extends RuptureComparison<E>> comps,
			Map<E, List<FaultSection>> rupSectMappings, Map<E, FaultSection> rupNuclSects, IMT[] imts,
			List<Site> highlightSites) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		HazardRatioType ratioType = HazardRatioType.PDIFF;
		double linearRatioFactor = 2d;
		
		lines.add("# Non-Ergodic Source & Site Maps");
		lines.add("");
		lines.add("This page compares non-ergodic simulation results against an ergodic GMPE, primarily through z-score"
				+ " analysis of residuals. The comparison fully ergodic GMPE is _"+gmpeRef.getName()+"_ and results are"
				+ " shown relative to the predicted mean and standard deviations from the model.");
		lines.add("");
		
		// header
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.getLast().isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		GeographicMapMaker siteMapMaker = new RupSetMapMaker(subSects, siteRegion);
		siteMapMaker.setWriteGeoJSON(false);
//		siteMapMaker.setScatterSymbol(PlotSymbol.FILLED_TRIANGLE, 7f, PlotSymbol.TRIANGLE, new Color(0, 0, 0, 127));
		siteMapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 10f, PlotSymbol.CIRCLE, new Color(0, 0, 0, 127));
		GeographicMapMaker sourceMapMaker = new RupSetMapMaker(subSects, sourceRegion);
		sourceMapMaker.setSkipNaNs(true);
		
		CPT zScoreCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1d, 1d);
		CPT sigmaCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0.5d, 1d);
		
		CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1d, 1d);
		if (ratioType == HazardRatioType.RATIO) {
			CPT leftCPT = new CPT();
			CPT rightCPT = new CPT();
			for (CPTVal val : ratioCPT) {
				if (val.start < 0f)
					leftCPT.add(val);
				else
					rightCPT.add(val);
			}
			leftCPT = leftCPT.rescale(1d/linearRatioFactor, 1d);
			rightCPT = rightCPT.rescale(1d, linearRatioFactor);
			ratioCPT = new CPT();
			ratioCPT.addAll(leftCPT);
			ratioCPT.addAll(rightCPT);
			ratioCPT.setBelowMinColor(ratioCPT.getMinColor());
			ratioCPT.setAboveMaxColor(ratioCPT.getMaxColor());
		} else if (ratioType == HazardRatioType.PDIFF) {
			ratioCPT = ratioCPT.rescale(-100d, 100d);
		} else {
			Preconditions.checkState(ratioType == HazardRatioType.LOG_RATIO);
		}
		CPT rawHazardCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
		
		List<Site> sites = new ArrayList<>(this.sites);
		for (int i=sites.size(); --i>=0;)
			if (sites.get(i) == null)
				sites.remove(i);
		List<Location> siteLocs = new ArrayList<>();
		for (Site site : sites)
			siteLocs.add(site.getLocation());
		
		List<Future<ZScoreResult[]>> scoreFutures = new ArrayList<>();
		System.out.println("Calculating site z-scores");
		for (Site site : sites) {
			List<Site> mySites = new ArrayList<>();
			mySites.add(site);
			
			MatchesSiteFilter<E> siteFilter = new RuptureComparisonFilter.MatchesSiteFilter<E>(site);
			
			List<? extends RuptureComparison<E>> siteComps = siteFilter.getMatches(comps, site);
			scoreFutures.add(exec.submit(new Callable<ZScoreResult[]>() {

				@Override
				public ZScoreResult[] call() throws Exception {
					return ZScoreHistPlot.calcZScores(simProv, siteComps, mySites, imts, siteFilter);
				}
			}));
		}
		List<ZScoreResult[]> siteScores = new ArrayList<>();
		for (Future<ZScoreResult[]> future : scoreFutures) {
			try {
				siteScores.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			System.out.println("Done with site "+siteScores.size()+"/"+scoreFutures.size());
		}
		
		System.out.println("Calculating sect z-scores");
		scoreFutures = new ArrayList<>();
		int numSectsInside = 0;
		boolean[] sectInsides = new boolean[subSects.size()];
		for (FaultSection sect : subSects) {
			for (Location loc : sect.getFaultTrace()) {
				if (sourceRegion.contains(loc)) {
					numSectsInside++;
					sectInsides[sect.getSectionId()] = true;
					break;
				}
			}
		}
		
		if (rupNuclSects != null) {
			lines.add("## Source Nucleation");
			lines.add(topLink); lines.add("");
			
			lines.add("The left plot gives the fraction of ruptures that nucleate on a given section, relative to the "
					+ "total number of ruptures for which that section participates. The right plot gives the primary "
					+ "propagation direction for ruptures involving each section. Colors are dark when all ruptures "
					+ "follow the plotted propagation direction, and white when each direction is equally utilized.");
			lines.add("");
			
			Map<FaultSection, List<E>> sectRups = new HashMap<>();
			for (E rup : rupSectMappings.keySet()) {
				for (FaultSection sect : rupSectMappings.get(rup)) {
					if (!sectRups.containsKey(sect))
						sectRups.put(sect, new ArrayList<>());
					sectRups.get(sect).add(rup);
				}
			}
			
			List<Color> directionColors = new ArrayList<>();
			List<Double> nuclFracts = new ArrayList<>();
			CPT nuclCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0d, 1d);
			CPT dirCPT = new CPT();
			Color prevColor = null;
			double bright = 0.9d;
			for (int i=0; i<=360; i++) {
				double hue = (double)i/360d;
				Color color = Color.getHSBColor((float)hue, 1f, (float)bright);
				if (prevColor != null)
					dirCPT.add(new CPTVal(i-1, prevColor, i, color));
				prevColor = color;
			}
			for (int s=0; s<subSects.size(); s++) {
				FaultSection sect = subSects.get(s);
				List<E> rups = sectRups.get(sect);
				if (!sectInsides[s] || rups == null) {
					// not in our map region, or no ruptures
					directionColors.add(null);
					nuclFracts.add(Double.NaN);
					continue;
				}
				Location first = sect.getFaultTrace().first();
				Location last = sect.getFaultTrace().last();
				
				int numParticipations = 0;
				int numNucleations = 0;
				int numForwards = 0;
				int numBackwards = 0;
				for (E rup : rups) {
					numParticipations++;
					FaultSection nuclSect = rupNuclSects.get(rup);
					Preconditions.checkNotNull(nuclSect, "not all ruptures have nucleation sections");
					if (nuclSect.equals(sect)) {
						numNucleations++;
					} else {
						// see what way it propagated
						double distToFirst = Double.POSITIVE_INFINITY;
						double distToLast = Double.POSITIVE_INFINITY;
						for (Location loc : nuclSect.getFaultTrace()) {
							distToFirst = Math.min(distToFirst, LocationUtils.horzDistanceFast(first, loc));
							distToLast = Math.min(distToLast, LocationUtils.horzDistanceFast(last, loc));
						}
						if (distToFirst < distToLast)
							// nucleation section is closer to the first point on this fault trace, propagated forwards
							numForwards++;
						else
							// opposite
							numBackwards++;
					}
				}
				nuclFracts.add((double)numNucleations/(double)numParticipations);
				Color dirColor;
				if (numForwards == numBackwards) {
					dirColor = Color.WHITE;
				} else {
					double az;
					int numExcessInDir;
					if (numForwards > numBackwards) {
						numExcessInDir = numForwards-numBackwards;
						az = LocationUtils.azimuth(first, last);
					} else {
						numExcessInDir = numBackwards-numForwards;
						az = LocationUtils.azimuth(last, first);
					}
					double blendFract = (double)numExcessInDir/(double)numParticipations;
					dirColor = Color.getHSBColor((float)(az/360d), (float)blendFract, (float)(bright*blendFract + (1-blendFract)));
				}
				directionColors.add(dirColor);
			}
			sourceMapMaker.plotSectScalars(nuclFracts, nuclCPT, "Nucleation/Participation Ratio");
			sourceMapMaker.plot(resourcesDir, "nucl_partic_ratios", " ", 900);
			sourceMapMaker.clearSectScalars();
			sourceMapMaker.plotSectColors(directionColors, dirCPT, "Primary Propagation Direction (deg)");
			
			PlotSpec spec = sourceMapMaker.buildPlot(" ");
			// build color wheel annotation
			
			// figure out aspect ratio of lat to lon
			Location centerLoc = new Location(0.5*(sourceRegion.getMinLat()+sourceRegion.getMaxLat()),
					0.5*(sourceRegion.getMinLon()+sourceRegion.getMaxLon()));
			double latLen = LocationUtils.horzDistance(centerLoc, new Location(centerLoc.getLatitude()+1, centerLoc.getLongitude()));
			double lonLen = LocationUtils.horzDistance(centerLoc, new Location(centerLoc.getLatitude(), centerLoc.getLongitude()+1));
			double aspect = latLen/lonLen;
			
			// center location of the wheel
			double centerOffset = 0.12*(sourceRegion.getMaxLon()-sourceRegion.getMinLon());
			double wheelCenterX = sourceRegion.getMinLon() + centerOffset;
			double wheelCenterY = sourceRegion.getMinLat() + centerOffset/aspect;
//			double widthX = 0.20*(sourceRegion.getMaxLat()-sourceRegion.getMinLat());
//			double aspect = PlotUtils.calcAspectRatio(
//					new Range(sourceRegion.getMinLon(), sourceRegion.getMaxLon()),
//					new Range(sourceRegion.getMinLat(), sourceRegion.getMaxLat()), true);
			
			double innerMult = 0.3;
			double outerMult = 0.4;
//			double sumY = Math.max(1d, hist.calcSumOfY_Vals());
			double halfDelta = 15d;
			for (double centerAz=0d; centerAz<=350d; centerAz+=30d) {
				double startAz = Math.toRadians(centerAz-halfDelta);
				double endAz = Math.toRadians(centerAz+halfDelta);
				
				List<Point2D> points = new ArrayList<>();
				
				double startX = Math.sin(startAz);
				double startY = Math.cos(startAz);
				double endX = Math.sin(endAz);
				double endY = Math.cos(endAz);
				
				points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
				points.add(new Point2D.Double(outerMult*startX, outerMult*startY));
				points.add(new Point2D.Double(outerMult*endX, outerMult*endY));
				points.add(new Point2D.Double(innerMult*endX, innerMult*endY));
				points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
				
				double[] polygon = new double[points.size()*2];
				int cnt = 0;
				for (Point2D pt : points) {
					polygon[cnt++] = pt.getX()+wheelCenterX;
					polygon[cnt++] = pt.getY()/aspect+wheelCenterY;
				}
				Color color = dirCPT.getColor((float)centerAz);
				
				Stroke stroke = PlotLineType.SOLID.buildStroke(1f);
				spec.addPlotAnnotation(new XYPolygonAnnotation(polygon, stroke, Color.DARK_GRAY, color));
			}
			sourceMapMaker.plot(resourcesDir, "sect_prop_az", spec, 900);
//			sourceMapMaker.plot(resourcesDir, "sect_prop_az", " ", 900);
			sourceMapMaker.clearSectColors();
//			System.out.println("ASPECT = "+latLen+" / "+lonLen+" = "+aspect);
//			System.out.flush();
//			System.exit(0);
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			table.initNewLine();
			table.addColumn("![Section Nucleation]("+resourcesDir.getName()+"/nucl_partic_ratios.png)");
			table.addColumn("![Section Propagation]("+resourcesDir.getName()+"/sect_prop_az.png)");
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		// pre-sort by section
		List<List<RuptureComparison<E>>> sectComps = new ArrayList<>();
		for (int s=0; s<sectInsides.length; s++) {
			if (sectInsides[s])
				sectComps.add(new ArrayList<>());
			else
				sectComps.add(null);
		}
		for (RuptureComparison<E> comp : comps) {
			for (FaultSection sect : rupSectMappings.get(comp.getRupture())) {
				if (sectInsides[sect.getSectionId()])
					sectComps.get(sect.getSectionId()).add(comp);
			}
		}
		
		for (FaultSection sect : subSects) {
			if (sectInsides[sect.getSectionId()]) {
				scoreFutures.add(exec.submit(new Callable<ZScoreResult[]>() {

					@Override
					public ZScoreResult[] call() throws Exception {
						return ZScoreHistPlot.calcZScores(simProv, sectComps.get(sect.getSectionId()), sites, imts, null);
					}
				}));
			} else {
				scoreFutures.add(null);
			}
		}
		
		List<ZScoreResult[]> sectScores = new ArrayList<>();
		int numSectsProcessed = 0;
		for (Future<ZScoreResult[]> future : scoreFutures) {
			if (future == null) {
				sectScores.add(null);
			} else {
				try {
					sectScores.add(future.get());
					numSectsProcessed++;
					System.out.println("Done with sect "+numSectsProcessed+"/"+numSectsInside);
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		}
		
		SimulationHazardCurveCalc<E> simCurveCalc = new SimulationHazardCurveCalc<>(simProv);
		
		for (int p=0; p<imts.length; p++) {
			IMT imt = imts[p];
			
			lines.add("## "+imt.getDisplayName());
			lines.add(topLink); lines.add("");
			
			lines.add("### "+imt.getDisplayName()+" z-scores");
			lines.add(topLink); lines.add("");
			
			List<Double> siteZMeans = new ArrayList<>();
			List<Double> siteZStdDevs = new ArrayList<>();
			for (int i=0; i<siteScores.size(); i++) {
				ZScoreResult score = siteScores.get(i)[p];
				siteZMeans.add(score.mean);
				siteZStdDevs.add(score.stdDevFract);
			}
			
			String prefix = imt.getPrefix();
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			siteMapMaker.plotScatterScalars(siteLocs, siteZMeans, zScoreCPT, imt.getDisplayName()+" site z-scores");
			siteMapMaker.plot(resourcesDir, prefix+"_site_z_means", " ", 900);
			siteMapMaker.plotScatterScalars(siteLocs, siteZStdDevs, sigmaCPT, imt.getDisplayName()+" site σ-fracts");
			siteMapMaker.plot(resourcesDir, prefix+"_site_z_sds", " ", 900);
			
			table.initNewLine();
			table.addColumn("![Site z-scores]("+resourcesDir.getName()+"/"+prefix+"_site_z_means.png)");
			table.addColumn("![Site z-scores]("+resourcesDir.getName()+"/"+prefix+"_site_z_sds.png)");
			table.finalizeLine();
			
			siteMapMaker.clearScatters();
			
			// now histograms
			ZScoreResult zSiteScores = buildZResult(siteZMeans, true, true);
			
			ZScoreHistPlot.plotStandardNormal(new ZScoreResult[] {zSiteScores}, "Site z-scores", new IMT[] {imt},
					gmpeRef, null, resourcesDir, prefix+"_site_z_hist", 0, true, -1);
			sigFractHistPlot(resourcesDir, prefix+"_site_z_sds_hist", siteZStdDevs, "Site σ-fracts", " ");
			
			table.initNewLine();
			table.addColumn("![Site histogram]("+resourcesDir.getName()+"/"+prefix+"_site_z_hist.png)");
			table.addColumn("![Site histogram]("+resourcesDir.getName()+"/"+prefix+"_site_z_sds_hist.png)");
			table.finalizeLine();
			
			List<Double> sourceZMeans = new ArrayList<>();
			List<Double> sourceZStdDevs = new ArrayList<>();
			for (int s=0; s<subSects.size(); s++) {
				ZScoreResult[] scores = sectScores.get(s);
				if (scores == null) {
					sourceZMeans.add(Double.NaN);
					sourceZStdDevs.add(Double.NaN);
				} else {
					sourceZMeans.add(scores[p].mean);
					sourceZStdDevs.add(scores[p].stdDevFract);
				}
			}
			
			sourceMapMaker.setScatterSymbol(PlotSymbol.CIRCLE, 2f);
			sourceMapMaker.plotScatters(siteLocs, new Color(0, 0, 0, 100));
			sourceMapMaker.plotSectScalars(sourceZMeans, zScoreCPT, imt.getDisplayName()+" source/path z-scores");
			sourceMapMaker.plot(resourcesDir, prefix+"_source_z_means", " ", 900);
			sourceMapMaker.plotSectScalars(sourceZStdDevs, sigmaCPT, imt.getDisplayName()+" source/path σ-fracts");
			sourceMapMaker.plot(resourcesDir, prefix+"_source_z_sds", " ", 900);
			sourceMapMaker.clearScatters();
			
			table.initNewLine();
			table.addColumn("![Source z-scores]("+resourcesDir.getName()+"/"+prefix+"_source_z_means.png)");
			table.addColumn("![Source z-scores]("+resourcesDir.getName()+"/"+prefix+"_source_z_sds.png)");
			table.finalizeLine();
			
			// now histograms
			ZScoreResult zSourceScores = buildZResult(sourceZMeans, true, true);

			ZScoreHistPlot.plotStandardNormal(new ZScoreResult[] {zSourceScores}, "Source/path z-scores", new IMT[] {imt},
					gmpeRef, null, resourcesDir, prefix+"_source_z_hist", 0, true, -1);
			sigFractHistPlot(resourcesDir, prefix+"_source_z_sds_hist", sourceZStdDevs, "Source/path σ-fracts", " ");

			table.initNewLine();
			table.addColumn("![Source histogram]("+resourcesDir.getName()+"/"+prefix+"_source_z_hist.png)");
			table.addColumn("![Source histogram]("+resourcesDir.getName()+"/"+prefix+"_source_z_sds_hist.png)");
			table.finalizeLine();
			
			zMeanVsStdDevScatterPlot(resourcesDir, prefix+"_site_z_scatter", siteZMeans, siteZStdDevs,
					"Site z/σ-fract relationship");
			zMeanVsStdDevScatterPlot(resourcesDir, prefix+"_source_z_scatter", sourceZMeans, sourceZStdDevs,
					"Source/path z/σ-fract relationship");
			
			table.initNewLine();
			table.addColumn("![Site scatter]("+resourcesDir.getName()+"/"+prefix+"_site_z_scatter.png)");
			table.addColumn("![Source scatter]("+resourcesDir.getName()+"/"+prefix+"_source_z_scatter.png)");
			table.finalizeLine();
			
			lines.addAll(table.build());
			
			lines.add("### "+imt.getDisplayName()+" Hazard Comparisons");
			lines.add(topLink); lines.add("");
			
			System.out.println("Calculating "+imt.getDisplayName()+" hazard curves");
			List<Future<DiscretizedFunc[]>> curveFutures = new ArrayList<>();
			for (Site site : sites) {
				curveFutures.add(exec.submit(new Callable<DiscretizedFunc[]>() {

					@Override
					public DiscretizedFunc[] call() throws Exception {
						SimulationHazardPlotter<E> curvePlotter = new SimulationHazardPlotter<>(
								simCurveCalc, comps, site, 1d, gmpeRef);
						DiscretizedFunc gmpeCurve = curvePlotter.getCalcGMPECurve(imt);
						DiscretizedFunc simCurve = curvePlotter.getCalcSimCurve(simCurveCalc, imt);
						return new DiscretizedFunc[] { simCurve, gmpeCurve };
					}
				}));
			}
			
			List<DiscretizedFunc> simCurves = new ArrayList<>();
			List<DiscretizedFunc> gmpeCurves = new ArrayList<>();
			
			for (Future<DiscretizedFunc[]> future : curveFutures) {
				try {
					DiscretizedFunc[] curves = future.get();
					simCurves.add(curves[0]);
					gmpeCurves.add(curves[1]);
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			
			table = MarkdownUtils.tableBuilder();
			
			List<String> hazScatterPrefixes = new ArrayList<>();
			
			for (boolean rtgm : new boolean [] {false, true}) {
				List<Double> simVals = new ArrayList<>();
				List<Double> ratios = new ArrayList<>();
				double maxVal = 0d;
				
				for (int s=0; s<sites.size(); s++) {
					DiscretizedFunc simCurve = simCurves.get(s);
					DiscretizedFunc gmpeCurve = gmpeCurves.get(s);
					
					double simVal, gmpeVal;
					if (rtgm) {
						simVal =  SimulationHazardPlotter.calcRTGM(simCurve, 1d);
						gmpeVal =  SimulationHazardPlotter.calcRTGM(gmpeCurve, 1d);
					} else {
						simVal = simCurve.getFirstInterpolatedX_inLogXLogYDomain(4e-4);
						gmpeVal = gmpeCurve.getFirstInterpolatedX_inLogXLogYDomain(4e-4);
					}
					maxVal = Math.max(maxVal, simVal);
					
					simVals.add(simVal);
					double ratio;
					switch (ratioType) {
					case LOG_RATIO:
						ratio = Math.log10(simVal/gmpeVal);
						break;
					case RATIO:
						ratio = simVal/gmpeVal;
						break;
					case PDIFF:
						ratio = 100d*(simVal - gmpeVal)/gmpeVal;
						break;

					default:
						throw new IllegalStateException();
					}
					ratios.add(ratio);
				}
				
//				String logStr = ratioType == HazardRatioType.LOG_RATIO ? "Log10 " : "";
//				String ratioStr = ratioType == HazardRatioType.PDIFF ? ", % Difference" : ", Ratio";
				String lableSuffix = "Simulation vs "+gmpeRef.getShortName();
				switch (ratioType) {
				case LOG_RATIO:
					lableSuffix += ", Log10 Ratio";
					break;
				case RATIO:
					lableSuffix += ", Ratio";
					break;
				case PDIFF:
					lableSuffix += ", % Difference";
					break;

				default:
					throw new IllegalStateException();
				}
				String label, ratioLabel, myPrefix;
				if (rtgm) {
					label = imt.getDisplayName()+", RTGM ("+imt.getUnits()+")";
					ratioLabel = imt.getDisplayName()+", RTGM, "+lableSuffix;
					myPrefix = prefix+"_rtgm";
				} else {
					label = imt.getDisplayName()+", 2% in 50yr Hazard ("+imt.getUnits()+")";
					ratioLabel = imt.getDisplayName()+", 2% in 50yr, "+lableSuffix;
					myPrefix = prefix+"_2in50";
				}
				double ceilScalar;
				if (maxVal > 0.5)
					ceilScalar = 0.2;
				else
					ceilScalar = 0.1;
				CPT hazardCPT = rawHazardCPT.rescale(0d, 0.25*Math.ceil(maxVal*4d));
				siteMapMaker.plotScatterScalars(siteLocs, simVals, hazardCPT, label);
				siteMapMaker.plot(resourcesDir, myPrefix, " ", 900);
				siteMapMaker.plotScatterScalars(siteLocs, ratios, ratioCPT, ratioLabel);
				siteMapMaker.plot(resourcesDir, myPrefix+"_ratios", " ", 900);
				
				table.initNewLine();
				table.addColumn("![Hazard]("+resourcesDir.getName()+"/"+myPrefix+".png)");
				table.addColumn("![Hazard]("+resourcesDir.getName()+"/"+myPrefix+"_ratios.png)");
				table.finalizeLine();
				
				zMeanVsStdDevHazardScatterPlot(resourcesDir, myPrefix+"_z_scatter", siteZMeans, siteZStdDevs, ratios, ratioType,
						ratioCPT, ratioLabel, "Hazard, z, and σ-fracts relationship");
				hazScatterPrefixes.add(myPrefix+"_z_scatter");
			}
			
			table.initNewLine();
			for (String myPrefix : hazScatterPrefixes)
				table.addColumn("![Scatter]("+resourcesDir.getName()+"/"+myPrefix+".png)");
			table.finalizeLine();
			
			lines.addAll(table.build());
			lines.add("");
		}
		
		if (highlightSites != null && !highlightSites.isEmpty()) {
			lines.add("## Highlight Sites");
			lines.add(topLink); lines.add("");
			
			for (Site site : highlightSites) {
				ZScoreResult[] siteZScores = null;
				for (int s=0; s<sites.size(); s++) {
					Site oSite = sites.get(s);
					if (site == oSite || site.equals(oSite)) {
						siteZScores = siteScores.get(s);
						break;
					}
				}
				Preconditions.checkNotNull(siteZScores, "highlight site not found in site list");
				
				System.out.println("Calculating source z-scores for "+site.getName());
				List<Site> mySites = new ArrayList<>();
				mySites.add(site);
				scoreFutures = new ArrayList<>();
				
				MatchesSiteFilter<E> siteFilter = new RuptureComparisonFilter.MatchesSiteFilter<E>(site);
				for (FaultSection sect : subSects) {
					if (sectInsides[sect.getSectionId()]) {
						List<RuptureComparison<E>> siteSectComps = siteFilter.getMatches(sectComps.get(sect.getSectionId()), site);
						scoreFutures.add(exec.submit(new Callable<ZScoreResult[]>() {

							@Override
							public ZScoreResult[] call() throws Exception {
								return ZScoreHistPlot.calcZScores(simProv, siteSectComps, mySites, imts, siteFilter);
							}
						}));
					} else {
						scoreFutures.add(null);
					}
				}
				List<ZScoreResult[]> siteSectScores = new ArrayList<>();
				numSectsProcessed = 0;
				for (Future<ZScoreResult[]> future : scoreFutures) {
					if (future == null) {
						siteSectScores.add(null);
					} else {
						try {
							siteSectScores.add(future.get());
							numSectsProcessed++;
							System.out.println("Done with sect "+numSectsProcessed+"/"+numSectsInside+" for "+site.getName());
						} catch (InterruptedException | ExecutionException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					}
				}
				TableBuilder table = MarkdownUtils.tableBuilder();
				for (int p=0; p<imts.length; p++) {
					IMT imt = imts[p];
					
					String prefix = site.getName().replaceAll("\\W+", "_")+"_"+imt.getPrefix();
					
					List<Color> sectZColors = new ArrayList<>();
					List<Color> sectSdColors = new ArrayList<>();
					for (int s=0; s<subSects.size(); s++) {
						ZScoreResult[] scores = siteSectScores.get(s);
						if (scores == null || Double.isNaN(scores[p].mean)) {
							sectZColors.add(null);
							sectSdColors.add(null);
						} else {
							sectZColors.add(zScoreCPT.getColor((float)scores[p].mean));
							sectSdColors.add(sigmaCPT.getColor((float)scores[p].stdDevFract));
						}
					}
					
//					siteMapMaker.plotScatterScalars(siteLocs, zMeans, zScoreCPT, imt.getDisplayName()+" Site z-scores");
//					siteMapMaker.plot(resourcesDir, prefix+"_site_z_means", " ", 900);
//					siteMapMaker.plotScatterScalars(siteLocs, zStdDevs, sigmaCPT, imt.getDisplayName()+" Site σ-fracts");
//					siteMapMaker.plot(resourcesDir, prefix+"_site_z_sds", " ", 900);
					
					sourceMapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 15f, PlotSymbol.CIRCLE, Color.BLACK);
					
					sourceMapMaker.clearSectColors();
					sourceMapMaker.clearSectScalars();
					sourceMapMaker.plotScatterScalars(
							List.of(site.getLocation()), List.of(siteZScores[p].mean), zScoreCPT,
							site.getName()+" "+imt.getDisplayName()+" source/path z-scores");
					sourceMapMaker.plotSectColors(sectZColors);
					sourceMapMaker.plot(resourcesDir, prefix+"_z_means", " ", 900);
					sourceMapMaker.clearSectColors();
					sourceMapMaker.clearSectScalars();
					sourceMapMaker.plotScatterScalars(
							List.of(site.getLocation()), List.of(siteZScores[p].stdDevFract), sigmaCPT,
							site.getName()+" "+imt.getDisplayName()+" source/path σ-fracts");
					sourceMapMaker.plotSectColors(sectSdColors);
					sourceMapMaker.plot(resourcesDir, prefix+"_z_sds", " ", 900);
					
					table.initNewLine();
					table.addColumn("![Source z-scores]("+resourcesDir.getName()+"/"+prefix+"_z_means.png)");
					table.addColumn("![Source z-scores]("+resourcesDir.getName()+"/"+prefix+"_z_sds.png)");
					table.finalizeLine();
				}
				
				lines.add("### "+site.getName());
				lines.add(topLink); lines.add("");
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static ZScoreResult buildZResult(List<Double> values, boolean zHistRange, boolean rescale) {
		double avgVal = 0d;
		double maxVal = Double.NEGATIVE_INFINITY;
		double minVal = Double.POSITIVE_INFINITY;
		StandardDeviation stdDev = new StandardDeviation();
		
		int numValid = 0;
		for (Double val : values) {
			if (val != null && Double.isFinite(val)) {
				numValid++;
				minVal = Math.min(minVal, val);
				maxVal = Math.max(maxVal, val);
				avgVal += val;
				stdDev.increment(val);
			}
		}
		
		Preconditions.checkState(numValid > 0);
		
		avgVal /= numValid;
		
		HistogramFunction hist;
		if (zHistRange) {
			int numBins;
			if (numValid < 100)
				numBins = 10;
			else if (numValid < 500)
				numBins = 40;
			else
				numBins = 100;
			hist = new HistogramFunction(-ZScoreHistPlot.numStdDev, ZScoreHistPlot.numStdDev, numBins);
		} else {
			double spread = maxVal - minVal;
			double delta;
			if (numValid == 1 || spread > 1d)
				delta = .1;
			else if (spread > 0.5)
				delta = 0.05;
			else
				delta = 0.02;
			hist = HistogramFunction.getEncompassingHistogram(minVal, maxVal, delta);
		}
		
		for (Double val : values)
			if (val != null && Double.isFinite(val))
				hist.add(hist.getClosestXIndex(val), 1d);
		
		if (rescale) {
			double area = ZScoreHistPlot.calcArea(hist);
			hist.scale(1d/area);
		}
		
		return new ZScoreResult(avgVal, stdDev.getResult(), hist);
	}
	
	private static void sigFractHistPlot(File resourcesDir, String prefix, List<Double> values, String name,
			String title) throws IOException {
		ZScoreResult fakeZ = buildZResult(values, false, false);
		double avgVal = fakeZ.mean;
		
		HistogramFunction hist = fakeZ.hist;
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		hist.setName(name);
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		
		Range yRange = new Range(0d, hist.calcSumOfY_Vals()*1.1);
		DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
		meanLine.set(avgVal, 0d);
		meanLine.set(avgVal, yRange.getUpperBound());
		
		meanLine.setName("Mean: "+twoSigFig.format(avgVal));
		
		funcs.add(meanLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE.darker()));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "σ-fract", "Count");
		spec.setLegendVisible(true);
		
		double minX = hist.getMinX();
		double maxX = hist.getMaxX();
		Range xRange = sigFractRange(minX, maxX);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		int width = 800;
		int height = 650;
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(resourcesDir, prefix, gp, width, height, true, true, false);
	}
	
	private static Range sigFractRange(double min, double max) {
		if (min < 0.25)
			min = 0d;
		else if (min < 0.5)
			min = 0.25;
		else
			min = 0.5;
		if (max > 1.5)
			max = 2;
		else if (max > 1.25)
			max = 1.5;
		else if (max > 1)
			max = 1.25;
		else
			max = 1;
			
		return new Range(min, max);
	}
	
	private static void zMeanVsStdDevScatterPlot(File resourcesDir, String prefix, List<Double> zMeans,
			List<Double> zStdDevs, String title) throws IOException {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		for (int i=0; i<zMeans.size(); i++) {
			Double zMean = zMeans.get(i);
			Double zStdDev = zStdDevs.get(i);
			if (zMean == null || !Double.isFinite(zMean))
				continue;
			xy.set(zMean, zStdDev);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_CROSS, 3f, Color.BLACK));
		
		Range xRange = new Range(-2, 2);
		Range yRange = sigFractRange(xy.getMinY(), xy.getMaxY());
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "z-score", "σ-fract");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		int width = 900;
		int height = 800;
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(resourcesDir, prefix, gp, width, height, true, true, false);
	}
	
	private static void zMeanVsStdDevHazardScatterPlot(File resourcesDir, String prefix, List<Double> zMeans,
			List<Double> zStdDevs, List<Double> hazardRatios, HazardRatioType ratioType, CPT cpt, String zLabel,
			String title) throws IOException {
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		
		for (int i=0; i<zMeans.size(); i++) {
			Double zMean = zMeans.get(i);
			if (zMean == null || !Double.isFinite(zMean))
				continue;
			Double zStdDev = zStdDevs.get(i);
			Double hazard = hazardRatios.get(i);
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			xy.set(zMean, zStdDev);
			Color color = cpt.getColor(hazard.floatValue());
			
			minY = Math.min(minY, zStdDev);
			maxY = Math.max(maxY, zStdDev);
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 8f, color));
			
			DefaultXY_DataSet outlineXY = new DefaultXY_DataSet();
			outlineXY.set(xy.get(0));
			
			funcs.add(outlineXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 8f, new Color(0, 0, 0, 127)));
		}
		
		Range xRange = new Range(-2, 2);
		Range yRange = sigFractRange(minY, maxY);
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "z-score", "σ-fract");
		
		PaintScaleLegend cptLabel = RupSetMapMaker.buildCPTLegend(cpt, zLabel);
		Font axisLabelFont = cptLabel.getAxis().getLabelFont();
//		cptLabel.getAxis().setLabelFont(new Font(axisLabelFont.getFontName(),axisLabelFont.getStyle(),20));
		spec.addSubtitle(cptLabel);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		
		int width = 900;
		int height = 900;
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(resourcesDir, prefix, gp, width, height, true, true, false);
	}
	
	public void shutdown() {
		exec.shutdown();
	}

}
