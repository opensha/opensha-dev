package scratch.kevin.simulators.ruptures.rotation;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.MathArrays;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.mod.impl.BaylessSomerville2013DirectivityModifier;
import org.opensha.sha.imr.mod.impl.BaylessSomerville2013DirectivityModifier.DirectivityParams;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceX_Parameter;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.ASK_EventData;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen.ValidationResult;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public abstract class RotatedRupVariabilityPageGen<E> {
	
	private static final boolean D = false;

	private FilterMethod filter;
	
	protected Map<Double, ? extends RotatedRupVariabilityConfig<E>> magConfigs;
	protected Map<Double, SimulationRotDProvider<RotationSpec>> magProvs;
	
	protected double[] calcPeriods;
	
	protected List<Site> sites;
	protected Map<Float, List<Site>> vs30SiteBundles;
	protected List<Float> sourceAzimuths;
	protected List<Float> siteSourceAzimuths;
	protected List<Float> distances;
	protected Map<Double, List<Integer>> magEventIDs;
	private int maxNumEvents = 0;
	private int minNumEvents = Integer.MAX_VALUE;
	private Table<Double, Quantity, List<?>> magQuantitiesTable;
	
	private boolean replotAzimuthDependence = true;
	
	private LoadingCache<E, EqkRupture> gmpeEventCache;
	
	private NGAW2_WrapperFullParam[] gmpes;
	
	protected static final int numExampleRotations = 5;
	
	private Color gmpeColor;
	private List<Color> siteColors;
	
	protected LoadingCache<VarGroupingKey, VariabilityResult> varResultCache;
	protected LoadingCache<VarGroupingKey, VariabilityResult[]> downsampledVarResultCache;
	protected LoadingCache<GMPE_GroupingKey, GMPE_Result> gmpeResultCache;
	protected LoadingCache<EventTermKey, EventTerm> eventTermCache;
	
	private Map<Integer, List<ASK_EventData>> realEventData;
	private int numRealDataSamples;
	
	private boolean hasMagDist;

	public RotatedRupVariabilityPageGen(RotatedRupVariabilityConfig<E> config,
			FilterMethod filter, double mag, SimulationRotDProvider<RotationSpec> prov, double[] calcPeriods) {
		this(filter, emptyMagMap(mag, config), emptyMagMap(mag, prov), calcPeriods);
	}
	
	private static <T> HashMap<Double, T> emptyMagMap(double mag, T value) {
		HashMap<Double, T> map = new HashMap<>();
		map.put(mag, value);
		return map;
	}

	public RotatedRupVariabilityPageGen(FilterMethod filter, Map<Double, ? extends RotatedRupVariabilityConfig<E>> magConfigs,
			Map<Double, SimulationRotDProvider<RotationSpec>> magProvs, double[] calcPeriods) {
		this.magConfigs = magConfigs;
		this.magProvs = magProvs;
		this.filter = filter;
		this.calcPeriods = calcPeriods;
		
		RotatedRupVariabilityConfig<E> config0 = magConfigs.values().iterator().next();
		
		sites = config0.getValues(Site.class, Quantity.SITE);
		
		// look for bundles of 2 or more sites with the same Vs30
		vs30SiteBundles = new HashMap<>();
		for (Site site : sites) {
			float vs30 = site.getParameter(Double.class, Vs30_Param.NAME).getValue().floatValue();
			List<Site> vsSites = vs30SiteBundles.get(vs30);
			if (vsSites == null) {
				vsSites = new ArrayList<>();
				vs30SiteBundles.put(vs30, vsSites);
			}
			vsSites.add(site);
		}
		for (Float vs30 : new ArrayList<Float>(vs30SiteBundles.keySet()))
			if (vs30SiteBundles.get(vs30).size() < 2)
				vs30SiteBundles.remove(vs30);
		
		sourceAzimuths = config0.getValues(Float.class, Quantity.SOURCE_AZIMUTH);
		siteSourceAzimuths = config0.getValues(Float.class, Quantity.SITE_TO_SOURTH_AZIMUTH);
		distances = config0.getValues(Float.class, Quantity.DISTANCE);
		magQuantitiesTable = HashBasedTable.create();
		
		magEventIDs = new HashMap<>();
		for (Double mag : magConfigs.keySet()) {
			RotatedRupVariabilityConfig<E> config = magConfigs.get(mag);
			List<Integer> eventIDs = config.getValues(Integer.class, Quantity.EVENT_ID);
			magEventIDs.put(mag, eventIDs);
			maxNumEvents = Integer.max(maxNumEvents, eventIDs.size());
			minNumEvents = Integer.min(minNumEvents, eventIDs.size());
			Preconditions.checkState(sites.size() == config.getValues(Site.class, Quantity.SITE).size());
			Preconditions.checkState(sourceAzimuths.size() == config.getValues(Float.class, Quantity.SOURCE_AZIMUTH).size());
			Preconditions.checkState(siteSourceAzimuths.size() == config.getValues(Float.class, Quantity.SITE_TO_SOURTH_AZIMUTH).size());
			Preconditions.checkState(distances.size() == config.getValues(Float.class, Quantity.DISTANCE).size());
			for (Quantity quantity : Quantity.values())
				magQuantitiesTable.put(mag, quantity, config.getQuantitiesMap().get(quantity));
		}
		
		CPT siteCPT = null;
		gmpeColor = Color.BLUE;
		if (sites.size() > 1) {
			try {
				siteCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, sites.size());
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			siteColors = new ArrayList<>();
			for (int i=0; i<sites.size(); i++)
				siteColors.add(siteCPT.getColor((float)i+1f).darker());
		}
		
		varResultCache = CacheBuilder.newBuilder().maximumSize(1000).build(new CacheLoader<VarGroupingKey, VariabilityResult>() {

			@Override
			public VariabilityResult load(VarGroupingKey key) throws Exception {
				return calcVarResult(key);
			}
			
		});
		
		downsampledVarResultCache = CacheBuilder.newBuilder().maximumSize(100).build(new CacheLoader<VarGroupingKey, VariabilityResult[]>() {

			@Override
			public VariabilityResult[] load(VarGroupingKey key) throws Exception {
				return calcDownsampledVarResults(key, numRealDataSamples);
			}
			
		});
		
		gmpeResultCache = CacheBuilder.newBuilder().build(new CacheLoader<GMPE_GroupingKey, GMPE_Result>() {

			@Override
			public GMPE_Result load(GMPE_GroupingKey key) throws Exception {
				return calcGMPE(key);
			}
			
		});
		
		gmpeEventCache = CacheBuilder.newBuilder().build(new CacheLoader<E, EqkRupture>() {

			@Override
			public EqkRupture load(E key) throws Exception {
				return buildGMPE_Rupture(key);
			}
			
		});
		
		eventTermCache = CacheBuilder.newBuilder().build(new CacheLoader<EventTermKey, EventTerm>() {

			@Override
			public EventTerm load(EventTermKey key) throws Exception {
				return calcEventTerm(key.eventID, key.mag, key.distance, key.sites);
			}
			
		});
		
		hasMagDist = magConfigs.size() > 3 && distances.size() > 3;
	}
	
	// TODO
	public static Map<Integer, RSQSimEvent> loadEvents(RSQSimCatalog catalog, Collection<Integer> ids) throws IOException {
		System.out.println("Loading "+ids.size()+" events...");
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		List<RSQSimEvent> events = catalog.loader().byIDs(Ints.toArray(ids));
		Preconditions.checkState(events.size() == ids.size(), "Loaded %s events, expected %s", events.size(), ids.size());
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		return eventsMap;
	}
	
	protected void clearCaches() {
		varResultCache.invalidateAll();
		downsampledVarResultCache.invalidateAll();
		gmpeResultCache.invalidateAll();
		gmpeEventCache.invalidateAll();
	}
	
	protected Collection<Integer> getAllEventIDs() {
		HashSet<Integer> idSet = new HashSet<>();
		for (List<Integer> ids : magEventIDs.values())
			idSet.addAll(ids);
		return idSet;
	}
	
	public void setRealEventData(Map<Integer, List<ASK_EventData>> realEventData, int numRealDataSamples) {
		this.realEventData = realEventData;
		this.numRealDataSamples = numRealDataSamples;
	}
	
	protected abstract E getEvent(int eventID);
	
	protected abstract EqkRupture buildGMPE_Rupture(E event);
	
	protected Scenario getBBP_PartB_Scenario(RotatedRupVariabilityConfig<E> config) {
		return null;
	}
	
	/**
	 * Set a GMPE for comparison.
	 * @param gmpe
	 */
	public void setGMPEs(NGAW2_WrapperFullParam... gmpes) {
		this.gmpes = gmpes;
	}
	
	private static final String al_atik = "Al Atik (2010)";
	private static final String aki_richards = "Aki & Richards (1980)";
	
	private enum GMPE_Var {
		TOTAL("Ïƒ", "&sigma;") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return gm.stdDev();
			}
		},
		PHI("ϕ", "&phi;") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return gm.phi();
			}
		},
		PHI_SS("ϕ_ss", "&phi;<sub>SS</sub>") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return Math.sqrt(gm.phi()*gm.phi() - 0.3*0.3);
			}
		},
		TAU("τ", "&tau;") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return gm.tau();
			}
		};

		String symbol;
		String htmlSymbol;
		private GMPE_Var(String symbol, String htmlSymbol) {
			this.symbol = symbol;
			this.htmlSymbol = htmlSymbol;
		}
		
		public abstract double calculate(ScalarGroundMotion gm);
	}
	
	private static Quantity[] qarr(Quantity... quantities) {
		return quantities;
	}
	
	public enum VariabilityType {
		PATH("Path-to-path", "path", "ϕ_p2p", "&phi;<sub>P2P</sub>", al_atik,
				null, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH), // group quantities
				qarr(Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false), // std dev of residuals
		SOUCE_STRIKE("Source-strike", "source_strike", "ϕ_s", "&phi;<sub>s</sub>", aki_richards,
//				GMPE_Var.PHI_SS, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				GMPE_Var.PHI_SS, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
//				Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH),
		WITHIN_EVENT_SS("Within-event, single-site", "within_event_ss", "ϕ_ss", "&phi;<sub>SS</sub>", al_atik,
//				GMPE_Var.PHI_SS, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				GMPE_Var.PHI_SS, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
		WITHIN_EVENT("Within-event", "within_event", "ϕ", "&phi;", al_atik,
				GMPE_Var.PHI, null,
				qarr(Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SITE, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(Quantity.SITE), // singleton quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
//		BETWEEN_EVENTS_SINGLE_PATH("Between-events, single-path", "between_events_single_path", "τ_0", "&tau;<sub>0</sub>", al_atik,
//				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH),
		BETWEEN_EVENTS("Between-events", "between_events", "τ", "&tau;", al_atik,
				GMPE_Var.TAU, null,
				qarr(Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SITE, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				true, "δB_e", "&delta;B<sub>e</sub>"); // medains
		// TODO phi_s2s is relative to the site classification, so leave it out for now
//		SITE_TO_SITE("Site Variability", "site_var", "ϕ_s2s", "&phi;<sub>S2S</sub>", al_atik,
//				false, true, null, null,
//				Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH);
		
		private String name;
		private String prefix;
		private String symbol;
		private String htmlSymbol;
		private String reference;
		private boolean separateSites;
		private boolean separateDists;
		private GMPE_Var gmpeStdDevType;
		private ScatterDisaggQuantity[] disaggScatters;
		private Quantity[] separateQuantities;
		private Quantity[] groupQuantities;
		private Quantity[] variedQuantities;
		private Quantity[] singletonQuantities;
		private boolean stdDevOfMedians; // otherwise std of residuals
		private String variedSymbol;
		private String variedSymbolHTML;
		
		private VariabilityType(String name, String prefix, String symbol, String htmlSymbol, String reference,
				GMPE_Var gmpeStdDevType, ScatterDisaggQuantity[] disaggScatters,
				Quantity[] separateQuantities, // compute sigma independently for all of these
				Quantity[] groupQuantities, // compute sigma of values for each unique combination of these
				Quantity[] variedQuantities, // across all unique combinations of these
				Quantity[] singletonQuantities, // but with only one instance of these per calc
				boolean stdDevOfMedians) { 
			this(name, prefix, symbol, htmlSymbol, reference, gmpeStdDevType, disaggScatters,
					separateQuantities, groupQuantities, variedQuantities, singletonQuantities, stdDevOfMedians, null, null);
		}

		private VariabilityType(String name, String prefix, String symbol, String htmlSymbol, String reference,
				GMPE_Var gmpeStdDevType, ScatterDisaggQuantity[] disaggScatters,
				Quantity[] separateQuantities, // compute sigma independently for all of these
				Quantity[] groupQuantities, // compute sigma of values for each unique combination of these
				Quantity[] variedQuantities, // across all unique combinations of these
				Quantity[] singletonQuantities, // but with only one instance of these per calc
				boolean stdDevOfMedians, String variedSymbol, String variedSymbolHTML) {
			this.name = name;
			this.prefix = prefix;
			this.symbol = symbol;
			this.htmlSymbol = htmlSymbol;
			this.reference = reference;
			this.gmpeStdDevType = gmpeStdDevType;
			this.disaggScatters = disaggScatters;
			this.separateQuantities = separateQuantities;
			this.groupQuantities = groupQuantities;
			this.variedQuantities = variedQuantities;
			this.singletonQuantities = singletonQuantities;
			this.stdDevOfMedians = stdDevOfMedians;
			this.variedSymbol = variedSymbol;
			this.variedSymbolHTML = variedSymbolHTML;
			
			for (Quantity separateQuantity : separateQuantities) {
				if (separateQuantity == Quantity.DISTANCE)
					separateDists = true;
				if (separateQuantity == Quantity.SITE)
					separateSites = true;
			}
			
			Preconditions.checkState(variedQuantities.length > 0, "must have varied quantities");
			Preconditions.checkState(groupQuantities.length > 0 || !stdDevOfMedians,
					"Must have group quantities if computing the std deviation of medians");
		}
		
		public int getRotationsPerStdDev(Map<Quantity, List<?>> quantitiesMap) {
			if (stdDevOfMedians)
				return getCombinations(quantitiesMap, groupQuantities);
			else
				return getCombinations(quantitiesMap, variedQuantities);
		}
		
		public int getCombinations(Map<Quantity, List<?>> quantitiesMap, Quantity[] quantities) {
			int variedCount = 1;
			for (Quantity quantity : quantities) {
				int num = quantitiesMap.get(quantity).size();
				variedCount *= num;
			}
			return variedCount;
		}
		
		public List<String> buildMethodologyLines(RotatedRupVariabilityPageGen<?> pageGen, File resourcesDir) throws IOException {
			List<String> lines = new ArrayList<>();
			String line = name+" variability, denoted "+htmlSymbol;
			if (reference != null)
				line += " in "+reference;
			line += ", is computed separately for each:";
			lines.add(line);
			lines.add("");
			Map<Quantity, List<?>> qMap = pageGen.magQuantitiesTable.rowMap().values().iterator().next();
			for (Quantity quantity : separateQuantities) {
				int num = qMap.get(quantity).size();
				lines.add("* "+quantity+" *["+num+" unique]*");
			}
			lines.add("");
			String symbolAdd = "";
			if (variedSymbolHTML != null)
				symbolAdd = ", "+variedSymbolHTML+",";
			if (stdDevOfMedians) {
				lines.add("We first compute the median natural-log ground motion"+symbolAdd+" for each combination of:");
				lines.add("");
				for (Quantity quantity : groupQuantities) {
					int num = qMap.get(quantity).size();
					lines.add("* "+quantity+" *["+num+" unique]*");
				}
				lines.add("");
				lines.add("That median"+symbolAdd+" is computed across all "+getCombinations(qMap, variedQuantities)+" combinations of:");
				lines.add("");
				for (Quantity quantity : variedQuantities) {
					int num = qMap.get(quantity).size();
					String str = quantity+" *["+num+" unique]*";
					if (singletonQuantities != null && Arrays.binarySearch(singletonQuantities, quantity) >= 0)
						str += " (Singleton)";
					lines.add("* "+str);
				}
				lines.add("");
				if (singletonQuantities != null && singletonQuantities.length > 0) {
					lines.add("We subdivide residual sets such that there is only a single unique value of each quantity "
							+ "marked '(Singleton)' above.");
					lines.add("");
				}
				if (variedSymbolHTML == null)
					line = "We take "+htmlSymbol+" to be the standard deviation of those medians.";
				else
					line = "We take "+htmlSymbol+" to be the standard deviation of all "+variedSymbolHTML+".";
				if (pageGen.sites.size() > 1) {
					if (separateSites) {
						line += " Finally, we compute the mean standard deviation across all sites. This"
								+ " total value is reported as **ALL SITES** and in summary plots/tables.";
					} else {
						line += " This is done separately for each set of sites with the same V<sub>S30</sub> value. Single-site "
								+ "values are also reported.";
					}
				}
				lines.add(line);
				lines.add("");
				if (separateDists && pageGen.distances.size() > 1) {
					lines.add("We also compute distance-independent "+htmlSymbol
							+", which we take to be the mean value across all distances.");
					lines.add("");
				}
			} else {
				lines.add("Then, for each unique combination of:");
				lines.add("");
				for (Quantity quantity : groupQuantities) {
					int num = qMap.get(quantity).size();
					lines.add("* "+quantity+" *["+num+" unique]*");
				}
				lines.add("");
				lines.add("we compute residuals"+symbolAdd+" of the natural-log ground motions "
						+ "(relative to the median), computed across all "+getCombinations(qMap, variedQuantities)+" combinations of:");
				lines.add("");
				for (Quantity quantity : variedQuantities) {
					int num = qMap.get(quantity).size();
					String str = quantity+" *["+num+" unique]*";
					if (singletonQuantities != null && Arrays.binarySearch(singletonQuantities, quantity) >= 0)
						str += " (Singleton)";
					lines.add("* "+str);
				}
				lines.add("");
				if (singletonQuantities != null && singletonQuantities.length > 0) {
					lines.add("We subdivide residual sets such that there is only a single unique value of each quantity "
							+ "marked '(Singleton)' above.");
					lines.add("");
				}
				line = "We take "+htmlSymbol+" to be the standard deviation of all residuals"+symbolAdd+" across each combination of "
						+Joiner.on(", ").join(groupQuantities)+".";
				if (pageGen.sites.size() > 1) {
					if (separateSites)
						line += " We also compute the total standard deviation across all residuals from all sites. This"
								+ " total value is reported as **ALL SITES** and in summary plots/tables.";
					else
						line += " This is done separately for each set of sites with the same V<sub>S30</sub> value.";
				}
				lines.add(line);
				lines.add("");
				if (separateDists && pageGen.distances.size() > 1) {
					lines.add("We also compute distance-independent "+htmlSymbol+", which is computed as the standard deviation of all residuals"
							+symbolAdd+" across all distances. Each residual is still computed relative to the log-median ground motion at "
							+"it's distance.");
					lines.add("");
				}
			}
			
			if (pageGen.numExampleRotations > 0) {
				boolean hasEvent = false;
				for (Quantity quantity : groupQuantities)
					if (quantity == Quantity.EVENT_ID)
						hasEvent = true;
				if (hasEvent) {
					System.out.println("Plotting example!");
					// we can plot an example as we are doing separate calcs for site/distance groups, and are grouping by event
					List<String> repeats = new ArrayList<>();
					for (Quantity quantity : groupQuantities)
						repeats.add(quantity.getName());
					lines.add("Here is an exmample with "+pageGen.numExampleRotations+" rotations, which would be repeated for "
							+ "each combination of ["+Joiner.on(", ").join(repeats)+"]. The site is shown with a blue square, and initially "
							+ "oriented rupture in bold with its hypocenter as a red star and centroid a green circle. Rotations of that "
							+ "rupture are in gray:");
					lines.add("");
					pageGen.plotExample(resourcesDir, "example_"+this.prefix, pageGen.distances.get(0), Lists.newArrayList(variedQuantities));
					lines.add("![Example](resources/example_"+this.prefix+".png)");
					lines.add("");
				} else {
					System.out.println("Not plotting example. SeparateSites? "+separateSites
							+". Quantities: "+Joiner.on("; ").join(groupQuantities));
				}
			}
			
			return lines;
		}
	}
	
	protected abstract double calcVprop(E event);
	protected abstract double getMag(E event);
	protected abstract double getArea(E event);
	protected abstract double getMaxSlip(E event);
	protected abstract double getMeanSlip(E event);
	protected abstract double getSlipStdDev(E event);
	protected abstract double getMeanMidSeisSlip(E event);
	
	private enum ScatterDisaggQuantity {
		V_PROP("Vprop", "v_prop") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return pageGen.calcVprop(pageGen.getEvent(rotation.eventID));
			}
		},
		MAG("Mag", "mag") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return pageGen.getMag(pageGen.getEvent(rotation.eventID));
			}
		},
		AREA("Area", "area") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return pageGen.getArea(pageGen.getEvent(rotation.eventID));
			}
		},
		MAX_SLIP("Max Slip", "max_slip") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return pageGen.getMaxSlip(pageGen.getEvent(rotation.eventID));
			}
		},
		MEAN_SLIP("Mean Slip", "mean_slip") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return pageGen.getMeanSlip(pageGen.getEvent(rotation.eventID));
			}
		},
		SOURCE_AZ("Source Azimuth", "src_az") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return rotation.sourceAz == null ? 0d : rotation.sourceAz;
			}
		},
		SITE_TO_SOURCE_AZ("Site-to-Source Az", "site_source_az") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation) {
				return rotation.siteToSourceAz == null ? 0d : rotation.siteToSourceAz;
			}
		};
		
		private String name;
		private String prefix;
		
		private ScatterDisaggQuantity(String name, String prefix) {
			this.name = name;
			this.prefix = prefix;
		}
		
		public abstract <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, RotationSpec rotation);
	}
	
	protected abstract String getScenarioName();
	protected abstract String getScenarioShortName();
	protected abstract String[] getScenarioMatchCriteria();
	
	protected List<Site> getSites() {
		return sites;
	}
	
	protected List<Integer> getEventIDs(Double mag) {
		return magEventIDs.get(mag);
	}
	
	public void generatePage(File outputDir, double[] periods, List<String> methodSpecificLines) throws IOException {
		generatePage(outputDir, periods, methodSpecificLines, null, null);
	}
	
	protected abstract String getModelName();
	
	public void generatePage(File outputDir, double[] periods, List<String> methodSpecificLines,
			double[] highlightMags, float[] highlightDists) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new LinkedList<>();
		
		String distName = BBP_PartBValidationConfig.DIST_JB ? "Joyner-Boore distance" : "3-dimensional distance";
		String distSymbol = BBP_PartBValidationConfig.DIST_JB ? "Rjb" : "Rrup";
		
		lines.add("# "+getModelName()+" Rotated Rupture Variability, "+getScenarioShortName());
		lines.add("");
		lines.add("This exercise uses translations and rotations to estimate ground motion variability from different "
				+ "sources. We begin by selecting a subset of similar ruptures which match a set of criteria (in this case, "
				+ getScenarioName()+"). Each rupture is then reoriented such that its strike (following the Aki & Richards "
				+ "1980 convention) is 0 degrees (due North, dipping to the right for normal or reverse ruptures). For each site, "
				+ "ruptures are translated such that their scalar seismic moment centroid is directly North of the site, and their "
				+ distName+" ("+distSymbol+") is as specified (we consider "+distances.size()+" distance[s] here).");
		lines.add("");
		lines.add("We then  perform various rotations. We rotate the rupture in place around its centroid, holding the site-to-source "
				+ "centroid path and "+distSymbol+" constant (henceforth '"+Quantity.SOURCE_AZIMUTH.getName()+"'). We also rotate ruptures around the site, "
				+ "holding "+distSymbol+" and source orientation relative to the site constant but sampling different various paths (henceforth "
				+ "'"+Quantity.SITE_TO_SOURTH_AZIMUTH.getName()+"'). We do this for each unique combination of "+Quantity.SOURCE_AZIMUTH+", "
				+ Quantity.SITE_TO_SOURTH_AZIMUTH.getName()+", "+Quantity.DISTANCE.getName()+", "+Quantity.SITE.getName()+", and "
				+ Quantity.EVENT_ID.getName()+".");
		lines.add("");
		
		if (methodSpecificLines != null && !methodSpecificLines.isEmpty()) {
			lines.addAll(methodSpecificLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Rupture Rotation Parameters");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Quantity", "Variations", "Description");
		if (maxNumEvents == minNumEvents)
			table.addLine(Quantity.EVENT_ID.getName(), maxNumEvents, Quantity.EVENT_ID.getDescription());
		else
			table.addLine(Quantity.EVENT_ID.getName(), minNumEvents+" - "+maxNumEvents, Quantity.EVENT_ID.getDescription());
		table.addLine(Quantity.SITE.getName(), sites.size(), Quantity.SITE.getDescription());
		table.addLine(Quantity.SOURCE_AZIMUTH.getName(), sourceAzimuths.size(), Quantity.SOURCE_AZIMUTH.getDescription());
		table.addLine(Quantity.SITE_TO_SOURTH_AZIMUTH.getName(), siteSourceAzimuths.size(), Quantity.SITE_TO_SOURTH_AZIMUTH.getDescription());
		table.addLine(Quantity.DISTANCE.getName(), Joiner.on(", ").join(distances)+" km", Quantity.DISTANCE.getDescription());
		table.addLine("**Total # Simulations**", "**"+magConfigs.values().iterator().next().getRotations().size()+"**",
				"Total number of combinations of the above.");
		lines.addAll(table.build());
		lines.add("");
		lines.add("## "+getScenarioShortName()+" Rupture Match Criteria");
		lines.add(topLink); lines.add("");
		String[] criteria = getScenarioMatchCriteria();
		lines.add("We condisder "+maxNumEvents+" events which match the following criteria:");
		lines.add("");
		for (String criterion : criteria)
			lines.add("* "+criterion);
		lines.add("");

		if (filter != null) {
			lines.add("In the case of more than "+maxNumEvents+" available matches, we use the "+filter.getName()+" filter method: "+filter.getDescription());
			lines.add("");
		}
		
		if (magConfigs.size() == 1 && this instanceof RSQSimRotatedRupVariabilityPageGen) {
			lines.add("### Fault Section Counts");
			lines.add(topLink); lines.add("");
			lines.add("This tables gives a list of all fault sections which participate in the ruptures matching the above "
					+ "criteria. Many ruptures include multiple sections, so the sum of counts may be larger than the number "
					+ "of events.");
			lines.add("");
			Map<String, Integer> parentCountsMap = new HashMap<>();
			List<Double> mags = new ArrayList<>();
			int totalCount = 0;
			for (int eventID : magConfigs.values().iterator().next().getValues(Integer.class, Quantity.EVENT_ID)) {
				RSQSimEvent event = (RSQSimEvent)getEvent(eventID);
				mags.add(event.getMagnitude());
				RSQSimSubSectEqkRupture rup;
				try {
					rup = (RSQSimSubSectEqkRupture)gmpeEventCache.get((E)event);
				} catch (ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				HashSet<String> parentNames = new HashSet<>();
				for (FaultSectionPrefData sect : rup.getSubSections())
					parentNames.add(sect.getParentSectionName());
				
				for (String parentName : parentNames) {
					Integer prevCount = parentCountsMap.get(parentName);
					if (prevCount == null)
						prevCount = 0;
					parentCountsMap.put(parentName, prevCount+1);
				}
				
				totalCount += parentNames.size();
			}
			List<String> parentNames = ComparablePairing.getSortedData(parentCountsMap);
			Collections.reverse(parentNames);
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Section Name", "Participation Count");
			for (String parentName : parentNames)
				table.addLine(parentName, parentCountsMap.get(parentName));
			table.addLine("SUM OF PARENT PARTICIPATIONS", totalCount);
			table.addLine("# UNIQUE PARENTS", parentNames);
			lines.addAll(table.build());
			lines.add("");
			double[] magsArray = Doubles.toArray(mags);
			lines.add("Actual magnitude range: ["+(float)StatUtils.min(magsArray)+","
					+(float)StatUtils.max(magsArray)
					+"], average: "+(float)StatUtils.mean(magsArray)
					+", stdDev: "+(float)Math.sqrt(StatUtils.variance(magsArray)));
			lines.add("");
		}
		
		lines.add("## Sites");
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine("Name", "Location", "Vs30 (m/s)", "Z1.0 (km)", "Z2.5 (km)");
		for (Site s : this.sites) {
			table.initNewLine();
			table.addColumn(s.getName());
			Location loc = s.getLocation();
			table.addColumn("*"+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
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
		lines.add("");
		
		lines.add("## Result Summary Table");
		lines.add("");
		int summaryTableIndex = lines.size();
		lines.add("");
		
		List<Double> plotMags = new ArrayList<>();
		if (highlightMags != null) {
			for (double mag : highlightMags)
				if (magConfigs.containsKey(mag))
					plotMags.add(mag);
		} else {
			plotMags.addAll(magConfigs.keySet());
		}
		Collections.sort(plotMags);
		List<Float> plotDists;
		if (highlightDists == null) {
			plotDists = distances;
		} else {
			for (float dist : highlightDists)
				Preconditions.checkState(distances.contains(dist), "Highlight distance not in dataset: %s", dist);
			plotDists = Floats.asList(highlightDists);
		}
		
		Map<VariabilityType, File[]> magDistPlots = new HashMap<>();
		
		HashSet<VariabilityType> computedTypes = new HashSet<>();
		
		for (VariabilityType type : VariabilityType.values()) {
			if (type.getRotationsPerStdDev(magQuantitiesTable.rowMap().values().iterator().next()) == 1)
				continue;
			if (!type.separateSites && vs30SiteBundles.isEmpty() && !type.stdDevOfMedians)
				continue;
			if (type.singletonQuantities != null && Arrays.binarySearch(type.singletonQuantities, Quantity.SITE) >= 0
					&& vs30SiteBundles.isEmpty())
				continue;
			computedTypes.add(type);
			
			System.out.println("*** Processing variability type: "+type.name+" ***");
			lines.add("## "+type.name+" Variability");
			lines.add(topLink); lines.add("");
			
			lines.add("### "+type.name+" Variability Methodology");
			lines.add(topLink); lines.add("");
			
			lines.addAll(type.buildMethodologyLines(this, resourcesDir));
			lines.add("");
			
			if (hasMagDist) {
				System.out.println("*** Plotting Mag-Dist for "+type.name+" ***");
				lines.add("### "+type.name+" Variability Mag-Distance Plots");
				lines.add(topLink); lines.add("");
				if (sites.size() > 1 && type.separateSites) {
					lines.add("#### "+type.name+" Site-Specific Variability Mag-Distance Plots");
					lines.add(topLink); lines.add("");
					table = MarkdownUtils.tableBuilder();
					table.initNewLine();
					table.addColumn("Site");
					for (double period : periods)
						table.addColumn(period == (int)period ? (int)period+"s" : (float)period+"s");
					table.finalizeLine();
					for (Site site : sites) {
						table.initNewLine();
						table.addColumn("**"+site.getName()+"**");
						String prefix = type.prefix+"_"+site.getName()+"mag_dist_std_dev";
						File[] files = plotMagDistCheckerboard(resourcesDir, prefix, type, site, periods).get(MagDistPlotType.SIM_STD_DEV);
						for (File file : files)
							table.addColumn("![Mag-Dist Plot](resources/"+file.getName()+")");
						table.finalizeLine();
					}
					lines.addAll(table.build());
					lines.add("");
					lines.add("#### "+type.name+" Total Variability Mag-Distance Plots");
					lines.add(topLink); lines.add("");
				}
				String prefix = type.prefix+"_mag_dist_std_dev";
				Map<MagDistPlotType, File[]> filesMap = plotMagDistCheckerboard(resourcesDir, prefix, type, null, periods);
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				if (filesMap.keySet().size() > 1)
					table.addColumn("Plot Type");
				for (double period : periods)
					table.addColumn(period == (int)period ? (int)period+"s" : (float)period+"s");
				table.finalizeLine();
				for (MagDistPlotType plotType : MagDistPlotType.values()) {
					if (!filesMap.containsKey(plotType))
						continue;
					table.initNewLine();
					if (filesMap.keySet().size() > 1)
						table.addColumn("**"+plotType.getName()+"**");
					for (File file : filesMap.get(plotType))
						table.addColumn("![Mag-Dist Plot](resources/"+file.getName()+")");
					table.finalizeLine();
				}
				magDistPlots.put(type, filesMap.get(MagDistPlotType.SIM_STD_DEV));
				lines.addAll(table.build());
				lines.add("");
				System.out.println("*** DONE Mag-Dist for "+type.name+" ***");
			}
			
			for (Double mag : plotMags) {
				System.out.println("*** Plotting M="+mag+" for "+type.name+" ***");
				RotatedRupVariabilityConfig<E> config = magConfigs.get(mag);
				String curHeading = "##";
				String uniqueTitleSection = mag == null ? type.name : "M"+optionalDigitDF.format(mag)+" "+type.name;
				if (!type.separateDists || distances.size() == 1) {
					curHeading += "#";
					lines.add(curHeading+" "+uniqueTitleSection+" Results");
					lines.add(topLink); lines.add("");
				}
				
				if (type.separateDists) {
					for (float distance : plotDists) {
						lines.addAll(generateVariabilityLines(config, type, mag, distance, periods, topLink, resourcesDir));
						lines.add("");
					}
					// also do total over all dists
					lines.addAll(generateVariabilityLines(config, type, mag, null, periods, topLink, resourcesDir));
					lines.add("");
				} else {
					lines.addAll(generateVariabilityLines(config, type, mag, null, periods, topLink, resourcesDir));
					lines.add("");
				}
			}
			varResultCache.invalidateAll();
			downsampledVarResultCache.invalidateAll();
			
			System.out.println("*** DONE with variability type: "+type.name+" ***");
		}
		
		// build summary table(s)
		List<String> summaryLines = new ArrayList<>();
		if (hasMagDist) {
			System.out.println("*** Building Mag-Dist summary table ***");
			summaryLines.add("### Mag-Dist Plots");
			summaryLines.add(topLink); summaryLines.add("");
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("Type");
			table.addColumn("Notation");
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+"s Std. Dev.");
			table.finalizeLine();
			for (VariabilityType type : VariabilityType.values()) {
				if (!computedTypes.contains(type))
					continue;
				table.initNewLine();
				table.addColumn("["+type.name+"](#"+MarkdownUtils.getAnchorName(type.name+" Variability")+")");
				table.addColumn(type.htmlSymbol);
				File[] plots = magDistPlots.get(type);
				for (File plot : plots)
					table.addColumn("![Mag-Dist Plot](resources/"+plot.getName()+")");
				table.finalizeLine();
			}
			summaryLines.addAll(table.build());
			System.out.println("*** DONE Mag-Dist summary table ***");
		}
		System.out.println("*** Building summary table ***");
		for (Double mag : plotMags) {
			if (mag != null && plotMags.size() > 1) {
				if (!summaryLines.isEmpty())
					summaryLines.add("");
				summaryLines.add("### M"+optionalDigitDF.format(mag)+" Result Summary Table");
				summaryLines.add(topLink); summaryLines.add("");
			}
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("Type");
			table.addColumn("Notation");
			table.addColumn("Distance");
			table.addColumn("T-independent Std. Dev.");
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+"s Std. Dev.");
			table.finalizeLine();
			TableBuilder gmpeTable = null;
			if (gmpes != null && gmpes.length > 0) {
				gmpeTable = MarkdownUtils.tableBuilder();
				gmpeTable.initNewLine();
				gmpeTable.addColumn("Type");
				gmpeTable.addColumn("Notation");
				gmpeTable.addColumn("Distance");
				for (NGAW2_WrapperFullParam gmpe : gmpes)
					for (double period : periods)
						gmpeTable.addColumn(gmpe.getShortName()+" "
								+optionalDigitDF.format(period)+"s");
				gmpeTable.finalizeLine();
			}
			for (VariabilityType type : VariabilityType.values()) {
				if (!computedTypes.contains(type))
					continue;
				List<Float> myDists;
				if (type.separateDists) {
					myDists = new ArrayList<>(this.distances);
					if (distances.size() > 1)
						myDists.add(null);
				} else {
					myDists = new ArrayList<>();
					myDists.add(null);
				}
				for (Float distance : myDists) {
					table.initNewLine();
					VarGroupingKey key = new VarGroupingKey(type, mag, distance, sites);
					VariabilityResult varResult;
					try {
						varResult = varResultCache.get(key);
					} catch (ExecutionException e) {
						if (e.getCause() instanceof IOException)
							throw (IOException)e.getCause();
						if (e.getCause() instanceof RuntimeException)
							throw (RuntimeException)e.getCause();
						throw ExceptionUtils.asRuntimeException(e);
					}
					table.addColumn(type.name);
					table.addColumn(type.htmlSymbol);
					if (distance == null && distances.size() == 1)
						table.addColumn(optionalDigitDF.format(distances.get(0))+" km");
					else if (distance == null && distances.size() > 1)
						table.addColumn("(all)");
					else
						table.addColumn(optionalDigitDF.format(distance)+" km");
					if (type.stdDevOfMedians) {
						table.addColumn(optionalDigitDF.format(varResult.getPeriodIndepMedianStdDevSet().stdDev));
						for (double period : periods)
							table.addColumn(optionalDigitDF.format(varResult.getMedianStdDevSet(period).stdDev));
					} else {
						table.addColumn(optionalDigitDF.format(varResult.getPeriodIndepResidualStdDevSet().total));
						for (double period : periods)
							table.addColumn(optionalDigitDF.format(varResult.getResidualStdDevSet(period).total));
					}
					table.finalizeLine();
					if (gmpeTable != null && distance != null && type.gmpeStdDevType != null) {
						gmpeTable.initNewLine();
						gmpeTable.addColumn(type.name);
						gmpeTable.addColumn(type.gmpeStdDevType.htmlSymbol);
						gmpeTable.addColumn(optionalDigitDF.format(distance)+" km");
						for (NGAW2_WrapperFullParam gmpe : gmpes) {
							GMPE_GroupingKey gKey = new GMPE_GroupingKey(gmpe, sites.get(0), mag, distance);
							GMPE_Result gVals;
							try {
								gVals = gmpeResultCache.get(gKey);
							} catch (ExecutionException e) {
								throw ExceptionUtils.asRuntimeException(e);
							}
							for (int p=0; p<periods.length; p++) {
								ScalarGroundMotion[] gms = gVals.gms[p];
								double[] vals = new double[gms.length];
								for (int i=0; i<vals.length; i++)
									vals[i] = type.gmpeStdDevType.calculate(gms[i]);
								gmpeTable.addColumn(optionalDigitDF.format(StatUtils.mean(vals)));
							}
						}
						gmpeTable.finalizeLine();
					}
				}
			}
			summaryLines.addAll(table.build());
			if (gmpeTable != null) {
				summaryLines.add("");
				if (mag != null && plotMags.size() > 1) {
					summaryLines.add("#### M"+optionalDigitDF.format(mag)+" GMPE Table");
					summaryLines.add(topLink); summaryLines.add("");
				} else {
					summaryLines.add("### GMPE Table");
					summaryLines.add(topLink); summaryLines.add("");
				}
				summaryLines.addAll(gmpeTable.build());
			}
			if (plotDists.size() > 1) {
				summaryLines.add("");
				if (mag != null && plotMags.size() > 1) {
					summaryLines.add("#### M"+optionalDigitDF.format(mag)+" Dist-Dependent Plot Table");
					summaryLines.add(topLink); summaryLines.add("");
				} else {
					summaryLines.add("### Dist-Dependent Plot Table");
					summaryLines.add(topLink); summaryLines.add("");
				}
				table = MarkdownUtils.tableBuilder();
				for (VariabilityType type : VariabilityType.values()) {
					if (!computedTypes.contains(type) || !type.separateDists)
						continue;
					String prefix = type.prefix+(mag == null ? "" : "_m"+optionalDigitDF.format(mag))+"_dist_periods";
					File plot = plotMultiDistPeriodDependentStdDevs(resourcesDir, prefix, type, mag, plotDists);
					table.addLine("**"+type.htmlSymbol+"**", "!["+type.htmlSymbol+"](resources/"+plot.getName()+")");
				}
				summaryLines.addAll(table.build());
				summaryLines.add("");
			}
		}
		System.out.println("*** DONE summary table ***");
		lines.addAll(summaryTableIndex, summaryLines);
		
		if (!hasMagDist) {
			System.out.println("*** doing event term scatters ***");
			Site[] scatterSites;
			lines.add("## Event Term Scatters");
			lines.add(topLink); lines.add("");
			if (sites.size() == 1 || vs30SiteBundles.isEmpty()) {
				scatterSites = new Site[] { sites.get(0) };
			} else {
				scatterSites = null;
				Float vs30 = null;
				for (Float myVs30 : vs30SiteBundles.keySet()) {
					List<Site> bundle = vs30SiteBundles.get(myVs30);
					if (scatterSites == null || bundle.size() > scatterSites.length) {
						scatterSites = bundle.toArray(new Site[0]);
						vs30 = myVs30;
					}
				}
				
				if (scatterSites.length < sites.size()) {
					lines.add("These results use the "+scatterSites.length+" sites with Vs30="+vs30);
					lines.add("");
				}
			}
			Double mag = plotMags.get(0);
			for (EventDisaggQuantity q : EventDisaggQuantity.values()) {
				lines.add("### "+q.name+" Event Term Scatters");
				lines.add(topLink); lines.add("");
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				table.addColumn("");
				for (double period : periods)
					table.addColumn(optionalDigitDF.format(period)+" s");
				table.finalizeLine();
				
				for (Float dist : distances) {
					table.initNewLine();
					table.addColumn("**"+optionalDigitDF.format(dist)+" km**");
					for (double period : periods) {
						File plot = plotEventTermScatter(q, mag, dist, period, scatterSites, resourcesDir);
						if (plot == null)
							table.addColumn("*N/A*");
						else
							table.addColumn("![plot]("+resourcesDir.getName()+"/"+plot.getName()+")");
					}
					table.finalizeLine();
				}
				lines.addAll(table.build());
			}
			System.out.println("*** DONE event term scatters ***");
		}
		
		if (!hasMagDist) {
			System.out.println("*** doing directivity ***");
			Site[] scatterSites;
			lines.add("## Directivity Comparisons");
			lines.add(topLink); lines.add("");
			if (sites.size() == 1 || vs30SiteBundles.isEmpty()) {
				scatterSites = new Site[] { sites.get(0) };
			} else {
				scatterSites = null;
				Float vs30 = null;
				for (Float myVs30 : vs30SiteBundles.keySet()) {
					List<Site> bundle = vs30SiteBundles.get(myVs30);
					if (scatterSites == null || bundle.size() > scatterSites.length) {
						scatterSites = bundle.toArray(new Site[0]);
						vs30 = myVs30;
					}
				}
				
				if (scatterSites.length < sites.size()) {
					lines.add("These results use the "+scatterSites.length+" sites with Vs30="+vs30);
					lines.add("");
				}
			}
			Double mag = plotMags.get(0);
			
			System.out.println("Calculating directivity...");
			Table<Float, Double, XY_DataSet> scatters =
					calcDirectivityComparisons(mag, periods, scatterSites, resourcesDir);
			System.out.println("DONE calculating directivity");
			if (DIRECTIVITY_DEBUG) {
				lines.add("Directivity comparisons for individual ruptures can be found [here]("
						+resourcesDir.getName()+"/"+DIRECTIVITY_DEBUG_DIRNAME+"/README.md).");
				lines.add("");
			}
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			table.addColumn("");
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+" s");
			table.finalizeLine();
			
			for (Float dist : distances) {
				table.initNewLine();
				table.addColumn("**"+optionalDigitDF.format(dist)+" km**");
				for (double period : periods) {
					String prefix = "directivity_"+optionalDigitDF.format(dist)
						+"km_"+optionalDigitDF.format(period)+"s";
					if (plotDirectivityComparison(resourcesDir, prefix, scatters.get(dist, period)))
						table.addColumn("![plot]("+resourcesDir.getName()+"/"+prefix+".png)");
					else
						table.addColumn("*N/A*");
				}
				table.finalizeLine();
			}
			lines.addAll(table.build());
			System.out.println("*** DONE directivity scatters ***");
		}
		
		if (sourceAzimuths.size() > 1 || siteSourceAzimuths.size() > 1 && !plotMags.isEmpty()) {
			lines.add("## Azumth Dependence");
			lines.add(topLink); lines.add("");
			
			List<Site> mySites = new ArrayList<>(sites);
			if (sites.size() > 1)
				mySites.add(null);
			
			for (Site site : mySites) {
				String sitePrefix = sites.size() > 1 && site != null ? site.getName()+"_" : "";
				if (sites.size() > 1) {
					if (site == null)
						lines.add("### All Sites Azumth Dependence");
					else
						lines.add("### "+site.getName()+" Azumth Dependence");
					lines.add(topLink); lines.add("");
				}
				for (Quantity quantity : new Quantity[] { Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH }) {
					if (magQuantitiesTable.get(plotMags.get(0), quantity).size() < 2)
						continue;
					
					Map<Quantity, List<?>> myQuantitiesMap = new HashMap<>();
					for (Quantity q2 : Quantity.values()) {
						if (q2 == quantity) {
							List<Object> list = new ArrayList<>();
							list.add(null); // single item;
							myQuantitiesMap.put(q2, list);
						} else {
							myQuantitiesMap.put(q2, magQuantitiesTable.get(plotMags.get(0), q2));
						}
					}
					List<VariabilityType> myTypes = new ArrayList<>();
					for (VariabilityType type : computedTypes) {
						if (type.getRotationsPerStdDev(myQuantitiesMap) > 1)
							myTypes.add(type);
						else
							System.out.println("Skipping type "+type.name+" for quantity "+quantity.getName());
					}
					myTypes.add(null);
					
					System.out.println("*** Doing "+quantity.getName()+" vs Dist ***");
					
					if (sites.size() > 1) {
						if (site == null)
							lines.add("#### All Sites "+quantity.getName()+" Dependence");
						else
							lines.add("#### "+site.getName()+" "+quantity.getName()+" Dependence");
						lines.add(topLink); lines.add("");
					} else {
						lines.add("### "+quantity.getName()+" Dependence");
						lines.add(topLink); lines.add("");
					}
					
					table = MarkdownUtils.tableBuilder();
					table.initNewLine();
					if (plotMags.size() > 1)
						table.addColumn("Mag");
					table.addColumn("Type");
					
					for (double period : periods)
						table.addColumn(optionalDigitDF.format(period)+"s");
					table.finalizeLine();
					
					Map<Double, Map<VariabilityType, File[]>> plotMaps2D = new HashMap<>();
					Map<Double, Map<VariabilityType, File[]>> plotMaps = new HashMap<>();
					
					for (double magnitude : plotMags) {
						String prefix = sitePrefix+"m"+optionalDigitDF.format(magnitude)+"_dist_"+quantity.name();
						if (hasMagDist) {
							Map<VariabilityType, File[]> plotsMap = plotQuantityDistCheckerboards(
									resourcesDir, prefix+"_2d", quantity, computedTypes, site, magnitude, periods, replotAzimuthDependence);
							plotMaps2D.put(magnitude, plotsMap);
						}
						
						Map<VariabilityType, File[]> plotsMap = plotQuantityDistFuncs(
									resourcesDir, prefix, quantity, computedTypes, site, magnitude, periods, plotDists, replotAzimuthDependence);
						plotMaps.put(magnitude, plotsMap);
					}
					
					boolean[] plot2Ds;
					if (hasMagDist)
						plot2Ds = new boolean[] { true, false };
					else
						plot2Ds = new boolean[] { false };
					
					for (VariabilityType type : myTypes) {
						for (double magnitude : plotMags) {
							for (boolean plot2D : plot2Ds) {
								table.initNewLine();
								if (plotMags.size() > 1)
									table.addColumn("**M"+optionalDigitDF.format(magnitude)+"**");
								if (type == null)
									table.addColumn("**Median SA**");
								else
									table.addColumn("**"+type.htmlSymbol+"**");
								File[] files = plot2D ? plotMaps2D.get(magnitude).get(type) : plotMaps.get(magnitude).get(type);
								for (File file : files)
									table.addColumn("!["+quantity.getName()+"](resources/"+file.getName()+")");
								table.finalizeLine();
							}
						}
					}
					
					lines.addAll(table.build());
					lines.add("");
					
					System.out.println("*** DONE Doing "+quantity.getName()+" vs Dist ***");
				}
			}
		}
		
		Map<Scenario, Double> partBConfigMap = new HashMap<>();
		for (Double mag : magConfigs.keySet()) {
			Scenario scenario = getBBP_PartB_Scenario(magConfigs.get(mag));
			if (scenario != null)
				partBConfigMap.put(scenario, mag);
		}
		
		if (!partBConfigMap.isEmpty()) {
			lines.add("## BBP PartB Comparison");
			lines.add(topLink); lines.add("");
			
			lines.add("Here we attempt to reproduce the SCEC BroadBand Platform \"Part B\" validation exercise as defined in:");
			lines.add("");
			lines.add("*Goulet, C. A., Abrahamson, N. A., Somerville, P. G., & Wooddell, K. E. (2014). The SCEC broadband platform "
					+ "validation exercise: Methodology for code validation in the context of seismicâ€�hazard analyses. "
					+ "Seismological Research Letters, 86(1), 17-26.* [(link)](https://pubs.geoscienceworld.org/ssa/srl/article/86/1/17/315438/"
					+ "the-scec-broadband-platform-validation-exercise)");
			lines.add("");
			lines.add("The BBP exercise positioned sites in a 'racetrack' around the ruptures. Here, we instead position and rotate the ruptures"
					+ " around the site in order to work in 3-D with CyberShake reciprical calculations. Results for official scenarios and distances"
					+ " are in **bold**, results for additional magnitudes or distances not defined in the Goulet et. al. (2014) are *italicised*.");
			lines.add("");
			lines.add("### BBP PartB Summary Table");
			lines.add(topLink); lines.add("");
			int summaryIndex = lines.size();
			lines.add("");
			
			Table<Scenario, Double, List<ValidationResult>> validationResultsTable = HashBasedTable.create();
			Table<Scenario, Double, String> plotsTable = HashBasedTable.create();
			
			TableBuilder summaryTable = MarkdownUtils.tableBuilder();
			
			summaryTable.initNewLine();
			summaryTable.addColumn("Scenario");
			if (sites.size() > 1)
				summaryTable.addColumn("Site");
			for (Float distance : distances)
				summaryTable.addColumn(distance+" km");
			summaryTable.finalizeLine();
			
			for (Scenario scenario : Scenario.values()) {
				Double mag = partBConfigMap.get(scenario);
				if (mag == null)
					continue;
				RotatedRupVariabilityConfig<E> config = magConfigs.get(mag);
				
				lines.add("### BBP PartB, "+scenario.getName());
				lines.add(topLink); lines.add("");
				
				int numEvents = magEventIDs.get(mag).size();
				
				SimulationRotDProvider<RotationSpec> simProv = magProvs.get(mag);
				
				TableBuilder plotTable = MarkdownUtils.tableBuilder();
				
				plotTable.initNewLine();
				if (sites.size() > 1)
					plotTable.addColumn("Site");
				for (Float distance : distances)
					plotTable.addColumn(distance+" km");
				plotTable.finalizeLine();
				
				boolean footwall = scenario.isFootwallOnly();
				
				Table<Float, Double, List<DiscretizedFunc>> distVsBunleTable = HashBasedTable.create();
				Map<Double, List<Site>> vsSiteMap = new HashMap<>();
				
				for (Site site : sites) {
					summaryTable.initNewLine();
					plotTable.initNewLine();
					summaryTable.addColumn("**"+scenario.getShortName()+"**");
					if (sites.size() > 1) {
						summaryTable.addColumn("**"+site.getName()+"**");
						plotTable.addColumn("**"+site.getName()+"**");
					}
					
					double vs30 = site.getParameter(Double.class, Vs30_Param.NAME).getValue();
					List<Site> otherSites = vsSiteMap.get(vs30);
					if (otherSites == null) {
						otherSites = new ArrayList<>();
						vsSiteMap.put(vs30, otherSites);
					}
					otherSites.add(site);
					
					for (Float distance : distances) {
						boolean official = scenario.isOfficialCriteria() && isOfficialPartB_Distance(distance);
						
						List<DiscretizedFunc> rd50s = new ArrayList<>();
						for (RotationSpec rotation : config.getRotationsForQuantities(Quantity.SITE, site, Quantity.DISTANCE, distance)) {
							Float sourceAz = (Float)rotation.getValue(Quantity.SOURCE_AZIMUTH);
							if (footwall && sourceAz < 180f)
								// only footwall
								continue;
							DiscretizedFunc spectrum = simProv.getRotD50(site, rotation, 0);
							rd50s.add(spectrum);
						}
						
						List<DiscretizedFunc> prevRDs = distVsBunleTable.get(distance, vs30);
						if (prevRDs == null) {
							prevRDs = new ArrayList<>();
							distVsBunleTable.put(distance, vs30, prevRDs);
						}
						prevRDs.addAll(rd50s);
						
						String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
						
						// plot it
						String prefix = "bbp_partB_"+scenario.getPrefix()+distPrefix+(sites.size() > 1 ? "_"+site.getName() : "");
//						
						List<ValidationResult> results = BBP_PartBValidationPageGen.calcPlotScenarioResults(
								scenario, distance, vs30, numEvents, rd50s, resourcesDir, prefix);
						
						validationResultsTable.put(scenario, distance.doubleValue(), results);
						plotsTable.put(scenario, distance.doubleValue(), prefix+".png");
						
						boolean passes = true;
						for (ValidationResult result : results)
							passes = passes && result.passes();
						
						summaryTable.addColumn(bbpTableStr(passes, official));
						plotTable.addColumn("![PartB Plot](resources/"+prefix+".png)");
					}
					
					summaryTable.finalizeLine();
					plotTable.finalizeLine();
				}
				
				// now do any combined ones
				for (Double vs30 : vsSiteMap.keySet()) {
					List<Site> sites = vsSiteMap.get(vs30);
					if (sites.size() == 1)
						continue;
					
					summaryTable.initNewLine();
					plotTable.initNewLine();
					summaryTable.addColumn("**"+scenario.getShortName()+"**");
					
					String siteName = sites.size()+"sites_Vs30_"+vs30.intValue();
					
					if (sites.size() > 1) {
						summaryTable.addColumn("**"+siteName+"**");
						plotTable.addColumn("**"+siteName+"**");
					}
					
					for (Float distance : distances) {
						boolean official = scenario.isOfficialCriteria() && isOfficialPartB_Distance(distance);
						List<DiscretizedFunc> rd50s = distVsBunleTable.get(distance, vs30);
						
						String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
						
						// plot it
						String prefix = "bbp_partB_"+scenario.getPrefix()+distPrefix+"_"+siteName;
//						
						List<ValidationResult> results = BBP_PartBValidationPageGen.calcPlotScenarioResults(
								scenario, distance, vs30, numEvents, rd50s, resourcesDir, prefix);
						
						validationResultsTable.put(scenario, distance.doubleValue(), results);
						plotsTable.put(scenario, distance.doubleValue(), prefix+".png");
						
						boolean passes = true;
						for (ValidationResult result : results)
							passes = passes && result.passes();
						
						summaryTable.addColumn(bbpTableStr(passes, official));
						plotTable.addColumn("![PartB Plot](resources/"+prefix+".png)");
					}
					
					summaryTable.finalizeLine();
					plotTable.finalizeLine();
				}
				
				lines.addAll(plotTable.wrap(5, 0).build());
				lines.add("");
			}
			
			lines.addAll(summaryIndex, summaryTable.build());
			
			BBP_PartBValidationPageGen.writeResultsCSV(new File(resourcesDir, "part_b_results.csv"), validationResultsTable, plotsTable);
		}
		
		System.out.println("*** Writing CSVs ***");
		lines.add("## CSV Files");
		lines.add(topLink); lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine("Magnitude", "Distance", "Site", "CSV File");
		for (double mag : plotMags) {
			for (Float dist : plotDists) {
				for (Site site : sites) {
					File file = writeSiteCSV(resourcesDir, site, mag, dist);
					table.addLine("M"+optionalDigitDF.format(mag), dist+" km", site.getName(),
							"["+file.getName()+"](resources/"+file.getName()+")");
				}
			}
		}
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private String bbpTableStr(boolean passes, boolean official) {
		if (official) {
			if (passes)
				return "**PASS**";
			return "**FAIL**";
		}
		if (passes)
			return "*(PASS)*";
		return "*(FAIL)*";
	}
	
	private boolean isOfficialPartB_Distance(float distance) {
		for (double d : BBP_PartBValidationConfig.OFFICIAL_DISTANCES)
			if ((float)d == distance)
				return true;
		return false;
	}
	
	private List<String> generateVariabilityLines(RotatedRupVariabilityConfig<E> config, VariabilityType type, Double mag, Float distance,
			double[] periods, String topLink, File resourcesDir) throws IOException {
		List<String> lines = new ArrayList<>();
//		if (!type.stdDevOfMedians)
//			return lines;
		String curHeading = "##";
		String uniqueTitleSection = mag == null ? type.name : "M"+optionalDigitDF.format(mag)+" "+type.name;
		if (distance != null && distances.size() > 1) {
			curHeading += "#";
			uniqueTitleSection = distance+" km "+uniqueTitleSection;
			lines.add(curHeading+" "+uniqueTitleSection+" Results");
			lines.add(topLink); lines.add("");
		} else if (distance == null && type.separateDists && distances.size() > 1) {
			// this is total across multiple distances
			curHeading += "#";
			uniqueTitleSection = "All Distances "+uniqueTitleSection;
			lines.add(curHeading+" "+uniqueTitleSection+" Results");
			lines.add(topLink); lines.add("");
		}
		curHeading += "#";
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		if (sites.size() > 1)
			table.addColumn("Site");
		for (double period : periods) {
			String periodStr = optionalDigitDF.format(period)+"s";
			table.addColumn(periodStr+" "+type.htmlSymbol);
			if (type.stdDevOfMedians) {
				if (type.variedSymbolHTML != null) {
					table.addColumn("Mean "+type.variedSymbolHTML);
					table.addColumn(type.variedSymbolHTML+" Range");
				} else {
					table.addColumn("Mean ln(Value)");
					table.addColumn("ln(Value) Range");
				}
			} else {
				table.addColumn("Total");
				table.addColumn("Mean");
				table.addColumn("Median");
				table.addColumn("Range");
			}
		}
		table.finalizeLine();
		
		List<String> siteBundleTableNames = new ArrayList<>();
		List<String> siteBundlePlotLabels = new ArrayList<>();
		List<List<Site>> siteBundles = new ArrayList<>();
		List<PlotCurveCharacterstics> bundleCurveChars = new ArrayList<>();
		if (sites.size() > 1) {
			if (type.separateSites) {
				siteBundleTableNames.add("**ALL SITES**");
				siteBundlePlotLabels.add("Total");
				siteBundles.add(sites);
				bundleCurveChars.add(new PlotCurveCharacterstics(
						PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
			} else {
				List<Float> vs30s = new ArrayList<>(vs30SiteBundles.keySet());
				Collections.sort(vs30s, new Comparator<Float>() {

					@Override
					public int compare(Float o1, Float o2) {
						int size1 = vs30SiteBundles.get(o1).size();
						int size2 = vs30SiteBundles.get(o2).size();
						if (size1 == size2)
							return Float.compare(o1, o2);
						return Integer.compare(size2, size1); // greatest first
					}
					
				});
				CPT grayCPT = new CPT(0, Integer.max(1, vs30s.size()-1), Color.BLACK, Color.GRAY);
				for (int i=0; i<vs30s.size(); i++) {
					float vs30 = vs30s.get(i);
					String name = vs30SiteBundles.get(vs30).size()+" sites, V<sub>S30</sub>="+optionalDigitDF.format(vs30);
					siteBundleTableNames.add("**"+name+"**");
					siteBundlePlotLabels.add("Vs30="+optionalDigitDF.format(vs30));
					siteBundles.add(vs30SiteBundles.get(vs30));
					bundleCurveChars.add(new PlotCurveCharacterstics(
							PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, grayCPT.getColor((float)i)));
				}
			}
			if (!(type.singletonQuantities != null && Arrays.binarySearch(type.singletonQuantities, Quantity.SITE) >= 0)) {
				// don't do individual if SITE is a singleton
				for (int i=0; i<sites.size(); i++) {
					Site site = sites.get(i);
					siteBundleTableNames.add(site.name);
					siteBundlePlotLabels.add(site.name);
					List<Site> siteList = new ArrayList<>();
					siteList.add(site);
					siteBundles.add(siteList);
					bundleCurveChars.add(new PlotCurveCharacterstics(
							PlotLineType.SOLID, 1.5f, siteColors.get(i % siteColors.size())));
				}
			}
		} else {
			siteBundleTableNames.add(null);
			siteBundlePlotLabels.add("Total");
			siteBundles.add(sites);
			bundleCurveChars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		}

		List<List<ResidualStdDevSet>> bundleTotResidualStdDevs = new ArrayList<>();
//		List<List<MedianStdDevSet>> bundleTotMedianStdDevs = new ArrayList<>();
		List<DiscretizedFunc> bundleStdDevFuncs = new ArrayList<>();
		List<StdDevPercentileFuncs> bundlePercentilesList = null;
		for (int i=0; i<periods.length; i++) {
			bundleTotResidualStdDevs.add(new ArrayList<>());
//			bundleTotMedianStdDevs.add(new ArrayList<>());
		}

		PeriodDepResidualsList downsampledList = null;
		List<List<ASK_EventData>> realDataList = null;
		List<VariabilityResult> allVarResults = new ArrayList<>();
		List<String> dataSizeDependenceLines = null;
		
		DiscretizedFunc gmpeFunc = null;
		GMPE_Result gmpeResult = null;
		if (gmpes != null && type.gmpeStdDevType != null) {
//			System.out.println("Doing GMPE!");
			List<GMPE_Result> gmpeResults = new ArrayList<>();
			List<Double> myDists = new ArrayList<>();
			if (distance == null)
				for (Float dist : distances)
					myDists.add(dist.doubleValue());
			else
				myDists.add(distance.doubleValue());
			for (Site site : sites) {
				for (Double dist : myDists) {
					for (NGAW2_WrapperFullParam gmpe : gmpes) {
						GMPE_GroupingKey key = new GMPE_GroupingKey(gmpe, site, mag, dist);
						try {
							gmpeResults.add(calcGMPE(key));
						} catch (ExecutionException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					}
				}
			}
			if (D) System.out.println("Have "+gmpeResults.size()+" gmpe results");
			// pack to a single GMPE result
			gmpeResult = packGMs(gmpeResults.toArray(new GMPE_Result[0]));
			gmpeFunc = new ArbitrarilyDiscretizedFunc();
			double[] gmpeVars = gmpeResult.getVariability(type.gmpeStdDevType);
			for (int p=0; p<calcPeriods.length; p++)
				gmpeFunc.set(calcPeriods[p], gmpeVars[p]);
			gmpeFunc.setName("GMPE "+type.gmpeStdDevType.symbol);
//			System.out.println(gmpeFunc);
		} else {
//			System.out.println("NOT doing GMPEs. Type: "+type.gmpeStdDevType+". GMPEs null? "+(gmpes == null));
		}

		String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
		String magPrefix = mag == null ? "" : "_m"+optionalDigitDF.format(mag);
		
		// plot it
		String prefix = type.prefix+magPrefix+distPrefix+"_std_dev";
		
		if (type.stdDevOfMedians && distance != null) {
			CSVFile<String> eventTermsCSV = new CSVFile<>(true);
			Map<Integer, List<Double>> eventTerms = new HashMap<>();
			List<String> header = new ArrayList<>();
			header.add("Event ID");
			for (String siteName : siteBundlePlotLabels)
				header.add(siteName);
			eventTermsCSV.addLine(header);
			
			System.out.println("writing debug table for "+type);
			for (int i=0; i<siteBundles.size(); i++) {
				Site[] mySites = siteBundles.get(i).toArray(new Site[0]);
				VarGroupingKey key = new VarGroupingKey(type, mag, distance, mySites);
				PeriodDepResidualsList periodResiduals = loadResiduals(key);
				List<ResidualSet> residuals = periodResiduals.get(0);
				System.out.println("Calculating for "+siteBundlePlotLabels.get(i));
				VariabilityResult var = calcVarResult(key);
				System.out.println("\t"+type.symbol+" = "+var.medianStdDevs[0].stdDev);
				System.out.println("\tmean median = "+var.medianStdDevs[0].meanMedian);
				System.out.println("\tmin median = "+var.medianStdDevs[0].minMedian);
				System.out.println("\tmax median = "+var.medianStdDevs[0].maxMedian);
				if (var.logVals != null)
					System.out.println("\tnum vals = "+var.logVals[0].rawVals.length);
				
				for (ResidualSet set : residuals) {
					Integer eventID = set.commonRotSpec.eventID;
					List<Double> terms;
					if (i == 0) {
						Preconditions.checkState(!eventTerms.containsKey(eventID));
						terms = new ArrayList<>();
						eventTerms.put(eventID, terms);
					} else {
						terms = eventTerms.get(eventID);
					}
					terms.add(set.median);
				}
			}
			for (int eventID : config.getValues(Integer.class, Quantity.EVENT_ID)) {
				List<String> line = new ArrayList<>();
				line.add(eventID+"");
				List<Double> terms = eventTerms.get(eventID);
				for (Double term : terms)
					line.add(term.floatValue()+"");
				eventTermsCSV.addLine(line);
			}
			File outputFile = new File(resourcesDir, prefix+"_"+optionalDigitDF.format(periods[0])+"s_medians_table.csv");
			System.out.println("writing "+outputFile.getName());
			eventTermsCSV.writeToFile(outputFile);
//			System.exit(0);
//		} else {
//			return new ArrayList<>();
		}
		
		for (int i=0; i<siteBundles.size(); i++) {
			Site[] mySites = siteBundles.get(i).toArray(new Site[0]);
			table.initNewLine();
			String siteTableName = siteBundleTableNames.get(i);
			String emphasis = "";
			if (sites.size() > 1) {
				for (int j=0; j<siteTableName.length() && siteTableName.charAt(j) == '*'; j++)
					emphasis += '*';
				table.addColumn(siteTableName);
			}
			DiscretizedFunc stdDevFunc = new ArbitrarilyDiscretizedFunc();
			StdDevPercentileFuncs stdDevPercentiles = null;
			VarGroupingKey key = new VarGroupingKey(type, mag, distance, mySites);
			VariabilityResult result; 
			try {
				result = varResultCache.get(key);
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (type.stdDevOfMedians) {
				System.out.println("Calculating for "+siteBundlePlotLabels.get(i));
				System.out.println("\t"+type.symbol+" = "+result.medianStdDevs[0].stdDev);
				System.out.println("\tmean median = "+result.medianStdDevs[0].meanMedian);
				System.out.println("\tmin median = "+result.medianStdDevs[0].minMedian);
				System.out.println("\tmax median = "+result.medianStdDevs[0].maxMedian);
				if (result.logVals != null)
					System.out.println("\tnum vals = "+result.logVals[0].rawVals.length);
			}
//			System.exit(0);
			for (int p=0; p<periods.length; p++) {
				if (type.stdDevOfMedians) {
					MedianStdDevSet set = result.getMedianStdDevSet(periods[p]);
					table.addColumn(emphasis+optionalDigitDF.format(set.stdDev)+emphasis);
					table.addColumn(emphasis+optionalDigitDF.format(set.meanMedian)+emphasis);
					table.addColumn(emphasis+"["+optionalDigitDF.format(set.minMedian)
						+" "+optionalDigitDF.format(set.maxMedian)+"]"+emphasis);
				} else {
					table.addColumn(""); // empty for period
					ResidualStdDevSet set = result.getResidualStdDevSet(periods[p]);
					table.addColumn(emphasis+optionalDigitDF.format(set.total)+emphasis);
					table.addColumn(emphasis+optionalDigitDF.format(set.mean)+emphasis);
					table.addColumn(emphasis+optionalDigitDF.format(set.median)+emphasis);
					table.addColumn(emphasis+"["+optionalDigitDF.format(set.min)
						+" "+optionalDigitDF.format(set.max)+"]"+emphasis);
				}
			}
			// now do the functions for all periods
			for (double period : calcPeriods) {
				if (type.stdDevOfMedians) {
					MedianStdDevSet set = result.getMedianStdDevSet(period);
					stdDevFunc.set(period, set.stdDev);
				} else {
					ResidualStdDevSet set = result.getResidualStdDevSet(period);
					stdDevFunc.set(period, set.total);
				}
			}
			// add independent
			if (type.stdDevOfMedians) {
				MedianStdDevSet set = result.getPeriodIndepMedianStdDevSet();
				stdDevFunc.set(0d, set.stdDev);
			} else {
				ResidualStdDevSet set = result.getPeriodIndepResidualStdDevSet();
				stdDevFunc.set(0d, set.total);

				for (int p=0; p<periods.length; p++) {
					bundleTotResidualStdDevs.get(p).add(result.getResidualStdDevSet(periods[p]));
				}
			}
			if (realEventData != null) {
				VariabilityResult[] downsampledResults;
				try {
					downsampledResults = downsampledVarResultCache.get(key);
				} catch (ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				if (downsampledResults.length > 0) {
					for (VariabilityResult v : downsampledResults)
						allVarResults.add(v);
					if (downsampledList == null) {
						VarGroupingKey sampleKey = new VarGroupingKey(type, mag, distance, mySites);
						realDataList = getRealDataListForKey(sampleKey);
						if (distance == null && type.separateDists)
							// we're doing this across all distances, so get the rotations for 1.
							// we'll scale to the number of distances later
							sampleKey = new VarGroupingKey(type, mag, distances.get(0), mySites);
						if (type.separateSites && mySites.length > 1)
							sampleKey = new VarGroupingKey(type, mag, sampleKey.distance, mySites[0]);
						sampleKey.realDataList = realDataList;
						downsampledList = loadResiduals(sampleKey);
						bundlePercentilesList = new ArrayList<>();
					}
					DiscretizedFunc lower95Func = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc upper95Func = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc lower68Func = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc upper68Func = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc medianFunc = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc sigmaFunc = new ArbitrarilyDiscretizedFunc();
					List<Double> myCalcPeriods = new ArrayList<>();
					myCalcPeriods.add(0d);
					for (double period : calcPeriods)
						myCalcPeriods.add(period);
					for (double period : myCalcPeriods) {
						double[] values = new double[downsampledResults.length];
						for (int j=0; j<values.length; j++) {
							if (type.stdDevOfMedians) {
								MedianStdDevSet set = period == 0d ? downsampledResults[j].getPeriodIndepMedianStdDevSet() :
									downsampledResults[j].getMedianStdDevSet(period);
								Preconditions.checkState(set.medians == null || set.medians.length > 1);
								values[j] = set.stdDev;
							} else {
								ResidualStdDevSet set = period == 0d ? downsampledResults[j].getPeriodIndepResidualStdDevSet() :
									downsampledResults[j].getResidualStdDevSet(period);
								values[j] = set.total;
							}
//							if (period == 1d) System.out.println("1s "+j+": "+values[j]);
						}
						lower95Func.set(period, StatUtils.percentile(values, 2.5));
						upper95Func.set(period, StatUtils.percentile(values, 97.5));
						lower68Func.set(period, StatUtils.percentile(values, 16));
						upper68Func.set(period, StatUtils.percentile(values, 84));
						medianFunc.set(period, StatUtils.percentile(values, 50d));
						sigmaFunc.set(period, Math.sqrt(StatUtils.variance(values)));
					}
//					for (double period : calcPeriods)
//						System.out.println((float)period+"\t"+optionalDigitDF.format(stdDevFunc.getY(period))+"\t["
//								+optionalDigitDF.format(lowerFunc.getY(period))+" "+optionalDigitDF.format(upperFunc.getY(period))+"]");
					
					stdDevPercentiles = new StdDevPercentileFuncs(new UncertainArbDiscDataset(medianFunc, lower95Func, upper95Func),
							new UncertainArbDiscDataset(medianFunc, lower68Func, upper68Func), sigmaFunc);
					
					if (type == VariabilityType.BETWEEN_EVENTS || type == VariabilityType.WITHIN_EVENT_SS)
						dataSizeDependenceLines = plotDataSizeDependence(key, resourcesDir, stdDevFunc, gmpeResult);
				}
			}
			stdDevFunc.setName(siteBundlePlotLabels.get(i));
			bundleStdDevFuncs.add(stdDevFunc);
			if (stdDevPercentiles != null)
				bundlePercentilesList.add(stdDevPercentiles);
			table.finalizeLine();
		}
		plotPeriodDependentStdDevs(resourcesDir, prefix, type.name+" ("+type.symbol+")",
				bundleStdDevFuncs, bundlePercentilesList, bundleCurveChars, gmpeFunc);
		lines.add("!["+type.name+" Variability](resources/"+prefix+".png)");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		if (!type.stdDevOfMedians) {
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+"s");
			table.finalizeLine();
			table.initNewLine();
			for (int p=0; p<periods.length; p++) {
				prefix = type.prefix+magPrefix+distPrefix+"_"+optionalDigitDF.format(periods[p])+"s_hist";
				plotStdDevsHistogram(resourcesDir, prefix, optionalDigitDF.format(periods[p])+"s "+type.name+" ("+type.symbol+")",
						periods[p], bundleTotResidualStdDevs.get(p).get(0), null);
				table.addColumn("!["+optionalDigitDF.format(periods[p])+"s](resources/"+prefix+".png)");
			}
			table.finalizeLine();
			table.wrap(periods.length > 4 ? 3 : 2, 0);
			lines.add("Here are plots of the histogram of "+type.htmlSymbol+" for each individual rupture, from "
					+ "which we compute a total "+type.htmlSymbol);
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
		}
		if (downsampledList != null && bundlePercentilesList != null && bundlePercentilesList.get(0) != null) {
			lines.add(curHeading+" "+uniqueTitleSection+" Downsampled Results");
			lines.add(topLink); lines.add("");
			
			DiscretizedFunc totStdDevFunc = bundleStdDevFuncs.get(0);
			StdDevPercentileFuncs totStdDevPercentiles = bundlePercentilesList.get(0);
			String distLine;
			if (distance == null) {
				distLine = "and all distances";
			} else {
				float distAdd = distance < 100f ? 10 : 20;
				distLine = "and distance within the range ["+(distance-distAdd)+" "+(distance+distAdd)+"] km";
			}
			int numEventsUsed = 0;
			int numRecordingsUsed = downsampledList.rotations.size();
			if (distance == null)
				numRecordingsUsed *= distances.size();
			int possibleRecordings = 0;
			for (List<ASK_EventData> data : realDataList) {
				if (data.size() > 1) {
					numEventsUsed++;
					possibleRecordings += data.size();
				}
			}
			String samplesPart = numRealDataSamples+" times";
			if (sites.size() > 1)
				samplesPart += " for each site";
			lines.add("We compute uncertainties on "+type.htmlSymbol+" through downsampling the rotational synthetic data to match the sample "
					+ "sizes used in the ASK 2014 regressions. We search the ASK dataset for ruptures with the same mechanism, magnitude in the "
					+ "range ["+(float)(mag-0.2)+" "+(float)(mag+0.2)+"], "+distLine+". We throw out any events with only 1 recording, leaving "
					+ "us with "+numEventsUsed+" events and a total of "+numRecordingsUsed+" recordings. We then downsample our "
					+ "simulated data "+samplesPart+", and compute "+type.htmlSymbol+" from each sample. The 95% confidence range "
					+ "from these samples is plotted as a shaded region above, and listed in the table below. Weighted standard deviations are "
					+ "calculated, weighted by the square-root of the number of recordings in each event.");
			lines.add("");
			if (possibleRecordings > numRecordingsUsed) {
				lines.add("*WARNING: Some real events had more recordings than we have rotations per event, so our dataset for this test "
						+ "is smaller. We are using "+(possibleRecordings-numRecordingsUsed)+" fewer data points.*");
				lines.add("");
			}
			table = MarkdownUtils.tableBuilder();
			table.addLine("Period (s)", "Full "+type.htmlSymbol, "Downsampled median "+type.htmlSymbol,
					"Downsampled "+type.htmlSymbol+" std. dev.", "Downsampled "+type.htmlSymbol+" 68% conf range",
					"Downsampled "+type.htmlSymbol+" 95% conf range");
			List<Double> myCalcPeriods = new ArrayList<>();
			myCalcPeriods.add(0d);
			for (double period : calcPeriods)
				myCalcPeriods.add(period);
			for (double period : myCalcPeriods) {
				String periodStr = period == 0d ? "T-independent" : optionalDigitDF.format(period);
				table.addLine(periodStr, optionalDigitDF.format(totStdDevFunc.getY(period)),
						optionalDigitDF.format(totStdDevPercentiles.bounds68.getY(period)),
						optionalDigitDF.format(totStdDevPercentiles.sigmaFunc.getY(period)),
						"["+optionalDigitDF.format(totStdDevPercentiles.bounds68.getLowerY(period))
						+" "+optionalDigitDF.format(totStdDevPercentiles.bounds68.getUpperY(period))+"]",
						"["+optionalDigitDF.format(totStdDevPercentiles.bounds95.getLowerY(period))
						+" "+optionalDigitDF.format(totStdDevPercentiles.bounds95.getUpperY(period))+"]");
			}
			lines.addAll(table.build());
			lines.add("");
			// now plot site histograms
			Double[] dsPeriods = { null, periods[0] };
			if (type.separateSites) {
				lines.add("These plots show the distribution of period-independent downsampled "+type.htmlSymbol
						+" for each site.");
				lines.add("");
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				table.addColumn("Period");
				if (sites.size() > 1 || sites.get(0) != null) {
					for (Site site : sites)
						table.addColumn("**"+site.getName()+"**");
				}
				table.finalizeLine();
				List<File[]> histFiles = new ArrayList<>();
				for (Site site : sites) {
					String sitePrefix = site == null ? "" : "_"+site.getName();
					prefix = type.prefix+magPrefix+distPrefix+sitePrefix+"_downsampled_hist";
					histFiles.add(plotSiteDownsampledHistogram(resourcesDir, prefix, type, mag, distance, dsPeriods, site));
				}
				for (int i=0; i<dsPeriods.length; i++) {
					Double period = dsPeriods[i];
					table.initNewLine();
					table.addColumn(period == null ? "Period-Indep" : optionalDigitDF.format(period)+"s");
					for (File[] files : histFiles)
						table.addColumn("![Dowmsampled Histogram](resources/"+files[i].getName()+")");
					table.finalizeLine();
				}
				lines.addAll(table.wrap(4, 1).build());
				lines.add("");
			} else {
				lines.add("This plot shows the distribution of period-independent downsampled "+type.htmlSymbol+".");
				lines.add("");
				prefix = type.prefix+magPrefix+distPrefix+"_downsampled_hist";
				File[] rets = plotSiteDownsampledHistogram(resourcesDir, prefix, type, mag, distance,
						dsPeriods, siteBundles.get(0).toArray(new Site[0]));
				table = MarkdownUtils.tableBuilder();
				for (int i=0; i<dsPeriods.length; i++) {
					Double period = dsPeriods[i];
					table.initNewLine();
					table.addColumn(period == null ? "Period-Indep" : optionalDigitDF.format(period)+"s");
					table.addColumn("![Dowmsampled Histogram](resources/"+rets[i].getName()+")");
					table.finalizeLine();
				}
				lines.addAll(table.build());
				lines.add("");
			}
			if (dataSizeDependenceLines != null) {
				lines.addAll(dataSizeDependenceLines);
				lines.add("");
			}
		}
		if (type.disaggScatters != null && type.disaggScatters.length > 0 && (distance != null || !type.separateDists)) {
			lines.add(curHeading+" "+uniqueTitleSection+" Disaggregations");
			lines.add(topLink); lines.add("");
			System.out.println("*** Plotting disagg scatters ***");
			VariabilityResult result;
			try {
				result = varResultCache.get(new VarGroupingKey(type, mag, distance, sites.toArray(new Site[0])));
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			lines.add("Here are plots of the "+type.htmlSymbol+" as a function of various parameters for disaggregation.");
			lines.add("");
			for (ScatterDisaggQuantity scatterQuantity : type.disaggScatters) {
				System.out.println("*** Plotting disagg scatter of "+scatterQuantity+" ***");
				File[][] plots = plotScatter(resourcesDir, type.prefix+"_scatter_", periods, type, result, scatterQuantity);
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				for (double period : periods)
					table.addColumn(optionalDigitDF.format(period)+"s");
				table.finalizeLine();
				for (File[] periodPlots : plots) {
					if (periodPlots == null || periodPlots[0] == null)
						continue;
					table.initNewLine();
					for (File plot : periodPlots)
						table.addColumn("![Scatter](resources/"+plot.getName()+")");
					table.finalizeLine();
				}
				lines.addAll(table.build());
				lines.add("");
				System.out.println("*** DONE disagg scatter of "+scatterQuantity+" ***");
			}
			System.out.println("*** DONE disagg scatters ***");
		}
		return lines;
	}
	
	private boolean replotDependence = false;
	private boolean dependUseSimCounts = true;
	
	private List<String> plotDataSizeDependence(VarGroupingKey key, File resourcesDir, DiscretizedFunc simFunc, GMPE_Result gmpeVal)
			throws IOException {
		List<List<ASK_EventData>> realData;

		double dm = 0.2;
		Map<Integer, List<ASK_EventData>> totalDataMatches = ASK_EventData.getMatches(realEventData,
				com.google.common.collect.Range.closed(key.magnitude-dm, key.magnitude+dm),
				null, null, 0d);
		List<List<ASK_EventData>> subDataMatches = getRealDataListForKey(key);
		
		if (dependUseSimCounts) {
			// we want actual num recordings and num events for all matching ruptures here,
			// regardless of distance
			
			realData = new ArrayList<>();
			for (List<ASK_EventData> values : totalDataMatches.values())
				realData.add(values);
		} else {
			realData = subDataMatches;
		}
		if (realData == null || realData.isEmpty())
			return null;
		int numEvents = realData.size();
		double aveNumRecordings = 0d;
		for (List<ASK_EventData> data : realData)
			aveNumRecordings += data.size();
		aveNumRecordings /= (double)numEvents;
		numEvents = Integer.max(numEvents, 2);
		int roundedNumRecordings = Integer.max(2, (int)Math.round(aveNumRecordings));
//		RotatedRupVariabilityConfig config = magConfigs.get(key.magnitude);
//		int simNumEvents = config.getValues(Integer.class, Quantity.EVENT_ID).size();
//		int simNumRecordings = 1;
		VariabilityType type = key.type;
		List<VarGroupingKey> keys = new ArrayList<>();
		keys.add(new VarGroupingKey(type, key.magnitude, key.distance, key.sites));
//		for (Quantity varyQ : type.variedQuantities)
//			simNumRecordings *= config.getQuantitiesMap().get(varyQ).size();
//		int maxNumEvents = Integer.max(numEvents, simNumEvents);
		int maxNumEvents = Integer.max(50, 2*numEvents);
		
		int numRealizations = 100;
		ArbitrarilyDiscretizedFunc numEventsFunc = new ArbitrarilyDiscretizedFunc();
		int num = 2;
		while (num<=maxNumEvents) {
			numEventsFunc.set((double)num, 0d);
			if (num>=100)
				num += 20;
			else if (num > 50)
				num += 10;
			else if (num > 20)
				num += 5;
			else if (num > 10)
				num += 2;
			else
				num++;
		}
		if (numEventsFunc.getMaxX() < maxNumEvents)
			numEventsFunc.set((double)maxNumEvents, 0d);
		
		RotatedRupVariabilityConfig<E> config = magConfigs.get(key.magnitude);
		Map<Quantity, List<?>> qMap = config.getQuantitiesMap();
		int fixedNumRecordings;
		if (dependUseSimCounts) {
			fixedNumRecordings = 1;
			for (Quantity q : type.variedQuantities)
				fixedNumRecordings *= qMap.get(q).size();
		} else {
			fixedNumRecordings = roundedNumRecordings;
		}
		int fixedNumEvents;
		if (dependUseSimCounts) {
			fixedNumEvents = qMap.get(Quantity.EVENT_ID).size();
		} else {
			fixedNumEvents = numEvents;
		}
		
//		int maxNumRecordings = Integer.max(roundedNumRecordings, simNumRecordings);
		int maxNumRecordings = Integer.max(50, 2*roundedNumRecordings);
		ArbitrarilyDiscretizedFunc numRecordingsFunc = new ArbitrarilyDiscretizedFunc();
		num = 2;
		while (num<=maxNumRecordings) {
			numRecordingsFunc.set((double)num, 0d);
			if (num>=100)
				num += 20;
			else if (num > 50)
				num += 10;
			else if (num > 20)
				num += 5;
			else if (num > 10)
				num += 2;
			else
				num++;
		}
		if (numRecordingsFunc.getMaxX() < maxNumRecordings)
			numRecordingsFunc.set((double)maxNumRecordings, 0d);
		
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Period", "Event Count Dependence", "Recordings/Event Dependence");
		Double[] periods = new Double[] { null, calcPeriods[0] };
		for (Double period : periods) {
			table.initNewLine();
			String prefix = type.prefix+"_event_count_dependence";
			if (key.distance == null)
				prefix += "_all_dists";
			else
				prefix += "_"+optionalDigitDF.format(key.distance)+"km";
			double simVal;
			if (period == null) {
				prefix += "_period_indep";
				table.addColumn("Period Indep.");
				simVal = simFunc.getY(0d);
			} else {
				prefix += "_"+optionalDigitDF.format(period)+"s";
				table.addColumn(optionalDigitDF.format(period)+"s");
				simVal = simFunc.getY(period);
			}
			if (replotDependence || !(new File(resourcesDir, prefix+".png").exists())) {
				System.out.println("Calculating event count dependence (max="+maxNumEvents+"), period="+period);
				System.out.println("\tNum recordings each for event count dependence: "+fixedNumRecordings);
				DiscretizedFunc lower95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc lower68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc sigmaFunc = new ArbitrarilyDiscretizedFunc();
				List<Future<double[]>> futures = new ArrayList<>();
				for (int i=0; i<numEventsFunc.size(); i++) {
					List<List<ASK_EventData>> realDataList = new ArrayList<>();
					int myNumEvents = (int)Math.round(numEventsFunc.getX(i));
					for (int j=0; j<myNumEvents; j++) {
						List<ASK_EventData> dataList = new ArrayList<>();
						for (int k=0; k<fixedNumRecordings; k++)
							dataList.add(null); // only the size of this list is used, not the contents. null is ok
						realDataList.add(dataList);
					}
					
					futures.add(exec.submit(new DependCalcCallable(key.type, numRealizations, keys, realDataList, period)));
				}
				for (int i=0; i<numEventsFunc.size(); i++) {
					int myNumEvents = (int)Math.round(numEventsFunc.getX(i));
					System.out.print(myNumEvents+" ");
					
					double[] results;
					try {
						results = futures.get(i).get();
					} catch (Exception e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					
					double x = numEventsFunc.getX(i);
					lower95Func.set(x, StatUtils.percentile(results, 2.5));
					upper95Func.set(x, StatUtils.percentile(results, 97.5));
					lower68Func.set(x, StatUtils.percentile(results, 16));
					upper68Func.set(x, StatUtils.percentile(results, 84));
					numEventsFunc.set(x, StatUtils.percentile(results, 50d));
					sigmaFunc.set(x, Math.sqrt(StatUtils.variance(results)));
				}
				System.out.println();
				StdDevPercentileFuncs numEventsPercentiles = new StdDevPercentileFuncs(new UncertainArbDiscDataset(numEventsFunc, lower95Func, upper95Func),
						new UncertainArbDiscDataset(numEventsFunc, lower68Func, upper68Func), sigmaFunc);
				Map<String, Double> dataVals = new HashMap<>();
				dataVals.put("ASK (2014) Num Events", (double)numEvents);
				plotDataSizeDependence(numEventsFunc, numEventsPercentiles, type.name+" ("+type.symbol+") Event Count Dependence",
						"Num Events", type.symbol, dataVals, simVal, gmpeVal, resourcesDir, prefix, type);
			}
			
			table.addColumn("![num events dependence]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			prefix = type.prefix+"_event_recordings_dependence";
			if (key.distance == null)
				prefix += "_all_dists";
			else
				prefix += "_"+optionalDigitDF.format(key.distance)+"km";
			if (period == null)
				prefix += "_period_indep";
			else
				prefix += "_"+optionalDigitDF.format(period);
			if (replotDependence || !(new File(resourcesDir, prefix+".png").exists())) {
				System.out.println("Calculating recording count dependence (max="+maxNumRecordings+"), period="+period);
				System.out.println("\tNum events each for recording count dependence: "+fixedNumEvents);
				DiscretizedFunc lower95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc lower68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc sigmaFunc = new ArbitrarilyDiscretizedFunc();
				List<Future<double[]>> futures = new ArrayList<>();
				for (int i=0; i<numRecordingsFunc.size(); i++) {
					List<List<ASK_EventData>> realDataList = new ArrayList<>();
					int myNumRecordings = (int)Math.round(numRecordingsFunc.getX(i));
					for (int j=0; j<fixedNumEvents; j++) {
						List<ASK_EventData> dataList = new ArrayList<>();
						for (int k=0; k<myNumRecordings; k++)
							dataList.add(null); // only the size of this list is used, not the contents. null is ok
						realDataList.add(dataList);
					}

					futures.add(exec.submit(new DependCalcCallable(key.type, numRealizations, keys, realDataList, period)));
				}
				for (int i=0; i<numRecordingsFunc.size(); i++) {
					int myNumRecordings = (int)Math.round(numRecordingsFunc.getX(i));
					System.out.print(myNumRecordings+" ");

					double[] results;
					try {
						results = futures.get(i).get();
					} catch (Exception e) {
						throw ExceptionUtils.asRuntimeException(e);
					}

					double x = numRecordingsFunc.getX(i);
					lower95Func.set(x, StatUtils.percentile(results, 2.5));
					upper95Func.set(x, StatUtils.percentile(results, 97.5));
					lower68Func.set(x, StatUtils.percentile(results, 16));
					upper68Func.set(x, StatUtils.percentile(results, 84));
					numRecordingsFunc.set(x, StatUtils.percentile(results, 50d));
					sigmaFunc.set(x, Math.sqrt(StatUtils.variance(results)));
				}
				System.out.println();
				StdDevPercentileFuncs numRecordingsPercentiles = new StdDevPercentileFuncs(new UncertainArbDiscDataset(numRecordingsFunc, lower95Func, upper95Func),
						new UncertainArbDiscDataset(numRecordingsFunc, lower68Func, upper68Func), sigmaFunc);
				Map<String, Double> dataVals = new HashMap<>();
				if (key.distance == null) {
					dataVals.put("ASK2014 Recs/Event", aveNumRecordings);
				} else {					
					double aveTotRecs = 0d;
					for (List<ASK_EventData> data : totalDataMatches.values())
						aveTotRecs += data.size();
					aveTotRecs /= totalDataMatches.size();
					double aveDistRecs = 0d;
					for (List<ASK_EventData> data : subDataMatches)
						aveDistRecs += data.size();
					aveDistRecs /= subDataMatches.size();
					dataVals.put("ASK2014 Râ‰…"+optionalDigitDF.format(key.distance)+" Recs/Event", aveDistRecs);
					dataVals.put("Total Recs/Event", aveTotRecs);
				}
				
				plotDataSizeDependence(numRecordingsFunc, numRecordingsPercentiles, type.name+" ("+type.symbol+") Event Recordings Dependence",
						"Num Recordings Per Event", type.symbol, dataVals, simVal, gmpeVal, resourcesDir, prefix, type);
			}
			table.addColumn("![num recordings dependence]("+resourcesDir.getName()+"/"+prefix+".png)");
			table.finalizeLine();
		}
		
		exec.shutdown();
		
		List<String> lines = new ArrayList<>();
		String line = "These plots show the dependence of "+type.htmlSymbol+" to the number of events included and the number of recordings "
				+ "per event. The left plot holds the number of recordings per event fixed at the ";
		if (dependUseSimCounts)
			line += "full set of simulated recordings";
		else
			line += "average data value";
		line += " ("+fixedNumRecordings+"), varying the number of events. The right plot holds the number of events fixed at the ";
		if (dependUseSimCounts)
			line += "full set of simulated events";
		else
			line += "data value";
		line += " ("+fixedNumEvents+"), varying the number of recordings per event.";
		lines.add(line);
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		float lowerMag = (float)(key.magnitude - dm);
		float upperMag = (float)(key.magnitude + dm);
		line = "This is a histogram of the number of recordings per event from ASK 2014 with M=["+lowerMag+","+upperMag+"].";
		String prefix = type.prefix+"_event_recordings_hist";
		if (key.distance == null) {
			prefix += "_all_dists";
			ASK_EventData.plotMultiCountHist(resourcesDir, prefix, "ASK 2014 Recordings Distribution",
					totalDataMatches.values(), null, null, null);
		} else {
			prefix += "_"+optionalDigitDF.format(key.distance)+"km";
			double dr = getDeltaDistance(key.distance);
			float lower = key.distance - (float)dr;
			float upper = key.distance + (float)dr;
			ASK_EventData.plotMultiCountHist(resourcesDir, prefix, "ASK 2014 Recordings Distribution",
					subDataMatches, "Dist âˆˆ ["+lower+","+upper+"] km", totalDataMatches.values(), "All Distances");
			line += " The top plot shows the subset with distance in the range ["+lower+","+upper+"], "
					+ "and the bottom the whole distribution at all distances.";
		}
		lines.add(line);
		lines.add("");
		lines.add("![Histogram]("+resourcesDir.getName()+"/"+prefix+".png)");
		
		return lines;
	}
	
	private class DependCalcCallable implements Callable<double[]> {
		
		private VariabilityType type;
		private int numRealizations;
		private List<VarGroupingKey> keys;
		private List<List<ASK_EventData>> realDataList;
		private Double period;

		public DependCalcCallable(VariabilityType type, int numRealizations, List<VarGroupingKey> keys,
				List<List<ASK_EventData>> realDataList, Double period) {
			this.type = type;
			this.numRealizations = numRealizations;
			this.keys = keys;
			this.realDataList = realDataList;
			this.period = period;
		}

		@Override
		public double[] call() throws Exception {
			List<Double> resultsList = new ArrayList<>();
			
			for (int n=0; n<numRealizations; n++) {
				// each call will be a new stochastic realization due to custom equals/hashCode implementation with realDataList set
				List<VariabilityResult> results = new ArrayList<>();
				for (VarGroupingKey myKey : keys) {
					VarGroupingKey randKey = new VarGroupingKey(myKey.type, myKey.magnitude, myKey.distance, myKey.sites);
					randKey.realDataList = realDataList;
					results.add(calcVarResult(randKey));
				}
				for (VariabilityResult result : results) {
					if (period == null) {
						if (type.stdDevOfMedians)
							resultsList.add(result.getPeriodIndepMedianStdDevSet().stdDev);
						else
							resultsList.add(result.getPeriodIndepResidualStdDevSet().total);
					} else {
						if (type.stdDevOfMedians)
							resultsList.add(result.getMedianStdDevSet(period).stdDev);
						else
							resultsList.add(result.getResidualStdDevSet(period).total);
					}
					
				}
			}
			return Doubles.toArray(resultsList);
		}
		
	}
	
	private void plotDataSizeDependence(DiscretizedFunc medianFunc, StdDevPercentileFuncs percentiles,
			String title, String xAxisLabel, String yAxisLabel, Map<String, Double> dataVals,
			double simVal, GMPE_Result gmpeVal, File resourcesDir, String prefix, VariabilityType type)
					throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

//		Color color95 = new Color(20, 20, 20);
//		Color color68 = new Color(60, 60, 60);
		Color color95 = Color.LIGHT_GRAY;
		Color color68 = Color.GRAY;
		
		percentiles.bounds95.setName("95% Range");
		funcs.add(percentiles.bounds95);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, color95));
//		percentiles.bounds68.setName("68% Range");
//		funcs.add(percentiles.bounds68);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, color68));
		medianFunc.setName("Median");
		funcs.add(medianFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		if (dataVals != null && !dataVals.isEmpty()) {
			PlotLineType[] dataTypes = { PlotLineType.DOTTED, PlotLineType.DASHED, PlotLineType.DOTTED_AND_DASHED };
			List<String> dataNames = ComparablePairing.getSortedData(dataVals);
			for (int i=0; i<dataNames.size(); i++) {
				String dataName = dataNames.get(i);
				double dataX = dataVals.get(dataName);
				XY_DataSet dataFunc = new DefaultXY_DataSet();
				dataFunc.setName(dataName);
				dataFunc.set(dataX, 0d);
				dataFunc.set(dataX, 1d);
				funcs.add(dataFunc);
				chars.add(new PlotCurveCharacterstics(dataTypes[i % dataTypes.length], 3f, Color.RED));
			}
		}
		
		if (gmpeVal != null) {
			double gmpeY = gmpeVal.medianVars.get(type.gmpeStdDevType)[0];
			XY_DataSet gmpeFunc = new DefaultXY_DataSet();
			if (gmpes.length > 1)
				gmpeFunc.setName("GMPE "+type.gmpeStdDevType.symbol);
			else
				gmpeFunc.setName(gmpes[0].getShortName()+" "+type.gmpeStdDevType.symbol);
			gmpeFunc.set(0d, gmpeY);
			gmpeFunc.set(medianFunc.getMaxX(), gmpeY);
			funcs.add(gmpeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLUE));
		}
		
		XY_DataSet simFunc = new DefaultXY_DataSet();
		simFunc.setName("Simulated "+type.symbol);
		simFunc.set(0d, simVal);
		simFunc.set(medianFunc.getMaxX(), simVal);
		funcs.add(simFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GREEN.darker()));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		Range xRange = new Range(0d, medianFunc.getMaxX());
		Range yRange = new Range(0d, 1d);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
	}

	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	protected static Float nullAsZero(Float value) {
		if (value == null)
			return 0f;
		return value;
	}
	
	private static class ResidualSet {
		double[] residuals;
		double median;
		double residualStdDev;
		
		private RotationSpec commonRotSpec;
		
		public ResidualSet(List<RotationSpec> rotations, double[] values) {
			this.residuals = values; // will subtract median in a sec
			this.median = DataUtils.median(values);
			Preconditions.checkState(Double.isFinite(median), "Non-finite median: %s", median);
			for (int i=0; i<values.length; i++)
				residuals[i] -= median;
			double variance = StatUtils.variance(residuals, 0d);
			if (variance < 0) {
				Preconditions.checkState(variance > -1e-10, "Negative variance: %s", variance);
				variance = 0;
			}
			this.residualStdDev = Math.sqrt(variance);
			Preconditions.checkState(Double.isFinite(residualStdDev),
					"Non-finite std dev: %s, var=%s, median=%s, size=%s", residualStdDev, variance, median, values.length);
			
			for (RotationSpec rotation : rotations) {
				if (commonRotSpec == null) {
					commonRotSpec = new RotationSpec(-1, rotation.site, rotation.eventID, rotation.distance,
							rotation.sourceAz, rotation.siteToSourceAz);
				} else {
					if (!Objects.equals(commonRotSpec.site, rotation.site))
						commonRotSpec = new RotationSpec(-1, null, commonRotSpec.eventID, commonRotSpec.distance,
								commonRotSpec.sourceAz, commonRotSpec.siteToSourceAz);
					if (commonRotSpec.eventID != rotation.eventID)
						commonRotSpec = new RotationSpec(-1, commonRotSpec.site, -1, commonRotSpec.distance,
								commonRotSpec.sourceAz, commonRotSpec.siteToSourceAz);
					if (!Objects.equals(commonRotSpec.distance, rotation.distance))
						commonRotSpec = new RotationSpec(-1, commonRotSpec.site, commonRotSpec.eventID, null,
								commonRotSpec.sourceAz, commonRotSpec.siteToSourceAz);
					if (!Objects.equals(commonRotSpec.sourceAz, rotation.sourceAz))
						commonRotSpec = new RotationSpec(-1, commonRotSpec.site, commonRotSpec.eventID, commonRotSpec.distance,
								null, commonRotSpec.siteToSourceAz);
					if (!Objects.equals(commonRotSpec.siteToSourceAz, rotation.siteToSourceAz))
						commonRotSpec = new RotationSpec(-1, commonRotSpec.site, commonRotSpec.eventID, commonRotSpec.distance,
								commonRotSpec.sourceAz, null);
				}
			}
		}
		
		public int size() {
			return residuals.length;
		}
		
		public double getMedian() {
			return median;
		}
		
		public double[] getResiduals() {
			return residuals;
		}
		
		public double[] getLogValues() {
			double[] ret = new double[residuals.length];
			for (int i=0; i<ret.length; i++)
				ret[i] = residuals[i]+median;
			return ret;
		}
		
		public double getStdDev() {
			return residualStdDev;
		}
	}
	
	private class PeriodDepResidualsList extends ArrayList<List<ResidualSet>> {
		
		private final List<RotationSpec> rotations;

		private PeriodDepResidualsList(List<RotationSpec> rotations) {
			this.rotations = rotations;
		}
	}
	
	private PeriodDepResidualsList loadResiduals(VarGroupingKey key) throws IOException {
		// do for all sites in the key
		if (!key.separateSites && (key.restrictTo == null || !key.restrictTo.containsKey(Quantity.SITE))) {
			// restrict it to just the specific sites
			Quantity[] separateQuantities = new Quantity[key.separateQuantities.length+1];
			separateQuantities[0] = Quantity.SITE;
			System.arraycopy(key.separateQuantities, 0, separateQuantities, 1, key.separateQuantities.length);
			VarGroupingKey newKey = new VarGroupingKey(key.type, separateQuantities, key.groupQuantities,
					key.magnitude, key.distance, key.sites);
			newKey.restrictTo = key.restrictTo;
			newKey.realDataList = key.realDataList;
			key = newKey;
		}
		Quantity[] constQuantities = key.separateQuantities;
		Object[] constValues = new Object[constQuantities.length];
		for (int i=0; i<constQuantities.length; i++) {
			if (constQuantities[i] == Quantity.SITE) {
//				Preconditions.checkState(key.sites.length == 1,
//						"Holding site constant, but VarGroupingKey has %s sites!", key.sites.length);
				if (key.sites.length == 1)
					constValues[i] = key.sites[0];
				else
					constValues[i] = Lists.newArrayList(key.sites);
			} else  if (constQuantities[i] == Quantity.DISTANCE) {
				constValues[i] = key.distance;
			} else {
				throw new IllegalStateException("Currently only set up for constant site or distance");
			}
		}
		
		if (key.restrictTo != null && !key.restrictTo.isEmpty()) {
			int prevNum = constQuantities.length;
			constQuantities = Arrays.copyOf(constQuantities, prevNum+key.restrictTo.size());
			constValues = Arrays.copyOf(constValues, prevNum+key.restrictTo.size());
			int index = prevNum;
			for (Quantity q : key.restrictTo.keySet()) {
				Object v = key.restrictTo.get(q);
				constQuantities[index] = q;
				constValues[index] = v;
				index++;
			}
		}
		
		PeriodDepResidualsList ret =  calcResiduals(key.magnitude, constQuantities, constValues,
				key.groupQuantities, key.type.singletonQuantities, key.realDataList);
		
		if (key.sites.length > 1 && !key.type.stdDevOfMedians && !key.type.separateSites) {
			// we have multiple sites and they're all bundled together, but this is a within-event type
			// variability, which means we should only have 1 per site in each residual set
		}
		
		return ret;
	}
	
	private VariabilityResult calcVarResult(VarGroupingKey key) throws IOException {
		boolean lightweight = key.realDataList != null;
		if (key.distance == null) {
			// combine across all distances
			VariabilityResult[] allVars = new VariabilityResult[distances.size()];
			for (int i=0; i<allVars.length; i++) {
				VarGroupingKey distKey = new VarGroupingKey(key.type,
						key.separateQuantities, key.groupQuantities, key.magnitude, distances.get(i), key.sites);
				distKey.restrictTo = key.restrictTo;
				distKey.realDataList = key.realDataList;
				try {
					if (key.realDataList == null)
						allVars[i] = varResultCache.get(distKey);
					else
						allVars[i] = calcVarResult(distKey);
				} catch (ExecutionException e) {
					ExceptionUtils.throwAsRuntimeException(e);
				}
			}
			if (allVars.length == 1)
				// only one distances
				return allVars[0];
			ResidualStdDevSet[] residualStdDevs;
			MedianStdDevSet[] medianStdDevs;
			if (key.type.stdDevOfMedians) {
				residualStdDevs = null;
				medianStdDevs = new MedianStdDevSet[calcPeriods.length];
				for (int p=0; p<calcPeriods.length; p++) {
					MedianStdDevSet[] allMedStdDevs = new MedianStdDevSet[allVars.length];
					for (int i=0; i<allVars.length; i++) {
						allMedStdDevs[i] = allVars[i].medianStdDevs[p];
					}
					medianStdDevs[p] = new MedianStdDevSet(allMedStdDevs);
				}
			} else {
				residualStdDevs = new ResidualStdDevSet[calcPeriods.length];
				medianStdDevs = null;
				for (int p=0; p<calcPeriods.length; p++) {
					ResidualStdDevSet[] allResStdDevs = new ResidualStdDevSet[allVars.length];
					for (int i=0; i<allVars.length; i++) {
						allResStdDevs[i] = allVars[i].residualStdDevs[p];
					}
					residualStdDevs[p] = new ResidualStdDevSet(allResStdDevs);
				}
			}
			return new VariabilityResult(null, null, residualStdDevs, medianStdDevs, calcPeriods, hasMagDist);
		}
		if (key.separateSites && key.sites.length > 1) {
			// do all sites at once
			List<List<ResidualSet>> totalPeriodResiduals = new ArrayList<>();
			for (int p=0; p<calcPeriods.length; p++)
				totalPeriodResiduals.add(new ArrayList<>());
			
			List<MedianStdDevSet[]> sitePeriodMediansList = new ArrayList<>();
			for (int p=0; p<calcPeriods.length; p++)
				sitePeriodMediansList.add(new MedianStdDevSet[key.sites.length]);
			
			for (int s=0; s<key.sites.length; s++) {
				Site site = key.sites[s];
				VarGroupingKey siteKey = new VarGroupingKey(key.type,
						key.separateQuantities, key.groupQuantities, key.magnitude, key.distance, site);
				siteKey.restrictTo = key.restrictTo;
				siteKey.realDataList = key.realDataList;
				
				if (D) System.out.println("Loading residuals for site "+site.getName()+", dist="+key.distance);
				List<List<ResidualSet>> residuals = loadResiduals(siteKey);

				LogValueSet[] logVals = lightweight ? null : new LogValueSet[calcPeriods.length];
				ResidualStdDevSet[] residualStdDevs;
				MedianStdDevSet[] medianStdDevs;
				RotationSpec[] commonRotSpecs = null;
				if (key.type.stdDevOfMedians) {
					residualStdDevs = null;
					medianStdDevs = new MedianStdDevSet[calcPeriods.length];
					for (int p=0; p<calcPeriods.length; p++) {
						List<ResidualSet> periodResiduals = residuals.get(p);
						if (!lightweight)
							logVals[p] = new LogValueSet(periodResiduals);
						medianStdDevs[p] = new MedianStdDevSet(lightweight, periodResiduals);
						if (p == 0) {
							commonRotSpecs = new RotationSpec[periodResiduals.size()];
							for (int i=0; i<periodResiduals.size(); i++)
								commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
						}
						totalPeriodResiduals.get(p).addAll(residuals.get(p));
						sitePeriodMediansList.get(p)[s] = medianStdDevs[p];
					}
				} else {
					residualStdDevs = new ResidualStdDevSet[calcPeriods.length];
					medianStdDevs = null;
					for (int p=0; p<calcPeriods.length; p++) {
						List<ResidualSet> periodResiduals = residuals.get(p);
						if (!lightweight)
							logVals[p] = new LogValueSet(periodResiduals);
						residualStdDevs[p] = new ResidualStdDevSet(lightweight, periodResiduals);
						if (p == 0) {
							commonRotSpecs = new RotationSpec[periodResiduals.size()];
							for (int i=0; i<periodResiduals.size(); i++)
								commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
						}
						totalPeriodResiduals.get(p).addAll(residuals.get(p));
					}
				}
				
				if (key.realDataList == null)
					varResultCache.put(siteKey, new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs,
							calcPeriods, hasMagDist));
			}

			// now combine sites
			if (D) System.out.println("Combining residuals for "+key.sites.length+" sites, dist="+key.distance);
			LogValueSet[] totLogVals = lightweight ? null : new LogValueSet[calcPeriods.length];
			ResidualStdDevSet[] totResidualStdDevs;
			MedianStdDevSet[] totMedianStdDevs;
			RotationSpec[] commonRotSpecs = null;
			if (key.type.stdDevOfMedians) {
				totResidualStdDevs = null;
				totMedianStdDevs = new MedianStdDevSet[calcPeriods.length];
				for (int p=0; p<calcPeriods.length; p++) {
					List<ResidualSet> periodResiduals = totalPeriodResiduals.get(p);
					if (!lightweight)
						totLogVals[p] = new LogValueSet(periodResiduals);
					if (p == 0) {
						commonRotSpecs = new RotationSpec[periodResiduals.size()];
						for (int i=0; i<periodResiduals.size(); i++)
							commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
					}
					totMedianStdDevs[p] = new MedianStdDevSet(sitePeriodMediansList.get(p));
				}
			} else {
				totResidualStdDevs = new ResidualStdDevSet[calcPeriods.length];
				totMedianStdDevs = null;
				for (int p=0; p<calcPeriods.length; p++) {
					List<ResidualSet> periodResiduals = totalPeriodResiduals.get(p);
					if (!lightweight)
						totLogVals[p] = new LogValueSet(periodResiduals);
					totResidualStdDevs[p] = new ResidualStdDevSet(lightweight, periodResiduals);
					if (p == 0) {
						commonRotSpecs = new RotationSpec[periodResiduals.size()];
						for (int i=0; i<periodResiduals.size(); i++)
							commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
					}
				}
			}
			
			return new VariabilityResult(commonRotSpecs, totLogVals, totResidualStdDevs, totMedianStdDevs, calcPeriods, hasMagDist);
		}
		List<List<ResidualSet>> residuals = loadResiduals(key);

//		System.out.println("have "+key.sites.length+" sites, "+residuals.get(0).size()+" residuals, separate="+key.separateSites);
		if (D) System.out.print("Combining residuals for "
				+(key.sites.length == 1 ? key.sites[0].getName() : key.sites.length+" sites")
				+", mag="+key.magnitude+", dist="+key.distance+"...");
		LogValueSet[] logVals = lightweight ? null : new LogValueSet[calcPeriods.length];
		ResidualStdDevSet[] residualStdDevs;
		MedianStdDevSet[] medianStdDevs;
		RotationSpec[] commonRotSpecs = null;
		if (key.type.stdDevOfMedians) {
			residualStdDevs = null;
			medianStdDevs = new MedianStdDevSet[calcPeriods.length];
			for (int p=0; p<calcPeriods.length; p++) {
				List<ResidualSet> periodResiduals = residuals.get(p);
				if (!lightweight)
					logVals[p] = new LogValueSet(periodResiduals);
				medianStdDevs[p] = new MedianStdDevSet(lightweight, periodResiduals);
				if (p == 0) {
					commonRotSpecs = new RotationSpec[periodResiduals.size()];
					for (int i=0; i<periodResiduals.size(); i++)
						commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
				}
			}
		} else {
			residualStdDevs = new ResidualStdDevSet[calcPeriods.length];
			medianStdDevs = null;
			for (int p=0; p<calcPeriods.length; p++) {
				List<ResidualSet> periodResiduals = residuals.get(p);
				if (!lightweight)
					logVals[p] = new LogValueSet(periodResiduals);
				residualStdDevs[p] = new ResidualStdDevSet(lightweight, periodResiduals);
				if (p == 0) {
					commonRotSpecs = new RotationSpec[periodResiduals.size()];
					for (int i=0; i<periodResiduals.size(); i++)
						commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
				}
			}
		}
		if (D) System.out.println("DONE");
		return new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs, calcPeriods, hasMagDist);
	}
	
	double getDeltaDistance(Float distance) {
		return distance == null || distance > 80 ? 20 : 10;
	}
	
	private List<List<ASK_EventData>> getRealDataListForKey(VarGroupingKey key) {
		double dm = 0.2;
		double dr = getDeltaDistance(key.distance);
		Map<Integer, List<ASK_EventData>> dataMatches = ASK_EventData.getMatches(realEventData,
				com.google.common.collect.Range.closed(key.magnitude-dm, key.magnitude+dm),
				key.distance == null ? null : com.google.common.collect.Range.closed(key.distance-dr, key.distance+dr), null, 0d);
		Preconditions.checkNotNull(dataMatches);
		// don't want to do that actually, dist-indep is dist-indep, so use all distances
//		if (key.distance == null && key.type.separateDists && distances.size() > 1) {
//			// we're calculating across all of our distances, to compare to distance independent recordings
//			// instead partition the number of recordings equally across each distance
//			Map<Integer, List<ASK_EventData>> filteredMatches = new HashMap<>();
//			for (int id : dataMatches.keySet()) {
//				List<ASK_EventData> recordings = dataMatches.get(id);
//				int targetNum = (int)Math.round((double)recordings.size()/(double)distances.size());
//				if (targetNum < 2 && recordings.size() >= 2)
//					targetNum = 2;
//				if (D) System.out.println("Target num="+targetNum+" for "+recordings.size()+" recordings and "+distances.size()+" dists");
//				if (targetNum > 1)
//					filteredMatches.put(id, recordings.subList(0, targetNum));
//			}
//			dataMatches = filteredMatches;
//		}
		List<List<ASK_EventData>> dataList = new ArrayList<>();
		for (int id : dataMatches.keySet())
			dataList.add(dataMatches.get(id));
		return dataList;
	}
	
	private VariabilityResult[] calcDownsampledVarResults(VarGroupingKey key, int numRealizations) throws IOException {
		List<List<ASK_EventData>> dataList = getRealDataListForKey(key);
		
		if (dataList.size() < 2)
			return new VariabilityResult[0];
		
		VariabilityResult[] ret = new VariabilityResult[numRealizations];
		
		VarGroupingKey randKey = new VarGroupingKey(key.type, key.separateQuantities, key.groupQuantities, key.magnitude, key.distance, key.sites);
		randKey.realDataList = dataList;
		
		for (int i=0; i<numRealizations; i++)
			// each call will be a new stochastic realization due to custom equals/hashCode implementation with realDataList set
			ret[i] = calcVarResult(randKey);
		
		return ret;
	}
	
	private static Random r = new Random();
	
	private class VarGroupingKey {
		private final VariabilityType type;
		private final Quantity[] separateQuantities;
		private final Quantity[] groupQuantities;
		private final Double magnitude;
		private final Float distance;
		private final Site[] sites;
		private Map<Quantity, Object> restrictTo;
		private List<List<ASK_EventData>> realDataList;
		
		// derrivative quantities
		private boolean separateSites;
		
		public VarGroupingKey(VariabilityType type, Double magnitude, Float distance, List<Site> sites) {
			this(type, magnitude, distance, sites.toArray(new Site[0]));
		}
		
		public VarGroupingKey(VariabilityType type, Double magnitude, Float distance, Site... sites) {
			this(type, type.separateQuantities, type.groupQuantities, magnitude, distance, sites);
		}
		
		public VarGroupingKey(VariabilityType type, Quantity[] separateQuantities, Quantity[] groupQuantities, Double magnitude,
				Float distance, Site... sites) {
			this.type = type;
			this.separateQuantities = separateQuantities;
			this.groupQuantities = groupQuantities;
			this.magnitude = magnitude;
			this.distance = distance;
			if (sites == null || sites.length == 0)
				sites = RotatedRupVariabilityPageGen.this.sites.toArray(new Site[0]);
			this.sites = sites;
			
			separateSites = false;
			for (Quantity quantity : separateQuantities)
				if (quantity == Quantity.SITE)
					separateSites = true;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + (separateSites ? 1231 : 1237);
			result = prime * result + ((type == null) ? 0 : type.hashCode());
			result = prime * result + ((distance == null) ? 0 : distance.hashCode());
			result = prime * result + Arrays.hashCode(groupQuantities);
			result = prime * result + ((magnitude == null) ? 0 : magnitude.hashCode());
			result = prime * result + Arrays.hashCode(separateQuantities);
			result = prime * result + Arrays.hashCode(sites);
			result = prime * result + ((restrictTo == null) ? 0 : restrictTo.hashCode());
			// never match, always generate new set stochastically
			result = prime * result + ((realDataList == null) ? 0 : r.nextInt());
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
			VarGroupingKey other = (VarGroupingKey) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (separateSites != other.separateSites)
				return false;
			if (type != other.type)
				return false;
			if (distance == null) {
				if (other.distance != null)
					return false;
			} else if (!distance.equals(other.distance))
				return false;
			if (!Arrays.equals(groupQuantities, other.groupQuantities))
				return false;
			if (magnitude == null) {
				if (other.magnitude != null)
					return false;
			} else if (!magnitude.equals(other.magnitude))
				return false;
			if (!Arrays.equals(separateQuantities, other.separateQuantities))
				return false;
			if (!Arrays.equals(sites, other.sites))
				return false;
			if (restrictTo == null) {
				if (other.restrictTo != null)
					return false;
			} else if (!restrictTo.equals(other.restrictTo))
				return false;
			if (realDataList != null)
				// always generate a new one with restrict to
				return false;
			return true;
		}

		private RotatedRupVariabilityPageGen getOuterType() {
			return RotatedRupVariabilityPageGen.this;
		}

	}
	
	private static class VariabilityResult {
		private final RotationSpec[] commonRotationSpecs;
		private final LogValueSet[] logVals;
		private final ResidualStdDevSet[] residualStdDevs;
		private final ResidualStdDevSet periodIndepResidualStdDev;
		private final MedianStdDevSet[] medianStdDevs;
		private final MedianStdDevSet periodIndepMedianStdDev;
		
		private double[] calcPeriods;
		private boolean hasMagDist;
		
		public VariabilityResult(RotationSpec[] commonRotationSpecs, LogValueSet[] logVals,
				ResidualStdDevSet[] residualStdDevs, MedianStdDevSet[] medianStdDevs,
				double[] calcPeriods, boolean hasMagDist) {
			this.calcPeriods = calcPeriods;
			this.hasMagDist = hasMagDist;
			this.commonRotationSpecs = commonRotationSpecs;
//			if (hasMagDist)
				this.logVals = logVals;
//			else
//				// only need them for mag-dist
//				this.logVals = null;
			this.residualStdDevs = residualStdDevs;
			this.medianStdDevs = medianStdDevs;
			
			if (residualStdDevs == null)
				periodIndepResidualStdDev = null;
			else
				periodIndepResidualStdDev = new ResidualStdDevSet(residualStdDevs);
			if (medianStdDevs == null)
				periodIndepMedianStdDev = null;
			else
				periodIndepMedianStdDev = new MedianStdDevSet(medianStdDevs);
		}
		
		public LogValueSet getLogValues(double period) {
			return logVals[Doubles.indexOf(calcPeriods, period)];
		}
		
		public ResidualStdDevSet getResidualStdDevSet(double period) {
			return residualStdDevs[Doubles.indexOf(calcPeriods, period)];
		}
		
		public MedianStdDevSet getMedianStdDevSet(double period) {
			return medianStdDevs[Doubles.indexOf(calcPeriods, period)];
		}
		
		public ResidualStdDevSet getPeriodIndepResidualStdDevSet() {
			return periodIndepResidualStdDev;
		}
		
		public MedianStdDevSet getPeriodIndepMedianStdDevSet() {
			return periodIndepMedianStdDev;
		}
	}
	
	private static class LogValueSet {
		private final double mean;
		private final double median;
		private final double min;
		private final double max;
		private final double[] residuals;
		private final double[] rawVals;
		
		public LogValueSet(List<ResidualSet> residuals) {
			int totSize = 0;
			for (ResidualSet set : residuals)
				totSize += set.size();
			this.residuals = new double[totSize];
			rawVals = new double[totSize];
			int index = 0;
			for (ResidualSet set : residuals) {
				double[] myRawVals = set.getLogValues();
				System.arraycopy(myRawVals, 0, rawVals, index, myRawVals.length);
				System.arraycopy(set.getResiduals(), 0, this.residuals, index, myRawVals.length);
				index += myRawVals.length;
			}
			Preconditions.checkState(index == totSize);
			mean = StatUtils.mean(rawVals);
			median = DataUtils.median(rawVals);
			min = StatUtils.min(rawVals);
			max = StatUtils.max(rawVals);
		}
	}
	
	private static List<ResidualSet> getCombinedResiduals(ResidualStdDevSet... sets) {
		List<ResidualSet> allResiduals = new ArrayList<>();
		for (ResidualStdDevSet set : sets)
			allResiduals.addAll(set.residuals);
		return allResiduals;
	}
	
	private static class ResidualStdDevSet {
		private final double total;
		private final double mean;
		private final double median;
		private final double min;
		private final double max;
		private final double[] stdDevs;
		private List<ResidualSet> residuals;
		
		public ResidualStdDevSet(ResidualStdDevSet... sets) {
			this(sets[0].stdDevs == null, getCombinedResiduals(sets));
		}
		
		public ResidualStdDevSet(boolean lightweight, List<ResidualSet> residuals) {
			this.residuals = residuals;
			double[] stdDevs = new double[residuals.size()];
			
			// see if equal weight
			int totLen = 0;
			int sizeEach = 0;
			for (int i=0; i<stdDevs.length; i++) {
				int size = residuals.get(i).size();
				totLen += size;
				if (i == 0) {
					sizeEach = size;
				} else {
					if (size != sizeEach)
						sizeEach = -1;
				}
			}
			double[] residualsArray = new double[totLen];
			double[] weights = sizeEach > 0 ? null : new double[totLen];
			int index = 0;
			for (int i=0; i<stdDevs.length; i++) {
				ResidualSet set = residuals.get(i);
				stdDevs[i] = set.getStdDev();
				int myLen = set.residuals.length;
				System.arraycopy(set.residuals, 0, residualsArray, index, myLen);
				if (weights != null) {
					double weight = Math.sqrt(myLen);
					for (int j=index; j<index+myLen; j++)
						weights[j] = weight;
				}
				index += myLen;
			}
			if (weights == null) {
				total = Math.sqrt(new Variance().evaluate(residualsArray));
			} else {
				weights = MathArrays.normalizeArray(weights, weights.length);
				total = Math.sqrt(new Variance().evaluate(residualsArray, weights));
			}
			mean = StatUtils.mean(stdDevs);
			median = DataUtils.median(stdDevs);
			min = StatUtils.min(stdDevs);
			max = StatUtils.max(stdDevs);
			
			if (lightweight)
				this.stdDevs = null;
			else
				this.stdDevs = stdDevs;
		}
	}
	
	private static class MedianStdDevSet {
		private final double stdDev;
		private final double[] medians;
		private final int num;
		private final double meanMedian;
		private final double minMedian;
		private final double maxMedian;
		
		public MedianStdDevSet(MedianStdDevSet... otherMedians) {
			boolean lightweight = otherMedians[0].medians == null;
			double[] stdDevs = new double[otherMedians.length];
			int totSize = 0;
			for (int i=0; i<otherMedians.length; i++) {
				MedianStdDevSet o = otherMedians[i];
				totSize += o.num;
				stdDevs[i] = o.stdDev;
			}
			// mean of the std dev's for each site (or distance if across multiple distances)
			stdDev = StatUtils.mean(stdDevs);
			if (lightweight) {
				this.medians = null;
				double minMedian = Double.POSITIVE_INFINITY;
				double maxMedian = Double.NEGATIVE_INFINITY;
				double meanMedian = 0d;
				for (int i=0; i<otherMedians.length; i++) {
					minMedian = Math.min(minMedian, otherMedians[i].minMedian);
					maxMedian = Math.max(maxMedian, otherMedians[i].maxMedian);
					meanMedian += otherMedians[i].meanMedian*(double)otherMedians[i].num;
				}
				meanMedian /= (double)totSize;
				this.minMedian = minMedian;
				this.maxMedian = maxMedian;
				this.meanMedian = meanMedian;
			} else {
				double[] medians = new double[totSize];
				int ind = 0;
				for (int i=0; i<otherMedians.length; i++) {
					System.arraycopy(otherMedians[i].medians, 0, medians, ind, otherMedians[i].medians.length);
					ind += otherMedians[i].medians.length;
				}
				meanMedian = StatUtils.mean(medians);
				minMedian = StatUtils.min(medians);
				maxMedian = StatUtils.max(medians);
				this.medians = medians;
			}
			
			num = totSize;
		}
		
		public MedianStdDevSet(boolean lightweight, List<ResidualSet> residuals) {
			double[] medians = new double[residuals.size()];
			Preconditions.checkState(medians.length > 1, "Only 1 event term, can't compute std dev");
			double[] weights = new double[residuals.size()];
			boolean equalWeight = true;
			for (int i=0; i<medians.length; i++) {
				ResidualSet set = residuals.get(i);
				medians[i] = set.median;
				weights[i] = Math.sqrt(Math.sqrt(set.size()));
				equalWeight = equalWeight && (i == 0 || (float)weights[i] == (float)weights[0]);
			}
			if (equalWeight)
				stdDev = Math.sqrt(StatUtils.variance(medians));
			else
				stdDev = Math.sqrt(new Variance().evaluate(medians, MathArrays.normalizeArray(weights, weights.length)));
//			System.out.println("Calculated stdDev from "+medians.length+" medians: "+stdDev);
			meanMedian = StatUtils.mean(medians);
			minMedian = StatUtils.min(medians);
			maxMedian = StatUtils.max(medians);
			
			if (lightweight)
				this.medians = null;
			else
				this.medians = medians;
			this.num = medians.length;
		}
	}
	
	protected static class EventTermKey {
		private final int eventID;
		private final double mag;
		private final float distance;
		private final Site[] sites;
		public EventTermKey(int eventID, double mag, float distance, Site[] sites) {
			super();
			this.eventID = eventID;
			this.mag = mag;
			this.distance = distance;
			this.sites = sites;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Float.floatToIntBits(distance);
			result = prime * result + eventID;
			long temp;
			temp = Double.doubleToLongBits(mag);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			result = prime * result + Arrays.hashCode(sites);
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
			EventTermKey other = (EventTermKey) obj;
			if (Float.floatToIntBits(distance) != Float.floatToIntBits(other.distance))
				return false;
			if (eventID != other.eventID)
				return false;
			if (Double.doubleToLongBits(mag) != Double.doubleToLongBits(other.mag))
				return false;
			if (!Arrays.equals(sites, other.sites))
				return false;
			return true;
		}
	}
	
	protected static class EventTerm {
		final double[] eventTerms;
		final List<RotationSpec> rotations;
		final List<double[]> rotLogVals;
		
		public EventTerm(double[] eventTerms, List<RotationSpec> rotations, List<double[]> rotLogVals) {
			super();
			this.eventTerms = eventTerms;
			this.rotations = rotations;
			this.rotLogVals = rotLogVals;
		}
	}
	
	private EventTerm calcEventTerm(Integer eventID, Double mag, Float distance,
			List<Site> sites) throws IOException {
		return calcEventTerm(eventID, mag, distance, sites.toArray(new Site[0]));
	}
	
	private EventTerm calcEventTerm(Integer eventID, Double mag, Float distance,
			Site... sites) throws IOException {
		Preconditions.checkState(sites.length > 0);
		HashSet<RotationSpec> rotationSet = new HashSet<>();
		RotatedRupVariabilityConfig<E> config = magConfigs.get(mag);
		for (Site site : sites)
			rotationSet.addAll(config.getRotationsForQuantities(Quantity.SITE, site,
					Quantity.EVENT_ID, eventID, Quantity.DISTANCE, distance));
		Mean[] means = new Mean[calcPeriods.length];
		for (int i=0; i<means.length; i++)
			means[i] = new Mean();
		SimulationRotDProvider<RotationSpec> prov = magProvs.get(mag);
		List<RotationSpec> rotations = new ArrayList<>(rotationSet);
		
		List<double[]> rotLogVals = new ArrayList<>();
		for (RotationSpec spec : rotations) {
			DiscretizedFunc rd50s = prov.getRotD50(spec.site, spec, 0);
			double[] vals = new double[calcPeriods.length];
			for (int i=0; i<calcPeriods.length; i++) {
				vals[i] = Math.log(rd50s.getInterpolatedY(calcPeriods[i]));
				means[i].increment(vals[i]);
			}
			rotLogVals.add(vals);
		}
		
		double[] ret = new double[means.length];
		for (int i=0; i<ret.length; i++)
			ret[i] = means[i].getResult();
		
		return new EventTerm(ret, rotations, rotLogVals);
	}
	
	private Random repeatableRandom;
	
	private PeriodDepResidualsList calcResiduals(Double magnitude, Quantity[] constQuantities, Object[] constValues,
			Quantity[] groupQuantities, Quantity[] singletons, List<List<ASK_EventData>> realDataList)
					throws IOException {
//		if (D) System.out.println("Rotations for mag: "+magnitude);
		RotatedRupVariabilityConfig<E> config = magConfigs.get(magnitude);
		List<RotationSpec> totalRotations;
		if (constQuantities != null && constQuantities.length > 0) {
//			for (int i=0; i<constQuantities.length; i++)
//				System.out.println("Constant: "+constQuantities[i]+": "+constValues[i]);
			totalRotations = config.getRotationsForQuantities(constQuantities, constValues);
			// remove any constants from the group quantities list
			List<Quantity> groupList = Lists.newArrayList(groupQuantities);
			for (Quantity quantity : constQuantities)
				groupList.remove(quantity);
			groupQuantities = groupList.toArray(new Quantity[0]);
		} else {
			totalRotations = config.getRotations();
		}
		if (realDataList != null) {
			// now downsample randomly
			
			if (repeatableRandom == null)
				repeatableRandom = new Random(magConfigs.get(magnitude).getRotations().size());
			
			List<List<ASK_EventData>> filteredRealDataList = new ArrayList<>();
			for (List<ASK_EventData> datas : realDataList)
				if (datas.size() > 1)
					filteredRealDataList.add(datas);
			realDataList = filteredRealDataList;
			
			Map<Integer, List<RotationSpec>> eventRotMap = new HashMap<>();
			List<Integer> eventIDs = new ArrayList<>();
			for (RotationSpec rot : totalRotations) {
				int eventID = rot.eventID;
				List<RotationSpec> eventRots = eventRotMap.get(eventID);
				if (eventRots == null) {
					eventRots = new ArrayList<>();
					eventRotMap.put(eventID, eventRots);
					eventIDs.add(eventID);
				}
				eventRots.add(rot);
			}
			Collections.shuffle(eventIDs, repeatableRandom);
			List<RotationSpec> downsampledRots = new ArrayList<>();
			int origNumRecordings = 0;
			for (List<ASK_EventData> data : realDataList)
				origNumRecordings += data.size();
			Preconditions.checkState(eventIDs.size() > 1, "Don't have enough events: %s", eventIDs.size());
			int eventsToInclude = Integer.min(eventIDs.size(), realDataList.size());
			Preconditions.checkState(eventsToInclude > 1, "eventsToInclude must be > 1: %s", eventsToInclude);
//			System.out.println("Doing "+eventsToInclude+" events");
			for (int i=0; i<eventsToInclude; i++) {
				int eventID = eventIDs.get(i);
				List<RotationSpec> allEventRotations = new ArrayList<>(eventRotMap.get(eventID));
				List<RotationSpec> eventGroupedRotations;
				if (groupQuantities != null && groupQuantities.length > 0) {
					// we have grouped quantities, and need to randomly choose one group
					Object[] randGroupValues = new Object[groupQuantities.length];
					for (int g=0; g<groupQuantities.length; g++) {
						if (groupQuantities[g] == Quantity.EVENT_ID) {
							randGroupValues[g] = eventID;
						} else {
							List<?> allValues = config.getQuantitiesMap().get(groupQuantities[g]);
							// randomly pick one
							randGroupValues[g] = allValues.get(repeatableRandom.nextInt(allValues.size()));
//							System.out.println("Random choice for group "+groupQuantities[g]+": "+randGroupValues[g]);
						}
					}
					eventGroupedRotations = RotatedRupVariabilityConfig.getRotationsForQuantities(
							allEventRotations, groupQuantities, randGroupValues);
//					System.out.println("Have "+eventGroupedRotations.size()+" for random group");
				} else {
					eventGroupedRotations = allEventRotations;
					System.out.println("Have "+eventGroupedRotations.size()+" (no groups)");
				}
				int rawNum = realDataList.get(i).size();
//				System.out.println("\tEvent "+i+" has "+rawNum+" recordings");
				int numRecordings = Integer.min(eventGroupedRotations.size(), rawNum);
				if (numRecordings < 2) {
//					System.out.println("\tWARNING: skipping recorded event with only 1 recording");
					continue;
				}
//				if (numRecordings < realDataList.get(i).size())
//					System.out.println("\tWARNING: only using "+numRecordings+" rotations for event, "
//							+ "even though we have "+realDataList.get(i).size()+" recordings");
				if (singletons != null && singletons.length > 0) {
					for (Quantity singleton : singletons) {
						List<List<RotationSpec>> temp = processSingletons(eventGroupedRotations, singleton);
						eventGroupedRotations = temp.get(repeatableRandom.nextInt(temp.size()));
					}
				}
				if (numRecordings < eventGroupedRotations.size()) {
					Collections.shuffle(eventGroupedRotations, repeatableRandom);
					for (int j=0; j<numRecordings; j++)
						downsampledRots.add(eventGroupedRotations.get(j));
				} else {
					downsampledRots.addAll(eventGroupedRotations);
				}
//				System.out.println("\tActually used "+numRecordings+" recordings. Cumulative# : "+downsampledRots.size());
			}
			if (D) System.out.println("Downsampled from "+totalRotations.size()+" to "+downsampledRots.size()+" rotations to match data. Have "
					+realDataList.size()+" real events, "+eventIDs.size()+" simulated events. Originally had "+origNumRecordings+" recordings.");
			totalRotations = downsampledRots;
		}
		int num = totalRotations.size();
		PeriodDepResidualsList ret = new PeriodDepResidualsList(totalRotations);
		for (int i=0; i<calcPeriods.length; i++)
			ret.add(new ArrayList<>());
		
		if (D) System.out.println("Computing residuals/medians for Constants ["
				+(constQuantities == null ? "" : Joiner.on(",").join(constQuantities))
				+"], Groups ["+Joiner.on(",").join(groupQuantities)+"]: "+num+" simulations");
		
//		System.out.println("Computing residuals from "+num+" simulations");
		
		calcGroupedResiduals(magnitude, totalRotations, calcPeriods, ret, groupQuantities, singletons);
		int count = 0;
		for (ResidualSet set : ret.get(new Random().nextInt(calcPeriods.length)))
			count += set.size();
		Preconditions.checkState((singletons != null && singletons.length > 0) || count == num,
				"Bad end index (and no singletons). Expected %s, have %s", num, count);
//		System.out.println(numThatVary+"/"+totNum+" have path variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
		
		return ret;
	}
	
	private void calcGroupedResiduals(double magnitude, List<RotationSpec> rotations, double[] periods,
			List<List<ResidualSet>> ret, Quantity[] groupQuantities, Quantity[] singletons) throws IOException {
//		Preconditions.checkState(!rotations.isEmpty()); // can be empty for downsampling
		if (rotations.isEmpty())
			return;
		if (groupQuantities.length == 0) {
			List<List<RotationSpec>> rotsLists = new ArrayList<>();
			rotsLists.add(rotations);
			if (singletons != null && singletons.length > 0) {
				// process singletons
				for (Quantity singleton : singletons) {
					List<List<RotationSpec>> newRotsLists = new ArrayList<>();
					for (List<RotationSpec> rots : rotsLists)
						newRotsLists.addAll(processSingletons(rots, singleton));
					rotsLists = newRotsLists;
				}
			}
//			System.out.println("Computing");
			for (List<RotationSpec> rotsList : rotsLists) {
				// we're at a set for which to calculate residuals, do it
				List<DiscretizedFunc> spectra = new ArrayList<>();
				SimulationRotDProvider<RotationSpec> prov = magProvs.get(magnitude);
				for (RotationSpec rotation : rotsList)
					spectra.add(prov.getRotD50(rotation.site, rotation, 0));
				for (int p=0; p<periods.length; p++) {
					double[] values = new double[spectra.size()];
					for (int i=0; i<values.length; i++)
						values[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
					ret.get(p).add(new ResidualSet(rotsList, values));
				}
			}
//			System.out.println("Done computing");
		} else {
//			System.out.println("have "+rotations.size()+" rotations before "+groupQuantities[0]);
			Quantity[] current = { groupQuantities[0] };
			Quantity[] downstream = Arrays.copyOfRange(groupQuantities, 1, groupQuantities.length);
			for (Object value : magQuantitiesTable.get(magnitude, groupQuantities[0])) {
				List<RotationSpec> myRotations = RotatedRupVariabilityConfig.getRotationsForQuantities(
						rotations, current, new Object[] {value});
//				System.out.println("\t"+value+": "+myRotations.size());
				calcGroupedResiduals(magnitude, myRotations, periods, ret, downstream, singletons);
			}
		}
	}
	
	private List<List<RotationSpec>> processSingletons(List<RotationSpec> rotations, Quantity singleton) {
		Map<Object, List<RotationSpec>> valueMap = new HashMap<>();
		for (RotationSpec rot : rotations) {
			Object value = rot.getValue(singleton);
			List<RotationSpec> valRots = valueMap.get(value);
			if (valRots == null) {
				valRots = new LinkedList<>();
				valueMap.put(value, valRots);
			}
			valRots.add(rot);
		}
		List<List<RotationSpec>> ret = new ArrayList<>();
		if (valueMap.size() == 1) {
			// only one unique value
			ret.add(rotations);
			return ret;
		}
		
		int minNum = Integer.MAX_VALUE;
		int maxNum = 0;
		for (Object value : valueMap.keySet()) {
			List<RotationSpec> list = valueMap.get(value);
			minNum = Integer.min(list.size(), minNum);
			maxNum = Integer.max(list.size(), maxNum);
			Collections.shuffle(list);
		}
		
		if (minNum != maxNum)
			System.err.println("WARNING: Processing singleton for "+singleton+" and have uneven list sizes, "
					+ "will skip some rotations. Have at least "+minNum+" rotations for each unique value, "
					+ "but the largest value list has "+maxNum);
		Preconditions.checkState(minNum > 0);
		
		for (int i=0; i<minNum; i++) {
			List<RotationSpec> list = new ArrayList<>();
			for (Object value : valueMap.keySet())
				list.add(valueMap.get(value).get(i));
			ret.add(list);
		}
		return ret;
	}
	
	private class GMPE_GroupingKey {
		private final NGAW2_WrapperFullParam gmpe;
		private final Site site;
		private final double magnitude;
		private final double distance;
		private Map<Quantity, Object> restrictTo;
		
		public GMPE_GroupingKey(NGAW2_WrapperFullParam gmpe, Site site, double magnitude, double distance) {
			this.gmpe = gmpe;
			this.site = site;
			this.magnitude = magnitude;
			this.distance = distance;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			long temp;
			temp = Double.doubleToLongBits(distance);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			result = prime * result + ((gmpe == null) ? 0 : gmpe.hashCode());
			temp = Double.doubleToLongBits(magnitude);
			result = prime * result + ((site == null) ? 0 : site.getLocation().hashCode());
			result = prime * result + ((restrictTo == null) ? 0 : restrictTo.hashCode());
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
			GMPE_GroupingKey other = (GMPE_GroupingKey) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (Double.doubleToLongBits(distance) != Double.doubleToLongBits(other.distance))
				return false;
			if (gmpe == null) {
				if (other.gmpe != null)
					return false;
			} else if (!gmpe.equals(other.gmpe))
				return false;
			if (Double.doubleToLongBits(magnitude) != Double.doubleToLongBits(other.magnitude))
				return false;
			if (site == null) {
				if (other.site != null)
					return false;
			} else if (!site.equals(other.site))
				return false;
			if (restrictTo == null) {
				if (other.restrictTo != null)
					return false;
			} else if (!restrictTo.equals(other.restrictTo))
				return false;
			return true;
		}

		private RotatedRupVariabilityPageGen getOuterType() {
			return RotatedRupVariabilityPageGen.this;
		}
	}
	
	private static GMPE_Result packGMs(GMPE_Result... results) {
		if (results.length == 1)
			return results[0];
		int lenEach = -1;
		for (GMPE_Result result : results) {
			if (lenEach < 0)
				lenEach = result.gms[0].length;
			else
				Preconditions.checkState(lenEach == result.gms[0].length);
		}
		ScalarGroundMotion[][] gms = new ScalarGroundMotion[results[0].gms.length][lenEach*results.length];
		for (int i=0; i<results.length; i++)
			for (int p=0; p<gms.length; p++)
				System.arraycopy(results[i].gms[p], 0, gms[p], i*lenEach, lenEach);
		return new GMPE_Result(gms);
	}
	
	private static class GMPE_Result {
		final ScalarGroundMotion[][] gms;
		final double[] logMedian;
		final Map<GMPE_Var, double[]> medianVars;
		
		public GMPE_Result(ScalarGroundMotion[][] gms) {
			this.gms = gms;
			logMedian = new double[gms.length];
			medianVars = new HashMap<>();
			for (GMPE_Var var : GMPE_Var.values())
				medianVars.put(var, new double[gms.length]);
			for (int p=0; p<gms.length; p++) {
				double[] medians = new double[gms[p].length];
				Map<GMPE_Var, double[]> varVals = new HashMap<>();
				for (GMPE_Var var : GMPE_Var.values())
					varVals.put(var, new double[gms[p].length]);
				for (int i=0; i<gms[p].length; i++) {
					ScalarGroundMotion gm = gms[p][i];
					medians[i] = gm.mean();
					for (GMPE_Var var : GMPE_Var.values())
						varVals.get(var)[i] = var.calculate(gm);
				}
				logMedian[p] = DataUtils.median(medians);
				for (GMPE_Var var : GMPE_Var.values())
					medianVars.get(var)[p] = DataUtils.median(varVals.get(var));
			}
		}
		
		public double[] getVariability(GMPE_Var type) {
			return medianVars.get(type);
		}
	}
	
	private GMPE_Result calcGMPE(GMPE_GroupingKey key) throws ExecutionException {
		List<Integer> eventIDs = magEventIDs.size() == 1 ?
				magEventIDs.values().iterator().next() : magEventIDs.get(key.magnitude);
		if (sites.size() > 1 && key.site == null) {
			// do for all sites
			GMPE_Result[] results = new GMPE_Result[sites.size()];
			for (int i=0; i<sites.size(); i++) {
				Site site = sites.get(i);
				GMPE_GroupingKey newKey = new GMPE_GroupingKey(key.gmpe, site, key.magnitude, key.distance);
				newKey.restrictTo = key.restrictTo;
				results[i] = gmpeResultCache.get(newKey);
			}
			return packGMs(results);
		}
		Site site = key.site;
		if (site == null) {
			Preconditions.checkState(sites.size() == 1);
			site = sites.get(0);
		}
		NGAW2_WrapperFullParam gmpe = key.gmpe;
		gmpe.setSite(site);
		
		DistanceJBParameter rJBParam = (DistanceJBParameter) gmpe.getParameter(DistanceJBParameter.NAME);
		DistanceRupParameter rRupParam = (DistanceRupParameter) gmpe.getParameter(DistanceRupParameter.NAME);
		DistanceX_Parameter rXParam = (DistanceX_Parameter) gmpe.getParameter(DistanceX_Parameter.NAME);
		
		Float restrictSourceAz = key.restrictTo == null ? null : (Float)key.restrictTo.get(Quantity.SOURCE_AZIMUTH);
		int numSourceAz = restrictSourceAz == null ? sourceAzimuths.size() : 1;
		ScalarGroundMotion[][] result = new ScalarGroundMotion[calcPeriods.length][eventIDs.size()*numSourceAz];
		
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param saParam = (SA_Param) gmpe.getIntensityMeasure();
		int ind = 0;
		
		for (int eventID : eventIDs) {
			E event = getEvent(eventID);
			EqkRupture rup = gmpeEventCache.get(event);
			gmpe.setEqkRupture(rup);
			
			double zTOR = rup.getRuptureSurface().getAveRupTopDepth();
			double rRup, rJB;
			if (BBP_PartBValidationConfig.DIST_JB) {
				rJB = key.distance;
				rRup = zTOR == 0d ? rJB : Math.sqrt(zTOR*zTOR + rJB*rJB);
			} else {
				rRup = key.distance;
				rJB = zTOR == 0d ? rRup : Math.sqrt(rRup*rRup - zTOR*zTOR);
			}
			
			// override distances
			rJBParam.setValueIgnoreWarning(rJB);
			rRupParam.setValueIgnoreWarning(rRup);
			
			for (int i=0; i<numSourceAz; i++) {
				double sourceAz = restrictSourceAz == null ? sourceAzimuths.get(i).doubleValue() : restrictSourceAz.doubleValue();
				double rX = rJB * Math.sin(Math.toRadians(sourceAz));
				
				rXParam.setValueIgnoreWarning(rX);
				
				for (int p=0; p<calcPeriods.length; p++) {
					SA_Param.setPeriodInSA_Param(saParam, calcPeriods[p]);
					result[p][ind] = gmpe.getGroundMotion();
					if (eventID == 165223 && p == 0 && D)
						System.out.println("GMPE\t"+eventID+"\tsAz="+(float)sourceAz+"\trstAz="+restrictSourceAz
								+"\t\trRup="+(float)rRup+"\trJB="+(float)rJB+"\trX="+(float)rX+"\tgm="+gmpe.getGroundMotion().mean());
				}
				ind++;
			}
		}
		return new GMPE_Result(result);
	}
	
	private class StdDevPercentileFuncs {
		private final UncertainArbDiscDataset bounds95;
		private final UncertainArbDiscDataset bounds68;
		private final DiscretizedFunc sigmaFunc;
		
		public StdDevPercentileFuncs(UncertainArbDiscDataset bounds95, UncertainArbDiscDataset bounds68, DiscretizedFunc sigmaFunc) {
			this.bounds95 = bounds95;
			this.bounds68 = bounds68;
			this.sigmaFunc = sigmaFunc;
		}
	}
	
	private void plotPeriodDependentStdDevs(File resourcesDir, String prefix, String title,
			List<DiscretizedFunc> stdDevFuncs, List<StdDevPercentileFuncs> stdDevPercentiles,
			List<PlotCurveCharacterstics> bundleChars, DiscretizedFunc gmpeFunc) throws IOException {
		Preconditions.checkState(stdDevFuncs.size() == bundleChars.size());
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double minPeriod = StatUtils.min(calcPeriods);
		double maxPeriod = StatUtils.max(calcPeriods);
		
		if (stdDevPercentiles != null && stdDevPercentiles.get(0) != null) {
			// only plot the first one
			Color color95 = new Color(0, 0, 0, 40);
			Color color68 = new Color(0, 0, 0, 60);
			
//			DiscretizedFunc indepUpper95 = new ArbitrarilyDiscretizedFunc();
//			indepUpper95.set(minPeriod, totalStdDevPercentiles.bounds95.getUpperY(0d));
//			indepUpper95.set(maxPeriod, totalStdDevPercentiles.bounds95.getUpperY(0d));
//			funcs.add(indepUpper95);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color95));
//			
//			DiscretizedFunc indepLower95 = new ArbitrarilyDiscretizedFunc();
//			indepLower95.set(minPeriod, totalStdDevPercentiles.bounds95.getLowerY(0d));
//			indepLower95.set(maxPeriod, totalStdDevPercentiles.bounds95.getLowerY(0d));
//			funcs.add(indepLower95);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color95));
//			
//			DiscretizedFunc indepUpper68 = new ArbitrarilyDiscretizedFunc();
//			indepUpper68.set(minPeriod, totalStdDevPercentiles.bounds68.getUpperY(0d));
//			indepUpper68.set(maxPeriod, totalStdDevPercentiles.bounds68.getUpperY(0d));
//			funcs.add(indepUpper68);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color68));
//			
//			DiscretizedFunc indepLower68 = new ArbitrarilyDiscretizedFunc();
//			indepLower68.set(minPeriod, totalStdDevPercentiles.bounds68.getLowerY(0d));
//			indepLower68.set(maxPeriod, totalStdDevPercentiles.bounds68.getLowerY(0d));
//			funcs.add(indepLower68);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color68));
			
			color95 = new Color(0, 0, 0, 20);
			color68 = new Color(0, 0, 0, 40);
			
			StdDevPercentileFuncs totalStdDevPercentiles = stdDevPercentiles.get(0);
			
			totalStdDevPercentiles.bounds95.setName(null);
			funcs.add(totalStdDevPercentiles.bounds95);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, color95));
			totalStdDevPercentiles.bounds68.setName(null);
			funcs.add(totalStdDevPercentiles.bounds68);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, color68));
		}
		
		for (int i=stdDevFuncs.size(); --i>=0;) {
			DiscretizedFunc stdDevFunc = stdDevFuncs.get(i);
			funcs.add(stdDevFunc);
			PlotCurveCharacterstics bChar = bundleChars.get(i);
			chars.add(bChar);
			if (bChar.getLineWidth() >= 2f) {
				double distPeriodIndepVal = stdDevFunc.getY(0d);
				DiscretizedFunc periodIndepFunc = new ArbitrarilyDiscretizedFunc(stdDevFunc.getName()+" Period-Indep.");
				periodIndepFunc.set(minPeriod, distPeriodIndepVal);
				periodIndepFunc.set(maxPeriod, distPeriodIndepVal);
				
				funcs.add(periodIndepFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, bChar.getColor()));
			}
		}
		
		if (gmpeFunc != null) {
			funcs.add(gmpeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, gmpeColor));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Period (s)", "Standard Deviation");
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		double minX = StatUtils.min(calcPeriods);
		double maxX = StatUtils.max(calcPeriods);
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(0d, 1d);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
	}
	
	private File plotMultiDistPeriodDependentStdDevs(File resourcesDir, String prefix,
			VariabilityType type, Double mag, List<Float> dists) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double[] distArray = Doubles.toArray(dists);
		CPT distCPT = new CPT(StatUtils.min(distArray), StatUtils.max(distArray), Color.BLUE, Color.RED);
		
		for (Float dist : dists) {
			VarGroupingKey key = new VarGroupingKey(type, mag, dist, sites);
			VariabilityResult result;
			try {
				result = varResultCache.get(key);
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			
			DiscretizedFunc spectrum = new ArbitrarilyDiscretizedFunc(optionalDigitDF.format(dist)+" km");
			for (double period : calcPeriods) {
				if (type.stdDevOfMedians)
					spectrum.set(period, result.getMedianStdDevSet(period).stdDev);
				else
					spectrum.set(period, result.getResidualStdDevSet(period).total);
			}
			
			funcs.add(spectrum);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, distCPT.getColor(dist)));
		}
		
		// now distance independent
		VarGroupingKey key = new VarGroupingKey(type, mag, null, sites);
		VariabilityResult result;
		try {
			result = varResultCache.get(key);
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		DiscretizedFunc spectrum = new ArbitrarilyDiscretizedFunc("Dist.-Indep.");
		for (double period : calcPeriods) {
			if (type.stdDevOfMedians)
				spectrum.set(period, result.getMedianStdDevSet(period).stdDev);
			else
				spectrum.set(period, result.getResidualStdDevSet(period).total);
		}
		
		funcs.add(spectrum);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		double distPeriodIndepVal = type.stdDevOfMedians ? result.getPeriodIndepMedianStdDevSet().stdDev :
			result.getPeriodIndepResidualStdDevSet().total;
		DiscretizedFunc distPeriodIndepFunc = new ArbitrarilyDiscretizedFunc("Dist.- & Period-Indep.");
		distPeriodIndepFunc.set(StatUtils.min(calcPeriods), distPeriodIndepVal);
		distPeriodIndepFunc.set(StatUtils.max(calcPeriods), distPeriodIndepVal);
		
		funcs.add(distPeriodIndepFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLACK));
		
		if (realEventData != null) {
			List<VariabilityResult> downsampledResults = new ArrayList<>();
			try {
				if (type.separateSites && sites.size() > 1) {
					for (Site oSite : sites)
						for (VariabilityResult r : downsampledVarResultCache.get(new VarGroupingKey(type, mag, null, oSite)))
							downsampledResults.add(r);
				} else {
					for (VariabilityResult r : downsampledVarResultCache.get(key))
						downsampledResults.add(r);
				}
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (!downsampledResults.isEmpty()) {
				DiscretizedFunc lower95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper95Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc lower68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upper68Func = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc medianFunc = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc sigmaFunc = new ArbitrarilyDiscretizedFunc();
				List<Double> myCalcPeriods = new ArrayList<>();
				myCalcPeriods.add(0d);
				for (double period : calcPeriods)
					myCalcPeriods.add(period);
				for (double period : myCalcPeriods) {
					double[] values = new double[downsampledResults.size()];
					for (int j=0; j<values.length; j++) {
						if (type.stdDevOfMedians) {
							MedianStdDevSet set = period == 0d ? downsampledResults.get(j).getPeriodIndepMedianStdDevSet() :
								downsampledResults.get(j).getMedianStdDevSet(period);
							Preconditions.checkState(set.medians == null || set.medians.length > 1);
							values[j] = set.stdDev;
						} else {
							ResidualStdDevSet set = period == 0d ? downsampledResults.get(j).getPeriodIndepResidualStdDevSet() :
								downsampledResults.get(j).getResidualStdDevSet(period);
							values[j] = set.total;
						}
//						if (period == 1d) System.out.println("1s "+j+": "+values[j]);
					}
					lower95Func.set(period, StatUtils.percentile(values, 2.5));
					upper95Func.set(period, StatUtils.percentile(values, 97.5));
					lower68Func.set(period, StatUtils.percentile(values, 16));
					upper68Func.set(period, StatUtils.percentile(values, 84));
					medianFunc.set(period, StatUtils.percentile(values, 50d));
					sigmaFunc.set(period, Math.sqrt(StatUtils.variance(values)));
				}
//				for (double period : calcPeriods)
//					System.out.println((float)period+"\t"+optionalDigitDF.format(stdDevFunc.getY(period))+"\t["
//							+optionalDigitDF.format(lowerFunc.getY(period))+" "+optionalDigitDF.format(upperFunc.getY(period))+"]");
				
				UncertainArbDiscDataset func95 = new UncertainArbDiscDataset(medianFunc, lower95Func, upper95Func);
				UncertainArbDiscDataset func68 = new UncertainArbDiscDataset(medianFunc, lower68Func, upper68Func);
				
				func95.setName(null);
				funcs.add(0, func95);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
						new Color(0, 0, 0, 20)));
				func68.setName(null);
				funcs.add(1, func68);
				chars.add(1, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
						new Color(0, 0, 0, 40)));
			}
		}
		
		String title = mag == null ? "" : "M"+optionalDigitDF.format(mag)+" ";
		title += type.name+", "+type.symbol;
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Period (s)", type.symbol);
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setLegendFontSize(24);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		Range xRange = new Range(StatUtils.min(calcPeriods), StatUtils.max(calcPeriods));
		Range yRange = new Range(0d, 1d);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}
	
	private void plotStdDevsHistogram(File resourcesDir, String prefix, String title, double period,
			ResidualStdDevSet totalStdDev, List<ResidualStdDevSet> siteStdDevs) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Preconditions.checkState(siteStdDevs == null || siteStdDevs.size() == siteColors.size());
		
		double delta = 0.05;
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 1d, delta);
		hist.setName("Total Histogram");
		
//		for (ResidualSet set : residuals)
//			hist.add(hist.getClosestXIndex(set.residualStdDev), 1d);
		for (double stdDev : totalStdDev.stdDevs)
			hist.add(hist.getClosestXIndex(stdDev), 1d);
		hist.normalizeBySumOfY_Vals();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		
		if (siteStdDevs != null) {
			for (int i=0; siteStdDevs != null && i<siteStdDevs.size(); i++) {
				ResidualStdDevSet siteStdDev = siteStdDevs.get(i);
				funcs.add(line(siteStdDev.total, 0d, siteStdDev.total, 1d, sites.get(i).getName()));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, siteColors.get(i)));
			}
		}
		
		funcs.add(line(totalStdDev.total, 0d, totalStdDev.total, 1d, "Total"));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Standard Deviation", "");
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		Range xRange = new Range(0d, 1d);
		Range yRange = new Range(0d, 1d);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 450);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
	}
	
	private File[] plotSiteDownsampledHistogram(File resourcesDir, String prefix,
			VariabilityType type, Double magnitude, Float distance, Double[] periods,
			Site... sites) throws IOException {
		
		VarGroupingKey key = new VarGroupingKey(type, magnitude, distance, sites);
		VariabilityResult[] dsResults;
		VariabilityResult totResult;
		try {
			dsResults = downsampledVarResultCache.get(key);
			totResult = varResultCache.get(key);
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		double delta = 0.05;

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		File[] ret = new File[periods.length];
		for (int i=0; i<periods.length; i++) {
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 1d, delta);
			hist.setName("Downsampled Realizations");
			
			Double period = periods[i];
			
			double meanDownsampled = 0d;
			for (VariabilityResult result : dsResults) {
				double stdDev;
				if (type.stdDevOfMedians) {
					stdDev = period == null ? result.getPeriodIndepMedianStdDevSet().stdDev
							: result.getMedianStdDevSet(period).stdDev;
				} else {
					stdDev = period == null ? result.getPeriodIndepResidualStdDevSet().total
							: result.getResidualStdDevSet(period).total;
				}
				meanDownsampled += stdDev;
				hist.add(hist.getClosestXIndex(stdDev), 1d);
			}
			hist.normalizeBySumOfY_Vals();
			meanDownsampled /= (double)dsResults.length;
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			funcs.add(line(meanDownsampled, 0d, meanDownsampled, 1d, "Mean of Downsampled"));
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
			
			double totStdDev;
			if (type.stdDevOfMedians) {
				totStdDev = period == null ? totResult.getPeriodIndepMedianStdDevSet().stdDev
						: totResult.getMedianStdDevSet(period).stdDev;
			} else {
				totStdDev = period == null ? totResult.getPeriodIndepResidualStdDevSet().total
						: totResult.getResidualStdDevSet(period).total;
			}
			funcs.add(line(totStdDev, 0d, totStdDev, 1d, "Total"));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, type.symbol+" Downsampled Distribution", "Standard Deviation", "");
			spec.setLegendVisible(true);
			
			Range xRange = new Range(0d, 1d);
			Range yRange = new Range(0d, 1d);

			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(800, 450);
			
			String myPrefix = prefix;
			if (period == null)
				myPrefix += "_period_indep";
			else
				myPrefix += "_"+optionalDigitDF.format(period)+"s";
			File pngFile = new File(resourcesDir, myPrefix+".png");
			gp.saveAsPNG(pngFile.getAbsolutePath());
			ret[i] = pngFile;
		}
		return ret;
	}
	
	private XY_DataSet line(double x0, double y0, double x1, double y1, String name) {
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.setName(name);
		line.set(x0, y0);
		line.set(x1, y1);
		return line;
	}
	
	private File[][] plotScatter(File resourcesDir, String prefix, double[] periods, VariabilityType type,
			VariabilityResult result, ScatterDisaggQuantity scatterQuantity) throws IOException {
		if (D) System.out.println("Scatter will have "+result.commonRotationSpecs.length+" x values");
		double[] xs = new double[result.commonRotationSpecs.length];
		for (int i=0; i<xs.length; i++)
			xs[i] = scatterQuantity.getValue(this, result.commonRotationSpecs[i]);

		String[] yLabels = { "Standard Deviation", "Median SA" };
		String[] yPrefixes = { "std_dev", "median", "residual" };
		Range[] yRanges = { new Range (0, 1), null };
		boolean[] yLogs = { false, true };
		
		File[][] ret = new File[yLabels.length][periods.length];
		
		for (int yI=0; yI<yLabels.length; yI++) {
			if (type.stdDevOfMedians && yI == 0)
				continue;
			if (!type.stdDevOfMedians && yI == 1)
				continue;
			for (int p=0; p<periods.length; p++) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				
				for (int i=0; i<xs.length; i++) {
					if (yI == 0)
						xy.set(xs[i], result.getResidualStdDevSet(periods[p]).stdDevs[i]);
					else if (yI == 1)
						xy.set(xs[i], Math.exp(result.getMedianStdDevSet(periods[p]).medians[i]));
				}
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
				
				String periodStr = optionalDigitDF.format(periods[p])+"s";
				
				PlotSpec spec = new PlotSpec(funcs, chars, periodStr+" "+scatterQuantity.name+" Scatter", scatterQuantity.name, yLabels[yI]);
				spec.setLegendVisible(false);

				PlotPreferences plotPrefs = PlotPreferences.getDefault();
				plotPrefs.setTickLabelFontSize(18);
				plotPrefs.setAxisLabelFontSize(20);
				plotPrefs.setPlotLabelFontSize(21);
				plotPrefs.setBackgroundColor(Color.WHITE);

				HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);

				gp.drawGraphPanel(spec, false, yLogs[yI], null, yRanges[yI]);
				gp.getChartPanel().setSize(800, 450);
				String myPrefix = prefix+"_"+scatterQuantity.prefix+"_"+periodStr+"_"+yPrefixes[yI];
				ret[yI][p] = new File(resourcesDir, myPrefix+".png");
				gp.saveAsPNG(ret[yI][p].getAbsolutePath());
			}
		}
		
		return ret;
	}
	
	private enum MagDistPlotType {
		SIM_STD_DEV("Simulated", "sim"),
		GMPE_STD_DEV("GMPE", "gmpe"),
		SIM_GMPE_STD_DEV_DIFF("Simulated - GMPE", "sim_gmpe_diff"),
		SIM_MEDIAN("Sim Median SA", "sim_median");
		
		private String name;
		private String prefix;

		private MagDistPlotType(String name, String prefix) {
			this.name = name;
			this.prefix = prefix;
		}
		
		public String getName() {
			return name;
		}
		
		public String getPrefix() {
			return prefix;
		}
		
		private boolean hasGMPE() {
			return this == MagDistPlotType.GMPE_STD_DEV || this == MagDistPlotType.SIM_GMPE_STD_DEV_DIFF;
		}
		
		private boolean hasMedian() {
			return this == MagDistPlotType.SIM_MEDIAN;
		}
		
		private boolean hasStdDev() {
			return this == MagDistPlotType.SIM_STD_DEV || this == MagDistPlotType.GMPE_STD_DEV
					|| this == MagDistPlotType.SIM_GMPE_STD_DEV_DIFF;
		}
		
		private boolean isDiff() {
			return this == MagDistPlotType.SIM_GMPE_STD_DEV_DIFF;
		}
	}
	
	private Map<MagDistPlotType, File[]> plotMagDistCheckerboard(File resourcesDir, String prefix, VariabilityType type, Site site, double[] periods)
			throws IOException {
		List<Double> magnitudes = new ArrayList<>(magConfigs.keySet());
		Collections.sort(magnitudes);
		
		double deltaMag = magnitudes.get(1) - magnitudes.get(0);
		double minMag = magnitudes.get(0);
		
		double deltaDist = distances.get(1) - distances.get(0);
		double minDist = distances.get(0);
		
		Map<MagDistPlotType, EvenlyDiscrXYZ_DataSet[]> xyzsMap = new HashMap<>();
		EvenlyDiscrXYZ_DataSet refXYZ = null;
		
		MagDistPlotType[] allTypes;
		if (site == null)
			allTypes = MagDistPlotType.values();
		else
			allTypes = new MagDistPlotType[] { MagDistPlotType.SIM_STD_DEV };
		
		for (MagDistPlotType plotType : allTypes) {
			if (plotType.hasGMPE()) {
				if (gmpes == null)
					continue;
				if (plotType.hasStdDev() && type.gmpeStdDevType == null)
					continue;
			}
			EvenlyDiscrXYZ_DataSet[] xyzs = new EvenlyDiscrXYZ_DataSet[periods.length];
			for (int p=0; p<periods.length; p++)
				xyzs[p] = new EvenlyDiscrXYZ_DataSet(distances.size(), magnitudes.size(), minDist, minMag, deltaDist, deltaMag);
			xyzsMap.put(plotType, xyzs);
			if (refXYZ == null)
				refXYZ = xyzs[0];
		}
		
		for (int xInd=0; xInd<distances.size(); xInd++) {
			float distance = distances.get(xInd);
			double x = refXYZ.getX(xInd);
			Preconditions.checkState(distance == (float)x, "Expected distance bin %s to be %s, have %s", xInd, distance, x);
			for (int yInd=0; yInd<magnitudes.size(); yInd++) {
				double magnitude = magnitudes.get(yInd);
				double y = refXYZ.getY(yInd);
				Preconditions.checkState((float)magnitude == (float)y, "Expected mag bin %s to be %s, have %s", yInd, magnitude, y);
				VarGroupingKey key = new VarGroupingKey(type, magnitude, distance, site);
				VariabilityResult varResult;
				try {
					varResult = varResultCache.get(key);
				} catch (ExecutionException e) {
					if (e.getCause() instanceof IOException)
						throw (IOException)e.getCause();
					if (e.getCause() instanceof RuntimeException)
						throw (RuntimeException)e.getCause();
					throw ExceptionUtils.asRuntimeException(e);
				}
				double[] stdDevs = new double[periods.length];
				for (int p=0; p<periods.length; p++) {
					if (type.stdDevOfMedians)
						stdDevs[p] = varResult.getMedianStdDevSet(periods[p]).stdDev;
					else
						stdDevs[p] = varResult.getResidualStdDevSet(periods[p]).total;
				}
				GMPE_Result gmpeResult = null;
				if (gmpes != null) {
					GMPE_Result[] gmpeResults = new GMPE_Result[gmpes.length];
					for (int g=0; g<gmpes.length; g++) {
						try {
							gmpeResults[g] = gmpeResultCache.get(
									new GMPE_GroupingKey(gmpes[g], site, magnitude, distance));
						} catch (ExecutionException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					}
					gmpeResult = packGMs(gmpeResults);
				}
				for (MagDistPlotType plotType : xyzsMap.keySet()) {
					EvenlyDiscrXYZ_DataSet[] xyzs = xyzsMap.get(plotType);
					double[] gmpeVars = gmpeResult == null || type.gmpeStdDevType == null ?
							null : gmpeResult.getVariability(type.gmpeStdDevType);
					for (int p=0; p<periods.length; p++) {
						int gmpeP = Doubles.indexOf(calcPeriods, periods[p]);
						double gmpeSD = gmpeVars == null ? -1 : gmpeVars[gmpeP];
						double val;
						switch (plotType) {
						case SIM_STD_DEV:
							val = stdDevs[p];
							break;
						case GMPE_STD_DEV:
							val = gmpeSD;
							break;
						case SIM_GMPE_STD_DEV_DIFF:
							val = stdDevs[p] - gmpeSD;
							break;
						case SIM_MEDIAN:
							val = Math.log10(Math.exp(varResult.getLogValues(periods[p]).median));
							break;

						default:
							throw new IllegalStateException("unknown: "+plotType);
						}
						xyzs[p].set(xInd, yInd, val);
					}
				}
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		Map<MagDistPlotType, File[]> retMap = new HashMap<>();
		
		for (MagDistPlotType plotType : xyzsMap.keySet()) {
			File[] ret = new File[periods.length];
			retMap.put(plotType, ret);
			
			CPT cpt;
			String zLabel;
			switch (plotType) {
			case SIM_STD_DEV:
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
				zLabel = plotType.getName()+" "+type.symbol;
				break;
			case GMPE_STD_DEV:
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
				zLabel = plotType.getName()+" "+type.gmpeStdDevType.symbol;
				break;
			case SIM_GMPE_STD_DEV_DIFF:
				cpt = GMT_CPT_Files.GMT_POLAR.instance().rescale(-0.5, 0.5);
				zLabel = plotType.getName()+" "+type.symbol;
				break;
			case SIM_MEDIAN:
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-4, 0d);
				zLabel = "Log10 "+plotType.getName();
				break;

			default:
				throw new IllegalStateException("unknown: "+plotType);
			}
			
			EvenlyDiscrXYZ_DataSet[] xyzs = xyzsMap.get(plotType);
			
			for (int p=0; p<periods.length; p++) {
				String periodStr;
				if (periods[p] == Math.round(periods[p]))
					periodStr = (int)periods[p]+"s";
				else
					periodStr = (float)periods[p]+"s";
				String title = periodStr+" "+type.name+" ("+type.symbol+"), "+plotType.getName();
				String myPrefix = prefix+"_"+periodStr+"_"+plotType.getPrefix();
				String xAxisLabel = BBP_PartBValidationConfig.DIST_JB ? "DistanceJB" : "Distance";
				XYZPlotSpec xyzSpec = new XYZPlotSpec(xyzs[p], cpt, title, xAxisLabel, "Magnitude", zLabel);
				XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
				xyzGP.drawPlot(xyzSpec, false, false, new Range(minDist-0.5*deltaDist, distances.get(distances.size()-1)+0.5*deltaDist),
						new Range(minMag-0.5*deltaMag, magnitudes.get(magnitudes.size()-1)+0.5*deltaMag));
				xyzGP.getChartPanel().getChart().setBackgroundPaint(Color.WHITE);
				xyzGP.getChartPanel().setSize(700, 550);
				File file = new File(resourcesDir, myPrefix+".png");
				xyzGP.saveAsPNG(file.getAbsolutePath());
				ret[p] = file;
			}
		}
		return retMap;
	}
	
	private Map<VariabilityType, File[]> plotQuantityDistCheckerboards(File resourcesDir, String prefix, Quantity quantity,
			Collection<VariabilityType> types, Site site, Double magnitude, double[] periods, boolean replot) throws IOException {
		Map<VariabilityType, File[]> retMap = new HashMap<>();
		
		List<VariabilityType> plotTypes = new ArrayList<>(types);
		plotTypes.add(null);
		
		for (VariabilityType type : plotTypes) {
			File[] ret = new File[periods.length];
			retMap.put(type, ret);
			
			for (int p=0; p<periods.length; p++) {
				String periodStr;
				if (periods[p] == Math.round(periods[p]))
					periodStr = (int)periods[p]+"s";
				else
					periodStr = (float)periods[p]+"s";
				String myPrefix;
				if (type == null)
					myPrefix = prefix+"_"+periodStr+"_median_sa";
				else
					myPrefix = prefix+"_"+periodStr+"_"+type.prefix;
				ret[p] = new File(resourcesDir, myPrefix+".png");
				replot = replot || !ret[p].exists();
			}
		}
		if (!replot)
			return retMap;
		
		RotatedRupVariabilityConfig<E> config = magConfigs.get(magnitude);
		List<Float> values = config.getValues(Float.class, quantity);
		// will already be sorted
		
		Preconditions.checkState(values.size() > 1, "%s doesn't vary!", quantity);
		
		double deltaValue = values.get(1) - values.get(0);
		double minValue = values.get(0);
		
		double deltaDist = distances.get(1) - distances.get(0);
		double minDist = distances.get(0);
		
		Map<VariabilityType, EvenlyDiscrXYZ_DataSet[]> xyzsMap = new HashMap<>();
		EvenlyDiscrXYZ_DataSet refXYZ = null;
		
		for (VariabilityType type : types) {
			EvenlyDiscrXYZ_DataSet[] xyzs = new EvenlyDiscrXYZ_DataSet[periods.length];
			for (int p=0; p<periods.length; p++)
				xyzs[p] = new EvenlyDiscrXYZ_DataSet(values.size(), distances.size(), minValue, minDist, deltaValue, deltaDist);
			xyzsMap.put(type, xyzs);
			if (refXYZ == null)
				refXYZ = xyzs[0];
		}
		EvenlyDiscrXYZ_DataSet[] medianXYZs = new EvenlyDiscrXYZ_DataSet[periods.length];
		for (int p=0; p<periods.length; p++)
			medianXYZs[p] = new EvenlyDiscrXYZ_DataSet(values.size(), distances.size(), minValue, minDist, deltaValue, deltaDist);
		xyzsMap.put(null, medianXYZs);
		
		for (int yInd=0; yInd<distances.size(); yInd++) {
			float distance = distances.get(yInd);
			double y = refXYZ.getY(yInd);
			Preconditions.checkState(distance == (float)y, "Expected distance bin %s to be %s, have %s", yInd, distance, y);
			for (int xInd=0; xInd<values.size(); xInd++) {
				Float value = values.get(xInd);
				double x = refXYZ.getX(xInd);
				Preconditions.checkState((float)value == (float)x, "Expected value bin %s to be %s, have %s", xInd, value, x);
				double[] medians = null;
				Map<Quantity, Object> restrictMap = new HashMap<>();
				restrictMap.put(quantity, value);
				for (VariabilityType type : types) {
					VarGroupingKey key = new VarGroupingKey(type, magnitude, distance,
							site == null ? sites.toArray(new Site[0]) : new Site[] {site});
					key.restrictTo = restrictMap;
					VariabilityResult varResult;
					try {
						varResult = varResultCache.get(key);
					} catch (ExecutionException e) {
						if (e.getCause() instanceof IOException)
							throw (IOException)e.getCause();
						if (e.getCause() instanceof RuntimeException)
							throw (RuntimeException)e.getCause();
						throw ExceptionUtils.asRuntimeException(e);
					}
					if (medians == null) {
						medians = new double[periods.length];
						for (int p=0; p<periods.length; p++)
							medians[p] = Math.log10(Math.exp(varResult.getLogValues(periods[p]).median));
					}
					double[] stdDevs = new double[periods.length];
					for (int p=0; p<periods.length; p++) {
						if (type.stdDevOfMedians)
							stdDevs[p] = varResult.getMedianStdDevSet(periods[p]).stdDev;
						else
							stdDevs[p] = varResult.getResidualStdDevSet(periods[p]).total;
					}
					EvenlyDiscrXYZ_DataSet[] xyzs = xyzsMap.get(type);
					for (int p=0; p<periods.length; p++)
						xyzs[p].set(xInd, yInd, stdDevs[p]);
				}
				// now do median
				EvenlyDiscrXYZ_DataSet[] xyzs = xyzsMap.get(null);
				for (int p=0; p<periods.length; p++)
					xyzs[p].set(xInd, yInd, medians[p]);
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(minValue-0.5*deltaValue, values.get(values.size()-1)+0.5*deltaValue);
		Range yRange = new Range(minDist-0.5*deltaDist, distances.get(distances.size()-1)+0.5*deltaDist);
		
		for (VariabilityType type : xyzsMap.keySet()) {
			File[] ret = retMap.get(type);
			
			CPT cpt;
			String zLabel;
			if (type == null) {
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-4, 0d);
				zLabel = "Log10 Median SA";
			} else {
				cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
				zLabel = type.symbol;
			}
			
			EvenlyDiscrXYZ_DataSet[] xyzs = xyzsMap.get(type);
			
			for (int p=0; p<periods.length; p++) {
				String periodStr;
				if (periods[p] == Math.round(periods[p]))
					periodStr = (int)periods[p]+"s";
				else
					periodStr = (float)periods[p]+"s";
				String title;
				if (type == null)
					title = periodStr+" "+quantity.getName()+", Median SA";
				else
					title = periodStr+" "+quantity.getName()+", "+type.name+" ("+type.symbol+")";
				
				String xAxisLabel = quantity.getName();
				String yAxisLabel = BBP_PartBValidationConfig.DIST_JB ? "DistanceJB" : "Distance";
				
				XYZPlotSpec xyzSpec = new XYZPlotSpec(xyzs[p], cpt, title, xAxisLabel, yAxisLabel, zLabel);
				XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
				xyzGP.drawPlot(xyzSpec, false, false, xRange, yRange);
				xyzGP.getChartPanel().getChart().setBackgroundPaint(Color.WHITE);
				xyzGP.getChartPanel().setSize(700, 550);
				xyzGP.saveAsPNG(ret[p].getAbsolutePath());
			}
		}
		return retMap;
	}
	
	private Map<VariabilityType, File[]> plotQuantityDistFuncs(File resourcesDir, String prefix, Quantity quantity,
			Collection<VariabilityType> types, Site site, Double magnitude, double[] periods, List<Float> distances, boolean replot)
					throws IOException {
		Map<VariabilityType, File[]> retMap = new HashMap<>();
		
		List<VariabilityType> plotTypes = new ArrayList<>(types);
		plotTypes.add(null);
		
		for (VariabilityType type : plotTypes) {
			File[] ret = new File[periods.length];
			retMap.put(type, ret);
			
			for (int p=0; p<periods.length; p++) {
				String periodStr;
				if (periods[p] == Math.round(periods[p]))
					periodStr = (int)periods[p]+"s";
				else
					periodStr = (float)periods[p]+"s";
				String myPrefix;
				if (type == null)
					myPrefix = prefix+"_"+periodStr+"_median_sa";
				else
					myPrefix = prefix+"_"+periodStr+"_"+type.prefix;
				ret[p] = new File(resourcesDir, myPrefix+".png");
				replot = replot || !ret[p].exists();
			}
		}
		if (!replot)
			return retMap;
		
		RotatedRupVariabilityConfig<E> config = magConfigs.get(magnitude);
		List<Float> values = config.getValues(Float.class, quantity);
		// will already be sorted
		
		Preconditions.checkState(values.size() > 1, "%s doesn't vary!", quantity);

		Map<VariabilityType, DiscretizedFunc[][]> xysMap = new HashMap<>();
		Map<VariabilityType, DiscretizedFunc[][]> gmpeXYsMap = new HashMap<>();
		
		for (boolean gmpe : new boolean[] { false, true }) {
			if (gmpe && this.gmpes == null)
				continue;
			for (VariabilityType type : types) {
				if (gmpe && type.gmpeStdDevType == null)
					continue;
				DiscretizedFunc[][] xys = new DiscretizedFunc[distances.size()][periods.length];
				for (int d=0; d<distances.size(); d++)
					for (int p=0; p<periods.length; p++)
						xys[d][p] = new ArbitrarilyDiscretizedFunc();
				if (gmpe)
					gmpeXYsMap.put(type, xys);
				else
					xysMap.put(type, xys);
			}
			DiscretizedFunc[][] medianXYs = new DiscretizedFunc[distances.size()][periods.length];
			for (int d=0; d<distances.size(); d++)
				for (int p=0; p<periods.length; p++)
					medianXYs[d][p] = new ArbitrarilyDiscretizedFunc();
			if (gmpe)
				gmpeXYsMap.put(null, medianXYs);
			else
				xysMap.put(null, medianXYs);
		}
		
		for (int d=0; d<distances.size(); d++) {
			Float distance = distances.get(d);
			for (int yInd=0; yInd<values.size(); yInd++) {
				Float value = values.get(yInd);
				double[] medians = null;
				Map<Quantity, Object> restrictMap = new HashMap<>();
				restrictMap.put(quantity, value);
				GMPE_Result gmpeResult = null;
				if (gmpes != null) {
					GMPE_Result[] gmpeResults = new GMPE_Result[gmpes.length];
					for (int g=0; g<gmpes.length; g++) {
						try {
							GMPE_GroupingKey gmpeKey = new GMPE_GroupingKey(gmpes[g], site, magnitude, distance);
							gmpeKey.restrictTo = restrictMap;
							gmpeResults[g] = gmpeResultCache.get(gmpeKey);
						} catch (ExecutionException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					}
					gmpeResult = packGMs(gmpeResults);
				}
				for (VariabilityType type : types) {
					VarGroupingKey key = new VarGroupingKey(type, magnitude, distance,
							site == null ? sites.toArray(new Site[0]) : new Site[] {site});
					key.restrictTo = restrictMap;
					VariabilityResult varResult;
					try {
						varResult = varResultCache.get(key);
					} catch (ExecutionException e) {
						if (e.getCause() instanceof IOException)
							throw (IOException)e.getCause();
						if (e.getCause() instanceof RuntimeException)
							throw (RuntimeException)e.getCause();
						throw ExceptionUtils.asRuntimeException(e);
					}
					if (medians == null && varResult.logVals != null) {
						medians = new double[periods.length];
						for (int p=0; p<periods.length; p++)
							medians[p] = Math.exp(varResult.getLogValues(periods[p]).median);
					}
					double[] stdDevs = new double[periods.length];
					for (int p=0; p<periods.length; p++) {
						if (type.stdDevOfMedians)
							stdDevs[p] = varResult.getMedianStdDevSet(periods[p]).stdDev;
						else
							stdDevs[p] = varResult.getResidualStdDevSet(periods[p]).total;
					}
					DiscretizedFunc[][] xys = xysMap.get(type);
					for (int p=0; p<periods.length; p++)
						xys[d][p].set(value, stdDevs[p]);
					if (type.gmpeStdDevType != null && gmpeResult != null) {
						DiscretizedFunc[][] gmpeXYs = gmpeXYsMap.get(type);
						double[] gmpeVars = gmpeResult.getVariability(type.gmpeStdDevType);
						for (int p=0; p<periods.length; p++) {
							int gmpeP = Doubles.indexOf(calcPeriods, periods[p]);
							gmpeXYs[d][p].set(value, gmpeVars[gmpeP]);
						}
					}
				}
				// now do median
				if (medians != null) {
					DiscretizedFunc[][] xys = xysMap.get(null);
					for (int p=0; p<periods.length; p++)
						xys[d][p].set(value, medians[p]);
				} else {
					xysMap.remove(null);
				}
				if (gmpeResult != null) {
					DiscretizedFunc[][] gmpeXYs = gmpeXYsMap.get(null);
					for (int p=0; p<periods.length; p++) {
						int gmpeP = Doubles.indexOf(calcPeriods, periods[p]);
						gmpeXYs[d][p].set(value, Math.exp(gmpeResult.logMedian[gmpeP]));
					}
				}
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		double[] distArray = Doubles.toArray(distances);
		double minDist = StatUtils.min(distArray);
		double maxDist = StatUtils.max(distArray);
		if (maxDist == minDist)
			maxDist++;
		CPT distCPT = new CPT(minDist, maxDist, Color.GRAY, Color.BLACK);
		
		for (VariabilityType type : xysMap.keySet()) {
			File[] ret = retMap.get(type);
			
			String xAxisLabel = quantity.getName();
			String yAxisLabel;
			if (type == null)
				yAxisLabel = "Log10 Median SA";
			else
				yAxisLabel = type.symbol;
			
			DiscretizedFunc[][] xys = xysMap.get(type);
			DiscretizedFunc[][] gmpeXYs = gmpeXYsMap.get(type);
			
			for (int p=0; p<periods.length; p++) {
				String periodStr;
				if (periods[p] == Math.round(periods[p]))
					periodStr = (int)periods[p]+"s";
				else
					periodStr = (float)periods[p]+"s";
				String title;
				if (type == null)
					title = periodStr+" "+quantity.getName()+", Median SA";
				else
					title = periodStr+" "+quantity.getName()+", "+type.name+" ("+type.symbol+")";
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				for (int d=0; d<distances.size(); d++) {
					xys[d][p].setName(optionalDigitDF.format(distances.get(d))+" km");
					Color color = distCPT.getColor(distances.get(d));
					
					funcs.add(xys[d][p]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
					
					if (gmpeXYs != null) {
						if (type == null)
							gmpeXYs[d][p].setName("GMPE");
						else
							gmpeXYs[d][p].setName("GMPE "+type.gmpeStdDevType.symbol);
						
						funcs.add(gmpeXYs[d][p]);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, color));
					}
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
				spec.setLegendVisible(distances.size() > 1);
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
				
				Range xRange = quantity == Quantity.SOURCE_AZIMUTH || quantity == Quantity.SITE_TO_SOURTH_AZIMUTH ?
						new Range(0, 360) : null;
				Range yRange;
				boolean yLog;
				if (type == null) {
					yRange = new Range(1e-4, 1e0);
					yLog = true;
				} else {
					yRange = new Range(0d, 1d);
					yLog = false;
				}

				gp.drawGraphPanel(spec, false, yLog, xRange, yRange);
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(ret[p].getAbsolutePath());
			}
		}
		return retMap;
	}
	
	protected abstract void plotExample(File resourcesDir, String prefix, double distance, List<Quantity> variedQuantities)
			throws IOException;
	
	private File writeSiteCSV(File resourcesDir, Site site, Double magnitude, Float distance) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = Lists.newArrayList(Quantity.EVENT_ID.getName(), Quantity.SOURCE_AZIMUTH.getName(),
				Quantity.SITE_TO_SOURTH_AZIMUTH.getName());
		for (double period : calcPeriods)
			header.add(optionalDigitDF.format(period)+"s SA");
		
		csv.addLine(header);
		
		SimulationRotDProvider<RotationSpec> prov = magProvs.get(magnitude);
		RotatedRupVariabilityConfig<E> config = magConfigs.get(magnitude);
		
		for (RotationSpec rotation : config.getRotationsForQuantities(Quantity.SITE, site, Quantity.DISTANCE, distance)) {
			DiscretizedFunc spectrum = prov.getRotD50(site, rotation, 0);
			List<String> line = new ArrayList<>();
			line.add(rotation.eventID+"");
			line.add(rotation.sourceAz == null ? "0.0" : rotation.sourceAz.toString());
			line.add(rotation.siteToSourceAz == null ? "0.0" : rotation.siteToSourceAz.toString());
			for (double period : calcPeriods)
				line.add((float)spectrum.getInterpolatedY(period)+"");
			csv.addLine(line);
		}
		
		File outFile = new File(resourcesDir, "sa_"+site.getName()+"_m"+magnitude.floatValue()+"_"+distance+"km.csv.gz");
		csv.writeToStream(new GZIPOutputStream(new FileOutputStream(outFile)));
		
		return outFile;
	}
	
	private File writeCombinedCSV(File resourcesDir, Double magnitude, double period) throws IOException {
		List<String> logValueHeader = Lists.newArrayList("Event ID", "Magnitude", "Source Azimuth (degrees)",
				"Site-To-Source Azimuth (degrees)", "Rrup Distance (km)");
		List<String> residualHeader = new ArrayList<>(logValueHeader);
		List<String> eventTermHeader = Lists.newArrayList("Event ID", "Magnitude", "Rrup Distance (km)");
		List<Site[]> siteBundles = new ArrayList<>();
		
		String valPeriodStr = "ln("+optionalDigitDF.format(period)+"s RotD50)";
		String residualPeriodStr = optionalDigitDF.format(period)+"s Residual";
		String eventTermPeriodStr = optionalDigitDF.format(period)+"s Event Term";
		if (sites.size() > 1) {
			if (!vs30SiteBundles.isEmpty()) {
				List<Float> vs30s = new ArrayList<>(vs30SiteBundles.keySet());
				Collections.sort(vs30s);
				for (Float vs30 : vs30s) {
					Site[] siteArray = vs30SiteBundles.get(vs30).toArray(new Site[0]);
					siteBundles.add(siteArray);
					String prefix = "Average of "+siteArray.length+" Sites w/ Vs30="+optionalDigitDF.format(vs30)+", ";
					logValueHeader.add(prefix+valPeriodStr);
					residualHeader.add(prefix+residualPeriodStr);
					eventTermHeader.add(prefix+eventTermPeriodStr);
				}
			}
			for (Site site : sites) {
				logValueHeader.add(site.getName()+", "+valPeriodStr);
				residualHeader.add(site.getName()+", "+residualPeriodStr);
				eventTermHeader.add(site.getName()+", "+eventTermPeriodStr);
				siteBundles.add(new Site[] {site});
			}
		} else {
			logValueHeader.add(valPeriodStr);
			residualHeader.add(residualPeriodStr);
			eventTermHeader.add(eventTermPeriodStr);
			siteBundles.add(new Site[] {sites.get(0)});
		}
		
		CSVFile<String> logValCSV = new CSVFile<>(true);
		CSVFile<String> residualCSV = new CSVFile<>(true);
		CSVFile<String> eventTermCSV = new CSVFile<>(true);
		logValCSV.addLine(logValueHeader);
		residualCSV.addLine(residualHeader);
		eventTermCSV.addLine(logValueHeader);
		
		SimulationRotDProvider<RotationSpec> prov = magProvs.get(magnitude);
		RotatedRupVariabilityConfig<E> config = magConfigs.get(magnitude);
		
		// TODO
		
//		Site site0 = sites.get(0);
//		for (RotationSpec rotation : config.getRotationsForQuantities(Quantity.SITE, site0)) {
//			DiscretizedFunc spectrum = prov.getRotD50(site, rotation, 0);
//			List<String> logValLine = new ArrayList<>();
//			List<String> residualLine = new ArrayList<>();
//			List<String> eventTermLine = new ArrayList<>();
//			line.add(rotation.eventID+"");
//			line.add(rotation.sourceAz == null ? "0.0" : rotation.sourceAz.toString());
//			line.add(rotation.siteToSourceAz == null ? "0.0" : rotation.siteToSourceAz.toString());
//			for (double period : calcPeriods)
//				line.add((float)spectrum.getInterpolatedY(period)+"");
//			csv.addLine(line);
//		}
//		
//		File outFile = new File(resourcesDir, "sa_"+site.getName()+"_m"+magnitude.floatValue()+"_"+distance+"km.csv.gz");
//		csv.writeToStream(new GZIPOutputStream(new FileOutputStream(outFile)));
		
		return null;
	}
	
	private enum EventDisaggQuantity {
		V_PROP("Propagation Velocity", "v_prop", "m/s") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.calcVprop(event);
			}
		},
		MAG("Mag", "mag", null) {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getMag(event);
			}
		},
		AREA("Log10(Area)", "area", "m^2") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getArea(event);
			}
		},
		MAX_SLIP("Max Slip", "max_slip", "m") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getMaxSlip(event);
			}
		},
		MEAN_SLIP("Mean Slip", "mean_slip", "m") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getMeanSlip(event);
			}
		},
		SLIP_STD_DEV("Slip Std Dev", "slip_std_dev", "m") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getSlipStdDev(event);
			}
		},
		MID_SEIS_MEAN_SLIP("Mid-Seismogenic Mean Slip", "mid_seis_mean_slip", "m") {
			@Override
			public <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event) {
				return pageGen.getMeanMidSeisSlip(event);
			}
		};
		
		private String name;
		private String prefix;
		private String units;
		
		private EventDisaggQuantity(String name, String prefix, String units) {
			this.name = name;
			this.prefix = prefix;
			this.units = units;
		}
		
		public abstract <R> double getValue(RotatedRupVariabilityPageGen<R> pageGen, R event);
	}
	
	private File plotEventTermScatter(EventDisaggQuantity q, Double mag,
			Float distance, double period, Site[] sites, File resourcesDir)
					throws IOException {
		RotatedRupVariabilityConfig<E> config = magConfigs.get(mag);
		
		List<Integer> eventIDs = config.getValues(Integer.class, Quantity.EVENT_ID);
		
		double[] eventTerms = new double[eventIDs.size()];
		double[] values = new double[eventIDs.size()];
		
		int periodIndex = Doubles.indexOf(calcPeriods, period);
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		
		boolean allNaN = true;
		try {
			for (int i=0; i<eventTerms.length; i++) {
				eventTerms[i] = eventTermCache.get(
						new EventTermKey(eventIDs.get(i), mag, distance, sites)).eventTerms[periodIndex];
				values[i] = q.getValue(this, getEvent(eventIDs.get(i)));
				
				allNaN = allNaN && !Double.isFinite(values[i]);
				
				scatter.set(values[i], eventTerms[i]);
			}
		} catch (ExecutionException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		if (allNaN)
			return null;

		double maxX = StatUtils.max(values);
		double minX = StatUtils.min(values);
		if (minX == maxX)
			return null;
		double spanX = maxX - minX;
		Range xRange = new Range(minX - 0.05*spanX, maxX + 0.05*spanX);
		
		HistogramFunction xHist = calcHist(values, xRange, 20);
		
		double maxY = StatUtils.max(eventTerms);
		double minY = StatUtils.min(eventTerms);
		double spanY = maxY - minY;
		Range yRange = new Range(minY - 0.05*spanY, maxY + 0.05*spanY);
		
		HistogramFunction yHist = calcHist(eventTerms, yRange, 20);

//		double maxHistY = Math.max(yHist.getMaxY(), xHist.getMaxY());
//		double scalarX = 0.25*yRange.getLength()/maxHistY;
//		double scalarY = 0.25*xRange.getLength()/maxHistY;
		double scalarX = 0.25*yRange.getLength()/xHist.getMaxY();
		double scalarY = 0.25*xRange.getLength()/yHist.getMaxY();
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 5f, Color.BLACK));
		
		SimpleRegression regression = new SimpleRegression();
		for (Point2D pt : scatter)
			regression.addData(pt.getX(), pt.getY());
		double b = regression.getIntercept();
		double m = regression.getSlope();
		DefaultXY_DataSet fit = new DefaultXY_DataSet();
		fit.set(xRange.getLowerBound(), m*xRange.getLowerBound()+b);
		fit.set(xRange.getUpperBound(), m*xRange.getUpperBound()+b);
		
		funcs.add(fit);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		String prefix = "event_term_scatter_"+q.prefix;
		if (magConfigs.size() > 0)
			prefix += "_m"+optionalDigitDF.format(mag);
		prefix += "_"+optionalDigitDF.format(distance)+"km";
		if (sites.length == 1)
			prefix += "_"+sites[0].getName();
		else
			prefix += "_"+sites.length+"sites";
		prefix += "_"+optionalDigitDF.format(period)+"s";
		
		String xAxisLabel = q.name;
		if (q.units != null)
			xAxisLabel += " ("+q.units+")";
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", xAxisLabel, "Event Term");
		List<XYAnnotation> anns = new ArrayList<>();
		anns.addAll(buildHistAnns(xHist, false, xRange, yRange, scalarX,
				new Color(255,0,0,80)));
		anns.addAll(buildHistAnns(yHist, true, xRange, yRange, scalarY,
				new Color(0,0,255,80)));
		spec.setPlotAnnotations(anns);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}
	
	private HistogramFunction calcHist(double[] values, Range range, int num) {
		HistogramFunction hist = new HistogramFunction(range.getLowerBound(),
				range.getUpperBound(), num);
		
		for (double value : values)
			hist.add(hist.getClosestXIndex(value), 1d);
		
		return hist;
	}
	
	private List<XYAnnotation> buildHistAnns(HistogramFunction hist, boolean vertical,
			Range xAxis, Range yAxis, double scalar, Color color) {
		List<XYAnnotation> anns = new ArrayList<>();
		
		double baselineValue = vertical ? xAxis.getLowerBound() : yAxis.getLowerBound();
		
		double delta = hist.getDelta();
		double halfDelta = 0.5*delta;
		
		for (int i=0; i<hist.size(); i++) {
			double middle = hist.getX(i);
			double start = middle - halfDelta;
			double end = middle + halfDelta;
			
			double rawHeight = hist.getY(i);
			if (rawHeight == 0d)
				continue;
			double top = baselineValue + scalar*rawHeight;
			
			double x0, y0, x1, y1;
			if (vertical) {
				x0 = baselineValue;
				x1 = top;
				y0 = start;
				y1 = end;
			} else {
				y0 = baselineValue;
				y1 = top;
				x0 = start;
				x1 = end;
			}
			
			anns.add(new XYBoxAnnotation(x0, y0, x1, y1, null, null, color));
		}
		
		return anns;
	}
	protected static boolean DIRECTIVITY_DEBUG = true;
	protected static String DIRECTIVITY_DEBUG_DIRNAME = "directivity_debug";
	
	protected abstract Table<Float, Double, XY_DataSet> calcDirectivityComparisons(Double mag, double[] periods,
			Site[] sites, File resourcesDir);
	
	private boolean plotDirectivityComparison(File resourcesDir, String prefix, XY_DataSet scatter)
			throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, Color.GRAY));
		
		double minX = scatter.getMinX();
		double maxX = scatter.getMaxX();
		
		if (minX == maxX && minX == 0d)
			return false;
		
		EvenlyDiscretizedFunc xFunc = new EvenlyDiscretizedFunc(minX, maxX, 21);
		// now shift so that the edges are at min and max
		double delta = xFunc.getDelta();
		xFunc = new EvenlyDiscretizedFunc(minX + 0.5*delta, xFunc.size()-1, delta);
		
		List<List<Double>> binnedVals = new ArrayList<>();
		for (int i=0; i<xFunc.size(); i++)
			binnedVals.add(new ArrayList<>());
		
		for (Point2D pt : scatter) {
			int index = xFunc.getClosestXIndex(pt.getX());
			binnedVals.get(index).add(pt.getY());
		}
		
		ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc medianFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xFunc.size(); i++) {
			List<Double> vals = binnedVals.get(i);
			if (vals.isEmpty())
				continue;
			double[] array = Doubles.toArray(vals);
			double x = xFunc.getX(i);
			meanFunc.set(x, StatUtils.mean(array));
			medianFunc.set(x, DataUtils.median(array));
		}
		
		meanFunc.setName("Mean");
		medianFunc.setName("Median");
		
		funcs.add(meanFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Color.BLACK));
		funcs.add(medianFunc);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, Color.BLUE));
		
		double absMaxX = Math.max(maxX, -minX);
		minX = -absMaxX - 0.1;
		maxX = absMaxX + 0.1;
		
		Range xRange = new Range(minX, maxX);
		Range yRange = new Range(scatter.getMinY()-0.1, scatter.getMaxY()+0.1);
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(minX, minX);
		oneToOne.set(maxX, maxX);
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ",
				"Bayless & Sommerville (2013) fD", "Within-Event Residual");
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		
		return true;
	}
	
	protected double[] getDirectivityPlotPeriods(double[] periods) {
		if (periods.length > 3) {
			List<Double> keep = new ArrayList<>();
			for (double period : periods)
				if ((float)period <= 3f || (float)period == (float)Math.floor(period))
					keep.add(period);
			periods = Doubles.toArray(keep);
		}
		return periods;
	}
	protected File[] plotIndvRupDirectivity(int eventID, RSQSimEvent oriented, EqkRupture rup,
			Table<Float, Float, double[]> azDistFdMap, Table<Float, Float, Location> azDistLocMap,
			Map<RotationSpec, double[]> diffsMap, double[] periods, File outputDir) throws IOException {
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		periods = getDirectivityPlotPeriods(periods);
		
		String commonPrefix = "event_"+eventID;

		Float chartDistance = distances.get(0);
		float maxMapDistance = 60f;
		
		PlotSpec mapSpec = RupturePlotGenerator.writeMapPlot(null, oriented, null, null, null,
				null, null, rup.getRuptureSurface(), null, null, null, null);
		if (oriented == null) {
			double hypoRadius = 0.02;
			BasicStroke hypoStroke = new BasicStroke(1f);
			Location hypoLoc = rup.getHypocenterLocation();
			XYPolygonAnnotation rectHypoPoly = new XYPolygonAnnotation(
					RupturePlotGenerator.star(hypoLoc.getLongitude(), hypoLoc.getLatitude(), hypoRadius), hypoStroke,
					Color.BLACK, RupturePlotGenerator.HYPO_COLOR);
			mapSpec.getPlotAnnotations().add(rectHypoPoly);
		}
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (Cell<Float, Float, Location> cell : azDistLocMap.cellSet()) {
			if (cell.getColumnKey() > maxMapDistance)
				continue;
			Location loc = cell.getValue();
			latTrack.addValue(loc.getLatitude());
			lonTrack.addValue(loc.getLongitude());
		}
		double minLat = latTrack.getMin();
		double maxLat = latTrack.getMax();
		double minLon = lonTrack.getMin();
		double maxLon = lonTrack.getMax();
		double maxSpan = Math.max(maxLat-minLat, maxLon-minLon);
		double centerLat = 0.5*(minLat+maxLat);
		double centerLon = 0.5*(minLon+maxLon);
//		System.out.println("orig lat range: "+minLat+" "+maxLat);
//		System.out.println("center lat: "+centerLat);
//		System.out.println("maxSpan: "+maxSpan);
		minLat = centerLat - 0.6*maxSpan;
		maxLat = centerLat + 0.6*maxSpan;
		minLon = centerLon - 0.6*maxSpan;
		maxLon = centerLon + 0.6*maxSpan;
		Region mapRegion = new Region(new Location(minLat, minLon), new Location(maxLat, maxLon));
		GriddedRegion gridReg = new GriddedRegion(mapRegion, 0.05, null);
		GriddedGeoDataSet[] xyzs = new GriddedGeoDataSet[periods.length];
		for (int p=0; p<periods.length; p++)
			xyzs[p] = new GriddedGeoDataSet(gridReg, false);
		BaylessSomerville2013DirectivityModifier bs = new BaylessSomerville2013DirectivityModifier();
		for (int i=0; i<xyzs[0].size(); i++) {
			Location loc = gridReg.getLocation(i);
			for (int p=0; p<periods.length; p++) {
				double fd = bs.getFd(rup, loc, periods[p]);
				xyzs[p].set(i, fd);
			}
		}
		
		Map<Float, List<RotationSpec>> azRotsMap = new HashMap<>();
		double minAz = Double.POSITIVE_INFINITY;
		double maxAz = Double.NEGATIVE_INFINITY;
		for (RotationSpec spec : diffsMap.keySet()) {
			Float az = nullAsZero(spec.sourceAz);
			List<RotationSpec> list = azRotsMap.get(az);
			if (list == null) {
				list = new ArrayList<>();
				azRotsMap.put(az, list);
			}
			list.add(spec);
			
			minAz = Math.min(minAz, az);
			maxAz = Math.max(maxAz, az);
		}
		
		DefaultXY_DataSet[] scatters = new DefaultXY_DataSet[periods.length];
		ArbitrarilyDiscretizedFunc[] meanFuncs = new ArbitrarilyDiscretizedFunc[periods.length];
		ArbitrarilyDiscretizedFunc[] medianFuncs = new ArbitrarilyDiscretizedFunc[periods.length];
		ArbitrarilyDiscretizedFunc[] fdFuncs = new ArbitrarilyDiscretizedFunc[periods.length];
		
		for (int p=0; p<scatters.length; p++) {
			scatters[p] = new DefaultXY_DataSet();
			meanFuncs[p] = new ArbitrarilyDiscretizedFunc();
			medianFuncs[p] = new ArbitrarilyDiscretizedFunc();
			fdFuncs[p] = new ArbitrarilyDiscretizedFunc();
		}
		
		boolean hasMultiple = false;
		
		List<XYAnnotation> anns = mapSpec.getPlotAnnotations();
		if (anns == null)
			anns = new ArrayList<>();
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 12);
		
		CPT cpt = GMT_CPT_Files.GMT_POLAR.instance();
		cpt = cpt.rescale(-1, 1d);
		
		for (Float dist : azDistFdMap.columnKeySet()) {
			if (dist > maxMapDistance)
				continue;
			
			for (Float az : azRotsMap.keySet()) {
				List<double[]> vals = new ArrayList<>();
				for (RotationSpec rot : azRotsMap.get(az)) {
					if (!rot.distance.equals(dist))
						continue;
					double[] diffs = diffsMap.get(rot);
					vals.add(diffs);
				}
				hasMultiple = hasMultiple || vals.size() > 1;
				for (int p=0; p<periods.length; p++) {
					double fd = azDistFdMap.get(az, dist)[p];
					double[] valsArray = new double[vals.size()];
					for (int i=0; i<vals.size(); i++) {
						valsArray[i] = vals.get(i)[p];
						scatters[p].set(az, valsArray[i]);
					}
					double mean = StatUtils.mean(valsArray);
					double median = DataUtils.median(valsArray);
					if (dist.equals(chartDistance)) {
						fdFuncs[p].set(az, fd);
						meanFuncs[p].set(az, mean);
						medianFuncs[p].set(az, median);
					}
					
					if (p == 0) {
						Location loc = azDistLocMap.get(az, dist);
						double lat = loc.getLatitude();
						double lon = loc.getLongitude();
						
						Color c = cpt.getColor((float)mean);
						double width = maxSpan/15d;
						double d = width * 0.5;
						double x0 = lon-d;
						double x1 = lon+d;
						double y0 = lat-d;
						double y1 = lat+d;
						XYBoxAnnotation boxAnn = new XYBoxAnnotation(x0, y0, x1, y1,
								PlotLineType.SOLID.buildStroke(1f), Color.BLACK, c);
						
						anns.add(boxAnn);
						XYTextAnnotation textAnn = new XYTextAnnotation(
								optionalDigitDF.format(az), lon, lat);
						textAnn.setFont(font);
						textAnn.setTextAnchor(TextAnchor.CENTER);
						anns.add(textAnn);
					}
				}
			}
		}
		
		Range yRange = new Range(-2d, 2d);
		Range xRange = new Range(minAz-5d, maxAz+5d);
		
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		List<Range> yRanges = new ArrayList<>();
		
		List<PlotSpec> specs = new ArrayList<>();
		for (int p=0; p<periods.length; p++) {
			yRanges.add(yRange);
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (hasMultiple) {
				scatters[p].setName("Simulated Data");
				funcs.add(scatters[p]);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, Color.GRAY));
				
				meanFuncs[p].setName("Simulated Mean");
				funcs.add(meanFuncs[p]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
				
				medianFuncs[p].setName("Simulated Median");
				funcs.add(medianFuncs[p]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
			} else {
				meanFuncs[p].setName("Simulated Data");
				funcs.add(meanFuncs[p]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			}
			
			fdFuncs[p].setName("B&S fD");
			funcs.add(fdFuncs[p]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Event "+eventID, "Source Azimuth",
					optionalDigitDF.format(periods[p])+"s Residual");
			spec.setLegendVisible(p == 0);
			specs.add(spec);
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);

		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(outputDir, commonPrefix+"_residuals.png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		
		File[] ret = new File[periods.length+1];
		ret[ret.length-1] = pngFile;
		
		for (int p=0; p<periods.length; p++) {
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyzs[p], cpt, " ", "Latitude", "Longitude",
					optionalDigitDF.format(periods[p])+"s Residual");
			xyzSpec.setXYElems((List<XY_DataSet>)mapSpec.getPlotElems());
			xyzSpec.setXYChars(mapSpec.getChars());
			xyzSpec.setPlotAnnotations(anns);
			
			XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
			double spacing = gridReg.getSpacing();
			xyzGP.drawPlot(xyzSpec, false, false, new Range(minLon-0.5*spacing, maxLon+0.5*spacing),
					new Range(minLat-0.5*spacing, maxLat+0.5*spacing));
			
			double tick;
			if (maxSpan > 3d)
				tick = 1d;
			else if (maxSpan > 1.5d)
				tick = 0.5;
			else if (maxSpan > 0.8)
				tick = 0.25;
			else
				tick = 0.1;
			TickUnits tus = new TickUnits();
			TickUnit tu = new NumberTickUnit(tick);
			tus.add(tu);
			xyzGP.getXAxis().setStandardTickUnits(tus);
			xyzGP.getYAxis().setStandardTickUnits(tus);
			
			// write plot
			xyzGP.getChartPanel().setSize(800, 800);
			String fName = commonPrefix+"_map_"+optionalDigitDF.format(periods[p])+"s.png";
			File file = new File(outputDir, fName);
			ret[p] = file;
			xyzGP.saveAsPNG(file.getAbsolutePath());
		}
		
		return ret;
	}
}
