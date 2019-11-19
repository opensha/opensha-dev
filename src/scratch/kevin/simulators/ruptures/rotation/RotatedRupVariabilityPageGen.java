package scratch.kevin.simulators.ruptures.rotation;

import java.awt.Color;
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
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jfree.data.Range;
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
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
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
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
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

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;

import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.ASK_EventData;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen.ValidationResult;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public abstract class RotatedRupVariabilityPageGen {
	
	private static final boolean D = false;

	private RSQSimCatalog catalog;
	private Map<Double, RotatedRupVariabilityConfig> magConfigs;
	private Map<Double, SimulationRotDProvider<RotationSpec>> magProvs;
	
	private double[] calcPeriods;
	
	private List<Site> sites;
	private Map<Float, List<Site>> vs30SiteBundles;
	private List<Float> sourceAzimuths;
	private List<Float> siteSourceAzimuths;
	private List<Float> distances;
	private Map<Double, List<Integer>> magEventIDs;
	private int maxNumEvents = 0;
	private int minNumEvents = Integer.MAX_VALUE;
	private Table<Double, Quantity, List<?>> magQuantitiesTable;
	
	private boolean replotAzimuthDependence = false;
	
	private Map<Integer, RSQSimEvent> eventsMap;
	private LoadingCache<RSQSimEvent, RSQSimSubSectEqkRupture> gmpeEventCache;
	
	private NGAW2_WrapperFullParam[] gmpes;
	
	private int numExampleRotations = 5;
	
	private Color gmpeColor;
	private List<Color> siteColors;
	
	private LoadingCache<VarGroupingKey, VariabilityResult> varResultCache;
	private LoadingCache<VarGroupingKey, VariabilityResult[]> downsampledVarResultCache;
	private LoadingCache<GMPE_GroupingKey, GMPE_Result> gmpeResultCache;
	
	private Map<Integer, List<ASK_EventData>> realEventData;
	private int numRealDataSamples;
	
	private boolean hasMagDist;

	public RotatedRupVariabilityPageGen(RSQSimCatalog catalog, RotatedRupVariabilityConfig config,
			double mag, SimulationRotDProvider<RotationSpec> prov, double[] calcPeriods) {
		this(catalog, emptyMagMap(mag, config), emptyMagMap(mag, prov), calcPeriods);
	}
	
	private static <T> HashMap<Double, T> emptyMagMap(double mag, T value) {
		HashMap<Double, T> map = new HashMap<>();
		map.put(mag, value);
		return map;
	}

	public RotatedRupVariabilityPageGen(RSQSimCatalog catalog, Map<Double, RotatedRupVariabilityConfig> magConfigs,
			Map<Double, SimulationRotDProvider<RotationSpec>> magProvs, double[] calcPeriods) {
		this.catalog = catalog;
		this.magConfigs = magConfigs;
		this.magProvs = magProvs;
		this.calcPeriods = calcPeriods;
		
		RotatedRupVariabilityConfig config0 = magConfigs.values().iterator().next();
		
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
			RotatedRupVariabilityConfig config = magConfigs.get(mag);
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
		
		varResultCache = CacheBuilder.newBuilder().maximumSize(500).build(new CacheLoader<VarGroupingKey, VariabilityResult>() {

			@Override
			public VariabilityResult load(VarGroupingKey key) throws Exception {
				return calcVarResult(key);
			}
			
		});
		
		downsampledVarResultCache = CacheBuilder.newBuilder().maximumSize(5).build(new CacheLoader<VarGroupingKey, VariabilityResult[]>() {

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
		
		gmpeEventCache = CacheBuilder.newBuilder().build(new CacheLoader<RSQSimEvent, RSQSimSubSectEqkRupture>() {

			@Override
			public RSQSimSubSectEqkRupture load(RSQSimEvent key) throws Exception {
				return buildGMPE_Rupture(key);
			}
			
		});
		
		hasMagDist = magConfigs.size() > 3 && distances.size() > 3;
	}
	
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
	
	protected void setEventsMap(Map<Integer, RSQSimEvent> eventsMap) {
		this.eventsMap = eventsMap;
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
	
	private synchronized RSQSimEvent getEvent(int eventID) {
		if (eventsMap == null) {
			try {
				eventsMap = loadEvents(catalog, getAllEventIDs());
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		return eventsMap.get(eventID);
	}
	
	private RSQSimSubSectEqkRupture buildGMPE_Rupture(RSQSimEvent event) {
		return catalog.getMappedSubSectRupture(event);
	}
	
	protected Scenario getBBP_PartB_Scenario(RotatedRupVariabilityConfig config) {
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
		TOTAL("σ", "&sigma;") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return gm.stdDev();
			}
		},
		PHI("φ", "&phi;") {
			@Override
			public double calculate(ScalarGroundMotion gm) {
				return gm.phi();
			}
		},
		PHI_SS("φ_ss", "&phi;<sub>SS</sub>") {
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
		PATH("Path-to-path", "path", "φ_p2p", "&phi;<sub>P2P</sub>", al_atik,
				null, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH), // group quantities
				qarr(Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false), // std dev of residuals
		SOUCE_STRIKE("Source-strike", "source_strike", "φ_s", "&phi;<sub>s</sub>", aki_richards,
//				GMPE_Var.PHI_SS, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				GMPE_Var.PHI_SS, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
//				Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH),
		WITHIN_EVENT_SS("Within-event, single-site", "within_event_ss", "φ_ss", "&phi;<sub>SS</sub>", al_atik,
//				GMPE_Var.PHI_SS, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				GMPE_Var.PHI_SS, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				qarr(), // singleton quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
		WITHIN_EVENT("Within-event", "within_event", "φ", "&phi;", al_atik,
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
//		SITE_TO_SITE("Site Variability", "site_var", "φ_s2s", "&phi;<sub>S2S</sub>", al_atik,
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
		
		public List<String> buildMethodologyLines(RotatedRupVariabilityPageGen pageGen, File resourcesDir) throws IOException {
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
	
	private Map<IDPairing, Double> elemDistCache = new HashMap<>();
	private double calcVprop(RSQSimEvent event) {
		double minTime = Double.POSITIVE_INFINITY;
		SimulatorElement hypo = null;
		for (EventRecord rec : event) {
			double[] times = rec.getElementTimeFirstSlips();
			Preconditions.checkNotNull(times, "Event doesn't have timing information");
			List<SimulatorElement> elems = rec.getElements();
			
			for (int i=0; i<elems.size(); i++) {
				if (times[i] < minTime) {
					minTime = times[i];
					hypo = elems.get(i);
				}
			}
		}

		int hypoID = hypo.getID();
		List<Double> vels = new ArrayList<>();
		for (EventRecord rec : event) {
			double[] times = rec.getElementTimeFirstSlips();
			Preconditions.checkNotNull(times, "Event doesn't have timing information");
			List<SimulatorElement> elems = rec.getElements();
			
			for (int i=0; i<elems.size(); i++) {
				SimulatorElement elem = elems.get(i);
				if (elem.getID() == hypoID)
					continue;
				int elemID = elem.getID();
				IDPairing pair = hypoID > elemID ? new IDPairing(elemID, hypoID) : new IDPairing(hypoID, elemID);
				Double dist = elemDistCache.get(pair);
				if (dist == null) {
					dist = LocationUtils.linearDistanceFast(hypo.getCenterLocation(), elem.getCenterLocation());
					elemDistCache.put(pair, dist);
				}
				double tDelta = times[i] - minTime;
				if (tDelta == 0)
					continue;
				double vel = dist/(tDelta);
				vels.add(vel);
			}
			
		}
		
		return DataUtils.median(Doubles.toArray(vels));
	}
	
	private enum ScatterDisaggQuantity {
		V_PROP("Vprop", "v_prop") {
			@Override
			public double getValue(RotatedRupVariabilityPageGen pageGen, RotationSpec rotation) {
				return pageGen.calcVprop(pageGen.getEvent(rotation.eventID));
			}
		},
		SOURCE_AZ("Source Azimuth", "src_az") {
			@Override
			public double getValue(RotatedRupVariabilityPageGen pageGen, RotationSpec rotation) {
				return rotation.sourceAz == null ? 0d : rotation.sourceAz;
			}
		},
		SITE_TO_SOURCE_AZ("Site-to-Source Az", "site_source_az") {
			@Override
			public double getValue(RotatedRupVariabilityPageGen pageGen, RotationSpec rotation) {
				return rotation.siteToSourceAz == null ? 0d : rotation.siteToSourceAz;
			}
		};
		
		private String name;
		private String prefix;
		
		private ScatterDisaggQuantity(String name, String prefix) {
			this.name = name;
			this.prefix = prefix;
		}
		
		public abstract double getValue(RotatedRupVariabilityPageGen pageGen, RotationSpec rotation);
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
	
	public void generatePage(File outputDir, double[] periods, List<String> methodSpecificLines,
			double[] highlightMags, float[] highlightDists) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new LinkedList<>();
		
		String distName = BBP_PartBValidationConfig.DIST_JB ? "Joyner-Boore distance" : "3-dimensional distance";
		String distSymbol = BBP_PartBValidationConfig.DIST_JB ? "Rjb" : "Rrup";
		
		lines.add("# "+catalog.getName()+" Rotated Rupture Variability, "+getScenarioShortName());
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
		lines.add("## "+getScenarioShortName()+" RSQSim Rupture Match Criteria");
		lines.add(topLink); lines.add("");
		String[] criteria = getScenarioMatchCriteria();
		lines.add("We condisder "+maxNumEvents+" events in the catalog which match the following criteria:");
		lines.add("");
		for (String criterion : criteria)
			lines.add("* "+criterion);
		lines.add("");
		
		if (magConfigs.size() == 1) {
			lines.add("### Fault Section Counts");
			lines.add(topLink); lines.add("");
			lines.add("This tables gives a list of all fault sections which participate in the ruptures matching the above "
					+ "criteria. Many ruptures include multiple sections, so the sum of counts may be larger than the number "
					+ "of events.");
			lines.add("");
			Map<String, Integer> parentCountsMap = new HashMap<>();
			for (int eventID : magConfigs.values().iterator().next().getValues(Integer.class, Quantity.EVENT_ID)) {
				RSQSimEvent event = getEvent(eventID);
				RSQSimSubSectEqkRupture rup;
				try {
					rup = gmpeEventCache.get(event);
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
			}
			List<String> parentNames = ComparablePairing.getSortedData(parentCountsMap);
			Collections.reverse(parentNames);
			
			table = MarkdownUtils.tableBuilder();
			table.addLine("Section Name", "Participation Count");
			for (String parentName : parentNames)
				table.addLine(parentName, parentCountsMap.get(parentName));
			lines.addAll(table.build());
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
				RotatedRupVariabilityConfig config = magConfigs.get(mag);
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
					+ "validation exercise: Methodology for code validation in the context of seismic‐hazard analyses. "
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
				RotatedRupVariabilityConfig config = magConfigs.get(mag);
				
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
					File file = writeCSV(resourcesDir, site, mag, dist);
					table.addLine("M"+optionalDigitDF.format(mag), dist+" km", site.getName(),
							"["+file.getName()+"](resources/"+file.getName()+")");
				}
			}
		}
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
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
	
	private List<String> generateVariabilityLines(RotatedRupVariabilityConfig config, VariabilityType type, Double mag, Float distance,
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
		List<String> siteBundlePlotLables = new ArrayList<>();
		List<List<Site>> siteBundles = new ArrayList<>();
		List<PlotCurveCharacterstics> bundleCurveChars = new ArrayList<>();
		if (sites.size() > 1) {
			if (type.separateSites) {
				siteBundleTableNames.add("**ALL SITES**");
				siteBundlePlotLables.add("Total");
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
					siteBundlePlotLables.add("Vs30="+optionalDigitDF.format(vs30));
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
					siteBundlePlotLables.add(site.name);
					List<Site> siteList = new ArrayList<>();
					siteList.add(site);
					siteBundles.add(siteList);
					bundleCurveChars.add(new PlotCurveCharacterstics(
							PlotLineType.SOLID, 1.5f, siteColors.get(i % siteColors.size())));
				}
			}
		} else {
			siteBundleTableNames.add(null);
			siteBundlePlotLables.add("Total");
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
			for (String siteName : siteBundleTableNames)
				header.add(siteName);
			eventTermsCSV.addLine(header);
			
			System.out.println("writing debug table for "+type);
			for (int i=0; i<siteBundles.size(); i++) {
				Site[] mySites = siteBundles.get(i).toArray(new Site[0]);
				VarGroupingKey key = new VarGroupingKey(type, mag, distance, mySites);
				PeriodDepResidualsList periodResiduals = loadResiduals(key);
				List<ResidualSet> residuals = periodResiduals.get(0);
				System.out.println("Calculating for "+siteBundlePlotLables.get(i));
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
				System.out.println("Calculating for "+siteBundlePlotLables.get(i));
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
					
					dataSizeDependenceLines = plotDataSizeDependence(key, resourcesDir, stdDevFunc, gmpeResult);
				}
			}
			stdDevFunc.setName(siteBundlePlotLables.get(i));
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
		if (downsampledList != null && bundlePercentilesList != null && bundlePercentilesList.get(0) != null) {
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
					+ "from these samples is plotted as a shaded region above, and listed in the table below.");
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
			if (dataSizeDependenceLines != null) {
				lines.addAll(dataSizeDependenceLines);
				lines.add("");
			}
			// now plot site histograms
			if (type.separateSites) {
				lines.add("These plots show the distribution of period-independent downsampled "+type.htmlSymbol
						+" for each site.");
				lines.add("");
				table = MarkdownUtils.tableBuilder();
				if (sites.size() > 1 || sites.get(0) != null) {
					table.initNewLine();
					for (Site site : sites)
						table.addColumn("**"+site.getName()+"**");
					table.finalizeLine();
				}
				table.initNewLine();
				for (Site site : sites) {
					String sitePrefix = site == null ? "" : "_"+site.getName();
					prefix = type.prefix+magPrefix+distPrefix+sitePrefix+"_downsampled_hist";
					plotSitePeriodIndepDownsampledHistogram(resourcesDir, prefix, type, mag, distance, site);
					table.addColumn("![Dowmsampled Histogram](resources/"+prefix+".png)");
				}
				table.finalizeLine();
				lines.addAll(table.wrap(4, 0).build());
				lines.add("");
			} else {
				lines.add("This plot shows the distribution of period-independent downsampled "+type.htmlSymbol+".");
				lines.add("");
				prefix = type.prefix+magPrefix+distPrefix+"_downsampled_hist";
				plotSitePeriodIndepDownsampledHistogram(resourcesDir, prefix, type, mag, distance,
						siteBundles.get(0).toArray(new Site[0]));
				lines.add("![Dowmsampled Histogram](resources/"+prefix+".png)");
				lines.add("");
			}
		}
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
		if (type.disaggScatters != null && type.disaggScatters.length > 0 && (distance != null || !type.separateDists)) {
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
		
		if (dependUseSimCounts) {
			double dm = 0.2;
			// we want actual num recordings and num events for all matching ruptures here,
			// regardless of distance
			Map<Integer, List<ASK_EventData>> dataMatches = ASK_EventData.getMatches(realEventData,
					com.google.common.collect.Range.closed(key.magnitude-dm, key.magnitude+dm),
					null, null, 0d);
			realData = new ArrayList<>();
			for (List<ASK_EventData> values : dataMatches.values())
				realData.add(values);
		} else {
			realData = getRealDataListForKey(key);
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
		
		RotatedRupVariabilityConfig config = magConfigs.get(key.magnitude);
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
				plotDataSizeDependence(numEventsFunc, numEventsPercentiles, type.name+" ("+type.symbol+") Event Count Dependence",
						"Num Events", type.symbol, numEvents, "ASK (2014) Num Events", simVal, gmpeVal, resourcesDir, prefix, type);
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
				plotDataSizeDependence(numRecordingsFunc, numRecordingsPercentiles, type.name+" ("+type.symbol+") Event Recordings Dependence",
						"Num Recordings Per Event", type.symbol, aveNumRecordings, "ASK (2014) Num Recs/Event", simVal, gmpeVal, resourcesDir, prefix, type);
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
			String title, String xAxisLabel, String yAxisLabel, double dataX, String dataName,
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
		percentiles.bounds68.setName("68% Range");
		funcs.add(percentiles.bounds68);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, color68));
		medianFunc.setName("Median");
		funcs.add(medianFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		XY_DataSet dataFunc = new DefaultXY_DataSet();
		dataFunc.setName(dataName);
		dataFunc.set(dataX, 0d);
		dataFunc.set(dataX, 1d);
		funcs.add(dataFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.RED));
		
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
	
	private static Float nullAsZero(Float value) {
		if (value == null)
			return 0f;
		return value;
	}
	
	private class ResidualSet {
		double[] residuals;
		double median;
		double residualStdDev;
		
		private RotationSpec commonRotSpec;
		
		public ResidualSet(List<RotationSpec> rotations, double[] values) {
			this.residuals = values; // will subract median in a sec
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
			return new VariabilityResult(null, null, residualStdDevs, medianStdDevs);
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
				
				if (key.realDataList != null)
					varResultCache.put(siteKey, new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs));
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
			
			return new VariabilityResult(commonRotSpecs, totLogVals, totResidualStdDevs, totMedianStdDevs);
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
		return new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs);
	}
	
	private List<List<ASK_EventData>> getRealDataListForKey(VarGroupingKey key) {
		double dm = 0.2;
		double dr = key.distance == null || key.distance > 80 ? 20 : 10;
		Map<Integer, List<ASK_EventData>> dataMatches = ASK_EventData.getMatches(realEventData,
				com.google.common.collect.Range.closed(key.magnitude-dm, key.magnitude+dm),
				key.distance == null ? null : com.google.common.collect.Range.closed(key.distance-dr, key.distance+dr), null, 0d);
		Preconditions.checkNotNull(dataMatches);
		if (key.distance == null && key.type.separateDists && distances.size() > 1) {
			// we're calculating across all of our distances, to compare to distance independent recordings
			// instead partition the number of recordings equally across each distance
			Map<Integer, List<ASK_EventData>> filteredMatches = new HashMap<>();
			for (int id : dataMatches.keySet()) {
				List<ASK_EventData> recordings = dataMatches.get(id);
				int targetNum = (int)Math.round((double)recordings.size()/(double)distances.size());
				if (targetNum < 2 && recordings.size() >= 2)
					targetNum = 2;
				if (D) System.out.println("Target num="+targetNum+" for "+recordings.size()+" recordings and "+distances.size()+" dists");
				if (targetNum > 1)
					filteredMatches.put(id, recordings.subList(0, targetNum));
			}
			dataMatches = filteredMatches;
		}
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
	
	private class VariabilityResult {
		private final RotationSpec[] commonRotationSpecs;
		private final LogValueSet[] logVals;
		private final ResidualStdDevSet[] residualStdDevs;
		private final ResidualStdDevSet periodIndepResidualStdDev;
		private final MedianStdDevSet[] medianStdDevs;
		private final MedianStdDevSet periodIndepMedianStdDev;
		
		public VariabilityResult(RotationSpec[] commonRotationSpecs, LogValueSet[] logVals,
				ResidualStdDevSet[] residualStdDevs, MedianStdDevSet[] medianStdDevs) {
			this.commonRotationSpecs = commonRotationSpecs;
			if (hasMagDist)
				this.logVals = logVals;
			else
				// only need them for mag-dist
				this.logVals = null;
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
	
	private class LogValueSet {
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
	
	private class ResidualStdDevSet {
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
			StandardDeviation totalCalc = new StandardDeviation();
			for (int i=0; i<stdDevs.length; i++) {
				ResidualSet set = residuals.get(i);
				stdDevs[i] = set.getStdDev();
				totalCalc.incrementAll(set.getResiduals());
			}
			total = totalCalc.getResult();
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
	
	private class MedianStdDevSet {
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
			for (int i=0; i<medians.length; i++) {
				ResidualSet set = residuals.get(i);
				medians[i] = set.median;
			}
			stdDev = Math.sqrt(StatUtils.variance(medians));
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
	
	private Random repeatableRandom;
	
	private PeriodDepResidualsList calcResiduals(Double magnitude, Quantity[] constQuantities, Object[] constValues,
			Quantity[] groupQuantities, Quantity[] singletons, List<List<ASK_EventData>> realDataList)
					throws IOException {
//		if (D) System.out.println("Rotations for mag: "+magnitude);
		RotatedRupVariabilityConfig config = magConfigs.get(magnitude);
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
				repeatableRandom = new Random(eventsMap.size());
			
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
			int eventsToInclude = Integer.min(eventIDs.size(), realDataList.size());
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
			RSQSimEvent event = getEvent(eventID);
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
	
	private void plotSitePeriodIndepDownsampledHistogram(File resourcesDir, String prefix,
			VariabilityType type, Double magnitude, Float distance, Site... sites) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
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
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 1d, delta);
		hist.setName("Downsampled Realizations");
		
		double meanDownsampled = 0d;
		for (VariabilityResult result : dsResults) {
			double stdDev = type.stdDevOfMedians ? result.getPeriodIndepMedianStdDevSet().stdDev
					: result.getPeriodIndepResidualStdDevSet().total;
			meanDownsampled += stdDev;
			hist.add(hist.getClosestXIndex(stdDev), 1d);
		}
		hist.normalizeBySumOfY_Vals();
		meanDownsampled /= (double)dsResults.length;
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
		
		funcs.add(line(meanDownsampled, 0d, meanDownsampled, 1d, "Mean of Downsampled"));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Color.BLACK));
		
		double totStdDev = type.stdDevOfMedians ? totResult.getPeriodIndepMedianStdDevSet().stdDev
				: totResult.getPeriodIndepResidualStdDevSet().total;
		funcs.add(line(totStdDev, 0d, totStdDev, 1d, "Total"));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, type.symbol+" Downsampled Distribution", "Standard Deviation", "");
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
						double gmpeSD = gmpeVars == null ? -1 : gmpeVars[p];
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
		
		RotatedRupVariabilityConfig config = magConfigs.get(magnitude);
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
		
		RotatedRupVariabilityConfig config = magConfigs.get(magnitude);
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
					if (medians == null) {
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
						for (int p=0; p<periods.length; p++)
							gmpeXYs[d][p].set(value, gmpeVars[p]);
					}
				}
				// now do median
				DiscretizedFunc[][] xys = xysMap.get(null);
				for (int p=0; p<periods.length; p++)
					xys[d][p].set(value, medians[p]);
				if (gmpeResult != null) {
					DiscretizedFunc[][] gmpeXYs = gmpeXYsMap.get(null);
					for (int p=0; p<periods.length; p++)
						gmpeXYs[d][p].set(value, Math.exp(gmpeResult.logMedian[p]));
				}
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		double[] distArray = Doubles.toArray(distances);
		CPT distCPT = new CPT(StatUtils.min(distArray), StatUtils.max(distArray), Color.GRAY, Color.BLACK);
		
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
	
	private void plotExample(File resourcesDir, String prefix, double distance, List<Quantity> variedQuantities)
			throws IOException {
		List<Site> sites = new ArrayList<>();
		
		Double minMag = Double.POSITIVE_INFINITY;
		for (Double mag : magEventIDs.keySet())
			minMag = Double.min(minMag, mag);
		RSQSimEvent exampleRupture = getEvent(magEventIDs.get(minMag).get(0));
		sites.add(this.sites.get(0));
		List<RSQSimEvent> ruptures = new ArrayList<>();
		ruptures.add(exampleRupture);
		int numSourceAz = variedQuantities.contains(Quantity.SOURCE_AZIMUTH) ? numExampleRotations : 1;
		int numSiteToSourceAz = variedQuantities.contains(Quantity.SITE_TO_SOURTH_AZIMUTH) ? numExampleRotations : 1;
		RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(catalog, sites, ruptures, new double[] {distance},
				numSourceAz, numSiteToSourceAz);
		
		config.plotRotations(resourcesDir, prefix, config.getRotations(), true);
	}
	
	private File writeCSV(File resourcesDir, Site site, Double magnitude, Float distance) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = Lists.newArrayList(Quantity.EVENT_ID.getName(), Quantity.SOURCE_AZIMUTH.getName(),
				Quantity.SITE_TO_SOURTH_AZIMUTH.getName());
		for (double period : calcPeriods)
			header.add(optionalDigitDF.format(period)+"s SA");
		
		csv.addLine(header);
		
		SimulationRotDProvider<RotationSpec> prov = magProvs.get(magnitude);
		RotatedRupVariabilityConfig config = magConfigs.get(magnitude);
		
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
}
