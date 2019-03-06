package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
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
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
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
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
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
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationPageGen.ValidationResult;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public abstract class RotatedRupVariabilityPageGen {

	private RSQSimCatalog catalog;
	private Map<Double, RotatedRupVariabilityConfig> magConfigs;
	private Map<Double, SimulationRotDProvider<RotationSpec>> magProvs;
	
	private double[] calcPeriods;
	
	private List<Site> sites;
	private List<Float> sourceAzimuths;
	private List<Float> siteSourceAzimuths;
	private List<Float> distances;
	private Map<Double, List<Integer>> magEventIDs;
	private int maxNumEvents = 0;
	private int minNumEvents = Integer.MAX_VALUE;
	private Table<Double, Quantity, List<?>> magQuantitiesTable;
	
	private Map<Integer, RSQSimEvent> eventsMap;
	private LoadingCache<RSQSimEvent, EqkRupture> gmpeEventCache;
	
	private NGAW2_WrapperFullParam[] gmpes;
	
	private int numExampleRotations = 5;
	
	private List<Color> siteColors;
	
	private LoadingCache<VarGroupingKey, VariabilityResult> varResultCache;
	private LoadingCache<GMPE_GroupingKey, GMPE_Result> gmpeResultCache;

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
		if (sites.size() > 1) {
			try {
				siteCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, sites.size()-1);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			siteColors = new ArrayList<>();
			for (int i=0; i<sites.size(); i++)
				siteColors.add(siteCPT.getColor((float)i).darker());
		}
		
		varResultCache = CacheBuilder.newBuilder().build(new CacheLoader<VarGroupingKey, VariabilityResult>() {

			@Override
			public VariabilityResult load(VarGroupingKey key) throws Exception {
				return calcVarResult(key);
			}
			
		});
		
		gmpeResultCache = CacheBuilder.newBuilder().build(new CacheLoader<GMPE_GroupingKey, GMPE_Result>() {

			@Override
			public GMPE_Result load(GMPE_GroupingKey key) throws Exception {
				return calcGMPE(key);
			}
			
		});
		
		gmpeEventCache = CacheBuilder.newBuilder().build(new CacheLoader<RSQSimEvent, EqkRupture>() {

			@Override
			public EqkRupture load(RSQSimEvent key) throws Exception {
				return buildGMPE_Rupture(key);
			}
			
		});
	}
	
	static Map<Integer, RSQSimEvent> loadEvents(RSQSimCatalog catalog, Collection<Integer> ids) throws IOException {
		System.out.println("Loading "+ids.size()+" events...");
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		List<RSQSimEvent> events = catalog.loader().byIDs(Ints.toArray(ids));
		Preconditions.checkState(events.size() == ids.size(), "Loaded %s events, expected %s", events.size(), ids.size());
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		return eventsMap;
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
	
	private EqkRupture buildGMPE_Rupture(RSQSimEvent event) {
		return catalog.getGMPE_Rupture(event, RSQSimBBP_Config.MIN_SUB_SECT_FRACT);
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
		TOTAL("σ"),
		PHI("φ"),
		TAU("τ");
		
		String symbol;
		private GMPE_Var(String symbol) {
			this.symbol = symbol;
		}
	}
	
	private static Quantity[] qarr(Quantity... quantities) {
		return quantities;
	}
	
	private enum VariabilityType {
		PATH("Path-to-path", "path", "φ_p2p", "&phi;<sub>P2P</sub>", al_atik,
				null, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH), // group quantities
				qarr(Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				false), // std dev of residuals
		SOUCE_STRIKE("Source-strike", "source_strike", "φ_s", "&phi;<sub>s</sub>", aki_richards,
				GMPE_Var.PHI, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH), // vary quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
//				Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH),
		WITHIN_EVENTS("Within-event, single-site", "within_event_ss", "φ_ss", "&phi;<sub>SS</sub>", al_atik,
				GMPE_Var.PHI, new ScatterDisaggQuantity[] {ScatterDisaggQuantity.V_PROP},
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
				false, "δw_es", "&delta;W<sub>es</sub>"), // std dev of residuals
//		BETWEEN_EVENTS_SINGLE_PATH("Between-events, single-path", "between_events_single_path", "τ_0", "&tau;<sub>0</sub>", al_atik,
//				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH),
		BETWEEN_EVENTS("Between-events", "between_events", "τ", "&tau;", al_atik,
				GMPE_Var.TAU, null,
				qarr(Quantity.SITE, Quantity.DISTANCE), // separate quantities
				qarr(Quantity.EVENT_ID), // group quantities
				qarr(Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH), // vary quantities
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
		private boolean stdDevOfMedians; // otherwise std of residuals
		private String variedSymbol;
		private String variedSymbolHTML;
		
		private VariabilityType(String name, String prefix, String symbol, String htmlSymbol, String reference,
				GMPE_Var gmpeStdDevType, ScatterDisaggQuantity[] disaggScatters,
				Quantity[] separateQuantities, // compute sigma independently for all of these
				Quantity[] groupQuantities, // compute sigma of values for each unique combination of these
				Quantity[] variedQuantities,// across all unique combinations of these
				boolean stdDevOfMedians) { 
			this(name, prefix, symbol, htmlSymbol, reference, gmpeStdDevType, disaggScatters,
					separateQuantities, groupQuantities, variedQuantities, stdDevOfMedians, null, null);
		}

		private VariabilityType(String name, String prefix, String symbol, String htmlSymbol, String reference,
				GMPE_Var gmpeStdDevType, ScatterDisaggQuantity[] disaggScatters,
				Quantity[] separateQuantities, // compute sigma independently for all of these
				Quantity[] groupQuantities, // compute sigma of values for each unique combination of these
				Quantity[] variedQuantities, // across all unique combinations of these
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
					lines.add("* "+quantity+" *["+num+" unique]*");
				}
				lines.add("");
				if (variedSymbolHTML == null)
					line = "We take "+htmlSymbol+" to be the standard deviation of those medians.";
				else
					line = "We take "+htmlSymbol+" to be the standard deviation of all "+variedSymbolHTML+".";
				if (separateSites && pageGen.sites.size() > 1)
					line += " Finally, we compute the median standard deviation across all sites. This"
							+ " total value is reported as **ALL SITES** and in summary plots/tables.";
				lines.add(line);
				lines.add("");
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
					lines.add("* "+quantity+" *["+num+" unique]*");
				}
				lines.add("");
				line = "We take "+htmlSymbol+" to be the standard deviation of all residuals"+symbolAdd+" across each combination of "
						+Joiner.on(", ").join(groupQuantities)+".";
				if (separateSites && pageGen.sites.size() > 1)
					line += " We also compute the total standard deviation across all residuals from all sites. This"
							+ " total value is reported as **ALL SITES** and in summary plots/tables.";
				lines.add(line);
				lines.add("");
			}
			
			if (pageGen.numExampleRotations > 0) {
				boolean hasEvent = false;
				for (Quantity quantity : groupQuantities)
					if (quantity == Quantity.EVENT_ID)
						hasEvent = true;
				if (separateSites && hasEvent) {
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
		
		boolean hasMagDist = magConfigs.size() > 3 && distances.size() > 3;
		Map<VariabilityType, File[]> magDistPlots = new HashMap<>();
		
		HashSet<VariabilityType> computedTypes = new HashSet<>();
		
		for (VariabilityType type : VariabilityType.values()) {
			if (type.getRotationsPerStdDev(magQuantitiesTable.rowMap().values().iterator().next()) == 1)
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
				} else {
					lines.addAll(generateVariabilityLines(config, type, mag, null, periods, topLink, resourcesDir));
					lines.add("");
				}
			}
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
			for (double period : periods)
				table.addColumn(optionalDigitDF.format(period)+"s Std. Dev.");
			table.finalizeLine();
			for (VariabilityType type : VariabilityType.values()) {
				if (!computedTypes.contains(type))
					continue;
				List<Float> myDists;
				if (type.separateDists) {
					myDists = this.distances;
				} else {
					myDists = new ArrayList<>();
					myDists.add(null);
				}
				for (Float distance : myDists) {
					table.initNewLine();
					VarGroupingKey key = new VarGroupingKey(type, mag, distance, null);
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
						for (double period : periods)
							table.addColumn(optionalDigitDF.format(varResult.getMedianStdDevSet(period).stdDev));
					} else {
						for (double period : periods)
							table.addColumn(optionalDigitDF.format(varResult.getResidualStdDevSet(period).total));
					}
					table.finalizeLine();
				}
			}
			summaryLines.addAll(table.build());
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
									resourcesDir, prefix+"_2d", quantity, computedTypes, site, magnitude, periods);
							plotMaps2D.put(magnitude, plotsMap);
						}
						
						Map<VariabilityType, File[]> plotsMap = plotQuantityDistFuncs(
									resourcesDir, prefix, quantity, computedTypes, site, magnitude, periods, plotDists);
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
				
				for (Site site : sites) {
					summaryTable.initNewLine();
					plotTable.initNewLine();
					summaryTable.addColumn("**"+scenario.getShortName()+"**");
					if (sites.size() > 1) {
						summaryTable.addColumn("**"+site.getName()+"**");
						plotTable.addColumn("**"+site.getName()+"**");
					}
					
					double vs30 = site.getParameter(Double.class, Vs30_Param.NAME).getValue();
					
					for (Float distance : distances) {
						boolean official = scenario.isOfficialCriteria() && isOfficialPartB_Distance(distance);
						
						List<DiscretizedFunc> rd50s = new ArrayList<>();
						for (RotationSpec rotation : config.getRotationsForQuantities(Quantity.SITE, site, Quantity.DISTANCE, distance)) {
							Float sourceAz = (Float)rotation.getValue(Quantity.SOURCE_AZIMUTH);
							if (footwall && sourceAz < 180f)
								// only footwall
								continue;
							rd50s.addAll(simProv.getRotD50s(site, rotation));
						}
						
						String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
						
						// plot it
						String prefix = "bbp_partB_"+scenario.getPrefix()+distPrefix+(sites.size() > 1 ? "_"+site.getName() : "");
//						
						List<ValidationResult> results = BBP_PartBValidationPageGen.calcPlotScenarioResults(
								scenario, distance, vs30, numEvents, rd50s, resourcesDir, prefix);
						
						boolean passes = true;
						for (ValidationResult result : results)
							passes = passes && result.passes();
						
						summaryTable.addColumn(bbpTableStr(passes, official));
						plotTable.addColumn("![PartB Plot](resources/"+prefix+".png)");
					}
					
					summaryTable.finalizeLine();
					plotTable.finalizeLine();
				}
				
				lines.addAll(plotTable.build());
				lines.add("");
			}
			
			lines.addAll(summaryIndex, summaryTable.build());
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
		String curHeading = "##";
		String uniqueTitleSection = mag == null ? type.name : "M"+optionalDigitDF.format(mag)+" "+type.name;
		if (type.separateDists && distances.size() > 1) {
			curHeading += "#";
			uniqueTitleSection = distance+" km "+uniqueTitleSection;
			lines.add(curHeading+" "+uniqueTitleSection+" Results");
			lines.add(topLink); lines.add("");
		}
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		if (type.separateSites && sites.size() > 1)
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
		DiscretizedFunc totStdDevFunc = null;
		ResidualStdDevSet[] totResidualStdDevs = new ResidualStdDevSet[periods.length];
		List<List<ResidualStdDevSet>> siteTotResidualStdDevs = null;
		List<List<MedianStdDevSet>> siteTotMedianStdDevs = null;
		if (type.separateSites && sites.size() > 1) {
			siteTotResidualStdDevs = new ArrayList<>();
			siteTotMedianStdDevs = new ArrayList<>();
			for (int i=0; i<periods.length; i++) {
				siteTotResidualStdDevs.add(new ArrayList<>());
				siteTotMedianStdDevs.add(new ArrayList<>());
			}
		}
		
		List<Site> mySites = new ArrayList<>();
		List<DiscretizedFunc> stdDevFuncs = null;
		if (type.separateSites && sites.size() > 1) {
			mySites.addAll(sites);
			stdDevFuncs = new ArrayList<>();
		}
		mySites.add(null);
		
		for (int i=0; i<mySites.size(); i++) {
			Site site = mySites.get(i);
			table.initNewLine();
			String emphasis = "";
			if (type.separateSites && sites.size() > 1) {
				if (site == null) {
					table.addColumn("**ALL SITES**");
					emphasis = "**";
				} else {
					table.addColumn(site.getName());
				}
			}
			DiscretizedFunc stdDevFunc = new ArbitrarilyDiscretizedFunc();
			if (site == null) {
				stdDevFunc.setName("Total");
				totStdDevFunc = stdDevFunc;
			} else {
				stdDevFunc.setName(site.getName());
				stdDevFuncs.add(stdDevFunc);
			}
			VariabilityResult result; 
			try {
				result = varResultCache.get(new VarGroupingKey(type, mag, distance, site));
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
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
			if (site == null) {
				if (!type.stdDevOfMedians) {
					for (int p=0; p<periods.length; p++)
						totResidualStdDevs[p] = result.getResidualStdDevSet(periods[p]);
				}
			} else {
				for (int p=0; p<periods.length; p++) {
					siteTotResidualStdDevs.get(p).add(result.getResidualStdDevSet(periods[p]));
					siteTotMedianStdDevs.get(p).add(result.getMedianStdDevSet(periods[p]));
				}
			}
			table.finalizeLine();
		}

		String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
		String magPrefix = mag == null ? "" : "_m"+optionalDigitDF.format(mag);
		
		// plot it
		String prefix = type.prefix+magPrefix+distPrefix+"_std_dev";
		plotPeriodDependentStdDevs(resourcesDir, prefix, type.name+" ("+type.symbol+")", stdDevFuncs, totStdDevFunc);
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
						periods[p], totResidualStdDevs[p], type.separateSites && sites.size() > 1 ? siteTotResidualStdDevs.get(p): null);
				table.addColumn("!["+optionalDigitDF.format(periods[p])+"s](resources/"+prefix+".png)");
			}
			table.finalizeLine();
			table.wrap(periods.length > 4 ? 3 : 2, 0);
			lines.addAll(table.build());
			lines.add("");
		}
		if (type.disaggScatters != null && type.disaggScatters.length > 0) {
			System.out.println("*** Plotting disagg scatters ***");
			Site site = sites.size() > 1 ? null : sites.get(0);
			VariabilityResult result;
			try {
				result = varResultCache.get(new VarGroupingKey(type, mag, distance, site));
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			for (ScatterDisaggQuantity scatterQuantity : type.disaggScatters) {
				System.out.println("*** Plotting disagg scatter of "+scatterQuantity+" ***");
				File[][] plots = plotScatter(resourcesDir, type.prefix+"_scatter_", periods, type, result, scatterQuantity);
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				for (double period : periods)
					table.addColumn(optionalDigitDF.format(period)+"s");
				table.finalizeLine();
				for (File[] periodPlots : plots) {
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
	
	private List<List<ResidualSet>> loadResiduals(VarGroupingKey key) throws IOException {
		Quantity[] constQuantities = key.separateQuantities;
		Object[] constValues = new Object[constQuantities.length];
		for (int i=0; i<constQuantities.length; i++) {
			if (constQuantities[i] == Quantity.SITE) {
				if (key.site == null && sites.size() == 1)
					constValues[i] = sites.get(0);
				else
					constValues[i] = key.site;
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
		
		return calcResiduals(key.magnitude, constQuantities, constValues, key.groupQuantities);
	}
	
	private VariabilityResult calcVarResult(VarGroupingKey key) throws IOException {
		if (key.separateSites && sites.size() > 1) {
			// do all sites at once
			List<List<ResidualSet>> totalPeriodResiduals = new ArrayList<>();
			for (int p=0; p<calcPeriods.length; p++)
				totalPeriodResiduals.add(new ArrayList<>());
			
			List<MedianStdDevSet[]> sitePeriodMediansList = new ArrayList<>();
			for (int p=0; p<calcPeriods.length; p++)
				sitePeriodMediansList.add(new MedianStdDevSet[sites.size()]);
			
			for (int s=0; s<sites.size(); s++) {
				Site site = sites.get(s);
				VarGroupingKey siteKey = new VarGroupingKey(
						key.separateQuantities, key.groupQuantities, key.magnitude, key.distance, site);
				siteKey.restrictTo = key.restrictTo;
				
				System.out.println("Loading residuals for site "+site.getName()+", dist="+key.distance);
				List<List<ResidualSet>> residuals = loadResiduals(siteKey);

				LogValueSet[] logVals = new LogValueSet[calcPeriods.length];
				ResidualStdDevSet[] residualStdDevs = new ResidualStdDevSet[calcPeriods.length];
				MedianStdDevSet[] medianStdDevs = new MedianStdDevSet[calcPeriods.length];
				RotationSpec[] commonRotSpecs = null;
				for (int p=0; p<calcPeriods.length; p++) {
					List<ResidualSet> periodResiduals = residuals.get(p);
					logVals[p] = new LogValueSet(periodResiduals);
					residualStdDevs[p] = new ResidualStdDevSet(periodResiduals);
					medianStdDevs[p] = new MedianStdDevSet(periodResiduals);
					if (p == 0) {
						commonRotSpecs = new RotationSpec[periodResiduals.size()];
						for (int i=0; i<periodResiduals.size(); i++)
							commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
					}
					totalPeriodResiduals.get(p).addAll(residuals.get(p));
					sitePeriodMediansList.get(p)[s] = medianStdDevs[p];
				}
				
				varResultCache.put(siteKey, new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs));
			}

			// now combine sites
			System.out.println("Combining residuals for "+sites.size()+" sites, dist="+key.distance);
			LogValueSet[] totLogVals = new LogValueSet[calcPeriods.length];
			ResidualStdDevSet[] totResidualStdDevs = new ResidualStdDevSet[calcPeriods.length];
			MedianStdDevSet[] totMedianStdDevs = new MedianStdDevSet[calcPeriods.length];
			RotationSpec[] commonRotSpecs = null;
			for (int p=0; p<calcPeriods.length; p++) {
				List<ResidualSet> periodResiduals = totalPeriodResiduals.get(p);
				totLogVals[p] = new LogValueSet(periodResiduals);
				totResidualStdDevs[p] = new ResidualStdDevSet(periodResiduals);
				if (p == 0) {
					commonRotSpecs = new RotationSpec[periodResiduals.size()];
					for (int i=0; i<periodResiduals.size(); i++)
						commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
				}
				totMedianStdDevs[p] = new MedianStdDevSet(sitePeriodMediansList.get(p));
			}
			VarGroupingKey noSiteKey = new VarGroupingKey(
					key.separateQuantities, key.groupQuantities, key.magnitude, key.distance, null);
			noSiteKey.restrictTo = key.restrictTo;
			varResultCache.put(noSiteKey, new VariabilityResult(commonRotSpecs, totLogVals, totResidualStdDevs, totMedianStdDevs));
			
			VariabilityResult ret = varResultCache.getIfPresent(key);
			Preconditions.checkState(ret != null, "Should have just been cached?");
			return ret;
		}
		// do for all sites
		List<List<ResidualSet>> residuals = loadResiduals(key);

		System.out.print("Combining residuals for "+(key.site == null ? "all sites" : key.site.getName())
				+", mag="+key.magnitude+", dist="+key.distance+"...");
		LogValueSet[] logVals = new LogValueSet[calcPeriods.length];
		ResidualStdDevSet[] residualStdDevs = new ResidualStdDevSet[calcPeriods.length];
		MedianStdDevSet[] medianStdDevs = new MedianStdDevSet[calcPeriods.length];
		RotationSpec[] commonRotSpecs = null;
		for (int p=0; p<calcPeriods.length; p++) {
			List<ResidualSet> periodResiduals = residuals.get(p);
			logVals[p] = new LogValueSet(periodResiduals);
			residualStdDevs[p] = new ResidualStdDevSet(periodResiduals);
			medianStdDevs[p] = new MedianStdDevSet(periodResiduals);
			if (p == 0) {
				commonRotSpecs = new RotationSpec[periodResiduals.size()];
				for (int i=0; i<periodResiduals.size(); i++)
					commonRotSpecs[i] = periodResiduals.get(i).commonRotSpec;
			}
		}
		System.out.println("DONE");
		return new VariabilityResult(commonRotSpecs, logVals, residualStdDevs, medianStdDevs);
	}
	
	private class VarGroupingKey {
		private final Quantity[] separateQuantities;
		private final Quantity[] groupQuantities;
		private final Double magnitude;
		private final Float distance;
		private final Site site;
		private Map<Quantity, Object> restrictTo;
		
		// derrivative quantities
		private boolean separateSites;
		
		public VarGroupingKey(VariabilityType type, Double magnitude, Float distance, Site site) {
			this(type.separateQuantities, type.groupQuantities, magnitude, distance, site);
		}
		public VarGroupingKey(Quantity[] separateQuantities, Quantity[] groupQuantities, Double magnitude,
				Float distance, Site site) {
			this.separateQuantities = separateQuantities;
			this.groupQuantities = groupQuantities;
			this.magnitude = magnitude;
			this.distance = distance;
			if (site == null && sites.size() == 1)
				site = sites.get(0);
			this.site = site;
			
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
			result = prime * result + ((distance == null) ? 0 : distance.hashCode());
			result = prime * result + Arrays.hashCode(groupQuantities);
			result = prime * result + ((magnitude == null) ? 0 : magnitude.hashCode());
			result = prime * result + Arrays.hashCode(separateQuantities);
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
			VarGroupingKey other = (VarGroupingKey) obj;
			if (!getOuterType().equals(other.getOuterType()))
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
	
	private class VariabilityResult {
		private final RotationSpec[] commonRotationSpecs;
		private final LogValueSet[] logVals;
		private final ResidualStdDevSet[] residualStdDevs;
		private final MedianStdDevSet[] medianStdDevs;
		
		public VariabilityResult(RotationSpec[] commonRotationSpecs, LogValueSet[] logVals,
				ResidualStdDevSet[] residualStdDevs, MedianStdDevSet[] medianStdDevs) {
			this.commonRotationSpecs = commonRotationSpecs;
			this.logVals = logVals;
			this.residualStdDevs = residualStdDevs;
			this.medianStdDevs = medianStdDevs;
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
	
	private class ResidualStdDevSet {
		private final double total;
		private final double mean;
		private final double median;
		private final double min;
		private final double max;
		private final double[] stdDevs;
		
		public ResidualStdDevSet(List<ResidualSet> residuals) {
			stdDevs = new double[residuals.size()];
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
		}
	}
	
	private class MedianStdDevSet {
		private final double stdDev;
		private final double[] medians;
		private final double meanMedian;
		private final double minMedian;
		private final double maxMedian;
		
		public MedianStdDevSet(MedianStdDevSet... otherMedians) {
			double[] stdDevs = new double[otherMedians.length];
			int totSize = 0;
			for (int i=0; i<otherMedians.length; i++) {
				MedianStdDevSet o = otherMedians[i];
				totSize += o.medians.length;
				stdDevs[i] = o.stdDev;
			}
			// median of medains
			stdDev = DataUtils.median(stdDevs);
			medians = new double[totSize];
			int ind = 0;
			for (int i=0; i<otherMedians.length; i++) {
				System.arraycopy(otherMedians[i].medians, 0, medians, ind, otherMedians[i].medians.length);
				ind += otherMedians[i].medians.length;
			}
			meanMedian = StatUtils.mean(medians);
			minMedian = StatUtils.min(medians);
			maxMedian = StatUtils.max(medians);
		}
		
		public MedianStdDevSet(List<ResidualSet> residuals) {
			medians = new double[residuals.size()];
			for (int i=0; i<medians.length; i++) {
				ResidualSet set = residuals.get(i);
				medians[i] = set.median;
			}
			stdDev = Math.sqrt(StatUtils.variance(medians));
			meanMedian = StatUtils.mean(medians);
			minMedian = StatUtils.min(medians);
			maxMedian = StatUtils.max(medians);
		}
	}
	
	private List<List<ResidualSet>> calcResiduals(Double magnitude,
			Quantity[] constQuantities, Object[] constValues, Quantity[] groupQuantities)
					throws IOException {
//		System.out.println("Rotations for mag: "+magnitude);
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
		int num = totalRotations.size();
		List<List<ResidualSet>> ret = new ArrayList<>();
		for (int i=0; i<calcPeriods.length; i++)
			ret.add(new ArrayList<>());
		
		System.out.println("Computing residuals/medians for Constants ["
				+(constQuantities == null ? "" : Joiner.on(",").join(constQuantities))
				+"], Groups ["+Joiner.on(",").join(groupQuantities)+"]: "+num+" simulations");
		
//		System.out.println("Computing resduals from "+num+" simulations");
		
		calcGroupedResiduals(magnitude, totalRotations, calcPeriods, ret, groupQuantities);
		int count = 0;
		for (ResidualSet set : ret.get(new Random().nextInt(calcPeriods.length)))
			count += set.size();
		Preconditions.checkState(count == num, "Bad end index. Expected %s, have %s", num, count);
//		System.out.println(numThatVary+"/"+totNum+" have path variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
		
		return ret;
	}
	
	private void calcGroupedResiduals(double magnitude, List<RotationSpec> rotations, double[] periods, List<List<ResidualSet>> ret,
			Quantity[] groupQuantities) throws IOException {
		Preconditions.checkState(!rotations.isEmpty());
		if (groupQuantities.length == 0) {
//			System.out.println("Computing with index="+index);
			// we're at a set for which to calculate residuals, do it
			List<DiscretizedFunc> spectra = new ArrayList<>();
			SimulationRotDProvider<RotationSpec> prov = magProvs.get(magnitude);
			for (RotationSpec rotation : rotations)
				spectra.add(prov.getRotD50(rotation.site, rotation, 0));
			for (int p=0; p<periods.length; p++) {
				double[] values = new double[spectra.size()];
				for (int i=0; i<values.length; i++)
					values[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
				ret.get(p).add(new ResidualSet(rotations, values));
			}
		} else {
//			System.out.println("have "+rotations.size()+" rotations before "+groupQuantities[0]);
			Quantity[] current = {groupQuantities[0]};
			Quantity[] downstream = Arrays.copyOfRange(groupQuantities, 1, groupQuantities.length);
			for (Object value : magQuantitiesTable.get(magnitude, groupQuantities[0])) {
				List<RotationSpec> myRotations = RotatedRupVariabilityConfig.getRotationsForQuantities(
						rotations, current, new Object[] {value});
//				System.out.println("\t"+value+": "+myRotations.size());
				calcGroupedResiduals(magnitude, myRotations, periods, ret, downstream);
			}
		}
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
		final double[] medianPhi;
		final double[] medianTau;
		final double[] medianTotal;
		
		public GMPE_Result(ScalarGroundMotion[][] gms) {
			this.gms = gms;
			logMedian = new double[gms.length];
			medianPhi = new double[gms.length];
			medianTau = new double[gms.length];
			medianTotal = new double[gms.length];
			for (int p=0; p<gms.length; p++) {
				double[] medians = new double[gms[p].length];
				double[] phis = new double[gms[p].length];
				double[] taus = new double[gms[p].length];
				double[] totals = new double[gms[p].length];
				for (int i=0; i<gms[p].length; i++) {
					ScalarGroundMotion gm = gms[p][i];
					medians[i] = gm.mean();
					phis[i] = gm.phi();
					taus[i] = gm.tau();
					totals[i] = gm.stdDev();
				}
				logMedian[p] = DataUtils.median(medians);
				medianPhi[p] = DataUtils.median(phis);
				medianTau[p] = DataUtils.median(taus);
				medianTotal[p] = DataUtils.median(totals);
			}
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
					if (eventID == 165223 && p == 0)
						System.out.println("GMPE\t"+eventID+"\tsAz="+(float)sourceAz+"\trstAz="+restrictSourceAz
								+"\t\trRup="+(float)rRup+"\trJB="+(float)rJB+"\trX="+(float)rX+"\tgm="+gmpe.getGroundMotion().mean());
				}
				ind++;
			}
		}
		return new GMPE_Result(result);
	}
	
	private void plotPeriodDependentStdDevs(File resourcesDir, String prefix, String title,
			List<DiscretizedFunc> siteStdDevFuncs, DiscretizedFunc totalStdDevFunc) throws IOException {
		Preconditions.checkState(siteStdDevFuncs == null || siteStdDevFuncs.size() == siteColors.size());
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; siteStdDevFuncs != null && i<siteStdDevFuncs.size(); i++) {
			funcs.add(siteStdDevFuncs.get(i));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, siteColors.get(i)));
		}
		
		funcs.add(totalStdDevFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Period (s)", "Standard Deviation");
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		double minX = Double.POSITIVE_INFINITY;
		double maxX = Double.NEGATIVE_INFINITY;
		for (DiscretizedFunc func : funcs) {
			minX = Double.min(minX, func.getMinX());
			maxX = Double.max(maxX, func.getMaxX());
		}
		
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
		CPT distCPT = new CPT(StatUtils.min(distArray), StatUtils.max(distArray), Color.GRAY, Color.BLACK);
		
		for (Float dist : dists) {
			VarGroupingKey key = new VarGroupingKey(type, mag, dist, sites.size() == 1 ? sites.get(0) : null);
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
	
	private XY_DataSet line(double x0, double y0, double x1, double y1, String name) {
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.setName(name);
		line.set(x0, y0);
		line.set(x1, y1);
		return line;
	}
	
	private File[][] plotScatter(File resourcesDir, String prefix, double[] periods, VariabilityType type,
			VariabilityResult result, ScatterDisaggQuantity scatterQuantity) throws IOException {
		System.out.println("Scatter will have "+result.commonRotationSpecs.length+" x values");
		double[] xs = new double[result.commonRotationSpecs.length];
		for (int i=0; i<xs.length; i++)
			xs[i] = scatterQuantity.getValue(this, result.commonRotationSpecs[i]);

		String[] yLabels = { "Standard Deviation", "Median SA" };
		String[] yPrefixes = { "std_dev", "residual", "median" };
		Range[] yRanges = { new Range (0, 1), null };
		boolean[] yLogs = { false, true };
		
		File[][] ret = new File[yLabels.length][periods.length];
		
		for (int yI=0; yI<yLabels.length; yI++) {
			if (type.stdDevOfMedians && yI == 0)
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
					for (int p=0; p<periods.length; p++) {
						double gmpeSD = -1;
						if (gmpeResult != null && type.gmpeStdDevType != null) {
							if (type.gmpeStdDevType == GMPE_Var.PHI)
								gmpeSD = gmpeResult.medianPhi[p];
							else if (type.gmpeStdDevType == GMPE_Var.TAU)
								gmpeSD = gmpeResult.medianTau[p];
							else
								gmpeSD = gmpeResult.medianTotal[p];
						}
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
			Collection<VariabilityType> types, Site site, Double magnitude, double[] periods) throws IOException {
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
					VarGroupingKey key = new VarGroupingKey(type, magnitude, distance, site);
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
		
		Map<VariabilityType, File[]> retMap = new HashMap<>();
		
		Range xRange = new Range(minValue-0.5*deltaValue, values.get(values.size()-1)+0.5*deltaValue);
		Range yRange = new Range(minDist-0.5*deltaDist, distances.get(distances.size()-1)+0.5*deltaDist);
		
		for (VariabilityType type : xyzsMap.keySet()) {
			File[] ret = new File[periods.length];
			retMap.put(type, ret);
			
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
				String myPrefix, title;
				if (type == null) {
					title = periodStr+" "+quantity.getName()+", Median SA";
					myPrefix = prefix+"_"+periodStr+"_median_sa";
				} else {
					title = periodStr+" "+quantity.getName()+", "+type.name+" ("+type.symbol+")";
					myPrefix = prefix+"_"+periodStr+"_"+type.prefix;
				}
				
				String xAxisLabel = quantity.getName();
				String yAxisLabel = BBP_PartBValidationConfig.DIST_JB ? "DistanceJB" : "Distance";
				
				XYZPlotSpec xyzSpec = new XYZPlotSpec(xyzs[p], cpt, title, xAxisLabel, yAxisLabel, zLabel);
				XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
				xyzGP.drawPlot(xyzSpec, false, false, xRange, yRange);
				xyzGP.getChartPanel().getChart().setBackgroundPaint(Color.WHITE);
				xyzGP.getChartPanel().setSize(700, 550);
				File file = new File(resourcesDir, myPrefix+".png");
				xyzGP.saveAsPNG(file.getAbsolutePath());
				ret[p] = file;
			}
		}
		return retMap;
	}
	
	private Map<VariabilityType, File[]> plotQuantityDistFuncs(File resourcesDir, String prefix, Quantity quantity,
			Collection<VariabilityType> types, Site site, Double magnitude, double[] periods, List<Float> distances) throws IOException {
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
					VarGroupingKey key = new VarGroupingKey(type, magnitude, distance, site);
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
						for (int p=0; p<periods.length; p++) {
							double gmpeSD;
							if (type.gmpeStdDevType == GMPE_Var.PHI)
								gmpeSD = gmpeResult.medianPhi[p];
							else if (type.gmpeStdDevType == GMPE_Var.TAU)
								gmpeSD = gmpeResult.medianTau[p];
							else
								gmpeSD = gmpeResult.medianTotal[p];
							gmpeXYs[d][p].set(value, gmpeSD);
						}
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
		
		Map<VariabilityType, File[]> retMap = new HashMap<>();
		
		double[] distArray = Doubles.toArray(distances);
		CPT distCPT = new CPT(StatUtils.min(distArray), StatUtils.max(distArray), Color.GRAY, Color.BLACK);
		
		for (VariabilityType type : xysMap.keySet()) {
			File[] ret = new File[periods.length];
			retMap.put(type, ret);
			
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
				String myPrefix, title;
				if (type == null) {
					title = periodStr+" "+quantity.getName()+", Median SA";
					myPrefix = prefix+"_"+periodStr+"_median_sa";
				} else {
					title = periodStr+" "+quantity.getName()+", "+type.name+" ("+type.symbol+")";
					myPrefix = prefix+"_"+periodStr+"_"+type.prefix;
				}
				
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
				File pngFile = new File(resourcesDir, myPrefix+".png");
				gp.saveAsPNG(pngFile.getAbsolutePath());
				ret[p] = pngFile;
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
