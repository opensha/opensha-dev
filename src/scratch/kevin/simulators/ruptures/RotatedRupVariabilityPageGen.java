package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.data.Range;
import org.opensha.commons.data.NamedComparator;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.ArbDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public abstract class RotatedRupVariabilityPageGen {

	private RSQSimCatalog catalog;
	private Map<Double, RotatedRupVariabilityConfig> magConfigs;
	private Map<Double, SimulationRotDProvider<RotationSpec>> magProvs;
	
	private List<Site> sites;
	private List<Float> sourceAzimuths;
	private List<Float> siteSourceAzimuths;
	private List<Float> distances;
	private Map<Double, List<Integer>> magEventIDs;
	private int maxNumEvents = 0;
	private int minNumEvents = Integer.MAX_VALUE;
	private Table<Double, Quantity, List<?>> magQuantitiesTable;
	
	private RSQSimEvent exampleRupture;
	private Site exampleSite;
	private int numExampleRotations;
	
	private List<Color> siteColors;
	
	private LoadingCache<StdDevSetKey, StdDevSet[]> stdDevCache;

	public RotatedRupVariabilityPageGen(RSQSimCatalog catalog, RotatedRupVariabilityConfig config,
			SimulationRotDProvider<RotationSpec> prov) {
		this(catalog, emptyMagMap(config), emptyMagMap(prov));
	}
	
	private static <T> HashMap<Double, T> emptyMagMap(T value) {
		HashMap<Double, T> map = new HashMap<>();
		map.put(null, value);
		return map;
	}

	public RotatedRupVariabilityPageGen(RSQSimCatalog catalog, Map<Double, RotatedRupVariabilityConfig> magConfigs,
			Map<Double, SimulationRotDProvider<RotationSpec>> magProvs) {
		this.catalog = catalog;
		this.magConfigs = magConfigs;
		this.magProvs = magProvs;
		
		RotatedRupVariabilityConfig config0 = magConfigs.values().iterator().next();
		
		sites = config0.getValues(Site.class, Quantity.SITE);
		sourceAzimuths = config0.getValues(Float.class, Quantity.SOURCE_AZIMUTH);
		siteSourceAzimuths = config0.getValues(Float.class, Quantity.SITE_TO_SOURTH_AZIMUTH);
		distances = config0.getValues(Float.class, Quantity.DISTANCE);
		magQuantitiesTable = HashBasedTable.create();
		
		magEventIDs = new HashMap<>();
		if (magConfigs.size() == 1) {
			List<Integer> eventIDs = config0.getValues(Integer.class, Quantity.EVENT_ID);
			magEventIDs.put(null, eventIDs);
			maxNumEvents = eventIDs.size();
			minNumEvents = eventIDs.size();
			for (Quantity quantity : Quantity.values())
				magQuantitiesTable.put(null, quantity, config0.getQuantitiesMap().get(quantity));
		} else {
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
		
		stdDevCache = CacheBuilder.newBuilder().build(new CacheLoader<StdDevSetKey, StdDevSet[]>() {

			@Override
			public StdDevSet[] load(StdDevSetKey key) throws Exception {
				return calcStdDevSet(key);
			}
			
		});
	}
	
	public void setExampleRupture(RSQSimEvent exampleRupture, Site exampleSite, int numRotations) {
		this.exampleRupture = exampleRupture;
		this.exampleSite = exampleSite;
		this.numExampleRotations = numRotations;
	}
	
	private static final String al_atik = "Al Atik (2010)";
	private static final String aki_richards = "Aki & Richards (1980)";
	
	private enum VariabilityType {
		PATH("Path-to-path", "path", "φ_p2p", "&phi;<sub>P2P</sub>", al_atik,
				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH),
		SOUCE_STRIKE("Source-strike", "source_strike", "φ_s", "&phi;<sub>s</sub>", aki_richards,
				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SITE_TO_SOURTH_AZIMUTH),
		WITHIN_EVENTS("Within-event", "within_event", "φ", "&phi;", al_atik,
				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.EVENT_ID),
		BETWEEN_EVENTS_SINGLE_PATH("Between-events, single-path", "between_events_single_path", "τ_0", "&tau;<sub>0</sub>", al_atik,
				true, true, Quantity.SITE, Quantity.DISTANCE, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH),
		BETWEEN_EVENTS("Between-events", "between_events", "τ", "&tau;", al_atik,
				true, true, Quantity.SITE, Quantity.DISTANCE),
		SITE_TO_SITE("Site-to-site", "site_to_site", "φ_s2s", "&phi;<sub>S2S</sub>", al_atik,
				false, true, Quantity.DISTANCE, Quantity.EVENT_ID, Quantity.SOURCE_AZIMUTH, Quantity.SITE_TO_SOURTH_AZIMUTH);
		
		private String name;
		private String prefix;
		private String symbol;
		private String htmlSymbol;
		private String reference;
		private boolean separateSites;
		private boolean separateDists;
		private Quantity[] groupQuantities;
		
		private List<Quantity> groupedList;
		private List<Quantity> variedList;

		private VariabilityType(String name, String prefix, String symbol, String htmlSymbol, String reference, boolean separateSites,
				boolean separateDists, Quantity... groupQuantities) {
			this.name = name;
			this.prefix = prefix;
			this.symbol = symbol;
			this.htmlSymbol = htmlSymbol;
			this.reference = reference;
			this.separateSites = separateSites;
			this.separateDists = separateDists;
			this.groupQuantities = groupQuantities;
			
			groupedList = new ArrayList<>();
			variedList = Lists.newArrayList(Quantity.values());
			
			for (Quantity quantity : groupQuantities) {
				groupedList.add(quantity);
				variedList.remove(quantity);
			}
		}
		
		public int getGroupedCount(Map<Quantity, List<?>> quantitiesMap) {
			int groupedCount = 1;
			for (Quantity quantity : groupedList) {
				int num = quantitiesMap.get(quantity).size();
				groupedCount *= num;
			}
			return groupedCount;
		}
		
		public int getVariedCount(Map<Quantity, List<?>> quantitiesMap) {
			int variedCount = 1;
			for (Quantity quantity : variedList) {
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
			line += ", is computed from ground motion residuals where the following quantities are held constant:";
			lines.add(line);
			lines.add("");
			
			for (Quantity quantity : groupedList) {
				int num = pageGen.magQuantitiesTable.rowMap().values().iterator().next().get(quantity).size();
				lines.add("* "+quantity+" *["+num+" unique]*");
			}
			lines.add("");
			lines.add("and the following quantities vary:");
			lines.add("");
			for (Quantity quantity : variedList) {
				int num = pageGen.magQuantitiesTable.rowMap().values().iterator().next().get(quantity).size();
				lines.add("* "+quantity+" *["+num+" unique]*");
			}
			lines.add("");
//					+ "orientation constant. Here is an exmample with "+pageGen.numExampleRotations+" rotations:");
			if (pageGen.exampleRupture != null) {
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
					pageGen.plotExample(resourcesDir, "example_"+this.prefix, pageGen.distances.get(0), variedList);
					lines.add("![Example](resources/example_"+this.prefix+".png)");
					lines.add("");
				} else {
					System.out.println("Not plotting example. SeparateSites? "+separateSites
							+". Quantities: "+Joiner.on("; ").join(groupQuantities));
				}
			}
			
			if (separateSites || separateDists) {
				line = "Standard deviation is computed and tabulated separately for each ";
				line += separateSites && separateDists ? "site and distance" : separateSites ? "site" : "distance";
				if (separateSites)
					line += ", then a total standard deviation is computed from all sitess and reported in the \"**ALL SITES**\" row.";
				else
					line += ".";
				if (separateDists)
					line += " Results are reported separately for each distance";
				lines.add(line);
			}
			
			return lines;
		}
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
		
		lines.add("# "+catalog.getName()+" Rotated Rupture Variability, "+getScenarioShortName());
		lines.add("");
		lines.add("This exercise uses translations and rotations to estimate ground motion variability from different "
				+ "sources. We begin by selecting a subset of similar ruptures which match a set of criteria (in this case, "
				+ getScenarioName()+"). Each rupture is then reoriented such that its strike (following the Aki & Richards "
				+ "1980 convention) is 0 degrees (due North, dipping to the right for normal or reverse ruptures). For each site, "
				+ "ruptures are translated such that their scalar seismic moment centroid is directly North of the site, and their "
				+ "Joyner-Boore distance (rJB) is as specified (we consider "+distances.size()+" distance[s] here).");
		lines.add("");
		lines.add("We then  perform various rotations. We rotate the rupture in place around its centroid, holding the site-to-source "
				+ "centroid path and rJB constant (henceforth '"+Quantity.SOURCE_AZIMUTH.getName()+"'). We also rotate ruptures around the site, "
				+ "holding rJB and source orientation relative to the site constant but sampling different various paths (henceforth "
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
		if (magConfigs.size() == 1) {
			plotMags.add(null);
		} else {
			if (highlightMags != null) {
				for (double mag : highlightMags)
					if (magConfigs.containsKey(mag))
						plotMags.add(mag);
			} else {
				plotMags.addAll(magConfigs.keySet());
			}
			Collections.sort(plotMags);
		}
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
			if (type.getVariedCount(magQuantitiesTable.rowMap().values().iterator().next()) == 1)
				continue;
			computedTypes.add(type);
			
			System.out.println("Processing variability type: "+type.name);
			lines.add("## "+type.name+" Variability");
			lines.add(topLink); lines.add("");
			
			lines.add("### "+type.name+" Variability Methodology");
			lines.add(topLink); lines.add("");
			
			lines.addAll(type.buildMethodologyLines(this, resourcesDir));
			lines.add("");
			
			if (hasMagDist) {
				lines.add("### "+type.name+" Variability Mag-Distance Plots");
				lines.add(topLink); lines.add("");
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				if (sites.size() > 1 && type.separateSites)
					table.addColumn("Site");
				for (double period : periods)
					table.addColumn(period == (int)period ? (int)period+"s" : (float)period+"s");
				table.finalizeLine();
				if (sites.size() > 1 && type.separateSites) {
					for (Site site : sites) {
						table.initNewLine();
						table.addColumn("**"+site.getName()+"**");
						String prefix = type.prefix+"mag_dist_std_dev";
						File[] files = plotMagDistCheckerboard(resourcesDir, prefix, type, site, periods);
						for (File file : files)
							table.addColumn("![Mag-Dist Plot](resources/"+file.getName()+")");
						table.finalizeLine();
					}
					table.initNewLine();
					table.addColumn("**ALL SITES***");
				} else {
					table.initNewLine();
				}
				String prefix = type.prefix+"_mag_dist_std_dev";
				File[] files = plotMagDistCheckerboard(resourcesDir, prefix, type, null, periods);
				for (File file : files)
					table.addColumn("![Mag-Dist Plot](resources/"+file.getName()+")");
				magDistPlots.put(type, files);
				table.finalizeLine();
				lines.addAll(table.build());
				lines.add("");
			}
			
			for (Double mag : plotMags) {
				RotatedRupVariabilityConfig config = magConfigs.get(mag);
				if (!type.separateDists || distances.size() == 1) {
					if (mag == null)
						lines.add("### "+type.name+" Variability Results");
					else
						lines.add("### M"+optionalDigitDF.format(mag)+" "+type.name+" Results");
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
		}
		
		// build summary table(s)
		List<String> summaryLines = new ArrayList<>();
		if (hasMagDist) {
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
				table.addColumn(type.name);
				table.addColumn(type.htmlSymbol);
				File[] plots = magDistPlots.get(type);
				for (File plot : plots)
					table.addColumn("![Mag-Dist Plot](resources/"+plot.getName()+")");
				table.finalizeLine();
			}
			summaryLines.addAll(table.build());
		}
		for (Double mag : plotMags) {
			if (mag != null) {
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
					StdDevSetKey key = new StdDevSetKey(type, mag, distance, null, periods);
					StdDevSet[] stdDevs;
					try {
						stdDevs = stdDevCache.get(key);
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
					for (StdDevSet stdDev : stdDevs)
						table.addColumn(optionalDigitDF.format(stdDev.total));
					table.finalizeLine();
				}
			}
			summaryLines.addAll(table.build());
		}
		lines.addAll(summaryTableIndex, summaryLines);
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private List<String> generateVariabilityLines(RotatedRupVariabilityConfig config, VariabilityType type, Double mag, Float distance,
			double[] periods, String topLink, File resourcesDir) throws IOException {
		List<String> lines = new ArrayList<>();
		if (type.separateDists && distances.size() > 1) {
			if (mag == null)
				lines.add("### "+distance+" km "+type.name+" Results");
			else
				lines.add("### "+distance+" km M"+optionalDigitDF.format(mag)+" "+type.name+" Results");
			lines.add(topLink); lines.add("");
		}
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		if (type.separateSites && sites.size() > 1)
			table.addColumn("Site");
		for (double period : periods) {
			String periodStr = optionalDigitDF.format(period)+"s";
			table.addColumn(periodStr+" "+type.htmlSymbol);
			table.addColumn("Total");
			table.addColumn("Mean");
			table.addColumn("Median");
			table.addColumn("Range");
		}
		table.finalizeLine();
		DiscretizedFunc totStdDevFunc = null;
		StdDevSet[] totStdDevs = null;
		List<List<StdDevSet>> siteTotStdDevs = null;
		if (type.separateSites && sites.size() > 1) {
			siteTotStdDevs = new ArrayList<>();
			for (int i=0; i<periods.length; i++)
				siteTotStdDevs.add(new ArrayList<>());
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
			StdDevSet[] myStdDevs;
			try {
				myStdDevs = stdDevCache.get(new StdDevSetKey(type, mag, distance, site, periods));
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			for (int p=0; p<periods.length; p++) {
				table.addColumn(""); // empty for period
				table.addColumn(emphasis+optionalDigitDF.format(myStdDevs[p].total)+emphasis);
				table.addColumn(emphasis+optionalDigitDF.format(myStdDevs[p].mean)+emphasis);
				table.addColumn(emphasis+optionalDigitDF.format(myStdDevs[p].median)+emphasis);
				table.addColumn(emphasis+"["+optionalDigitDF.format(myStdDevs[p].min)+" "+optionalDigitDF.format(myStdDevs[p].max)+"]"+emphasis);
				stdDevFunc.set(periods[p], myStdDevs[p].total);
			}
			if (site == null)
				totStdDevs = myStdDevs;
			else
				for (int p=0; p<periods.length; p++)
					siteTotStdDevs.get(p).add(myStdDevs[p]);
			table.finalizeLine();
		}

		String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
		String magPrefix = mag == null ? "" : "_m"+optionalDigitDF.format(mag);
		
		// plot it
		String prefix = type.prefix+magPrefix+distPrefix+"_std_dev";
		plotPeriodDependentStdDevs(resourcesDir, prefix, type.name+" ("+type.symbol+")", periods, stdDevFuncs, totStdDevFunc);
		lines.add("!["+type.name+" Variability](resources/"+prefix+".png)");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		for (double period : periods)
			table.addColumn(optionalDigitDF.format(period)+"s");
		table.finalizeLine();
		table.initNewLine();
		for (int p=0; p<periods.length; p++) {
			prefix = type.prefix+magPrefix+distPrefix+"_"+optionalDigitDF.format(periods[p])+"s_hist";
			plotStdDevsHistogram(resourcesDir, prefix, optionalDigitDF.format(periods[p])+"s "+type.name+" ("+type.symbol+")",
					periods[p], totStdDevs[p], type.separateSites && sites.size() > 1 ? siteTotStdDevs.get(p): null);
			table.addColumn("!["+optionalDigitDF.format(periods[p])+"s](resources/"+prefix+".png)");
		}
		table.finalizeLine();
		table.wrap(periods.length > 4 ? 3 : 2, 0);
		lines.addAll(table.build());
		lines.add("");
		return lines;
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	private static Float nullAsZero(Float value) {
		if (value == null)
			return 0f;
		return value;
	}
	
	private class ResidualSet {
		double[] zeroCenteredValues;
		double meanValue;
		double residualStdDev;
		
		public ResidualSet(double[] values) {
			this.zeroCenteredValues = values;
			this.meanValue = StatUtils.mean(values);
			Preconditions.checkState(Double.isFinite(meanValue), "Non-finite mean: %s", meanValue);
			for (int i=0; i<values.length; i++)
				zeroCenteredValues[i] -= meanValue;
			double variance = StatUtils.variance(zeroCenteredValues, meanValue);
			if (variance < 0) {
				Preconditions.checkState(variance > -1e-10, "Negative variance: %s", variance);
				variance = 0;
			}
			this.residualStdDev = Math.sqrt(variance);
			Preconditions.checkState(Double.isFinite(residualStdDev),
					"Non-finite std dev: %s, var=%s, mean=%s, size=%s", residualStdDev, variance, meanValue, values.length);
		}
		
		public int size() {
			return zeroCenteredValues.length;
		}
		
		public double getMean() {
			return meanValue;
		}
		
		public double[] getZeroCentered() {
			return zeroCenteredValues;
		}
		
		public double[] getRawValues() {
			double[] ret = new double[zeroCenteredValues.length];
			for (int i=0; i<ret.length; i++)
				ret[i] = zeroCenteredValues[i]+meanValue;
			return ret;
		}
		
		public double getStdDev() {
			return residualStdDev;
		}
	}
	
	private List<List<ResidualSet>> loadResiduals(StdDevSetKey key) throws IOException {
		Quantity[] constQuantities = null;
		Object[] constValues = null;
		if (key.type.separateSites) {
			if (key.type.separateDists) {
				constQuantities = new Quantity[] { Quantity.DISTANCE, Quantity.SITE };
				if (sites.size() == 1 && key.site == null)
					constValues = new Object[] { key.distance, sites.get(0) };
				else
					constValues = new Object[] { key.distance, key.site };
			} else {
				constQuantities = new Quantity[] { Quantity.SITE };
				if (sites.size() == 1 && key.site == null)
					constValues = new Object[] { sites.get(0) };
				else
					constValues = new Object[] { key.site };
			}
		} else {
			if (key.type.separateDists) {
				constQuantities = new Quantity[] { Quantity.DISTANCE };
				constValues = new Object[] { key.distance };
			}
		}
		
		return calcResiduals(key.magnitude, key.periods, constQuantities, constValues, key.type.groupQuantities);
	}
		
	
	private StdDevSet[] calcStdDevSet(StdDevSetKey key) throws IOException {
		if (key.type.separateSites && sites.size() > 1) {
			// do all sites at once
			List<List<ResidualSet>> totalPeriodResiduals = new ArrayList<>();
			for (int p=0; p<key.periods.length; p++)
				totalPeriodResiduals.add(new ArrayList<>());
			
			for (Site site : sites) {
				StdDevSetKey siteKey = new StdDevSetKey(key.type, key.magnitude, key.distance, site, key.periods);
				
				List<List<ResidualSet>> residuals = loadResiduals(siteKey);
				
				StdDevSet[] stdDevs = new StdDevSet[key.periods.length];
				for (int p=0; p<key.periods.length; p++) {
					stdDevs[p] = new StdDevSet(residuals.get(p));
					totalPeriodResiduals.get(p).addAll(residuals.get(p));
				}
				
				stdDevCache.put(siteKey, stdDevs);
			}
			
			StdDevSet[] totStdDevs = new StdDevSet[key.periods.length];
			for (int p=0; p<key.periods.length; p++)
				totStdDevs[p] = new StdDevSet(totalPeriodResiduals.get(p));
			StdDevSetKey noSiteKey = new StdDevSetKey(key.type, key.magnitude, key.distance, null, key.periods);
			stdDevCache.put(noSiteKey, totStdDevs);
			
			StdDevSet[] ret = stdDevCache.getIfPresent(key);
			Preconditions.checkState(ret != null, "Should have just been cached?");
			return ret;
		}
		// do for all sites
		List<List<ResidualSet>> residuals = loadResiduals(key);
		
		StdDevSet[] stdDevs = new StdDevSet[key.periods.length];
		for (int p=0; p<key.periods.length; p++)
			stdDevs[p] = new StdDevSet(residuals.get(p));
		return stdDevs;
	}
	
	private class StdDevSetKey {
		private final VariabilityType type;
		private final Double magnitude;
		private final Float distance;
		private final Site site;
		private final double[] periods;
		
		public StdDevSetKey(VariabilityType type, Double magnitude, Float distance, Site site, double[] periods) {
			this.type = type;
			this.magnitude = magnitude;
			this.distance = distance;
			this.site = site;
			this.periods = periods;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + ((distance == null) ? 0 : distance.hashCode());
			result = prime * result + ((magnitude == null) ? 0 : magnitude.hashCode());
			result = prime * result + Arrays.hashCode(periods);
			result = prime * result + ((site == null) ? 0 : site.getLocation().hashCode());
			result = prime * result + ((type == null) ? 0 : type.hashCode());
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
			StdDevSetKey other = (StdDevSetKey) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (distance == null) {
				if (other.distance != null)
					return false;
			} else if (!distance.equals(other.distance))
				return false;
			if (magnitude == null) {
				if (other.magnitude != null)
					return false;
			} else if (!magnitude.equals(other.magnitude))
				return false;
			if (!Arrays.equals(periods, other.periods))
				return false;
			if (site == null) {
				if (other.site != null)
					return false;
			} else if (!site.equals(other.site))
				return false;
			if (type != other.type)
				return false;
			return true;
		}

		private RotatedRupVariabilityPageGen getOuterType() {
			return RotatedRupVariabilityPageGen.this;
		}

	}
	
	private class StdDevSet {
		private final double total;
		private final double mean;
		private final double median;
		private final double min;
		private final double max;
		private final double[] stdDevs;
		
		public StdDevSet(List<ResidualSet> residuals) {
			stdDevs = new double[residuals.size()];
			StandardDeviation totalCalc = new StandardDeviation();
			for (int i=0; i<stdDevs.length; i++) {
				ResidualSet set = residuals.get(i);
				stdDevs[i] = set.getStdDev();
				totalCalc.incrementAll(set.getZeroCentered());
			}
			total = totalCalc.getResult();
			mean = StatUtils.mean(stdDevs);
			median = DataUtils.median(stdDevs);
			min = StatUtils.min(stdDevs);
			max = StatUtils.max(stdDevs);
		}
	}
	
	private List<List<ResidualSet>> calcResiduals(Double magnitude, double[] periods,
			Quantity[] constQuantities, Object[] constValues, Quantity[] groupQuantities) throws IOException {
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
		for (int i=0; i<periods.length; i++)
			ret.add(new ArrayList<>());
		
		System.out.println("Computing residuals for Constants ["
				+(constQuantities == null ? "" : Joiner.on(",").join(constQuantities))
						+"] and Groups ["+Joiner.on(",").join(groupQuantities)+"]: "+num+" simulations");
		
//		System.out.println("Computing resduals from "+num+" simulations");
		
		calcGroupedResiduals(magnitude, totalRotations, periods, ret, groupQuantities);
		int count = 0;
		for (ResidualSet set : ret.get(new Random().nextInt(periods.length)))
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
				ret.get(p).add(new ResidualSet(values));
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
	
//	private List<double[]> calcSourceAzResiduals(double[] periods, Site site, float distance) throws IOException {
//		int num = sourceAzimuths.size()*siteSourceAzimuths.size()*eventIDs.size();
//		List<double[]> ret = new ArrayList<>();
//		for (int i=0; i<periods.length; i++)
//			ret.add(new double[num]);
//		
//		int index = 0;
//		double[] workingArray = new double[sourceAzimuths.size()];
//		System.out.println("Calculating path residuals with "+workingArray.length+" azimuths per site");
//		int numThatVary = 0;
//		int totNum = 0;
//		for (int eventID : eventIDs) {
//			for (float siteSourceAzimuth : siteSourceAzimuths) {
//				List<RotationSpec> rotations = config.getSourceRotations(site, eventID, distance, siteSourceAzimuth);
//				Preconditions.checkState(rotations.size() == workingArray.length,
//						"Unexpected number of configurations for event %s, site %s, distance %s, az %s. Expected %s, have %s",
//						eventID, site.getName(), distance, siteSourceAzimuth, workingArray.length, rotations.size());
//				List<DiscretizedFunc> spectra = new ArrayList<>();
//				for (RotationSpec rotation : rotations)
//					spectra.add(prov.getRotD50(site, rotation, 0));
//				boolean varies = false;
//				for (int p=0; p<periods.length; p++) {
//					for (int i=0; i<workingArray.length; i++)
//						workingArray[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
//					// now detrend to remove any source effects
//					double mean = StatUtils.mean(workingArray);
//					for (int i=0; i<workingArray.length; i++) {
//						ret.get(p)[index+i] = workingArray[i]-mean;
//						if (Math.abs(ret.get(p)[index+i]) > 1e-2)
//							varies = true;
//					}
//				}
//				index += workingArray.length;
//				totNum++;
//				if (varies)
//					numThatVary++;
//			}
//		}
//		Preconditions.checkState(index == num, "Bad end index. Expected %s, have %s", num, index);
//		System.out.println(numThatVary+"/"+totNum+" have source az variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
//		
//		return ret;
//	}
//	
//	private List<double[]> calcEventResiduals(double[] periods, Site site, float distance) throws IOException {
//		int num = sourceAzimuths.size()*siteSourceAzimuths.size()*eventIDs.size();
//		List<double[]> ret = new ArrayList<>();
//		for (int i=0; i<periods.length; i++)
//			ret.add(new double[num]);
//		
//		int index = 0;
//		double[] workingArray = new double[eventIDs.size()];
//		System.out.println("Calculating path residuals with "+workingArray.length+" azimuths per site");
//		int numThatVary = 0;
//		int totNum = 0;
//		for (float sourceAzimuth : sourceAzimuths) {
//			for (float siteSourceAzimuth : siteSourceAzimuths) {
//				List<RotationSpec> rotations = config.getSourcesForRotation(site, distance, sourceAzimuth, siteSourceAzimuth);
//				Preconditions.checkState(rotations.size() == workingArray.length,
//						"Unexpected number of configurations for site %s, distance %s, sourceAz %s, siteSourceAz %s. Expected %s, have %s",
//						site.getName(), distance, sourceAzimuth, siteSourceAzimuth, workingArray.length, rotations.size());
//				List<DiscretizedFunc> spectra = new ArrayList<>();
//				for (RotationSpec rotation : rotations)
//					spectra.add(prov.getRotD50(site, rotation, 0));
//				boolean varies = false;
//				for (int p=0; p<periods.length; p++) {
//					for (int i=0; i<workingArray.length; i++)
//						workingArray[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
//					// now detrend to remove any source effects
//					double mean = StatUtils.mean(workingArray);
//					for (int i=0; i<workingArray.length; i++) {
//						ret.get(p)[index+i] = workingArray[i]-mean;
//						if (Math.abs(ret.get(p)[index+i]) > 1e-2)
//							varies = true;
//					}
//				}
//				index += workingArray.length;
//				totNum++;
//				if (varies)
//					numThatVary++;
//			}
//		}
//		Preconditions.checkState(index == num, "Bad end index. Expected %s, have %s", num, index);
//		System.out.println(numThatVary+"/"+totNum+" have event variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
//		
//		return ret;
//	}
	
	private void plotPeriodDependentStdDevs(File resourcesDir, String prefix, String title, double[] periods,
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
		
		Range xRange = new Range(StatUtils.min(periods), StatUtils.max(periods));
		Range yRange = new Range(0d, 1d);

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
	}
	
	private void plotStdDevsHistogram(File resourcesDir, String prefix, String title, double period,
			StdDevSet totalStdDev, List<StdDevSet> siteStdDevs) throws IOException {
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
				StdDevSet siteStdDev = siteStdDevs.get(i);
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
	
	private File[] plotMagDistCheckerboard(File resourcesDir, String prefix, VariabilityType type, Site site, double[] periods)
			throws IOException {
		List<Double> magnitudes = new ArrayList<>(magConfigs.keySet());
		Collections.sort(magnitudes);
		
		double deltaMag = magnitudes.get(1) - magnitudes.get(0);
		double minMag = magnitudes.get(0);
		
		double deltaDist = distances.get(1) - distances.get(0);
		double minDist = distances.get(0);
		
		EvenlyDiscrXYZ_DataSet[] xyzs = new EvenlyDiscrXYZ_DataSet[periods.length];
		for (int p=0; p<periods.length; p++)
			xyzs[p] = new EvenlyDiscrXYZ_DataSet(distances.size(), magnitudes.size(), minDist, minMag, deltaDist, deltaMag);
		
		for (int xInd=0; xInd<distances.size(); xInd++) {
			float distance = distances.get(xInd);
			double x = xyzs[0].getX(xInd);
			Preconditions.checkState(distance == (float)x, "Expected distance bin %s to be %s, have %s", xInd, distance, x);
			for (int yInd=0; yInd<magnitudes.size(); yInd++) {
				double magnitude = magnitudes.get(yInd);
				double y = xyzs[0].getY(yInd);
				Preconditions.checkState((float)magnitude == (float)y, "Expected mag bin %s to be %s, have %s", yInd, magnitude, y);
				StdDevSetKey key = new StdDevSetKey(type, magnitude, distance, site, periods);
				StdDevSet[] stdDevs;
				try {
					stdDevs = stdDevCache.get(key);
				} catch (ExecutionException e) {
					if (e.getCause() instanceof IOException)
						throw (IOException)e.getCause();
					if (e.getCause() instanceof RuntimeException)
						throw (RuntimeException)e.getCause();
					throw ExceptionUtils.asRuntimeException(e);
				}
				for (int p=0; p<periods.length; p++) {
					xyzs[p].set(xInd, yInd, stdDevs[p].total);
				}
			}
		}
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		File[] ret = new File[periods.length];
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		
		for (int p=0; p<periods.length; p++) {
			String periodStr;
			if (periods[p] == Math.round(periods[p]))
				periodStr = (int)periods[p]+"s";
			else
				periodStr = (float)periods[p]+"s";
			String title = periodStr+" "+type.name+" ("+type.symbol+")";
			String myPrefix = prefix+"_"+periodStr;
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyzs[p], cpt, title, "DistanceJB", "Magnitude", type.symbol);
			XYZGraphPanel xyzGP = new XYZGraphPanel(plotPrefs);
			xyzGP.drawPlot(xyzSpec, false, false, new Range(minDist-0.5*deltaDist, distances.get(distances.size()-1)+0.5*deltaDist),
					new Range(minMag-0.5*deltaMag, magnitudes.get(magnitudes.size()-1)+0.5*deltaMag));
			xyzGP.getChartPanel().getChart().setBackgroundPaint(Color.WHITE);
			xyzGP.getChartPanel().setSize(700, 550);
			File file = new File(resourcesDir, myPrefix+".png");
			xyzGP.saveAsPNG(file.getAbsolutePath());
			ret[p] = file;
		}
		return ret;
	}
	
	private void plotExample(File resourcesDir, String prefix, double distance, List<Quantity> variedQuantities)
			throws IOException {
		List<Site> sites = new ArrayList<>();
		sites.add(exampleSite);
		List<RSQSimEvent> ruptures = new ArrayList<>();
		ruptures.add(exampleRupture);
		int numSourceAz = variedQuantities.contains(Quantity.SOURCE_AZIMUTH) ? numExampleRotations : 1;
		int numSiteToSourceAz = variedQuantities.contains(Quantity.SITE_TO_SOURTH_AZIMUTH) ? numExampleRotations : 1;
		RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(catalog, sites, ruptures, new double[] {distance},
				numSourceAz, numSiteToSourceAz);
		
		config.plotRotations(resourcesDir, prefix, config.getRotations(), true);
	}
}
