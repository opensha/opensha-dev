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
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
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
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVariabilityPageGen {

	private RSQSimCatalog catalog;
	private Scenario scenario;
	private RotatedRupVariabilityConfig config;
	private SimulationRotDProvider<RotationSpec> prov;
	
	private List<Site> sites;
	private List<Float> sourceAzimuths;
	private List<Float> siteSourceAzimuths;
	private List<Float> distances;
	private List<Integer> eventIDs;
	private Map<Quantity, List<?>> quantitiesMap;
	
	private RSQSimEvent exampleRupture;
	private Site exampleSite;
	private int numExampleRotations;
	
	private List<Color> siteColors;

	public RotatedRupVariabilityPageGen(RSQSimCatalog catalog, Scenario scenario, RotatedRupVariabilityConfig config,
			SimulationRotDProvider<RotationSpec> prov) {
		this.catalog = catalog;
		this.scenario = scenario;
		this.config = config;
		this.prov = prov;
		
		HashSet<Site> sitesSet = new HashSet<>();
		HashSet<Float> sourceAzSet = new HashSet<>();
		HashSet<Float> siteSourceAzSet = new HashSet<>();
		HashSet<Float> distancesSet = new HashSet<>();
		HashSet<Integer> eventIDsSet = new HashSet<>();
		for (RotationSpec rotation : config.getRotations()) {
			sitesSet.add(rotation.site);
			sourceAzSet.add(nullAsZero(rotation.sourceAz));
			siteSourceAzSet.add(nullAsZero(rotation.siteToSourceAz));
			distancesSet.add(nullAsZero(rotation.distance));
			eventIDsSet.add(rotation.eventID);
		}
		sites = new ArrayList<>(sitesSet);
		sourceAzimuths = new ArrayList<>(sourceAzSet);
		siteSourceAzimuths = new ArrayList<>(siteSourceAzSet);
		distances = new ArrayList<>(distancesSet);
		eventIDs = new ArrayList<>(eventIDsSet);
		Collections.sort(sites, new NamedComparator());
		Collections.sort(sourceAzimuths);
		Collections.sort(siteSourceAzimuths);
		Collections.sort(distances);
		Collections.sort(eventIDs);
		quantitiesMap = new HashMap<>();
		quantitiesMap.put(Quantity.SITE, sites);
		quantitiesMap.put(Quantity.EVENT_ID, eventIDs);
		quantitiesMap.put(Quantity.DISTANCE, distances);
		quantitiesMap.put(Quantity.SOURCE_AZIMUTH, sourceAzimuths);
		quantitiesMap.put(Quantity.SITE_TO_SOURTH_AZIMUTH, siteSourceAzimuths);
		
		CPT siteCPT;
		try {
			siteCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, sites.size()-1);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		siteColors = new ArrayList<>();
		for (int i=0; i<sites.size(); i++)
			siteColors.add(siteCPT.getColor((float)i).darker());
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
				int num = pageGen.quantitiesMap.get(quantity).size();
				lines.add("* "+quantity+" *["+num+" unique]*");
			}
			lines.add("");
			lines.add("and the following quantities vary:");
			lines.add("");
			for (Quantity quantity : variedList) {
				int num = pageGen.quantitiesMap.get(quantity).size();
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
			
//			lines.add("Individual standard deviations are computed for each unique combination of") // TODO
			
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
	
	public void generatePage(File outputDir, double[] periods, List<String> methodSpecificLines) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new LinkedList<>();
		
		lines.add("# "+catalog.getName()+" Rotated Rupture Variability, "+scenario.getShortName());
		lines.add("");
		lines.add("This exercise uses translations and rotations to estimate ground motion variability from different "
				+ "sources. We begin by selecting a subset of similar ruptures which match a set of criteria (in this case, "
				+ scenario.getName()+"). Each rupture is then reoriented such that its strike (following the Aki & Richards "
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
		table.addLine(Quantity.EVENT_ID.getName(), eventIDs.size(), Quantity.EVENT_ID.getDescription());
		table.addLine(Quantity.SITE.getName(), sites.size(), Quantity.SITE.getDescription());
		table.addLine(Quantity.SOURCE_AZIMUTH.getName(), sourceAzimuths.size(), Quantity.SOURCE_AZIMUTH.getDescription());
		table.addLine(Quantity.SITE_TO_SOURTH_AZIMUTH.getName(), siteSourceAzimuths.size(), Quantity.SITE_TO_SOURTH_AZIMUTH.getDescription());
		table.addLine(Quantity.DISTANCE.getName(), Joiner.on(", ").join(distances)+" km", Quantity.DISTANCE.getDescription());
		table.addLine("**Total # Simulations**", "**"+config.getRotations().size()+"**", "Total number of combinations of the above.");
		lines.addAll(table.build());
		lines.add("");
		lines.add("## "+scenario.getShortName()+" RSQSim Rupture Match Criteria");
		lines.add(topLink); lines.add("");
		String[] criteria = scenario.getMatchCriteria();
		lines.add("We condisder "+eventIDs.size()+" events in the catalog match the following criteria:");
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
		
		Table<VariabilityType, Float, List<StdDevSet>> resultsTable = HashBasedTable.create();
		
		for (VariabilityType type : VariabilityType.values()) {
			if (type.getVariedCount(quantitiesMap) == 1)
				continue;
			System.out.println("Processing variability type: "+type.name);
			lines.add("## "+type.name+" Variability");
			lines.add(topLink); lines.add("");
			
			lines.add("### "+type.name+" Variability Methodology");
			lines.add(topLink); lines.add("");
			
			lines.addAll(type.buildMethodologyLines(this, resourcesDir));
			lines.add("");
			
			if (!type.separateDists || distances.size() == 1) {
				lines.add("### "+type.name+" Variability Results");
				lines.add(topLink); lines.add("");
			}
			
			if (type.separateDists) {
				for (float distance : distances) {
					lines.addAll(generateVariabilityLines(type, distance, periods, topLink, resourcesDir, resultsTable));
					lines.add("");
				}
			} else {
				lines.addAll(generateVariabilityLines(type, null, periods, topLink, resourcesDir, resultsTable));
				lines.add("");
			}
		}
		
		// build summary table
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Type");
		table.addColumn("Notation");
		table.addColumn("Distance");
		for (double period : periods)
			table.addColumn(optionalDigitDF.format(period)+"s Std. Dev.");
		table.finalizeLine();
		for (VariabilityType type : VariabilityType.values()) {
			if (!resultsTable.containsRow(type))
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
				List<StdDevSet> stdDevs = resultsTable.get(type, distance);
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
		lines.addAll(summaryTableIndex, table.build());
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private List<String> generateVariabilityLines(VariabilityType type, Float distance, double[] periods, String topLink, File resourcesDir,
			Table<VariabilityType, Float, List<StdDevSet>> resultsTable) throws IOException {
		List<String> lines = new ArrayList<>();
		if (type.separateDists && distances.size() > 1) {
			lines.add("### "+distance+" km "+type.name+" Results");
			lines.add(topLink); lines.add("");
		}
		
//		lines.add("Standard deviations are computed from ");
		
		List<Site> mySites = new ArrayList<>();
		List<List<List<ResidualSet>>> myResiduals = new ArrayList<>();
		List<DiscretizedFunc> stdDevFuncs = null;
		
		if (type.separateSites) {
			Quantity[] constQuantities;
			
			if (type.separateDists)
				constQuantities = new Quantity[] { Quantity.DISTANCE, Quantity.SITE };
			else
				constQuantities = new Quantity[] { Quantity.SITE };
			
			stdDevFuncs = new ArrayList<>();
			
			List<List<ResidualSet>> totalPeriodResiduals = new ArrayList<>();
			for (int p=0; p<periods.length; p++)
				totalPeriodResiduals.add(new ArrayList<>());
			
			for (Site site : sites) {
				Object[] constValues;
				if (type.separateDists)
					constValues = new Object[] { distance, site };
				else
					constValues = new Object[] { site };
				
				List<List<ResidualSet>> residuals = calcResiduals(periods, constQuantities, constValues, type.groupQuantities);
				mySites.add(site);
				myResiduals.add(residuals);
				for (int i=0; i<periods.length; i++)
					totalPeriodResiduals.get(i).addAll(residuals.get(i));
			}
			mySites.add(null);
			myResiduals.add(totalPeriodResiduals);
		} else {
			Quantity[] constQuantities = null;
			Object[] constValues = null;
			
			if (type.separateDists) {
				constQuantities = new Quantity[] { Quantity.DISTANCE };
				constValues = new Object[] { distance };
			}
			List<List<ResidualSet>> residuals = calcResiduals(periods, constQuantities, constValues, type.groupQuantities);
			mySites.add(null);
			myResiduals.add(residuals);
		}
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		if (type.separateSites)
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
		List<StdDevSet> totStdDevs = new ArrayList<>();
		List<List<StdDevSet>> siteTotStdDevs = null;
		if (type.separateSites) {
			siteTotStdDevs = new ArrayList<>();
			for (int i=0; i<periods.length; i++)
				siteTotStdDevs.add(new ArrayList<>());
		}
		
		for (int i=0; i<mySites.size(); i++) {
			Site site = mySites.get(i);
			table.initNewLine();
			String emphasis = "";
			if (type.separateSites) {
				if (site == null) {
					table.addColumn("**ALL SITES**");
					emphasis = "**";
				} else {
					table.addColumn(site.getName());
				}
			}
			List<List<ResidualSet>> residuals = myResiduals.get(i);
			DiscretizedFunc stdDevFunc = new ArbitrarilyDiscretizedFunc();
			if (site == null) {
				stdDevFunc.setName("Total");
				totStdDevFunc = stdDevFunc;
			} else {
				stdDevFunc.setName(site.getName());
				stdDevFuncs.add(stdDevFunc);
			}
			List<StdDevSet> myStdDevs = new ArrayList<>();
			for (int p=0; p<periods.length; p++) {
				List<ResidualSet> periodResiduals = residuals.get(p);
				StdDevSet stat = new StdDevSet(periodResiduals);
				myStdDevs.add(stat);
				table.addColumn(""); // empty for period
				table.addColumn(emphasis+optionalDigitDF.format(stat.total)+emphasis);
				table.addColumn(emphasis+optionalDigitDF.format(stat.mean)+emphasis);
				table.addColumn(emphasis+optionalDigitDF.format(stat.median)+emphasis);
				table.addColumn(emphasis+"["+optionalDigitDF.format(stat.min)+" "+optionalDigitDF.format(stat.max)+"]"+emphasis);
				stdDevFunc.set(periods[p], stat.total);
			}
			if (site == null)
				totStdDevs = myStdDevs;
			else
				for (int p=0; p<periods.length; p++)
					siteTotStdDevs.get(p).add(myStdDevs.get(p));
			table.finalizeLine();
		}
		
		String distPrefix = distance == null ? "" : "_"+optionalDigitDF.format(distance)+"km";
		
		// plot it
		String prefix = type.prefix+distPrefix+"_std_dev";
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
			prefix = type.prefix+distPrefix+"_"+optionalDigitDF.format(periods[p])+"s_hist";
			plotStdDevsHistogram(resourcesDir, prefix, optionalDigitDF.format(periods[p])+"s "+type.name+" ("+type.symbol+")",
					periods[p], myResiduals.get(myResiduals.size()-1).get(p), totStdDevs.get(p), type.separateSites ? siteTotStdDevs.get(p): null);
			table.addColumn("!["+optionalDigitDF.format(periods[p])+"s](resources/"+prefix+".png)");
		}
		resultsTable.put(type, distance, totStdDevs);
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
	
	private class StdDevSet {
		private final double total;
		private final double mean;
		private final double median;
		private final double min;
		private final double max;
		
		public StdDevSet(List<ResidualSet> residuals) {
			double[] stdDevs = new double[residuals.size()];
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
	
	private List<List<ResidualSet>> calcResiduals(double[] periods, Quantity[] constQuantities, Object[] constValues, Quantity[] groupQuantities)
			throws IOException {
		List<RotationSpec> totalRotations;
		if (constQuantities != null && constQuantities.length > 0) {
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
		
		calcGroupedResiduals(totalRotations, periods, ret, groupQuantities);
		int count = 0;
		for (ResidualSet set : ret.get(new Random().nextInt(periods.length)))
			count += set.size();
		Preconditions.checkState(count == num, "Bad end index. Expected %s, have %s", num, count);
//		System.out.println(numThatVary+"/"+totNum+" have path variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
		
		return ret;
	}
	
	private void calcGroupedResiduals(List<RotationSpec> rotations, double[] periods, List<List<ResidualSet>> ret,
			Quantity[] groupQuantities) throws IOException {
		Preconditions.checkState(!rotations.isEmpty());
		if (groupQuantities.length == 0) {
//			System.out.println("Computing with index="+index);
			// we're at a set for which to calculate residuals, do it
			List<DiscretizedFunc> spectra = new ArrayList<>();
			for (RotationSpec rotation : rotations)
				spectra.add(prov.getRotD50(rotation.site, rotation, 0));
			for (int p=0; p<periods.length; p++) {
				double[] values = new double[spectra.size()];
				for (int i=0; i<values.length; i++)
					values[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
				ret.get(p).add(new ResidualSet(values));
			}
		} else {
//			System.out.println("have "+rotations.size()+" rotations before "+groupQuantities[0]+", index="+index);
			Quantity[] current = {groupQuantities[0]};
			Quantity[] downstream = Arrays.copyOfRange(groupQuantities, 1, groupQuantities.length);
			for (Object value : quantitiesMap.get(groupQuantities[0])) {
				List<RotationSpec> myRotations = RotatedRupVariabilityConfig.getRotationsForQuantities(
						rotations, current, new Object[] {value});
//				System.out.println("\t"+value+": "+myRotations.size());
				calcGroupedResiduals(myRotations, periods, ret, downstream);
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
			List<ResidualSet> residuals, StdDevSet totalStdDev, List<StdDevSet> siteStdDevs) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Preconditions.checkState(siteStdDevs == null || siteStdDevs.size() == siteColors.size());
		
		double delta = 0.05;
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, 1d, delta);
		hist.setName("Total Histogram");
		
		for (ResidualSet set : residuals)
			hist.add(hist.getClosestXIndex(set.residualStdDev), 1d);
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

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		
		double[] periods = {3d, 5d, 7.5, 10d};
		boolean doExample = true;
		
		System.out.println("Catalog: "+catalog.getName());
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
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-rotatedRups-")) {
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
		System.out.println("\tInput file: "+bbpZipFile.getName());
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site(RSQSimBBP_Config.VM));
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		Map<Scenario, RotatedRupVariabilityPageGen> pageGensMap = new HashMap<>();
		for (Scenario scenario : Scenario.values()) {
			File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			if (rotConfFile.exists()) {
				RotatedRupVariabilityConfig config = RotatedRupVariabilityConfig.loadCSV(catalog, rotConfFile, null, sites);
				
				BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario);
				
				RotatedRupVariabilityPageGen pageGen = new RotatedRupVariabilityPageGen(catalog, scenario, config, bbpLoader);
				
				pageGensMap.put(scenario, pageGen);
			}
		}
		
		Map<Scenario, RSQSimEvent> exampleEventsMap = null;
		if (doExample) {
			Map<Scenario, Integer> exampleIDs = new HashMap<>();
			for (RotatedRupVariabilityPageGen pageGen : pageGensMap.values())
				exampleIDs.put(pageGen.scenario, pageGen.eventIDs.get(0));
			System.out.println("Loading "+exampleIDs.size()+" example ruptures");
			List<RSQSimEvent> events = catalog.loader().byIDs(Ints.toArray(exampleIDs.values()));
			exampleEventsMap = new HashMap<>();
			for (Scenario scenario : exampleIDs.keySet()) {
				int exampleID = exampleIDs.get(scenario);
				for (RSQSimEvent event : events)
					if (event.getID() == exampleID)
						exampleEventsMap.put(scenario, event);
			}
		}
		
		for (Scenario scenario : pageGensMap.keySet()) {
			System.out.println("Doing scenario: "+scenario);
			
			RotatedRupVariabilityPageGen pageGen = pageGensMap.get(scenario);
			
			if (doExample) {
				RSQSimEvent exampleRup = exampleEventsMap.get(scenario);
				pageGen.setExampleRupture(exampleRup, pageGen.sites.get(0), 5);
			}
			
			File rotDir = new File(catalogOutputDir, "rotated_ruptures_"+scenario.getPrefix());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			methodSpecificLines.add("");
			methodSpecificLines.add("[RSQSim Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines);
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}
}
