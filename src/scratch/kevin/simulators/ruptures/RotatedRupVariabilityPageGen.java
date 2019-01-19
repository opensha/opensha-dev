package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.data.Range;
import org.opensha.commons.data.NamedComparator;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
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

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
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
	
	public void generatePage(File outputDir, double[] periods, List<String> headerLines) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new LinkedList<>();
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Rupture Rotation Parameters");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Events", eventIDs.size());
		table.addLine("Sites", sites.size());
		table.addLine("Source Rotation Azimuths", sourceAzimuths.size());
		table.addLine("Site-To-Source Path Azimuths", siteSourceAzimuths.size());
		table.addLine("Source-Site Distance[s]", Joiner.on(",").join(distances)+" km");
		table.addLine("Total # Simulations", +config.getRotations().size());
		lines.addAll(table.build());
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
		
		if (siteSourceAzimuths.size() > 1) {
			lines.add("## Path Variability");
			lines.add(topLink); lines.add("");
			
			lines.add("### Path Variability Methodology");
			lines.add(topLink); lines.add("");
			
			lines.add("Path variability is computed by spinning each rupture around each site, holding distance and relative rupture "
					+ "orientation constant. Here is an exmample with "+numExampleRotations+" rotations:");
			lines.add("");
			plotExample(resourcesDir, "example_path", distances.get(0), 1, numExampleRotations);
			lines.add("![Example](resources/example_path.png)");
			lines.add("");
			lines.add("This is done separately for each site, rupture, distance, and any rupture rotations (about its centroid). "
					+ "This calculation uses "+siteSourceAzimuths.size()+" rotations for each such site/rupture pair, where residuals "
					+ "are computed (relative to the mean intensity) for each path. A total path standard devaiation is then computed "
					+ "from these path residuals for all ruptures (including centroid rotations) for a given site and distance. The "
					+ "\"**ALL SITES**\" standard deviation is computed from all residuals across all sites for a given distance.");
			
			if (distances.size() == 1) {
				lines.add("### Path Variability Results");
				lines.add(topLink); lines.add("");
			}
			
			for (float distance : distances) {
				if (distances.size() > 1) {
					lines.add("### "+distance+" km Path Residuals");
					lines.add(topLink); lines.add("");
				}
				
				table = MarkdownUtils.tableBuilder();
				table.initNewLine();
				table.addColumn("Site");
				for (double period : periods) {
					String periodStr = optionalDigitDF.format(period)+"s";
					table.addColumn(periodStr+" Std. Dev.");
					table.addColumn(periodStr+" Range");
				}
				table.finalizeLine();
				StandardDeviation[] stdDevStat = new StandardDeviation[periods.length];
				double[] overallMin = new double[periods.length];
				double[] overallMax = new double[periods.length];
				for (int p=0; p<periods.length; p++) {
					stdDevStat[p] = new StandardDeviation();
					overallMin[p] = Double.POSITIVE_INFINITY;
					overallMax[p] = Double.NEGATIVE_INFINITY;
				}
				List<DiscretizedFunc> stdDevFuncs = new ArrayList<>(); 
				for (Site site : sites) {
					List<double[]> pathResiduals = calcPathResiduals(periods, site, distance);
					table.initNewLine();
					table.addColumn(site.getName());
					DiscretizedFunc stdDevFunc = new ArbitrarilyDiscretizedFunc(site.getName());
					stdDevFuncs.add(stdDevFunc);
					for (int i=0; i<periods.length; i++) {
						double[] periodResiduals = pathResiduals.get(i);
						double variance = StatUtils.variance(periodResiduals);
						double stdDev = Math.sqrt(variance);
						double min = StatUtils.min(periodResiduals);
						double max = StatUtils.max(periodResiduals);
						table.addColumn(optionalDigitDF.format(stdDev));
						table.addColumn("["+optionalDigitDF.format(min)+" "+optionalDigitDF.format(max)+"]");
						stdDevStat[i].incrementAll(periodResiduals);
						overallMin[i] = Double.min(min, overallMin[i]);
						overallMax[i] = Double.max(max, overallMax[i]);
						stdDevFunc.set(periods[i], stdDev);
					}
					table.finalizeLine();
				}
				table.initNewLine();
				table.addColumn("**ALL SITES**");
				DiscretizedFunc totStdDevFunc = new ArbitrarilyDiscretizedFunc("Total");
				for (int i=0; i<periods.length; i++) {
					table.addColumn("**"+optionalDigitDF.format(stdDevStat[i].getResult())+"**");
					table.addColumn("**["+optionalDigitDF.format(overallMin[i])+" "+optionalDigitDF.format(overallMax[i])+"]**");
					totStdDevFunc.set(periods[i], stdDevStat[i].getResult());
				}
				table.finalizeLine();
				// plot it
				plotSiteStdDevs(resourcesDir, "path_std_dev", "Path Variability", periods, stdDevFuncs, totStdDevFunc);
				lines.add("![Path Variability](resources/path_std_dev.png)");
				lines.add("");
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		if (sourceAzimuths.size() > 1) {
			lines.add("## Source Orientation (Directivity) Variability");
			lines.add(topLink); lines.add("");
			
			lines.add("### Source Orientation Methodology");
			lines.add(topLink); lines.add("");
			
			lines.add("Source orientation variability is computed by spinning each rupture around its centroid, holding distance and "
					+ "path constant. Here is an exmample with "+numExampleRotations+" rotations:");
			lines.add("");
			plotExample(resourcesDir, "example_source_orientation", distances.get(0), numExampleRotations, 1);
			lines.add("![Example](resources/example_source_orientation.png)");
			lines.add("");
			lines.add("This is done separately for each site, rupture, distance, and path. This calculation uses "+sourceAzimuths.size()
					+" centroid rotations for each rupture (for each site/distance/path), where residuals are computed (relative to the mean "
					+ "intensity) for each rupture orientation. A total rupture orientation standard devaiation is then computed "
					+ "from these residuals for all ruptures (including multiple paths) for a given site and distance. The "
					+ "\"**ALL SITES**\" standard deviation is computed from all residuals across all sites for a given distance.");
			
			if (distances.size() == 1) {
				lines.add("### Source Orientation Results");
				lines.add(topLink); lines.add("");
			}
			
//			for (float distance : distances) {
//				if (distances.size() > 1) {
//					lines.add("### "+distance+" km Source Orientation Residuals");
//					lines.add(topLink); lines.add("");
//				}
//				
//				table = MarkdownUtils.tableBuilder();
//				table.initNewLine();
//				table.addColumn("Site");
//				for (double period : periods) {
//					String periodStr = optionalDigitDF.format(period)+"s";
//					table.addColumn(periodStr+" Std. Dev.");
//					table.addColumn(periodStr+" Range");
//				}
//				table.finalizeLine();
//				StandardDeviation[] stdDevStat = new StandardDeviation[periods.length];
//				double[] overallMin = new double[periods.length];
//				double[] overallMax = new double[periods.length];
//				for (int p=0; p<periods.length; p++) {
//					stdDevStat[p] = new StandardDeviation();
//					overallMin[p] = Double.POSITIVE_INFINITY;
//					overallMax[p] = Double.NEGATIVE_INFINITY;
//				}
//				List<DiscretizedFunc> stdDevFuncs = new ArrayList<>(); 
//				for (Site site : sites) {
//					List<double[]> pathResiduals = calcPathResiduals(periods, site, distance);
//					table.initNewLine();
//					table.addColumn(site.getName());
//					DiscretizedFunc stdDevFunc = new ArbitrarilyDiscretizedFunc(site.getName());
//					stdDevFuncs.add(stdDevFunc);
//					for (int i=0; i<periods.length; i++) {
//						double[] periodResiduals = pathResiduals.get(i);
//						double variance = StatUtils.variance(periodResiduals);
//						double stdDev = Math.sqrt(variance);
//						double min = StatUtils.min(periodResiduals);
//						double max = StatUtils.max(periodResiduals);
//						table.addColumn(optionalDigitDF.format(stdDev));
//						table.addColumn("["+optionalDigitDF.format(min)+" "+optionalDigitDF.format(max)+"]");
//						stdDevStat[i].incrementAll(periodResiduals);
//						overallMin[i] = Double.min(min, overallMin[i]);
//						overallMax[i] = Double.max(max, overallMax[i]);
//						stdDevFunc.set(periods[i], stdDev);
//					}
//					table.finalizeLine();
//				}
//				table.initNewLine();
//				table.addColumn("**ALL SITES**");
//				DiscretizedFunc totStdDevFunc = new ArbitrarilyDiscretizedFunc("Total");
//				for (int i=0; i<periods.length; i++) {
//					table.addColumn("**"+optionalDigitDF.format(stdDevStat[i].getResult())+"**");
//					table.addColumn("**["+optionalDigitDF.format(overallMin[i])+" "+optionalDigitDF.format(overallMax[i])+"]**");
//					totStdDevFunc.set(periods[i], stdDevStat[i].getResult());
//				}
//				table.finalizeLine();
//				// plot it
//				plotSiteStdDevs(resourcesDir, "path_std_dev", "Path Variability", periods, stdDevFuncs, totStdDevFunc);
//				lines.add("![Path Variability](resources/path_std_dev.png)");
//				lines.add("");
//				lines.addAll(table.build());
//				lines.add("");
//			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	private static Float nullAsZero(Float value) {
		if (value == null)
			return 0f;
		return value;
	}
	
	private List<double[]> calcPathResiduals(double[] periods, Site site, float distance) throws IOException {
		int num = sourceAzimuths.size()*siteSourceAzimuths.size()*eventIDs.size();
		List<double[]> ret = new ArrayList<>();
		for (int i=0; i<periods.length; i++)
			ret.add(new double[num]);
		
		int index = 0;
		double[] workingArray = new double[siteSourceAzimuths.size()];
		System.out.println("Calculating path residuals with "+workingArray.length+" azimuths per site");
		int numThatVary = 0;
		int totNum = 0;
		for (int eventID : eventIDs) {
			for (float sourceAzimuth : sourceAzimuths) {
				List<RotationSpec> rotations = config.getSiteToSourceRotations(site, eventID, distance, sourceAzimuth);
				Preconditions.checkState(rotations.size() == workingArray.length,
						"Unexpected number of configurations for event %s, site %s, distance %s, az %s. Expected %s, have %s",
						eventID, site.getName(), distance, sourceAzimuth, workingArray.length, rotations.size());
				List<DiscretizedFunc> spectra = new ArrayList<>();
				for (RotationSpec rotation : rotations)
					spectra.add(prov.getRotD50(site, rotation, 0));
				boolean varies = false;
				for (int p=0; p<periods.length; p++) {
					for (int i=0; i<workingArray.length; i++)
						workingArray[i] = Math.log(spectra.get(i).getInterpolatedY(periods[p]));
					// now detrend to remove any source effects
					double mean = StatUtils.mean(workingArray);
					for (int i=0; i<workingArray.length; i++) {
						ret.get(p)[index+i] = workingArray[i]-mean;
						if (Math.abs(ret.get(p)[index+i]) > 1e-2)
							varies = true;
					}
				}
				index += workingArray.length;
				totNum++;
				if (varies)
					numThatVary++;
			}
		}
		Preconditions.checkState(index == num, "Bad end index. Expected %s, have %s", num, index);
		System.out.println(numThatVary+"/"+totNum+" have path variability ("+optionalDigitDF.format(100d*numThatVary/totNum)+" %)");
		
		return ret;
	}
	
	private void plotSiteStdDevs(File resourcesDir, String prefix, String title, double[] periods,
			List<DiscretizedFunc> siteStdDevFuncs, DiscretizedFunc totalStdDevFunc) throws IOException {
		Preconditions.checkState(siteStdDevFuncs.size() == siteColors.size());
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<siteStdDevFuncs.size(); i++) {
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
	
	private void plotExample(File resourcesDir, String prefix, double distance, int numSourceAz, int numSiteToSourceAz)
			throws IOException {
		List<Site> sites = new ArrayList<>();
		sites.add(exampleSite);
		List<RSQSimEvent> ruptures = new ArrayList<>();
		ruptures.add(exampleRupture);
		RotatedRupVariabilityConfig config = new RotatedRupVariabilityConfig(catalog, sites, ruptures, new double[] {distance},
				numSourceAz, numSiteToSourceAz);
		
		config.plotRotations(resourcesDir, prefix, config.getRotations());
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		
		double[] periods = {3d, 5d, 7.5, 10d};
		
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
		
		for (Scenario scenario : Scenario.values()) {
			File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			if (rotConfFile.exists()) {
				System.out.println("Doing scenario: "+scenario);
				
				RotatedRupVariabilityConfig config = RotatedRupVariabilityConfig.loadCSV(catalog, rotConfFile, null, sites);
				
				BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario);
				
				RotatedRupVariabilityPageGen pageGen = new RotatedRupVariabilityPageGen(catalog, scenario, config, bbpLoader);
				
				int exampleEventID = pageGen.eventIDs.get(0);
				System.out.println("Loading event "+exampleEventID+" for example plots");
				RSQSimEvent exampleRup = catalog.loader().byID(exampleEventID);
				pageGen.setExampleRupture(exampleRup, pageGen.sites.get(0), 5);
				
				File rotDir = new File(catalogOutputDir, "rotated_ruptures_"+scenario.getPrefix());
				Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
				
				List<String> lines = new ArrayList<>();
				lines.add("# "+catalog.getName()+" Rotated Rupture Variability, "+scenario.getShortName());
				lines.add("");
				lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
				
				pageGen.generatePage(rotDir, periods, lines);
			}
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}
}
