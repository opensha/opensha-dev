package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.laughTest.AbstractLaughTest;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter.CleanupFilter;
import scratch.UCERF3.inversion.laughTest.MinSectsPerParentFilter.ContinualFilter;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class RSQSimU3RupturePageGen {

	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException, RuntimeException {
		File catalogsBaseDir = new File("/data/kevin/simulators/catalogs");
		File mainOutputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");

//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(catalogsBaseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(catalogsBaseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_3062.instance(catalogsBaseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(catalogsBaseDir);
		
		boolean rebuildSol = false;
		
		File catalogDir = catalog.getCatalogDir();
		
		File catalogOutputDir = new File(mainOutputDir, catalogDir.getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File multiFaultDir = new File(catalogOutputDir, "multi_fault");
		Preconditions.checkState(multiFaultDir.exists() || multiFaultDir.mkdir());
		
		File resourcesDir = new File(multiFaultDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		File fssDir = new File(catalogDir, "fss");
		Preconditions.checkState(fssDir.exists() || fssDir.mkdir());
		
		String catalogName = catalog.getName();
		String catalogType = "RSQSim";
		String catalogTypeFileName = "rsqsim";
		
		double minMag = 6.5;
		int skipYears = 5000;
		boolean includeNumSects = false;
		
		List<String> lines = new ArrayList<>();
		
		FaultBasedMapGen.MAP_LABEL_SIZE = 24;
		FaultBasedMapGen.MAP_LABEL_TICK_SIZE = 20;
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.FAULT_THICKNESS = 4d;
		
		double minFractForInclusion = catalog.getMinSubSectFractForInclusion();
		
		// header
		lines.add("# Multi Fault Rupture Comparisons");
		lines.add("");
		lines.add("*Subsections participate in a rupture if at least "
				+(float)(minFractForInclusion*100d)+" % of its area ruptures*");
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		FaultModels fm = catalog.getFaultModel();
		DeformationModels dm = catalog.getDeformationModel();
		
		String catParams = "m"+(float)minMag+"_skip"+skipYears+"_sectArea"+(float)minFractForInclusion;
		File solFile = new File(fssDir, "rsqsim_sol_"+catParams+".zip");
		FaultSystemSolution sol;
		if (solFile.exists() && !rebuildSol) {
			System.out.println("Loading solution from: "+solFile.getAbsolutePath());
			sol = FaultSystemIO.loadSol(solFile);
		} else {
			System.out.println("Loading events from: "+catalogDir.getAbsolutePath());
			Loader loader = catalog.loader().minMag(minMag).skipYears(skipYears);
			sol = catalog.buildSolution(loader, minMag);
			System.out.println("Writing solution to: "+solFile.getAbsolutePath());
			FaultSystemIO.writeSol(sol, solFile);
		}
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		System.out.println(rupSet.getNumRuptures()+" unique ruptures");
		
		LaughTestFilter filter = LaughTestFilter.getDefault();
		if (!includeNumSects)
			filter.setMinNumSectInRup(0);
		
		File scratchDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/");
		DeformationModelFetcher fetch = new DeformationModelFetcher(fm, dm, scratchDir, 0.1);
		
		List<? extends FaultSection> datas = fetch.getSubSectionList();
		
		Map<IDPairing, Double> distances = fetch.getSubSectionDistanceMap(1000d);
		Map<IDPairing, Double> azimuths = fetch.getSubSectionAzimuthMap(distances.keySet());
		Map<Integer, Double> rakesMap = new HashMap<Integer, Double>();
		for (FaultSection data : rupSet.getFaultSectionDataList())
			rakesMap.put(data.getSectionId(), data.getAveRake());
		boolean applyGarlockPintoMtnFix = true;
		
		CoulombRates coulombRates = null;
		if (filter.getCoulombFilter() != null) {
			try {
				coulombRates = CoulombRates.loadUCERF3CoulombRates(fm);
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
		
		// we want actual catalog rupture counts before binning into U3 style ruptures
		// find the smallest rate, which will be 1/catLen, then numRups = solRate/minRate
		double minRate = StatUtils.min(sol.getRateForAllRups());
		
		List<List<Integer>> sectionConnectionsListList = SectionClusterList.computeCloseSubSectionsListList(
				datas, distances, filter.getMaxJumpDist(), coulombRates);
		
		List<AbstractLaughTest> tests = filter.buildLaughTests(azimuths, distances, rakesMap, coulombRates, applyGarlockPintoMtnFix,
				sectionConnectionsListList, rupSet.getFaultSectionDataList());
		
		// doesn't come with jump dist filter by default (it is included explicitly in generation for UCERF3)
		tests.add(0, new JumpDistFilter(distances, 5d));
		
		if (includeNumSects) {
			// replace the separate min sects per parent filters with a single one
			tests.add(1, new CombinedMinSectsFilter(removeByClass(tests, ContinualFilter.class),
				removeByClass(tests, CleanupFilter.class)));
		}
		
		Color[] colors = { Color.DARK_GRAY, Color.RED, Color.BLUE, Color.GREEN.darker(), Color.CYAN, Color.ORANGE };
		
		int allPassCount = 0;
		int[] failCounts = new int[tests.size()];
		int[] onlyFailCounts = new int[tests.size()];
		int[] erredCounts = new int[tests.size()];
		
		int tot = 0;
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			int numCatalogOccurances = (int)Math.round(sol.getRateForRup(r)/minRate);
			tot += numCatalogOccurances;
			Preconditions.checkState(numCatalogOccurances >= 1);
			List<? extends FaultSection> rupture = rupSet.getFaultSectionDataForRupture(r);
			boolean allPass = true;
			int onlyFailureIndex = -1;
			for (int t=0; t<tests.size(); t++) {
				AbstractLaughTest test = tests.get(t);
				boolean subPass;
				try {
					subPass = test.doesRupturePass(rupture);
				} catch (Exception e) {
					if (erredCounts[t] == 0) {
						System.err.println("First exception for "+test.getName()+":");
						e.printStackTrace();
					}
					erredCounts[t] += numCatalogOccurances;
					subPass = true; // do not fail on error
				}
				if (!subPass && allPass) {
					// this is the first failure
					onlyFailureIndex = t;
				} else if (!subPass) {
					// failed more than 1
					onlyFailureIndex = -1;
				}
				allPass = subPass && allPass;
				if (!subPass)
					failCounts[t] += numCatalogOccurances;
			}
			if (allPass)
				allPassCount += numCatalogOccurances;
			if (onlyFailureIndex >= 0)
				onlyFailCounts[onlyFailureIndex] += numCatalogOccurances;
		}
		System.out.println("Passed all filters: "+countStats(allPassCount, tot));
		for (int t=0; t<tests.size(); t++) {
			System.out.println(tests.get(t).getName());
			System.out.println("\tFailed: "+countStats(failCounts[t], tot));
			System.out.println("\tOnly Failure: "+countStats(onlyFailCounts[t], tot));
			System.out.println("\tErred: "+countStats(erredCounts[t], tot));
		}
		
		// now plot
		double dx = 1d;
		double buffer = 0.2*dx;
		double deltaEachSide = (dx - buffer)/2d;
		float thickness = 80f;
		double maxY = 50;
		
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		
		List<PlotElement> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(new DefaultXY_DataSet(new double[] {0d, 1d}, new double[] {0d, 0d}));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, Color.WHITE));
		
		List<XYAnnotation> anns = new ArrayList<>();
		
		for (int i=0; i<tests.size(); i++) {
			double x = i*dx + 0.5*dx;
			double percentFailed = 100d*failCounts[i]/tot;
			double percentOnly = 100d*onlyFailCounts[i]/tot;
			double percentErred = 100d*erredCounts[i]/tot;
			
			Color c = colors[i % colors.length];
			
//			funcs.add(vertLine(x, 0, percentFailed));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, c));
			anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentFailed, c));
			
			if (percentOnly > 0) {
//				funcs.add(vertLine(x, 0, percentOnly));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, darker(c)));
				anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentOnly, darker(c)));
			}
			
			String title = tests.get(i).getShortName();
			
			if (percentErred > 0) {
//				funcs.add(vertLine(x, percentFailed, percentFailed + percentErred));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, Color.LIGHT_GRAY));
				anns.add(emptyBox(x-deltaEachSide, percentFailed, x+deltaEachSide, percentFailed + percentErred,
						PlotLineType.DASHED, Color.LIGHT_GRAY, 2f));
				title += "*";
			}
			
			XYTextAnnotation ann = new XYTextAnnotation(title, x, maxY*0.95);
			ann.setTextAnchor(TextAnchor.TOP_CENTER);
			ann.setPaint(c);
			ann.setFont(font);
			
			anns.add(ann);
			
			ann = new XYTextAnnotation(percentDF.format(percentFailed/100d), x, percentFailed+0.6);
			ann.setTextAnchor(TextAnchor.BOTTOM_CENTER);
			ann.setPaint(Color.BLACK);
			ann.setFont(font);
			
			anns.add(ann);
		}
		
		XYTextAnnotation ann = new XYTextAnnotation(
				percentDF.format((double)allPassCount/tot)+" passed all", dx*0.25, maxY*0.88);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setPaint(Color.BLACK);
		ann.setFont(font);
		
		anns.add(ann);
		
		String title = "Rupture Plausibility Filters, M≥"+(float)minMag+", SectArea≥"+(float)minFractForInclusion;
		PlotSpec spec = new PlotSpec(funcs, chars, title, " ", "Percent Failed");
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(resourcesDir, "filters_"+catParams).getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, false, new Range(0, tests.size()*dx), new Range(0, maxY));
		gp.getXAxis().setTickLabelsVisible(false);
//		gp.getXAxis().setvisi
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		
		lines.add("## Plausibility Filter Comparisons");
//		lines.add(topLink); lines.add("");
		lines.add("");
		
		lines.add("### Rupture Failure Percentages");
		lines.add(topLink); lines.add("");
		lines.add("");
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/filters_"+catParams+".png)");
		
		// now jumps
		lines.add("## 1km Jump Count");
		lines.add(topLink); lines.add("");
		lines.add("");
		System.out.println("Plotting num jumps");
		RSQSimRupJumpCompare.plotFixedJumpDist(catalog.getU3CompareSol(), distances, sol, catalogName, minMag, 1d, resourcesDir);
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/jumps_1.0km.png)");
		lines.add("");
		
		// cumulant mag
		System.out.println("Plotting cumulant mag");
		plotCumulantMags(catalog.getU3CompareSol(), sol, catalogType, resourcesDir);
		
		lines.add("## Cumulant Magnitude");
		lines.add(topLink); lines.add("");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3", "Difference");
		File rsPlot = new File(resourcesDir, "mag_cumulant_medians_"+catalogType.toLowerCase()+".png");
		Preconditions.checkState(rsPlot.exists());
		File u3Plot = new File(resourcesDir, "mag_cumulant_medians_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		File diffPlot = new File(resourcesDir, "mag_cumulant_medians_diff.png");
		Preconditions.checkState(diffPlot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")",
				"![Difference]("+resourcesDir.getName()+"/"+diffPlot.getName()+")");
		rsPlot = new File(resourcesDir, "mag_cumulant_iqr_"+catalogType.toLowerCase()+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "mag_cumulant_iqr_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		diffPlot = new File(resourcesDir, "mag_cumulant_iqr_diff.png");
		Preconditions.checkState(diffPlot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")",
				"![Difference]("+resourcesDir.getName()+"/"+diffPlot.getName()+")");
		lines.addAll(table.build());
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine("![Median Scatter]("+resourcesDir.getName()+"/mag_cumulant_medians_scatter.png)",
				"![IQR Scatter]("+resourcesDir.getName()+"/mag_cumulant_iqr_scatter.png)");
		lines.addAll(table.build());
		lines.add("");
		
		System.out.println("Plotting section connections");
		Map<IDPairing, Double> rsConnections = calcConnecitons(sol, distances, 100d);
		System.out.println("Detected "+rsConnections.size()+" connections for "+prefix);
		FaultSystemSolution u3Sol = catalog.getU3CompareSol();
		Map<IDPairing, Double> u3Connections = calcConnecitons(u3Sol, distances, 25d);
		System.out.println("Detected "+u3Connections.size()+" U3 connections");
//		Color rsColor = Color.RED.darker();
//		Color u3Color = Color.BLUE.darker();
		Color rsColor = Color.RED;
		Color u3Color = Color.BLUE;
		Region fullReg = new CaliforniaRegions.RELM_TESTING();
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_"+catalogTypeFileName,
				catalogName+" Connectivity", rsConnections.keySet(), rsColor, fullReg);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_ucerf3",
				"UCERF3 Connectivity", u3Connections.keySet(), u3Color, fullReg);
		Map<IDPairing, Double> rsUniqueConnections = new HashMap<>(rsConnections);
		for (IDPairing pair : u3Connections.keySet())
			if (rsUniqueConnections.containsKey(pair))
				rsUniqueConnections.remove(pair);
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName,
				catalogName+" Unique Connectivity", rsUniqueConnections.keySet(), rsColor, fullReg);
		Map<IDPairing, Double> u3UniqueConnections = new HashMap<>(u3Connections);
		for (IDPairing pair : rsConnections.keySet())
			if (u3UniqueConnections.containsKey(pair))
				u3UniqueConnections.remove(pair);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_unique_ucerf3",
				"UCERF3 Unique Connectivity", u3UniqueConnections.keySet(), u3Color, fullReg);
		
		double maxConnDist = 0d;
		for (IDPairing pair : rsConnections.keySet())
			maxConnDist = Math.max(maxConnDist, distances.get(pair));
		for (IDPairing pair : u3Connections.keySet())
			maxConnDist = Math.max(maxConnDist, distances.get(pair));
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_"+catalogTypeFileName,
				catalogName+" Connectivity", rsConnections, rsUniqueConnections, maxConnDist,
				distances, rsColor, false);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_ucerf3",
				"UCERF3 Connectivity", u3Connections, u3UniqueConnections, maxConnDist,
				distances, u3Color, false);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName,
				catalogName+" Connectivity", rsConnections, rsUniqueConnections, maxConnDist,
				distances, rsColor, true);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_ucerf3",
				"UCERF3 Connectivity", u3Connections, u3UniqueConnections, maxConnDist,
				distances, u3Color, true);
		
		lines.add("## Fault Section Connections");
		lines.add(topLink); lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3");
		rsPlot = new File(resourcesDir, "sect_connectivity_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		rsPlot = new File(resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_unique_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		rsPlot = new File(resourcesDir, "sect_connectivity_hist_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_hist_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		rsPlot = new File(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_hist_rates_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("### Zoomed Connection Maps");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3");

		Region northZoomReg = new Region(new Location(35.5, -120.5), new Location(39.5, -123));
		northZoomReg.setName("zoom_north");
		Region southZoomReg = new Region(new Location(32.5, -115), new Location(35.5, -121));
		southZoomReg.setName("south_north");
		Region[] zoomRegions = { northZoomReg, southZoomReg };
		for (Region zoomReg : zoomRegions) {
			String zoomPrefix = "sect_connectivity_"+zoomReg.getName();
			plotConnectivityLines(sol.getRupSet(), resourcesDir, zoomPrefix+"_"+catalogTypeFileName,
					catalogName+" Connectivity", rsConnections.keySet(), rsColor, zoomReg);
			plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, zoomPrefix+"_ucerf3",
					"UCERF3 Connectivity", u3Connections.keySet(), u3Color, zoomReg);
			rsPlot = new File(resourcesDir, zoomPrefix+"_"+catalogTypeFileName+".png");
			Preconditions.checkState(rsPlot.exists());
			u3Plot = new File(resourcesDir, zoomPrefix+"_ucerf3.png");
			Preconditions.checkState(u3Plot.exists());
			table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
					"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		System.out.println("Plotting connectivity clusters");
		plotConnectivity(catalog.getU3CompareSol().getRupSet(), resourcesDir, "connectivity_ucerf3", "UCERF3 Connectivity");
		plotConnectivity(sol.getRupSet(), resourcesDir, "connectivity_"+catalogTypeFileName, catalogName+" Connectivity");
		
		lines.add("## Fault Connectivity Clusters");
		lines.add(topLink); lines.add("");
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3");
		rsPlot = new File(resourcesDir, "connectivity_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "connectivity_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+")");
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, multiFaultDir);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
	}
	
	private static Color darker(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
//		r += (255-r)/2;
//		g += (255-g)/2;
//		b += (255-b)/2;
		r /= 2;
		g /= 2;
		b /= 2;
		return new Color(r, g, b);
	}
	
	private static DefaultXY_DataSet vertLine(double x, double y0, double y1) {
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.set(x, y0);
		line.set(x, y1);
		return line;
	}
	
	private static XYBoxAnnotation filledBox(double x0, double y0, double x1, double y1, Color c) {
		XYBoxAnnotation ann = new XYBoxAnnotation(x0, y0, x1, y1, null, null, c);
		return ann;
	}
	
	private static XYBoxAnnotation emptyBox(double x0, double y0, double x1, double y1,
			PlotLineType lineType, Color c, float thickness) {
		Stroke stroke = lineType.buildStroke(thickness);
		XYBoxAnnotation ann = new XYBoxAnnotation(x0, y0, x1, y1, stroke, c, null);
		return ann;
	}
	
	private static final DecimalFormat percentDF = new DecimalFormat("0.00%");
	private static String countStats(int count, int tot) {
		return count+"/"+tot+" ("+percentDF.format((double)count/(double)tot)+")";
	}
	
	private static class JumpDistFilter extends AbstractLaughTest {
		
		private Map<IDPairing, Double> distances;
		private double maxJumpDist;
		
		public JumpDistFilter(Map<IDPairing, Double> distances, double maxJumpDist) {
			this.distances = distances;
			this.maxJumpDist = maxJumpDist;
		}

		@Override
		public String getShortName() {
			return "JumpDist";
		}

		@Override
		public String getName() {
			return "Maximum Jump Dist";
		}

		@Override
		public boolean doesLastSectionPass(List<? extends FaultSection> rupture, List<IDPairing> pairings,
				List<Integer> junctionIndexes) {
			if (junctionIndexes.isEmpty())
				return true;
			IDPairing pair = pairings.get(junctionIndexes.size()-1);
			return distances.get(pair) <= maxJumpDist;
		}

		@Override
		public boolean isContinueOnFaulure() {
			return false;
		}

		@Override
		public boolean isApplyJunctionsOnly() {
			return true;
		}
	}
	
	private static class CombinedMinSectsFilter extends AbstractLaughTest {
		
		private ContinualFilter continualFilter;
		private CleanupFilter cleanupFilter;

		public CombinedMinSectsFilter(MinSectsPerParentFilter.ContinualFilter continualFilter,
				MinSectsPerParentFilter.CleanupFilter cleanupFilter) {
			this.continualFilter = continualFilter;
			this.cleanupFilter = cleanupFilter;
		}

		@Override
		public String getShortName() {
			return "SectsPerParent";
		}

		@Override
		public String getName() {
			return "Min Sects Per Parent";
		}

		@Override
		public boolean doesLastSectionPass(List<? extends FaultSection> rupture, List<IDPairing> pairings,
				List<Integer> junctionIndexes) {
			boolean passContinual = continualFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
			boolean passCleanup = cleanupFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
//			if (!junctionIndexes.isEmpty() && junctionIndexes.get(junctionIndexes.size()-1) == pairings.size()-1)
//				return cleanupFilter.doesLastSectionPass(rupture, pairings, junctionIndexes);
			int i = rupture.size()-1;
			boolean junction = i > 0 &&
					rupture.get(i).getParentSectionId() != rupture.get(i-1).getParentSectionId();
			return passContinual && (!junction || passCleanup);
		}

		@Override
		public boolean isContinueOnFaulure() {
			return true;
		}

		@Override
		public boolean isApplyJunctionsOnly() {
			return false;
		}
		
	}
	
	private static <E extends AbstractLaughTest> E removeByClass(List<AbstractLaughTest> tests, Class<E> clazz) {
		for (int i=tests.size(); --i>=0;) {
			AbstractLaughTest test = tests.get(i);
			if (clazz.isInstance(test))
				return (E)tests.remove(i);
		}
		throw new IllegalStateException();
	}
	
	static void plotCumulantMags(FaultSystemSolution u3Sol, FaultSystemSolution rsSol,
			String catalogName, File outputDir) throws IOException, GMT_MapException, RuntimeException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Sect Index", "Sect Name", catalogName+" Median", catalogName+" IQR",
				"UCERF3 Median", "UCERF3 IQR");
		
		FaultSystemRupSet u3RupSet = u3Sol.getRupSet();
		FaultSystemRupSet rsRupSet = rsSol.getRupSet();
		
		Preconditions.checkState(u3RupSet.getNumSections() == rsRupSet.getNumSections());
		
		List<Double> rsMedians = new ArrayList<>();
		List<Double> u3Medians = new ArrayList<>();
		
		List<Double> rsIQRs = new ArrayList<>();
		List<Double> u3IQRs = new ArrayList<>();
		
		DefaultXY_DataSet medianScatter = new DefaultXY_DataSet();
		DefaultXY_DataSet iqrScatter = new DefaultXY_DataSet();
		
		for (int s=0; s<u3RupSet.getNumSections(); s++) {
			List<String> line = new ArrayList<>();
			line.add(s+"");
			line.add(u3RupSet.getFaultSectionData(s).getName());
			EvenlyDiscretizedFunc rsFunc = calcCumulantMagFunc(rsSol, s);
			double rsMedian, rsIQR;
			if (rsFunc == null) {
				rsMedian = Double.NaN;
				rsIQR = Double.NaN;
			} else {
				rsMedian = rsFunc.getFirstInterpolatedX(0.5);
				rsIQR = rsFunc.getFirstInterpolatedX(0.75) - rsFunc.getFirstInterpolatedX(0.25);
			}
			rsMedians.add(rsMedian);
			rsIQRs.add(rsIQR);
			line.add(rsMedian+"");
			line.add(rsIQR+"");
			EvenlyDiscretizedFunc u3Func = calcCumulantMagFunc(u3Sol, s);
			double u3Median, u3IQR;
			if (u3Func == null) {
				u3Median = Double.NaN;
				u3IQR = Double.NaN;
			} else {
				u3Median = u3Func.getFirstInterpolatedX(0.5);
				u3IQR = u3Func.getFirstInterpolatedX(0.75) - u3Func.getFirstInterpolatedX(0.25);
				if (!Double.isNaN(rsMedian)) {
					medianScatter.set(u3Median, rsMedian);
					iqrScatter.set(u3IQR, rsIQR);
				}
			}
			u3Medians.add(u3Median);
			u3IQRs.add(u3IQR);
			line.add(u3Median+"");
			line.add(u3IQR+"");
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, "mag_cumulant_medians.csv"));
		
		List<LocationList> faults = new ArrayList<>();
		for (FaultSection fault : u3RupSet.getFaultSectionDataList())
			faults.add(fault.getFaultTrace());
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(6d,  8.5d);
		Region reg = new CaliforniaRegions.RELM_TESTING();
		FaultBasedMapGen.makeFaultPlot(cpt, faults, Doubles.toArray(u3Medians), reg, outputDir,
				"mag_cumulant_medians_ucerf3", false, false, "UCERF3 Mag Cumulant Median");
		FaultBasedMapGen.makeFaultPlot(cpt, faults, Doubles.toArray(rsMedians), reg, outputDir,
				"mag_cumulant_medians_"+catalogName.toLowerCase(), false, false, catalogName+" Mag Cumulant Median");
		
//		CPT diffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1d, 1d);
		CPT diffCPT = new CPT(-1, 1d,
				new Color(0, 0, 140), new Color(0, 60, 200 ), new Color(0, 120, 255),
				Color.WHITE,
				new Color(255, 120, 0), new Color(200, 60, 0), new Color(140, 0, 0));
		diffCPT.setBelowMinColor(diffCPT.getMinColor());
		diffCPT.setAboveMaxColor(diffCPT.getMaxColor());
		double[] diffVals = new double[faults.size()];
		
		for (int i=0; i<diffVals.length; i++)
			diffVals[i] = rsMedians.get(i) - u3Medians.get(i);
		
		FaultBasedMapGen.makeFaultPlot(diffCPT, faults, diffVals, reg, outputDir,
				"mag_cumulant_medians_diff", false, true, catalogName+"-U3 Mag Cumulant Median");
		
		// IQRs
		double maxIQR = Math.ceil(2*Math.max(iqrScatter.getMaxX(), iqrScatter.getMaxY()))*0.5d;
//		CPT iqrCPT = new CPT(0d, 1d, new Color(40, 0, 0), new Color(80, 0, 0), 	new Color(140, 0, 0),
//				new Color(200, 60, 0), new Color(255, 120, 0), new Color(255, 200, 0),
//				new Color(255, 225, 0), Color.WHITE);
		CPT iqrCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().rescale(0d, 0.7);
		Color lastIQRCOlor = new Color(255, 255, 200);
		iqrCPT.add(new CPTVal(iqrCPT.getMaxValue(), iqrCPT.getMaxColor(), 1f, lastIQRCOlor));
		iqrCPT.setAboveMaxColor(lastIQRCOlor);
		iqrCPT.setNanColor(Color.GRAY);
		FaultBasedMapGen.makeFaultPlot(iqrCPT, faults, Doubles.toArray(u3IQRs), reg, outputDir,
				"mag_cumulant_iqr_ucerf3", false, false, "UCERF3 Mag Cumulant IQR");
		FaultBasedMapGen.makeFaultPlot(iqrCPT, faults, Doubles.toArray(rsIQRs), reg, outputDir,
				"mag_cumulant_iqr_"+catalogName.toLowerCase(), false, false, catalogName+" Mag Cumulant IQR");
		
		diffCPT = diffCPT.rescale(-0.5, 0.5);
		diffVals = new double[faults.size()];
		
		for (int i=0; i<diffVals.length; i++)
			diffVals[i] = rsIQRs.get(i) - u3IQRs.get(i);
		
		FaultBasedMapGen.makeFaultPlot(diffCPT, faults, diffVals, reg, outputDir,
				"mag_cumulant_iqr_diff", false, true, catalogName+"-U3 Mag Cumulant IQR");
		
		// scatters
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		double minMedian = 0.5*Math.floor(2d*Math.min(medianScatter.getMinX(), medianScatter.getMinY()));
		double maxMedian = 0.5*Math.ceil(2d*Math.max(medianScatter.getMaxX(), medianScatter.getMaxY()));
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(minMedian, minMedian);
		oneToOne.set(maxMedian, maxMedian);
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		funcs.add(medianScatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Median Cumulant Mag Scatter", "UCERF3", catalogName);
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		DecimalFormat magDF = new DecimalFormat("0.0");
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 22);
		double[] rsMedianArray = getNoNans(rsMedians);
		double rsMeanMedian = StatUtils.mean(rsMedianArray);
		double rsMedianMedian = DataUtils.median(rsMedianArray);
		XYTextAnnotation ann = new XYTextAnnotation(
				"  "+catalogName+": mean="+magDF.format(rsMeanMedian)+", mdn.="+magDF.format(rsMedianMedian),
				minMedian, minMedian + 0.95*(maxMedian-minMedian));
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(annFont);
		anns.add(ann);
		double[] u3MedianArray = getNoNans(u3Medians);
		double u3MeanMedian = StatUtils.mean(u3MedianArray);
		double u3MedianMedian = DataUtils.median(u3MedianArray);
		ann = new XYTextAnnotation(
				"  UCERF3: mean="+magDF.format(u3MeanMedian)+", mdn.="+magDF.format(u3MedianMedian),
				minMedian, minMedian + 0.9*(maxMedian-minMedian));
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(annFont);
		anns.add(ann);
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		File pngFile = new File(outputDir, "mag_cumulant_medians_scatter.png");
		
		gp.drawGraphPanel(spec, false, false, new Range(minMedian, maxMedian), new Range(minMedian, maxMedian));
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(pngFile.getAbsolutePath());
		
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		
		oneToOne = new DefaultXY_DataSet();
		oneToOne.set(0d, 0d);
		oneToOne.set(maxIQR, maxIQR);
		
		funcs.add(oneToOne);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		
		funcs.add(iqrScatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		spec = new PlotSpec(funcs, chars, "Cumulant Mag IQR Scatter", "UCERF3", catalogName);
		
		anns = new ArrayList<>();
		DecimalFormat iqrDF = new DecimalFormat("0.00");
		double[] rsIQRArray = getNoNans(rsIQRs);
		double rsMeanIQR = StatUtils.mean(rsIQRArray);
		double rsMedianIQR = DataUtils.median(rsIQRArray);
		ann = new XYTextAnnotation(
				"  "+catalogName+": mean="+iqrDF.format(rsMeanIQR)+", mdn.="+iqrDF.format(rsMedianIQR),
				0d, 0.95*maxIQR);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(annFont);
		anns.add(ann);
		double[] u3IQRArray = getNoNans(u3IQRs);
		double u3MeanIQR = StatUtils.mean(u3IQRArray);
		double u3MedianIQR = DataUtils.median(u3IQRArray);
		ann = new XYTextAnnotation(
				"  UCERF3: mean="+iqrDF.format(u3MeanIQR)+", mdn.="+iqrDF.format(u3MedianIQR),
				0d, 0.9*maxIQR);
		ann.setTextAnchor(TextAnchor.TOP_LEFT);
		ann.setFont(annFont);
		anns.add(ann);
		spec.setPlotAnnotations(anns);
		
		pngFile = new File(outputDir, "mag_cumulant_iqr_scatter.png");
		
		gp.drawGraphPanel(spec, false, false, new Range(0, maxIQR), new Range(0, maxIQR));
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(pngFile.getAbsolutePath());
	}
	
	private static double[] getNoNans(List<Double> vals) {
		List<Double> ret = new ArrayList<>();
		for (double val : vals)
			if (Double.isFinite(val))
				ret.add(val);
		return Doubles.toArray(ret);
	}
	
	private static EvenlyDiscretizedFunc calcCumulantMagFunc(FaultSystemSolution sol, int s) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(5d, 9d, (int)((9d-5d)/0.01) + 1);
		for (int r : sol.getRupSet().getRupturesForSection(s)) {
			double mag = sol.getRupSet().getMagForRup(r);
			int i = func.getClosestXIndex(mag);
			for (int x=i; x<func.size(); x++)
				func.add(x, MagUtils.magToMoment(mag));
		}
		if (func.calcSumOfY_Vals() == 0)
			return null;
		func.scale(1d/func.getMaxY());
		return func;
	}
	
	static void plotConnectivity(FaultSystemRupSet rupSet, File outputDir, String prefix, String title)
			throws IOException, GMT_MapException, RuntimeException {
		List<HashSet<Integer>> clusters = new ArrayList<>();
		Map<Integer, Integer> sectIndexToClusterIndexMap = new HashMap<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++)
			processClusterRecursive(rupSet, s, clusters.size(), clusters, sectIndexToClusterIndexMap);
		
		System.out.println("Detected "+clusters.size()+" clusters for "+prefix);
		
		CPT refCPT = GMT_CPT_Files.MAX_SPECTRUM.instance();
		refCPT = refCPT.rescale(0d, 1d);
		// list of values for each discrete color, initially sorted from first color to last
		List<Double> colorValues = Lists.newArrayList();
//		for (CPTVal cptVal : refCPT)
//			colorValues.add((double)cptVal.start);
//		colorValues.add((double)refCPT.get(refCPT.size()-1).end);
		for (double v=0; v<=1d; v+=0.1)
			colorValues.add(v);
		// now sorted from last color to first
		Collections.reverse(colorValues);
		refCPT.setNanColor(Color.GRAY);
		
		// sort from smallest to largest
		List<Integer> sizes = Lists.newArrayList();
		for (HashSet<Integer> cluster : clusters)
			sizes.add(cluster.size());
		clusters = ComparablePairing.getSortedData(sizes, clusters);
		// now reverse, largest to smallest
		Collections.reverse(clusters);
		
		if (clusters.size() > colorValues.size())
			clusters = clusters.subList(0, colorValues.size());
		
		List<LocationList> faults = new ArrayList<>();
		List<Double> values = new ArrayList<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			faults.add(rupSet.getFaultSectionData(s).getFaultTrace());
			double val = Double.NaN;
			for (int i=0; i<clusters.size(); i++) {
				if (clusters.get(i).contains(s))
					val = colorValues.get(i);
			}
			values.add(val);
		}
		
		Region reg = new CaliforniaRegions.RELM_TESTING();
		FaultBasedMapGen.makeFaultPlot(refCPT, faults, Doubles.toArray(values), reg, outputDir,
				prefix, false, false, title+" ("+colorValues.size()+" largest)");
	}
	
	private static void processClusterRecursive(FaultSystemRupSet rupSet, int sect, int clusterIndex, List<HashSet<Integer>> clusters,
			Map<Integer, Integer> sectIndexToClusterIndexMap) {
		if (sectIndexToClusterIndexMap.containsKey(sect))
			// we've already done this one
			return;
		if (clusters.size() == clusterIndex)
			clusters.add(new HashSet<>());
		clusters.get(clusterIndex).add(sect);
		sectIndexToClusterIndexMap.put(sect, clusterIndex);
		for (int r : rupSet.getRupturesForSection(sect)) {
			for (int sect2 : rupSet.getSectionsIndicesForRup(r)) {
				processClusterRecursive(rupSet, sect2, clusterIndex, clusters, sectIndexToClusterIndexMap);
			}
		}
	}
	
	static void plotConnectivityLines(FaultSystemRupSet rupSet, File outputDir, String prefix, String title,
			Set<IDPairing> connections, Color connectedColor, Region reg) throws IOException {
		Color faultColor = Color.DARK_GRAY;
		Color faultOutlineColor = Color.LIGHT_GRAY;
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		XY_DataSet[] outlines = PoliticalBoundariesData.loadCAOutlines();
		PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(PlotLineType.SOLID, (float)1d, Color.GRAY);
		
		for (XY_DataSet outline : outlines) {
			funcs.add(outline);
			chars.add(outlineChar);
		}
		
		List<Location> middles = new ArrayList<>();
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			RuptureSurface surf = sect.getFaultSurface(1d);
			
			XY_DataSet trace = new DefaultXY_DataSet();
			for (Location loc : surf.getEvenlyDiscritizedUpperEdge())
				trace.set(loc.getLongitude(), loc.getLatitude());
			
			if (sect.getAveDip() != 90d) {
				XY_DataSet outline = new DefaultXY_DataSet();
				LocationList perimeter = surf.getPerimeter();
				for (Location loc : perimeter)
					outline.set(loc.getLongitude(), loc.getLatitude());
				Location first = perimeter.first();
				outline.set(first.getLongitude(), first.getLatitude());
				
				funcs.add(0, outline);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, faultOutlineColor));
			}
			
			middles.add(GriddedSurfaceUtils.getSurfaceMiddleLoc(surf));
			
			if (s == 0)
				trace.setName("Fault Sections");
			
			funcs.add(trace);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, faultColor));
		}
		
		boolean first = true;
		for (IDPairing connection : connections) {
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			
			if (first) {
				xy.setName("Connections");
				first = false;
			}
			
			Location loc1 = middles.get(connection.getID1());
			Location loc2 = middles.get(connection.getID2());
			
			xy.set(loc1.getLongitude(), loc1.getLatitude());
			xy.set(loc2.getLongitude(), loc2.getLatitude());
			
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, connectedColor));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Longitude", "Latitude");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(reg.getMinLon(), reg.getMaxLon());
		Range yRange = new Range(reg.getMinLat(), reg.getMaxLat());
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		double tick = 2d;
		TickUnits tus = new TickUnits();
		TickUnit tu = new NumberTickUnit(tick);
		tus.add(tu);
		gp.getXAxis().setStandardTickUnits(tus);
		gp.getYAxis().setStandardTickUnits(tus);
		
		File file = new File(outputDir, prefix);
		double aspectRatio = yRange.getLength() / xRange.getLength();
		gp.getChartPanel().setSize(800, 200 + (int)(600d*aspectRatio));
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	static void plotConnectivityHistogram(File outputDir, String prefix, String title,
			Map<IDPairing, Double> connections, Map<IDPairing, Double> uniqueConnections,
			double maxDist,	Map<IDPairing, Double> distances, Color connectedColor, boolean rateWeighted)
					throws IOException {
		double delta = 1d;
//		if (maxDist > 90)
//			delta = 5d;
//		else if (maxDist > 40)
//			delta = 2;
//		else if (maxDist > 20)
//			delta = 1d;
//		else
//			delta = 0.5d;
		
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, maxDist, delta);
		hist.setName("All Connections");
		HistogramFunction uniqueHist = HistogramFunction.getEncompassingHistogram(0d, maxDist, delta);
		uniqueHist.setName("Unique To Model");
		
		double myMax = 0d;
		double mean = 0d;
		double sumWeights = 0d;
		double meanAbove = 0d;
		double sumWeightsAbove = 0d;
		
		for (IDPairing pair : connections.keySet()) {
			double dist = distances.get(pair);
			double weight = rateWeighted ? connections.get(pair) : 1d;
			
			myMax = Math.max(myMax, dist);
			mean += dist*weight;
			sumWeights += weight;
			if (dist >= 0.1) {
				meanAbove += dist*weight;
				sumWeightsAbove += weight;
			}
			
			int xIndex = hist.getClosestXIndex(dist);
			hist.add(xIndex, weight);
			if (uniqueConnections.containsKey(pair))
				uniqueHist.add(xIndex, weight);
		}

		mean /= sumWeights;
		meanAbove /= sumWeightsAbove;
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Color uniqueColor = new Color(connectedColor.getRed()/4,
				connectedColor.getGreen()/4, connectedColor.getBlue()/4);
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, connectedColor));
		
		funcs.add(uniqueHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, uniqueColor));
		
		String yAxisLabel = rateWeighted ? "Annual Rate" : "Count";
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Jump Distance (km)", yAxisLabel);
		spec.setLegendVisible(true);
		
		Range xRange = new Range(0d, maxDist);
		Range yRange = new Range(0d, 1.05*hist.getMaxY());
		
		DecimalFormat distDF = new DecimalFormat("0.0");
		double annX = 0.975*maxDist;
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		
		double annYScalar = 0.975;
		double annYDelta = 0.05;
		
		double annY = annYScalar*yRange.getUpperBound();
		XYTextAnnotation maxAnn = new XYTextAnnotation(
				"Max: "+distDF.format(myMax), annX, annY);
		maxAnn.setFont(annFont);
		maxAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(maxAnn);
		
		annYScalar -= annYDelta;
		annY = annYScalar*yRange.getUpperBound();
		XYTextAnnotation meanAnn = new XYTextAnnotation(
				"Mean: "+distDF.format(mean), annX, annY);
		meanAnn.setFont(annFont);
		meanAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(meanAnn);
		
		annYScalar -= annYDelta;
		annY = annYScalar*yRange.getUpperBound();
		if (rateWeighted) {
			XYTextAnnotation rateAnn = new XYTextAnnotation(
					"Total Rate: "+distDF.format(sumWeights), annX, annY);
			rateAnn.setFont(annFont);
			rateAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			spec.addPlotAnnotation(rateAnn);
		} else {
			XYTextAnnotation countAnn = new XYTextAnnotation(
					"Total Count: "+(int)sumWeights, annX, annY);
			countAnn.setFont(annFont);
			countAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			spec.addPlotAnnotation(countAnn);
		}
		
		annYScalar -= annYDelta;
		annY = annYScalar*yRange.getUpperBound();
		XYTextAnnotation meanAboveAnn = new XYTextAnnotation(
				"Mean >0.1: "+distDF.format(meanAbove), annX, annY);
		meanAboveAnn.setFont(annFont);
		meanAboveAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		spec.addPlotAnnotation(meanAboveAnn);
		
		annYScalar -= annYDelta;
		annY = annYScalar*yRange.getUpperBound();
		if (rateWeighted) {
			XYTextAnnotation rateAnn = new XYTextAnnotation(
					"Total Rate >0.1: "+distDF.format(sumWeightsAbove), annX, annY);
			rateAnn.setFont(annFont);
			rateAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			spec.addPlotAnnotation(rateAnn);
		} else {
			XYTextAnnotation countAnn = new XYTextAnnotation(
					"Total Count >0.1: "+(int)sumWeightsAbove, annX, annY);
			countAnn.setFont(annFont);
			countAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
			spec.addPlotAnnotation(countAnn);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 650);
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	private static Map<IDPairing, Double> calcConnecitons(FaultSystemSolution sol,
			Map<IDPairing, Double> distCache, double maxPossilbeJump) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		RupSetConnectionSearch search = new RupSetConnectionSearch(
				rupSet, distCache, maxPossilbeJump, RupSetConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		
		List<Future<HashSet<IDPairing>>> futures = new ArrayList<>();
		for (int r=0; r<rupSet.getNumRuptures(); r++)
			futures.add(exec.submit(new ConnectionCalc(search, r)));
		
		Map<IDPairing, Double> connections = new HashMap<>();
		
		for (int r=0; r<futures.size(); r++) {
			double rate = sol.getRateForRup(r);
			if (r % 1000 == 0)
				System.out.println("Calculating for rupture "+r+"/"+rupSet.getNumRuptures()
					+" ("+connections.size()+" connections found so far)");
			Future<HashSet<IDPairing>> future = futures.get(r);
			try {
				HashSet<IDPairing> pairings = future.get();
				for (IDPairing pair : pairings) {
					Double prevRate = connections.get(pair);
					if (prevRate == null)
						prevRate = 0d;
					connections.put(pair, prevRate + rate);
				}
			} catch (InterruptedException | ExecutionException e) {
				exec.shutdown();
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		System.out.println("Found "+connections.size()+" total connections");
		
		exec.shutdown();
		
		return connections;
	}
	
	private static class ConnectionCalc implements Callable<HashSet<IDPairing>> {
		
		private RupSetConnectionSearch search;
		private int rupIndex;

		public ConnectionCalc(RupSetConnectionSearch search, int rupIndex) {
			this.search = search;
			this.rupIndex = rupIndex;
		}

		@Override
		public HashSet<IDPairing> call() throws Exception {
			return search.calcConnections(rupIndex);
		}
		
	}

}
