package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
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
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRuptureBuilder;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.CoulombJunctionFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.CumulativeAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpDistFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpAzimuthChangeFilter.AzimuthCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.MinSectsPerParentFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.SingleClusterPerParentFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.SplayCountFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.TotalAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.U3CompatibleCumulativeRakeChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.DistCutoffClosestSectClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.UCERF3ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.UCERF3SectionConnectionStrategy;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester.TestType;
import scratch.UCERF3.inversion.laughTest.PlausibilityResult;
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
		System.out.println(rupSet.getNumRuptures()+" ruptures");
		
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(
				rupSet.getFaultSectionDataList());
		File distAzCacheFile = new File(fssDir, "dist_az_cache.csv");
		if (distAzCacheFile.exists())
			distAzCalc.loadCacheFile(distAzCacheFile);
		RupSetConnectionSearch rsConnSearch = new RupSetConnectionSearch(rupSet, distAzCalc,
				new DistCutoffClosestSectClusterConnectionStrategy(1000d),
				RupSetConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		
		System.out.println("Building ClusterRuptures for RSQSim");
		List<ClusterRupture> rsClusterRups = buildClusterRups(rupSet, rsConnSearch);
		Map<Jump, List<Integer>> rsJumpsToRupsMap = new HashMap<>();
		Map<Jump, Double> rsJumps = getJumps(sol, rsClusterRups, rsJumpsToRupsMap);
		HashSet<UniqueRupture> rsUniqueRups = new HashSet<>();
		for (ClusterRupture rup : rsClusterRups)
			rsUniqueRups.add(rup.unique);
		System.out.println("Found "+rsUniqueRups.size()+" unique ruptures");
		System.out.println("Detected "+rsJumps.size()+" RSQSim connections");
		
		System.out.println("Building ClusterRuptures for UCERF3");
		FaultSystemSolution u3Sol = catalog.getComparisonSolution();
		FaultSystemRupSet u3RupSet = u3Sol.getRupSet();
		RupSetConnectionSearch u3ConnSearch = new RupSetConnectionSearch(u3RupSet, distAzCalc,
				new DistCutoffClosestSectClusterConnectionStrategy(15d),
				RupSetConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		List<ClusterRupture> u3ClusterRups = new ArrayList<>();
		for (int r=0; r<u3RupSet.getNumRuptures(); r++)
			u3ClusterRups.add(ClusterRupture.forOrderedSingleStrandRupture(
					u3RupSet.getFaultSectionDataForRupture(r), distAzCalc));
		Map<Jump, List<Integer>> u3JumpsToRupsMap = new HashMap<>();
		Map<Jump, Double> u3Jumps = getJumps(u3Sol, u3ClusterRups, u3JumpsToRupsMap);
		System.out.println("Detected "+u3Jumps.size()+" U3 connections");
		
		// we want actual catalog rupture counts before binning into U3 style ruptures
		// find the smallest rate, which will be 1/catLen, then numRups = solRate/minRate
		double minRate = StatUtils.min(sol.getRateForAllRups());
		
		List<PlausibilityFilter> filters = new ArrayList<>();
		double maxJumpDist = 5d;
		// these are sort of by-construction filters
		filters.add(new JumpDistFilter(maxJumpDist));
		CoulombRates coulombRates = CoulombRates.loadUCERF3CoulombRates(fm);
		List<FaultSubsectionCluster> clusters = ClusterRuptureBuilder.buildClusters(
				rupSet.getFaultSectionDataList(),
				new UCERF3ClusterConnectionStrategy(maxJumpDist, coulombRates), distAzCalc);
		filters.add(new MinSectsPerParentFilter(2, true, clusters));
		filters.add(new SingleClusterPerParentFilter());
		filters.add(new SplayCountFilter(0));
		
		// these are rupture properties themselves
		AzimuthCalc u3AzCalc = new JumpAzimuthChangeFilter.UCERF3LeftLateralFlipAzimuthCalc(distAzCalc);
		filters.add(new FailOnCantEvalAzFilter(new JumpAzimuthChangeFilter(u3AzCalc, 60f), false));
		filters.add(new FailOnCantEvalAzFilter(new TotalAzimuthChangeFilter(u3AzCalc, 60f, true, true), true));
		filters.add(new CumulativeAzimuthChangeFilter(
				new JumpAzimuthChangeFilter.SimpleAzimuthCalc(distAzCalc), 560f));
//		filters.add(new CumulativeRakeChangeFilter(180f));
//		filters.add(new JumpCumulativeRakeChangeFilter(180f));
		filters.add(new U3CompatibleCumulativeRakeChangeFilter(180d));
		CoulombRatesTester coulombTester = new CoulombRatesTester(
				TestType.COULOMB_STRESS, 0.04, 0.04, 1.25d, true, true);
		filters.add(new CoulombJunctionFilter(coulombTester, coulombRates));
		
		Color[] colors = { Color.DARK_GRAY, new Color(102, 51, 0), Color.RED, Color.BLUE,
				Color.GREEN.darker(), Color.CYAN, Color.PINK, Color.ORANGE.darker(), Color.MAGENTA };
		
		int allPassCount = 0;
		int[] failCounts = new int[filters.size()];
		int[] failCanContinueCounts = new int[filters.size()];
		int[] onlyFailCounts = new int[filters.size()];
		int[] erredCounts = new int[filters.size()];
		
		int tot = 0;
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			int numCatalogOccurances = (int)Math.round(sol.getRateForRup(r)/minRate);
			tot += numCatalogOccurances;
			Preconditions.checkState(numCatalogOccurances >= 1);
			ClusterRupture rupture = rsClusterRups.get(r);
			boolean allPass = true;
			int onlyFailureIndex = -1;
			for (int t=0; t<filters.size(); t++) {
				PlausibilityFilter test = filters.get(t);
				PlausibilityResult result;
				boolean subPass;
				try {
					result = test.apply(rupture, false);
					if (result == PlausibilityResult.FAIL_FUTURE_POSSIBLE)
						failCanContinueCounts[t]++;
					subPass = result.isPass();
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
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Filter", "Failed", "Only Failure", "Erred");
		System.out.println("Passed all filters: "+countStats(allPassCount, tot));
		for (int t=0; t<filters.size(); t++) {
			table.initNewLine();
			table.addColumn("**"+filters.get(t).getName()+"**");
			table.addColumn(countStats(failCounts[t], tot));
			table.addColumn(countStats(onlyFailCounts[t], tot));
			table.addColumn(countStats(erredCounts[t], tot));
			table.finalizeLine();
		}
		
		// now plot
		double dx = 1d;
		double buffer = 0.2*dx;
		double deltaEachSide = (dx - buffer)/2d;
		float thickness = 80f;
		double maxY = 50;

		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 22);
		Font allFont = new Font(Font.SANS_SERIF, Font.BOLD, 26);
		
		List<PlotElement> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(new DefaultXY_DataSet(new double[] {0d, 1d}, new double[] {0d, 0d}));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, Color.WHITE));
		
		List<XYAnnotation> anns = new ArrayList<>();
		
		double topRowY = maxY*0.95;
		double secondRowY = maxY*0.91;
		double thirdRowY = maxY*0.85;
		
		for (int i=0; i<filters.size(); i++) {
			double x = i*dx + 0.5*dx;
			double percentFailed = 100d*failCounts[i]/tot;
			double percentOnly = 100d*onlyFailCounts[i]/tot;
			double percentErred = 100d*erredCounts[i]/tot;
			
			Color c = colors[i % colors.length];
			
			String title = filters.get(i).getShortName();
			
			if (percentErred > 0) {
//				funcs.add(vertLine(x, percentFailed, percentFailed + percentErred));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, Color.LIGHT_GRAY));
				anns.add(emptyBox(x-deltaEachSide, 0d, x+deltaEachSide, percentFailed + percentErred,
						PlotLineType.DASHED, Color.LIGHT_GRAY, 2f));
				title += "*";
			}
			
//			funcs.add(vertLine(x, 0, percentFailed));
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, c));
			anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentFailed, c));
			
			if (percentOnly > 0) {
//				funcs.add(vertLine(x, 0, percentOnly));
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, thickness, darker(c)));
				anns.add(filledBox(x-deltaEachSide, 0, x+deltaEachSide, percentOnly, darker(c)));
			}
			
			XYTextAnnotation ann = new XYTextAnnotation(title, x, i % 2 == 0 ? secondRowY : thirdRowY);
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
		
		Range xRange = new Range(-0.15*dx, (filters.size()+0.15)*dx);
		
		XYTextAnnotation ann = new XYTextAnnotation(
				percentDF.format((double)allPassCount/tot)+" passed all", xRange.getCentralValue(), topRowY);
		ann.setTextAnchor(TextAnchor.CENTER);
		ann.setPaint(Color.BLACK);
		ann.setFont(allFont);
		
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
		
		gp.drawGraphPanel(spec, false, false, xRange, new Range(0, maxY));
		gp.getXAxis().setTickLabelsVisible(false);
//		gp.getXAxis().setvisi
		gp.getChartPanel().setSize(1200, 600);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		
		lines.add("## Plausibility Filter Comparisons");
		lines.add(topLink); lines.add("");
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/filters_"+catParams+".png)");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		float[] maxJumpDists = { 0.1f, 1f, 3f };
		
		// now jumps
		lines.add("## Jump Counts Over Distance");
		lines.add(topLink); lines.add("");
		for (float jumpDist : maxJumpDists) {
			lines.add("");
			System.out.println("Plotting num jumps");
			File plotFile = plotFixedJumpDist(u3Sol, u3ClusterRups, sol, rsClusterRups, distAzCalc, catalogName,
					minMag, jumpDist, resourcesDir);
			lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/"+plotFile.getName()+")");
		}
		lines.add("");
		
		// now azimuths
		List<RakeType> rakeTypes = new ArrayList<>();
		rakeTypes.add(null);
		for (RakeType type : RakeType.values())
			rakeTypes.add(type);
		lines.add("## Jump Azimuths");
		lines.add(topLink); lines.add("");
		
		Table<RakeType, RakeType, List<Double>> rsRakeAzTable = calcJumpAzimuths(rsClusterRups, distAzCalc);
		Table<RakeType, RakeType, List<Double>> u3RakeAzTable = calcJumpAzimuths(u3ClusterRups, distAzCalc);
		
		for (RakeType sourceType : rakeTypes) {
			if (sourceType == null) {
				prefix = "jump_az_any";
				title = "Jumps from Any";
				lines.add("### Jump Azimuths From Any");
			} else {
				prefix = "jump_az_"+sourceType.prefix;
				title = "Jumps from "+sourceType.name;
				lines.add("### Jump Azimuths From "+sourceType.name);
			}
			
			System.out.println("Plotting "+title);

			lines.add(topLink); lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.addLine(catalogName, "UCERF3");
			
			table.initNewLine();
			File plotFile = plotJumpAzimuths(sourceType, rakeTypes, rsRakeAzTable,
					resourcesDir, prefix, title);
			table.addColumn("!["+title+"](resources/"+plotFile.getName()+")");
			plotFile = plotJumpAzimuths(sourceType, rakeTypes, u3RakeAzTable,
					resourcesDir, "u3_"+prefix, title);
			table.addColumn("!["+title+"](resources/"+plotFile.getName()+")");
			table.finalizeLine();
			lines.addAll(table.build());
			lines.add("");
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			
			for (RakeType destType : rakeTypes) {
				String myPrefix = prefix+"_";
				String myTitle = title+" to ";
				if (destType == null) {
					myPrefix += "any";
					myTitle += "Any";
				} else {
					myPrefix += destType.prefix;
					myTitle += destType.name;
				}
				
				plotFile = plotJumpAzimuthsRadial(sourceType, destType, rsRakeAzTable,
						resourcesDir, myPrefix, myTitle);
				table.addColumn("!["+title+"](resources/"+plotFile.getName()+")");
			}
			table.finalizeLine();
			lines.addAll(table.wrap(3, 0).build());
			lines.add("");
		}
		
		System.out.println("Plotting section connections");
//		Color rsColor = Color.RED.darker();
//		Color u3Color = Color.BLUE.darker();
		Color rsColor = Color.RED;
		Color u3Color = Color.BLUE;
		Color bothColor = Color.GREEN.darker();
		Region fullReg = new CaliforniaRegions.RELM_TESTING();
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps.keySet(), rsColor, fullReg, 800);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_ucerf3",
				"UCERF3 Connectivity", u3Jumps.keySet(), u3Color, fullReg, 800);
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_"+catalogTypeFileName+"_hires",
				catalogName+" Connectivity", rsJumps.keySet(), rsColor, fullReg, 3000);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_ucerf3_hires",
				"UCERF3 Connectivity", u3Jumps.keySet(), u3Color, fullReg, 3000);
		Map<Jump, Double> rsUniqueJumps = new HashMap<>(rsJumps);
		Set<Jump> commonJumps = new HashSet<>();
		for (Jump jump : u3Jumps.keySet()) {
			if (rsUniqueJumps.containsKey(jump)) {
				rsUniqueJumps.remove(jump);
				commonJumps.add(jump);
			}
		}
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName,
				catalogName+" Unique Connectivity", rsUniqueJumps.keySet(), rsColor, fullReg, 800);
		plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName+"_hires",
				catalogName+" Unique Connectivity", rsUniqueJumps.keySet(), rsColor, fullReg, 3000);
		Map<Jump, Double> u3UniqueJumps = new HashMap<>(u3Jumps);
		for (Jump jump : rsJumps.keySet())
			if (u3UniqueJumps.containsKey(jump))
				u3UniqueJumps.remove(jump);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_unique_ucerf3",
				"UCERF3 Unique Connectivity", u3UniqueJumps.keySet(), u3Color, fullReg, 800);
		plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_unique_ucerf3_hires",
				"UCERF3 Unique Connectivity", u3UniqueJumps.keySet(), u3Color, fullReg, 3000);
		
		double maxConnDist = 0d;
		for (Jump jump : rsJumps.keySet())
			maxConnDist = Math.max(maxConnDist, jump.distance);
		for (Jump jump : u3Jumps.keySet())
			maxConnDist = Math.max(maxConnDist, jump.distance);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps, rsUniqueJumps, maxConnDist,
				rsColor, false);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_ucerf3",
				"UCERF3 Connectivity", u3Jumps, u3UniqueJumps, maxConnDist,
				u3Color, false);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps, rsUniqueJumps, maxConnDist,
				rsColor, true);
		plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_ucerf3",
				"UCERF3 Connectivity", u3Jumps, u3UniqueJumps, maxConnDist,
				u3Color, true);
		
		lines.add("## Fault Section Connections");
		lines.add(topLink); lines.add("");
		
		List<Set<Jump>> connectionsList = new ArrayList<>();
		List<Color> connectedColors = new ArrayList<>();
		List<String> connNames = new ArrayList<>();
		
		connectionsList.add(rsUniqueJumps.keySet());
		connectedColors.add(rsColor);
		connNames.add("RSQSim Only");
		
		connectionsList.add(u3UniqueJumps.keySet());
		connectedColors.add(u3Color);
		connNames.add("UCERF3 Only");
		
		connectionsList.add(commonJumps);
		connectedColors.add(bothColor);
		connNames.add("Common Connections");
		
		String combConnPrefix = "sect_connectivity_combined";
		plotConnectivityLines(rupSet, resourcesDir, combConnPrefix, "Combined Connectivity",
				connectionsList, connectedColors, connNames, fullReg, 800);
		plotConnectivityLines(rupSet, resourcesDir, combConnPrefix+"_hires", "Combined Connectivity",
				connectionsList, connectedColors, connNames, fullReg, 3000);
		lines.add("![Combined]("+resourcesDir.getName()+"/"+combConnPrefix+".png)");
		lines.add("");
		lines.add("[View high resolution]("+resourcesDir.getName()+"/"+combConnPrefix+"_hires.png)");
		lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3");
		File rsPlot = new File(resourcesDir, "sect_connectivity_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		File u3Plot = new File(resourcesDir, "sect_connectivity_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+") "
				+ "[View high resolution]("+resourcesDir.getName()
				+"/"+rsPlot.getName().replace(".png", "_hires.png")+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+") "
				+ "[View high resolution]("+resourcesDir.getName()
				+"/"+u3Plot.getName().replace(".png", "_hires.png")+")");
		rsPlot = new File(resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_unique_ucerf3.png");
		Preconditions.checkState(u3Plot.exists());
		table.addLine("!["+catalogName+"]("+resourcesDir.getName()+"/"+rsPlot.getName()+") "
				+ "[View high resolution]("+resourcesDir.getName()
				+"/"+rsPlot.getName().replace(".png", "_hires.png")+")",
				"![UCERF3]("+resourcesDir.getName()+"/"+u3Plot.getName()+") "
				+ "[View high resolution]("+resourcesDir.getName()
				+"/"+u3Plot.getName().replace(".png", "_hires.png")+")");
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
		
		System.out.println("Plotting connection examples");
		lines.add("### Unique Connection Example Ruptures");
		lines.add(topLink); lines.add("");
		
		lines.add("**RSQSim Ruptures with Unique Connections**");
		int maxRups = 20;
		int maxCols = 5;
		table = plotConnRupExamples(rsConnSearch, rsUniqueJumps.keySet(),
				rsJumpsToRupsMap, maxRups, maxCols, resourcesDir, "rs_conn_example");
		lines.add("");
		if (table == null)
			lines.add("*N/A*");
		else
			lines.addAll(table.build());
		lines.add("");
		lines.add("**UCERF3 Ruptures with Unique Connections**");
		table = plotConnRupExamples(u3ConnSearch, u3UniqueJumps.keySet(),
				u3JumpsToRupsMap, maxRups, maxCols, resourcesDir, "u3_conn_example");
		lines.add("");
		if (table == null)
			lines.add("*N/A*");
		else
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
		
		// cumulant mag
		System.out.println("Plotting cumulant mag");
		plotCumulantMags(catalog.getU3CompareSol(), sol, catalogType, resourcesDir);

		lines.add("## Cumulant Magnitude");
		lines.add(topLink); lines.add("");
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine(catalogName, "UCERF3", "Difference");
		rsPlot = new File(resourcesDir, "mag_cumulant_medians_"+catalogType.toLowerCase()+".png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "mag_cumulant_medians_ucerf3.png");
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
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, multiFaultDir);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(mainOutputDir);
		
		distAzCalc.writeCacheFile(distAzCacheFile);
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
	
	private static class FailOnCantEvalAzFilter implements PlausibilityFilter {
		
		private PlausibilityFilter filter;
		private boolean endsOnly;

		private FailOnCantEvalAzFilter(PlausibilityFilter filter, boolean endsOnly) {
			this.filter = filter;
			this.endsOnly = endsOnly;
		}

		@Override
		public String getShortName() {
			return filter.getShortName();
		}

		@Override
		public String getName() {
			return filter.getName();
		}

		@Override
		public PlausibilityResult apply(ClusterRupture rupture, boolean verbose) {
			PlausibilityResult result = filter.apply(rupture, verbose);
			if (result.isPass() || !result.canContinue())
				return result;
			if (endsOnly) {
				checkStrandEndsRecursive(rupture);
			} else {
				for (Jump jump : rupture.getJumpsIterable()) {
					if (rupture.sectDescendantsMap.get(jump.toSection).isEmpty())
						throw new IllegalStateException("Jump to single section. Rupture:\n"+rupture);
				}
			}
			return result;
		}
		
		private void checkStrandEndsRecursive(ClusterRupture rupture) {
			for (Jump jump : rupture.internalJumps)
				if (jump.toCluster == rupture.clusters[rupture.clusters.length-1])
					if (rupture.sectDescendantsMap.get(jump.toSection).isEmpty())
						throw new IllegalStateException("Jump to single section. Rupture:\n"+rupture);
			for (ClusterRupture splay : rupture.splays.values())
				checkStrandEndsRecursive(splay);
		}

		@Override
		public PlausibilityResult testJump(ClusterRupture rupture, Jump newJump, boolean verbose) {
			return apply(rupture.take(newJump), verbose);
		}
		
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
			Set<Jump> connections, Color connectedColor, Region reg, int width) throws IOException {
		List<Set<Jump>> connectionsList = new ArrayList<>();
		List<Color> connectedColors = new ArrayList<>();
		List<String> connNames = new ArrayList<>();
		
		connectionsList.add(connections);
		connectedColors.add(connectedColor);
		connNames.add("Connections");
		
		plotConnectivityLines(rupSet, outputDir, prefix, title, connectionsList, connectedColors, connNames, reg, width);
	}
	
	static void plotConnectivityLines(FaultSystemRupSet rupSet, File outputDir, String prefix, String title,
			List<Set<Jump>> connectionsList, List<Color> connectedColors, List<String> connNames,
			Region reg, int width) throws IOException {
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
		
		for (int i=0; i<connectionsList.size(); i++) {
			Set<Jump> connections = connectionsList.get(i);
			Color connectedColor = connectedColors.get(i);
			String connName = connNames.get(i);
			
			boolean first = true;
			for (Jump connection : connections) {
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				
				if (first) {
					xy.setName(connName);
					first = false;
				}
				
				Location loc1 = middles.get(connection.fromSection.getSectionId());
				Location loc2 = middles.get(connection.toSection.getSectionId());
				
				xy.set(loc1.getLongitude(), loc1.getLatitude());
				xy.set(loc2.getLongitude(), loc2.getLatitude());
				
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, connectedColor));
			}
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
		gp.getChartPanel().setSize(width, 200 + (int)((width-200d)*aspectRatio));
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	static void plotConnectivityHistogram(File outputDir, String prefix, String title,
			Map<Jump, Double> connections, Map<Jump, Double> uniqueConnections,
			double maxDist, Color connectedColor, boolean rateWeighted)
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
		
		for (Jump pair : connections.keySet()) {
			double dist = pair.distance;
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
	
	private static List<ClusterRupture> buildClusterRups(FaultSystemRupSet rupSet,
			RupSetConnectionSearch search) {
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		
		List<Future<ClusterRupture>> futures = new ArrayList<>();
		for (int r=0; r<rupSet.getNumRuptures(); r++)
			futures.add(exec.submit(new ClusterRupCalc(search, r)));
		
		List<ClusterRupture> ruptures = new ArrayList<>();
		
		for (int r=0; r<futures.size(); r++) {
			if (r % 1000 == 0)
				System.out.println("Calculating for rupture "+r+"/"+rupSet.getNumRuptures());
			Future<ClusterRupture> future = futures.get(r);
			try {
//				HashSet<IDPairing> pairings = future.get();
//				for (IDPairing pair : pairings) {
//					Double prevRate = connections.get(pair);
//					if (prevRate == null)
//						prevRate = 0d;
//					connections.put(pair, prevRate + rate);
//					List<Integer> prevRups = pairToRupsMap.get(pair);
//					if (prevRups == null) {
//						prevRups = new ArrayList<>();
//						pairToRupsMap.put(pair, prevRups);
//					}
//					prevRups.add(r);
//				}
				ruptures.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				exec.shutdown();
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		System.out.println("Built "+ruptures.size()+" ruptures");
		
		exec.shutdown();
		Preconditions.checkState(ruptures.size() == rupSet.getNumRuptures());
		
		return ruptures;
	}
	
	private static Map<Jump, Double> getJumps(FaultSystemSolution sol, List<ClusterRupture> ruptures,
			Map<Jump, List<Integer>> jumpToRupsMap) {
		Map<Jump, Double> jumpRateMap = new HashMap<>();
		for (int r=0; r<ruptures.size(); r++) {
			double rate = sol.getRateForRup(r);
			ClusterRupture rupture = ruptures.get(r);
			for (Jump jump : rupture.getJumpsIterable()) {
				if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
					jump = jump.reverse();
				Double prevRate = jumpRateMap.get(jump);
				if (prevRate == null)
					prevRate = 0d;
				jumpRateMap.put(jump, prevRate + rate);
				if (jumpToRupsMap != null) {
					List<Integer> prevRups = jumpToRupsMap.get(jump);
					if (prevRups == null) {
						prevRups = new ArrayList<>();
						jumpToRupsMap.put(jump, prevRups);
					}
					prevRups.add(r);
				}
			}
		}
		return jumpRateMap;
	}
	
	private static class ClusterRupCalc implements Callable<ClusterRupture> {
		
		private RupSetConnectionSearch search;
		private int rupIndex;

		public ClusterRupCalc(RupSetConnectionSearch search, int rupIndex) {
			this.search = search;
			this.rupIndex = rupIndex;
		}

		@Override
		public ClusterRupture call() throws Exception {
			return search.buildClusterRupture(rupIndex, false);
		}
		
	}
	
	private static TableBuilder plotConnRupExamples(RupSetConnectionSearch search, Set<Jump> pairings,
			Map<Jump, List<Integer>> pairRupsMap, int maxRups, int maxCols,
			File resourcesDir, String prefix) throws IOException {
		List<Jump> sortedPairings = new ArrayList<>(pairings);
		Collections.sort(sortedPairings, Jump.id_comparator);
		
		Random r = new Random(sortedPairings.size()*maxRups);
		Collections.shuffle(sortedPairings, r);
		
		int possibleRups = 0;
		for (Jump pair : pairings)
			possibleRups += pairRupsMap.get(pair).size();
		if (possibleRups < maxRups)
			maxRups = possibleRups;
		if (maxRups == 0)
			return null;
		
		int indInPairing = 0;
		List<Integer> rupsToPlot = new ArrayList<>();
		while (rupsToPlot.size() < maxRups) {
			for (Jump pair : sortedPairings) {
				List<Integer> rups = pairRupsMap.get(pair);
				if (rups.size() > indInPairing) {
					rupsToPlot.add(rups.get(indInPairing));
					if (rupsToPlot.size() == maxRups)
						break;
				}
			}
			indInPairing++;
		}
		
		System.out.println("Plotting "+rupsToPlot+" ruptures");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		for (int rupIndex : rupsToPlot) {
			String rupPrefix = prefix+"_"+rupIndex;
			search.plotConnections(resourcesDir, rupPrefix, rupIndex, pairings, "Unique Connections");
			table.addColumn("![Rupture "+rupIndex+"]("
					+resourcesDir.getName()+"/"+rupPrefix+".png)");
		}
		table.finalizeLine();
		return table.wrap(maxCols, 0);
	}
	
	private static File plotFixedJumpDist(FaultSystemSolution u3Sol, List<ClusterRupture> u3ClusterRups,
			FaultSystemSolution rsSol, List<ClusterRupture> rsClusterRups, SectionDistanceAzimuthCalculator distAzCalc,
			String rsName, double minMag, float jumpDist, File outputDir) throws IOException {
		DiscretizedFunc u3Func = calcJumpDistFunc(u3Sol, u3ClusterRups, minMag, jumpDist);
		u3Func.scale(1d/u3Func.calcSumOfY_Vals());
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(u3Func);
		u3Func.setName("UCERF3");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		
		DiscretizedFunc rsFunc = calcJumpDistFunc(rsSol, rsClusterRups, minMag, jumpDist);
		rsFunc.scale(1d/rsFunc.calcSumOfY_Vals());
		rsFunc.setName(rsName);
		funcs.add(rsFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.DARK_GRAY));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "M≥"+(float)minMag+" Jump Comparison",
				"Num Jumps ≥"+(float)jumpDist+"km", "Fraction (Rate-Weighted)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(outputDir, "jumps_"+(float)jumpDist+"km").getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, false, null, new Range(0d, 1d));
		TickUnits tus = new TickUnits();
		TickUnit tu = new NumberTickUnit(1d);
		tus.add(tu);
		gp.getXAxis().setStandardTickUnits(tus);
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		gp.saveAsTXT(prefix+".txt");
		return new File(prefix+".png");
	}
	
	private static DiscretizedFunc calcJumpDistFunc(FaultSystemSolution sol, List<ClusterRupture> clusterRups,
			double minMag, float jumpDist) {
		EvenlyDiscretizedFunc solFunc = new EvenlyDiscretizedFunc(0d, 5, 1d);
		FaultSystemRupSet rupSet = sol.getRupSet();

		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double mag = rupSet.getMagForRup(r);

			if (mag < minMag)
				continue;
			
			ClusterRupture rup = clusterRups.get(r);
			int jumpsOverDist = 0;
			for (Jump jump : rup.getJumpsIterable()) {
				if ((float)jump.distance > jumpDist)
					jumpsOverDist++;
			}

			double rate = sol.getRateForRup(r);
			
			// indexes are fine to use here since it starts at zero with a delta of one 
			if (jumpsOverDist < solFunc.size())
				solFunc.set(jumpsOverDist, solFunc.getY(jumpsOverDist) + rate);
		}
		
		return solFunc;
	}
	
	private enum RakeType {
		RIGHT_LATERAL("Right-Lateral SS", "rl", Color.RED.darker()) {
			@Override
			public boolean isMatch(double rake) {
				return (float)rake >= -180f && (float)rake <= -170f
						|| (float)rake <= 180f && (float)rake >= 170f;
			}
		},
		LEFT_LATERAL("Left-Lateral SS", "ll", Color.GREEN.darker()) {
			@Override
			public boolean isMatch(double rake) {
				return (float)rake >= -10f && (float)rake <= 10f;
			}
		},
		REVERSE("Reverse", "rev", Color.BLUE.darker()) {
			@Override
			public boolean isMatch(double rake) {
				return (float)rake >= 80f && (float)rake <= 100f;
			}
		},
		NORMAL("Normal", "norm", Color.YELLOW.darker()) {
			@Override
			public boolean isMatch(double rake) {
				return (float)rake >= -100f && (float)rake <= -80f;
			}
		},
		OBLIQUE("Oblique", "oblique", Color.MAGENTA.darker()) {
			@Override
			public boolean isMatch(double rake) {
				for (RakeType type : values())
					if (type != this && type.isMatch(rake))
						return false;
				return true;
			}
		};
		
		private String name;
		private String prefix;
		private Color color;

		private RakeType(String name, String prefix, Color color) {
			this.name = name;
			this.prefix = prefix;
			this.color = color;
		}
		
		public abstract boolean isMatch(double rake);
	}
	
	private static Table<RakeType, RakeType, List<Double>> calcJumpAzimuths(
			List<ClusterRupture> rups, SectionDistanceAzimuthCalculator distAzCalc) {
		AzimuthCalc azCalc = new JumpAzimuthChangeFilter.SimpleAzimuthCalc(distAzCalc);
		Table<RakeType, RakeType, List<Double>> ret = HashBasedTable.create();
		for (RakeType r1 : RakeType.values())
			for (RakeType r2 : RakeType.values())
				ret.put(r1, r2, new ArrayList<>());
		for (ClusterRupture rup : rups) {
			for (Jump jump : rup.getJumpsIterable()) {
				RakeType sourceRake = null, destRake = null;
				for (RakeType type : RakeType.values()) {
					if (type.isMatch(jump.fromSection.getAveRake()))
						sourceRake = type;
					if (type.isMatch(jump.toSection.getAveRake()))
						destRake = type;
				}
				Preconditions.checkNotNull(sourceRake);
				Preconditions.checkNotNull(destRake);
				FaultSection before1 = rup.sectPredecessorsMap.get(jump.fromSection);
				if (before1 == null)
					continue;
				FaultSection before2 = jump.fromSection;
				double beforeAz = azCalc.calcAzimuth(before1, before2);
				FaultSection after1 = jump.toSection;
				for (FaultSection after2 : rup.sectDescendantsMap.get(after1)) {
					double afterAz = azCalc.calcAzimuth(after1, after2);
					double rawDiff = JumpAzimuthChangeFilter.getAzimuthDifference(beforeAz, afterAz);
					Preconditions.checkState(rawDiff >= -180 && rawDiff <= 180);
					double[] azDiffs;
					if ((float)before2.getAveDip() == 90f) {
						// strike slip, include both directions
						azDiffs = new double[] { rawDiff, -rawDiff };
					} else {
						// follow the aki & richards convention
						double dipDir = before2.getDipDirection();
						double dipDirDiff = JumpAzimuthChangeFilter.getAzimuthDifference(dipDir, beforeAz);
						if (dipDirDiff < 0)
							// this means that the fault dips to the right of beforeAz, we're good
							azDiffs = new double[] { rawDiff };
						else
							// this means that the fault dips to the left of beforeAz, flip it
							azDiffs = new double[] { -rawDiff };
					}
					for (double azDiff : azDiffs)
						ret.get(sourceRake, destRake).add(azDiff);
				}
			}
		}
		return ret;
	}
	
	private static Map<RakeType, List<Double>> getAzimuthsFrom (RakeType sourceRake,
			Table<RakeType, RakeType, List<Double>> azTable) {
		Map<RakeType, List<Double>> azMap;
		if (sourceRake == null) {
			azMap = new HashMap<>();
			for (RakeType type : RakeType.values())
				azMap.put(type, new ArrayList<>());
			for (RakeType source : RakeType.values()) {
				Map<RakeType, List<Double>> row = azTable.row(source);
				for (RakeType dest : row.keySet()) 
					azMap.get(dest).addAll(row.get(dest));
			}
		} else {
			azMap = azTable.row(sourceRake);
		}
		return azMap;
	}
	
	private static File plotJumpAzimuths(RakeType sourceRake, List<RakeType> destRakes,
			Table<RakeType, RakeType, List<Double>> azTable,
			File outputDir, String prefix, String title) throws IOException {
		Map<RakeType, List<Double>> azMap = getAzimuthsFrom(sourceRake, azTable);
		
		Range xRange = new Range(-180d, 180d);
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		
		List<Range> yRanges = new ArrayList<>();
		List<PlotSpec> specs = new ArrayList<>();
		
		for (int i=0; i<destRakes.size(); i++) {
			RakeType destRake = destRakes.get(i);
			
			HistogramFunction hist = HistogramFunction.getEncompassingHistogram(-179d, 179d, 15d);
			for (RakeType oRake : azMap.keySet()) {
				if (destRake != null && destRake != oRake)
					continue;
				for (double azDiff : azMap.get(oRake)) {
					hist.add(hist.getClosestXIndex(azDiff), 1d);
				}
			}

			Color color;
			String label;
			if (destRake == null) {
				color = Color.DARK_GRAY;
				label = "Any";
			} else {
				color = destRake.color;
				label = destRake.name;
			}
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(hist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, color));
			
			double maxY = Math.max(1.1*hist.getMaxY(), 1d);
			Range yRange = new Range(0d, maxY);
			
			PlotSpec spec = new PlotSpec(funcs, chars, title, "Azimuthal Difference", "Count");
			
			XYTextAnnotation ann = new XYTextAnnotation("To "+label, 175, maxY*0.975);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			spec.addPlotAnnotation(ann);
			
			specs.add(spec);
			yRanges.add(yRange);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(outputDir, prefix+".png");
		gp.getChartPanel().setSize(700, 1000);
		gp.saveAsPNG(file.getAbsolutePath());
		return file;
	}
	
	private static double azDiffDegreesToAngleRad(double azDiff) {
		// we want zero to be up, 90 to be right, 180 to be down, -90 to be left
		// sin/cos convention is zero at the right, 90 up, 180 left, -90 down
		
		Preconditions.checkState((float)azDiff >= (float)-180f && (float)azDiff <= 180f,
				"Bad azDiff: %s", azDiff);
		// first mirror it
		azDiff *= -1;
		// now rotate 90 degrees
		azDiff += 90d;
		
		return Math.toRadians(azDiff);
	}
	
	private static File plotJumpAzimuthsRadial(RakeType sourceRake, RakeType destRake,
			Table<RakeType, RakeType, List<Double>> azTable,
			File outputDir, String prefix, String title) throws IOException {
		System.out.println("Plotting "+title);
		Map<RakeType, List<Double>> azMap = getAzimuthsFrom(sourceRake, azTable);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Map<Float, List<Color>> azColorMap = new HashMap<>();
		
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(-179d, 179d, 15d);
		long totCount = 0;
		for (RakeType oRake : azMap.keySet()) {
			if (destRake != null && destRake != oRake)
				continue;
			for (double azDiff : azMap.get(oRake)) {
				hist.add(hist.getClosestXIndex(azDiff), 1d);
				
				Float azFloat = (float)azDiff;
				List<Color> colors = azColorMap.get(azFloat);
				if (colors == null) {
					colors = new ArrayList<>();
					azColorMap.put(azFloat, colors);
				}
				colors.add(oRake.color);
				totCount++;
			}
		}
		
		System.out.println("Have "+azColorMap.size()+" unique azimuths, "+totCount+" total");
//		Random r = new Random(azColorMap.keySet().size());
		double alphaEach = 0.025;
		if (totCount > 0)
			alphaEach = Math.max(alphaEach, 1d/totCount);
		for (Float azFloat : azColorMap.keySet()) {
			double sumRed = 0d;
			double sumGreen = 0d;
			double sumBlue = 0d;
			double sumAlpha = 0;
			int count = 0;
			for (Color azColor : azColorMap.get(azFloat)) {
				sumRed += azColor.getRed();
				sumGreen += azColor.getGreen();
				sumBlue += azColor.getBlue();
				if (sumAlpha < 1d)
					sumAlpha += alphaEach;
				count++;
			}
			double red = sumRed/(double)count;
			double green = sumGreen/(double)count;
			double blue = sumBlue/(double)count;
			if (red > 1d)
				red = 1d;
			if (green > 1d)
				green = 1d;
			if (blue > 1d)
				blue = 1d;
			if (sumAlpha > 1d)
				sumAlpha = 1d;
			Color color = new Color((float)red, (float)green, (float)blue, (float)sumAlpha);
//			if (destRake == null) {
//				// multipe types, choose a random color sampled from the actual colors
//				// for this azimuth
//				List<Color> colorList = azColorMap.get(azFloat);
//				color = colorList.get(r.nextInt(colorList.size()));
//			} else {
//				color = destRake.color;
//			}
			
			DefaultXY_DataSet line = new DefaultXY_DataSet();
			line.set(0d, 0d);
			double azRad = azDiffDegreesToAngleRad(azFloat);
			double x = Math.cos(azRad);
			double y = Math.sin(azRad);
			line.set(x, y);
			
			funcs.add(line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, color));
		}
		
		double dip;
		if (sourceRake == RakeType.LEFT_LATERAL || sourceRake == RakeType.RIGHT_LATERAL)
			dip = 90d;
		else if (sourceRake == RakeType.NORMAL || sourceRake == RakeType.REVERSE)
			dip = 60d;
		else
			dip = 75d;
		
		double traceLen = 0.5d;
		double lowerDepth = 0.25d;
		if (dip < 90d) {
			// add surface
			
			double horzWidth = lowerDepth/Math.tan(Math.toRadians(dip));
			DefaultXY_DataSet outline = new DefaultXY_DataSet();
			outline.set(0d, 0d);
			outline.set(horzWidth, 0d);
			outline.set(horzWidth, -traceLen);
			outline.set(0d, -traceLen);
			
			funcs.add(outline);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		}
		
		DefaultXY_DataSet trace = new DefaultXY_DataSet();
		trace.set(0d, 0d);
		trace.set(0d, -traceLen);
		
		funcs.add(trace);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Color.BLACK));
		PlotSpec spec = new PlotSpec(funcs, chars, title, "", " ");
		
		CPT cpt = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse();
		cpt = cpt.rescale(2d*Float.MIN_VALUE, 0.25d);
		cpt.setBelowMinColor(Color.WHITE);
		double halfDelta = 0.5*hist.getDelta();
		double innerMult = 0.95;
		double outerMult = 1.05;
		double sumY = Math.max(1d, hist.calcSumOfY_Vals());
		for (int i=0; i<hist.size(); i++) {
			double centerAz = hist.getX(i);
			double startAz = azDiffDegreesToAngleRad(centerAz-halfDelta);
			double endAz = azDiffDegreesToAngleRad(centerAz+halfDelta);
			
			List<Point2D> points = new ArrayList<>();
			
			double startX = Math.cos(startAz);
			double startY = Math.sin(startAz);
			double endX = Math.cos(endAz);
			double endY = Math.sin(endAz);
			
			points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
			points.add(new Point2D.Double(outerMult*startX, outerMult*startY));
			points.add(new Point2D.Double(outerMult*endX, outerMult*endY));
			points.add(new Point2D.Double(innerMult*endX, innerMult*endY));
			points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
			
			double[] polygon = new double[points.size()*2];
			int cnt = 0;
			for (Point2D pt : points) {
				polygon[cnt++] = pt.getX();
				polygon[cnt++] = pt.getY();
			}
			Color color = cpt.getColor((float)(hist.getY(i)/sumY));
			
			Stroke stroke = PlotLineType.SOLID.buildStroke(2f);
			spec.addPlotAnnotation(new XYPolygonAnnotation(polygon, stroke, Color.DARK_GRAY, color));
		}
		
		PaintScaleLegend cptBar = XYZGraphPanel.getLegendForCPT(cpt, "Fraction",
				24, 18, 0.05d, RectangleEdge.BOTTOM);
		spec.addSubtitle(cptBar);
		
		Range xRange = new Range(-1.1d, 1.1d);
		Range yRange = new Range(-1.1d, 1.1d);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(22);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		gp.getXAxis().setTickLabelsVisible(false);
		gp.getYAxis().setTickLabelsVisible(false);
		
		File file = new File(outputDir, prefix+".png");
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath());
		return file;
	}
	
	static LocationVector calcSlipVector(FaultSection sect, boolean footwall) {
		double strike = sect.getFaultTrace().getAveStrike();
		double rake = sect.getAveRake();
		double dip = sect.getAveDip();
		if (footwall)
			// we want the footwall motion
			rake = FaultUtils.getInRakeRange(rake + 180d);
		double azimuth = strike + rake;
		if ((float)rake == -180f || (float)rake == 0f || (float)rake == 180f) {
			// strike-slip motion, don't need to account for dip
			return new LocationVector(azimuth, 1d, 0d);
		}
		// at least some dip-slip motion
		
		 // this is the component of slip that is in the down-dip dir
		// it will be positive if rake ~ +90 (down-dip) and negative if rake ~ -90 (up-dip)
		double fractInDipDir = Math.sin(Math.toRadians(rake));
//		System.out.println("fractInDipDir: "+fractInDipDir);
		// sin(dip) = vertical/fractInDipDir
		double vertical = fractInDipDir*Math.sin(Math.toRadians(dip));
		// vert^2 + horiz^2 = 1^2
		double horizontal = Math.sqrt(1 - vertical*vertical);
		return new LocationVector(azimuth, horizontal, vertical);
	}
	
	static boolean isOnFootwall(RuptureSurface surf, Location loc) {
		return surf.getDistanceX(loc) < 0;
	}

}
