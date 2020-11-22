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
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
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
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRuptureBuilder;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.U3CoulombJunctionFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.ClusterCoulombCompatibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.ClusterPathCoulombCompatibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.CumulativeAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.GapWithinSectFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpAzimuthChangeFilter.AzimuthCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.JumpDistFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.MinSectsPerParentFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.MultiDirectionalPlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.ParentCoulombCompatibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.ParentCoulombCompatibilityFilter.Directionality;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.SingleClusterPerParentFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.SplayCountFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.TotalAzimuthChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.U3CompatibleCumulativeRakeChangeFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.DistCutoffClosestSectClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.UCERF3ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetDiagnosticsPageGen;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetDiagnosticsPageGen.RakeType;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetDiagnosticsPageGen.RupSetPlausibilityResult;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;
import org.opensha.sha.simulators.stiffness.RuptureCoulombResult;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.RuptureCoulombResult.RupCoulombQuantity;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessAggregationMethod;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessResult;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.FaultModels;
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
		double minMag = 6.5;
		int skipYears = 65000;
		double minFractForInclusion = 0.5;
		catalog.setFractForInclusion(minFractForInclusion);
		
		boolean rebuildSol = false;
		
		File catalogDir = catalog.getCatalogDir();
		
		File catalogOutputDir = new File(mainOutputDir, catalogDir.getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File multiFaultDir;
		if ((float)minFractForInclusion != (float)RSQSimCatalog.MIN_SUB_SECT_FRACT_DEFAULT)
			multiFaultDir = new File(catalogOutputDir, "multi_fault_area_fract_"+(float)minFractForInclusion);
		else
			multiFaultDir = new File(catalogOutputDir, "multi_fault");
		Preconditions.checkState(multiFaultDir.exists() || multiFaultDir.mkdir());
		
		File resourcesDir = new File(multiFaultDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		File fssDir = new File(catalogDir, "fss");
		Preconditions.checkState(fssDir.exists() || fssDir.mkdir());
		
		String catalogName = catalog.getName();
		String catalogType = "RSQSim";
		String catalogTypeFileName = "rsqsim";
		
		List<String> lines = new ArrayList<>();
		
		FaultBasedMapGen.MAP_LABEL_SIZE = 24;
		FaultBasedMapGen.MAP_LABEL_TICK_SIZE = 20;
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		FaultBasedMapGen.FAULT_THICKNESS = 4d;
		
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
		RuptureConnectionSearch rsConnSearch = new RuptureConnectionSearch(rupSet, distAzCalc,
				1000d, RuptureConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		
		System.out.println("Building ClusterRuptures for RSQSim");
//		List<ClusterRupture> rsClusterRups = RupSetDiagnosticsPageGen.buildClusterRups(rupSet, rsConnSearch);
		List<ClusterRupture> rsClusterRups = rupSet.getClusterRuptures();
		if (rsClusterRups == null) {
			rupSet.buildClusterRups(rsConnSearch);
			rsClusterRups = rupSet.getClusterRuptures();
		}
		Map<Jump, List<Integer>> rsJumpsToRupsMap = new HashMap<>();
		Map<Jump, Double> rsJumps = RupSetDiagnosticsPageGen.getJumps(sol, rsClusterRups, rsJumpsToRupsMap);
		HashSet<UniqueRupture> rsUniqueRups = new HashSet<>();
		for (ClusterRupture rup : rsClusterRups)
			rsUniqueRups.add(rup.unique);
		System.out.println("Found "+rsUniqueRups.size()+" unique ruptures");
		System.out.println("Detected "+rsJumps.size()+" RSQSim connections");
		
		System.out.println("Building ClusterRuptures for UCERF3");
		FaultSystemSolution u3Sol = catalog.getComparisonSolution();
		FaultSystemRupSet u3RupSet = u3Sol.getRupSet();
		RuptureConnectionSearch u3ConnSearch = new RuptureConnectionSearch(u3RupSet, distAzCalc,
				15d, RuptureConnectionSearch.CUMULATIVE_JUMPS_DEFAULT);
		List<ClusterRupture> u3ClusterRups = new ArrayList<>();
		for (int r=0; r<u3RupSet.getNumRuptures(); r++)
			u3ClusterRups.add(ClusterRupture.forOrderedSingleStrandRupture(
					u3RupSet.getFaultSectionDataForRupture(r), distAzCalc));
		Map<Jump, List<Integer>> u3JumpsToRupsMap = new HashMap<>();
		Map<Jump, Double> u3Jumps = RupSetDiagnosticsPageGen.getJumps(u3Sol, u3ClusterRups, u3JumpsToRupsMap);
		System.out.println("Detected "+u3Jumps.size()+" U3 connections");
		
		// we want actual catalog rupture counts before binning into U3 style ruptures
		// find the smallest rate, which will be 1/catLen, then numRups = solRate/minRate
		double minRate = StatUtils.min(sol.getRateForAllRups());
		
		List<PlausibilityFilter> filters = new ArrayList<>();
		double maxJumpDist = 5d;
		// these are sort of by-construction filters
		filters.add(new JumpDistFilter(maxJumpDist));
		CoulombRates coulombRates = CoulombRates.loadUCERF3CoulombRates(fm);
		ClusterConnectionStrategy u3ConnStrategy = new UCERF3ClusterConnectionStrategy(
				u3RupSet.getFaultSectionDataList(), distAzCalc, maxJumpDist, coulombRates);
		filters.add(new MinSectsPerParentFilter(2, true, true, u3ConnStrategy));
		filters.add(new GapWithinSectFilter());
		filters.add(new SplayCountFilter(0));
		
		// these are rupture properties themselves
		AzimuthCalc u3AzCalc = new JumpAzimuthChangeFilter.UCERF3LeftLateralFlipAzimuthCalc(distAzCalc);
		JumpAzimuthChangeFilter jumpAz = new JumpAzimuthChangeFilter(u3AzCalc, 60f);
		jumpAz.setErrOnCantEvaluate(true);
		filters.add(jumpAz);
		filters.add(new TotalAzimuthChangeFilter(u3AzCalc, 60f, true, true));
		filters.add(new CumulativeAzimuthChangeFilter(
				new JumpAzimuthChangeFilter.SimpleAzimuthCalc(distAzCalc), 560f));
//		filters.add(new CumulativeRakeChangeFilter(180f));
//		filters.add(new JumpCumulativeRakeChangeFilter(180f));
		filters.add(new U3CompatibleCumulativeRakeChangeFilter(180d));
		CoulombRatesTester coulombTester = new CoulombRatesTester(
				TestType.COULOMB_STRESS, 0.04, 0.04, 1.25d, true, true);
		filters.add(new U3CoulombJunctionFilter(coulombTester, coulombRates));
		
		RupSetPlausibilityResult result = RupSetDiagnosticsPageGen.testRupSetPlausibility(rsClusterRups, filters, null, rsConnSearch);
		TableBuilder table = RupSetDiagnosticsPageGen.getRupSetPlausibilityTable(result);
		
		// now plot
		String title = "Rupture Plausibility Filters, M≥"+(float)minMag+", SectArea≥"+(float)minFractForInclusion;
		String prefix = "filters_"+catParams;
		
		File filterPlot = RupSetDiagnosticsPageGen.plotRupSetPlausibility(result, resourcesDir, prefix, title);
		
		lines.add("## Plausibility Filter Comparisons");
		lines.add(topLink); lines.add("");
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/"+filterPlot.getName()+")");
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
			File plotFile = RupSetDiagnosticsPageGen.plotFixedJumpDist(rupSet, sol, rsClusterRups,
					catalogName, u3RupSet, u3Sol, u3ClusterRups, "UCERF3", distAzCalc, minMag,
					jumpDist, resourcesDir);
//			File plotFile = RupSetDiagnosticsPageGen.plotFixedJumpDist(u3Sol, u3ClusterRups, sol, rsClusterRups, distAzCalc, catalogName,
//					minMag, jumpDist, resourcesDir);
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
		
		Table<RakeType, RakeType, List<Double>> rsRakeAzTable =
				RupSetDiagnosticsPageGen.calcJumpAzimuths(rsClusterRups, distAzCalc);
		Table<RakeType, RakeType, List<Double>> u3RakeAzTable =
				RupSetDiagnosticsPageGen.calcJumpAzimuths(u3ClusterRups, distAzCalc);
		
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
			File plotFile = RupSetDiagnosticsPageGen.plotJumpAzimuths(sourceType, rakeTypes, rsRakeAzTable,
					resourcesDir, prefix, title);
			table.addColumn("!["+title+"](resources/"+plotFile.getName()+")");
			plotFile = RupSetDiagnosticsPageGen.plotJumpAzimuths(sourceType, rakeTypes, u3RakeAzTable,
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
				
				plotFile = RupSetDiagnosticsPageGen.plotJumpAzimuthsRadial(sourceType, destType, rsRakeAzTable,
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
		RupSetDiagnosticsPageGen.plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps.keySet(), rsColor, fullReg, 800);
		RupSetDiagnosticsPageGen.plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_ucerf3",
				"UCERF3 Connectivity", u3Jumps.keySet(), u3Color, fullReg, 800);
		RupSetDiagnosticsPageGen.plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_"+catalogTypeFileName+"_hires",
				catalogName+" Connectivity", rsJumps.keySet(), rsColor, fullReg, 3000);
		RupSetDiagnosticsPageGen.plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_ucerf3_hires",
				"UCERF3 Connectivity", u3Jumps.keySet(), u3Color, fullReg, 3000);
		Map<Jump, Double> rsUniqueJumps = new HashMap<>(rsJumps);
		Set<Jump> commonJumps = new HashSet<>();
		for (Jump jump : u3Jumps.keySet()) {
			if (rsUniqueJumps.containsKey(jump)) {
				rsUniqueJumps.remove(jump);
				commonJumps.add(jump);
			}
		}
		RupSetDiagnosticsPageGen.plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName,
				catalogName+" Unique Connectivity", rsUniqueJumps.keySet(), rsColor, fullReg, 800);
		RupSetDiagnosticsPageGen.plotConnectivityLines(sol.getRupSet(), resourcesDir, "sect_connectivity_unique_"+catalogTypeFileName+"_hires",
				catalogName+" Unique Connectivity", rsUniqueJumps.keySet(), rsColor, fullReg, 3000);
		Map<Jump, Double> u3UniqueJumps = new HashMap<>(u3Jumps);
		for (Jump jump : rsJumps.keySet())
			if (u3UniqueJumps.containsKey(jump))
				u3UniqueJumps.remove(jump);
		RupSetDiagnosticsPageGen.plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_unique_ucerf3",
				"UCERF3 Unique Connectivity", u3UniqueJumps.keySet(), u3Color, fullReg, 800);
		RupSetDiagnosticsPageGen.plotConnectivityLines(u3Sol.getRupSet(), resourcesDir, "sect_connectivity_unique_ucerf3_hires",
				"UCERF3 Unique Connectivity", u3UniqueJumps.keySet(), u3Color, fullReg, 3000);
		
		double maxConnDist = 0d;
		for (Jump jump : rsJumps.keySet())
			maxConnDist = Math.max(maxConnDist, jump.distance);
		for (Jump jump : u3Jumps.keySet())
			maxConnDist = Math.max(maxConnDist, jump.distance);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps, rsUniqueJumps, maxConnDist,
				rsColor, false, false);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_ucerf3",
				"UCERF3 Connectivity", u3Jumps, u3UniqueJumps, maxConnDist,
				u3Color, false, false);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName,
				catalogName+" Connectivity", rsJumps, rsUniqueJumps, maxConnDist,
				rsColor, true, false);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_ucerf3",
				"UCERF3 Connectivity", u3Jumps, u3UniqueJumps, maxConnDist,
				u3Color, true, false);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName+"_log",
				catalogName+" Connectivity", rsJumps, rsUniqueJumps, maxConnDist,
				rsColor, true, true);
		RupSetDiagnosticsPageGen.plotConnectivityHistogram(resourcesDir, "sect_connectivity_hist_rates_ucerf3_log",
				"UCERF3 Connectivity", u3Jumps, u3UniqueJumps, maxConnDist,
				u3Color, true, true);
		
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
		RupSetDiagnosticsPageGen.plotConnectivityLines(rupSet, resourcesDir, combConnPrefix, "Combined Connectivity",
				connectionsList, connectedColors, connNames, fullReg, 800);
		RupSetDiagnosticsPageGen.plotConnectivityLines(rupSet, resourcesDir, combConnPrefix+"_hires", "Combined Connectivity",
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
		rsPlot = new File(resourcesDir, "sect_connectivity_hist_rates_"+catalogTypeFileName+"_log.png");
		Preconditions.checkState(rsPlot.exists());
		u3Plot = new File(resourcesDir, "sect_connectivity_hist_rates_ucerf3_log.png");
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
		table = RupSetDiagnosticsPageGen.plotConnRupExamples(rsConnSearch, rupSet, rsUniqueJumps.keySet(),
				rsJumpsToRupsMap, maxRups, maxCols, resourcesDir, "rs_conn_example");
		lines.add("");
		if (table == null)
			lines.add("*N/A*");
		else
			lines.addAll(table.build());
		lines.add("");
		lines.add("**UCERF3 Ruptures with Unique Connections**");
		table = RupSetDiagnosticsPageGen.plotConnRupExamples(u3ConnSearch, u3RupSet, u3UniqueJumps.keySet(),
				u3JumpsToRupsMap, maxRups, maxCols, resourcesDir, "u3_conn_example");
		lines.add("");
		if (table == null)
			lines.add("*N/A*");
		else
			lines.addAll(table.build());
		lines.add("");
		
		// now stiffness comparisons
		Map<String, String> params = catalog.getParams();
		double lambda = 30000;
		double mu = 30000;
		if (params.containsKey("lameLambda"))
			lambda = Double.parseDouble(params.get("lameLambda"));
		if (params.containsKey("lameMu"))
			lambda = Double.parseDouble(params.get("lameMu"));
		double coeffOfFriction = 0.5;
		SubSectStiffnessCalculator stiffnessCalc = new SubSectStiffnessCalculator(
				rupSet.getFaultSectionDataList(), 2d, lambda, mu, coeffOfFriction);
		StiffnessType type = StiffnessType.CFF;
		
		System.out.println("Plotting Jump Coulomb");
		
		lines.add("## Coulomb Comparisons");
		lines.add(topLink); lines.add("");
		lines.add("These calculations compute "+type.getHTML()+" between each subsection. This is done by "
				+ "discretizing each subsection into "+(float)stiffnessCalc.getGridSpacing()+" x "
				+ (float)stiffnessCalc.getGridSpacing()+" km rectangular patches, and them computing "
				+ "stiffness between each patch on each surface. So for two subsections with N and M "
				+ "patches each, we get N x M total stiffness calculations. We then use various aggregation "
				+ "methods to to compute compatibility between those subsections:");
		lines.add("");
		StiffnessAggregationMethod[] plotQuantities = { StiffnessAggregationMethod.MEDIAN,
				StiffnessAggregationMethod.MAX, StiffnessAggregationMethod.FRACT_POSITIVE };
		for (StiffnessAggregationMethod q : plotQuantities)
			lines.add("* **"+q.toString()+"**: "+q.getDescription());
		lines.add("");
		
		List<PlausibilityFilter> coulombFilters = new ArrayList<>();
		coulombFilters.add(new ParentCoulombCompatibilityFilter(
				stiffnessCalc, StiffnessAggregationMethod.MEDIAN, 0f, Directionality.EITHER));
		coulombFilters.add(new ClusterCoulombCompatibilityFilter(
				stiffnessCalc, StiffnessAggregationMethod.MEDIAN, 0f));
		coulombFilters.add(new ClusterPathCoulombCompatibilityFilter(
				stiffnessCalc, StiffnessAggregationMethod.MEDIAN, 0f));
		
		System.out.println("Testing new coulomb filters");
		result = RupSetDiagnosticsPageGen.testRupSetPlausibility(rsClusterRups, coulombFilters, null, rsConnSearch);
		table = RupSetDiagnosticsPageGen.getRupSetPlausibilityTable(result);
		
		// now plot
		title = "New Coulomb Plausibility Filters, M≥"+(float)minMag+", SectArea≥"+(float)minFractForInclusion;
		prefix = "filters_coulomb_"+catParams;
		
		filterPlot = RupSetDiagnosticsPageGen.plotRupSetPlausibility(result, resourcesDir, prefix, title);
		
		lines.add("### Coulomb Plausibility Filter Comparisons");
		lines.add(topLink); lines.add("");
		lines.add("![Plausibility Filter]("+resourcesDir.getName()+"/"+filterPlot.getName()+")");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		for (boolean fullClusters : new boolean[] {false, true}) {
			table = MarkdownUtils.tableBuilder();
			table.addLine("RSQSim", "UCERF3");
			
			prefix = "jump_cff";
			if (fullClusters)
				prefix += "_clusters";
			else
				prefix += "_subsections";

			System.out.println("Plotting RSQSim, fullClusters="+fullClusters);
			File[][] rsPlots = plotJumpStiffnessHistogram(rsClusterRups, sol, stiffnessCalc, type,
					plotQuantities, fullClusters, resourcesDir, prefix+"_rsqsim");
			System.out.println("Plotting UCERF3, fullClusters="+fullClusters);
			File[][] u3Plots = plotJumpStiffnessHistogram(u3ClusterRups, u3Sol, stiffnessCalc, type,
					plotQuantities, fullClusters, resourcesDir, prefix+"_ucerf3");
			for (int q=0; q<plotQuantities.length; q++) {
				table.initNewLine();
				if (rsPlots[q] == null)
					table.addColumn("N/A");
				else
					table.addColumn("![plot](resources/"+rsPlots[q][q].getName()+")");
				if (u3Plots[q] == null)
					table.addColumn("N/A");
				else
					table.addColumn("![plot](resources/"+u3Plots[q][q].getName()+")");
				table.finalizeLine();
			}
			// now scatters
			for (int q1=0; q1<plotQuantities.length; q1++) {
				for (int q2=q1+1; q2<plotQuantities.length; q2++) {
					table.initNewLine();
					if (rsPlots[q1] == null)
						table.addColumn("N/A");
					else
						table.addColumn("![plot](resources/"+rsPlots[q1][q2].getName()+")");
					if (u3Plots[q1] == null)
						table.addColumn("N/A");
					else
						table.addColumn("![plot](resources/"+u3Plots[q1][q2].getName()+")");
					table.finalizeLine();
				}
			}
			if (fullClusters) {
				lines.add("### Jump Coulomb Comparisons, Full Clusters");
				lines.add(topLink); lines.add("");
				lines.add("This gives the aggregated coulomb compatibility of every jump, computed between "
						+ "the full subsection clusters on either side of that jump.");
			} else {
				lines.add("### Jump Coulomb Comparisons, Jump Subsections Only");
				lines.add(topLink); lines.add("");
				lines.add("This gives the aggregated coulomb compatibility of every jump, computed between "
						+ "only the closest subsections on either side of that jump.");
			}
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
		}
		
		table = MarkdownUtils.tableBuilder();
		table.addLine("RSQSim", "UCERF3");
		
		prefix = "rup_to_cluster_cff";

		System.out.println("Plotting Rup to Cluster, RSQSim");
		File[][] rsPlots = plotRupToClusterCoulombHistogram(rsClusterRups, sol, stiffnessCalc,
				type, plotQuantities, resourcesDir, prefix+"_rsqsim");
		System.out.println("Plotting Rup to Cluster, UCERF3");
		File[][] u3Plots = plotRupToClusterCoulombHistogram(u3ClusterRups, u3Sol, stiffnessCalc,
				type, plotQuantities, resourcesDir, prefix+"_ucerf3");
		for (int q=0; q<plotQuantities.length; q++) {
			table.initNewLine();
			if (rsPlots[q] == null)
				table.addColumn("N/A");
			else
				table.addColumn("![plot](resources/"+rsPlots[q][q].getName()+")");
			if (u3Plots[q] == null)
				table.addColumn("N/A");
			else
				table.addColumn("![plot](resources/"+u3Plots[q][q].getName()+")");
			table.finalizeLine();
		}
		// now scatters
		for (int q1=0; q1<plotQuantities.length; q1++) {
			for (int q2=q1+1; q2<plotQuantities.length; q2++) {
				table.initNewLine();
				if (rsPlots[q1] == null)
					table.addColumn("N/A");
				else
					table.addColumn("![plot](resources/"+rsPlots[q1][q2].getName()+")");
				if (u3Plots[q1] == null)
					table.addColumn("N/A");
				else
					table.addColumn("![plot](resources/"+u3Plots[q1][q2].getName()+")");
				table.finalizeLine();
			}
		}
		lines.add("### Jump Coulomb Comparisons, Full Rupture to Cluster");
		lines.add(topLink); lines.add("");
		lines.add("This gives the aggregated coulomb compatibility of each individual subsection cluster, "
				+ "computed to the rest of the rupture");
		lines.add("");
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("### Rupture Aggregate Section Coulomb Comparisons");
		lines.add(topLink); lines.add("");
		lines.add("This gives aggregated coulomb compatibilty of entire ruptures. We first calculate a net "
				+ "coulomb stress change for each subsection, conditioned on unit slip on all other subsection "
				+ "which participate in the rupture. We then plot the mean/min of those subsection net coulomb "
				+ "stress changes to ge a single number for the full rupture. We also keep track of the number "
				+ "of subsections in each rupture which are net negative.");
		lines.add("");

		System.out.println("Calculating RSQSim Rupture Coulomb");
		List<RuptureCoulombResult> rsRupCoulombs = calcRupCoulombs(
				rsClusterRups, stiffnessCalc, StiffnessAggregationMethod.MEDIAN);
		System.out.println("Calculating UCERF3 Rupture Coulomb");
		List<RuptureCoulombResult> u3RupCoulombs = calcRupCoulombs(
				u3ClusterRups, stiffnessCalc, StiffnessAggregationMethod.MEDIAN);
		
		System.out.println("Plotting Rupture Coulomb");
		
		table = MarkdownUtils.tableBuilder();
		table.addLine("RSQSim", "UCERF3");
		StiffnessAggregationMethod stiffAggMethod = StiffnessAggregationMethod.MEDIAN;
		for (RupCoulombQuantity quantity : RupCoulombQuantity.values()) {
			System.out.println("Plotting RSQSim, "+quantity);
			prefix = "rup_cff_"+quantity.name();
			File plot = plotRupCoulombHistogram(rsRupCoulombs, sol, quantity,
					resourcesDir, prefix+"_rsqsim");
			table.initNewLine();
			if (plot == null)
				table.addColumn("N/A");
			else
				table.addColumn("![plot](resources/"+plot.getName()+")");
			System.out.println("Plotting UCERF3, "+quantity);
			plot = plotRupCoulombHistogram(u3RupCoulombs, u3Sol, quantity,
					resourcesDir, prefix+"_ucerf3");
			if (plot == null)
				table.addColumn("N/A");
			else
				table.addColumn("![plot](resources/"+plot.getName()+")");
			table.finalizeLine();
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
	
	private static File[][] plotJumpStiffnessHistogram(List<ClusterRupture> ruptures, FaultSystemSolution sol,
			SubSectStiffnessCalculator calc, StiffnessType type, StiffnessAggregationMethod[] quantities, boolean fullClusters,
			File resourcesDir, String prefix) throws IOException {
		List<List<Double>> valuesList = new ArrayList<>();
		double[] allMaxAbs = new double[quantities.length];
		for (int q=0; q<quantities.length; q++)
			valuesList.add(new ArrayList<>());
		List<Double> rates = new ArrayList<>();
		
		for (int r=0; r<ruptures.size(); r++) {
			ClusterRupture rup = ruptures.get(r);
			double rate = sol.getRateForRup(r);
			for (Jump rawJump : rup.getJumpsIterable()) {
				for (Jump jump : new Jump[] { rawJump, rawJump.reverse() }) {
//					System.out.println("Calculating for jump "+jump);
					StiffnessResult stiffness;
					if (fullClusters)
						stiffness = calc.calcClusterStiffness(type, jump.fromCluster, jump.toCluster);
					else
						stiffness = calc.calcStiffness(type, jump.fromSection, jump.toSection);
					if (stiffness == null)
						continue;
					for (int q=0; q<quantities.length; q++) {
						double val = stiffness.getValue(quantities[q]);
						allMaxAbs[q] = Math.max(allMaxAbs[q], Math.abs(val));
						valuesList.get(q).add(val);
					}
					rates.add(0.5*rate); // 0.5 since we're including it both forwards and backwards
				}
			}
		}
		if (valuesList.get(0).isEmpty())
			return null;
		File[][] ret = new File[quantities.length][quantities.length];
		for (int q=0; q<quantities.length; q++) {
			double maxAbs = allMaxAbs[q];
			String title = "Jump "+quantities[q].toString();
			if (fullClusters)
				title += ", Full Clusters";
			else
				title += ", Closest Subsections";
			
			String xAxisLabel = type+", "+quantities[q];
			List<Double> values = valuesList.get(q);
			
			ret[q][q] = plotValueHist(resourcesDir, prefix+"_"+quantities[q].name(), values, rates,
					maxAbs, title, xAxisLabel, quantities[q]);
			for (int q2=q+1; q2<quantities.length; q2++) {
				String scatterTitle = "Jump "+type.toString()+" Scatter";
				String scatterPrefix = prefix+"_"+quantities[q].name()+"_vs_"+quantities[q2].name();
				ret[q][q2] = plotValueScatter(resourcesDir, scatterPrefix, values, valuesList.get(q2),
						maxAbs, allMaxAbs[q2], scatterTitle, quantities[q], quantities[q2]);
				ret[q2][q] = ret[q][q2];
			}
		}
		
		return ret;
	}

	static File plotValueHist(File resourcesDir, String myPrefix, List<Double> values, List<Double> rates, double maxAbs,
			String title, String xAxisLabel, StiffnessAggregationMethod quantity) throws IOException {
		maxAbs = Math.ceil(10d*maxAbs)/10d;
		
		HistogramFunction hist;
		Range xRange;
		if (quantity == StiffnessAggregationMethod.FRACT_POSITIVE) {
			hist = HistogramFunction.getEncompassingHistogram(0d, 1d, 0.02);
			xRange = new Range(0d, hist.getMaxX()+0.5*hist.getDelta());
		} else {
			hist = new HistogramFunction(-maxAbs, maxAbs, 50);
			xRange = new Range(-maxAbs-0.5*hist.getDelta(), maxAbs+0.5*hist.getDelta());
		}
		for (int i=0; i<values.size(); i++) {
			double val = values.get(i);
			hist.add(hist.getClosestXIndex(val), rates.get(i));
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, "Jump Annual Rate");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		double minNonZero = Double.POSITIVE_INFINITY;
		for (Point2D pt : hist)
			if (pt.getY() > 0)
				minNonZero = Math.min(minNonZero, pt.getY());
		
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
				Math.pow(10, Math.ceil(Math.log10(hist.getMaxY()))));
		
		File pngFile = new File(resourcesDir, myPrefix+".png");
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}

	static File plotValueScatter(File resourcesDir, String myPrefix, List<Double> values1, List<Double> values2,
			double maxAbs1, double maxAbs2, String title, StiffnessAggregationMethod quantity1,
			StiffnessAggregationMethod quantity2) throws IOException {
		maxAbs1 = Math.ceil(10d*maxAbs1)/10d;
		maxAbs2 = Math.ceil(10d*maxAbs2)/10d;
		
		Range xRange, yRange;
		if (quantity1 == StiffnessAggregationMethod.FRACT_POSITIVE)
			xRange = new Range(0d, 1d);
		else
			xRange = new Range(-maxAbs1, maxAbs1);
		if (quantity2 == StiffnessAggregationMethod.FRACT_POSITIVE)
			yRange = new Range(0d, 1d);
		else
			yRange = new Range(-maxAbs2, maxAbs2);
		Preconditions.checkState(values1.size() == values2.size());
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		for (int i=0; i<values1.size(); i++)
			scatter.set(values1.get(i), values2.get(i));
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, quantity1.toString(), quantity2.toString());
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		File pngFile = new File(resourcesDir, myPrefix+".png");
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}
	
	private static List<RuptureCoulombResult> calcRupCoulombs(List<ClusterRupture> ruptures,
			SubSectStiffnessCalculator calc, StiffnessAggregationMethod stiffAggMethod) {
		List<Future<RuptureCoulombResult>> futures = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		for (ClusterRupture rup : ruptures)
			futures.add(exec.submit(new RupCoulombCalcCall(rup, calc, stiffAggMethod)));
		List<RuptureCoulombResult> results = new ArrayList<>();
		for (Future<RuptureCoulombResult> future : futures) {
			try {
				results.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		exec.shutdown();
		return results;
	}
	
	private static class RupCoulombCalcCall implements Callable<RuptureCoulombResult> {
		
		private ClusterRupture rup;
		private SubSectStiffnessCalculator calc;
		private StiffnessAggregationMethod stiffAggMethod;

		public RupCoulombCalcCall(ClusterRupture rup, SubSectStiffnessCalculator calc,
				StiffnessAggregationMethod stiffAggMethod) {
			this.rup = rup;
			this.calc = calc;
			this.stiffAggMethod = stiffAggMethod;
		}

		@Override
		public RuptureCoulombResult call() throws Exception {
			if (rup.unique.size() < 2)
				return null;
			else
				return new RuptureCoulombResult(rup, calc, stiffAggMethod);
		}
		
	}
	
	private static File plotRupCoulombHistogram(List<RuptureCoulombResult> rupCoulombs, FaultSystemSolution sol,
			RupCoulombQuantity quantity, File resourcesDir, String prefix) throws IOException {
		List<Double> values = new ArrayList<>();
		List<Double> rates = new ArrayList<>();
		double maxAbs = 0d;
		for (int r=0; r<rupCoulombs.size(); r++) {
			RuptureCoulombResult result = rupCoulombs.get(r);
			if (result == null)
				continue;
			double rate = sol.getRateForRup(r);
			double val = result.getValue(quantity);
			maxAbs = Math.max(Math.abs(val), maxAbs);
			values.add(val);
			rates.add(rate);
		}
		if (values.isEmpty())
			return null;
		maxAbs = Math.ceil(10d*maxAbs)/10d;
		
		HistogramFunction hist;
		Range xRange;
		if (quantity == RupCoulombQuantity.MEAN_SECT_FRACT_POSITIVES
				|| quantity == RupCoulombQuantity.MIN_SECT_FRACT_POSITIVES) {
			hist = HistogramFunction.getEncompassingHistogram(0d, 1d, 0.02);
			xRange = new Range(0d, hist.getMaxX()+0.5*hist.getDelta());
		} else if (quantity == RupCoulombQuantity.NUM_NET_NEGATIVE_SECTS) {
			hist = new HistogramFunction(0d, 10d, 11);
			xRange = new Range(-0.5, hist.getMaxX()+0.5*hist.getDelta());
		} else {
			hist = new HistogramFunction(-maxAbs, maxAbs, 50);
			xRange = new Range(-maxAbs-0.5*hist.getDelta(), maxAbs+0.5*hist.getDelta());
		}
		for (int i=0; i<values.size(); i++) {
			double val = values.get(i);
			hist.add(hist.getClosestXIndex(val), rates.get(i));
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		String title = "Full Rupture Aggregated Coulomb";
		
		String xAxisLabel = quantity.toString();
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, "Annual Rate");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		double minNonZero = Double.POSITIVE_INFINITY;
		for (Point2D pt : hist)
			if (pt.getY() > 0)
				minNonZero = Math.min(minNonZero, pt.getY());
		
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
				Math.pow(10, Math.ceil(Math.log10(hist.getMaxY()))));
		
		File pngFile = new File(resourcesDir, prefix+".png");
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}
	
	private static File[][] plotRupToClusterCoulombHistogram(List<ClusterRupture> ruptures, FaultSystemSolution sol,
			SubSectStiffnessCalculator calc, StiffnessType type, StiffnessAggregationMethod[] quantities,
			File resourcesDir, String prefix) throws IOException {
		List<List<Double>> valuesList = new ArrayList<>();
		double[] allMaxAbs = new double[quantities.length];
		for (int q=0; q<quantities.length; q++)
			valuesList.add(new ArrayList<>());
		List<Double> rates = new ArrayList<>();
		for (int r=0; r<ruptures.size(); r++) {
			ClusterRupture rupture = ruptures.get(r);
			if (rupture.unique.size() < 2 || rupture.getTotalNumJumps() == 0)
				continue;
			double rate = sol.getRateForRup(r);
			for (FaultSubsectionCluster cluster : rupture.getClustersIterable()) {
				StiffnessResult stiffness = calc.calcAggRupToClusterStiffness(type, rupture, cluster);
				if (stiffness == null)
					continue;
				for (int q=0; q<quantities.length; q++) {
					double val = stiffness.getValue(quantities[q]);
					allMaxAbs[q] = Math.max(allMaxAbs[q], Math.abs(val));
					valuesList.get(q).add(val);
				}
				rates.add(rate);
			}
		}
		if (valuesList.get(0).isEmpty())
			return null;
		File[][] ret = new File[quantities.length][quantities.length];
		for (int q=0; q<quantities.length; q++) {
			double maxAbs = allMaxAbs[q];
			String title = "Aggregate Full Rup to Cluster "+quantities[q].toString();
			
			String xAxisLabel = type+", "+quantities[q];
			List<Double> values = valuesList.get(q);
			
			ret[q][q] = plotValueHist(resourcesDir, prefix+"_"+quantities[q].name(), values, rates,
					maxAbs, title, xAxisLabel, quantities[q]);
			for (int q2=q+1; q2<quantities.length; q2++) {
				String scatterTitle = "Aggregate Full Rup to Cluster "+type.toString()+" Scatter";
				String scatterPrefix = prefix+"_"+quantities[q].name()+"_vs_"+quantities[q2].name();
				ret[q][q2] = plotValueScatter(resourcesDir, scatterPrefix, values, valuesList.get(q2),
						maxAbs, allMaxAbs[q2], scatterTitle, quantities[q], quantities[q2]);
				ret[q2][q] = ret[q][q2];
			}
		}
		
		return ret;
	}

}
