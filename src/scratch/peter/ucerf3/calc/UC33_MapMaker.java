package scratch.peter.ucerf3.calc;

import static com.google.common.base.Charsets.US_ASCII;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.CPT_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.GRID_SPACING_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.LOG_PLOT_NAME;
import static org.opensha.nshmp2.tmp.TestGrid.CA_NSHMP;
import static org.opensha.nshmp2.tmp.TestGrid.CA_RELM;
import static org.opensha.nshmp2.tmp.TestGrid.LOS_ANGELES;
import static org.opensha.nshmp2.tmp.TestGrid.SAN_FRANCISCO;
import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.nshmp2.util.Period.GM2P00;
import static org.opensha.nshmp2.util.Period.GM3P00;
import static scratch.peter.curves.ProbOfExceed.PE10IN50;
import static scratch.peter.curves.ProbOfExceed.PE2IN50;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.nshmp2.calc.HazardResultWriterLocal;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.peter.curves.ProbOfExceed;
import scratch.peter.nshmp.BinaryCurves;
import scratch.peter.nshmp.CurveContainer;
import scratch.peter.nshmp.NSHMP_DataUtils;
import scratch.peter.nshmp.NSHMP_PlotUtils;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

/**
 * Add comments here
 * 
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC33_MapMaker {

	private static final Splitter SPLIT = Splitter.on('_');
	private static final Splitter SPLIT_COMMA = Splitter.on(',');
	private static final Splitter SPLIT_TAB = Splitter.on('\t');
	private static final Joiner JOIN = Joiner.on(',');
	private static String format = "%.2f";
	private static final String S = File.separator;
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static final String ROOT = "tmp/UC33/maps/";
	private static final String SRC = ROOT + "src/";

	public static void main(String[] args) throws IOException {
//		makeBrAvgRatioMapUC33();
//		makePrelimBrAvgHazardMaps();
//		doInversionRunAnalysis("UC33-10runs-3sec", GM3P00, PE2IN50);
//		buildLogicTreeVarMaps();
//		buildBrAvgMaps();
//		buildBrAvgHazardMaps();
		
		// these two methods can be reconfigured to compute different map pairs
		// currently used for GMPE and newGridSrc comparisons
//		buildNSHMPratioDiff(); 
//		buildUC3ratioDiff();
		
//		buildBgFltMapsTestCPT();
		
//		buildUC3gmpe13();
//		buildBrAvgHazardMaps();
		
//		makeLocalRatioMaps();
//		makeLocalHazardMaps();
		
//		createCAxyzFiles();
		
//		buildGMPE14changes();
//		buildGMPE14_AS_bugfix();
		
//		buildGMEP13comparisons();
//		buildGMEP13comparisonsBinaries();
//		buildIdrissComparison();

//		buildUC3_NSHMP_binaries();
//		buildUC3_NSHMP_binaries40();
//		buildUC3_NSHMP_determ();
		
		// for Morgan:
//		buildUC3_NSHMP_BG_binariesMorgan();
//		buildUC3_NSHMP_binariesMorgan();
//		buildUC3_NSHMP_binariesMorgan2();
//		buildUC2_NSHMP_binaries();
//		build_UC3_GMM08_binary();
		
//		finalBSSCcheck();
//		finalMapsDebug();
//		finalMapsSpacingComparison();
//		tmp();
		
//		consolidateFinalDeterm();
		
		// time dependence
//		buildUCtimeDependentMaps();
						
//		testNSHMP();
//		buildUC3_NSHMP_versionMaps();
//		build_NSHMP_comp1();
//		build_NSHMP_comp2();
//		build_NSHMP_comp3();
//		build_NSHMP_comp4();
//		build_NSHMP_BgFltMaps();
//		build_NSHMP_LocalRatioMaps();
//		build_NSHMP_TotalHazardMaps();
		
//		build_NSHMP13B_TotalHazardMap();
		
		build_NSHMP_LocalRatioMaps_forEQS();
		
//		resamplingTest();
		
		// for Nico:
		//  -- 2013 model with 2008 NGAs
		//  -- 2013 model with 2013 NGAs
//		combineCurves();
		
		//  -- rtgm & risk comparisons with derivitave data sets of above curves
//		riskRatioMapsCA();
//		riskRatioMapsCA2();
//		riskRatioMapsUS();
//		riskRatioMapsCA2_CrBeta0p8();
//		riskRatioMapsUS_CrBeta0p8();
		
	}

	private static void testNSHMP() {
		// build brAvg curve containers
		Period p = GM0P20;
		ProbOfExceed pe = PE2IN50;
		TestGrid gr = TestGrid.CA_RELM;
		double spacing = 0.1;
		String srcDir = SRC + "NSHMP13-det-flt-8/UC33brAvg_FM31_ZENGBB";
//		CurveContainer cc_FM = buildBrAvgCurveContainer(srcDir, p, gr, spacing);
//		GeoDataSet xyz_FM = NSHMP_DataUtils.extractPE(cc_FM, gr.grid(spacing), pe);
		GeoDataSet xyz = loadSingle(srcDir, pe, gr, p, spacing);
		makeHazardMap(xyz, spacing, p, pe, gr, ROOT + "test/NSHMP-" + p.getLabel());
	}
	
	/*
	 * Used to make ratio maps for UCERF3.2 report 
	 */
	private static void buildLogicTreeVarMaps() throws IOException {
		List<String> brOverList = null;
		String brUnder = null;
		String branches = null;
		String srcDir = SRC + "UC33/";
		String outDir = ROOT + "LogicTreeRatios/";
		
//		// Mmax and Mgt5 comparisons
//		brOverList = Lists.newArrayList(
//			"M565-MX73", "M579-MX73", "M596-MX73",
//			"M565-MX76", "M579-MX76", "M596-MX76",
//			"M565-MX79", "M579-MX79", "M596-MX79");
//		brUnder = "M579-MX76";
//		makeRatioMaps(srcDir, outDir, brOverList, brUnder);
		
//		// other branch node comparisons
//		brOverList = Lists.newArrayList(
//			"FM31", "FM32", 
//			"DM_ZBB", "DM_ABM", "DM_GEOL", "DM_NEOK",
//			"MS_SH09M", "MS_ELLB", "MS_ELLBSL", "MS_HB08", "MS_SHCSD",
//			"DSR_TAP", "DSR_UNI",
//			"M565", "M579", "M596",
//			"MX73", "MX76", "MX79",
//			"UC2", "UC3");
//		brUnder = "all";
//		makeRatioMaps(srcDir, outDir, brOverList, brUnder);
		
//		// ratio map UC33 over NSMP (1440 branches)
//		branches = "all";
//		makeNSHMPratioMap(srcDir, outDir, branches);
//		branches = "UC3";
//		makeNSHMPratioMap(srcDir, outDir, branches);
//		branches = "UC2";
//		makeNSHMPratioMap(srcDir, outDir, branches);
	}
	
	private static void build_NSHMP_LocalRatioMaps_forEQS() {
		
		double spacing = 0.02;

		List<TestGrid> localGrids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		for (Period p : periods) {
						
			for (ProbOfExceed pe : PEs) {
				for (TestGrid locGrid : localGrids) {
					
					// NSHMP14 hi res local
					String overPath = SRC + "NSHMP14-SF-LA" + S + "mean_ucerf3_sol" + S + locGrid.name() + S + p + S + "curves.csv";
					File overFile = new File(overPath);
					CurveContainer ccOver = CurveContainer.create(overFile, locGrid, spacing);
					
					// NSHMP08 hi res local
					String underPath = SRC + "NSHMP08-SF-LA" + S + "ALL" + S + locGrid.name() + S + p + S + "curves.csv";
					File underFile = new File(underPath);
					CurveContainer ccUnder = CurveContainer.create(underFile, locGrid, spacing);
				
					// data
					GriddedRegion gr = locGrid.grid(spacing);
					GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(ccOver, gr, pe);
					GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(ccUnder, gr, pe);
					
					GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
					GeoDataSet xyzDiff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
					String id = locGrid + "-" + p.getLabel() +  "-" + pe;
					String ratioFile = ROOT + "NSHMP14-SF-LA" + S + "FM2" + S + id + "-ratio";
					String diffFile = ROOT + "NSHMP14-SF-LA" + S + "FM2" + S + id + "-diff";
					makeRatioPlotNSHMP(xyzRatio, spacing, locGrid.bounds(), ratioFile, "GM ratio", true);
					makeDiffPlotNSHMP(xyzDiff, spacing, locGrid.bounds(), diffFile, "GM diff", true);
				}
			}
		}
	}
	
	private static void build_NSHMP_LocalRatioMaps() {
		
		double spacing = 0.05;

		List<TestGrid> localGrids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50);
		
		for (Period p : periods) {
			// NSHMP13 hi res 40 branchAvg curves CA_NSHMP
			File fltDir = new File(SRC + "NSHMP13-40-FLT");
			CurveContainer ccOver = buildBrAvgCurveContainer(fltDir, p, CA_NSHMP, spacing);
			File bgDir = new File(SRC + "NSHMP13-2-BG");
			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, CA_NSHMP, spacing);
			ccOver.add(bgCC);
			
			for (ProbOfExceed pe : PEs) {
				for (TestGrid locGrid : localGrids) {
					// NSHMP08 hi res local
					String underPath = SRC + "nshmp_SF-LA" + S + locGrid.name() + S + p + S + "curves.csv";
					File underFile = new File(underPath);
					CurveContainer ccUnder = CurveContainer.create(underFile, locGrid, spacing);
				
					// data
					GriddedRegion gr = locGrid.grid(spacing);
					GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(ccOver, gr, pe);
					GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(ccUnder, gr, pe);
					
					GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
					GeoDataSet xyzDiff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
					String id = locGrid + "-" + p.getLabel() +  "-" + pe;
					String ratioFile = ROOT + "NSHMP13" + S + "Local_13over8" + S + id + "-ratio";
					String diffFile = ROOT + "NSHMP13" + S + "Local_13over8" + S + id + "-diff";
					makeRatioPlotNSHMP(xyzRatio, spacing, locGrid.bounds(), ratioFile, "GM ratio", true);
					makeDiffPlotNSHMP(xyzDiff, spacing, locGrid.bounds(), diffFile, "GM diff", true);
				}
			}
		}
	}

	
	private static void makeLocalRatioMaps() {
		makeLocalRatioMaps("nshmp_SF-LA-13", true, "nshmp_SF-LA", true, "NSHM13-NSHM", 0.1);
//		makeLocalRatioMaps("UC33-SF-LA-08", false, "nshmp_SF-LA", true, "UC3_08-NSHM", 0.1);
//		makeLocalRatioMaps("UC33-SF-LA-13", false, "nshmp_SF-LA", true, "UC3_13-NSHM", 0.1);
	}
	
	// creates higher-res hazard ratio maps of the SF and LA regions using 8
	// branch averaged solutions for UCERF3
	private static void makeLocalRatioMaps(String overSrc, boolean overIsNSHM,
			String underSrc, boolean underIsNSHM, String outName, double scale) {
	
		double spacing = 0.05;
		String outDir = ROOT + "Hazard-Local" + S;
		boolean showFaults = true;
		
		List<TestGrid> grids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		for (ProbOfExceed pe : PEs) {
			for (Period p : periods) {
				for (TestGrid grid : grids) {

					CurveContainer over = (overIsNSHM) ? 
						loadLocalNSHMP(overSrc, grid, p, spacing) :
						loadLocalUC3(overSrc, grid, p, spacing);
					CurveContainer under = (underIsNSHM) ? 
						loadLocalNSHMP(underSrc, grid, p, spacing) :
						loadLocalUC3(underSrc, grid, p, spacing);

					// data
					GriddedRegion gr = grid.grid(spacing);
					GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(over, gr, pe);
					GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(under, gr, pe);

					// ratio map
					GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);
					String label = "Ratio " + pe + " " + p.getLabel() + " (g)";
					String dlDir = outDir + dlDirNameRatio(grid, outName, pe, "ratio", p);
					makeRatioPlot(xyz, spacing, grid.bounds(), dlDir, label,
						true, scale, true, showFaults);
					
					// diff map
					xyz = GeoDataSetMath.subtract(xyzOver, xyzUnder);
					label = "Diff " + pe + " " + p.getLabel() + " (g)";
					dlDir = outDir + dlDirNameRatio(grid, outName, pe, "diff", p);
					double diffScale = p.equals(GM0P20) ? 0.4 : 0.2; // larger scale for 5Hz
					makeDiffPlot(xyz, spacing, grid.bounds(), dlDir, label,
						diffScale, true, showFaults);
				}
			}
		}
	}
	
	private static void makeLocalHazardMaps() {
		makeLocalHazardMaps("nshmp_SF-LA-13", true, "NSHM13");
		makeLocalHazardMaps("nshmp_SF-LA", true, "NSHM");
		makeLocalHazardMaps("UC33-SF-LA-08", false, "UC3-08");
		makeLocalHazardMaps("UC33-SF-LA-13", false, "UC3-13");
	}
	
	// creates higher-res hazard ratio maps of the SF and LA regions using 8
	// branch averaged solutions for UCERF3
	private static void makeLocalHazardMaps(String src, boolean isNSHM,
			String outName) {
	
		double spacing = 0.05;
		String outDir = ROOT + "Hazard-Local" + S;
		boolean showFaults = true;
		
		List<TestGrid> grids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		for (ProbOfExceed pe : PEs) {
			for (Period p : periods) {
				for (TestGrid grid : grids) {

					CurveContainer cc = (isNSHM) ? 
						loadLocalNSHMP(src, grid, p, spacing) :
						loadLocalUC3(src, grid, p, spacing);

					// data
					GriddedRegion gr = grid.grid(spacing);
					GeoDataSet xyz = NSHMP_DataUtils.extractPE(cc, gr, pe);

					// ratio map
//					GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);
					double[] minmax = NSHMP_PlotUtils.getRange(p);
					GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
					String label = "GM " + pe + " " + p.getLabel() + " (g)";
					String dlDir = outDir + dlDirNameRatio(grid, outName, pe, null, p);

					makeMapPlot(xyz, spacing, grid.bounds(), dlDir, label,
						minmax[0], minmax[1], cpt, true, showFaults);
				}
			}
		}
	}

	// moves/renames pdf one level up and deletes other map files
	private static void savePDF(String dir) {
		try {
			File pdfFrom = new File(dir, "map.pdf");
			File pdfTo = new File(dir + ".pdf");
			Files.copy(pdfFrom, pdfTo);
			FileUtils.deleteDirectory(new File(dir));
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	private static String dlDirNameRatio(TestGrid tg, String outName,
			ProbOfExceed pe, String type, Period p) {
		StringBuffer sb = new StringBuffer();
		sb.append(tg == SAN_FRANCISCO ? "SF" : "LA").append("-");
		sb.append(outName).append("-");
		if (type != null) sb.append(type).append("-");
		sb.append(p.getLabel()).append("-");
		sb.append(pe == PE2IN50 ? "2in50" : "10in50");
		return sb.toString();
	}

	
	private static CurveContainer loadLocalUC3(String dir, TestGrid grid, Period p, double spacing) {
		File srcDir = new File(SRC + dir);
		return buildBrAvgCurveContainer(srcDir, p, grid, spacing);
	}
	
	private static CurveContainer loadLocalNSHMP(String dir, TestGrid grid, Period p, double spacing) {
		String srcPath = SRC + dir + S + grid.name() + S + p + S + "curves.csv";
		File srcFile = new File(srcPath);
		return CurveContainer.create(srcFile, grid, spacing);
	}
	
	
	// builds UC2 and UC3 ratio and difference maps for timeDep and Indep data
	// the curves sources for this are 50 year forecasts and are probabilistic,
	// not rate based.
	private static void buildUCtimeDependentMaps() {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
//		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		String path = SRC + "UC-TimeDep/";
		double prob = 0.02;
		
		// UC2
		for (Period p : periods) {
			GeoDataSet xyzOver = loadSingle(path + "UC2-TD", prob, grid, p, spacing);
			GeoDataSet xyzUnder = loadSingle(path + "UC2-TI", prob, grid, p, spacing);
			String id = p.getLabel() + "-" + String.format("%.0fin50", prob * 100);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "TimeDep/UC2_ratio-" + id;
			String diffDir = ROOT + "TimeDep/UC2_diff-" + id;
			makeRatioPlotNSHMP(ratio, spacing, grid.bounds(), ratioDir, "UC2 ratio " + id, false);
			makeDiffPlotNSHMP(diff, spacing, grid.bounds(), diffDir, "UC2 diff " + id, false);
		}

		// UC3
		for (Period p : periods) {
			GeoDataSet xyzOver = loadSingle(path + "UC3-FM31-TD", prob, grid, p, spacing);
			GeoDataSet xyzUnder = loadSingle(path + "UC3-FM31-TI", prob, grid, p, spacing);
			String id = p.getLabel() + "-" + String.format("%.0fin50", prob * 100);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "TimeDep/UC3_ratio-" + id;
			String diffDir = ROOT + "TimeDep/UC3_diff-" + id;
			makeRatioPlotNSHMP(ratio, spacing, grid.bounds(), ratioDir, "UC3 ratio " + id, false);
			makeDiffPlotNSHMP(diff, spacing, grid.bounds(), diffDir, "UC3 diff " + id, false);
		}

	}
	
	
	// SRC13-GMM13 over SRC08-GMM08
	private static void build_NSHMP_comp1() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// nshmp data
				GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca", pe, grid, p, spacing);
				String mapID = p.getLabel() + "-" + pe.name();
				nshmpXYZ.put(mapID,  xyzNSHMP);
				// nshmp hazard map
//				String dir = ROOT + "Hazard/nshmp-" + mapID;
//				makeHazardMap(xyzNSHMP, spacing, p, pe, grid, dir);
			}
		}

		// UC3 8brAvgMaps
		TestGrid nshmpGrid = TestGrid.CA_NSHMP;
		Map<String, GeoDataSet> brAvgXYZ = Maps.newHashMap();
		String srcBase = SRC + "NSHMP13-8-prob";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, nshmpGrid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "NSHMP13/Hazard-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 / NSHMP
		for (String id : brAvgXYZ.keySet()) {
			GeoDataSet xyzUnder = nshmpXYZ.get(id);
			GeoDataSet xyzOver = brAvgXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/SRC13-GMM13_sup_SRC08-GMM08-" + id;
			String diffDir = ROOT + "NSHMP13/SRC13-GMM13_diff_SRC08-GMM08-" + id;
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", false);
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "GM diff", false);
		}

	}
	
	// SRC13-GMM08 over SRC08-GMM08
	private static void build_NSHMP_comp2() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// nshmp data
				GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca", pe, grid, p, spacing);
				String mapID = p.getLabel() + "-" + pe.name();
				nshmpXYZ.put(mapID,  xyzNSHMP);
				// nshmp hazard map
//				String dir = ROOT + "Hazard/nshmp-" + mapID;
//				makeHazardMap(xyzNSHMP, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 8brAvgMaps
		Map<String, GeoDataSet> brAvgXYZ = Maps.newHashMap();
		String srcBase = SRC + "UC33brAvg-FM-DM-08" + S + "all";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "Hazard/brAvg08-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 / NSHMP
		for (String id : brAvgXYZ.keySet()) {
			GeoDataSet xyzNSHMP = nshmpXYZ.get(id);
			GeoDataSet xyzBrAvg = brAvgXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzBrAvg, xyzNSHMP);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzBrAvg, xyzNSHMP);
//			String ratioDir = ROOT + "Hazard/brAvg13_sup_nshmp-" + id;
//			String diffDir = ROOT + "Hazard/brAvg13_diff_nshmp-" + id;
			String ratioDir = ROOT + "NSHMP13/SRC13-GMM08_sup_SRC08-GMM08-" + id;
			String diffDir = ROOT + "NSHMP13/SRC13-GMM08_diff_SRC08-GMM08-" + id;
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", true);
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "GM diff", true);
		}
		
	}


	// SRC08-GMM13 over SRC08-GMM08
	private static void build_NSHMP_comp3() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				GeoDataSet over = loadSingle(SRC + "nshmp_ca-13/", pe, grid, p, spacing);
				GeoDataSet under = loadSingle(SRC + "nshmp_ca/", pe, grid, p, spacing);

				GeoDataSet xyzRatio = GeoDataSetMath.divide(over, under);
				GeoDataSet xyzDiff = GeoDataSetMath.subtract(over, under);
				
				String id = p.getLabel() + "-" + pe.name();
				
				String ratioDir = ROOT + "NSHMP13/SRC08-GMM13_sup_SRC08-GMM08-" + id;
				String diffDir = ROOT + "NSHMP13/SRC08-GMM13_diff_SRC08-GMM08-" + id;
				makeRatioPlotNSHMP(xyzRatio, 0.1, grid.bounds(), ratioDir, "GM ratio", true);
				makeDiffPlotNSHMP(xyzDiff, 0.1, grid.bounds(), diffDir, "GM diff", true);
			}
		}
	}
	
	// SRC13-GMM13 over SRC13-GMM08
	private static void build_NSHMP_comp4() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		// SRC13 GMM13
		TestGrid nshmpGrid = TestGrid.CA_NSHMP;
		Map<String, GeoDataSet> overXYZ = Maps.newHashMap();
		String srcBase = SRC + "NSHMP13-8-prob";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, nshmpGrid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				overXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "NSHMP13/Hazard-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}

		// SRC13 GMM08
		Map<String, GeoDataSet> underXYZ = Maps.newHashMap();
		srcBase = SRC + "UC33brAvg-FM-DM-08" + S + "all";
		gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				underXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "Hazard/brAvg08-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// GMPE diff & ratio with 14 source model
		for (String id : overXYZ.keySet()) {
			GeoDataSet xyzOver = overXYZ.get(id);
			GeoDataSet xyzUnder = underXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/SRC13-GMM13_sup_SRC13-GMM08-" + id;
			String diffDir = ROOT + "NSHMP13/SRC13-GMM13_diff_SRC13-GMM08-" + id;
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", true);
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "GM diff", true);
		}

	}

	// build ratio and diff maps of NSHMP-13 and nshmp-08 for bg and flt sources
	private static void build_NSHMP_BgFltMaps() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		List<String> srcs = Lists.newArrayList("bg", "flt");
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmp08XYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				for (String src : srcs) {
					// nshmp data
					GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca_" + src, pe, grid, p, spacing);
					String mapID = p.getLabel() + "-" + pe.name() + "-" + src;
					nshmp08XYZ.put(mapID,  xyzNSHMP);
				}
			}
		}
		
		// UC3 faults 8brAvgMaps FM-DM
		Map<String, GeoDataSet> nshmp13XYZ = Maps.newHashMap();
		TestGrid nshmpGrid = CA_NSHMP;
		String detBase = SRC + "NSHMP13-8-det";
		String probBase = SRC + "NSHMP13-8-prob";
		for (Period p : periods) {
			File detSrcDir = new File(detBase);
			File probSrcDir = new File(probBase);
			CurveContainer fltCC = buildBrAvgCurveContainer(detSrcDir, p, nshmpGrid, spacing);
			CurveContainer bgCC = buildBrAvgCurveContainer(probSrcDir, p, nshmpGrid, spacing);
			bgCC.subtract(fltCC);
			for (ProbOfExceed pe : PEs) {
				GeoDataSet fltXYZ = NSHMP_DataUtils.extractPE(fltCC, grid.grid(spacing), pe);
				String mapID = p.getLabel() + "-" + pe.name() + "-flt";
				nshmp13XYZ.put(mapID,  fltXYZ);
				GeoDataSet bgXYZ = NSHMP_DataUtils.extractPE(bgCC, grid.grid(spacing), pe);
				mapID = p.getLabel() + "-" + pe.name() + "-bg";
				nshmp13XYZ.put(mapID,  bgXYZ);
			}
		}
		
		// NSHMP13 / NSHMP08
		for (String id : nshmp08XYZ.keySet()) {
			GeoDataSet xyzNSHMP08 = nshmp08XYZ.get(id);
			GeoDataSet xyzNSHMP13 = nshmp13XYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzNSHMP13, xyzNSHMP08);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzNSHMP13, xyzNSHMP08);
			String ratioDir = ROOT + "NSHMP13" + S + "BG-FLT" + S + id + "-ratio";
			String diffDir = ROOT + "NSHMP13" + S + "BG-FLT" + S + id + "-diff";
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", false);
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "GM diff", false);
		}
		
	}
	
	// for comparing 'release' UC3-NSHMp runs to earlier, 'practice' runs
	private static void buildUC3_NSHMP_versionMaps() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		// UC3 8brAvgMaps
		TestGrid nshmpGrid = TestGrid.CA_NSHMP;
		Map<String, GeoDataSet> brAvgXYZnew = Maps.newHashMap();
//		String srcBase = SRC + "UC33brAvg-FM-DM-13" + S + "all";
		String srcBase = SRC + "NSHMP13-8-prob";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, nshmpGrid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZnew.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "Hazard/brAvg13-" + mapID;
//				String dir = ROOT + "HazardNSHMP/UC3-NSHMP-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 8brAvgMaps
		Map<String, GeoDataSet> brAvgXYZold = Maps.newHashMap();
		srcBase = SRC + "UC33brAvg-FM-DM-13" + S + "all";
//		srcBase = SRC + "UC33brAvg-FM-DM-08" + S + "all";
		gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZold.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "Hazard/brAvg13-" + mapID;
//				String dir = ROOT + "Hazard/brAvg08-" + mapID;
//				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		

		
		// UC3 / NSHMP
		for (String id : brAvgXYZold.keySet()) {
			GeoDataSet xyzBrAvgOld = brAvgXYZold.get(id);
			GeoDataSet xyzBrAvgNew = brAvgXYZnew.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzBrAvgNew, xyzBrAvgOld);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzBrAvgNew, xyzBrAvgOld);
//			String ratioDir = ROOT + "Hazard/brAvg13_sup_nshmp-" + id;
//			String diffDir = ROOT + "Hazard/brAvg13_diff_nshmp-" + id;
			String ratioDir = ROOT + "HazardNSHMP/UC3-NSHMP_final_sup_beta-" + id;
			String diffDir = ROOT + "HazardNSHMP/UC3-NSHMP_final_diff_beta-" + id;
			makeRatioPlot(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.3, true, false);
			double diffScale = id.contains(GM0P20.getLabel()) ? 0.4 : 0.2; // larger scale for 5Hz
			makeDiffPlot(diff, 0.1, grid.bounds(), diffDir, "GM diff", diffScale, true, false);
		}

	}

	
	
	// combines 8 branch averaged solutions for faults (computed with
	// overlapping WUS fault removed, eg.g Carson) with 2 bravg solutions for
	// background (computed masking all those grid nodes outside CA)
	private static void buildUC3_NSHMP_binaries8() throws IOException {
		TestGrid grid = TestGrid.CA_NSHMP;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			// separate fault and bg files
//			File fltDir = new File(SRC + "NSHMP13-det-flt-8");
//			CurveContainer fltCC = buildBrAvgCurveContainer(fltDir, p, grid, spacing);
//			File bgDir = new File(SRC + "NSHMP13-bg");
//			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, grid, spacing);
//			fltCC.add(bgCC);
			
			// single files
			File curveDir = new File(SRC + "NSHMP13-8-prob");
			CurveContainer fltCC = buildBrAvgCurveContainer(curveDir, p, grid, spacing);

//			GeoDataSet xyz = NSHMP_DataUtils.extractPE(fltCC, grid.grid(spacing), PE2IN50);
//			makeHazardMap(xyz, spacing, p, PE2IN50, grid, ROOT + "test/NSHMP-" + p.getLabel());

			File out = new File("tmp/NSHMP-CA/binaries_09-06-2013/CA-" + p.getLabel() + ".curves");
			Files.createParentDirs(out);
			
			BinaryCurves.writeUC3(fltCC, p, spacing, "UCERF3.3 Sources", out);
		}
	}
	
	// converts curves build using MeanUCERF3 to NSHMP binary format
	private static void buildUC3_NSHMP_binaries40() throws IOException {
		TestGrid grid = TestGrid.CA_NSHMP;
		double spacing = 0.05;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			// separate fault and bg files
			File fltDir = new File(SRC + "NSHMP13-40-FLT");
			CurveContainer fltCC = buildBrAvgCurveContainer(fltDir, p, grid, spacing);
			File bgDir = new File(SRC + "NSHMP13-2-BG");
			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, grid, spacing);

			// make a map on the way if desired
//			GeoDataSet xyz = NSHMP_DataUtils.extractPE(fltCC, grid.grid(spacing), PE2IN50);
//			makeHazardMap(xyz, spacing, p, PE2IN50, TestGrid.CA_RELM, ROOT + "test/NSHMP-" + p.getLabel());

			// bg and flt separate
//			File outFlt = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-FLT-0p05.curves");
//			File outGrd = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-GRD-0p1.curves");
//			Files.createParentDirs(outFlt);
//			BinaryCurves.writeUC3(fltCC, p, 0.05, "UCERF3.3 Fault Sources", outFlt);
//			BinaryCurves.writeUC3(bgCC, p, 0.1, "UCERF3.3 Grid Sources", outGrd);
			
			// all together now
//			fltCC.add(bgCC);

//			GeoDataSet xyz = NSHMP_DataUtils.extractPE(fltCC, grid.grid(spacing), PE2IN50);
//			makeHazardMap(xyz, spacing, p, PE2IN50, TestGrid.CA_RELM, ROOT + "test/NSHMP-" + p.getLabel());

//			File outAll = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-ALL-0p05.curves");
			File outFlt = new File("tmp/morgan/curves/CA13-" + p.getLabel() + "-FLT-0p05.curves");
			BinaryCurves.writeUC3(fltCC, p, spacing, "UCERF3.3 FLT Sources", outFlt);
			File outBg = new File("tmp/morgan/curves/CA13-" + p.getLabel() + "-BG-0p05.curves");
			BinaryCurves.writeUC3(bgCC, p, spacing, "UCERF3.3 BG Sources", outBg);
		}
	}
	
	// flt and bg curves for combination (for figures only)
	private static void buildUC3_NSHMP_binariesMorgan() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			// separate fault and bg files
			File fltDir = new File(SRC + "NSHMP13-det-flt-8");
			CurveContainer fltCC = buildBrAvgCurveContainer(fltDir, p, grid, spacing);
			File bgDir = new File(SRC + "NSHMP13-bg");
			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, grid, spacing);

			File outFlt = new File("tmp/morgan/curves/CA13-" + p.getLabel() + "-FLT-0p1.curves");
			BinaryCurves.writeUC3(fltCC, p, spacing, "UCERF3.3 FLT Sources", outFlt);
			File outBg = new File("tmp/morgan/curves/CA13-" + p.getLabel() + "-BG-0p1.curves");
			BinaryCurves.writeUC3(bgCC, p, spacing, "UCERF3.3 BG Sources", outBg);
		}
	}

	// flt and bg curves for combination (for figures only)
	private static void buildUC3_NSHMP_binariesMorgan2() throws IOException {
		String sol = "mean_ucerf3_sol";
		TestGrid grid = TestGrid.CA_NSHMP;
		double fltSpacing = 0.05;
		double bgSpacing = 0.1;
		List<Period> periods = Lists.newArrayList(Period.GM0P30); //GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			// separate fault and bg files
			File fltDir = new File(SRC + "NSHMP14" + S + "FLT_0.05" + S + sol + S + grid + S + p + S + "curves.csv");
			CurveContainer fltCC = CurveContainer.create(fltDir, grid, fltSpacing);
			File bgDir = new File(SRC + "NSHMP14" + S + "BG" + S + sol + S + grid + S + p + S + "curves.csv");
			CurveContainer bgCC = CurveContainer.create(bgDir, grid, bgSpacing);

			File outFlt = new File("tmp/morgan/curves/CA14-" + p.getLabel() + "-FLT-0p05.curves");
			BinaryCurves.writeUC3(fltCC, p, fltSpacing, "UCERF3.3 FLT Sources 0.05deg", outFlt);
			File outBg = new File("tmp/morgan/curves/CA14-" + p.getLabel() + "-BG-0p1.curves");
			BinaryCurves.writeUC3(bgCC, p, bgSpacing, "UCERF3.3 BG Sources 0.1deg", outBg);
		}
	}
	
	// binaries of the two different background models (UC2 vs UC3) using the 'reference' branch
	private static void buildUC3_NSHMP_BG_binariesMorgan() throws IOException {
		TestGrid grid = TestGrid.CA_NSHMP;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		List<String> bgModels = Lists.newArrayList("UC2", "UC3");
		
		for (String bgModel : bgModels) {
			for (Period p : periods) {
				// separate fault and bg files
				File fltDir = new File(SRC + "NSHMP14-BG-Morgan" + S + bgModel + S + grid + S + p + S + "curves.csv");
				CurveContainer fltCC = CurveContainer.create(fltDir, grid, spacing);
	
				File out = new File("tmp/morgan/curves/CA14-BG-" + bgModel + "-" + p.getLabel() + "-0p1.curves");
				BinaryCurves.writeUC3(fltCC, p, spacing, "UCERF3.3 " + bgModel + " smooth seis", out);
			}
		}
	}


	
	private static void buildUC2_NSHMP_binaries() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			// separate fault and bg files
			CurveContainer fltCC = loadLocalNSHMP("nshmp_ca_flt", grid, p, spacing);
			File outFlt = new File("tmp/morgan/curves/CA08-" + p.getLabel() + "-FLT-0p1.curves");
			BinaryCurves.writeUC3(fltCC, p, spacing, "UCERF2 FLT Sources", outFlt);

			CurveContainer bgCC = loadLocalNSHMP("nshmp_ca_bg", grid, p, spacing);
			File outBg = new File("tmp/morgan/curves/CA08-" + p.getLabel() + "-BG-0p1.curves");
			BinaryCurves.writeUC3(bgCC, p, spacing, "UCERF2 BG Sources", outBg);
		}
	}
	
	// single brAveraged fault model computed using 2008 GMMs
	// for Morgan M. for WUS combination
	private static void build_UC3_GMM08_binary() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		for (Period p : periods) {
			// separate fault and bg files
			String dir = "NSHMP13-GMPE08" + S + "ALL" + S + "UC33brAvg_FM31";
			CurveContainer cc = loadLocalNSHMP(dir, grid, p, spacing);
			File outFlt = new File("tmp/morgan/curves/CA-SRC14-GMM08-" + p.getLabel() + "-ALL-0p1.curves");
			BinaryCurves.writeUC3(cc, p, spacing, "NSHMP14 with GMM08", outFlt);
		}
	}

	
	

	// combines 40 branch averaged solutions for faults (computed with
	// overlapping WUS fault removed, eg. Carson) with 2 bravg solutions for
	// background (computed masking all those grid nodes outside CA)
	private static void buildUC3_NSHMP_binaries() throws IOException {
		TestGrid grid = TestGrid.CA_NSHMP;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM3P00); //GM0P00, GM0P20, GM1P00, GM1P50, GM2P00, GM3P00);
		
		for (Period p : periods) {
			File f = new File(SRC + "NSHMP13B-epi" + S + "mean_ucerf3_sol" + S + grid + S + p + S + "curves.csv");
			CurveContainer cc = CurveContainer.create(f, grid, spacing);
			File out = new File("tmp/NSHMP-CA/binaries_02-07-2014/CA-" + p.getLabel() + "-0p1.curves");
			BinaryCurves.writeUC3(cc, p, 0.1, "UCERF3.3 All Sources", out);
		}
	}
	

	private static void finalBSSCcheck() throws IOException {
		double spacing = 0.1;
		TestGrid dataGrid = TestGrid.CA_NSHMP;
		TestGrid mapGrid = TestGrid.CA_RELM;
		GriddedRegion gr = mapGrid.grid(spacing);
		ProbOfExceed pe =  ProbOfExceed.PE2IN50;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			
			// feb 2014 vs sept 2013
			// early run, sept 13 -- separate fault and bg files
			File fltDir = new File(SRC + "NSHMP13-40-FLT");
			CurveContainer fltCC = buildBrAvgCurveContainer(fltDir, p, dataGrid, 0.05);
			File bgDir = new File(SRC + "NSHMP13-2-BG");
			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, dataGrid, 0.05);
			fltCC.add(bgCC);
			GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(fltCC, gr, pe);
	
			// early run, feb 14
			File dir = new File(SRC + "NSHMP13B-epi" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
			CurveContainer ccOver = CurveContainer.create(dir, dataGrid, spacing);
			GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(ccOver, gr, pe);

			String id = p.getLabel() + "-" + pe.name();
			String ratioDir = ROOT + "NSHMP14-BSSC" + S + "NSHMP13-FebVsSept_ratio_" + id;
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), ratioDir, "NSHMP13-FebVsSept_" + id, true, 0.02, true, false);

			// april 2014 vs feb 2014
//			// early run, feb 14
//			File dir = new File(SRC + "NSHMP13B-epi" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
//			CurveContainer ccUnder = CurveContainer.create(dir, dataGrid, spacing);
//			GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(ccUnder, gr, pe);
//					
//			// recent run, April
//			File bgDir = new File(SRC + "NSHMP14" + S + "BG" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
//			CurveContainer ccOver = CurveContainer.create(bgDir, dataGrid, spacing);
//			File fltDir = new File(SRC + "NSHMP14" + S + "FLT" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
//			CurveContainer fltCC = CurveContainer.create(fltDir, dataGrid, spacing);
//			ccOver.add(fltCC);
//			GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(ccOver, gr, pe);
//			
//			String id = p.getLabel() + "-" + pe.name();
//			String ratioDir = ROOT + "NSHMP14-BSSC" + S + "NSHMP14sup13_ratio_" + id;
//			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), ratioDir, "NSHMP14_sup_13_" + id, true, 0.01, true, false);
		}

	}
	
	/*
	 * Tracking down changes to the maps 
	 */
	private static void finalMapsDebug() throws IOException {
		double spacing = 0.1;
		TestGrid dataGrid = TestGrid.CA_NSHMP;
		TestGrid mapGrid = TestGrid.CA_RELM;
		GriddedRegion gr = mapGrid.grid(spacing);
		ProbOfExceed pe =  ProbOfExceed.PE2IN50;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);

		for (Period p : periods) {
			
			String id = p.getLabel() + "-" + pe.name();
			String baseDir = ROOT + "NSHMP14-BSSC-debug" + S;
			String binBase = baseDir + "bin" + S;

			// 1: 9/6/13 -- only 8 branch averaged solutions
			File curveDir = new File(SRC + "NSHMP13-8-prob");
			CurveContainer cc_9_6_13 = buildBrAvgCurveContainer(curveDir, p, dataGrid, spacing);
			GeoDataSet xyz_9_6_13 = NSHMP_DataUtils.extractPE(cc_9_6_13, gr, pe);
			File binOut = new File(binBase + "1_9-6-13_" + p.getLabel() + "_0p1.curves");
			BinaryCurves.writeUC3(cc_9_6_13, p, spacing, "UCERF3.3 9-6-13", binOut);
			
			// 2: 9/13/13 -- separate fault and bg files, 0.05 grid
			File f_9_13_13_flt = new File(SRC + "NSHMP13-40-FLT");
			CurveContainer cc_9_13_13 = buildBrAvgCurveContainer(f_9_13_13_flt, p, dataGrid, 0.05);
			File f_9_13_13_bg = new File(SRC + "NSHMP13-2-BG");
			CurveContainer cc_9_13_13_bg = buildBrAvgCurveContainer(f_9_13_13_bg, p, dataGrid, 0.05);
			cc_9_13_13.add(cc_9_13_13_bg);
			GeoDataSet xyz_9_13_13 = NSHMP_DataUtils.extractPE(cc_9_13_13, gr, pe);
			binOut = new File(binBase + "2_9-13-13_" + p.getLabel() + "_0p1.curves");
			BinaryCurves.writeUC3(cc_9_6_13, p, spacing, "UCERF3.3 9-13-13", binOut);
			
			// 3: 2/7/14 -- switch to mean solution,  0.1 grid
			File f_2_7_14 = new File(SRC + "NSHMP13B-epi" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
			CurveContainer cc_2_7_14 = CurveContainer.create(f_2_7_14, dataGrid, spacing);
			GeoDataSet xyz_2_7_14 = NSHMP_DataUtils.extractPE(cc_2_7_14, gr, pe);
			binOut = new File(binBase + "2_2-7-14_" + p.getLabel() + "_0p1.curves");
			BinaryCurves.writeUC3(cc_2_7_14, p, spacing, "UCERF3.3 2-7-14", binOut);
			
			// 4: 5/2/14 -- most recent run
			File f_5_2_14_flt = new File(SRC + "NSHMP14" + S + "FLT" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
			CurveContainer cc_5_2_14 = CurveContainer.create(f_5_2_14_flt, dataGrid, spacing);
			File f_5_2_14_bg = new File(SRC + "NSHMP14" + S + "BG" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
			CurveContainer cc_5_2_14_bg = CurveContainer.create(f_5_2_14_bg, dataGrid, spacing);
			cc_5_2_14.add(cc_5_2_14_bg);
			GeoDataSet xyz_5_2_14 = NSHMP_DataUtils.extractPE(cc_5_2_14, gr, pe);
			binOut = new File(binBase + "2_5-2-14_" + p.getLabel() + "_0p1.curves");
			BinaryCurves.writeUC3(cc_5_2_14, p, spacing, "UCERF3.3 5-2-14", binOut);

						
			// 4/1
			String name = "4over1_" + id;
			GeoDataSet ratio = GeoDataSetMath.divide(xyz_5_2_14, xyz_9_6_13);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
			
			// 4/2
			name = "4over2_" + id;
			ratio = GeoDataSetMath.divide(xyz_5_2_14, xyz_9_13_13);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);

			// 4/3
			name = "4over3_" + id;
			ratio = GeoDataSetMath.divide(xyz_5_2_14, xyz_2_7_14);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
			
			// 2/1
			name = "2over1_" + id;
			ratio = GeoDataSetMath.divide(xyz_9_13_13, xyz_9_6_13);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
			
			// 3/2
			name = "3over2_" + id;
			ratio = GeoDataSetMath.divide(xyz_2_7_14, xyz_9_13_13);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
			
			// 4/3
			name = "4over3_" + id;
			ratio = GeoDataSetMath.divide(xyz_5_2_14, xyz_2_7_14);
			makeRatioPlot(ratio, 0.1, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
		}
	}
	
	// comparing new 0.05 spaced fault curves to those from Sept 2013
	private static void finalMapsSpacingComparison() throws IOException {
		double spacing = 0.05;
		TestGrid dataGrid = TestGrid.CA_NSHMP;
		TestGrid mapGrid = TestGrid.CA_RELM;
		GriddedRegion gr = dataGrid.grid(spacing);
		ProbOfExceed pe =  ProbOfExceed.PE2IN50;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
	
		for (Period p : periods) {
			String id = p.getLabel() + "-" + pe.name();
			String baseDir = ROOT + "NSHMP14-BSSC-final" + S;
			String binBase = baseDir + "bin" + S;

			// 1: 9/13/13 -- separate fault and bg files, 0.05 grid
			File f_9_13_13_flt = new File(SRC + "NSHMP13-40-FLT");
			CurveContainer cc_9_13_13 = buildBrAvgCurveContainer(f_9_13_13_flt, p, dataGrid, spacing);
			GeoDataSet xyz_9_13_13 = NSHMP_DataUtils.extractPE(cc_9_13_13, gr, pe);

			// 2: 5/21/14 -- most recent run
			File f_5_21_14_flt = new File(SRC + "NSHMP14" + S + "FLT_0.05" + S + "mean_ucerf3_sol" + S + dataGrid + S + p + S + "curves.csv");
			CurveContainer cc_5_21_14 = CurveContainer.create(f_5_21_14_flt, dataGrid, spacing);
			GeoDataSet xyz_5_21_14 = NSHMP_DataUtils.extractPE(cc_5_21_14, gr, pe);
			File binOut = new File(binBase + "FLT_5-21-14_" + p.getLabel() + "_0p05.curves");
			BinaryCurves.writeUC3(cc_5_21_14, p, spacing, "UCERF3.3 5-21-14", binOut);

			// 2/1
			String name = "current/Sept" + id;
			GeoDataSet ratio = GeoDataSetMath.divide(xyz_5_21_14, xyz_9_13_13);
			makeRatioPlot(ratio, spacing, mapGrid.bounds(), baseDir + name, name, true, 0.02, true, false);
			
		}
	}
		
	
//	private static void resamplingTest() throws IOException {
//		TestGrid grid = CA_RELM;
//		double spacing = Double.NaN;
//		List<Period> periods = Lists.newArrayList(GM0P20); //GM0P00, GM0P20, GM1P00);
//		
//		for (Period p : periods) {
//			// separate fault and bg files
////			File fltDir = new File(SRC + "NSHMP13-40-FLT");
////			CurveContainer fltCC = buildBrAvgCurveContainer(fltDir, p, grid, spacing);
//			
//			// build bg container from 0.05 spacing file
//			spacing = 0.05;
//			File bgDir = new File(SRC + "NSHMP13-2-BG");
//			CurveContainer bgCC = buildBrAvgCurveContainer(bgDir, p, CA_NSHMP, spacing);
//			GeoDataSet xyz = NSHMP_DataUtils.extractPE(bgCC, grid.grid(spacing), PE2IN50);
//			makeHazardMap(xyz, spacing, p, PE2IN50, grid, ROOT + "test/resampTest-BGnoResamp-" + p.getLabel());
//
//			// resample container - upsizing grid to NSHMP
//			spacing = 0.1;
//			// bg coming in at high resolution but being downsampled when
//			// added to 
//			CurveContainer ccResamp = BinaryCurves.createNSHMP(bgCC, spacing, p);
//			xyz = NSHMP_DataUtils.extractPE(ccResamp, grid.grid(spacing), PE2IN50);
//			makeHazardMap(xyz, 0.1, p, PE2IN50, grid, ROOT + "test/resampTest-BGresamp-" + p.getLabel());
//
//			// bg and flt separate
////			File outFlt = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-FLT-0p05.curves");
////			File outGrd = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-GRD-0p1.curves");
////			Files.createParentDirs(outFlt);
////			BinaryCurves.writeUC3(fltCC, p, 0.05, "UCERF3.3 Fault Sources", outFlt);
////			BinaryCurves.writeUC3(bgCC, p, 0.1, "UCERF3.3 Grid Sources", outGrd);
//			
//			// all together now
////			fltCC.add(bgCC);
//
////			GeoDataSet xyz = NSHMP_DataUtils.extractPE(fltCC, grid.grid(spacing), PE2IN50);
////			makeHazardMap(xyz, spacing, p, PE2IN50, TestGrid.CA_RELM, ROOT + "test/NSHMP-" + p.getLabel());
//
////			File outAll = new File("tmp/NSHMP-CA/binaries_09-20-2013/CA-" + p.getLabel() + "-ALL-0p05.curves");
////			BinaryCurves.writeUC3(fltCC, p, 0.05, "UCERF3.3 All Sources", outAll);
//			
//		}
//	}

	
	// combine different fault based deterministic outputs
	private static void buildUC3_NSHMP_determ() throws IOException {
		TestGrid grid = TestGrid.CA_NSHMP;
		
//		String root = SRC + "NSHMP13-det-flt-40" + S;
//		List<String> dirs = Lists.newArrayList(
//			"UC33brAvg_FM31_ZENGBB_SH09M",
//			"UC33brAvg_FM32_ZENGBB_SH09M");
//		List<Period> periods =  Lists.newArrayList(GM0P00, GM0P20, GM1P00);
//		
//		for (Period p : periods) {
//			List<File> files = Lists.newArrayList();
//			for (String dir : dirs) {
//				files.add(new File(root + dir + S + grid + S + p + S + "determ.txt"));
//			}
//			File out = new File(root + "RefBranchDet-" + p.getLabel() + ".txt");
//			createMaxDeterm(files, out);
//		}
		
//		String root = SRC + "NSHMP13-det-flt-8" + S;
//		List<String> dirs = Lists.newArrayList(
//			"UC33brAvg_FM31_GEOL",
//			"UC33brAvg_FM31_ZENGBB",
//			"UC33brAvg_FM31_NEOK",
//			"UC33brAvg_FM32_GEOL",
//			"UC33brAvg_FM32_ZENGBB",
//			"UC33brAvg_FM32_NEOK");

		String root = SRC + "NSHMP13-40-FLT" + S;
		List<String> dirs = Lists.newArrayList(
			"UC33brAvg_FM31_GEOL_ELLB",
			"UC33brAvg_FM31_GEOL_ELLBSL",
			"UC33brAvg_FM31_GEOL_HB08",
			"UC33brAvg_FM31_GEOL_SH09M",
			"UC33brAvg_FM31_GEOL_SHCSD",
			"UC33brAvg_FM31_NEOK_ELLB",
			"UC33brAvg_FM31_NEOK_ELLBSL",
			"UC33brAvg_FM31_NEOK_HB08",
			"UC33brAvg_FM31_NEOK_SH09M",
			"UC33brAvg_FM31_NEOK_SHCSD",
			"UC33brAvg_FM31_ZENGBB_ELLB",
			"UC33brAvg_FM31_ZENGBB_ELLBSL",
			"UC33brAvg_FM31_ZENGBB_HB08",
			"UC33brAvg_FM31_ZENGBB_SH09M",
			"UC33brAvg_FM31_ZENGBB_SHCSD",
			"UC33brAvg_FM32_GEOL_ELLB",
			"UC33brAvg_FM32_GEOL_ELLBSL",
			"UC33brAvg_FM32_GEOL_HB08",
			"UC33brAvg_FM32_GEOL_SH09M",
			"UC33brAvg_FM32_GEOL_SHCSD",
			"UC33brAvg_FM32_NEOK_ELLB",
			"UC33brAvg_FM32_NEOK_ELLBSL",
			"UC33brAvg_FM32_NEOK_HB08",
			"UC33brAvg_FM32_NEOK_SH09M",
			"UC33brAvg_FM32_NEOK_SHCSD",
			"UC33brAvg_FM32_ZENGBB_ELLB",
			"UC33brAvg_FM32_ZENGBB_ELLBSL",
			"UC33brAvg_FM32_ZENGBB_HB08",
			"UC33brAvg_FM32_ZENGBB_SH09M",
			"UC33brAvg_FM32_ZENGBB_SHCSD");
		
		List<Period> periods =  Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		
		for (Period p : periods) {
			List<File> files = Lists.newArrayList();
			for (String dir : dirs) {
				files.add(new File(root + dir + S + grid + S + p + S + "determ.txt"));
			}
			File out = new File(root + "MaxDet-" + p.getLabel() + ".txt");
			createMaxDeterm(files, out);
		}
	}
	
	private static void createMaxDeterm(List<File> files, File out) throws IOException {
		// build maps with first file
		
		Map<String, String> lineMap = Maps.newHashMap();
		Map<String, Double> valMap = Maps.newHashMap();

		for (File f : files) {
			List<String> lines = Files.readLines(f, US_ASCII);
			for (String line : lines) {
				Iterable<String> it = SPLIT_TAB.split(line);
				String key = Iterables.get(it, 0) + "_" + Iterables.get(it, 1);
				double median = Double.valueOf(Iterables.get(it, 2));
				if (valMap.containsKey(key)) {
					double refMedian = valMap.get(key);
					if (median > refMedian) {
						lineMap.put(key, line);
						valMap.put(key, median);
					}
				} else {
					lineMap.put(key, line);
					valMap.put(key, median);
				}
			}
		}
		
		Files.write("", out, US_ASCII);
		for (String key : lineMap.keySet()) {
			Files.append(lineMap.get(key) + LF, out, US_ASCII);
		}
	}
	
	// compare pairs of nshmp maps
	// used to check ratio of optimized NSHMP to new grid implementation
	// 2008 GMPEs to 2013 GMPEs
	private static void buildNSHMPratioDiff() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
//		String outDir = ROOT + "NewGMPE" + S;
//		String outDir = ROOT + "NewGrid" + S;
		String outDir = ROOT + "NewGridGMPE" + S;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
//				GeoDataSet over = loadSingle(SRC + "nshmp_ca-13/", pe, grid, p, spacing);
//				GeoDataSet under = loadSingle(SRC + "nshmp_ca-08/", pe, grid, p, spacing);
			
//				GeoDataSet over = loadSingle(SRC + "nshmp_ca-08/", pe, grid, p, spacing);
//				GeoDataSet under = loadSingle(SRC + "nshmp_ca/", pe, grid, p, spacing);

				GeoDataSet over = loadSingle(SRC + "nshmp_ca-13/", pe, grid, p, spacing);
				GeoDataSet under = loadSingle(SRC + "nshmp_ca/", pe, grid, p, spacing);

				GeoDataSet xyzRatio = GeoDataSetMath.divide(over, under);
				GeoDataSet xyzDiff = GeoDataSetMath.subtract(over, under);
			
//				String id = "NSHMP-13-08";
//				String id = "NSHMP-newGrid";
				String id = "NSHMP-13-opt";
				
				String ratioDir = outDir + id + "-ratio-" + p.getLabel() + "-" + pe;
				makeRatioPlot(xyzRatio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.3, true, false);
				String diffDir = outDir +  id + "-diff-" + p.getLabel() + "-" + pe;
				double diffScale = p.equals(GM0P20) ? 0.4 : 0.2; // larger scale for 5Hz
				makeDiffPlot(xyzDiff, 0.1, grid.bounds(), diffDir, "GM diff", diffScale, true, false);
			}
		}
	}
	
	// compare pairs of UC3 maps
	// used to check ratio of original UC3 maps to new grid implementation
	// (originals used pt2vertPoiss)
	// 2008 GMPEs to 2013 GMPEs
	private static void buildUC3ratioDiff() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
//		String outDir = ROOT + "NewGMPE" + S;
//		String outDir = ROOT + "NewGrid" + S;
		String outDir = ROOT + "NewGridGMPE" + S;

		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		GriddedRegion gr = grid.grid(spacing);

//		String srcOver = SRC + "UC33brAvg-FM-DM-13" + S + "all";
//		String srcUnder = SRC + "UC33brAvg-FM-DM-08" + S + "all";
		
//		String srcOver = SRC + "UC33brAvg-FM-DM-08" + S + "all";
//		String srcUnder = SRC + "UC33brAvg-FM-DM";

		String srcOver = SRC + "UC33brAvg-FM-DM-13" + S + "all";
		String srcUnder = SRC + "UC33brAvg-FM-DM";

		for (Period p : periods) {
			CurveContainer brAvgCCover = buildBrAvgCurveContainer(new File(srcOver), p, grid, spacing);
			CurveContainer brAvgCCunder = buildBrAvgCurveContainer(new File(srcUnder), p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet over = NSHMP_DataUtils.extractPE(brAvgCCover, gr, pe);
				GeoDataSet under = NSHMP_DataUtils.extractPE(brAvgCCunder, gr, pe);
				GeoDataSet xyzRatio = GeoDataSetMath.divide(over, under);
				GeoDataSet xyzDiff = GeoDataSetMath.subtract(over, under);
				
//				String id = "UC3-13-08";
//				String id = "UC3-newGrid";
				String id = "UC3-13-pt2v";
				
				String ratioDir = outDir + id + "-ratio-" + p.getLabel() + "-" + pe;
				makeRatioPlot(xyzRatio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.3, true, false);
				String diffDir = outDir +  id + "-diff-" + p.getLabel() + "-" + pe;
				double diffScale = p.equals(GM0P20) ? 0.4 : 0.2; // larger scale for 5Hz
				makeDiffPlot(xyzDiff, 0.1, grid.bounds(), diffDir, "GM diff", diffScale, true, false);
			}
		}
	}

	// build ratio maps of the 'final' (or close to it) NGAW2 ground motion
	// models relative to those that were used to date (3/10/2014) for the
	// NSHMP; these only used the single FM3.1 brAvg solution (720 branches)
	private static void buildGMPE14changes() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		String outDir = ROOT + "NSHMP13B" + S + "GMPE13-14Change" + S;
		String brSol = "UC33brAvg_FM31";
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		List<String> gmpes = Lists.newArrayList("AS", "BS","CB","CY","ID");
		
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {				
				for (String gmpe : gmpes) {
					String srcOver = SRC + "NSHMP13-GMPE14" + S + gmpe + S + brSol;
					GeoDataSet over = loadSingle(srcOver, pe, grid, p, spacing);

					String srcUnder = SRC + "NSHMP13-GMPE13" + S + gmpe + S + brSol;
					GeoDataSet under = loadSingle(srcUnder, pe, grid, p, spacing);
					
					GeoDataSet ratio = GeoDataSetMath.divide(over, under);
					String ratioDir = outDir + "ratio_" + p.getLabel() + S + gmpe;
//					makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "2% in 50 " + p.getLabel() + " ratio " + gmpe + " 14/13", true);
					makeRatioPlot(ratio, 0.1, grid.bounds(), ratioDir, "2% in 50 " + p.getLabel() + " ratio " + gmpe + " 14/13", true, 0.01, true, false);

//					GeoDataSet diff = GeoDataSetMath.subtract(over, under);
//					String diffDir = outDir + "diff_" + p.getLabel() + S + gmpe;
//					makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "2% in 50 " + p.getLabel() + " diff " + gmpe + " 14-13", true);
				}
			}
		}
	}

	// build ratio maps of the 'final' (or close to it) NGAW2 ground motion
	// models relative to those that were used to date (3/10/2014) for the
	// NSHMP; these only used the single FM3.1 brAvg solution (720 branches)
	private static void buildGMPE14_AS_bugfix() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		String outDir = ROOT + "NSHMP13B" + S + "GMPE13-14Change" + S;
		String brSol = "UC33brAvg_FM31";
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50);
		
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				String srcOver = SRC + "NSHMP13-GMPE14" + S + "AS-bugfix" + S + brSol;
				GeoDataSet over = loadSingle(srcOver, pe, grid, p, spacing);

				String srcUnder = SRC + "NSHMP13-GMPE14" + S + "AS" + S + brSol;
				GeoDataSet under = loadSingle(srcUnder, pe, grid, p, spacing);
					
				GeoDataSet ratio = GeoDataSetMath.divide(over, under);
				String ratioDir = outDir + "ratio_" + p.getLabel() + S + "AS-bugfix";
				makeRatioPlot(ratio, 0.1, grid.bounds(), ratioDir, "2% in 50 " + p.getLabel() + " ratio AS 14bf/14", true, 0.05, true, false);
			}
		}
	}

	
	// build ratio maps of the each individual 2013 GMPE to the weighted combo
	// these only used the single FM3.1 brAvg solution (720 branches)
	private static void buildGMEP13comparisons() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		String outDir = ROOT + "NSHMP13B" + S + "GMPE13ComparisonTMP" + S;
		String brSol = "UC33brAvg_FM31";
		
		List<Period> periods = Lists.newArrayList(GM0P00); //, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		List<String> gmpes = Lists.newArrayList("AS"); //,"BS","CB","CY","ID");
		
		
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// String srcUnder = SRC + "NSHMP13-GMPE13" + S + "ALL" + S + brSol;
				// GeoDataSet under = loadSingle(srcUnder, pe, grid, p, spacing);
				GeoDataSet under = combineGMPE13(p, pe);
				
				for (String gmpe : gmpes) {
					String srcOver = SRC + "NSHMP13-GMPE13" + S + gmpe + S + brSol;
					GeoDataSet over = loadSingle(srcOver, pe, grid, p, spacing);
					
					GeoDataSet ratio = GeoDataSetMath.divide(over, under);
					String ratioDir = outDir + "ratio_" + p.getLabel() + S + gmpe;
					makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "2% in 50 " + p.getLabel() + " ratio " + gmpe + "/ALL", true);

					GeoDataSet diff = GeoDataSetMath.subtract(over, under);
					String diffDir = outDir + "diff_" + p.getLabel() + S + gmpe;
					makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "2% in 50 " + p.getLabel() + " diff " + gmpe + "-ALL", true);
				}
			}
		}
	}
	
	private static GeoDataSet combineGMPE13(Period p, ProbOfExceed pe) {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		GriddedRegion gr = grid.grid(spacing);
		String brSol = "UC33brAvg_FM31";
		List<String> gmpes = Lists.newArrayList("AS","BS","CB","CY","ID");
		List<Double> weights = Lists.newArrayList(0.22, 0.22, 0.22, 0.22, 0.12);
		CurveContainer cc = null;
		int index = 0;
		for (String gmpe : gmpes) {
			String path = SRC + "NSHMP13-GMPE13" + S + gmpe + S + brSol;
			File cFile = new File(path + S + grid + S + p + S + "curves.csv");
			CurveContainer ccTmp = CurveContainer.create(cFile, grid, spacing);
			if (cc == null) {
				cc = ccTmp.scale(weights.get(index++));
				continue;
			}
			cc.add(ccTmp.scale(weights.get(index++)));
		}
		return NSHMP_DataUtils.extractPE(cc, gr, pe);
	}
	
	private static void buildIdrissComparison() {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		GriddedRegion gr = grid.grid(spacing);
		String outDir = ROOT + "NSHMP13B" + S + "IdrissTest" + S;
		String brSol = "UC33brAvg_FM31";
		List<Period> periods = Lists.newArrayList(GM0P00); //GM0P00, GM0P20, GM1P00);
		String pathOver = SRC + "NSHMP13-GMPE13" + S + "ID" + S + brSol;
		String pathUnder = SRC + "NSHMP13-GMPE13" + S + "~ID" + S + brSol;
		for (Period p : periods) {
			File cFileOver = new File(pathOver + S + grid + S + p + S + "curves.csv");
			File cFileUnder = new File(pathUnder + S + grid + S + p + S + "curves.csv");
			CurveContainer ccOver = CurveContainer.create(cFileOver, grid, spacing);
			CurveContainer ccUnder = CurveContainer.create(cFileUnder, grid, spacing);
			
			System.out.println(ccOver.getCurve(new Location(34.300,-117.50)));
			GeoDataSet geoOver = NSHMP_DataUtils.extractPE(ccOver, gr, PE2IN50);
			GeoDataSet geoUnder = NSHMP_DataUtils.extractPE(ccUnder, gr, PE2IN50);
			
			GeoDataSet ratio = GeoDataSetMath.divide(geoOver, geoUnder);
			String ratioDir = outDir + "ratio_" + p.getLabel();
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "2% in 50 " + p.getLabel() + " Idriss New/Old", true);

			GeoDataSet diff = GeoDataSetMath.subtract(geoOver, geoUnder);
			String diffDir = outDir + "diff_" + p.getLabel();
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "2% in 50 " + p.getLabel() + " Idriss New-Old", true);
			
		}
		
	}
	
	
	// combines 8 branch averaged solutions for faults (computed with
	// overlapping WUS fault removed, eg.g Carson) with 2 bravg solutions for
	// background (computed masking all those grid nodes outside CA)
	private static void buildGMEP13comparisonsBinaries() throws IOException {
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.1;
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<String> gmpes = Lists.newArrayList("ID"); //"AS","BS","CB","CY","ID","ALL");
		
		for (Period p : periods) {
			for (String gmpe : gmpes) {
			// single files
			File cFile = new File(SRC + "NSHMP13-GMPE13" + S + gmpe + S + 
				"UC33brAvg_FM31" + S + grid + S + p + S + "curves.csv");
			
			CurveContainer gmpeCC = CurveContainer.create(cFile, grid, spacing);

			File out = new File("tmp/NSHMP-CA/gmm-binaries_2-4-2014/" + gmpe + "-" + p.getLabel() + ".curves");
			Files.createParentDirs(out);
			
			BinaryCurves.writeUC3(gmpeCC, p, spacing, "UCERF3.3 GMM13:" + gmpe, out);
			}
		}
	}

	
	
	
//	// NGAW2 2013 comparsions
//	private static void buildUC3gmpe13() {
//		TestGrid grid = CA_RELM;
//		double spacing = 0.1;
//		String outDir = ROOT + "NewGMPE/";
//		
//		List<Period> periods = Lists.newArrayList(GM1P00); //, GM0P20, GM1P00);
//		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
//		GriddedRegion gr = grid.grid(spacing);
//
//		String srcOver = SRC + "UC33brAvg-FM-DM-13";
//		String srcUnder = SRC + "UC33brAvg-FM-DM-08";
//		
//		for (Period p : periods) {
//			CurveContainer brAvgCCover = buildBrAvgCurveContainer(new File(srcOver), p, grid, spacing);
//			CurveContainer brAvgCCunder = buildBrAvgCurveContainer(new File(srcUnder), p, grid, spacing);
//			for (ProbOfExceed pe : PEs) {
//				// brAvg data
//				GeoDataSet over = NSHMP_DataUtils.extractPE(brAvgCCover, gr, pe);
//				GeoDataSet under = NSHMP_DataUtils.extractPE(brAvgCCunder, gr, pe);
//				GeoDataSet xyzRatio = GeoDataSetMath.divide(over, under);
//				GeoDataSet xyzDiff = GeoDataSetMath.subtract(over, under);
//				
//				String id = "UC3-13-08";
//				String ratioDir = outDir + id + "-ratio-" + p + "-" + pe;
//				makeRatioPlot(xyzRatio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.1, true, false);
//				String diffDir = outDir +  id + "-diff-" + p + "-" + pe;
//				makeDiffPlot(xyzDiff, 0.1, grid.bounds(), diffDir, "GM diff", 0.2, true, false);
//			}
//		}
//	}
//	
	
	// ground motion maps and ratios for NSMP and UC3 using the 8brAvg solutions
	private static void buildBrAvgHazardMaps() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00);//, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// nshmp data
				GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca", pe, grid, p, spacing);
				String mapID = p.getLabel() + "-" + pe.name();
				nshmpXYZ.put(mapID,  xyzNSHMP);
				// nshmp hazard map
//				String dir = ROOT + "Hazard/nshmp-" + mapID;
//				makeHazardMap(xyzNSHMP, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 8brAvgMaps
		Map<String, GeoDataSet> brAvgXYZ = Maps.newHashMap();
//		String srcBase = SRC + "UC33brAvg-FM-DM-13" + S + "all";
		String srcBase = SRC + "UC33brAvg-FM-DM-08" + S + "all";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
//				String dir = ROOT + "Hazard/brAvg13-" + mapID;
				String dir = ROOT + "Hazard/brAvg08-" + mapID;
				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 / NSHMP
		for (String id : brAvgXYZ.keySet()) {
			GeoDataSet xyzNSHMP = nshmpXYZ.get(id);
			GeoDataSet xyzBrAvg = brAvgXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzBrAvg, xyzNSHMP);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzBrAvg, xyzNSHMP);
//			String ratioDir = ROOT + "Hazard/brAvg13_sup_nshmp-" + id;
//			String diffDir = ROOT + "Hazard/brAvg13_diff_nshmp-" + id;
			String ratioDir = ROOT + "Hazard/brAvg08_sup_nshmp-" + id;
			String diffDir = ROOT + "Hazard/brAvg08_diff_nshmp-" + id;
			makeRatioPlot(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.3, true, false);
			double diffScale = id.contains(GM0P20.getLabel()) ? 0.4 : 0.2; // larger scale for 5Hz
			makeDiffPlot(diff, 0.1, grid.bounds(), diffDir, "GM diff", diffScale, true, false);
		}
		
	}
	
	private static void build_NSHMP13B_TotalHazardMap() throws IOException {
		TestGrid grid = CA_NSHMP;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00, GM2P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// nshmp data
				GeoDataSet xyzNSHMP = loadSingle(SRC + "NSHMP13B" + S + 
					"mean_ucerf3_sol", pe, grid, p, spacing);
				String mapID = p.getLabel() + "-" + pe.name();
				nshmpXYZ.put(mapID,  xyzNSHMP);
				// nshmp hazard map
				String dir = ROOT + "NSHMP13B" + S + mapID;
				makeHazardMap(xyzNSHMP, spacing, p, pe, grid, dir);
			}
		}

	}
	
	// ground motion maps and ratios for UC3wGMPE13 / NSHM, both including CASC
	// and CAdeep; UC3 uses 8brAvg solutions
	private static void build_NSHMP_TotalHazardMaps() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				// nshmp data
				GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca_nw", pe, grid, p, spacing);
				String mapID = p.getLabel() + "-" + pe.name();
				nshmpXYZ.put(mapID,  xyzNSHMP);
				// nshmp hazard map
				String dir = ROOT + "NSHMP13" + S + "Total" + S + "nshmp_nw-" + mapID;
				makeHazardMap(xyzNSHMP, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 8brAvgMaps
		Map<String, GeoDataSet> brAvgXYZ = Maps.newHashMap();
		String srcBase = SRC + "NSHMP13-8-prob";
		GriddedRegion gr = grid.grid(spacing);
		TestGrid nshmpGrid = CA_NSHMP;
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, nshmpGrid, spacing);
			
			String deepFile1 = "/Users/pmpowers/projects/svn/NSHMP14/tmp/out/CAdeep_2014_Mmax8." + p.getLabel().toLowerCase();
			CurveContainer ccDeep1 = CurveContainer.create(new File(deepFile1).toURI().toURL());
			ccDeep1.scale(0.1);
			brAvgCC.addNSHMP(ccDeep1);
			
			String deepFile2 = "/Users/pmpowers/projects/svn/NSHMP14/tmp/out/CAdeep_2014_Mmax75." + p.getLabel().toLowerCase();
			CurveContainer ccDeep2 = CurveContainer.create(new File(deepFile2).toURI().toURL());
			ccDeep2.scale(0.9);
			brAvgCC.addNSHMP(ccDeep2);

			String deepFile3 = "/Users/pmpowers/projects/svn/NSHMP14/tmp/out/CAdeep_2014." + p.getLabel().toLowerCase();
			CurveContainer ccDeep3 = CurveContainer.create(new File(deepFile3).toURI().toURL());
			brAvgCC.addNSHMP(ccDeep3);

			String cascFile = "/Users/pmpowers/projects/svn/NSHMP14/tmp/out/sub.2014." + p.getLabel();
			CurveContainer ccCasc = CurveContainer.create(new File(cascFile).toURI().toURL());
			brAvgCC.addNSHMP(ccCasc);

			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name();
				brAvgXYZ.put(mapID,  xyzBrAvg);
				// braAvg hazard map
				String dir = ROOT + "NSHMP13" + S + "Total" + S + "NSHMP13nw-" + mapID;
				makeHazardMap(xyzBrAvg, spacing, p, pe, grid, dir);
			}
		}
		
		// UC3 / NSHMP
		for (String id : brAvgXYZ.keySet()) {
			GeoDataSet xyzNSHMP = nshmpXYZ.get(id);
			GeoDataSet xyzBrAvg = brAvgXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzBrAvg, xyzNSHMP);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzBrAvg, xyzNSHMP);
			String ratioDir = ROOT + "NSHMP13" + S + "Total" + S + "NSHMP13nw_sup_nshmp_nw-" + id;
			String diffDir = ROOT + "NSHMP13" + S + "Total" + S + "NSHMP13nw_diff_nshmp_nw-" + id;
			makeRatioPlotNSHMP(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", false);
			makeDiffPlotNSHMP(diff, 0.1, grid.bounds(), diffDir, "GM diff", false);
		}
		
	}
	
	// build ratio and diff maps of UC3-13 and nshmp for bg and flt sources
	private static void buildBgFltMaps() throws IOException {
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		
		List<Period> periods = Lists.newArrayList(GM0P20); //GM0P00, GM0P20, GM1P00);
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50); //, PE10IN50);
		List<String> srcs = Lists.newArrayList("bg", "flt");
		
		
		// NSHMP maps
		Map<String, GeoDataSet> nshmpXYZ = Maps.newHashMap();
		for (Period p : periods) {
			for (ProbOfExceed pe : PEs) {
				for (String src : srcs) {
					// nshmp data
					GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca_" + src, pe, grid, p, spacing);
					String mapID = p.getLabel() + "-" + pe.name() + "-" + src;
					nshmpXYZ.put(mapID,  xyzNSHMP);
				}
			}
		}
		
		Map<String, GeoDataSet> brAvgXYZ = Maps.newHashMap();
		GriddedRegion gr = grid.grid(spacing);

		// UC3 faults 8brAvgMaps FM-DM
		String fltBase = SRC + "UC33brAvg-FM-DM-13" + S + "flt";
		for (Period p : periods) {
			File srcDir = new File(fltBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name() + "-flt";
				brAvgXYZ.put(mapID,  xyzBrAvg);
			}
		}
		
		// UC3 bg 2brAvgMaps FM
		String bgBase = SRC + "UC33brAvg-FM-DM-13" + S + "bg";
		for (Period p : periods) {
			File srcDir = new File(bgBase);
			CurveContainer brAvgCC = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
			for (ProbOfExceed pe : PEs) {
				// brAvg data
				GeoDataSet xyzBrAvg = NSHMP_DataUtils.extractPE(brAvgCC, gr, pe);
				String mapID = p.getLabel() + "-" + pe.name() + "-bg";
				brAvgXYZ.put(mapID,  xyzBrAvg);
			}
		}

		// UC3 / NSHMP
		for (String id : brAvgXYZ.keySet()) {
			GeoDataSet xyzNSHMP = nshmpXYZ.get(id);
			GeoDataSet xyzBrAvg = brAvgXYZ.get(id);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzBrAvg, xyzNSHMP);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzBrAvg, xyzNSHMP);
			String ratioDir = ROOT + "Hazard-BgFlt/brAvg13_sup_nshmp_" + id;
			String diffDir = ROOT + "Hazard-BgFlt/brAvg13_diff_nshmp_" + id;
			makeRatioPlot(ratio, 0.1, grid.bounds(), ratioDir, "GM ratio", true, 0.3, true, false);
			double diffScale = id.contains(GM0P20.getLabel()) ? 0.4 : 0.2; // larger scale for 5Hz
			makeDiffPlot(diff, 0.1, grid.bounds(), diffDir, "GM diff", diffScale, true, false);
		}
		
	}
	

	// builds the maps that compare FM-DM and FM-DM-MS brAVg maps to full tree
	// these were done prior to the grid source update
	private static void buildBrAvgMaps() throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		double spacing = 0.1;
		String dlDir = ROOT + "BranchAvgComparison/";
		GriddedRegion gr = grid.grid(spacing);
		File srcDir = null;
		
		// build brAvg curve containers
		srcDir = new File(SRC + "UC33brAvg-FM");
		CurveContainer cc_FM = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
		GeoDataSet xyz_FM = NSHMP_DataUtils.extractPE(cc_FM, gr, pe);
		makeHazardMap(xyz_FM, spacing, p, pe, grid, dlDir + "map-FM/");

//		srcDir = new File(SRC + "UC33brAvg-FM-DM");
//		CurveContainer cc_FM_DM = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
//		GeoDataSet xyz_FM_DM = NSHMP_DataUtils.extractPE(cc_FM_DM, gr, pe);
//		makeHazardMap(xyz_FM_DM, spacing, p, pe, grid, dlDir + "map-FM-DM/");
//		
//		srcDir = new File(SRC + "UC33brAvg-FM-DM-MS");
//		CurveContainer cc_FM_DM_MS = buildBrAvgCurveContainer(srcDir, p, grid, spacing);
//		GeoDataSet xyz_FM_DM_MS = NSHMP_DataUtils.extractPE(cc_FM_DM_MS, gr, pe);
//		makeHazardMap(xyz_FM_DM_MS, spacing, p, pe, grid, dlDir + "map-FM-DM-MS/");
		
		
		// load full tree denominator
		File brUnderFile = new File(dlDir + "../LogicTreeRatios/branchsets/all.txt");
		GeoDataSet xyzUnder = loadMulti(SRC + "UC33/", brUnderFile, pe, grid, p, suffix);
		
		String dlPath = null;
		GeoDataSet xyzOver = null;
		GeoDataSet xyzRatio = null;
		
		// FM map
		xyzOver = NSHMP_DataUtils.extractPE(cc_FM, gr, pe);
		
		dlPath = dlDir + "FM-0.1/";
		// ratio data sets are corrupted by gmt so recreate when using for multiple maps
		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.1, true, false);
		dlPath = dlDir + "FM-0.05/";
		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.05, true, false);
		dlPath = dlDir + "FM-0.01/";
		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.01, true, false);

//		// FM-DM map
//		xyzOver = NSHMP_DataUtils.extractPE(cc_FM_DM, gr, pe);
//		
//		dlPath = dlDir + "FM-DM-0.1/";
//		// ratio data sets are corrupted by gmt so recreate when using for multiple maps
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.1, true, false);
//		dlPath = dlDir + "FM-DM-0.05/";
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.05, true, false);
//		dlPath = dlDir + "FM-DM-0.01/";
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.01, true, false);
//
//		// FM-DM-MS map
//		xyzOver = NSHMP_DataUtils.extractPE(cc_FM_DM_MS, gr, pe);
//
//		dlPath = dlDir + "FM-DM-MS-0.1/";
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.1, true, false);
//		dlPath = dlDir + "FM-DM-MS-0.05/";
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.05, true, false);
//		dlPath = dlDir + "FM-DM-MS-0.01/";
//		xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlPath, "brAvg/fullTree", true, 0.01, true, false);
		
	}
	
	private static void makeHazardMap(GeoDataSet xyz, double spacing, Period p, ProbOfExceed pe, TestGrid grid, String outDir) {
		double[] minmax = NSHMP_PlotUtils.getRange(p);
		GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
		String label = pe + " " + p.getLabel() + " (g)";
		makeMapPlot(xyz, spacing, grid.bounds(), outDir, label, minmax[0], minmax[1], cpt, true, false);
	}
	
	private static CurveContainer buildBrAvgCurveContainer(File srcDir, Period p,
			TestGrid grid, double spacing) {
		
		Map<String, Double> fmWts = Maps.newHashMap();
		fmWts.put("FM31", 0.5);
		fmWts.put("FM32", 0.5);
		Map<String, Double> dmWts = Maps.newHashMap();
		dmWts.put("ABM", 0.1);
		dmWts.put("GEOL", 0.3);
		dmWts.put("NEOK", 0.3);
		dmWts.put("ZENGBB", 0.3);
		Map<String, Double> msWts = Maps.newHashMap();
		msWts.put("ELLB", 0.2);
		msWts.put("ELLBSL", 0.2);
		msWts.put("HB08", 0.2);
		msWts.put("SH09M", 0.2);
		msWts.put("SHCSD", 0.2);

		CurveContainer brAvgCC = null;

		for (File brAvgDir : srcDir.listFiles()) {
			if (!brAvgDir.isDirectory()) continue;
			String brAvgName = brAvgDir.getName();
			System.out.println(brAvgDir.getName());
			String brAvgPath = brAvgDir.getPath() + S + grid + S + p + S + "curves.csv";
			File brAvgFile = new File(brAvgPath);
			System.out.println("path: " + brAvgPath);
			
			// brAvg weight
			Iterator<String> ids = SPLIT.split(brAvgName).iterator();
			ids.next(); // skip first part
			Double fm = fmWts.get(ids.next());
			Double dm = ids.hasNext() ? dmWts.get(ids.next()) : 1.0;
			Double ms = ids.hasNext() ? msWts.get(ids.next()) : 1.0;
			double wt = fm * dm * ms;
			System.out.println("wt: " + wt);
			
			// create and weight curve container
			CurveContainer cc = CurveContainer.create(brAvgFile, grid, spacing);
			cc.scale(wt);
			if (brAvgCC == null) {
				brAvgCC = cc;
			} else {
				brAvgCC.add(cc);
			}
		}
		return brAvgCC;
	}
	

	
	// UCERF3.3 node ratio maps
	private static void makeRatioMaps(String srcDir, String outDir,
			List<String> brOverList, String brUnder) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		
		File brUnderFile = new File(outDir + "branchsets", brUnder + ".txt");
		System.out.println(brUnderFile);
		GeoDataSet xyzUnder = loadMulti(srcDir, brUnderFile, pe, grid, p, suffix);

		for (String brOver : brOverList) {
			File brOverFile = new File(outDir + "branchsets", brOver + ".txt");
			GeoDataSet xyzOver = loadMulti(srcDir, brOverFile, pe, grid, p, suffix);
		
			GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);

			String dlDir = outDir + brOver + "_sup_" + brUnder;
			makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "GM ratio", true, 0.1, true, false);
		}
	}
	
	// UCERF3.3 NSHMP ratio maps
	private static void makeNSHMPratioMap(String srcDir, String outDir,
			String branches) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		
		GeoDataSet xyzNSHMP = loadSingle(SRC + "nshmp_ca/", pe, grid, p, 0.1);
		File branchFile = new File(outDir + "branchsets", branches + ".txt");
		GeoDataSet xyzUC32 = loadMulti(srcDir, branchFile, pe, grid, p, "");
		GeoDataSet xyz = GeoDataSetMath.divide(xyzUC32, xyzNSHMP);
		String dlDir = outDir + branches + "_sup_NSHMP";
		makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "GM ratio", true, 0.3, false, false);
	}

	
	public static void makeRatioPlotNSHMP(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, 0.5);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 2.0);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = GMT_CPT_Files.NSHMP_RATIO;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(false);
			if (showFaults) {
				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
//				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
//				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public static void makeDiffPlotNSHMP(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, -0.5);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 0.5);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = GMT_CPT_Files.NSHMP_DIFF;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(true);
			if (showFaults) {
				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
//				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
//				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public static void makeRatioPlot(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title, boolean log,
			double logScale, boolean smooth, boolean showFaults) {
		double scale = log ? logScale : 0.2;
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, log ? -scale : 1-scale);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, log ? scale : 1+scale);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = log ? GMT_CPT_Files.UCERF3_HAZ_RATIO_P3 : GMT_CPT_Files.GMT_POLAR;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, log ? true : false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(smooth);
			if (showFaults) {
//				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static void makeDiffPlot(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title,
			double scale, boolean smooth, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, -scale);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, scale);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = GMT_CPT_Files.UCERF3_HAZ_RATIO_P3;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(smooth);
			if (showFaults) {
//				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}


	private static void addFaultTraces(FaultModels fm, GMT_Map map, Color c) {
		List<FaultSectionPrefData> faults = fm.fetchFaultSections();
		for (FaultSectionPrefData fspd : faults) {
			PSXYPolygon poly = new PSXYPolygon(fspd.getFaultTrace());
			poly.setPenColor(c);
			poly.setPenWidth(2);
			map.addPolys(poly);
		}
	}

	private static GeoDataSet loadSingle(String dir, ProbOfExceed pe,
			TestGrid grid, Period p, double spacing) {
		File curves = new File(dir + S + grid + S + p + S + "curves.csv");
		CurveContainer cc = CurveContainer.create(curves, grid, spacing);
		GriddedRegion gr = grid.grid(spacing);
		return NSHMP_DataUtils.extractPE(cc, gr, pe);
	}
	
	private static GeoDataSet loadSingle(String dir, double prob,
			TestGrid grid, Period p, double spacing) {
		File curves = new File(dir + S + grid + S + p + S + "curves.csv");
		CurveContainer cc = CurveContainer.create(curves, grid, spacing);
		GriddedRegion gr = grid.grid(spacing);
		return NSHMP_DataUtils.extractPE(cc, gr, prob);
	}

	private static GeoDataSet loadMulti(String srcDir, File branchListFile,
			ProbOfExceed pe, TestGrid grid, Period p, String suffix)
			throws IOException {

		List<String> branchNames = Files.readLines(branchListFile, US_ASCII);
		System.out.println("Loading: " + branchListFile.getName());
		
		// create wt list
		List<Double> wtList = Lists.newArrayList();
		for (String brName : branchNames) {
			LogicTreeBranch branch = LogicTreeBranch.fromFileName(brName);
			wtList.add(branch.getAprioriBranchWt());
		}
		DataUtils.asWeights(wtList);
		
		String cPath = grid + S + p + S + "curves.csv";
		GriddedRegion gr = grid.grid(0.1);
		CurveContainer mapcc = null;
		
		int idx = 0;
		for (String brName : branchNames) {
			if (idx % 100 == 0) System.out.print(idx + " ");
			String brID = brName + suffix;
			String brPath = srcDir + S + brID + S + cPath;
			File brFile = new File(brPath);
			CurveContainer cc = CurveContainer.create(brFile, grid, 0.1);
			cc.scale(wtList.get(idx++));

			if (mapcc == null) {
				mapcc = cc;
			} else {
				mapcc.add(cc);
			}
		}
		return NSHMP_DataUtils.extractPE(mapcc, gr, pe);
	}
	
	
	
	
	
	// UCERF3.3 prelim branchAvg over comparable UC32 branch avg.
	private static void makeBrAvgRatioMapUC33() throws IOException {

		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		double spacing = 0.1;

		String over = SRC + "UC33brAvg5x_fm32";
		String under = SRC + "UC32brAvg5x_fm32";
		String out = ROOT + "UC-33-32-brAvg-fm32";

		GeoDataSet xyzOver = loadSingle(over, pe, grid, p, spacing);
		GeoDataSet xyzUnder = loadSingle(under, pe, grid, p, spacing);
		GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);

		makeRatioPlot(xyz, 0.1, grid.bounds(), out, "GM ratio",
			true, 0.1, true, true);
	}
	
	// UCERF3.3 prelim ground motion maps
	private static void makePrelimBrAvgHazardMaps() {
		
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		double spacing = 0.1;
		
		String src33 = SRC + "UC33brAvg_prelim";
		String src32 = SRC + "UC32brAvg_cf33";

		GeoDataSet xyz33 = loadSingle(src33, pe, grid, p, spacing);
		GeoDataSet xyz32 = loadSingle(src32, pe, grid, p, spacing);

		String out33 = ROOT + "UC33brAvgMap";
		String out32 = ROOT + "UC32brAvgMap";

		// map
		double[] minmax = NSHMP_PlotUtils.getRange(p);
		GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
		String label = pe + " " + p.getLabel() + " (g)";
		
		makeMapPlot(xyz33, spacing, grid.bounds(), out33, label,
			minmax[0], minmax[1], cpt, true, true);
		makeMapPlot(xyz32, spacing, grid.bounds(), out32, label,
			minmax[0], minmax[1], cpt, true, true);
	}
	
	// UCERF3.3 tests of convergence; investigates the distribution of mean
	// hazard/ground motion across multiple inversion runs; sample = 10 runs
	// and we have maps for each
	private static void doInversionRunAnalysis(String dir, Period period, ProbOfExceed pe) throws IOException {
//		String dir = "UC33-10runs-PGA";
		String path1 = "UC33-brAvg-fm31-run";
		String path2 = "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_WITH_IND_RUNS_FM3_1_run";
		String path3 = "_MEAN_BRANCH_AVG_SOL";
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		int runCount = 10;
		GriddedRegion gr = grid.grid(spacing);
		
//		String under = SRC + "UC33brAvg5x_fm31";
//		GeoDataSet xyzUnder = loadSingle(under, pe, grid, period, spacing);
		
		// init stat arrays
		List<double[]> solMeans = Lists.newArrayList();
		for (int i=0; i<gr.getNodeCount(); i++) {
			solMeans.add(new double[10]);
		}
		
		// aggregate stat data and output ratio maps
		for (int i=0; i<runCount; i++) {
			String solDir = SRC + dir + S + path1 + i + S + path2 + i + path3;
			System.out.println("Processing sol: " + i);
			GeoDataSet xyzSol = loadSingle(solDir, pe, grid, period, spacing);
			for (int j=0; j<gr.getNodeCount(); j++) {
				solMeans.get(j)[i] = xyzSol.get(j);
			}
			
//			// maps - NOTE denominator is only first 5 inv runs, not 10
//			GeoDataSet xyz = GeoDataSetMath.divide(xyzSol, xyzUnder);
//			String out = ROOT + "UC33invTest-" + period + "-" + pe + "-" + i + S;
//			makeRatioPlot(xyz, 0.1, grid.bounds(), out, "GM ratio",
//				true, true, false);
		}
		
		LocationList locs = gr.getNodeList();
		
		// compute stats
		List<Double> means = Lists.newArrayList();
		List<Double> stds = Lists.newArrayList();
		
		Mean meanCalc = new Mean();
		StandardDeviation stdCalc = new StandardDeviation();
		for (int i=0; i<gr.getNodeCount(); i++) {
			double[] runData = solMeans.get(i);
			double mean = meanCalc.evaluate(runData);
			double std = stdCalc.evaluate(runData, mean);
			means.add(mean);
			stds.add(std);
		}
		
		// set mean and meanOverStd data
//		GriddedGeoDataSet xyzMean = new GriddedGeoDataSet(gr, true);
//		for (int i=0; i<gr.getNodeCount(); i++) {
//			xyzMean.set(i, means.get(i));
//		}
		
//		// compare mean(mean of 10 runs) to mean (5 runs - combined sol) 
//		GeoDataSet xyz = GeoDataSetMath.divide(xyzMean, xyzUnder);
//		String out = ROOT + "UC33invTest-10xMean-5xMean-" + period + "-" + pe;
//		makeRatioPlot(xyz, 0.1, grid.bounds(), out, "GM ratio",
//			true, true, false);

		// mean over std data
		String outPath = ROOT + "MapConvTests/stats/";
		String periodHead = (period.equals(GM0P00) ? "PGA" : "3sec");
		String peHead = (pe.equals(PE2IN50) ? "2in50" : (pe.equals(PE10IN50) ? "10in50" : "1in100"));
		String label = periodHead + "-" + peHead;
		File datFile = new File(outPath, "invStats-" + label + ".csv");
		List<String> headDat = Lists.newArrayList("lat", "lon", label+"-mean", label+"-std");
		String header = JOIN.join(headDat) + "\n";
		Files.write(header, datFile, US_ASCII);
		for (int i=0; i<locs.size(); i++) {
			Location loc = locs.get(i);
			List<String> locDat = Lists.newArrayList(
				String.format(format, loc.getLatitude()),
				String.format(format, loc.getLongitude()),
				Double.toString(means.get(i)),
				Double.toString(stds.get(i)));
			Files.append(JOIN.join(locDat) + "\n", datFile, Charsets.US_ASCII);
		}
		
	}


	
	public static void makeMapPlot(GeoDataSet xyz, double spacing, double[] bounds,
			String dlDir, String title, double scaleMin, double scaleMax,
			GMT_CPT_Files cpt, boolean smooth, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, scaleMin);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, scaleMax);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(smooth);
			if (showFaults) {
//				addFaultTraces(FaultModels.FM2_1, map, Color.BLUE);
				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	
	// output of 2in50 lat,lon,val data for K. Rukstales
	private static void createCAxyzFiles() throws IOException {
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		ProbOfExceed pe = PE10IN50;
		TestGrid grid = CA_RELM;
		double spacing = 0.1;
		String outDir = "tmp/forKR/";
		
		Region extRegion = new Region(
			new Location(50.0, -125.0),
			new Location(24.6, -100.0));
		
		for (Period p : periods) {
			// load new map data UC3 + NGA13
			File caDir = new File(SRC + "UC33brAvg-FM-DM-13" + S + "all");
			CurveContainer ccCA = buildBrAvgCurveContainer(caDir, p, grid, spacing);
			
			String extDir = "/Users/pmpowers/projects/NSHMP/tmp/out/";
			String pString = (p == GM0P00) ? "pga" : (p == GM0P20) ? "5hz" : "1hz";
			String cascPath = extDir + "sub.2014." + pString;
			String deepPath = extDir + "CAdeep_2014." + pString;
			URL cascUrl = new File(cascPath).toURI().toURL();
			URL deepUrl = new File(deepPath).toURI().toURL();
			CurveContainer ccCasc = CurveContainer.create(cascUrl);
			CurveContainer ccDeep = CurveContainer.create(deepUrl);
			
			File outFile = new File(outDir, "CA-" + p.getLabel() + "-" + pe + ".csv");
			Files.createParentDirs(outFile);
			Files.write("", outFile, US_ASCII); // clears preexisting file
			for (Location loc : ccCA) {
				if (!extRegion.contains(loc)) continue;
				DiscretizedFunc ca = ccCA.getCurve(loc);
				DiscretizedFunc casc = ccCasc.getCurve(loc);
				DiscretizedFunc deep = ccDeep.getCurve(loc);
				
				addCurves(ca, casc);
				addCurves(ca, deep);
				
				double val = ca.getFirstInterpolatedX_inLogXLogYDomain(pe.annualRate());
				
				List<String> locDat = Lists.newArrayList(
					String.format(format, loc.getLatitude()),
					String.format(format, loc.getLongitude()),
					Double.toString(val));
				String line = JOIN.join(locDat) + LF;
				Files.append(line, outFile, US_ASCII); 
			}
		}
	}
	
	// adds c2 to c1 in place - for the time being this shouldn't have problems
	// as we loop over the base OPENSHA calculated curve. PGA curves computed
	// by SH have an additional value and the last value (2.13) is instead 2.2.
	// We're just putting the 2.2 value in the 2.13 bin for now as this
	// shouldn't affect 2in50.
	// By looping sha, the last value, if present, will be ignored. Although
	// hazpoint generates curves that are short when they contain trailing
	// zero values, a curve container keeps the zeros (or very low values)
	private static void addCurves(DiscretizedFunc c1, DiscretizedFunc c2) {
		for (int i=0; i<c1.size(); i++) {
			c1.set(i, c1.getY(i) + c2.getY(i));
		}
	}
	

	private static void riskRatioMapsCA() throws IOException {
		
		// GMM08-combined-1Hz.2475
		String srcDir = SRC + "NSHMP13-risk/";
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		
//		// 2%in50 maps
//		for (Period p : periods) {
//			String overPath = srcDir + "GMM13-combined-"+p.getLabel()+".2475.csv";
//			String underPath = srcDir + "GMM08-combined-"+p.getLabel()+".2475.csv";
//			GeoDataSet xyzOver = loadLatLonValCA(new File(overPath), false);
//			GeoDataSet xyzUnder = loadLatLonValCA(new File(underPath), false);
//			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
//			String ratioDir = ROOT + "NSHMP13/Risk-CA/GMM13_sup_GMM08-2in50-" + p.getLabel();
//			String diffDir = ROOT + "NSHMP13/Risk-CA/GMM13_diff_GMM08-2in50-" + p.getLabel();
//			makeRatioPlotNSHMP(ratio, 0.1, CA_RELM.bounds(), ratioDir, "GM ratio " + p.getLabel(), false);
//			makeDiffPlotNSHMP(diff, 0.1, CA_RELM.bounds(), diffDir, "GM diff " + p.getLabel(), false);
//		}
//		
//		// RTGM maps
//		for (Period p : periods) {
//			String overPath = srcDir + "GMM13-combined-"+p.getLabel()+".C_Rs.csv";
//			String underPath = srcDir + "GMM08-combined-"+p.getLabel()+".C_Rs.csv";
//			GeoDataSet xyzOver = loadLatLonValCA(new File(overPath), true);
//			GeoDataSet xyzUnder = loadLatLonValCA(new File(underPath), true);
//			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
//			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
//			String ratioDir = ROOT + "NSHMP13/Risk-CA/GMM13_sup_GMM08-RTGM-" + p.getLabel();
//			String diffDir = ROOT + "NSHMP13/Risk-CA/GMM13_diff_GMM08-RTGM-" + p.getLabel();
//			makeRatioPlotNSHMP(ratio, 0.1, CA_RELM.bounds(), ratioDir, "GM ratio " + p.getLabel(), false);
//			makeDiffPlotNSHMP(diff, 0.1, CA_RELM.bounds(), diffDir, "GM diff " + p.getLabel(), false);
//		}

		// Cr maps
		for (Period p : periods) {
			String overPath = srcDir + "GMM13-combined-"+p.getLabel()+".C_Rs.csv";
			String underPath = srcDir + "GMM08-combined-"+p.getLabel()+".C_Rs.csv";
			GeoDataSet xyzOver = loadLatLonValCA(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValCA(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-CA/SRC13_GMM13_sup_GMM08-Cr-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-CA/SRC13-GMM13_diff_GMM08-Cr-" + p.getLabel();
			makeCrRatioPlot(ratio, 0.1, CA_RELM.bounds(), ratioDir, "Cr ratio " + p.getLabel(), false);
			makeCrDiffPlot(diff, 0.1, CA_RELM.bounds(), diffDir, "Cr diff " + p.getLabel(), false);
		}
	}
	
	private static void riskRatioMapsCA2() throws IOException {
		// total model for CA; zoom of US maps (below)
		String srcDir = SRC + "NSHMP13-risk/";
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.05;
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		// Uniform hazard ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.UH.csv";
			String underPath = srcDir + "UShazard.20081229."+p.getLabel()+".UH.csv";
			GeoDataSet xyzOver = loadLatLonValUSforCA(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUSforCA(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-CA/CA-UH-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-CA/CA-UH-diff-" + p.getLabel();
			makeRatioPlotNSHMP(ratio, spacing, grid.bounds(), ratioDir, "GM ratio " + p.getLabel(), false);
			makeDiffPlotNSHMP(diff, spacing, grid.bounds(), diffDir, "GM diff " + p.getLabel(), false);
		}
		// Cr ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.C_Rs.csv";
			String underPath = srcDir + "USriskCoefficients.20081229."+p.getLabel()+".C_Rs.csv";
			GeoDataSet xyzOver = loadLatLonValUSforCA(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUSforCA(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-CA/CA-Cr-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-CA/CA-Cr-diff-" + p.getLabel();
			makeCrRatioPlot(ratio, spacing, grid.bounds(), ratioDir, "Cr ratio " + p.getLabel(), false);
			makeCrDiffPlot(diff, spacing, grid.bounds(), diffDir, "Cr diff " + p.getLabel(), false);
		}
	}

	private static void riskRatioMapsCA2_CrBeta0p8() throws IOException {
		// total model for CA; zoom of US maps (below)
		String srcDir = SRC + "NSHMP13-risk/";
		TestGrid grid = TestGrid.CA_RELM;
		double spacing = 0.05;
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		// Cr ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.beta0p8.C_Rs.csv";
			String underPath = srcDir + "USriskCoefficients.20081229."+p.getLabel()+".C_Rs.csv";
			GeoDataSet xyzOver = loadLatLonValUSforCA(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUSforCA(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-CA/CA-Cr-b0p8-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-CA/CA-Cr-b0p8-diff-" + p.getLabel();
			makeCrRatioPlot(ratio, spacing, grid.bounds(), ratioDir, "Cr b0.8 ratio " + p.getLabel(), false);
			makeCrDiffPlot(diff, spacing, grid.bounds(), diffDir, "Cr b0.8 diff " + p.getLabel(), false);
		}
	}

	private static void riskRatioMapsUS() throws IOException {
		
		String srcDir = SRC + "NSHMP13-risk/";
		TestGrid grid = TestGrid.NATIONAL;
		double spacing = 0.05;
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		
		// Uniform hazard ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.UH.csv";
			String underPath = srcDir + "UShazard.20081229."+p.getLabel()+".UH.csv";
			GeoDataSet xyzOver = loadLatLonValUS(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUS(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-US/US-UH-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-US/US-UH-diff-" + p.getLabel();
			makeRatioPlotNSHMP(ratio, spacing, grid.bounds(), ratioDir, "GM ratio " + p.getLabel(), false);
			makeDiffPlotNSHMP(diff, spacing, grid.bounds(), diffDir, "GM diff " + p.getLabel(), false);
		}
		// Cr ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.C_Rs.csv";
			String underPath = srcDir + "USriskCoefficients.20081229."+p.getLabel()+".C_Rs.csv";
			GeoDataSet xyzOver = loadLatLonValUS(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUS(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-US/US-Cr-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-US/US-Cr-diff-" + p.getLabel();
			makeCrRatioPlot(ratio, spacing, grid.bounds(), ratioDir, "Cr ratio " + p.getLabel(), false);
			makeCrDiffPlot(diff, spacing, grid.bounds(), diffDir, "Cr diff " + p.getLabel(), false);
		}

	}
	
	private static void riskRatioMapsUS_CrBeta0p8() throws IOException {
		String srcDir = SRC + "NSHMP13-risk/";
		TestGrid grid = TestGrid.NATIONAL;
		double spacing = 0.05;
		List<Period> periods = Lists.newArrayList(GM0P20, GM1P00);
		// Cr ratio-diff maps
		for (Period p : periods) {
			String overPath = srcDir + "us_epimerge."+p.getLabel()+".v131018.beta0p8.C_Rs.csv";
			String underPath = srcDir + "USriskCoefficients.20081229."+p.getLabel()+".C_Rs.csv";
			GeoDataSet xyzOver = loadLatLonValUS(new File(overPath), false);
			GeoDataSet xyzUnder = loadLatLonValUS(new File(underPath), false);
			GeoDataSet ratio = GeoDataSetMath.divide(xyzOver, xyzUnder);
			GeoDataSet diff = GeoDataSetMath.subtract(xyzOver, xyzUnder);
			String ratioDir = ROOT + "NSHMP13/Risk-US/US-Cr-b0p8-ratio-" + p.getLabel();
			String diffDir = ROOT + "NSHMP13/Risk-US/US-Cr-b0p8-diff-" + p.getLabel();
			makeCrRatioPlot(ratio, spacing, grid.bounds(), ratioDir, "Cr b0.8 ratio " + p.getLabel(), false);
			makeCrDiffPlot(diff, spacing, grid.bounds(), diffDir, "Cr b0.8 diff " + p.getLabel(), false);
		}

	}


	private static GeoDataSet loadLatLonValUS(File f, boolean valInCol4) throws IOException {
		GriddedRegion grUS = TestGrid.NATIONAL_POLY.grid(0.05);
		GeoDataSet xyz = new GriddedGeoDataSet(grUS, true);
		for (String line : Files.readLines(f, Charsets.US_ASCII)) {
			List<String> parts = Lists.newArrayList(SPLIT_COMMA.split(line));
			double lat = Double.valueOf(parts.get(0));
			double lon = Double.valueOf(parts.get(1));
			double val = Double.valueOf(parts.get(valInCol4 ? 3 : 2));
			Location loc = new Location(lat, lon);
			if (grUS.contains(loc)) {
				xyz.set(loc, val);
			}
		}
		return xyz;
	}

	private static GeoDataSet loadLatLonValUSforCA(File f, boolean valInCol4) throws IOException {
		GriddedRegion grCA = new GriddedRegion(getCA_mask(), 0.05, GriddedRegion.ANCHOR_0_0);
		Region rUS = TestGrid.NATIONAL_POLY.grid(1.0);
		GeoDataSet xyz = new GriddedGeoDataSet(grCA, true);
		for (String line : Files.readLines(f, Charsets.US_ASCII)) {
			List<String> parts = Lists.newArrayList(SPLIT_COMMA.split(line));
			double lat = Double.valueOf(parts.get(0));
			double lon = Double.valueOf(parts.get(1));
			double val = Double.valueOf(parts.get(valInCol4 ? 3 : 2));
			Location loc = new Location(lat, lon);
			if (grCA.contains(loc) && rUS.contains(loc)) {
				xyz.set(loc, val);
			}
		}
		return xyz;
	}

	private static GeoDataSet loadLatLonValCA(File f, boolean valInCol4) throws IOException {
		GriddedRegion grCA = new GriddedRegion(getCA_mask(), 0.1, GriddedRegion.ANCHOR_0_0);
		GeoDataSet xyz = new GriddedGeoDataSet(grCA, true);
		for (String line : Files.readLines(f, Charsets.US_ASCII)) {
			List<String> parts = Lists.newArrayList(SPLIT_COMMA.split(line));
			double lat = Double.valueOf(parts.get(0));
			double lon = Double.valueOf(parts.get(1));
			double val = Double.valueOf(parts.get(valInCol4 ? 3 : 2));
			Location loc = new Location(lat, lon);
			if (grCA.contains(loc)) {
				xyz.set(loc, val);
			}
		}
		return xyz;
	}
	
	private static Region getCA_mask() {
			LocationList locs = new LocationList();
			locs.add(new Location(39.000, -119.999));
			locs.add(new Location(35.000, -114.635));
			locs.add(new Location(34.848, -114.616));
			locs.add(new Location(34.719, -114.482));
			locs.add(new Location(34.464, -114.371));
			locs.add(new Location(34.285, -114.122));
			locs.add(new Location(34.097, -114.413));
			locs.add(new Location(33.934, -114.519));
			locs.add(new Location(33.616, -114.511));
			locs.add(new Location(33.426, -114.636));
			locs.add(new Location(33.401, -114.710));
			locs.add(new Location(33.055, -114.676));
			locs.add(new Location(33.020, -114.501));
			locs.add(new Location(32.861, -114.455));
			locs.add(new Location(32.741, -114.575));
			locs.add(new Location(32.718, -114.719));
			locs.add(new Location(32.151, -120.861));
			locs.add(new Location(39.000, -126.000));
			locs.add(new Location(42.001, -126.000));
			locs.add(new Location(42.001, -119.999));
			locs.add(new Location(39.000, -119.999));
			return new Region(locs, null);
	}

	public static void makeCrRatioPlot(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, 0.8);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 1.2);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = GMT_CPT_Files.NSHMP_DIFF;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(true);
			if (showFaults) {
//				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static void makeCrDiffPlot(GeoDataSet xyz, double spacing,
			double[] bounds, String dlDir, String title, boolean showFaults) {
		GMT_MapGenerator mapGen = NSHMP_PlotUtils.create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, -0.2);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 0.2);
		mapGen.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		GMT_CPT_Files cpt = GMT_CPT_Files.NSHMP_DIFF;
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			map.setRescaleCPT(true);
			if (showFaults) {
//				addFaultTraces(FaultModels.FM2_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_1, map, Color.BLACK);
				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			savePDF(dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}


	// combines curves and outputs curve container as original curves.csv file
	private static void combineCurves() throws FileNotFoundException, IOException {
		List<Period> periods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		TestGrid grid = CA_RELM;
		double spacing = 0.1;

//		for (Period p : periods) {
//			File dir = new File(SRC + "UC33brAvg-FM-DM-13" + S + "all");
//			CurveContainer cc = buildBrAvgCurveContainer(dir, p, grid, spacing);
//			
//			File out = new File(SRC + "UC33brAvg-FM-DM-13" + S + "combined-" + 
//					p.getLabel() + ".csv");
//			BufferedWriter writer = Files.newWriter(out, Charsets.US_ASCII);
//			HazardResultWriterLocal.writeCurveHeader(writer, p);
//			for (Location loc : cc) {
//				String resultStr = formatCurve(loc, cc.getValues(loc));
//				writer.write(resultStr);
//				writer.newLine();
//			}
//		}
		
		// UC3 8brAvgMaps
		TestGrid nshmpGrid = TestGrid.CA_NSHMP;
		String srcBase = SRC + "NSHMP13-8-prob";
		GriddedRegion gr = grid.grid(spacing);
		for (Period p : periods) {
			File srcDir = new File(srcBase);
			CurveContainer cc = buildBrAvgCurveContainer(srcDir, p, nshmpGrid, spacing);
			File out = new File(SRC + "NSHMP13-8-prob" + S + "combined-" + 
					p.getLabel() + ".csv");
			BufferedWriter writer = Files.newWriter(out, Charsets.US_ASCII);
			HazardResultWriterLocal.writeCurveHeader(writer, p);
			for (Location loc : gr) {
				String resultStr = formatCurve(loc, cc.getValues(loc));
				writer.write(resultStr);
				writer.newLine();
			}
		}
	}
	
	
	private static final Joiner J = Joiner.on(',').useForNull(" ");
	private static String formatCurve(Location loc, List<Double> values) {
		List<String> dat = Lists.newArrayList();
		dat.add(String.format("%.1f", loc.getLatitude()));
		dat.add(String.format("%.1f", loc.getLongitude()));
		for (double v : values) {
			dat.add(Double.toString(v));
		}
		return J.join(dat);
	}


	


}
