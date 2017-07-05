package scratch.peter.ucerf3.calc;

import static com.google.common.base.Charsets.US_ASCII;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.CPT_PARAM_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.LOG_PLOT_NAME;
import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.GRID_SPACING_PARAM_NAME;
import static org.opensha.nshmp2.tmp.TestGrid.*;
import static org.opensha.nshmp2.util.Period.*;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.*;
import static scratch.UCERF3.enumTreeBranches.FaultModels.*;
import static scratch.UCERF3.enumTreeBranches.InversionModels.*;
import static scratch.UCERF3.enumTreeBranches.MaxMagOffFault.*;
import static scratch.UCERF3.enumTreeBranches.MomentRateFixes.*;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.*;
import static scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels.*;
import static scratch.UCERF3.enumTreeBranches.SpatialSeisPDF.*;
import static scratch.UCERF3.enumTreeBranches.TotalMag5Rate.*;
import static scratch.peter.curves.ProbOfExceed.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.collections.MultiMap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.peter.curves.ProbOfExceed;
import scratch.peter.nshmp.CurveContainer;
import scratch.peter.nshmp.NSHMP_DataUtils;
import scratch.peter.nshmp.NSHMP_PlotUtils;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.io.Files;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class UC3_MapMaker {
	
	private static final String S = File.separator;
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static final String ROOT = "tmp/UC3maps/";
	private static final String SUFFIX = "_0p1";

	public static void main(String[] args) throws IOException {
//		generateBranchList();
//		buildMaps();
//		makeMultiBranchMap();
//		buildMapsUC32();
//		makeCoulombTestMaps();
//		makeLocalHazardMaps();
		makeLocalRatioMaps();
//		makeLocalMapsLoRes();
//		makeDefModelSpatialSeisMap();
//		makeLittleSalmonPlots();
//		makeLittleSalmonMultiPlots();
//		makeAleatoryMagScalingaRatio();
//		makeBrAvgFMDMratioMap();
		
//		TestGrid grid = CA_RELM;
//		ProbOfExceed pe = PE2IN50;
//		Period p = GM0P00;
//		String suffix = "";
//		double spacing = 0.1;
//		
//		
//		CurveContainer u3p2cc = CurveContainer.create(new File("/tmp/curves.csv"), grid, spacing);
//		
//		GriddedRegion gr = grid.grid(spacing);
//		GeoDataSet xyz = NSHMP_DataUtils.extractPE(u3p2cc, gr, pe);
//
//		String outDir = "/tmp/peter";
////		makeMapPlot(xyz, spacing, grid.bounds(), dlFile.getAbsolutePath(), "U3.2 Br AVg", false, true, false);
//		double[] minmax = NSHMP_PlotUtils.getRange(p);
//		GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
//		String label = pe + " " + p.getLabel() + " (g)";
//		String dlDir = outDir + "map-" + dlDirNameRatio(grid, pe, p, "all");
//		makeMapPlot(xyz, spacing, grid.bounds(), dlDir, label,
//			minmax[0], minmax[1], cpt, false, false);
	}
	
	private static void makeMultiBranchMap() throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;

//		String treeSrcDir = "FM-DM-MS-DSR-UV";
//		String branchList = "test.txt";
//		GeoDataSet over = loadMulti(treeSrcDir, branchList, pe, grid, p, SUFFIX);
		
//		String treeSrcDir = "FM-DM-MS-DSR-UV";
//		String branchList = "ref-U3-80.txt";
//		GeoDataSet over = mergeBranches(treeSrcDir, branchList, pe, grid, p);

//		String treeSrcDir2 = "FM-DM-MS-DSR-UV";
//		String branchList2 = "ref-U2-80.txt";
//		GeoDataSet under = mergeBranches(treeSrcDir2, branchList2, pe, grid, p);

//		String nshmpSrcDir = "nshmp_ca";
//		GeoDataSet under = nshmpReference(nshmpSrcDir, pe, grid, p);
		
		// =====
		
//		String treeSrcDir = "tree_refNoBG";
//		String branchList = "ref-U3-80.txt";
//		GeoDataSet over = mergeBranches(treeSrcDir, branchList, pe, grid, p);

//		String nshmpSrcDir = "nshmp_ca_nobg";
//		GeoDataSet under = nshmpReference(nshmpSrcDir, pe, grid, p);

		// =====
		
//		String treeSrcDir = "FM-DM-MS-DSR-UV";
//		String branchList = "ref-160.txt";
//		GeoDataSet over = mergeBranches(treeSrcDir, branchList, pe, grid, p);

//		String treeSrcDir2 = "FM-DM-MS-DSR-UV-M576-MX72";
//		String branchList2 = "ref-M576-MX72-160.txt";
//		GeoDataSet under = mergeBranches(treeSrcDir2, branchList2, pe, grid, p);

		// =====
		
//		String treeSrcDir2 = "FM-DM-MS-DSR-UV-M576-MX76";
//		String branchList2 = "ref-M576-MX76-U2-80.txt";
//		GeoDataSet over = mergeBranches(treeSrcDir2, branchList2, pe, grid, p);

//		String nshmpSrcDir = "nshmp_ca";
//		GeoDataSet under = nshmpReference(nshmpSrcDir, pe, grid, p);

		// =====
		
//		String uc3srcDir = "uc3uc2mapTAP_0p1";
//		GeoDataSet over = loadSingle(uc3srcDir, pe, grid, p);
//
//		String nshmpSrcDir = "muc2up_fm2p1_nobg_0p1"; //"nshmp_ca_nobg";
//		GeoDataSet under = loadSingle(nshmpSrcDir, pe, grid, p);
//
//		GeoDataSet xyzout = GeoDataSetMath.divide(over, under);
//		
//		String dlDir = ROOT + "maps/UC3UC2MAPTAP-UC2FM2P1-nobg/PGA-10p50-log";
//		makeRatioPlot(xyzout, grid.bounds(), dlDir, "10% in 50 : PGA", true);

		// ===== rate (glenn biasi) test comparison
		
//		// glenn biasi rates using 20 runs
//		String srcDir = "uc3rateTest";
//		GeoDataSet rateTest = loadSingle(srcDir, pe, grid, p);
//
//		// ref branch using 7 runs
//		srcDir = "FM-DM-MS-DSR-UV/FM3_2_ZENG_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_0p1";
//		GeoDataSet refBranch = loadSingle(srcDir, pe, grid, p);
//		
//		// NSHMP reference
//		srcDir = "nshmp_ca";
//		GeoDataSet nshmp = loadSingle(srcDir, pe, grid, p);
//
//		GeoDataSet rateOverRef = GeoDataSetMath.divide(rateTest, refBranch);
//		GeoDataSet rateOverNSHMP = GeoDataSetMath.divide(rateTest, nshmp);
//		GeoDataSet refOverNSHMP = GeoDataSetMath.divide(rateTest, nshmp);
//		
//		String dlDir = ROOT + "maps/cptTest3/rate_over_ref-PGA-2p50-log";
////		makeRatioPlot(rateOverRef, grid.bounds(), dlDir, "2% in 50 : PGA", true);
////		dlDir = ROOT + "maps/cptTest/rate_over_nshmp-PGA-2p50-log";
////		makeRatioPlot(rateOverNSHMP, grid.bounds(), dlDir, "2% in 50 : PGA", true);
//		dlDir = ROOT + "maps/cptTest4/ref_over_nshmp-PGA-2p50-log";
//		makeRatioPlot(refOverNSHMP, grid.bounds(), dlDir, "test", true, true);

		
		// base hazard maps
//		String srcDir = "nshmp_ca";
//		GeoDataSet nshmp = loadSingle(srcDir, pe, grid, p);
//		String outDir = "NSHMP";
		
//		String srcDir = ROOT + "src" + S + "UC32";
//		File branchListFile = new File(ROOT + "branchsetsUC32", "all.txt");
//		GeoDataSet nshmp = loadMulti(srcDir, branchListFile, pe, grid, p, "");
//		String outDir = "ALL";
//		
//		GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(GM0P00);
//		double[] minmax = NSHMP_PlotUtils.getRange(GM0P00);
//		String dlDir = "/Users/pmpowers/projects/OpenSHA/tmp/UC3maps/mapsUC32b/" + outDir;
//		makeMapPlot(nshmp, grid.bounds(), dlDir, "2% in 50 PGA (g)", minmax[0],
//			minmax[1], cpt, true);

	}
	
	private static GeoDataSet loadSingle(String src, ProbOfExceed pe, 
			TestGrid grid, Period p) {
		return loadSingle( src, pe, grid, p, 0.1);
	}

	private static GeoDataSet loadSingle(String src, ProbOfExceed pe, 
			TestGrid grid, Period p, double spacing) {
		File nshmp = new File(ROOT + "src" + S + src + S + grid + S + p + S + "curves.csv");
		CurveContainer nshmpcc = CurveContainer.create(nshmp, grid, spacing);
		GriddedRegion gr = grid.grid(spacing);
		return NSHMP_DataUtils.extractPE(nshmpcc, gr, pe);
	}

	/*
	 * Creates a geo data set by combining multiple weighted curve sets.
	 * 
	 * This method fixed a problem in UC31 maps wherein GeoDataSets for a
	 * PE were being created, weighted and added. Although the differences are
	 * minor, with a few exceptions, the correct calculation is to combine
	 * weighted curves, then compute a PE GeoDataSet.
	 */
	private static GeoDataSet loadMulti(String srcDir, File branchListFile,
			ProbOfExceed pe, TestGrid grid, Period p, String suffix)
			throws IOException {

		List<String> branchNames = Files.readLines(branchListFile, US_ASCII);
		System.out.println("Loading: " + branchListFile.getName());
		
		// create wt list (subbing ZENGBB [wt=0.3] for ZENG [wt=0.0])
		List<Double> wtList = Lists.newArrayList();
		for (String brName : branchNames) {
			LogicTreeBranch branch = LogicTreeBranch.fromFileName(brName);
			if (branch.getValue(DeformationModels.class).equals(ZENG)) {
				branch.setValue(ZENGBB);
			}
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
			// catch when referencing UC3.1 branches which use ZENG, not ZENGBB
			if (!brFile.exists()) brFile = new File(replaceZeng(brPath));
			CurveContainer cc = CurveContainer.create(brFile, grid, 0.1);
			cc.scale(wtList.get(idx++));

			if (mapcc == null) {
				mapcc = cc;
			} else {
				mapcc.add(cc);
			}
		}
		System.out.println();
		return NSHMP_DataUtils.extractPE(mapcc, gr, pe);
	}
	
	private static GeoDataSet loadMultiNoBG(String srcDir, String bgSrcDir,
			File branchListFile, ProbOfExceed pe, TestGrid grid, Period p) 
					throws IOException {

		List<String> branchNames = Files.readLines(branchListFile, US_ASCII);
		System.out.println("Loading: " + branchListFile.getName());
		
		// create wt list (subbing ZENGBB [wt=0.3] for ZENG [wt=0.0])
		List<Double> wtList = Lists.newArrayList();
		for (String brName : branchNames) {
			LogicTreeBranch branch = LogicTreeBranch.fromFileName(brName);
			if (branch.getValue(DeformationModels.class).equals(ZENG)) {
				branch.setValue(ZENGBB);
			}
			wtList.add(branch.getAprioriBranchWt());
		}
		DataUtils.asWeights(wtList);
		
		String cPath = grid + S + p + S + "curves.csv";
		GriddedRegion gr = grid.grid(0.1);
		CurveContainer mapcc = null;
		
		int idx = 0;
		for (String brName : branchNames) {
			if (idx % 100 == 0) System.out.print(idx + " ");
			String brID = brName;
			String totBrPath = srcDir + S + brID + S + cPath;
			String bgBrPath = bgSrcDir + S + brID + S + cPath;
			File totBrFile = new File(totBrPath);
			File bgBrFile = new File(bgBrPath);
			CurveContainer ccTot = CurveContainer.create(totBrFile, grid, 0.1);
			CurveContainer ccBg = CurveContainer.create(bgBrFile, grid, 0.1);
			ccTot.scale(wtList.get(idx));
			ccBg.scale(wtList.get(idx++));
			ccTot.subtract(ccBg);

			if (mapcc == null) {
				mapcc = ccTot;
			} else {
				mapcc.add(ccTot);
			}
		}
		System.out.println();
		return NSHMP_DataUtils.extractPE(mapcc, gr, pe);
	}
	
	private static GeoDataSet loadNSHMPnoFlt(ProbOfExceed pe, TestGrid grid, Period p) {
		File nshmpTotal = new File(ROOT + "src" + S + "nshmp_ca" + S + grid + S + p + S + "curves.csv");
		CurveContainer nshmpTotCC = CurveContainer.create(nshmpTotal, grid, 0.1);
		File nshmpFlt = new File(ROOT + "src" + S + "nshmp_ca_nobg" + S + grid + S + p + S + "curves.csv");
		CurveContainer nshmpFltCC = CurveContainer.create(nshmpFlt, grid, 0.1);
		nshmpTotCC.subtract(nshmpFltCC);
		GeoDataSet xyz = NSHMP_DataUtils.extractPE(nshmpTotCC, grid.grid(0.1), pe);
		return xyz;
	}
	

	private static String replaceZeng(String branchID) {
		StringBuffer sb = new StringBuffer(branchID);
		int pos = sb.indexOf("ZENGBB");
		if (pos >= 0) sb.delete(pos+4, pos+6);
		return sb.toString();
	}
	
//	private static GeoDataSet loadMultiOLD(String srcDir, String branchList,
//			ProbOfExceed pe, TestGrid grid, Period p) {
//		File branchListFile = new File(ROOT + "branchsets", branchList);
//		List<String> branchNames = null;
//		try {
//			branchNames = Files.readLines(branchListFile, US_ASCII);
//		} catch (IOException ioe) {
//			ioe.printStackTrace();
//		}
//		
//		// create wt list
//		
//		List<Double> wtList = Lists.newArrayList();
//		for (String brName : branchNames) {
//			LogicTreeBranch branch = LogicTreeBranch.fromFileName(brName);
//			wtList.add(branch.getAprioriBranchWt());
//		}
//		double[] wts = DataUtils.asWeights(Doubles.toArray(wtList));
//		
//		
//		String cPath = grid + S + p + S + "curves.csv";
//		GriddedRegion gr = grid.grid(0.1);
//		GeoDataSet mapxyz = null;
//		
//		int idx = 0;
//		for (String brName : branchNames) {
//			System.out.println(brName);
//			String brID = brName + SUFFIX;
//			String brPath = ROOT + "src" + S + srcDir + S + brID + S + cPath;
//			File brFile = new File(brPath);
//			CurveContainer cc = CurveContainer.create(brFile, grid);
//			GeoDataSet xyz = NSHMP_DataUtils.extractPE(cc, gr, pe);
//			xyz = NSHMP_GeoDataUtils.multiply(xyz, wts[idx++]);
//			if (mapxyz == null) {
//				mapxyz = xyz;
//			} else {
//				mapxyz = GeoDataSetMath.add(mapxyz, xyz);
//			}
//		}
//		
//		return mapxyz;
//	}

	
	/*
	 * Used to make ratio maps for UCERF3.2 report 
	 */
	private static void buildMapsUC32() throws IOException {
		List<String> brOverList = null;
		String brUnder = null;
		String branches = null;
		String uc31srcDir = null;
		String srcDir = ROOT + "src/UC32/";
		String outDir = ROOT + "mapsUC32b/";
		
		// Mmax and Mgt5 comparisons
		brOverList = Lists.newArrayList(
			"M576-MX72", "M587-MX72", "M510-MX72",
			"M576-MX76", "M587-MX76", "M510-MX76",
			"M576-MX80", "M587-MX80", "M510-MX80");
		brUnder = "M587-MX76";
		makeRatioMap(srcDir, outDir, brOverList, brUnder);
		
		
		// other branch node comparisons
		
		
//		brOverList = Lists.newArrayList(
//			"FM31", "FM32", 
//			"DM_ZBB", "DM_ABM", "DM_GEOL", "DM_NEOK",
//			"MS_SH09M", "MS_ELLB", "MS_ELLBSL", "MS_HB08", "MS_SHCSD",
//			"DSR_TAP", "DSR_UNI",
//			"M576", "M587", "M510",
//			"MX72", "MX76", "MX80",
//			"UC2", "UC3");
//		makeRatioMap(srcDir, outDir, brOverList);
		
//		// comparisons with UC3.1 logic tree subsets; M5 and MX variants
//		// 160 branches each
//		branches = "UC31cf-M587-MX76";
//		uc31srcDir = ROOT + "src/FM-DM-MS-DSR-UV/";
//		makeUC31cfMap(srcDir, uc31srcDir, outDir, branches);
//		
//		branches = "UC31cf-M576-MX72";
//		uc31srcDir = ROOT + "src/FM-DM-MS-DSR-UV-M576-MX72/";
//		makeUC31cfMap(srcDir, uc31srcDir, outDir, branches);
//
//		branches = "UC31cf-M576-MX76";
//		uc31srcDir = ROOT + "src/FM-DM-MS-DSR-UV-M576-MX76/";
//		makeUC31cfMap(srcDir, uc31srcDir, outDir, branches);
		
		
		// ratio maps of UC3.2 over NHSMP, faults only
		// subsets
//		String bgDir = ROOT + "src/UC32-FM-DM-MS-bgOnly";
//		brOverList = Lists.newArrayList("NoBG-DM-ABM", "NoBG-DM-GEOL",
//			"NoBG-DM-NEOK", "NoBG-DM-ZENGBB", "NoBG-FM31", "NoBG-FM32",
//			"NoBG-MS-ELLB", "NoBG-MS-ELLBSL", "NoBG-MS-HB08", "NoBG-MS-SH09M",
//			"NoBG-MS-SHCSD");
//		makeRatioMapNoBG(srcDir, bgDir, outDir, brOverList);
		
		// all 40 branches run
//		String bgDir = ROOT + "src/UC32-FM-DM-MS-bgOnly";
//		brOverList = Lists.newArrayList("NoBG-FM-DM-MS");
//		makeRatioMapNoBG(srcDir, bgDir, outDir, brOverList);
		
		// ratio maps just UC3.2 background models over NSHMP
		// UC32_noBg_combo / nshmp_ca - nshmap_ca_nobg
		
//		brOverList = Lists.newArrayList(
//			"NoFLT-FM-DM-MS-U2",
//			"NoFLT-FM-DM-MS-U3",
//			"NoFLT-FM-DM-MS-SS");
//		srcDir = ROOT + "src/UC32-FM-DM-MS-SS-bg/";
//		makeRatioMapsNoFLT(srcDir, outDir, brOverList);
		
		
//		// ratio map UC32 over NSMP (1440 branches)
//		branches = "all";
//		makeRatioMap(srcDir, outDir, branches);

//		branches = "UC3";
//		makeRatioMap(srcDir, outDir, branches);
		
//		// ratio map of UC3.2 over NSHMP at 1-sec for
//		// [FM]-[DM]-[MS]-TAP-M587-MX76-U2
//		// for comparison (to PGA) must make a PGA equivalent
//		String branchesFile = "UC32-FM-DM-MS-U2";
//		brOverList = Lists.newArrayList(branchesFile);
//		make1secRatioMap(
//			ROOT + "src/UC32-FM-DM-MS-U2-1sec/", outDir, brOverList);		
//		makeRatioMap(srcDir, outDir, brnchesFile);
		
		// ratio branchAvg over FM31 and FM32
		branches = "FM32"; //"testSingle"; //"FM31";
		String brAvgSrcDir = "UC32branchAvg/FM32";
		
		makeBrAvgRatioMap(brAvgSrcDir, srcDir, outDir, branches);
	}
	
	// creates hi-res hazard ratio maps of the SF and LA regions;NOTE these only
	// consider the branch avg solutions for FM3.2 as there were errors running
	// some FM3.1 models
	private static void makeLocalRatioMaps() {
	
		double spacing = 0.02;
		String srcDir = ROOT + "src/SF-LA-0p02/";
		String outDir = ROOT + "mapsUC32local/ratiosWithFaults/";
				
//		List<Period> allPeriods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
//		Multimap<ProbOfExceed, Period> pePeriodMap = ArrayListMultimap.create();
//		pePeriodMap.putAll(PE2IN50, allPeriods);
//		pePeriodMap.putAll(PE10IN50, allPeriods);
//		pePeriodMap.put(PE40IN50, GM0P00);
//		List<String> srcTypes = Lists.newArrayList("all", "flt", "bg");
//		List<TestGrid> grids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		
//		for (ProbOfExceed pe : pePeriodMap.keySet()) {
//			for (Period p : pePeriodMap.get(pe)) {
//				for (String srcType : srcTypes) {
//					for (TestGrid grid : grids) {
//						makeLocalRatioMap(srcDir, srcType, grid, spacing, pe, p, outDir, true);
////						System.out.println(srcType + " " + grid + " " + pe + " " + p);
//					}
//				}
//			}
//		}
		makeLocalRatioMap(srcDir, "all", SAN_FRANCISCO, spacing, PE2IN50, GM0P00, outDir, true);

	}
	
	private static void makeLocalRatioMap(String srcDir,
			String srcType, TestGrid grid, double spacing, ProbOfExceed pe,
			Period p, String dlDir, boolean showFaults) {
		
		Lists.newArrayList("UC32brAvg-FM32", "nshmp_ca");
		
		// hazard curve source
		String UC32fm31srcPath = srcDir + "UC32brAvg-FM31-" + srcType + S + grid + S + p + S + "curves.csv";
		String UC32fm32srcPath = srcDir + "UC32brAvg-FM32-" + srcType + S + grid + S + p + S + "curves.csv";
		String NSHMsrcPath = srcDir + "nshmp_ca-" + srcType + S + grid + S + p + S + "curves.csv";
		File UC32fm31src = new File(UC32fm31srcPath);
		File UC32fm32src = new File(UC32fm32srcPath);
		File NSHMsrc = new File(NSHMsrcPath);
		
		// for two fault models
		CurveContainer fm31cc = CurveContainer.create(UC32fm31src, grid, spacing);
		CurveContainer fm32cc = CurveContainer.create(UC32fm32src, grid, spacing);
		fm32cc.add(fm31cc);
		fm32cc.scale(0.5);
		
		// curves
		// for one fault model
		// CurveContainer UC32cc = CurveContainer.create(UC32src, grid, spacing);
		CurveContainer NSHMcc = CurveContainer.create(NSHMsrc, grid, spacing);
		
		// data
		GriddedRegion gr = grid.grid(spacing);
		GeoDataSet UC32xyz = NSHMP_DataUtils.extractPE(fm32cc, gr, pe);
		GeoDataSet NSHMxyz = NSHMP_DataUtils.extractPE(NSHMcc, gr, pe);

		// ratio map
		GeoDataSet xyz = GeoDataSetMath.divide(UC32xyz, NSHMxyz);
		String label = "Ratio " + pe + " " + p.getLabel() + " (g)";
		dlDir += dlDirNameRatio(grid, pe, p, srcType) + "-p1";
		
		makeRatioPlotNSHMP(xyz, spacing, grid.bounds(), dlDir, label, showFaults);

	}

	// local loRes maps for comparison with hiRes versions
	private static void makeLocalMapsLoRes() {
		
		double spacing = 0.1;
		Period p = GM0P20;
		ProbOfExceed pe = PE2IN50;
		
		String srcDir = ROOT + "src/";
		String outDir = ROOT + "mapsUC32localLoRes/";
				
		String UC32path = "UC32branchAvg/FM32/CA_RELM/" + p + "/curves.csv";
		String NSHMpath = "nshmp_ca/CA_RELM/" + p + "/curves.csv";
		
		File UC32src = new File(srcDir + UC32path);
		File NSHMsrc = new File(srcDir + NSHMpath);

		List<TestGrid> grids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		
		for (TestGrid grid : grids) {
			System.out.println();
			
			CurveContainer UC32cc = CurveContainer.create(UC32src, grid, spacing);
			CurveContainer NSHMcc = CurveContainer.create(NSHMsrc, grid, spacing);
			
			// data
			GriddedRegion gr = grid.grid(spacing);
			GeoDataSet UC32xyz = NSHMP_DataUtils.extractPE(UC32cc, gr, pe);
			GeoDataSet NSHMxyz = NSHMP_DataUtils.extractPE(NSHMcc, gr, pe);
			
			// hazard map
			double[] minmax = NSHMP_PlotUtils.getRange(p);
			GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
			String label = pe + " " + p.getLabel() + " (g)";
			String dlDir = outDir + "map-" + dlDirNameRatio(grid, pe, p, "all");
			makeMapPlot(UC32xyz, spacing, grid.bounds(), dlDir, label,
				minmax[0], minmax[1], cpt, false, false);
			
			// ratio map
			GeoDataSet xyz = GeoDataSetMath.divide(UC32xyz, NSHMxyz);
			label = "Ratio " + pe + " " + p.getLabel() + " (g)";
			dlDir = outDir + "ratio-" + dlDirNameRatio(grid, pe, p, "all-p1");
			makeRatioPlot(xyz, spacing, grid.bounds(), dlDir, label, true, true, false);
		}
	}
	
	// creates hi-res hazard maps of the SF and LA regions; NOTE these only
	// consider the branch avg solutions for FM3.2 as there were errors running
	// some FM3.1 models
	private static void makeLocalHazardMaps() {

		double spacing = 0.02;
		String srcDir = ROOT + "src/SF-LA-0p02/";
		String outDir = ROOT + "mapsUC32localMapTest/";
				
		List<Period> allPeriods = Lists.newArrayList(GM0P00, GM0P20, GM1P00);
		Multimap<ProbOfExceed, Period> pePeriodMap = ArrayListMultimap.create();
		pePeriodMap.putAll(PE2IN50, allPeriods);
		pePeriodMap.putAll(PE10IN50, allPeriods);
		pePeriodMap.put(PE40IN50, GM0P00);
		List<String> srcIDs = Lists.newArrayList("UC32brAvg", "nshmp_ca");
		List<String> srcTypes = Lists.newArrayList("all", "flt", "bg");
		List<TestGrid> grids = Lists.newArrayList(SAN_FRANCISCO, LOS_ANGELES);
		
//		for (ProbOfExceed pe : pePeriodMap.keySet()) {
//			for (Period p : pePeriodMap.get(pe)) {
//				for (String srcID : srcIDs) {
//					for (String srcType : srcTypes) {
//						for (TestGrid grid : grids) {
//							System.out.println(srcID + "\t" + srcType + "\t" + grid + "\t" + pe + "\t" + p);
//							makeLocalHazardMap(srcDir, srcID, srcType, grid, spacing, pe, p, outDir);
//						}
//					}
//				}
//			}
//		}
		
		makeLocalHazardMap(srcDir, "UC32brAvg", "all", SAN_FRANCISCO, spacing, PE2IN50, GM0P00, outDir, true);
	}
	
	private static void makeLocalHazardMap(String srcDir, String srcID,
			String srcType, TestGrid grid, double spacing, ProbOfExceed pe,
			Period p, String dlDir, boolean showFaults) {
		
		// hazard curve source
		CurveContainer cc = null;
		
		// if brAvg, combine Fault Models
		if (srcID.equals("UC32brAvg")) {
			String fm31srcPath = srcDir + srcID + "-FM31-" + srcType + S + grid + S + p + S + "curves.csv";
			String fm32srcPath = srcDir + srcID + "-FM32-" + srcType + S + grid + S + p + S + "curves.csv";
			File fm31src = new File(fm31srcPath);
			File fm32src = new File(fm32srcPath);
			CurveContainer fm31cc = CurveContainer.create(fm31src, grid, spacing);
			CurveContainer fm32cc = CurveContainer.create(fm32src, grid, spacing);
			fm32cc.add(fm31cc);
			fm32cc.scale(0.5);
			cc = fm32cc;
		} else {
			String srcPath = srcDir + srcID + "-" + srcType + S;
			srcPath += grid + S + p + S + "curves.csv";
			File src = new File(srcPath);
			cc = CurveContainer.create(src, grid, spacing);
		}

		// data
		GriddedRegion gr = grid.grid(spacing);
		GeoDataSet xyz = NSHMP_DataUtils.extractPE(cc, gr, pe);

		// map
		double[] minmax = NSHMP_PlotUtils.getRange(p);
		GMT_CPT_Files cpt = NSHMP_PlotUtils.getCPT(p);
		String label = pe + " " + p.getLabel() + " (g)";
		dlDir += dlDirNameMap(srcID, grid, pe, p, srcType) + S;
		makeMapPlot(xyz, spacing, grid.bounds(), dlDir, label,
			minmax[0], minmax[1], cpt, false, showFaults);
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
	
	private static String dlDirNameMap(String srcID, TestGrid tg, ProbOfExceed pe, Period p, String srcType) {
		StringBuffer sb = new StringBuffer();
		sb.append(srcID.startsWith("UC32brAvg") ? "UC32brAvg" : "NSHMP");
		sb.append("-");
		sb.append(tg == SAN_FRANCISCO ? "SF" : "LA");
		sb.append("-");
		sb.append(pe == PE2IN50 ? "2in50" : pe == PE10IN50 ? "10in50" : "40in50");
		sb.append("-");
		sb.append(p == GM0P00 ? "PGA" : p == GM0P20 ? "5Hz" : "1Hz");
		sb.append("-");
		sb.append(srcType);
		return sb.toString();
	}
	
	private static String dlDirNameRatio(TestGrid tg, ProbOfExceed pe, Period p, String srcType) {
		StringBuffer sb = new StringBuffer();
		sb.append(tg == SAN_FRANCISCO ? "SF" : "LA");
		sb.append("-");
		sb.append(pe == PE2IN50 ? "2in50" : pe == PE10IN50 ? "10in50" : "40in50");
		sb.append("-");
		sb.append(p == GM0P00 ? "PGA" : p == GM0P20 ? "5Hz" : "1Hz");
		sb.append("-");
		sb.append(srcType);
		return sb.toString();
	}

	// ratio map testing if UC3 halos might contaminate Deformation Model
	// based gridded seis sources.
	private static void makeDefModelSpatialSeisMap() {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE10IN50;
		Period p = GM0P00;
		String dlDir = ROOT + "mapsUC32b/DM_spatialSeisTest/p3";
		String srcRoot = "UC32-DMspatialTest" + S;
		String srcU2 = srcRoot + "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU2";
		String srcU3 = srcRoot + "FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3";
		GeoDataSet xyzU2 = loadSingle(srcU2, pe, grid, p);
		GeoDataSet xyzU3 = loadSingle(srcU3, pe, grid, p);
		GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzU3, xyzU2);
		makeRatioPlot(xyzRatio, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
	}
	
	// Ratio of two conv run ground motion maps for Little Salmon onshore.
	// First pass at this used run01 (low rate) and run83 (high rate) for
	// rupture 73204. Second pass expanded to run00 to run22 in the
	// numerator. These tests used UC31 conv runs.
	private static void makeLittleSalmonPlots() {
		TestGrid grid = LITTLE_SALMON;
		Period p = GM0P00;
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50, PE1IN100);
		String srcRoot = "UC32convLittleSalmonSingle" + S;
		String srcUnder = srcRoot + "sol83";
		for (ProbOfExceed pe : PEs) {
			GeoDataSet xyzUnder = loadSingle(srcUnder, pe, grid, p, 0.02);
			for (int i=0; i<23; i++) {
				String solNum = String.format("%02d", i);
				String id  = "sol" + solNum + "_sup_sol83";
				String dlDir = ROOT + "mapsUC32b/UC32littleSalmon/single-UC31/" + pe + S +
					id + S;
				String srcOver = srcRoot + "sol" + solNum;
				GeoDataSet xyzOver = loadSingle(srcOver, pe, grid, p, 0.02);
				GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
				makeRatioPlot(xyzRatio, 0.02, grid.bounds(), dlDir,
					pe.toString() + " ratio " + solNum + "/83", true, true,
					true);
			}
		}
	}

	// Similar to above but using sets of 5 runs (i.e. run0 = average of
	// run00 to run04).  These tests used UC32 conv runs.
	private static void makeLittleSalmonMultiPlots() {
		TestGrid grid = LITTLE_SALMON;
		Period p = GM0P00;
		List<ProbOfExceed> PEs = Lists.newArrayList(PE2IN50, PE1IN100);
		String srcRoot = "UC32convLittleSalmonMulti" + S;
		String srcUnder = srcRoot + "sol9";
		int[] indices = new int[] {0,3,6,12,15,18};
		for (ProbOfExceed pe : PEs) {
			GeoDataSet xyzUnder = loadSingle(srcUnder, pe, grid, p, 0.02);
			for (int index : indices) {
				String id  = "sol" + index + "_sup_sol9";
				String dlDir = ROOT + "mapsUC32b/UC32littleSalmon/multi-UC32/" + pe + S +
					id + S;
				String srcOver = srcRoot + "sol" + index;
				GeoDataSet xyzOver = loadSingle(srcOver, pe, grid, p, 0.02);
				GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
				makeRatioPlot(xyzRatio, 0.02, grid.bounds(), dlDir,
					pe.toString() + " ratio sol" + index + "/sol9", true, true,
					true);
			}
		}
	}
	
	// Creates a ratio plot of branch averaged UC32 with aleatory uncertainty
	// on magnitude enabled over disabled; both fault models
	private static void makeAleatoryMagScalingaRatio() {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE10IN50;
		Period p = GM0P00;
		double spacing = 0.1;
		String root = ROOT + "src/UC32branchAvg/";
		String dlDir = ROOT + "mapsUC32b/UC32aleaMagTest";

		// load +alea sources
		String fm31on = root + "FM31alea/" + grid + S + p + S + "curves.csv";
		String fm32on = root + "FM32alea/" + grid + S + p + S + "curves.csv";
		String fm31off = root + "FM31/" + grid + S + p + S + "curves.csv";
		String fm32off = root + "FM32/" + grid + S + p + S + "curves.csv";

		File fm31onSrc = new File(fm31on);
		File fm32onSrc = new File(fm32on);
		File fm31offSrc = new File(fm31off);
		File fm32offSrc = new File(fm32off);

		CurveContainer fm31ccOn = CurveContainer.create(fm31onSrc, grid, spacing);
		CurveContainer fm32ccOn = CurveContainer.create(fm32onSrc, grid, spacing);
		CurveContainer fm31ccOff = CurveContainer.create(fm31offSrc, grid, spacing);
		CurveContainer fm32ccOff = CurveContainer.create(fm32offSrc, grid, spacing);

		fm32ccOn.add(fm31ccOn);
		fm32ccOn.scale(0.5);
		fm32ccOff.add(fm31ccOff);
		fm32ccOff.scale(0.5);

		GriddedRegion gr = grid.grid(spacing);
		GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(fm32ccOn, gr, pe);
		GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(fm32ccOff, gr, pe);
		GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
		
		makeRatioPlot(xyzRatio, 0.1, grid.bounds(), dlDir, "aleaOn/aleaOff", 
			true, true, false);
	}
		
		
	// creates ratios maps of coulomb variants to reference branch
	private static void makeCoulombTestMaps() {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
	
		// load single refBranch
//		String refBrSrc = "UC32/FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3";
//		GeoDataSet xyzRef = loadSingle(refBrSrc, pe, grid, p);
		
		String name = null;
		String cID = null;
		String dlDir = null;
		GeoDataSet xyz = null;
		GeoDataSet xyzRatio = null;
		
		name = "coulomb0.0";
		cID = "UC32coulombTest/" + name;
		dlDir = ROOT + "mapsUC32b/coulombTest/" + name;
		GeoDataSet xyzRef = loadSingle(cID, pe, grid, p);
//		xyzRatio = GeoDataSetMath.divide(xyz, xyzRef);
//		makeRatioPlot(xyzRatio, grid.bounds(), dlDir, "hazard ratio", true, true);
		
		name = "coulomb0.1";
		cID = "UC32coulombTest/" + name;
		dlDir = ROOT + "mapsUC32b/coulombTest/" + name + "_sup_0.0";
		xyz = loadSingle(cID, pe, grid, p);
		xyzRatio = GeoDataSetMath.divide(xyz, xyzRef);
		makeRatioPlot(xyzRatio, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
		
		name ="coulomb0.05";
		cID = "UC32coulombTest/" + name;
		dlDir = ROOT + "mapsUC32b/coulombTest/" + name + "_sup_0.0";
		xyz = loadSingle(cID, pe, grid, p);
		xyzRatio = GeoDataSetMath.divide(xyz, xyzRef);
		makeRatioPlot(xyzRatio, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);

	}
	
	// UCERF3.2 branchAvg over branches
	private static void makeBrAvgRatioMap(String srcDirForAvg,
			String SrcDirForTree, String outDir, String brUnder)
			throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		
		GeoDataSet xyzOver = loadSingle(srcDirForAvg, pe, grid, p);
	
		File brUnderFile = new File(ROOT + "branchsetsUC32", brUnder + ".txt");
		GeoDataSet xyzUnder = loadMulti(SrcDirForTree, brUnderFile, pe, grid, p, suffix);

		GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);

		String dlDir = outDir + srcDirForAvg;
		File dlFile = new File(dlDir);
		dlFile.mkdirs();
		makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
	}

	
	// UCERF3.2 branchAvg (DM and FM variants combined) over all branches
	private static void makeBrAvgFMDMratioMap() throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		double spacing = 0.1;
		String dlDir = ROOT + "mapsUC32b/UC32branchAvg/allFM/";

		
		// load Dm-FM branch avg solutions
		CurveContainer brAvgCC = null;
		String root = ROOT + "src/UC32branchAvg/";
		Set<DeformationModels> DMs = EnumSet.of(ABM, GEOLOGIC, NEOKINEMA, ZENGBB);
//		for (DeformationModels dm : DMs) {
					
			// load sources
//			String fm31 = root + "FM31-" + dm.getShortName() + S + grid + S + p + S + "curves.csv";
//			String fm32 = root + "FM32-" + dm.getShortName() + S + grid + S + p + S + "curves.csv";
			String fm31 = root + "FM31" + S + grid + S + p + S + "curves.csv";
			String fm32 = root + "FM32" + S + grid + S + p + S + "curves.csv";
			File fm31src = new File(fm31);
			File fm32src = new File(fm32);
			CurveContainer fm31cc = CurveContainer.create(fm31src, grid, spacing);
			CurveContainer fm32cc = CurveContainer.create(fm32src, grid, spacing);

			// sum, scale (FM), and scale (DM)
			fm32cc.add(fm31cc);
			fm32cc.scale(0.5); // * dm.getRelativeWeight(null));
			
//			if (brAvgCC == null) {
//				brAvgCC = fm32cc;
//			} else {
//				brAvgCC.add(fm32cc);
//			}
//		}
		GriddedRegion gr = grid.grid(spacing);
		GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(fm32cc, gr, pe);
		
		// full logic tree
		File underFile = new File(ROOT + "branchsetsUC32", "all.txt");
		GeoDataSet xyzUnder = loadMulti(ROOT + "src/UC32/", underFile, pe, grid, p, suffix);

		GeoDataSet xyzRatio = GeoDataSetMath.divide(xyzOver, xyzUnder);
		File dlFile = new File(dlDir);
		dlFile.mkdirs();
		makeRatioPlot(xyzRatio, spacing, grid.bounds(), dlDir, "brAvg/fullTree", true, true, false);
	}

	private static GeoDataSet UC32xyz;
	
	// UCERF3.2 node ratio maps multithreaded
	private static void makeRatioMap(String srcDir, String outDir,
			List<String> brOverList) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		
		if (UC32xyz == null) {
			File brFile = new File(ROOT + "branchsetsUC32", "all.txt");
			UC32xyz = loadMulti(srcDir, brFile, pe, grid, p, suffix);
		}
		
		int numProc = Runtime.getRuntime().availableProcessors();
		ExecutorService ex = Executors.newFixedThreadPool(numProc);

		for (String brOver : brOverList) {
			RatioMapMaker rmm = new RatioMapMaker(UC32xyz, brOver, srcDir,
				outDir, pe, grid, p, suffix);
			ex.submit(rmm);
		}
		ex.shutdown();
	}

	// UCERF3.2 node ratio mapsx
	private static void makeRatioMap(String srcDir, String outDir,
			List<String> brOverList, String brUnder) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		String suffix = "";
		
		File brUnderFile = new File(ROOT + "branchsetsUC32", brUnder + ".txt");
		GeoDataSet xyzUnder = loadMulti(srcDir, brUnderFile, pe, grid, p, suffix);

		for (String brOver : brOverList) {
			File brOverFile = new File(ROOT + "branchsetsUC32", brOver + ".txt");
			GeoDataSet xyzOver = loadMulti(srcDir, brOverFile, pe, grid, p, suffix);
		
			GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);

			String dlDir = outDir + brOver + "_sup_" + brUnder;
			makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
		}
	}
	
	private static class RatioMapMaker implements Callable<Void> {

		GeoDataSet xyzRef;
		String branches;
		String srcDir;
		String outDir;
		ProbOfExceed pe;
		TestGrid grid;
		Period p;
		String suffix;
		
		RatioMapMaker(GeoDataSet xyzRef, String branches, String srcDir, String outDir,
			ProbOfExceed pe, TestGrid grid, Period p, String suffix) {
			this.xyzRef = xyzRef;
			this.branches = branches;
			this.srcDir = srcDir;
			this.outDir = outDir;
			this.pe = pe;
			this.grid = grid;
			this.p = p;
			this.suffix = suffix;
		}
		
		@Override
		public Void call() throws IOException {
			File brFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
			GeoDataSet xyz = loadMulti(srcDir, brFile, pe, grid, p, suffix);
			GeoDataSet xyzRatio = GeoDataSetMath.divide(xyz, xyzRef);
			String dlDir = outDir + branches + "_sup_REF";
			makeRatioPlot(xyzRatio, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
			return null;
		}
		
	}
	
	// UCERF3.1 to 3.2 comparison maps
	private static void makeUC31cfMap(String srcDirUC32, String srcDirUC31, 
			String outDir, String branches) throws IOException {
		
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		
		File branchFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
		GeoDataSet xyz31 = loadMulti(srcDirUC31, branchFile, pe, grid, p, SUFFIX);
		GeoDataSet xyz32 = loadMulti(srcDirUC32, branchFile, pe, grid, p, "");

		GeoDataSet xyz = GeoDataSetMath.divide(xyz32, xyz31);
		String dlDir = outDir + branches;
		makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, true, false);
	}

	// UCERF3.2 NSHMP ratio maps
	private static void makeRatioMap(String srcDir, String outDir,
			String branches) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		
		GeoDataSet xyzNSHMP = loadSingle("nshmp_ca", pe, grid, p);
		File branchFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
		GeoDataSet xyzUC32 = loadMulti(srcDir, branchFile, pe, grid, p, "");
		GeoDataSet xyz = GeoDataSetMath.divide(xyzUC32, xyzNSHMP);
		String dlDir = outDir + branches + "-LTwtVars_sup_NSHMP";
		makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, false, false);
	}

	// UCERF3.2 NSHMP fault ratio maps (no bg)
	private static void makeRatioMapNoBG(String srcDir, String bgDir, String outDir,
			List<String> branchFileList) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		
		GeoDataSet xyzNSHMP = loadSingle("nshmp_ca_nobg", pe, grid, p);

		for (String branches : branchFileList) {
			File branchFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
			GeoDataSet xyzUC32 = loadMultiNoBG(srcDir, bgDir,branchFile, pe, grid, p);
			GeoDataSet xyz = GeoDataSetMath.divide(xyzUC32, xyzNSHMP);
			String dlDir = outDir + branches + "_sup_NSHMP_PGA";
			makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, false, false);
		}
	}
	
	// UCERF3.2 NSHMP bg ratio maps (no flt)
	private static void makeRatioMapsNoFLT(String srcDir, String outDir,
			List<String> branchFileList) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM0P00;
		
		GeoDataSet xyzNSHMP = loadNSHMPnoFlt(pe, grid, p);

		for (String branches : branchFileList) {
			File branchFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
			GeoDataSet xyzUC32 = loadMulti(srcDir, branchFile, pe, grid, p, "");
			GeoDataSet xyz = GeoDataSetMath.divide(xyzUC32, xyzNSHMP);
			String dlDir = outDir + branches + "_sup_NSHMP_BG";
			makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, false, false);
		}
	}

	// UCERF3.2 NSHMP 1sec ratio map
	private static void make1secRatioMap(String srcDir, String outDir,
			List<String> branchFileList) throws IOException {
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		Period p = GM1P00;
		
		GeoDataSet xyzNSHMP = loadSingle("nshmp_ca",pe, grid, p);

		for (String branches : branchFileList) {
			File branchFile = new File(ROOT + "branchsetsUC32", branches + ".txt");
			GeoDataSet xyzUC32 = loadMulti(srcDir, branchFile, pe, grid, p, "");
			GeoDataSet xyz = GeoDataSetMath.divide(xyzUC32, xyzNSHMP);
			String dlDir = outDir + branches + "_sup_NSHMP_1HZ";
			makeRatioPlot(xyz, 0.1, grid.bounds(), dlDir, "hazard ratio", true, false, false);
		}
	}

//	/*
//	 * Used to make ratio maps for UCERF3.1 report 
//	 */
//	private static void buildMaps() {
//		String mapGrpName = null;
//		String refDir = null;
//		String brDir = null;
//		boolean log = true;
//		try {
//			mapGrpName = "FM31-[DM]-[MS]-TAP-87-76-U3";
//			refDir = "FM-DM-MS-DSR";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM31-[DM]-[MS]-UNI-87-76-U3";
//			refDir = "FM-DM-MS-DSR";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM32-[DM]-[MS]-TAP-87-76-U3";
//			refDir = "FM-DM-MS-DSR";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM32-[DM]-[MS]-UNI-87-76-U3";
//			refDir = "FM-DM-MS-DSR";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "[FM]-ZENG-S09M-[DSR]-87-76-U3";
//			refDir = "FM-DM-MS-DSR";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM31-[DM]-S09M-TAP-87-76-U2";
//			refDir = "FM-DM-MS-UV";
//			brDir = "FM-DM-MS-UV";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM31-[DM]-S09M-TAP-87-76-U3";
//			refDir = "FM-DM-MS-UV";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM32-[DM]-S09M-TAP-87-76-U2";
//			refDir = "FM-DM-MS-UV";
//			brDir = "FM-DM-MS-UV";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM32-[DM]-S09M-TAP-87-76-U3";
//			refDir = "FM-DM-MS-UV";
//			brDir = "FM-DM-MS-DSR";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM31-ZENG-S09M-TAP-[M5]-[MM]-U3";
//			refDir = "FM-M5-MM";
//			brDir = "FM-M5-MM";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//			mapGrpName = "FM32-ZENG-S09M-TAP-[M5]-[MM]-U3";
//			refDir = "FM-M5-MM";
//			brDir = "FM-M5-MM";
//			makeRatioMaps(refDir, brDir, mapGrpName, log);
//
//		} catch (IOException ioe) {
//			ioe.printStackTrace();
//		}
//
//	}

	// make ratio maps from list of logic tree branches
	
	// for each map set need list and reference
	
	// dir = FM-DM-MS-DSR
	// FM31-[DM]-[MS]-TAP-87-76-U3 / FM31-ZENG-S09M-TAP-87-76-U3
	// FM31-[DM]-[MS]-UNI-87-76-U3 / FM31-ZENG-S09M-UNI-87-76-U3
	// FM32-[DM]-[MS]-TAP-87-76-U3 / FM32-ZENG-S09M-TAP-87-76-U3
	// FM32-[DM]-[MS]-UNI-87-76-U3 / FM32-ZENG-S09M-UNI-87-76-U3
	
	// dir = FM-DM-MS-DSR
	// [FM]-ZENG-S09M-[DSR]-87-76-U3 / FM31-ZENG-S09M-TAP-87-76-U3
	
	// dir = FM-DM-MS-UV (for U2) FM-DM-MS-DSR (for U3)
	// FM31-[DM]-S09M-TAP-87-76-U2 / FM31-ZENG-S09M-TAP-87-76-U2
	// FM31-[DM]-S09M-TAP-87-76-U3 / FM31-ZENG-S09M-TAP-87-76-U2
	// FM32-[DM]-S09M-TAP-87-76-U2 / FM32-ZENG-S09M-TAP-87-76-U2
	// FM32-[DM]-S09M-TAP-87-76-U3 / FM32-ZENG-S09M-TAP-87-76-U2
	
	// dir = FM-M5-MM
	// FM31-ZENG-S09M-TAP-[M5]-[MM]-U3 / FM31-ZENG-S09M-TAP-87-76-U3
	// FM32-ZENG-S09M-TAP-[M5]-[MM]-U3 / FM32-ZENG-S09M-TAP-87-76-U3
	
	// for UC3.1
	private static void makeRatioMaps(String refDir, String brDir,
			String listName, boolean log) throws IOException {
		
		File branchListFile = new File(ROOT + "branchsets", listName + ".txt");
		List<String> branchNames = Files.readLines(branchListFile, US_ASCII);
		TestGrid grid = CA_RELM;
		ProbOfExceed pe = PE2IN50;
		
		// reference branch is first in list
		// make ratios of all subsequent branches
		String refID = branchNames.get(0) + SUFFIX;
		String cPath = grid + S + GM0P00 + S + "curves.csv";
		String refPath = ROOT + "src" + S + refDir + S + refID + S + cPath;
		File refFile = new File(refPath);
		
		for (String brName : Iterables.skip(branchNames, 1)) {
			String dlDir = ROOT + "mapsUC32" + S + listName + (log ? "-log" : "") + S + brName + S;
			String brID = brName + SUFFIX;
			String brPath = ROOT + "src" + S + brDir + S + brID + S + cPath;
			File brFile = new File(brPath);
			
			makeRatioMap(dlDir, brFile, refFile, grid, pe, "Ratio", log, true, false);
		}
		
	}
	
	private static void makeRatioMap(String dlDir, File fOver, File fUnder,
			TestGrid grid, ProbOfExceed pe, String title,
			boolean log, boolean smooth, boolean showFaults) {
		
		double spacing = 0.1;
		GriddedRegion gr = grid.grid(spacing);
		CurveContainer cc = null;
		cc = CurveContainer.create(fOver, grid, spacing);
		GeoDataSet xyzOver = NSHMP_DataUtils.extractPE(cc, gr, pe);
		cc = CurveContainer.create(fUnder, grid, spacing);
		GeoDataSet xyzUnder = NSHMP_DataUtils.extractPE(cc, gr, pe);
		
		GeoDataSet xyz = GeoDataSetMath.divide(xyzOver, xyzUnder);
		makeRatioPlot(xyz, spacing, grid.bounds(), dlDir, title, log, smooth, showFaults);
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

	public static void makeRatioPlot(GeoDataSet xyz, double spacing, double[] bounds,
			String dlDir, String title, boolean log, boolean smooth, boolean showFaults) {
		double scale = log ? 0.1 : 0.2;
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
//				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
			
			// copy map.pdf to dl dir
			File pdfFrom = new File(dlDir, "map.pdf");
			File pdfTo = new File(dlDir + ".pdf");
			Files.copy(pdfFrom, pdfTo);
		} catch (IOException ioe) {
			ioe.printStackTrace();
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
//				addFaultTraces(FaultModels.FM3_2, map, Color.BLACK);
			}
			NSHMP_PlotUtils.makeMap(map, mapGen, "No metadata", dlDir);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}



	private static void generateBranchList() {
		
		String fileName = "tree1440";
		Set<FaultModels> fltModels = EnumSet.of(
			FM3_1, FM3_2); //FM3_2); // FM3_1, FM3_2);
		Set<DeformationModels> defModels = EnumSet.of(
//			ZENGBB);
			ABM, GEOLOGIC, NEOKINEMA, ZENGBB);
		Set<ScalingRelationships> scalingRel = EnumSet.of(
//			SHAW_2009_MOD);
			ELLSWORTH_B, ELLB_SQRT_LENGTH, HANKS_BAKUN_08,
			SHAW_CONST_STRESS_DROP, SHAW_2009_MOD);
		Set<SlipAlongRuptureModels> slipRup = EnumSet.of(
			UNIFORM, TAPERED); //UNIFORM, TAPERED);
		Set<InversionModels> invModels = EnumSet.of(
			CHAR_CONSTRAINED);
		Set<TotalMag5Rate> totM5rate = EnumSet.of(
			RATE_6p5, RATE_7p9, RATE_9p6);
		Set<MaxMagOffFault> mMaxOff = EnumSet.of(
			MAG_7p3, MAG_7p6, MAG_7p9); // MAG_7p2, MAG_7p6, MAG_8p0);
		Set<MomentRateFixes> momentFix = EnumSet.of(
			NONE);
		Set<SpatialSeisPDF> spatialSeis = EnumSet.of(
			UCERF2, UCERF3); // UCERF2, UCERF3);

		List<Set<? extends LogicTreeBranchNode<?>>> branchSets = Lists.newArrayList();
		branchSets.add(fltModels);
		branchSets.add(defModels);
		branchSets.add(scalingRel);
		branchSets.add(slipRup);
		branchSets.add(invModels);
		branchSets.add(totM5rate);
		branchSets.add(mMaxOff);
		branchSets.add(momentFix);
		branchSets.add(spatialSeis);
		
		int count = 0;
		Set<List<LogicTreeBranchNode<?>>> branches = Sets.cartesianProduct(branchSets);
		try {
		File out = new File("tmp/UC33/curvejobs/branches", fileName + ".txt");
		Files.write("", out, US_ASCII);
		for (List<LogicTreeBranchNode<?>> branch : branches) {
			LogicTreeBranch ltb = LogicTreeBranch.fromValues(branch);
			Files.append(ltb.buildFileName() + LF, out, US_ASCII);
			System.out.println((count++) + " " + ltb.buildFileName());
		}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	
	// make composite map across FM-DM-MS-DSR
	// need to get weight for each map and normalize wprior to summing
	
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

}
