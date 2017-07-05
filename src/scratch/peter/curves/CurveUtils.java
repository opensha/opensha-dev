package scratch.peter.curves;

import static com.google.common.base.Charsets.US_ASCII;
import static com.google.common.base.Preconditions.checkArgument;
import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.sra.rtgm.RTGM.Frequency.SA_0P20;
import static org.opensha.sra.rtgm.RTGM.Frequency.SA_1P00;
import static scratch.peter.curves.ProbOfExceed.PE10IN50;
import static scratch.peter.curves.ProbOfExceed.PE2IN50;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils;
import org.opensha.nshmp.NEHRP_TestCity;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.nshmp2.util.Period;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.SimpleFaultData;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import scratch.UCERF3.AverageFaultSystemSolution;
import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.peter.nshmp.CurveContainer;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.ArrayTable;
import com.google.common.collect.BiMap;
import com.google.common.collect.Collections2;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

/**
 * Utilities for organizing UC2 & UC3 logic tree branch results (hazard curves).
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class CurveUtils {

	private static final Splitter SPLIT = Splitter.on(',');
	private static final Joiner JOIN = Joiner.on(',');
	private static final Joiner TAB_JOIN = Joiner.on('\t');
	
	
	private static final Mean MEAN = new Mean();
	private static final String S = File.separator;
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static final List<String> STAT_FIELDS = Lists.newArrayList("stat",
		"2in50", "10in50", "rtgm");
	private static final List<String> SUMMARY_FIELDS = Lists.newArrayList("wt",
		"2in50", "10in50", "rtgm");
	private static final List<String> CITY_FIELDS = Lists.newArrayList("city",
		"2in50", "10in50", "rtgm");
	private static final String UC3_ROOT = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/";
	private static final String CURVE_FILE = "NSHMP08_WUS_curves";
	private static final String PARAM_FILE = "NSHMP08_WUS_params";

	public static void main(String[] args) throws IOException {

//		 generateFortranCityData();

//		 generateBranchSummaries2();
//		String srcPath = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/tree/SRP1440";
//		String locFile = srcPath + "/SRPsites.txt";
//		String curveDir = srcPath + "/reduce";
//		generateBranchSummaries2(locFile, curveDir, true);
		
//		String srcPath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		String locFile = srcPath + "sites.txt";
//		String curveDir = srcPath + "NEHRP-PBR-SRP/UC2-TimeIndep";
//		generateBranchSummaries2(locFile, curveDir, false);

//		String srcPath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		String locFile = srcPath + "sites.txt";
//		String curveDir = srcPath + "NEHRP-PBR-SRP/UC3.2-bg";
//		generateBranchSummaries2(locFile, curveDir, true);
		
//		String treePath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		File srcDir = new File(treePath + "NEHRP-PBR-SRP/UC32tree1440-bg/");
//		File outDir = new File(treePath + "NEHRP-PBR-SRP/UC32tree1440-bg-reduce");
//		File locFile = new File(treePath + "sites.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile.getPath(), false);
//		generateBranchSummaries2(locFile.getPath(), outDir.getPath(), true);

//		String treePath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		File srcDir = new File(treePath + "NEHRP-PBR-SRP/UC32treeNewSites/");
//		File outDir = new File(treePath + "NEHRP-PBR-SRP/UC32treeNewSites-reduce");
//		File locFile = new File(treePath + "newsites.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile.getPath(), false);
//		generateBranchSummaries2(locFile.getPath(), outDir.getPath(), true);

//		String treePath = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/tree/SRP1440";
//		File srcDir = new File(treePath + "/src");
//		File outDir = new File(treePath + "/reduce");
//		File locFile = new File(treePath + "/SRPsites.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile, false);

//		String treePath = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/tree/PBR1440";
//		File srcDir = new File(treePath + "/src");
//		File outDir = new File(treePath + "/reduce");
//		File locFile = new File(treePath + "/PBRsites.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile, false);

//		String treePath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		File srcDir = new File(treePath + "NEHRP-PBR-SRP/UC3.2/");
//		File locFile = new File(treePath + "/test.txt");
//		generateBranchSummaries2(locFile.getPath(), srcDir.getPath(), true);
		
		

//		String srcPath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		String locFile = srcPath + "sites.txt";
//		String curveDir = srcPath + "NEHRP-PBR-SRP/UC3.2-bg";
//		generateBranchSummaries2(locFile, curveDir, true);

//		File srcDir = new File(UC3_ROOT + "UC3.2-conv-src2");
//		File outDir = new File(UC3_ROOT + "UC3.2-conv");
//		File locFile = new File("/Users/pmpowers/projects/OpenSHA/tmp/curves/sites/all.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile.getPath(), true);
//		generateBranchSummaries2(locFile.getPath(), outDir.getPath(), false);

//		File srcDir = new File(UC3_ROOT + "UC3.2-vars-src");
//		File outDir = new File(UC3_ROOT + "UC3.2-vars");
//		File locFile = new File("/Users/pmpowers/projects/OpenSHA/tmp/curves/sites/all.txt");
//		reorganizeUC3branchResults(srcDir, outDir, locFile.getPath(), true);
//		generateBranchSummaries2(locFile.getPath(), outDir.getPath(), false);

		File srcDir = new File("/Users/pmpowers/projects/OpenSHA/tmp/hazard/PALO_VERDE/minitree");
		File outDir = new File("/Users/pmpowers/projects/OpenSHA/tmp/hazard/PALO_VERDE/minitreeOut");
		File locFile = new File("/Users/pmpowers/projects/OpenSHA/tmp/curves/sites/palo-verde.txt");
		reorganizeUC3branchResults(srcDir, outDir, locFile.getPath(), false);
		generateBranchSummaries2(locFile.getPath(), outDir.getPath(), true);

		
//		generateBranchList();

//		fix10in50s();
		
//		listAllBranches();
//		 writeConvTestMags();
		
//		String treePath = "/Users/pmpowers/projects/OpenSHA/tmp/hazard/";
//		String totSrc = treePath + "NEHRP-PBR-SRP/UC3.2/";
//		String bgSrc = treePath + "NEHRP-PBR-SRP/UC3.2-bg/";
//		String fltDest = treePath + "NEHRP-PBR-SRP/UC3.2-flt/";
//		String locFile = treePath + "sites.txt";
//		extractFaultCurves(totSrc, bgSrc, fltDest, locFile);
		
//		morganSurface();
		
	}
	
	

	/**
	 * UCERF3 logic tree reorganizer. Currently UC3 logic tree hazard curves
	 * computed for NEHRP test cities are grouped by logic tree branch; it is
	 * better to initialize a branch erf and loop location (cities) than to loop
	 * branch erf's at each location (city). This utility method groups curves
	 * for each city in a single file and writes statistical curve summaries.
	 * 
	 * @param src
	 * @param out
	 * @param ignoreWts
	 * @throws IOException
	 */
	public static void reorganizeUC3branchResults(File srcDir, File outDir,
			String locPath, boolean ignoreWts) throws IOException {
		
		// convert solutions grouped by branch to solutions grouped by city
		Set<Period> periods = EnumSet.of(GM0P00, GM1P00); //GM0P00, GM0P20, GM1P00);
		
		// create location list
		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locPath);
		Set<String> locNames = locMap.keySet();
		
		BiMap<String, Integer> indexMap = HashBiMap.create();
		Map<Integer, Double> wtMap = Maps.newHashMap();

		File[] branchDirs = srcDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File f) {
				return !f.getName().startsWith(".");
			}
		});

		for (int i=0; i<branchDirs.length; i++) {
			System.out.println(i + " " + branchDirs[i].getName());
		}

		APrioriBranchWeightProvider wtProvider = new APrioriBranchWeightProvider();
		int index = 0;
		for (File branch : branchDirs) {
			if (!branch.isDirectory()) continue;
			String branchName = branch.getName();
			indexMap.put(branchName, index);
			LogicTreeBranch ltb = LogicTreeBranch.fromFileName(branchName);
			double wt = ignoreWts ? 1.0 : wtProvider.getWeight(ltb);
			wtMap.put(index, wt);
			index++;
		}

		// normalize weights
		if (!ignoreWts) {
			Collection<Double> wts = wtMap.values();
			double sum = DataUtils.sum(Doubles.toArray(wts));
			System.out.println("Weight sum: " + sum);
			for (int idx : wtMap.keySet()) {
				double wt = wtMap.get(idx);
				wtMap.put(idx, wt / sum);
			}
		}

		Map<Period, Table<Integer, String, String>> curveMap = Maps
			.newHashMap();
		Table<Integer, String, String> table = null;
		for (Period period : periods) {
			table = ArrayTable.create(indexMap.values(), locNames);
			curveMap.put(period, table);
		}

		for (File branch : branchDirs) {
			if (!branch.isDirectory()) continue;
			String branchName = branch.getName();
			int branchIdx = indexMap.get(branchName);
			double branchWt = wtMap.get(branchIdx);
			for (Period period : periods) {
				table = curveMap.get(period);
				File periodDir = new File(branch, period.name());
				File curveFile = new File(periodDir, CURVE_FILE + ".csv");
				List<String> lines = Files.readLines(curveFile, US_ASCII);
				for (String line : Iterables.skip(lines, 1)) {
					String locStr = StringUtils.substringBefore(line, ",");
//					String curveStr = StringUtils.substringAfter(StringUtils
//						.substringAfter(
//							StringUtils.substringAfter(
//								StringUtils.substringAfter(line, ","), ","),
//							","), ",");
					String curveStr = StringUtils.substringAfter(line, ",");
					String curveStrOut = branchIdx + "," + branchWt + "," +
						curveStr;
					table.put(branchIdx, locStr, curveStrOut);
				}
			}
		}

		for (Period period : periods) {
			table = curveMap.get(period);
			for (String name : locNames) {

				File curvesOut = new File(outDir, period.name() + S +
					name + S + CURVE_FILE + ".csv");
				Files.createParentDirs(curvesOut);
				// header
				Iterable<String> gmVals = Collections2.transform(period
					.getFunction().xValues(), Functions.toStringFunction());
				List<String> headers = Lists.newArrayList("ERF#", "wt");
				Iterable<String> cityFields = Iterables.concat(headers, gmVals);
				String cityHeader = JOIN.join(cityFields) + LF;
				Files.write(cityHeader, curvesOut, US_ASCII);

				File paramsOut = new File(outDir, period.name() + S +
					name + S + PARAM_FILE + ".csv");
				Files.createParentDirs(paramsOut);
				Files.write("ERF#,BranchName" + LF, paramsOut, US_ASCII);

				// data
				Map<Integer, String> indexMapInverse = indexMap.inverse();
				for (int i = 0; i < indexMap.size(); i++) {
					String curveLine = table.get(i, name) + LF;
					String paramLine = i + "," + indexMapInverse.get(i) + LF;
					Files.append(curveLine, curvesOut, US_ASCII);
					Files.append(paramLine, paramsOut, US_ASCII);
				}
			}
		}
	}

	/**
	 * UCERF3 logic tree reorganizer. Currently UC3 logic tree hazard curves
	 * computed for NEHRP test cities are grouped by logic tree branch; it is
	 * better to initialize a branch erf and loop location (cities) than to loop
	 * branch erf's at each location (city). This utility method groups curves
	 * for each city in a single file and writes statistical curve summaries.
	 * 
	 * @param src
	 * @param out
	 * @param ignoreWts
	 * @throws IOException
	 */
	public static void reorganizeUC3branchResults(File srcDir, File outDir,
			boolean ignoreWts) throws IOException {
		// convert solutions grouped by branch to solutions grouped by city
		Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);

		BiMap<String, Integer> indexMap = HashBiMap.create();
		Map<Integer, Double> wtMap = Maps.newHashMap();

		File[] branchDirs = srcDir.listFiles();

		APrioriBranchWeightProvider wtProvider = new APrioriBranchWeightProvider();
		int index = 0;
		for (File branch : branchDirs) {
			if (!branch.isDirectory()) continue;
			String branchName = branch.getName();
			indexMap.put(branchName, index);
			LogicTreeBranch ltb = LogicTreeBranch.fromFileName(branchName);
			double wt = ignoreWts ? 1.0 : wtProvider.getWeight(ltb);
			wtMap.put(index, wt);
			index++;
		}

		// normalize weights
		if (!ignoreWts) {
			Collection<Double> wts = wtMap.values();
			double sum = DataUtils.sum(Doubles.toArray(wts));
			System.out.println("Weight sum: " + sum);
			for (int idx : wtMap.keySet()) {
				double wt = wtMap.get(idx);
				wtMap.put(idx, wt / sum);
			}
		}

		Map<Period, Table<Integer, NEHRP_TestCity, String>> curveMap = Maps
			.newHashMap();
		Table<Integer, NEHRP_TestCity, String> table = null;
		table = ArrayTable.create(indexMap.values(), NEHRP_TestCity.getCA());
		curveMap.put(GM0P00, table);
		table = ArrayTable.create(indexMap.values(), NEHRP_TestCity.getCA());
		curveMap.put(GM0P20, table);
		table = ArrayTable.create(indexMap.values(), NEHRP_TestCity.getCA());
		curveMap.put(GM1P00, table);

		for (File branch : branchDirs) {
			if (!branch.isDirectory()) continue;
			String branchName = branch.getName();
			int branchIdx = indexMap.get(branchName);
			double branchWt = wtMap.get(branchIdx);
			for (Period period : periods) {
				table = curveMap.get(period);
				File periodDir = new File(branch, period.name());
				File curveFile = new File(periodDir, CURVE_FILE + ".csv");
				List<String> lines = Files.readLines(curveFile, US_ASCII);
				for (String line : Iterables.skip(lines, 1)) {
					String cityStr = StringUtils.substringBefore(line, ",");
					NEHRP_TestCity city = NEHRP_TestCity.valueOf(cityStr);
					String curveStr = StringUtils.substringAfter(StringUtils
						.substringAfter(
							StringUtils.substringAfter(
								StringUtils.substringAfter(line, ","), ","),
							","), ",");
					String curveStrOut = branchIdx + "," + branchWt + "," +
						curveStr;
					table.put(branchIdx, city, curveStrOut);
				}
			}
		}

		for (Period period : periods) {
			table = curveMap.get(period);
			for (NEHRP_TestCity city : NEHRP_TestCity.getCA()) {

				File curvesOut = new File(outDir, period.name() + S +
					city.name() + S + CURVE_FILE + ".csv");
				Files.createParentDirs(curvesOut);
				// header
				Iterable<String> gmVals = Collections2.transform(period
					.getFunction().xValues(), Functions.toStringFunction());
				List<String> headers = Lists.newArrayList("ERF#", "wt");
				Iterable<String> cityFields = Iterables.concat(headers, gmVals);
				String cityHeader = JOIN.join(cityFields) + LF;
				Files.write(cityHeader, curvesOut, US_ASCII);

				File paramsOut = new File(outDir, period.name() + S +
					city.name() + S + PARAM_FILE + ".csv");
				Files.createParentDirs(paramsOut);
				Files.write("ERF#,BranchName" + LF, paramsOut, US_ASCII);

				// data
				Map<Integer, String> indexMapInverse = indexMap.inverse();
				for (int i = 0; i < indexMap.size(); i++) {
					String curveLine = table.get(i, city) + LF;
					String paramLine = i + "," + indexMapInverse.get(i) + LF;
					Files.append(curveLine, curvesOut, US_ASCII);
					Files.append(paramLine, paramsOut, US_ASCII);
				}
			}
		}
	}

	/**
	 * Utility method to create summaries of logic tree branch hazard curves.
	 */
	public static void generateBranchSummaries() {
		Iterable<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);
		Iterable<NEHRP_TestCity> cities = NEHRP_TestCity.getCA(); // EnumSet.of(VENTURA);
		String imrID = NSHMP08_WUS.SHORT_NAME;
//		 String dir = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC2-TimeIndep";
//		String dir = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/convABM";
		String dir = "/Users/pmpowers/Documents/OpenSHA/RTGM/data/UC3/treeFM32single";
		try {
			// boolean is tornado
			runBranchSummaries(dir, imrID, periods, cities, true);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static void generateBranchSummaries2(String locPath, String curveDir,
			boolean tornado) throws IOException {
		
		Iterable<Period> periods = EnumSet.of(GM0P00, GM1P00); //GM0P00, GM0P20, GM1P00);
		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locPath);
		String imrID = NSHMP08_WUS.SHORT_NAME;
		try {
			// boolean is tornado
			runBranchSummaries2(curveDir, imrID, periods, locMap.keySet(), tornado);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	/*
	 */
	private static void runBranchSummaries2(String dir, String imrID,
			Iterable<Period> periods, Iterable<String> locNames,
			boolean tornado) throws IOException {
		for (Period p : periods) {
			for (String name : locNames) {
				System.out.println(name + " " + p);
				File srcDir = new File(dir + S + p.name() + S + name);
				File srcFile = new File(srcDir, imrID + "_curves.csv");
				File branchFile = new File(srcDir, imrID + "_params.csv");
				File statFile = new File(srcDir, imrID + "_stats.csv");
				File sumFile = new File(srcDir, imrID + "_summary.csv");
				File torRTGM_File = null, tor2in50_File = null, tor10in50_File = null;
				if (tornado) {
					torRTGM_File = new File(srcDir, imrID + "_tornado_rtgm.csv");
					tor2in50_File = new File(srcDir, imrID + "_tornado_2in50.csv");
					tor10in50_File = new File(srcDir, imrID + "_tornado_10in50.csv");
				}
				summarizeBranches(srcFile, branchFile, statFile, sumFile,
					torRTGM_File, tor2in50_File, tor10in50_File, p);
			}
		}
	}


	/*
	 * Create summaries of logic tree branch hazard curves.
	 */
	private static void runBranchSummaries(String dir, String imrID,
			Iterable<Period> periods, Iterable<NEHRP_TestCity> cities,
			boolean tornado) throws IOException {
		for (Period p : periods) {
			for (NEHRP_TestCity c : cities) {
				File srcDir = new File(dir + S + p.name() + S + c.name());
				File srcFile = new File(srcDir, imrID + "_curves.csv");
				File branchFile = new File(srcDir, imrID + "_params.csv");
				File statFile = new File(srcDir, imrID + "_stats.csv");
				File sumFile = new File(srcDir, imrID + "_summary.csv");
				File torRTGM_File = null, tor2in50_File = null, tor10in50_File = null;
				if (tornado) {
					torRTGM_File = new File(srcDir, imrID + "_tornado_rtgm.csv");
					tor2in50_File = new File(srcDir, imrID + "_tornado_2in50.csv");
					tor10in50_File = new File(srcDir, imrID + "_tornado_10in50.csv");
				}
				summarizeBranches(srcFile, branchFile, statFile, sumFile,
					torRTGM_File, tor2in50_File, tor10in50_File, p);
			}
		}
	}

	/*
	 * Reads an erf branch file that has columns with wts, rtgm values and a
	 * hazard curve values and outputs one summary file with min, max, and mean
	 * data and another with just wt, pe2in50, and pe10in50.
	 */
	private static void summarizeBranches(File curveFile, File branchFile,
			File stat, File summary, File torRTGM, File tor2in50,
			File tor10in50, Period period)
					throws IOException {
		List<String> curveLines = Files.readLines(curveFile, US_ASCII);
		List<Double> weights = Lists.newArrayList();
		List<Double> rtgms = Lists.newArrayList();
		List<Double> pe2in50s = Lists.newArrayList();
		List<Double> pe10in50s = Lists.newArrayList();
		XY_DataSetList curves = new XY_DataSetList();
		
		// branch list and index reverse lookup 
		List<String> branchLines = Files.readLines(branchFile, US_ASCII);
		List<String> branchList = Lists.newArrayList();

		// create model function: first line has ERF#, wt, rtgm, gm-vals ...
		Iterable<String> firstLine = SPLIT.split(curveLines.get(0));
		DiscretizedFunc curveModel = new ArbitrarilyDiscretizedFunc();
		for (String gmStr : Iterables.skip(firstLine, 2)) {
			double gmVal = Double.parseDouble(gmStr);
			curveModel.set(gmVal, 0.0);
		}

		// fill curves, pe intercepts and rtgm lists
		for (String line : Iterables.skip(curveLines, 1)) {
			Iterable<String> vals = SPLIT.split(line);

			double weight = Double.parseDouble(Iterables.get(vals, 1));
			weights.add(weight);

			DiscretizedFunc curve = curveModel.deepClone();
			int idx = 0;
			for (String val : Iterables.skip(vals, 2)) {
				double annRate = Double.parseDouble(val);
				if (annRate < 0) annRate = 0;
				curve.set(idx++, annRate);
			}
			curves.add(curve);

			double pe2in50 = getPE(curve, PE2IN50);
			pe2in50s.add(pe2in50);

			double pe10in50 = getPE(curve, PE10IN50);
			pe10in50s.add(pe10in50);

			double rtgm = getRTGM(curve, period);
			rtgms.add(rtgm);
		}
		
		// fill branch list and index lookup map
		for (String line : Iterables.skip(branchLines, 1)) {
			Iterable<String> vals = SPLIT.split(line);
			String branchName = Iterables.get(vals, 1);
			branchList.add(branchName);
		}

		// write PEs
		Iterable<String> summaryFields = Iterables.concat(SUMMARY_FIELDS);
		String summaryHeader = JOIN.join(summaryFields) + LF;
		Files.write(summaryHeader, summary, US_ASCII);
		for (int i = 0; i < weights.size(); i++) {
			Iterable<String> lineData = createSummaryData(weights.get(i),
				pe2in50s.get(i), pe10in50s.get(i), rtgms.get(i));
			String summaryLine = JOIN.join(lineData) + LF;
			Files.append(summaryLine, summary, US_ASCII);
		}
		
		// write tornado data
		if (torRTGM != null) {
			TornadoData tdRTGM = new UC3_TornadoBuilder(branchList, rtgms).build();
			Files.write(tdRTGM.toSortedString(), torRTGM, US_ASCII);
			TornadoData td2in50 = new UC3_TornadoBuilder(branchList, pe2in50s).build();
			Files.write(td2in50.toSortedString(), tor2in50, US_ASCII);
			TornadoData td10in50 = new UC3_TornadoBuilder(branchList, pe10in50s).build();
			Files.write(td10in50.toSortedString(), tor10in50, US_ASCII);
		}
		
		// calc and write stats
		FractileCurveCalculator fcc = new FractileCurveCalculator(curves,
			weights);
		XY_DataSet minCurve = fcc.getMinimumCurve();
		XY_DataSet maxCurve = fcc.getMaximumCurve();
		XY_DataSet meanCurve = fcc.getMeanCurve();

		ArbitrarilyDiscretizedFunc f;

		f = xyDataToFunc(meanCurve);
		double mean2in50 = getPE(f, PE2IN50);
		double mean10in50 = getPE(f, PE10IN50);
		double meanRTGM = getRTGM(f, period);

		f = xyDataToFunc(minCurve);
		double min2in50 = getPE(f, PE2IN50);
		double min10in50 = getPE(f, PE10IN50);
		double minRTGM = getRTGM(f, period);

		f = xyDataToFunc(maxCurve);
		double max2in50 = getPE(f, PE2IN50);
		double max10in50 = getPE(f, PE10IN50);
		double maxRTGM = getRTGM(f, period);

		// header
		Iterable<String> statFields = Iterables.concat(STAT_FIELDS,
			Iterables.skip(firstLine, 2));
		String statHeader = JOIN.join(statFields) + LF;
		Files.write(statHeader, stat, US_ASCII);

		// mean
		Iterable<String> meanDat = createData("mean", mean2in50, mean10in50,
			meanRTGM, meanCurve);
		String meanLine = JOIN.join(meanDat) + LF;
		Files.append(meanLine, stat, US_ASCII);

		// min
		Iterable<String> minDat = createData("min", min2in50, min10in50,
			minRTGM, minCurve);
		String minLine = JOIN.join(minDat) + LF;
		Files.append(minLine, stat, US_ASCII);

		// max
		Iterable<String> maxDat = createData("max", max2in50, max10in50,
			maxRTGM, maxCurve);
		String maxLine = JOIN.join(maxDat) + LF;
		Files.append(maxLine, stat, US_ASCII);

	}

	private static double getPE(DiscretizedFunc f, ProbOfExceed pe) {
		return f.getFirstInterpolatedX_inLogXLogYDomain(pe.annualRate());
	}

	private static double getRTGM(DiscretizedFunc f, Period p) {
		if (!(p == GM0P20 || p == GM1P00)) return 0;
		Frequency freq = p.equals(GM0P20) ? SA_0P20 : SA_1P00;
		RTGM rtgm = RTGM.create(f, freq, 0.8).call();
		return rtgm.get();
	}

	private static Iterable<String> createData(String label, double pe2in50,
			double pe10in50, double rtgm, XY_DataSet curve) {

		Iterable<String> intercepts = Lists.newArrayList(label,
			Double.toString(pe2in50), Double.toString(pe10in50),
			Double.toString(rtgm));

		Iterable<String> values = Collections2.transform(curve.yValues(),
			Functions.toStringFunction());

		return Iterables.concat(intercepts, values);
	}

	private static Iterable<String> createSummaryData(double wt,
			double pe2in50, double pe10in50, double rtgm) {
		Iterable<String> values = Lists.newArrayList(Double.toString(wt),
			Double.toString(pe2in50), Double.toString(pe10in50),
			Double.toString(rtgm));
		return values;
	}

	/*
	 * Mine fortran results for curves and exceedance ground motions at
	 * specified NEHRP test cities.
	 */

	public static void generateFortranCityData() throws IOException {
		String srcPath = "/Volumes/Scratch/nshmp-sources/FortranUpdate";
		String outPath = "/Volumes/Scratch/rtgm/FortranUpdate2";

		Iterable<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);
		Map<String, Location> siteMap = 
				UC3_CalcUtils.readSiteFile("tmp/UC3sites/all.txt");

		for (Period p : periods) {
			File src = new File(srcPath + S + p.name() + S + "curves.dat");
			File out = new File(outPath + S + p.name() + S + "FORT_curves.csv");
			try {
				Files.createParentDirs(out);
				deriveCityData(siteMap, src, out, p);
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
	}

	private static void deriveCityData(Map<String, Location> siteMap,
			File src, File out, Period period) throws IOException {

		// header
		Iterable<String> gmVals = Collections2.transform(period.getFunction()
			.xValues(), Functions.toStringFunction());
		Iterable<String> cityFields = Iterables.concat(CITY_FIELDS, gmVals);
		String cityHeader = JOIN.join(cityFields) + LF;
		Files.write(cityHeader, out, US_ASCII);

		// extract data
		CurveContainer cc = CurveContainer.create(src);
		for (String siteName : siteMap.keySet()) {
			Location loc = siteMap.get(siteName);
			DiscretizedFunc curve = cc.getCurve(loc);
			double pe2in50 = getPE(curve, PE2IN50);
			double pe10in50 = getPE(curve, PE10IN50);
			double rtgm = getRTGM(curve, period);
			Iterable<String> cityData = createData(siteName, pe2in50,
				pe10in50, rtgm, curve);
			String cityLine = JOIN.join(cityData) + LF;
			Files.append(cityLine, out, US_ASCII);
		}
	}
	

	static class UC3_TornadoBuilder {
		
		private List<Double> values;
		private List<String> branchNames;
		private Map<String, Integer> branchIdxMap = Maps.newHashMap();
		
		UC3_TornadoBuilder(List<String> branchNames, List<Double> values) {
			checkArgument(branchNames.size() == values.size());
			this.branchNames = branchNames;
			this.values = values;
			int idx = 0;
			for (String name : branchNames) {
				branchIdxMap.put(name, idx++);
			}
		}
		
		public TornadoData build() {			

			List<Class> classList = Lists.newArrayList(); // for ordering
			Map<Class, Set<LogicTreeBranchNode>> nodeMap = Maps.newHashMap();

			// init class to node variant map
			for (Class clazz : LogicTreeBranch.getLogicTreeNodeClasses()) {
				Set<LogicTreeBranchNode> nodeSet = Sets.newHashSet();
				classList.add(clazz);
				nodeMap.put(clazz, nodeSet);
			}

			// fill nodeMap wtih valid logic tree variants
			for (String name : branchNames) {
				LogicTreeBranch ltb = LogicTreeBranch.fromFileName(name);
				for (LogicTreeBranchNode node : ltb) {
					Class clazz = LogicTreeBranch.getEnumEnclosingClass(node.getClass());
					Set<LogicTreeBranchNode> nodeSet = nodeMap.get(clazz);
					nodeSet.add(node);
				}
			}

			// median value and logic tree branch
			int medIdx = medianIndex(values);
			double medVal = values.get(medIdx);
			String medBrName = branchNames.get(medIdx);
			LogicTreeBranch medLTB = LogicTreeBranch.fromFileName(medBrName);
			
			// loop all valid nodes gathering values for branch variants
			TornadoData td = new TornadoData(medVal);
			for (Class clazz : classList) {
				Set<LogicTreeBranchNode> nodeSet = nodeMap.get(clazz);
				for (LogicTreeBranchNode node : nodeSet) {
					LogicTreeBranch ltb = (LogicTreeBranch) medLTB.clone();
					ltb.setValue(node);
					String brName = ltb.buildFileName();
					int brIdx = branchIdxMap.get(brName);
					double brVal = values.get(brIdx);
					td.add(clazz, (Enum) node, brVal);
				}
			}

			return td;
		}
		
	}
	
	/*
	 * If values.size() is odd, method return the index of the median value. If
	 * values.size() is even the index of values.size()/2 is returned.
	 */
	private static int medianIndex(List<Double> values) {
		double[] sortedVals = Doubles.toArray(values);
		Arrays.sort(sortedVals);
		int idx = (sortedVals.length - 1) / 2;
		double median = sortedVals[idx];
		return values.indexOf(median);
	}
	
	
	
	private static final String fix10in50path = "/Users/pmpowers/Documents/OpenSHA/RTGM/data";
	private static final String[] fixList = {
		"FortranUpdate", "FSS_UC2map", "MeanUCERF2", "MeanUCERF2update",
		"MeanUCERF2update_FM2P1", "ModMeanUCERF2update_FM2P1", "NSHMP_CA_SHA", 
		"NSHMP_CA_SHA-epi", "NSHMP_SHA", "NSHMP_SHA-epi"
	};
	
	private static void fix10in50s() {
		try {
			Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);
			for (String dir : fixList) {
				for (Period p : periods) {
					String path = fix10in50path + S + dir + S + p + S;
					File curveFile = new File(path + "NSHMP08_WUS_curves.csv");
					if (!curveFile.exists()) {
						curveFile = new File(path + "NSHMP08_curves.csv");
					}
					if (!curveFile.exists()) {
						curveFile = new File(path + "FORT_curves.csv");
					}
					
					List<String> linesIn = Files.readLines(curveFile, US_ASCII);
					Iterable<String> firstLine = SPLIT.split(linesIn.get(0));
					DiscretizedFunc curveModel = new ArbitrarilyDiscretizedFunc();
					for (String gmStr : Iterables.skip(firstLine, 4)) {
						double gmVal = Double.parseDouble(gmStr);
						curveModel.set(gmVal, 0.0);
					}
					Files.write(linesIn.get(0) + LF, curveFile, US_ASCII);
	
					// fill curves, pe intercepts and rtgm lists
					for (String line : Iterables.skip(linesIn, 1)) {
						Iterable<String> lineIter = SPLIT.split(line);
						List<String> lineList = Lists.newArrayList(lineIter);
						
						DiscretizedFunc curve = curveModel.deepClone();
						int idx = 0;
						for (String val : Iterables.skip(lineIter, 4)) {
							double annRate = Double.parseDouble(val);
							curve.set(idx++, annRate);
						}
						double pe10in50 = getPE(curve, PE10IN50);
						lineList.set(2, Double.toString(pe10in50));
						String fixedLine = JOIN.join(lineList) + LF;
						Files.append(fixedLine, curveFile, US_ASCII);
					}
				}
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	private static void listAllBranches() {
		try {
			String path1440 = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/tree/2012_10_29-tree-fm31_x7-fm32_x1_COMPOUND_SOL.zip";
			CompoundFaultSystemSolution cfss = UC3_CalcUtils
				.getCompoundSolution(path1440);
			List<LogicTreeBranch> branches = Lists.newArrayList(cfss
				.getBranches());
			File out = new File("tmp/branchlist1440.txt");
			Files.write("", out, US_ASCII);
			int idx = 0;

			for (LogicTreeBranch branch : branches) {
				Files.append((idx++) + " " + branch.buildFileName() + LF, out,
					US_ASCII);
			}

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	
	private static void writeConvTestMags() throws IOException {
		// rupture of interest: 29749
		int fssRupIdx = 29749;
		int maxIdx = 1;
		
		File out = new File("tmp/SRPconvTestRupRates.txt");
		String header = TAB_JOIN.join("FSSidx","fssRate","fssMag","erfRate","erfMag") + LF;
		Files.write(header, out, US_ASCII);
		
		// load conv fss
		String convSolPath = "/Users/pmpowers/projects/OpenSHA/tmp/invSols/conv/FM3_1_ZENG_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip";
		AverageFaultSystemSolution afss = UC3_CalcUtils.getAvgSolution(convSolPath);
		for (int i=0; i<maxIdx; i++) {
			InversionFaultSystemSolution fss = afss.getSolution(i);
			double fssRupRate = fss.getRateForRup(fssRupIdx);
			double fssRupMag = fss.getRupSet().getMagForRup(fssRupIdx);
			System.out.println(fssRupRate + "\t" + fssRupMag);
		
//			UCERF3_FaultSysSol_ERF erf = new UCERF3_FaultSysSol_ERF(fss);
//			UC3_CalcUtils.initUC3(erf, IncludeBackgroundOption.EXCLUDE, false, true, 1.0);
			FaultSystemSolutionERF erf = UC3_CalcUtils.getUC3_ERF(convSolPath, IncludeBackgroundOption.EXCLUDE, false, true, 1.0);
			erf.updateForecast();
			int srcIdx = -1;
			for (int j=0; j<erf.getNumFaultSystemSources(); j++) {
				int rupIdx = erf.getFltSysRupIndexForSource(j);
				if (rupIdx == fssRupIdx) {
					srcIdx = j;
					break;
				}
			}
			checkArgument(srcIdx != -1);
			ProbEqkSource src = erf.getSource(srcIdx);
			System.out.println(src.getSourceMetadata());
			ProbEqkRupture rup = src.getRupture(0);
			double erfRupRate = rup.getMeanAnnualRate(1d);
			double erfRupMag = rup.getMag();
			String outLine = TAB_JOIN.join(i, fssRupRate, fssRupMag, erfRupRate, erfRupMag) + LF;
			System.out.println(outLine);
			Files.append(outLine, out, US_ASCII);
		}
	}
	
	private static ArbitrarilyDiscretizedFunc xyDataToFunc(XY_DataSet xy) {
		ArbitrarilyDiscretizedFunc f = new ArbitrarilyDiscretizedFunc();
		for (Point2D p : xy) {
			f.set(p);
		}
		return f;
	}
	
	/*
	 * subtracts background from total curves to get fault source only curves
	 */
	private static void extractFaultCurves(String totalDir, String bgDir,
			String fltDir, String locFile) throws IOException {
		
		Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00);
		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locFile);
		String paramFile = PARAM_FILE + ".csv";
		String curveFile = CURVE_FILE + ".csv";
		
		for (Period p : periods) {
			for (String locName : locMap.keySet()) {
					
				File totSrc = new File(totalDir + p + S + locName + S + curveFile);
				File bgSrc = new File(bgDir + p + S + locName + S + curveFile);
				File fltDest = new File(fltDir + p + S + locName + S + curveFile);

				Iterator<String> totCurves = Files.readLines(totSrc, US_ASCII).iterator();
				Iterator<String> bgCurves = Files.readLines(bgSrc, US_ASCII).iterator();
				
				// header
				Files.createParentDirs(fltDest);
				Files.write(totCurves.next() + LF, fltDest, US_ASCII);
				bgCurves.next();
				
				// curves
				while (totCurves.hasNext()) {
					Iterable<String> totDat = SPLIT.split(totCurves.next());
					Iterable<String> bgDat = SPLIT.split(bgCurves.next());					
					Iterable<String> totValStr = Iterables.skip(totDat, 2);
					Iterable<String> bgValStr = Iterables.skip(bgDat, 2);
					List<Double> totVals = stringToDouble(totValStr);
					List<Double> bgVals = stringToDouble(bgValStr);
					List<Double> fltVals = DataUtils.subtract(totVals, bgVals);
					Iterable<?> dat = Iterables.concat(Iterables.limit(totDat, 2), fltVals);
					String fltStr = JOIN.join(dat);
					Files.append(fltStr + LF, fltDest, US_ASCII);
				}
				
				// copy params / logic tree branch list
				File paramSrc = new File(totalDir + p + S + locName + S + paramFile);
				File paramDest = new File(fltDir + p + S + locName + S + paramFile);
				Files.createParentDirs(paramDest);
				Files.copy(paramSrc, paramDest);
			}
		}
	}
	
	private static List<Double> stringToDouble(Iterable<String> strVals) {
		List<Double> vals = Lists.newArrayList();
		for (String s : strVals) {
			vals.add(Double.parseDouble(s));
		}
		return vals;
	}
	
	private static List<String> doubleToString(Iterable<Double> dblVals) {
		List<String> strs = Lists.newArrayList();
		for (Double d : dblVals) {
			strs.add(Double.toString(d));
		}
		return strs;
	}
	
	
	private static void morganSurface() throws IOException {
		// 20 km fault width @ 50 deg = 15.321 km depth
		File srcDat = new File("tmp/morgan/segment.txt");
		File outDat = new File("tmp/morgan/segment-surface.txt");
		Iterable<String> locLines = Files.readLines(srcDat, US_ASCII);
		FaultTrace trace = new FaultTrace("test");
		Splitter split =  Splitter.on(" ").omitEmptyStrings();
		for (String locLine : locLines) {
			Iterator<String> it = split.split(locLine).iterator();
			double lon = Double.parseDouble(it.next());
			double lat = Double.parseDouble(it.next());
			trace.add(new Location(lat, lon));
		}
		SimpleFaultData sfd = new SimpleFaultData(50, 15.321, 0, trace);
		double spacing = trace.getTraceLength()/255;
		StirlingGriddedSurface sgs = new StirlingGriddedSurface(sfd, spacing, 0.158); // CW_EB 256 x 128
		System.out.println(sgs);
		System.out.println(sgs.getGridSpacingDownDip());
		System.out.println(sgs.getGridSpacingAlongStrike());
		
		Files.write(sgs.toString(), outDat, US_ASCII);
		Files.append("Along strike spacing: " + sgs.getGridSpacingAlongStrike() + "\n", outDat, US_ASCII);
		Files.append("Down dip spacing: " + sgs.getGridSpacingDownDip() + "\n", outDat, US_ASCII);

	}
	
}
