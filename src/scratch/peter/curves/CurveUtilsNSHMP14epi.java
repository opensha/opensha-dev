package scratch.peter.curves;

import static com.google.common.base.Charsets.US_ASCII;
import static com.google.common.base.Preconditions.checkArgument;
import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.nshmp2.util.Period.GM4P00;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.ABM;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.GEOLOGIC;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.NEOKINEMA;
import static scratch.UCERF3.enumTreeBranches.DeformationModels.ZENGBB;
import static scratch.UCERF3.enumTreeBranches.FaultModels.FM3_1;
import static scratch.UCERF3.enumTreeBranches.FaultModels.FM3_2;
import static scratch.UCERF3.enumTreeBranches.MaxMagOffFault.MAG_7p3;
import static scratch.UCERF3.enumTreeBranches.MaxMagOffFault.MAG_7p6;
import static scratch.UCERF3.enumTreeBranches.MaxMagOffFault.MAG_7p9;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.ELLB_SQRT_LENGTH;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.ELLSWORTH_B;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.HANKS_BAKUN_08;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.SHAW_2009_MOD;
import static scratch.UCERF3.enumTreeBranches.ScalingRelationships.SHAW_CONST_STRESS_DROP;
import static scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels.TAPERED;
import static scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels.UNIFORM;
import static scratch.UCERF3.enumTreeBranches.SpatialSeisPDF.UCERF2;
import static scratch.UCERF3.enumTreeBranches.SpatialSeisPDF.UCERF3;
import static scratch.UCERF3.enumTreeBranches.TotalMag5Rate.RATE_6p5;
import static scratch.UCERF3.enumTreeBranches.TotalMag5Rate.RATE_7p9;
import static scratch.UCERF3.enumTreeBranches.TotalMag5Rate.RATE_9p6;
import static scratch.peter.curves.ProbOfExceed.PE10IN50;
import static scratch.peter.curves.ProbOfExceed.PE1IN100;
import static scratch.peter.curves.ProbOfExceed.PE2IN50;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.exceptions.InvalidRangeException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.Interpolate;
import org.opensha.nshmp2.util.Period;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.peter.curves.CurveUtilsUC33.UC3_TornadoBuilder;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Enums;
import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.ArrayTable;
import com.google.common.collect.Collections2;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

@SuppressWarnings("javadoc")
public class CurveUtilsNSHMP14epi {

	private static final Joiner JOIN = Joiner.on(',');
	private static final Splitter SPLIT = Splitter.on(',').trimResults();
	private static final String S = File.separator;
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static final String CURVE_FILE = "curves.csv";
	private static final String PARAM_FILE = "params.csv";
	private static final String ROOT = "tmp/UC33/";
	private static final Map<String, Double> GMM_ID_WT_MAP;
	private static final Map<String, Double> EPI_ID_WT_MAP;

	static {
		GMM_ID_WT_MAP = Maps.newLinkedHashMap();
		GMM_ID_WT_MAP.put("AS", 0.22);
		GMM_ID_WT_MAP.put("BS", 0.22);
		GMM_ID_WT_MAP.put("CB", 0.22);
		GMM_ID_WT_MAP.put("CY", 0.22);
		GMM_ID_WT_MAP.put("ID", 0.12);

		EPI_ID_WT_MAP = Maps.newLinkedHashMap();
		EPI_ID_WT_MAP.put("HI", 0.185);
		EPI_ID_WT_MAP.put("OFF", 0.630);
		EPI_ID_WT_MAP.put("LO", 0.185);
	}

	public static void main(String[] args) throws IOException {
		Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00, GM4P00);

		String srcPath = ROOT + "curves/src/tree-14-gmm-epi/";
		String outPath = ROOT + "curves/src-reduce/NSHMP14-GMM-EPI/";
		String locsPath = ROOT + "curvejobs/sites/noPBR.txt";
		String flagsPath = outPath + "nodeFlags.csv";

//		 reorganize(periods, srcPath, outPath, locsPath);

//		 if (!new File(flagsPath).exists()) {
//		 createFlagFile(outPath + "GM0P00/LOS_ANGELES/params.csv", flagsPath);
//		 }

//		 runBranchSummaries2(outPath, flagsPath, periods, locsPath, PE2IN50);

		buildTornados(outPath, periods, locsPath, PE2IN50);
	}

	private static void runBranchSummaries2(String dir, String flagsPath, Iterable<Period> periods,
			String locsPath, ProbOfExceed pe) throws IOException {

		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locsPath);
		Iterable<String> locNames = locMap.keySet();

		File flagFile = new File(flagsPath);

		for (Period p : periods) {
			System.out.println(p);
			for (String name : locNames) {
				System.out.println("    " + name);
				String srcPath = dir + S + p.name() + S + name;
				summarizePE(srcPath, flagFile, pe);
			}
		}
	}

	private static void summarizePE(String srcDir, File flagFile, ProbOfExceed pe)
			throws IOException {

		File curveFile = new File(srcDir, "curves.csv");

		XY_DataSetList curves = initCurves(curveFile);
		List<Double> weights = initWeights(curveFile);

		// derive and write fractile curves - file will be rewritten for each pe
		Map<Fractile, XY_DataSet> fractileMap = computeFractiles(curves, weights);
		File fracFile = new File(srcDir, "fractiles.csv");
		writeFractiles(fractileMap, fracFile);

		// compute true weighted mean
		double tMean = getPE(XYtoFunc(fractileMap.get(Fractile.MEAN)), pe);

		// compute ground motion at PE and te PE (or rate) at the mean gm
		List<Double> gms = Lists.newArrayList();
		List<Double> rates = Lists.newArrayList();
		for (XY_DataSet curve : curves) {
			DiscretizedFunc f = (DiscretizedFunc) curve;
			gms.add(getPE(f, pe));
			rates.add(getRate(f, tMean));
		}

		// write distribution data
		File distrosFile = new File(srcDir, "distros_" + pe.name() + ".csv");
		writeDistros(gms, rates, weights, distrosFile);

		// compute and write mean for each branch choice and tMean with branch
		// removed
		Map<String, List<Boolean>> flagMap = readFlags(flagFile);
		Map<String, double[]> leafDataMap = createLeafDataMap(curves, weights, flagMap);
		File leafFile = new File(srcDir, "leaves_" + pe.name() + ".csv");
		writeLeafData(leafDataMap, tMean, leafFile);

	}

	private static XY_DataSetList initCurves(File curveFile) throws IOException {
		List<String> curveLines = Files.readLines(curveFile, US_ASCII);
		XY_DataSetList curves = new XY_DataSetList();

		// create model function: first line has BranchIndex, BranchWt, gm-vals
		// ...
		Iterable<String> firstLine = SPLIT.split(curveLines.get(0));
		DiscretizedFunc curveModel = new ArbitrarilyDiscretizedFunc();
		for (String gmStr : Iterables.skip(firstLine, 2)) {
			double gmVal = Double.parseDouble(gmStr);
			curveModel.set(gmVal, 0.0);
		}

		// fill curves
		for (String line : Iterables.skip(curveLines, 1)) {
			DiscretizedFunc curve = curveModel.deepClone();
			int idx = 0;
			for (String val : Iterables.skip(SPLIT.split(line), 2)) {
				double annRate = Double.parseDouble(val);
				if (annRate < 0) annRate = 0;
				curve.set(idx++, annRate);
			}
			curves.add(curve);
		}
		return curves;
	}

	private static List<Double> initWeights(File curveFile) throws IOException {
		List<Double> weights = Lists.newArrayList();
		for (String line : Iterables.skip(Files.readLines(curveFile, US_ASCII), 1)) {
			weights.add(FluentIterable.from(SPLIT.split(line)).transform(DoubleValueOfFn.INSTANCE)
				.get(1));
		}
		return weights;
	}

	private static Map<Fractile, XY_DataSet> computeFractiles(XY_DataSetList curves,
			List<Double> weights) {
		Map<Fractile, XY_DataSet> fMap = Maps.newEnumMap(Fractile.class);
		FractileCurveCalculator fcc = new FractileCurveCalculator(curves, weights);
		fMap.put(Fractile.MIN, fcc.getMinimumCurve());
		fMap.put(Fractile.F02, fcc.getFractile(0.02));
		fMap.put(Fractile.F16, fcc.getFractile(0.16));
		fMap.put(Fractile.MEAN, fcc.getMeanCurve());
		fMap.put(Fractile.F84, fcc.getFractile(0.84));
		fMap.put(Fractile.F98, fcc.getFractile(0.98));
		fMap.put(Fractile.MAX, fcc.getMaximumCurve());
		return fMap;
	}

	private static void writeFractiles(Map<Fractile, XY_DataSet> fracMap, File fracFile)
			throws IOException {
		// header
		XY_DataSet model = fracMap.get(Fractile.MIN);
		Iterable<String> fracFields = Iterables.concat(
			Lists.newArrayList("fractile", "IML_At2in50", "RateAtMeanIML"),
			Iterables.transform(model.xValues(), Functions.toStringFunction()));
		Files.write(JOIN.join(fracFields) + LF, fracFile, US_ASCII);
		
		// true mean which is used to get fractile rates
		double tMean = getPE(XYtoFunc(fracMap.get(Fractile.MEAN)), PE2IN50);

		// body
		for (Entry<Fractile, XY_DataSet> fracEntry : fracMap.entrySet()) {

			ArbitrarilyDiscretizedFunc f = XYtoFunc(fracEntry.getValue());
			double IMLat2in50 = getPE(f, PE2IN50);
			double rateAtMeanIML = getRate(f, tMean);

			Iterable<String> fracData = Iterables.concat(
				Lists.newArrayList(fracEntry.getKey().name()),
				Iterables.transform(Lists.newArrayList(IMLat2in50, rateAtMeanIML),
					Functions.toStringFunction()),
				Iterables.transform(fracEntry.getValue().yValues(), Functions.toStringFunction()));
			Files.append(JOIN.join(fracData) + LF, fracFile, US_ASCII);
		}
	}

	private static void writeDistros(List<Double> gms, List<Double> rates, List<Double> weights,
			File distroFile) throws IOException {
		Files.write("wt,gm,rate" + LF, distroFile, US_ASCII);
		for (int i = 0; i < gms.size(); i++) {
			Files.append(weights.get(i) + "," + gms.get(i) + "," + rates.get(i) + LF, distroFile,
				US_ASCII);
		}
	}

	private static void writeLeafData(Map<String, double[]> leafMap, double mean, File leafFile)
			throws IOException {
		Files.write("leaf,mean,noMean" + LF, leafFile, US_ASCII);
		Files.append("NONE," + mean + "," + mean + LF, leafFile, US_ASCII);
		for (Entry<String, double[]> leafEntry : leafMap.entrySet()) {
			Files.append(
				leafEntry.getKey() + "," + leafEntry.getValue()[0] + "," + leafEntry.getValue()[1] +
					LF, leafFile, US_ASCII);
		}
	}

	private static Map<String, double[]> createLeafDataMap(XY_DataSetList curves,
			List<Double> weights, Map<String, List<Boolean>> flagMap) {
		Map<String, double[]> leafMap = Maps.newLinkedHashMap();
		for (Entry<String, List<Boolean>> leafEntry : flagMap.entrySet()) {
			double leafMean = meanForLeaf(curves, weights, leafEntry.getValue());
			double noLeafMean = meanWithoutLeaf(curves, weights, leafEntry.getValue());
			leafMap.put(leafEntry.getKey(), new double[] { leafMean, noLeafMean });
		}
		return leafMap;
	}

	private static double meanForLeaf(XY_DataSetList curves, List<Double> weights,
			List<Boolean> flags) {
		XY_DataSetList leafCurves = new XY_DataSetList();
		List<Double> leafWeights = Lists.newArrayList();
		for (int i = 0; i < flags.size(); i++) {
			if (!flags.get(i)) continue;
			leafCurves.add(curves.get(i));
			leafWeights.add(weights.get(i));
		}
		FractileCurveCalculator fcc = new FractileCurveCalculator(leafCurves, leafWeights);
		return getPE(XYtoFunc(fcc.getMeanCurve()), PE2IN50);
	}

	private static double meanWithoutLeaf(XY_DataSetList curves, List<Double> weights,
			List<Boolean> flags) {
		XY_DataSetList otherCurves = new XY_DataSetList();
		List<Double> otherWeights = Lists.newArrayList();
		for (int i = 0; i < flags.size(); i++) {
			if (flags.get(i)) continue;
			otherCurves.add(curves.get(i));
			otherWeights.add(weights.get(i));
		}
		FractileCurveCalculator fcc = new FractileCurveCalculator(otherCurves, otherWeights);
		return getPE(XYtoFunc(fcc.getMeanCurve()), PE2IN50);
	}

	private static void reorganize(Set<Period> periods, String srcPath, String outPath,
			String locPath) throws IOException {

		// init reference lists and output files
		List<String> branchIDs = Lists.newArrayList();
		List<Double> branchWTs = Lists.newArrayList();
		initBranchData(srcPath, branchIDs, branchWTs);
		Set<String> locIDs = UC3_CalcUtils.readSiteFile(locPath).keySet();
		initOutputDirsAndFiles(outPath, periods, locIDs);

		for (int i = 0; i < branchIDs.size(); i++) {
			String branchID = branchIDs.get(i);
			double branchWT = branchWTs.get(i);

			for (Period p : periods) {

				// base index per branch
				int branchIndex = i * GMM_ID_WT_MAP.size() * EPI_ID_WT_MAP.size();

				int gmmCount = 0;
				for (Entry<String, Double> gmmEntry : GMM_ID_WT_MAP.entrySet()) {
					String gmmID = gmmEntry.getKey();
					double gmmWT = gmmEntry.getValue();

					int gmmIndex = branchIndex + gmmCount * EPI_ID_WT_MAP.size();

					int epiCount = 0;
					for (Entry<String, Double> epiEntry : EPI_ID_WT_MAP.entrySet()) {
						String epiID = epiEntry.getKey();
						double epiWT = epiEntry.getValue();

						String gmmEpiDir = gmmID + "-" + epiID;
						File curveFile = new File(srcPath + gmmEpiDir + S + branchID + S + p + S +
							CURVE_FILE);

						int epiIndex = gmmIndex + epiCount;
						double wt = branchWT * gmmWT * epiWT;

						List<String> lines = Files.readLines(curveFile, US_ASCII);
						for (String line : Iterables.skip(lines, 1)) {
							String loc = StringUtils.substringBefore(line, ",");
							String curve = StringUtils.substringAfter(line, ",");

							String curveStr = epiIndex + "," + wt + "," + curve + LF;
							File cOut = new File(outPath + p + S + loc + S + CURVE_FILE);
							Files.append(curveStr, cOut, US_ASCII);

							String paramStr = epiIndex + "," + branchID + "," + gmmID + "," +
								epiID + LF;
							File pOut = new File(outPath + p + S + loc + S + PARAM_FILE);
							Files.append(paramStr, pOut, US_ASCII);
						}
						epiCount++;
					}
					gmmCount++;
				}
			}

		}
	}

	// init logic tree branch and weight lists
	private static void initBranchData(String srcPath, List<String> branchIDs,
			List<Double> branchWTs) {

		File firstSrcDir = new File(srcPath, "AS-OFF");
		File[] branchDirs = firstSrcDir.listFiles(new FileFilter() {
			@Override public boolean accept(File f) {
				return !f.getName().startsWith(".");
			}
		});

		APrioriBranchWeightProvider wtProvider = new APrioriBranchWeightProvider();

		for (File branchDir : branchDirs) {
			String branchName = branchDir.getName();
			branchIDs.add(branchName);
			LogicTreeBranch ltb = LogicTreeBranch.fromFileName(branchName);
			double wt = wtProvider.getWeight(ltb);
			branchWTs.add(wt);
		}
		DataUtils.asWeights(branchWTs); // normalize
	}

	private static void initOutputDirsAndFiles(String outPath, Set<Period> periods,
			Set<String> locIDs) throws IOException {
		for (Period p : periods) {
			for (String locID : locIDs) {

				File cFile = new File(outPath + p + S + locID + S + CURVE_FILE);
				Files.createParentDirs(cFile);
				// write header
				Iterable<String> gmVals = Collections2.transform(p.getFunction().xValues(),
					Functions.toStringFunction());
				List<String> headers = Lists.newArrayList("BranchIndex", "BranchWt");
				Iterable<String> cityFields = Iterables.concat(headers, gmVals);
				String cityHeader = JOIN.join(cityFields) + LF;
				Files.write(cityHeader, cFile, US_ASCII);

				File pFile = new File(outPath + p.name() + S + locID + S + PARAM_FILE);
				Files.createParentDirs(pFile);
				// write header
				Files.write("BranchIndex,BranchName,GMM,EPI" + LF, pFile, US_ASCII);
			}
		}
	}

	//@formatter:off
	
	@SuppressWarnings("unchecked")
	private static List<? extends LogicTreeBranchNode<?>> UC33nodeList = 
			Lists.newArrayList(
				FM3_1, FM3_2,
				ABM, GEOLOGIC, NEOKINEMA, ZENGBB,
				SHAW_2009_MOD, HANKS_BAKUN_08, ELLSWORTH_B, ELLB_SQRT_LENGTH, SHAW_CONST_STRESS_DROP,
				UNIFORM, TAPERED,
				RATE_6p5, RATE_7p9, RATE_9p6,
				MAG_7p3, MAG_7p6, MAG_7p9,
				UCERF2, UCERF3);
	
	@SuppressWarnings("unchecked")
	private static List<? extends Class<? extends LogicTreeBranchNode<?>>> UC33classList =
			Lists.newArrayList(
				FaultModels.class,
				DeformationModels.class,
				ScalingRelationships.class,
				SlipAlongRuptureModels.class,
				TotalMag5Rate.class,
				MaxMagOffFault.class,
				SpatialSeisPDF.class);

	
	private static Map<String, List<Boolean>> readFlags(File flagFile) throws IOException {
		List<String> flagLines = Files.readLines(flagFile, US_ASCII);
		// iteration ordered by leaf name keys
		Map<String, List<Boolean>> flagMap = Maps.newLinkedHashMap();
		// init lists
		Iterable<String> ID_Line = SPLIT.split(flagLines.get(0));
		for (String brID : ID_Line) {
			flagMap.put(brID, new ArrayList<Boolean>());
		}
		//load flags
		for (String flags : Iterables.skip(flagLines, 1)) {
			Iterator<Boolean> flagIt = Iterables.transform(
				SPLIT.split(flags), StrToBoolFn.INSTANCE).iterator();
			for (Entry<String, List<Boolean>> flagEntry : flagMap.entrySet()) {
				flagEntry.getValue().add(flagIt.next());
			}
		}
		return flagMap;
	}
	
	
	// creates and writes a table of flags indicating which nodes are used
	// in a branch; this simplifies process of building histograms in Matlab
	private static void createFlagFile(String branchPath, String flagPath)
			throws IOException {
		
		// branch list and index reverse lookup 
		List<String> branchLines = Files.readLines(new File(branchPath), US_ASCII);
		List<String> branchNames = Lists.newArrayList();
		
		// fill branch list and index lookup map
		for (String line : Iterables.skip(branchLines, 1)) {
			Iterable<String> vals = SPLIT.split(line);
			branchNames.add(Iterables.get(vals, 1) + "_" + Iterables.get(vals, 2) + "_" + 
					Iterables.get(vals, 3));
		}

		// want to redraw total GM distribution as stacked histograms for each
		// node - need a table of all branch nodes and flags [0,1] for each 
		// indicating if used on the branch at index.
		
		// intialize table of branches and nodes to false, first converting
		// logic tree identifiers to strings so we can incorporate gmm and epi 
		// branches
		
		Iterable<String> leafIDs = Iterables.concat(
			FluentIterable.from(UC33nodeList)
				.transform(LTB_ToStringFn.INSTANCE), 
				GMM_ID_WT_MAP.keySet(),
				EPI_ID_WT_MAP.keySet());
				
		ArrayTable<String, String, Boolean> nodeTable = ArrayTable.create(branchNames, leafIDs);
		for (String branchName : nodeTable.rowKeyList()) {
			for (String branchID : nodeTable.columnKeyList()) {
				nodeTable.put(branchName, branchID, false);
			}
		}
		
		// set trues
		for (String branchNameGmmEpi : branchNames) {
			int _idx1 = branchNameGmmEpi.lastIndexOf('_');
			String branchNameGmm = branchNameGmmEpi.substring(0, _idx1); // epi stripped
			int _idx2 = branchNameGmm.lastIndexOf('_');
			String branchID = branchNameGmm.substring(0, _idx2); // gmm stripped
					
			// UC3 logic tree branches
			LogicTreeBranch branch = LogicTreeBranch.fromFileName(branchID);
			for (Class<? extends LogicTreeBranchNode<?>> clazz : UC33classList) {
				String node = branch.getValueUnchecked(clazz).name();
				nodeTable.put(branchNameGmmEpi, node, true);
			}
			
			// GMMs
			String gmmID = branchNameGmm.substring(_idx2 + 1);
			nodeTable.put(branchNameGmmEpi, gmmID, true);
			
			// EPI
			String epiID = branchNameGmmEpi.substring(_idx1 + 1);
			nodeTable.put(branchNameGmmEpi, epiID, true);
		}
		
		// write file
		String nodeHeader = JOIN.join(leafIDs) + LF;
		File flagFile = new File(flagPath);
		Files.write(nodeHeader, flagFile, US_ASCII);
		for (String branchName : nodeTable.rowKeyList()) {
			List<Integer> flags = Lists.newArrayList();
			Map<String, Boolean>  row = nodeTable.row(branchName);
			for (String nodeID : nodeTable.columnKeyList()) {
				flags.add(row.get(nodeID) ? 1 : 0);
			}
			String flagLine = JOIN.join(flags) + LF;
			Files.append(flagLine, flagFile, US_ASCII);
		}
	}

	private enum LTB_ToStringFn implements Function<LogicTreeBranchNode<?>, String> {
		INSTANCE;
		@Override
		public String apply(LogicTreeBranchNode<?> ltbn) {
			return ltbn.name();
		}
	}
	
	private enum DoubleValueOfFn implements Function<String, Double> {
		INSTANCE;
		@Override public Double apply(String s) {
			return Double.valueOf(s.trim());
		}
	}

	private enum StrToBoolFn implements Function<String, Boolean> {
		INSTANCE;
		@Override public Boolean apply(String s) {
			return !s.equals("0");
		}
	}
	
	private enum Fractile {
		MIN,
		F02,
		F16,
		MEAN,
		F84,
		F98,
		MAX;
	}
		
	//@formatter:on

	private static ArbitrarilyDiscretizedFunc XYtoFunc(XY_DataSet xy) {
		ArbitrarilyDiscretizedFunc f = new ArbitrarilyDiscretizedFunc();
		for (Point2D p : xy) {
			f.set(p);
		}
		return f;
	}

	private static double getPE(DiscretizedFunc f, ProbOfExceed pe) {
		try {
			return f.getFirstInterpolatedX_inLogXLogYDomain(pe.annualRate());
		} catch (InvalidRangeException ire) {
			// we're probably below the last value on the curve so extrapolate
			// to use the current Interpolate implementation, xy values must be
			// reveresd (so that y is ascending) and flipped: supply x for y
			double[] Xs = Doubles.toArray(Lists.reverse(f.xValues()));
			double[] Ys = Doubles.toArray(Lists.reverse(f.yValues()));
			return Interpolate.findLogLogY(Ys, Xs, pe.annualRate());
		}
	}

	private static double getRate(DiscretizedFunc f, double iml) {
		try {
			return f.getInterpolatedY_inLogXLogYDomain(iml);
		} catch (InvalidRangeException ire) {
			// we're probably below the last value on the curve so extrapolate
			double[] Xs = Doubles.toArray(f.xValues());
			double[] Ys = Doubles.toArray(f.yValues());
			return Interpolate.findLogLogY(Xs, Ys, iml);
		}
	}

	
	// convert leaf data to from for matlab plotting
	private static void buildTornados(String dir, Iterable<Period> periods, String locsPath,
			ProbOfExceed pe) throws IOException {

		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locsPath);
		Iterable<String> locNames = locMap.keySet();

		List<Set<String>> keySetList = createTornadoKeyLists();

		// storing tornado rank counts: occurences of each branch
		// at each rank or level 0-8
		// will want to expand this to look at individual periods
		List<Multiset<String>> tornadoRanks = Lists.newArrayList();
		for (int i=0; i<keySetList.size(); i++) {
			Multiset<String> multiset = HashMultiset.create();
			tornadoRanks.add(multiset);
		}
		
		for (Period p : periods) {
			
			System.out.println(p);
			for (String name : locNames) {
				System.out.println("    " + name);
				String srcPath = dir + p.name() + S + name + S;

				// read noMean column leaves file
				File leafFile = new File(srcPath + "leaves_" + pe.name() + ".csv");
				List<String> leafLines = Files.readLines(leafFile, US_ASCII);
				Map<String, Double> leafMap = Maps.newHashMap();
				double mean = Double.valueOf(SPLIT.splitToList(leafLines.get(1)).get(1));
				for (String line : Iterables.skip(leafLines, 2)) {
					// skip header and NONE row
					List<String> parts = SPLIT.splitToList(line);
					double noMean = Double.valueOf(parts.get(2));
					leafMap.put(parts.get(0), mean / noMean);
				}

				// build strings of min and max for each branch and add
				// to sorted delta map
				TreeMap<Double, String> deltaMap = Maps.newTreeMap();
				for (Set<String> keySet : keySetList) {
					double minVal = Double.MAX_VALUE;
					String minKey = "";
					double maxVal = Double.MIN_VALUE;
					String maxKey = "";
					for (String key : keySet) {
						double val = leafMap.get(key);
						if (val < minVal) {
							minVal = val;
							minKey = key;
						}
						if (val > maxVal) {
							maxVal = val;
							maxKey = key;
						}
					}
					double delta = maxVal - minVal;
					List<?> outParts = Lists.newArrayList(minKey, minVal, maxKey, maxVal);
					String out = JOIN.join(outParts);
					deltaMap.put(delta, out);
				}

				// save map - while building leaf rank counts
				File tornadoFile = new File(srcPath, "tornado_" + pe.name() + ".csv");
				Files.write("", tornadoFile, US_ASCII);
				int count = 0;
				for (String line : deltaMap.descendingMap().values()) {
					Files.append(line + LF, tornadoFile, US_ASCII);
					
					String leaf = line.substring(0, line.indexOf(','));
					int leafIndex = lookupLeafIndex(keySetList, leaf);
					String branchName = branchNames.get(leafIndex);
					tornadoRanks.get(count).add(branchName);
					count++;
				}
			}
		}
		
		// write tornado ranks
		File tornadoRankFile = new File(dir, "tornadoRanks.csv");
		String tornadoRankHeader = "rank," + JOIN.join(branchNames) + LF;
		Files.write(tornadoRankHeader, tornadoRankFile, US_ASCII);
		int rankCount = 0;
		for (Multiset<String> rankData : tornadoRanks) {
			List<Integer> counts = Lists.newArrayList(rankCount);
			for (String branchName : branchNames) {
				counts.add(rankData.count(branchName));
			}
			String countStr = JOIN.join(counts) + LF;
			Files.append(countStr, tornadoRankFile, US_ASCII);
			rankCount++;
		}
	}

	private static int lookupLeafIndex(List<Set<String>> leafSetList, String leaf) {
		for (int i=0; i<leafSetList.size(); i++) {
			Set<String> leafSet = leafSetList.get(i);
			if (leafSet.contains(leaf)) return i;
		}
		throw new IllegalArgumentException("Bad leaf id: " + leaf);
	}
	
	// branch names in same order as tornado key lists
	private static List<String> branchNames = Lists.newArrayList("Fault Model",
		"Deformation Model", "Scaling Relationship", "Slip Along Rupture", "M\u22655 Rate",
		"Off-fault Mmax", "Smoothed Seismicity Model", "Ground Motion Model",
		"GMM Epistemic Uncertainty");

	private static List<Set<String>> createTornadoKeyLists() {
		List<Set<String>> setList = new ArrayList<Set<String>>();

		for (Class<? extends LogicTreeBranchNode<?>> clazz : UC33classList) {
			LogicTreeBranchNode<?>[] nodes = clazz.getEnumConstants();
			Set<String> ltbSet = Sets.newHashSet();
			Function<LogicTreeBranchNode<?>, String> fn = LTB_ToStringFn.INSTANCE;
			for (LogicTreeBranchNode<?> node : nodes) {
				if (!UC33nodeList.contains(node)) continue;
				ltbSet.add(fn.apply(node));
			}
			setList.add(ltbSet);
		}

		setList.add(GMM_ID_WT_MAP.keySet());
		setList.add(EPI_ID_WT_MAP.keySet());
		return setList;
	}

}
