package scratch.peter.curves;

import static com.google.common.base.Charsets.US_ASCII;
import static org.opensha.nshmp2.util.Period.GM0P00;
import static org.opensha.nshmp2.util.Period.GM0P20;
import static org.opensha.nshmp2.util.Period.GM1P00;
import static org.opensha.nshmp2.util.Period.GM4P00;
import static scratch.peter.curves.ProbOfExceed.PE10IN50;
import static scratch.peter.curves.ProbOfExceed.PE2IN50;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.io.IOUtils;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.exceptions.InvalidRangeException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.Interpolate;
import org.opensha.nshmp2.util.Period;

import scratch.peter.ucerf3.calc.UC3_CalcUtils;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

@SuppressWarnings("javadoc")
public class CurveUtilsNSHMP08 {

	private static final Joiner JOIN = Joiner.on(',');
	private static final Splitter SPLIT = Splitter.on(',').trimResults();
	private static final String S = File.separator;
	private static final String LF = IOUtils.LINE_SEPARATOR;
	private static final String ROOT = "tmp/UC33/";
	
	public static void main(String[] args) throws IOException {
		Set<Period> periods = EnumSet.of(GM0P00, GM0P20, GM1P00, GM4P00);
		
		String outPath = ROOT + "curves/src-reduce/NSHMP08-GMM/";
		String locsPath = ROOT + "curvejobs/sites/noPBR.txt";

		runBranchSummaries(outPath, periods, locsPath, PE2IN50);
	}
	
	private static void runBranchSummaries(String dir, 
			Iterable<Period> periods, String locsPath, ProbOfExceed pe)
			throws IOException {

		Map<String, Location> locMap = UC3_CalcUtils.readSiteFile(locsPath);
		Iterable<String> locNames = locMap.keySet();
		
		for (Period p : periods) {
			System.out.println(p);
			for (String name : locNames) {
				System.out.println("    " + name);
				String srcPath = dir + S + p.name() + S + name;
				summarizePE(srcPath, pe);
			}
		}
	}
	
	private static void summarizePE(String srcDir, ProbOfExceed pe) throws IOException {

		File curveFile1 = new File(srcDir, "Boore2008_curves.csv");
		File curveFile2 = new File(srcDir, "CB2008_curves.csv");
		File curveFile3 = new File(srcDir, "CY2008_curves.csv");

		XY_DataSetList curves = initCurves(curveFile1);
		curves.addAll(initCurves(curveFile2));
		curves.addAll(initCurves(curveFile3));

		List<Double> wtVals = initWeights(curveFile1);
		List<Double> wtCopy1 = Lists.newArrayList(wtVals);
		DataUtils.scale(0.3334, wtCopy1);
		List<Double> wtCopy2 = Lists.newArrayList(wtVals);
		DataUtils.scale(0.3333, wtCopy2);

		List<Double> weights = Lists.newArrayList();
		weights.addAll(wtCopy1);
		weights.addAll(wtCopy2);
		weights.addAll(wtCopy2);
		
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
		
	}
	
	private static XY_DataSetList initCurves(File curveFile) throws IOException {
		List<String> curveLines = Files.readLines(curveFile, US_ASCII);
		XY_DataSetList curves = new XY_DataSetList();
		
		// create model function: first line has BranchIndex, BranchWt, gm-vals ...
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
		for (String line : Iterables.skip(
				Files.readLines(curveFile, US_ASCII), 1)) {
			weights.add(FluentIterable
				.from(SPLIT.split(line))
				.transform(DoubleValueOfFn.INSTANCE)
				.get(1));
		}
		return weights;
	}
	
	private static Map<Fractile, XY_DataSet> computeFractiles(
			XY_DataSetList curves, List<Double> weights) {
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
	
	private static void writeFractiles(Map<Fractile, XY_DataSet> fracMap,
			File fracFile) throws IOException {
		// header
		XY_DataSet model = fracMap.get(Fractile.MIN);
		Iterable<String> fracFields = Iterables.concat(
			Lists.newArrayList("fractile", "IML_At2in50", "RateAtMeanIML"),
			Iterables.transform(model.xValues(),
				Functions.toStringFunction()));
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
				Iterables.transform(fracEntry.getValue().yValues(),
					Functions.toStringFunction()));
			Files.append(JOIN.join(fracData) + LF, fracFile, US_ASCII);
		}
	}
	
	private static void writeDistros(List<Double> gms, List<Double> rates,
			List<Double> weights, File distroFile) throws IOException {
		Files.write("wt,gm,rate" + LF, distroFile, US_ASCII);
		for (int i = 0; i < gms.size(); i++) {
			Files.append(weights.get(i) + "," + gms.get(i) + "," +
				rates.get(i) + LF, distroFile, US_ASCII);
		}
	}
	
		
	private enum DoubleValueOfFn implements Function<String, Double> {
		INSTANCE;
		@Override public Double apply(String s) {
			return Double.valueOf(s.trim());
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
		
	private static ArbitrarilyDiscretizedFunc XYtoFunc(XY_DataSet xy) {
		ArbitrarilyDiscretizedFunc f = new ArbitrarilyDiscretizedFunc();
		for (Point2D p : xy) {
			f.set(p);
		}
		return f;
	}
	
	// returns the ground motion at a specified prob of exceedance
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

	// returns the prob of exceedance (or rate) at a specified ground motion
	private static double getRate(DiscretizedFunc f,double iml) {
		try {
			return f.getInterpolatedY_inLogXLogYDomain(iml);
		} catch (InvalidRangeException ire) {
			// we're probably below the last value on the curve so extrapolate
			double[] Xs = Doubles.toArray(f.xValues());
			double[] Ys = Doubles.toArray(f.yValues());
			return Interpolate.findLogLogY(Xs, Ys, iml);
		}
	}

}
