package scratch.kevin.bbp;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.util.component.ShahiBaker2014Trans;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SpectraPlotter {
	
	private static DiscretizedFunc[] loadAll(File file) throws IOException {
		return loadAll(Files.readLines(file, Charset.defaultCharset()));
	}
	
	private static DiscretizedFunc[] loadAll(List<String> lines) {
		DiscretizedFunc[] ret = null;
		String[] names = null;
		for (String line : lines) {
			line = line.trim();
			if (names == null && line.startsWith("#")) {
				String nameStr = line.substring(1).trim();
				if (line.contains("\t"))
					names = nameStr.split("\t");
				else
					names = nameStr.split("\\s+");
//				System.out.println("Names:");
//				for (String name : names)
//					System.out.println("\t"+name);
			}
			if (line.startsWith("#") || line.isEmpty())
				continue;
			String[] split = line.split("\\s+");
			if (ret == null) {
				Preconditions.checkState(split.length > 1);
				ret = new DiscretizedFunc[split.length - 1];
				for (int i=0; i<ret.length; i++) {
					ret[i] = new ArbitrarilyDiscretizedFunc();
					if (names != null) {
						if (names.length == ret.length)
							ret[i].setName(names[i].trim());
						else if (names.length == ret.length+1) // has x axis
							ret[i].setName(names[i+1].trim());
					}
				}
			} else {
				Preconditions.checkState(ret.length == split.length-1);
			}
			double x = Double.parseDouble(split[0]);
			for (int i=0; i<ret.length; i++) {
				// *'s mean zero seimograms, which should never happen
//				if (split[i+1].startsWith("***"))
//					ret[i].set(x, 0d);
//				else
					ret[i].set(x, Double.parseDouble(split[i+1]));
			}
		}
		return ret;
	}
	
	public static DiscretizedFunc loadRotD50(File file) throws IOException {
		DiscretizedFunc[] ret = loadAll(file);
		Preconditions.checkState(ret.length == 3 || ret.length == 4);
		return new LightFixedXFunc(ret[2]);
	}
	
	public static DiscretizedFunc loadRotD50(List<String> lines) {
		DiscretizedFunc[] ret = loadAll(lines);
		Preconditions.checkState(ret.length == 3 || ret.length == 4);
		return new LightFixedXFunc(ret[2]);
	}
	
	public static DiscretizedFunc loadRotD100(File file) throws IOException {
		DiscretizedFunc[] ret = loadAll(file);
		Preconditions.checkState(ret.length == 4);
		return new LightFixedXFunc(ret[3]);
	}
	
	public static DiscretizedFunc loadRotD100(List<String> lines) {
		DiscretizedFunc[] ret = loadAll(lines);
		Preconditions.checkState(ret.length == 4);
		return new LightFixedXFunc(ret[3]);
	}
	
	/**
	 * 
	 * @param file
	 * @return array of [RotD50, RotD100]
	 * @throws IOException
	 */
	public static DiscretizedFunc[] loadRotD(File file) throws IOException {
		DiscretizedFunc[] ret = loadAll(file);
		Preconditions.checkState(ret.length == 4);
		return new DiscretizedFunc[] { new LightFixedXFunc(ret[2]), new LightFixedXFunc(ret[3]) };
	}
	
	/**
	 * 
	 * @param lines
	 * @return array of [RotD50, RotD100]
	 */
	public static DiscretizedFunc[] loadRotD(List<String> lines) {
		DiscretizedFunc[] ret = loadAll(lines);
		Preconditions.checkState(ret.length == 4);
		return new DiscretizedFunc[] { new LightFixedXFunc(ret[2]), new LightFixedXFunc(ret[3]) };
	}

	public static void plotRotD(File file, boolean rd50, boolean rd100) throws IOException {
		plotRotD(file, null, null, rd50, rd100, null);
	}
	
	public static void plotRotD(File file, File outputDir, String prefix, boolean rd50, boolean rd100, UncertainArbDiscDataset[] gmpes)
			throws IOException {
		Preconditions.checkState(rd50 || rd100);
		if (file.isDirectory()) {
			for (File sub : file.listFiles()) {
				if (!rd100 && sub.getName().endsWith(".rd50"))
					plotRotD(sub, outputDir, prefix, rd50, rd100, gmpes);
				else if (sub.getName().endsWith(".rd100"))
					plotRotD(sub, outputDir, prefix, rd50, rd100, gmpes);
			}
			return;
		}
		
		System.out.println("Plotting "+file.getAbsolutePath());
		DiscretizedFunc[] funcArray = loadAll(file);
		
		if (outputDir == null)
			outputDir = file.getParentFile();
		if (prefix == null)
			prefix = file.getName();
		
		plotSpectra(funcArray, file.getName(), "Period (s)", "PSA (g)", true, true, new Range(1e-3, 1e1),
				outputDir, prefix, gmpes, false, rd50, rd100);
	}
	
	public static DiscretizedFunc loadFAS(List<String> lines) {
		DiscretizedFunc[] ret = loadAll(lines);
		return ret[ret.length-1];
	}
	
	public static DiscretizedFunc loadFAS(File file) throws IOException {
		DiscretizedFunc[] ret = loadAll(file);
		return ret[ret.length-1];
	}

	public static void plotFAS(File file) throws IOException {
		plotFAS(file, null, null);
	}
	public static void plotFAS(File file, File outputDir, String prefix) throws IOException {
		if (file.isDirectory()) {
			for (File sub : file.listFiles()) {
				if (sub.getName().endsWith(".fs.col") || (sub.isDirectory() && sub.getName().equals("FAS")))
					plotFAS(sub, outputDir, prefix);
			}
			return;
		}
		
		System.out.println("Plotting "+file.getAbsolutePath());
		DiscretizedFunc[] funcArray = loadAll(file);
		
		if (outputDir == null)
			outputDir = file.getParentFile();
		if (prefix == null)
			prefix = file.getName();
		
		plotSpectra(funcArray, file.getName(), "Frequency (hz)", "Fourier Amplitude (cm/s)",
				true, true, null, outputDir, prefix, null, true, false, false);
	}
	
	private static Color[] GMPE_COLORS = { Color.RED.brighter(), Color.YELLOW.brighter(),
			Color.CYAN.brighter(), Color.MAGENTA.brighter() };
	
	private static void plotSpectra(DiscretizedFunc[] funcArray, String title, String xAxisLable, String yAxisLabel,
			boolean xLog, boolean yLog, Range yRange, File outputDir, String prefix, UncertainArbDiscDataset[] gmpes,
			boolean fas, boolean rd50, boolean rd100) throws IOException {
		Preconditions.checkState(fas || rd50 || rd100);
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (gmpes != null) {
			for (int i=0; i<gmpes.length; i++) {
				Color color = GMPE_COLORS[i % GMPE_COLORS.length];
				PlotCurveCharacterstics outsideChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color);
				PlotCurveCharacterstics meanChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, color);
				
				funcs.add(gmpes[i]);
				chars.add(meanChar);
				
				funcs.add(gmpes[i].getUpper());
				chars.add(outsideChar);
				
				funcs.add(gmpes[i].getLower());
				chars.add(outsideChar);
			}
		}
		
		for (int i=0; i<funcArray.length; i++) {
			funcs.add(funcArray[i]);
			if (i == funcArray.length -1 && (fas || (rd50 && !rd100))
					|| i == funcArray.length -2 && rd50 && rd100) {
				// RotD50 or FAS
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			} else if (i == funcArray.length -1 && rd100) {
				// RotD100
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK));
			} else {
				chars.add(new PlotCurveCharacterstics(
						individualCompLineTypes[i % individualCompLineTypes.length],1f, Color.GRAY));
			}
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLable, yAxisLabel);
		spec.setLegendVisible(funcArray[0].getName() != null);
		
		HeadlessGraphPanel gp = buildGP();
		
		gp.drawGraphPanel(spec, xLog, yLog, getXRange(funcs, xLog), yRange);
		
		gp.getChartPanel().setSize(800, 800);
		File file = new File(outputDir, prefix);
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static final Range getXRange(List<? extends XY_DataSet> funcs, boolean log) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (XY_DataSet func : funcs) {
			min = Math.min(min, func.getMinX());
			max = Math.max(max, func.getMaxX());
		}
		if (log)
			return getNiceLogRange(min, max);
		return new Range(min, max);
	}
	
	private static Range getNiceLogRange(double min, double max) {
		double lowerLog = Math.log10(min);
		double lowerPower = Math.floor(lowerLog);
		if (lowerLog - lowerPower > 0.9)
			lowerPower = lowerLog;
		min = Math.pow(10, lowerPower);
		double upperLog = Math.log10(max);
		double upperPower = Math.ceil(upperLog);
		if (upperPower - upperLog > 0.9)
			upperPower = upperLog;
		max = Math.pow(10, upperPower);
		return new Range(min, max);
	}
	
	private static final Range getYRange(List<? extends XY_DataSet> funcs, boolean log) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		List<XY_DataSet> allFuncs = new ArrayList<>();
		for (XY_DataSet func : funcs) {
			if (func instanceof UncertainArbDiscDataset) {
				UncertainArbDiscDataset rangeFunc = (UncertainArbDiscDataset)func;
				allFuncs.add(rangeFunc.getUpper());
				allFuncs.add(rangeFunc.getLower());
			}
			allFuncs.add(func);
		}
		for (XY_DataSet func : allFuncs) {
			min = Math.min(min, func.getMinY());
			max = Math.max(max, func.getMaxY());
		}
		if (log)
			return getNiceLogRange(min, max);
		return new Range(min, max);
	}
	
	private static final PlotLineType[] individualCompLineTypes = {PlotLineType.DASHED, PlotLineType.DOTTED,
			PlotLineType.DOTTED_AND_DASHED};
	
	private static final DecimalFormat percentDF = new DecimalFormat("0.##%");
	
	public static void plotMultiRotD50(List<File> refFiles, String refName, File dataFile, String dataName, String title,
			File outputDir, String prefix, UncertainArbDiscDataset[] gmpes) throws IOException {
		List<DiscretizedFunc> refSpectra = new ArrayList<>();
		for (File refFile : refFiles)
			refSpectra.add(loadRotD50(refFile));
		DiscretizedFunc dataSpectra = loadRotD50(dataFile);
		plotMultiRotD50(refSpectra, refName, dataSpectra, dataName, title, outputDir, prefix, gmpes);
	}
	
	public static void plotMultiRotD50(List<DiscretizedFunc> refFiles, String refName, DiscretizedFunc dataSpectra, String dataName,
			String title, File outputDir, String prefix, UncertainArbDiscDataset[] gmpes) throws IOException {
		plotMultiSpectra(refFiles, refName, dataSpectra, dataName, title, outputDir, prefix, true, gmpes);
	}
	
	public static void plotMultiFAS(List<File> refFiles, String refName, File dataFile, String dataName, String title,
			File outputDir, String prefix) throws IOException {
		List<DiscretizedFunc> refSpectra = new ArrayList<>();
		for (File refFile : refFiles)
			refSpectra.add(loadFAS(refFile));
		DiscretizedFunc dataSpectra = loadFAS(dataFile);
		plotMultiFAS(refSpectra, refName, dataSpectra, dataName, title, outputDir, prefix);
	}
	
	public static void plotMultiFAS(List<DiscretizedFunc> refFiles, String refName, DiscretizedFunc dataSpectra,
			String dataName, String title, File outputDir, String prefix) throws IOException {
		plotMultiSpectra(refFiles, refName, dataSpectra, dataName, title, outputDir, prefix, false, null);
	}
	
	private static void plotMultiSpectra(List<DiscretizedFunc> refSpectra, String refName, DiscretizedFunc dataSpectra,
			String dataName, String title, File outputDir, String prefix, boolean rotD50, UncertainArbDiscDataset[] gmpes)
					throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		populateRefRangeFuncs(refSpectra, refName, funcs, chars);
		
		if (gmpes != null) {
			for (int i=0; i<gmpes.length; i++) {
				Color color = GMPE_COLORS[i % GMPE_COLORS.length];
				PlotCurveCharacterstics outsideChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color);
				PlotCurveCharacterstics meanChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, color);
				
				funcs.add(gmpes[i]);
				chars.add(meanChar);
				
				funcs.add(gmpes[i].getUpper());
				chars.add(outsideChar);
				
				funcs.add(gmpes[i].getLower());
				chars.add(outsideChar);
			}
		}
		
		dataSpectra.setName(dataName);
		funcs.add(dataSpectra);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		String xAxisLabel, yAxisLabel;
		if (rotD50) {
			xAxisLabel = "Period (s)";
			yAxisLabel =  "RotD50 (g)";
		} else {
			xAxisLabel = "Frequency (hz)";
			yAxisLabel =  "Fourier Amplitude (cm/s)";
		}
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = buildGP();
		
		Range yRange;
		if (rotD50)
			yRange = new Range(1e-3, 1e1);
		else
			yRange = getYRange(funcs, true);
		gp.drawGraphPanel(spec, true, true, getXRange(funcs, true), yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 800);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static Range populateRefRangeFuncs(List<DiscretizedFunc> refSpectra, String refName,
			List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars) {
		XY_DataSetList refFuncs = new XY_DataSetList();
		List<Double> relativeWts = new ArrayList<>();
		for (DiscretizedFunc refFunc : refSpectra){
			relativeWts.add(1d);
			refFuncs.add(refFunc);
		}
		FractileCurveCalculator refFractCalc = new FractileCurveCalculator(refFuncs, relativeWts);
		
		List<double[]> fractiles = new ArrayList<>();
		List<Color> fractileColors = new ArrayList<>();
		
//		0.025, 0.16, 0.84, 0.975
		fractiles.add(new double[] { 0d, 1d });
		fractileColors.add(new Color(220, 220, 220)); // very light gray
		
		fractiles.add(new double[] { 0.025, 0.975 });
		fractileColors.add(Color.LIGHT_GRAY);
		
		fractiles.add(new double[] { 0.16, 0.84 });
		fractileColors.add(Color.GREEN.brighter());
		
		DiscretizedFunc refMean = (DiscretizedFunc) refFractCalc.getMeanCurve();
		refMean.setName(refName+" Mean");
		
		double maxY = 0d;
		double minY = Double.POSITIVE_INFINITY;
		
		for (int i=0; i<fractiles.size(); i++) {
			double[] range = fractiles.get(i);
			DiscretizedFunc lowerFunc = (DiscretizedFunc)refFractCalc.getFractile(range[0]);
			DiscretizedFunc upperFunc = (DiscretizedFunc)refFractCalc.getFractile(range[1]);
			UncertainArbDiscDataset rangeFunc = new UncertainArbDiscDataset(refMean, lowerFunc, upperFunc);
			rangeFunc.setName(percentDF.format(range[0])+"-"+percentDF.format(range[1]));
			funcs.add(rangeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, fractileColors.get(i)));
			maxY = Math.max(maxY, upperFunc.getMaxY());
//			System.out.println("Fractiles for "+rangeFunc.getName());
//			System.out.println(rangeFunc);
			for (Point2D pt : lowerFunc)
				if (pt.getY() > 0)
					minY = Math.min(minY, pt.getY());
		}
		
		funcs.add(refMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		
		return new Range(minY, maxY);
	}
	
	private static Range rot_d_ratio_y_range = new Range(1d, 1.5d);
	
	public static void plotRotDRatio(File dataFile, String dataName,
			String title, File outputDir, String prefix) throws IOException {
		plotRotDRatio(loadRotD(dataFile), dataName, null, null, title, outputDir, prefix);
	}
	
	public static void plotRotDRatio(DiscretizedFunc[] dataRD, String dataName,
			String title, File outputDir, String prefix) throws IOException {
		plotRotDRatio(dataRD, dataName, null, null, title, outputDir, prefix);
	}
	
	public static void plotRotDRatio(List<DiscretizedFunc[]> refRDs, String refName, String title, File outputDir, String prefix)
			throws IOException {
		plotRotDRatio(null, null, refRDs, refName, title, outputDir, prefix);
	}
	
	public static void plotRotDRatio(DiscretizedFunc[] dataRD, String dataName, List<DiscretizedFunc[]> refRDs, String refName,
			String title, File outputDir, String prefix) throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DiscretizedFunc refPeriods = null;
		
		if (refRDs != null && !refRDs.isEmpty()) {
			System.out.println("Calculating ratios for "+refRDs.size()+" inputs");
			List<DiscretizedFunc> refRatios = new ArrayList<>();
			for (DiscretizedFunc[] refRD : refRDs)
				refRatios.add(SimulationRotDProvider.calcRotDRatio(refRD));
			System.out.println("Done calculating ratios, populating range funcs");
			populateRefRangeFuncs(refRatios, refName, funcs, chars);
			System.out.println("Done populating range funcs");
			refPeriods = refRatios.get(0);
		}
		
		if (dataRD != null) {
			DiscretizedFunc dataRatio = SimulationRotDProvider.calcRotDRatio(dataRD);
			dataRatio.setName(dataName);
			funcs.add(dataRatio);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			refPeriods = dataRatio;
		}
		
		// Shahi & Baker
		ShahiBaker2014Trans shahi = new ShahiBaker2014Trans();
		DiscretizedFunc shahiFunc = new ArbitrarilyDiscretizedFunc();
		Preconditions.checkNotNull(refPeriods, "Must supply either single or range of data (or both)");
		for (Point2D pt : refPeriods) {
			double period = pt.getX();
			if (period < shahi.getMinPeriod() || period > shahi.getMaxPeriod())
				continue;
			shahiFunc.set(period, shahi.getScalingFactor(period));
		}
		shahiFunc.setName(shahi.getName());
		funcs.add(shahiFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
		
		String xAxisLabel = "Period (s)";
		String yAxisLabel = "RotD100/RotD50 Ratio";
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = buildGP();
		gp.drawGraphPanel(spec, true, false, getXRange(funcs, true), rot_d_ratio_y_range);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	public static void plotRotDRatioDependence(List<DiscretizedFunc[]> rds, List<Double> scalars, String scalarLabel,
			int numBins, double[] periods, String dataName, String title, File outputDir, String prefix, boolean logX)
					throws IOException {
		List<List<Double>> periodScalars = new ArrayList<>();
		for (int p=0; p<periods.length; p++)
			periodScalars.add(scalars);
		plotRotDRatioPeriodDependence(rds, periodScalars, scalarLabel, numBins, periods, dataName, title, outputDir, prefix, logX);
	}
	
	public static void plotRotDRatioPeriodDependence(List<DiscretizedFunc[]> rds, List<List<Double>> periodScalars, String scalarLabel,
			int numBins, double[] periods, String dataName, String title, File outputDir, String prefix, boolean logX)
					throws IOException {
		Preconditions.checkArgument(periodScalars.size() == periods.length,
				"Have %s periods but %s lists", periods.length, periodScalars.size());
		Preconditions.checkArgument(periodScalars.get(0).size() == rds.size(),
				"Each list should be %s long (first is %s)", rds.size(), periodScalars.get(0).size());
		
		List<DiscretizedFunc> ratios = new ArrayList<>();
		for (DiscretizedFunc[] rd : rds)
			ratios.add(SimulationRotDProvider.calcRotDRatio(rd));
		
		double minScalar = Double.POSITIVE_INFINITY;
		double maxScalar = 0;
		for (List<Double> myScalars : periodScalars) {
			for (double scalar : myScalars) {
				Preconditions.checkState(!logX || scalar > 0);
				minScalar = Math.min(minScalar, scalar);
				maxScalar = Math.max(maxScalar, scalar);
			}
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		List<XY_DataSet> rangeFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> rangeChars = new ArrayList<>();
		
		CPT periodCPT = new CPT(Doubles.min(periods), Doubles.max(periods), Color.DARK_GRAY, Color.RED.darker());
		
		int minBinCount = Integer.MAX_VALUE;
		int maxBinCount = 0;
		
		for (int p=0; p<periods.length; p++) {
			double period = periods[p];
			EvenlyDiscretizedFunc xValsFunc;
			if (logX)
				xValsFunc = new EvenlyDiscretizedFunc(Math.log10(minScalar), Math.log10(maxScalar), numBins);
			else
				xValsFunc = new EvenlyDiscretizedFunc(minScalar, maxScalar, numBins);
			ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
			List<List<Double>> vals = new ArrayList<>();
			for (int i=0; i<numBins; i++)
				vals.add(new ArrayList<>());
			
			List<Double> myScalars = periodScalars.get(p);
			
			for (int i=0; i<ratios.size(); i++) {
				int binIndex;
				if (logX)
					binIndex = xValsFunc.getClosestXIndex(Math.log10(myScalars.get(i)));
				else
					binIndex = xValsFunc.getClosestXIndex(myScalars.get(i));
				vals.get(binIndex).add(ratios.get(i).getInterpolatedY(period));
			}
			for (int i=0; i<numBins; i++) {
				double[] binVals = Doubles.toArray(vals.get(i));
				int count = binVals.length;
				if (count < minBinCount)
					minBinCount = count;
				if (count > maxBinCount)
					maxBinCount = count;
				double mean, stdDev;
				if (count > 0) {
					mean = StatUtils.mean(binVals);
					stdDev = Math.sqrt(StatUtils.variance(binVals));
					double x;
					if (logX)
						x = Math.pow(10, xValsFunc.getX(i));
					else
						x = xValsFunc.getX(i);
					meanFunc.set(x, mean);
					lowerFunc.set(x, mean - stdDev);
					upperFunc.set(x, mean + stdDev);
				}
			}
			
			UncertainArbDiscDataset rangeFunc = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
			Color c = periodCPT.getColor((float)period);
			
			rangeFuncs.add(rangeFunc);
			rangeChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(c.getRed(), c.getGreen(), c.getBlue(), 50)));
			
			meanFunc.setName((float)period+"s");
			funcs.add(meanFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 3.5f, c));
		}
		
		funcs.addAll(0, rangeFuncs);
		chars.addAll(0, rangeChars);
		
		String xAxisLabel = scalarLabel;
		String yAxisLabel = "RotD100/RotD50 Ratio";
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = buildGP();
		
		Range xRange;
		if (logX)
			xRange = getNiceLogRange(minScalar, maxScalar);
		else
			xRange = new Range(minScalar, maxScalar);
		
		List<XYTextAnnotation> anns = new ArrayList<>();
		double annY = rot_d_ratio_y_range.getUpperBound();
		XYTextAnnotation countAnn = new XYTextAnnotation("  "+rds.size()+" Spectra", xRange.getLowerBound(), annY);
		countAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 22));
		countAnn.setTextAnchor(TextAnchor.TOP_LEFT);
		anns.add(countAnn);
		XYTextAnnotation binsAnn = new XYTextAnnotation("["+minBinCount+" - "+maxBinCount+"] Per Bin  ", xRange.getUpperBound(), annY);
		binsAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 22));
		binsAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		anns.add(binsAnn);
		spec.setPlotAnnotations(anns);
		
		gp.drawGraphPanel(spec, logX, false, xRange, rot_d_ratio_y_range);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	public static void plotRotDRatioVsRd50Scatter(List<DiscretizedFunc[]> rds, double period,
			String dataName, String title, File outputDir, String prefix, int numMeanBins) throws IOException {
		List<Double> scalars = new ArrayList<>();
		for (DiscretizedFunc[] rd : rds)
			scalars.add(rd[0].getInterpolatedY(period));
		plotRotDRatioScatter(rds, scalars, (float)period+"s RotD50", period, dataName, title, outputDir, prefix, true, numMeanBins);
	}
	
	public static void plotRotDRatioScatter(List<DiscretizedFunc[]> rds, List<Double> scalars, String scalarLabel,
			double period, String dataName, String title, File outputDir, String prefix, boolean logX, int numMeanBins)
					throws IOException {
		Preconditions.checkArgument(scalars.size() == rds.size());
		
		List<DiscretizedFunc> ratios = new ArrayList<>();
		for (DiscretizedFunc[] rd : rds)
			ratios.add(SimulationRotDProvider.calcRotDRatio(rd));
		
		double minScalar = Double.POSITIVE_INFINITY;
		double maxScalar = 0;
		for (double dist : scalars) {
			minScalar = Math.min(minScalar, dist);
			maxScalar = Math.max(maxScalar, dist);
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		
		for (int i=0; i<ratios.size(); i++)
			scatter.set(scalars.get(i), ratios.get(i).getInterpolatedY(period));
		
		funcs.add(scatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		
		// regression
//		SimpleRegression regression = new SimpleRegression();
//		for (Point2D pt : scatter) {
//			if (logX)
//				regression.addData(Math.log10(pt.getX()), pt.getY());
//			else
//				regression.addData(pt.getX(), pt.getY());
//		}
//		double b = regression.getIntercept();
//		double m = regression.getSlope();
//		DefaultXY_DataSet fit = new DefaultXY_DataSet();
//		double[] regXs = { minScalar, maxScalar };
//		for (double origX : regXs) {
//			double x = origX;
//			if (logX)
//				x = Math.log10(origX);
//			double y = m*x + b;
//			fit.set(origX, y);
//		}
//		funcs.add(fit);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
		if (numMeanBins > 0) {
			EvenlyDiscretizedFunc xVals;
			if (logX) {
				Preconditions.checkState(minScalar > 0);
				// TODO this puts a bin center on the edges of the data, fix
				xVals = new EvenlyDiscretizedFunc(Math.log10(minScalar), Math.log10(maxScalar), numMeanBins);
			} else {
				xVals = new EvenlyDiscretizedFunc(minScalar, maxScalar, numMeanBins);
			}
			List<List<Double>> valsList = new ArrayList<>();
			for (int i=0; i<numMeanBins; i++)
				valsList.add(new ArrayList<>());
			for (Point2D pt : scatter) {
				int binIndex;
				if (logX)
					binIndex = xVals.getClosestXIndex(Math.log10(pt.getX()));
				else
					binIndex = xVals.getClosestXIndex(pt.getX());
				valsList.get(binIndex).add(pt.getY());
			}
			
			ArbitrarilyDiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
			
			for (int i=0; i<numMeanBins; i++) {
				double[] binVals = Doubles.toArray(valsList.get(i));
				int count = binVals.length;
				double mean, stdDev;
				if (count > 0) {
					mean = StatUtils.mean(binVals);
					stdDev = Math.sqrt(StatUtils.variance(binVals));
					double x;
					if (logX)
						x = Math.pow(10, xVals.getX(i));
					else
						x = xVals.getX(i);
					meanFunc.set(x, mean);
					lowerFunc.set(x, mean - stdDev);
					upperFunc.set(x, mean + stdDev);
				}
			}
			
			UncertainArbDiscDataset rangeFunc = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
			Color c = new Color(0, 200, 0);
			
			funcs.add(rangeFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(c.getRed(), c.getGreen(), c.getBlue(), 50)));
			
			meanFunc.setName((float)period+"s");
			funcs.add(meanFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 3.5f, c));
		}
		
		String xAxisLabel = scalarLabel;
		String yAxisLabel = (period)+"s RotD100/RotD50 Ratio";
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = buildGP();
		
		Range xRange;
		if (logX)
			xRange = getNiceLogRange(minScalar, maxScalar);
		else
			xRange = new Range(minScalar, maxScalar);
		
		gp.drawGraphPanel(spec, logX, false, xRange, rot_d_ratio_y_range);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 500);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
//		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static HeadlessGraphPanel buildGP() {
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(20);
		gp.setBackgroundColor(Color.WHITE);
		return gp;
	}
	
	public static File findRotD50File(File dir, String siteName) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (name.endsWith(".rd50") && name.contains(siteName))
				return file;
		}
		throw new FileNotFoundException("No .rd50 files found for site "+siteName+" in "+dir.getAbsolutePath());
	}
	
	public static File findRotD100File(File dir, String siteName) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (name.endsWith(".rd100") && name.contains(siteName))
				return file;
		}
		throw new FileNotFoundException("No .rd100 files found for site "+siteName+" in "+dir.getAbsolutePath());
	}
	
	public static File findFASFile(File dir, String siteName) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (file.isDirectory() && name.equals("FAS"))
				return findFASFile(file, siteName);
			if (name.endsWith(".fs.col") && name.contains(siteName))
				return file;
		}
		throw new FileNotFoundException("No .fs.col files found for site "+siteName+" in "+dir.getAbsolutePath());
	}
	
	public static List<File> findRunFiles(File dir, String siteName, boolean rotD50) throws FileNotFoundException {
		Preconditions.checkState(dir.exists());
		File resultsDir = new File(dir, "results");
		if (!resultsDir.exists())
			resultsDir = dir;
		
		List<File> files = new ArrayList<>();
		for (File runDir : resultsDir.listFiles()) {
			if (runDir.getName().startsWith("run_")) {
				if (rotD50)
					files.add(findRotD50File(runDir, siteName));
				else
					files.add(findFASFile(runDir, siteName));
			}
		}
		
		return files;
	}
	
	public static UncertainArbDiscDataset calcGMPE_RotD50(EqkRupture rupture, BBP_Site bbpSite, VelocityModel vm, ScalarIMR gmpe) {
		Site site = bbpSite.buildGMPE_Site(vm);
		
		gmpe.setSite(site);
		gmpe.setEqkRupture(rupture);
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param saParam = (SA_Param)gmpe.getIntensityMeasure();
		List<Double> periods = saParam.getPeriodParam().getAllowedDoubles();
		
		DiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
		
//		System.out.println("Mag: "+rupture.getMag());
//		System.out.println("Rake: "+rupture.getAveRake());
//		System.out.println("rRup: "+gmpe.getParameter(DistanceRupParameter.NAME).getValue());
		
		for (double period : periods) {
			SA_Param.setPeriodInSA_Param(saParam, period);
			
			double logMean = gmpe.getMean();
			double stdDev = gmpe.getStdDev();
//			System.out.println("Mean: "+logMean+"\tStdDev: "+stdDev);
			double logUpper = logMean + stdDev;
			double logLower = logMean - stdDev;
			
			upperFunc.set(period, Math.exp(logUpper));
			meanFunc.set(period, Math.exp(logMean));
			lowerFunc.set(period, Math.exp(logLower));
		}
		
		UncertainArbDiscDataset func = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
		func.setName(gmpe.getShortName()+" ±σ");
		return func;
	}

	public static void main(String[] args) throws IOException {
		plotRotDRatio(new File("/tmp/bbp_test1/774284145569850.SBSM.rd100"), "Test Data", "Test Ratio",
				new File("/tmp/bbp_test1"), "rot_d_ratio");
//		plotRotD(new File("/tmp/bbp_test1"), true, true);
//		plotRotD50(new File("/home/kevin/bbp/bbp_data/outdata/6462338"));
//		plotRotD50(new File("/home/kevin/bbp/bbp_data/outdata/6450680"));
//		plotRotD50(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/"
//				+ "event_4099020_0.05s_ADJ_VEL_seis"));
//		plotRotD50(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/event_4099020_0.05s_ADJ_VEL_bbp"));
//		plotFAS(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/event_4099020_0.05s_ADJ_VEL_bbp"));
		System.exit(0);
		
//		ScalarIMR[] gmpes = null;
		ScalarIMR[] gmpes = { AttenRelRef.ASK_2014.instance(null), AttenRelRef.BSSA_2014.instance(null),
				AttenRelRef.CB_2014.instance(null), AttenRelRef.CY_2014.instance(null) };
		
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
		File refDir = new File("/home/kevin/bbp/parallel/2017_10_04-rundir2194_long-event136704-dx1.16-noHF/results_orig");
		int eventID = 136704;
		
//		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
//		File refDir = new File("/home/kevin/bbp/parallel/2017_10_04-JG_UCERF3_millionElement-event4099020-dx0.48-noHF/results");
//		int eventID = 4099020;
		
		File srfDir = new File(catalog.getCatalogDir(), "event_srfs");
		File rsDir = new File(srfDir, "event_"+eventID+"_0.05s_ADJ_VEL_bbp");
		
		List<BBP_Site> sites = BBP_Site.readFile(refDir);
		
		int numRefRuns = 200;
		VelocityModel vm = VelocityModel.LA_BASIN;
		double minFractForInclusion = 0.2;
		
		EqkRupture gmpeRup = null;
		if (gmpes != null) {
			for (ScalarIMR gmpe : gmpes)
				gmpe.setParamDefaults();
			System.out.println("Loading event...");
			RSQSimEvent event = catalog.loader().byID(eventID);
			gmpeRup = catalog.getGMPE_Rupture(event, minFractForInclusion);
			System.out.println("DONE");
		}
		
		for (BBP_Site site : sites) {
			String siteName = site.getName();
			File rsRD50File = findRotD50File(rsDir, siteName);
			UncertainArbDiscDataset[] gmpeSpectra = null;
			if (gmpes != null) {
				gmpeSpectra = new UncertainArbDiscDataset[gmpes.length];
				for (int i=0; i<gmpes.length; i++) {
					System.out.println("Calculating spectra for "+gmpes[i].getShortName());
					gmpeSpectra[i] = calcGMPE_RotD50(gmpeRup, site, vm, gmpes[i]);
				}
				System.out.println("DONE spectra");
			}
			plotRotD(rsRD50File, null, null, true, false, gmpeSpectra);
			List<File> refRD50Files = new ArrayList<>();
			File rsFASFile = findFASFile(rsDir, siteName);
			plotFAS(rsFASFile);
			List<File> refFASFiles = new ArrayList<>();
			for (int i=0; i<numRefRuns; i++) {
				File subDir = new File(refDir, "run_"+i);
				refRD50Files.add(findRotD50File(subDir, siteName));
				refFASFiles.add(findFASFile(subDir, siteName));
			}
			plotMultiRotD50(refRD50Files, "Graves & Pitarka", rsRD50File, "RSQSim", "Event "+eventID+" "+siteName+" Spectra",
					rsDir, siteName+"_RotD50_compare_event_"+eventID, gmpeSpectra);
			plotMultiFAS(refFASFiles, "Graves & Pitarka", rsFASFile, "RSQSim", "Event "+eventID+" "+siteName+" Spectra",
					rsDir, siteName+"_FAS_compare_event_"+eventID);
		}
	}

}
