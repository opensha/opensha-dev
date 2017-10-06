package scratch.kevin.bbp;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
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

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SpectraPlotter {
	
	private static DiscretizedFunc[] loadAll(File file) throws IOException {
		DiscretizedFunc[] ret = null;
		String[] names = null;
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
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
			for (int i=0; i<ret.length; i++)
				ret[i].set(x, Double.parseDouble(split[i+1]));
		}
		return ret;
	}
	
	public static DiscretizedFunc loadRotD50(File file) throws IOException {
		DiscretizedFunc[] ret = loadAll(file);
		return ret[ret.length-1];
	}

	public static void plotRotD50(File file) throws IOException {
		plotRotD50(file, null, null, null);
	}
	
	public static void plotRotD50(File file, File outputDir, String prefix, UncertainArbDiscDataset[] gmpes)
			throws IOException {
		if (file.isDirectory()) {
			for (File sub : file.listFiles()) {
				if (sub.getName().endsWith(".rd50"))
					plotRotD50(sub, outputDir, prefix, gmpes);
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
				outputDir, prefix, gmpes);
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
				true, true, null, outputDir, prefix, null);
	}
	
	private static Color[] GMPE_COLORS = { Color.RED.brighter(), Color.YELLOW.brighter(),
			Color.CYAN.brighter(), Color.MAGENTA.brighter() };
	
	private static void plotSpectra(DiscretizedFunc[] funcArray, String title, String xAxisLable, String yAxisLabel,
			boolean xLog, boolean yLog, Range yRange, File outputDir, String prefix, UncertainArbDiscDataset[] gmpes) throws IOException {
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
			if (i == funcArray.length -1) {
				// RotD50
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
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
		plotMultiSpectra(refFiles, refName, dataFile, dataName, title, outputDir, prefix, true, gmpes);
	}
	
	public static void plotMultiFAS(List<File> refFiles, String refName, File dataFile, String dataName, String title,
			File outputDir, String prefix) throws IOException {
		plotMultiSpectra(refFiles, refName, dataFile, dataName, title, outputDir, prefix, false, null);
	}
	
	private static void plotMultiSpectra(List<File> refFiles, String refName, File dataFile, String dataName, String title,
			File outputDir, String prefix, boolean rotD50, UncertainArbDiscDataset[] gmpes) throws IOException {
		XY_DataSetList refFuncs = new XY_DataSetList();
		List<Double> relativeWts = new ArrayList<>();
		for (File refFile : refFiles) {
			if (rotD50)
				refFuncs.add(loadRotD50(refFile));
			else
				refFuncs.add(loadFAS(refFile));
			relativeWts.add(1d);
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
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
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
			for (Point2D pt : lowerFunc)
				if (pt.getY() > 0)
					minY = Math.min(minY, pt.getY());
		}
		
		funcs.add(refMean);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		
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
		
		DiscretizedFunc dataFunc;
		if (rotD50)
			dataFunc = loadRotD50(dataFile);
		else
			dataFunc = loadFAS(dataFile);
		for (Point2D pt : dataFunc)
			if (pt.getY() > 0)
				minY = Math.min(minY, pt.getY());
		dataFunc.setName(dataName);
		funcs.add(dataFunc);
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
		
		minY = Math.pow(10, Math.floor(Math.log10(minY)));
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
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
	
	public static UncertainArbDiscDataset calcGMPE_RotD50(EqkRupture rupture, Location loc, VelocityModel vm, ScalarIMR gmpe) {
		Site site = new Site(loc);
		
		Vs30_Param vs30Param = new Vs30_Param(vm.getVs30());
		vs30Param.setValueAsDefault();
		site.addParameter(vs30Param);
		
		Vs30_TypeParam vs30TypeParam = new Vs30_TypeParam();
		vs30TypeParam.setValue(Vs30_TypeParam.VS30_TYPE_MEASURED); // TODO
		site.addParameter(vs30TypeParam);
		
		DepthTo1pt0kmPerSecParam z10Param = new DepthTo1pt0kmPerSecParam();
		z10Param.setValue(null);
		site.addParameter(z10Param);
		
		DepthTo2pt5kmPerSecParam z25Param = new DepthTo2pt5kmPerSecParam();
		z25Param.setValue(null);
		site.addParameter(z25Param);
		
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
//		plotRotD50(new File("/home/kevin/bbp/bbp_data/outdata/6462338"));
//		plotRotD50(new File("/home/kevin/bbp/bbp_data/outdata/6450680"));
//		plotRotD50(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/"
//				+ "event_4099020_0.05s_ADJ_VEL_seis"));
//		plotRotD50(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/event_4099020_0.05s_ADJ_VEL_bbp"));
//		plotFAS(new File("/data/kevin/simulators/catalogs/JG_UCERF3_millionElement/event_srfs/event_4099020_0.05s_ADJ_VEL_bbp"));
//		System.exit(0);
		
//		ScalarIMR[] gmpes = null;
		ScalarIMR[] gmpes = { AttenRelRef.ASK_2014.instance(null), AttenRelRef.BSSA_2014.instance(null),
				AttenRelRef.CB_2014.instance(null), AttenRelRef.CY_2014.instance(null) };
		
		String[] siteNames = { "USC", "SBSM" };
		Location[] locs = { new Location(34.0192, -118.286), new Location(34.064986, -117.29201) };
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
		File refDir = new File("/home/kevin/bbp/parallel/2017_10_04-rundir2194_long-event136704-dx1.16-noHF/results_orig");
		int eventID = 136704;
		
//		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
//		File refDir = new File("/home/kevin/bbp/parallel/2017_10_04-JG_UCERF3_millionElement-event4099020-dx0.48-noHF/results");
//		int eventID = 4099020;
		
		File srfDir = new File(catalog.getCatalogDir(), "event_srfs");
		File rsDir = new File(srfDir, "event_"+eventID+"_0.05s_ADJ_VEL_bbp");
		
		int numRefRuns = 200;
		VelocityModel vm = VelocityModel.LA_BASIN;
		double minFractForInclusion = 0.2;
		
		EqkRupture gmpeRup = null;
		if (gmpes != null) {
			for (ScalarIMR gmpe : gmpes)
				gmpe.setParamDefaults();
			System.out.println("Loading event...");
			RSQSimEvent event = catalog.loadEventByID(eventID);
			gmpeRup = catalog.getGMPE_Rupture(event, minFractForInclusion);
			System.out.println("DONE");
		}
		
		for (int s=0; s<siteNames.length; s++) {
			String siteName = siteNames[s];
			Location loc = locs[s];
			File rsRD50File = findRotD50File(rsDir, siteName);
			UncertainArbDiscDataset[] gmpeSpectra = null;
			if (gmpes != null) {
				gmpeSpectra = new UncertainArbDiscDataset[gmpes.length];
				for (int i=0; i<gmpes.length; i++) {
					System.out.println("Calculating spectra for "+gmpes[i].getShortName());
					gmpeSpectra[i] = calcGMPE_RotD50(gmpeRup, loc, vm, gmpes[i]);
				}
				System.out.println("DONE spectra");
			}
			plotRotD50(rsRD50File, null, null, gmpeSpectra);
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
