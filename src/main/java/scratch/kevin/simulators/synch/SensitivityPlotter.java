package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.synch.MPJSynchSensTest.SensTestType;

public class SensitivityPlotter {
	
	private enum Fault {
		COACHELLA("Coachella"),
		CARRIZO("Carrizo"),
		CHOLAME("Cholame"),
		MOJAVE("Mojave"),
		GARLOCK("Garlock"),
		SAN_JACINTO("San Jacinto");
		
		private String name;
		
		private Fault(String name) {
			this.name = name;
		}
	}
	
	private SensTestType testType;
	private Table<Fault, Fault, SensitivityResult> results;
	
	public SensitivityPlotter(File dir, SensTestType testType) throws IOException {
		this.testType = testType;
		
		File subDir = new File(dir, testType.name());
		Preconditions.checkState(subDir.exists(), "Doesn't exist: "+subDir.getAbsolutePath());
		
		int maxTrials = -1;
		
		for (File file : subDir.listFiles()) {
			String name = file.getName();
			if (!name.endsWith(".csv") || !name.contains("trials"))
				continue;
			
			// parse num trials
			String subName = name.substring(name.lastIndexOf("_")+1, name.indexOf("trials"));
			int numTrials = Integer.parseInt(subName);
			if (numTrials > maxTrials) {
				results = HashBasedTable.create();
				maxTrials = numTrials;
			} else if (numTrials < maxTrials) {
				continue;
			}
			
			Fault[] faults = getFaults(name);
			CSVFile<String> csv = CSVFile.readFile(file, true);
			results.put(faults[0], faults[1], new SensitivityResult(faults[0], faults[1], csv));
		}
	}
	
	private static Fault[] getFaults(String fileName) {
		fileName = fileName.replaceAll("_", " ");
		Fault[] ret = new Fault[2];
		
		for (Fault fault : Fault.values()) {
			if (fileName.contains(fault.name)) {
				if (ret[0] == null)
					ret[0] = fault;
				else if (ret[1] == null)
					ret[1] = fault;
				else
					throw new IllegalStateException("Name contains 3 faults???");
			}
		}
		Preconditions.checkState(ret[0] != null, "Couldn't parse names for: "+fileName);
		Preconditions.checkState(ret[1] != null, "Couldn't parse names for: "+fileName);
		return ret;
	}
	
	private static DecimalFormat catDF = new DecimalFormat("0.00");
	
	private class SensitivityResult {
		private Fault fault1;
		private Fault fault2;
		private List<Double> values;
		private List<Double> gBars;
		
		private List<List<Double>> trials;
		private List<Double> trialMeans;
		private List<Range> confidenceInvervals;
		
		private DiscretizedFunc gBarFunc;
		private UncertainArbDiscDataset trialsFunc;
		
		public SensitivityResult(Fault fault1, Fault fault2, CSVFile<String> csv) {
			this.fault1 = fault1;
			this.fault2 = fault2;
			
			values = praseLine(csv.getLine(0));
			gBars = praseLine(csv.getLine(1));
			
			if (csv.getNumRows() > 2) {
				// we have trials
				trials = Lists.newArrayList();
				for (int row=2; row<csv.getNumRows(); row++)
					trials.add(praseLine(csv.getLine(row)));
				
				trialMeans = Lists.newArrayList();
				confidenceInvervals = Lists.newArrayList();
				
				for (int i=0; i<values.size(); i++) {
					double[] vals = new double[trials.size()];
					for (int j=0; j<trials.size(); j++)
						vals[j] = trials.get(j).get(i);
					double mean = StatUtils.mean(vals);
					double upper = StatUtils.percentile(vals, 97.5);
					double lower = StatUtils.percentile(vals, 2.5);
					if (lower == 0)
						lower = Math.exp(-10);
					if (upper == 0)
						upper = Math.exp(10);
					if (mean == 0)
						mean = Math.exp(-10);
					trialMeans.add(mean);
					confidenceInvervals.add(new Range(lower, upper));
				}
			}
			
			if (testType == SensTestType.CAT_FRACT) {
				// special case
				
				Map<String, List<Double>> mappedVals = Maps.newHashMap();
				Map<String, List<Double>> mappedTrialVals = Maps.newHashMap();
				
				for (int i=0; i<values.size(); i++) {
					double val = values.get(i);
					if (val > 1d)
						// remove sequence identifier
						val = val - Math.floor(val);
					String key = catDF.format(val);
					List<Double> mapped = mappedVals.get(key);
					if (mapped == null) {
						mapped = Lists.newArrayList();
						mappedVals.put(key, mapped);
					}
					mapped.add(gBars.get(i));
					
					if (trials != null) {
						mapped = mappedTrialVals.get(key);
						if (mapped == null) {
							mapped = Lists.newArrayList();
							mappedTrialVals.put(key, mapped);
						}
						for (List<Double> trialVals : trials)
							mapped.add(trialVals.get(i));
					}
				}
				
				// gBarFunc is itself uncertain
				DiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
				meanFunc.setName(testType.getDisplayName()+" sensitivity, "
						+fault1.name+" vs "+fault2.name);
				DiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
				DiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
				
				for (String key : mappedVals.keySet()) {
					double val = Double.parseDouble(key);
					List<Double> vals = mappedVals.get(key);
					int expectedNum = (int)(1d/val);
					Preconditions.checkState(expectedNum == vals.size(),
							"Expected "+expectedNum+" got "+vals.size()+" for val="+val);
					double[] valArray = Doubles.toArray(vals);
					meanFunc.set(val, Math.log(StatUtils.mean(valArray)));
					lowerFunc.set(val, Math.log(StatUtils.min(valArray)));
					upperFunc.set(val, Math.log(StatUtils.max(valArray)));
				}
				
				gBarFunc = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
				
				if (trials != null) {
					meanFunc = new ArbitrarilyDiscretizedFunc();
					meanFunc.setName(testType.getDisplayName()+" sensitivity random trials, "
							+fault1.name+" vs "+fault2.name);
					lowerFunc = new ArbitrarilyDiscretizedFunc();
					upperFunc = new ArbitrarilyDiscretizedFunc();
					
					for (String key : mappedTrialVals.keySet()) {
						double val = Double.parseDouble(key);
						List<Double> vals = mappedTrialVals.get(key);
						int expectedNum = (int)(1d/val) * trials.size();
						Preconditions.checkState(expectedNum == vals.size(),
								"Expected "+expectedNum+" got "+vals.size()+" for val="+val);
						double[] valArray = Doubles.toArray(vals);
						meanFunc.set(val, Math.log(StatUtils.mean(valArray)));
						lowerFunc.set(val, Math.log(StatUtils.percentile(valArray, 2.5)));
						upperFunc.set(val, Math.log(StatUtils.percentile(valArray, 97.5)));
					}
					
					trialsFunc = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
				}
			} else {
				gBarFunc = new ArbitrarilyDiscretizedFunc();
				gBarFunc.setName(testType.getDisplayName()+" sensitivity, "+fault1.name+" vs "+fault2.name);
				
				for (int i=0; i<values.size(); i++) {
					gBarFunc.set(values.get(i), Math.log(gBars.get(i)));
				}
				
				if (trials != null) {
					DiscretizedFunc meanFunc = new ArbitrarilyDiscretizedFunc();
					meanFunc.setName(testType.getDisplayName()+" sensitivity random trials, "
							+fault1.name+" vs "+fault2.name);
					DiscretizedFunc lowerFunc = new ArbitrarilyDiscretizedFunc();
					DiscretizedFunc upperFunc = new ArbitrarilyDiscretizedFunc();
					
					for (int i=0; i<values.size(); i++) {
						double x = values.get(i);
						meanFunc.set(x, Math.log(trialMeans.get(i)));
						lowerFunc.set(x, Math.log(confidenceInvervals.get(i).getLowerBound()));
						upperFunc.set(x, Math.log(confidenceInvervals.get(i).getUpperBound()));
					}
					
					trialsFunc = new UncertainArbDiscDataset(meanFunc, lowerFunc, upperFunc);
				}
			}
		}
		
	}
	
	private static List<Double> praseLine(List<String> line) {
		// first value is header
		List<Double> vals = Lists.newArrayList();
		for (int i=1; i<line.size(); i++)
			vals.add(Double.parseDouble(line.get(i)));
		return vals;
	}
	
	public void plotSinglePairing(Fault fault1, Fault fault2, File outputDir) throws IOException {
		SensitivityResult result = results.get(fault1, fault2);
		if (result == null)
			result = results.get(fault2, fault1);
		Preconditions.checkNotNull(result);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		if (result.trials != null) {
			// first uncertainties
			funcs.add(result.trialsFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.BLUE));
			// then plot mean
			funcs.add(result.trialsFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLUE));
		}
		
		// now main func
		if (result.gBarFunc instanceof UncertainArbDiscDataset) {
			funcs.add(result.gBarFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.RED));
		}
		funcs.add(result.gBarFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
		
		String title = testType.getDisplayName()+" Sensitivity Study: "+fault1.name+" vs "+fault2.name;
		String xAxisLabel = testType.getDisplayName();
		String yAxisLabel = "Ln(Gbar)";
		
		double maxVal = 0d;
		for (DiscretizedFunc func : funcs) {
			double funcMax = Math.max(Math.abs(func.getMaxY()), Math.abs(func.getMinY()));
			if (func instanceof UncertainArbDiscDataset) {
				UncertainArbDiscDataset ufunc = (UncertainArbDiscDataset)func;
				funcMax = Math.max(funcMax, Math.max(Math.abs(ufunc.getUpperMaxY()), Math.abs(ufunc.getUpperMinY())));
				funcMax = Math.max(funcMax, Math.max(Math.abs(ufunc.getLowerMaxY()), Math.abs(ufunc.getLowerMinY())));
			}
			if (funcMax > maxVal)
				maxVal = funcMax;
		}
		Preconditions.checkState(Doubles.isFinite(maxVal));
		double plotMax = 2;
		while (maxVal > plotMax)
			plotMax += 1d;
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(18);

		Range xRange = new Range(result.gBarFunc.getMinX(), result.gBarFunc.getMaxX());
		Range yRange = new Range(-plotMax, plotMax);
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		outputDir = new File(outputDir, testType.name());
		outputDir = new File(outputDir, "plots");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		File outputFile = new File(outputDir, testType.name()+"_"
				+fault1.name.replaceAll(" ", "_")+"_"+fault2.name.replaceAll(" ", "_"));
		
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		gp.getChartPanel().setSize(500, 400);
		File smallPlotDir = new File(outputDir, "small_plots");
		if (!smallPlotDir.exists())
			smallPlotDir.mkdir();
		outputFile = new File(smallPlotDir, outputFile.getName());
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
	}
	
	public void plotAllSinglePairings(File dir) throws IOException {
		for (Cell<Fault, Fault, SensitivityResult> cell : results.cellSet())
			plotSinglePairing(cell.getRowKey(), cell.getColumnKey(), dir);
	}
	
	public void plotAllOnePlot(File dir) throws IOException {
		List<Color> defaultColors = GraphWindow.generateDefaultColors();
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<DiscretizedFunc> uncertainFuncs = Lists.newArrayList();
		List<PlotCurveCharacterstics> uncertainChars = Lists.newArrayList();
		
		List<Cell<Fault, Fault, SensitivityResult>> cells = Lists.newArrayList(results.cellSet());
		// sort for repeatable coloring
		Collections.sort(cells, new Comparator<Cell<Fault, Fault, SensitivityResult>>() {

			@Override
			public int compare(Cell<Fault, Fault, SensitivityResult> o1,
					Cell<Fault, Fault, SensitivityResult> o2) {
				int ret = o1.getRowKey().name.compareTo(o2.getRowKey().name);
				if (ret != 0)
					return ret;
				return o1.getColumnKey().name.compareTo(o2.getColumnKey().name);
			}
		});
		
		int colorIndex = 0;
		for (Cell<Fault, Fault, SensitivityResult> cell : cells) {
			if (colorIndex == defaultColors.size())
				colorIndex = 0;
			
			Color color = defaultColors.get(colorIndex);
			
			DiscretizedFunc func = cell.getValue().gBarFunc;
			if (func instanceof UncertainArbDiscDataset) {
				uncertainFuncs.add(func);
				uncertainChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, color));
			}
			funcs.add(func);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, color));
			
			colorIndex++;
		}
		
		funcs.addAll(0, uncertainFuncs);
		chars.addAll(0, uncertainChars);
		
		String title = testType.getDisplayName()+" Sensitivity Study, All Faults";
		String xAxisLabel = testType.getDisplayName();
		String yAxisLabel = "Ln(Gbar)";
		
		double maxVal = 0d;
		for (DiscretizedFunc func : funcs) {
			double funcMax = Math.max(Math.abs(func.getMaxY()), Math.abs(func.getMinY()));
			if (func instanceof UncertainArbDiscDataset) {
				UncertainArbDiscDataset ufunc = (UncertainArbDiscDataset)func;
				funcMax = Math.max(funcMax, Math.max(Math.abs(ufunc.getUpperMaxY()), Math.abs(ufunc.getUpperMinY())));
				funcMax = Math.max(funcMax, Math.max(Math.abs(ufunc.getLowerMaxY()), Math.abs(ufunc.getLowerMinY())));
			}
			if (funcMax > maxVal)
				maxVal = funcMax;
		}
		double plotMax = 2;
		while (maxVal > plotMax)
			plotMax += 1d;
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(18);

		Range xRange = new Range(funcs.get(0).getMinX(), funcs.get(0).getMaxX());
		Range yRange = new Range(-plotMax, plotMax);
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		File outputDir = new File(dir, testType.name());
		outputDir = new File(outputDir, "plots");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		File outputFile = new File(outputDir, testType.name()+"_ALL");
		
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
		gp.getChartPanel().setSize(500, 400);
		File smallPlotDir = new File(outputDir, "small_plots");
		if (!smallPlotDir.exists())
			smallPlotDir.mkdir();
		outputFile = new File(smallPlotDir, testType.name()+"_ALL");
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsPDF(outputFile.getAbsolutePath()+".pdf");
	}
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators/synch/weight_CATALOG_G_indep_shiftLag/sens_tests/");
		SensTestType testType = SensTestType.ELEM_VOL;
		
		SensitivityPlotter plot = new SensitivityPlotter(dir, testType);
		plot.plotAllSinglePairings(dir);
		plot.plotAllOnePlot(dir);
	}

}
