package scratch.kevin.sampling;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.sampling.SamplingMethod;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceMetric;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceSummary;

/** Builds paper-oriented plots from the compact convergence summary CSV files. */
public class HazardConvergencePlots {

	private static final Range Y_RANGE = new Range(1e-2, 2e1);
	private static final Range SIGNED_Y_RANGE = new Range(-1d, 1d);
	private static final String SOBOL_REFERENCE = HazardConvergenceCalcs.LOO_SOBOL_REFERENCE_NAME;
	private static final String POOLED_SOBOL_REFERENCE = HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME;
	private static final String MCS_REFERENCE = HazardConvergenceCalcs.MCS_REFERENCE_NAME;
	private static final SamplingMethod MCS = SamplingMethod.MONTE_CARLO;
	private static final SamplingMethod SOBOL = SamplingMethod.OWEN_SCRAMBLED_SOBOL;


	private static final Map<ConvergenceMetric, Color> METRIC_COLORS = Map.of(
			ConvergenceMetric.MEAN_HAZARD, Colors.tab_blue,
			ConvergenceMetric.STANDARD_DEVIATION, Colors.tab_orange,
			ConvergenceMetric.IQR, Colors.tab_green,
			ConvergenceMetric.CENTRAL_68_RANGE, Colors.tab_red,
			ConvergenceMetric.CENTRAL_95_RANGE, Colors.tab_purple);

	private static final Map<ConvergenceMetric, PlotSymbol> METRIC_SYMBOLS = Map.of(
			ConvergenceMetric.MEAN_HAZARD, PlotSymbol.FILLED_CIRCLE,
			ConvergenceMetric.STANDARD_DEVIATION, PlotSymbol.FILLED_INV_TRIANGLE,
			ConvergenceMetric.IQR, PlotSymbol.FILLED_TRIANGLE,
			ConvergenceMetric.CENTRAL_68_RANGE, PlotSymbol.FILLED_SQUARE,
			ConvergenceMetric.CENTRAL_95_RANGE, PlotSymbol.FILLED_DIAMOND);
	
	private static String getMethodName(SamplingMethod method) {
		if (method == SamplingMethod.OWEN_SCRAMBLED_SOBOL)
			return SamplingMethod.SOBOL.getShortName();
		return method.getShortName();
	}
	
	private static final Map<SamplingMethod, String> METHOD_FILE_PREFIXES;
	static {
		Map<SamplingMethod, String> prefixes = new HashMap<>();
		for (SamplingMethod method : SamplingMethod.values())
			prefixes.put(method, getMethodName(method).toLowerCase().replaceAll("-", "_"));
		METHOD_FILE_PREFIXES = prefixes;
	}

	public static void main(String[] args) throws IOException {
		File convergenceDir = new File(PaperPaths.FIGURES_DIR, "hazard_convergence");
		plotPeriod(new File(convergenceDir, "pga_two_in_50"), "PGA");
		plotPeriod(new File(convergenceDir, "1s_sa_two_in_50"), "1 s SA");
	}

	static void plotPeriod(File outputDir, String periodName) throws IOException {
		List<ReferenceSummary> references = loadReferenceSummaries(
				new File(outputDir, "reference_comparisons.csv"));
		List<RealizationPairSummary> realizationPairs = loadRealizationPairSummaries(
				new File(outputDir, "realization_pair_comparisons.csv"));

		for (SamplingMethod method : references.stream().map(ReferenceSummary::method).distinct().sorted().toList()) {
			for (boolean sobolPool : new boolean[] {true, false}) {
				plotReference(outputDir, references, method, sobolPool, ConvergenceSummary.MEAN_ABSOLUTE);
				plotReference(outputDir, references, method, sobolPool, ConvergenceSummary.MEAN_SIGNED);
//				// Retain the original Sobol convergence maximum plots as standalone figures.
//				if (method == SOBOL)
//					plotReference(outputDir, references, method, sobolPool, ConvergenceSummary.MAXIMUM_ABSOLUTE);
			}
		}
		for (boolean sobolPool : new boolean[] {true, false})
			plotMethodReference(outputDir, references, sobolPool, ConvergenceSummary.MEAN_ABSOLUTE,
					"Absolute difference (%)");
		plotRealizationPairs(outputDir, realizationPairs, ConvergenceSummary.MEAN_ABSOLUTE,
				"Absolute difference (%)");
	}

	private static String referenceFor(SamplingMethod method, boolean sobolPool) {
		return sobolPool ? (method == SOBOL ? SOBOL_REFERENCE : POOLED_SOBOL_REFERENCE)
				: (method == MCS ? HazardConvergenceCalcs.LOO_MCS_REFERENCE_NAME : MCS_REFERENCE);
	}

	private static boolean includeSummary(ConvergenceSummary actual, ConvergenceSummary requested) {
		return actual == requested || requested == ConvergenceSummary.MEAN_ABSOLUTE
				&& actual == ConvergenceSummary.MAXIMUM_ABSOLUTE;
	}

	private static void plotReference(File outputDir, List<ReferenceSummary> rows,
			SamplingMethod method, boolean sobolPool, ConvergenceSummary summary) throws IOException {
		List<ReferenceSummary> matching = rows.stream().filter(row -> row.method() == method
				&& row.reference().equals(referenceFor(method, sobolPool))
				&& includeSummary(row.spatialSummary(), summary)).toList();
		int[] counts = matching.stream().mapToInt(ReferenceSummary::sampleCount).distinct().sorted().toArray();
		if (counts.length < 2)
			return;
		String pool = sobolPool ? "pooled_sobol" : "pooled_mcs";
		String prefix = "convergence_"+METHOD_FILE_PREFIXES.get(method)+"_vs_"
				+pool+"_"+summaryPrefix(summary);
		String yLabel = summary == ConvergenceSummary.MEAN_SIGNED ? "Signed bias (%)"
				: summary == ConvergenceSummary.MEAN_ABSOLUTE ? "Absolute difference (%)"
						: "Maximum absolute difference (%)";
		writePlot(outputDir, prefix, getMethodName(method)+" versus "+(sobolPool ? "Sobol pool" : "MCS pool"),
				"Sample count", yLabel, counts,
				Arrays.stream(counts).mapToObj(Integer::toString).toArray(String[]::new), matching, summary);
	}

	static void plotMethodReference(File outputDir, List<ReferenceSummary> rows,
			boolean sobolConsensus, ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		int[] counts = rows.stream().mapToInt(ReferenceSummary::sampleCount).distinct().sorted().toArray();
		for (int count : counts)
			plotMethodReference(outputDir, rows, sobolConsensus, spatialSummary, yLabel, count);
	}

	private static void plotMethodReference(File outputDir, List<ReferenceSummary> rows,
			boolean sobolConsensus, ConvergenceSummary spatialSummary, String yLabel, int sampleCount) throws IOException {
		List<MethodSummary> matching = new ArrayList<>();
		List<String> labels = new ArrayList<>();
		for (SamplingMethod method : rows.stream().map(row -> row.method()).distinct().sorted().toList()) {
			String reference = referenceFor(method, sobolConsensus);
			int index = labels.size();
			for (ReferenceSummary row : rows) {
				if (row.method() == method && row.sampleCount() == sampleCount
						&& row.reference().equals(reference)
						&& includeSummary(row.spatialSummary(), spatialSummary)) {
					matching.add(new MethodSummary(index, row.metric(), row.spatialSummary(), row.realizations(),
							row.minimum(), row.median(), row.maximum()));
				}
			}
			if (matching.stream().anyMatch(row -> row.count() == index))
				labels.add(getMethodName(method));
		}
		if (matching.isEmpty())
			return;
		String referencePrefix = sobolConsensus ? "pooled_sobol" : "pooled_mcs";
		String title = sampleCount+" samples versus "
				+(sobolConsensus ? POOLED_SOBOL_REFERENCE : MCS_REFERENCE);
		writePlot(outputDir, "method_comparison_"+sampleCount+"_"+referencePrefix+"_"+summaryPrefix(spatialSummary),
				title, "Sampling method", yLabel,
				IntStream.range(0, labels.size()).toArray(), labels.toArray(String[]::new), matching, spatialSummary);
	}

	private static void plotRealizationPairs(File outputDir, List<RealizationPairSummary> rows,
			ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		int[] counts = rows.stream().mapToInt(RealizationPairSummary::sampleCount).distinct().sorted().toArray();
		for (int count : counts)
			plotRealizationPairs(outputDir, rows, spatialSummary, yLabel, count);
	}

	private static void plotRealizationPairs(File outputDir, List<RealizationPairSummary> rows,
			ConvergenceSummary spatialSummary, String yLabel, int sampleCount) throws IOException {
		List<MethodSummary> matching = new ArrayList<>();
		List<String> labels = new ArrayList<>();
		for (SamplingMethod method : rows.stream().map(row -> row.method()).distinct().sorted().toList()) {
			int index = labels.size();
			for (RealizationPairSummary row : rows) {
				if (row.method() == method && row.sampleCount() == sampleCount
						&& includeSummary(row.spatialSummary(), spatialSummary)) {
					matching.add(new MethodSummary(index, row.metric(), row.spatialSummary(), row.realizationPairs(),
							row.minimum(), row.median(), row.maximum()));
				}
			}
			if (matching.stream().anyMatch(row -> row.count() == index))
				labels.add(getMethodName(method));
		}
		if (matching.isEmpty())
			return;
		writePlot(outputDir, "method_comparison_"+sampleCount+"_realization_pairs_"+summaryPrefix(spatialSummary),
				"Differences between "+sampleCount+"-sample realizations", "Sampling method", yLabel,
				IntStream.range(0, labels.size()).toArray(), labels.toArray(String[]::new), matching, spatialSummary);
	}

	private static void writePlot(File outputDir, String prefix, String title, String xLabel, String yLabel,
			int[] counts, String[] countLabels, List<? extends SummaryRow> rows,
			ConvergenceSummary primary) throws IOException {
		boolean signed = primary == ConvergenceSummary.MEAN_SIGNED;
		Map<ConvergenceMetric, List<? extends SummaryRow>> byMetric = new LinkedHashMap<>();
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			List<? extends SummaryRow> metricRows = rows.stream().filter(row -> row.metric().equals(metric)
					&& row.spatialSummary() == primary).toList();
			if (!metricRows.isEmpty())
				byMetric.put(metric, metricRows);
		}

		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<XY_DataSet> maxFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> maxChars = new ArrayList<>();
		List<XY_DataSet> medianFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> medianChars = new ArrayList<>();
		if (signed) {
			ArbitrarilyDiscretizedFunc zero = new ArbitrarilyDiscretizedFunc();
			zero.set(-0.3d, 0d);
			zero.set(counts.length-0.7d, 0d);
			funcs.add(zero);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.GRAY));
		}
		for (Map.Entry<ConvergenceMetric, List<? extends SummaryRow>> entry : byMetric.entrySet()) {
			ArbitrarilyDiscretizedFunc median = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lower = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upper = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<counts.length; i++) {
				int count = counts[i];
				SummaryRow row = entry.getValue().stream().filter(candidate -> candidate.count() == count)
						.findFirst().orElseThrow();
				median.set((double)i, row.median());
				lower.set((double)i, row.minimum());
				upper.set((double)i, row.maximum());
			}
			Color color = METRIC_COLORS.get(entry.getKey());
			// Signed ranges overlap heavily, so only show the realization envelope for mean hazard.
			if (!signed || entry.getKey() == ConvergenceMetric.MEAN_HAZARD) {
				UncertainArbDiscFunc uncertainty = new UncertainArbDiscFunc(median, lower, upper);
				funcs.add(uncertainty);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
						new Color(color.getRed(), color.getGreen(), color.getBlue(), 70)));
			}
			median.setName(entry.getKey().shortLabel);
			medianFuncs.add(median);
			PlotSymbol sym = METRIC_SYMBOLS.get(entry.getKey());
			if (primary == ConvergenceSummary.MEAN_ABSOLUTE) {
				ArbitrarilyDiscretizedFunc worst = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<counts.length; i++) {
					int count = counts[i];
					rows.stream().filter(row -> row.metric() == entry.getKey() && row.count() == count
							&& row.spatialSummary() == ConvergenceSummary.MAXIMUM_ABSOLUTE)
							.findFirst().ifPresent(row -> worst.set((double)Arrays.binarySearch(counts, count), row.maximum()));
				}
				if (worst.size() > 0) {
					maxFuncs.add(worst); // Unnamed, so it adds no legend entry.
					maxChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 0.5f,
							PlotSymbol.getOutlineSymbol(sym), 3f, color));
				}
			}
			medianChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sym, 5f, color));
			PlotSymbol outlineSym = PlotSymbol.getOutlineSymbol(sym);
			if (outlineSym != null) {
				// add slightly darker outline overlay
				median = median.deepClone();
				median.setName(null);
				medianFuncs.add(median);
				medianChars.add(new PlotCurveCharacterstics(outlineSym, 5f, color.darker().darker()));
			}
		}
		// Put every envelope in the dataset first so all median lines render above all shading.
		funcs.addAll(maxFuncs);
		chars.addAll(maxChars);
		funcs.addAll(medianFuncs);
		chars.addAll(medianChars);

		PlotSpec plot = new PlotSpec(funcs, chars, title, xLabel, yLabel);
//		plot.setLegendVisible(true);
//		plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
//		plot.setLegendInset(RectangleAnchor.TOP);
//		plot.setLegendInset(RectangleAnchor.TOP, 0.5, 0.975, 0.9, false);
		plot.setLegendInset(RectangleAnchor.BOTTOM, 0.5, 0.025, 0.9, false);
		HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
		gp.getPlotPrefs().setPlotLabelFontSize(10);
		gp.getPlotPrefs().setLegendFontSize(8);
		gp.getPlotPrefs().setLegendLineLength(6d);
		gp.getPlotPrefs().setTickLabelFontSize(8);
//		gp.getPlotPrefs().setPlotPadding(new RectangleInsets(5, 0, 0, 12));
		Range xRange = counts.length == 1 ? new Range(-0.5, 0.5) : new Range(-0.2d, counts.length-0.8d);
//		Range xRange = counts.length == 1 ? new Range(-0.5, 0.5) : new Range(-0.3d, counts.length-0.7d);
		gp.drawGraphPanel(plot, false, !signed, xRange, signed ? SIGNED_Y_RANGE : Y_RANGE);
		PlotUtils.setXTick(gp, 1d);
		((NumberAxis)gp.getXAxis()).setNumberFormatOverride(categoryFormat(countLabels));
		PlotUtils.writePrintPlots(outputDir, prefix, gp, PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d,
				3.4d, 300, true, true, false);
	}

	private static NumberFormat categoryFormat(String[] labels) {
		return new NumberFormat() {
			@Override public StringBuffer format(double value, StringBuffer buffer, FieldPosition pos) {
				int index = (int)Math.round(value);
				if (index >= 0 && index < labels.length && Math.abs(value-index) < 1e-6)
					buffer.append(labels[index]);
				return buffer;
			}
			@Override public StringBuffer format(long value, StringBuffer buffer, FieldPosition pos) {
				return format((double)value, buffer, pos);
			}
			@Override public Number parse(String source, ParsePosition pos) {
				pos.setErrorIndex(pos.getIndex());
				return null;
			}
		};
	}

	private static List<ReferenceSummary> loadReferenceSummaries(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		Map<ReferenceGroup, List<Double>> groups = new LinkedHashMap<>();
		for (int row=1; row<csv.getNumRows(); row++) {
			SamplingMethod method = SamplingMethod.valueOf(csv.get(row, 1));
			int sampleCount = Integer.parseInt(csv.get(row, 4));
			String reference = csv.get(row, 5);
			ConvergenceMetric metric = ConvergenceMetric.fromLabel(csv.get(row, 7));
			addValue(groups, new ReferenceGroup(method, sampleCount, reference, metric,
					ConvergenceSummary.MEAN_SIGNED), Double.parseDouble(csv.get(row, 8)));
			addValue(groups, new ReferenceGroup(method, sampleCount, reference, metric,
					ConvergenceSummary.MEAN_ABSOLUTE), Double.parseDouble(csv.get(row, 9)));
			addValue(groups, new ReferenceGroup(method, sampleCount, reference, metric,
					ConvergenceSummary.MAXIMUM_ABSOLUTE), Double.parseDouble(csv.get(row, 11)));
		}
		List<ReferenceSummary> rows = new ArrayList<>();
		for (Map.Entry<ReferenceGroup, List<Double>> entry : groups.entrySet()) {
			ReferenceGroup group = entry.getKey();
			double[] values = entry.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			rows.add(new ReferenceSummary(group.method(), group.sampleCount(), group.reference(), group.metric(),
					group.spatialSummary(), values.length, StatUtils.min(values),
					StatUtils.percentile(values, 50d), StatUtils.max(values)));
		}
		return rows;
	}

	private static List<RealizationPairSummary> loadRealizationPairSummaries(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		Map<RealizationPairGroup, List<Double>> groups = new LinkedHashMap<>();
		for (int row=1; row<csv.getNumRows(); row++) {
			SamplingMethod method = SamplingMethod.valueOf(csv.get(row, 0));
			int sampleCount = Integer.parseInt(csv.get(row, 1));
			ConvergenceMetric metric = ConvergenceMetric.fromLabel(csv.get(row, 6));
			addValue(groups, new RealizationPairGroup(method, sampleCount, metric,
					ConvergenceSummary.MEAN_ABSOLUTE), Double.parseDouble(csv.get(row, 8)));
			addValue(groups, new RealizationPairGroup(method, sampleCount, metric,
					ConvergenceSummary.MAXIMUM_ABSOLUTE), Double.parseDouble(csv.get(row, 10)));
		}
		List<RealizationPairSummary> rows = new ArrayList<>();
		for (Map.Entry<RealizationPairGroup, List<Double>> entry : groups.entrySet()) {
			RealizationPairGroup group = entry.getKey();
			double[] values = entry.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			rows.add(new RealizationPairSummary(group.method(), group.sampleCount(), group.metric(),
					group.spatialSummary(), values.length, StatUtils.min(values),
					StatUtils.percentile(values, 50d), StatUtils.max(values)));
		}
		return rows;
	}

	private static <K> void addValue(Map<K, List<Double>> groups, K key, double value) {
		groups.computeIfAbsent(key, unused -> new ArrayList<>()).add(value);
	}

	private static String summaryPrefix(ConvergenceSummary summary) {
		return switch (summary) {
		case MEAN_SIGNED -> "signed_bias";
		case MEAN_ABSOLUTE -> "mean_abs";
		case P95_ABSOLUTE -> "p95_abs";
		case MAXIMUM_ABSOLUTE -> "max_abs";
		};
	}

	interface SummaryRow {
		int count();
		ConvergenceSummary spatialSummary();
		ConvergenceMetric metric();
		double minimum();
		double median();
		double maximum();
	}

	record ReferenceSummary(SamplingMethod method, int sampleCount, String reference, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary,
			int realizations, double minimum, double median, double maximum) implements SummaryRow {
		@Override public int count() { return sampleCount; }
	}

	private record MethodSummary(int methodIndex, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary, int realizations,
			double minimum, double median, double maximum) implements SummaryRow {
		@Override public int count() { return methodIndex; }
	}


	private record RealizationPairSummary(SamplingMethod method, int sampleCount, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary, int realizationPairs,
			double minimum, double median, double maximum) {}

	private record ReferenceGroup(SamplingMethod method, int sampleCount, String reference,
			ConvergenceMetric metric,
			ConvergenceSummary spatialSummary) {}

	private record RealizationPairGroup(SamplingMethod method, int sampleCount,
			ConvergenceMetric metric, ConvergenceSummary spatialSummary) {}
}
