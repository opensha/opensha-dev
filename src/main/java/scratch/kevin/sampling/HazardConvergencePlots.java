package scratch.kevin.sampling;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.ui.RectangleAnchor;
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

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceMetric;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceSummary;

/** Builds paper-oriented plots from the compact convergence summary CSV files. */
public class HazardConvergencePlots {

	private static final Range Y_RANGE = new Range(1e-2, 2e1);
	private static final String SOBOL_REFERENCE = HazardConvergenceCalcs.LOO_SOBOL_REFERENCE_NAME;
	private static final String POOLED_SOBOL_REFERENCE = HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME;
	private static final String MCS_REFERENCE = HazardConvergenceCalcs.MCS_REFERENCE_NAME;
	private static final SamplingMethod SOBOL = SamplingMethod.OWEN_SCRAMBLED_SOBOL;
	private static final SamplingMethod PAIRWISE_LHS = SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE;
	private static final List<SamplingMethod> COMPARISON_METHODS = List.of(SOBOL, PAIRWISE_LHS);

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

	public static void main(String[] args) throws IOException {
		File convergenceDir = new File(PaperPaths.FIGURES_DIR, "hazard_convergence/sobol_convergence");
		plotPeriod(new File(convergenceDir, "pga_two_in_50"), "PGA");
		plotPeriod(new File(convergenceDir, "1s_sa_two_in_50"), "1 s SA");
	}

	static void plotPeriod(File outputDir, String periodName) throws IOException {
		List<ReferenceSummary> references = loadReferenceSummaries(
				new File(outputDir, "reference_comparisons.csv"));
		List<DoublingSummary> doublings = loadDoublingSummaries(
				new File(outputDir, "paired_doubling_comparisons.csv"));
		List<RealizationPairSummary> realizationPairs = loadRealizationPairSummaries(
				new File(outputDir, "realization_pair_comparisons.csv"));

		plotReference(outputDir, periodName, references, SOBOL_REFERENCE, "sobol_consensus",
				ConvergenceSummary.MEAN_ABSOLUTE, "Spatial mean absolute difference (%)");
		plotReference(outputDir, periodName, references, SOBOL_REFERENCE, "sobol_consensus",
				ConvergenceSummary.MAXIMUM_ABSOLUTE, "Maximum absolute difference (%)");
		plotReference(outputDir, periodName, references, MCS_REFERENCE, "mcs_reference",
				ConvergenceSummary.MEAN_ABSOLUTE, "Spatial mean absolute difference (%)");
		plotReference(outputDir, periodName, references, MCS_REFERENCE, "mcs_reference",
				ConvergenceSummary.MAXIMUM_ABSOLUTE, "Maximum absolute difference (%)");
		plotDoubling(outputDir, periodName, doublings, ConvergenceSummary.MEAN_ABSOLUTE,
				"Spatial mean absolute difference (%)");
		plotDoubling(outputDir, periodName, doublings, ConvergenceSummary.MAXIMUM_ABSOLUTE,
				"Maximum absolute difference (%)");
		plotMethodReference(outputDir, references, false, ConvergenceSummary.MEAN_ABSOLUTE,
				"Spatial mean absolute difference (%)");
		plotMethodReference(outputDir, references, false, ConvergenceSummary.MAXIMUM_ABSOLUTE,
				"Maximum absolute difference (%)");
		plotMethodReference(outputDir, references, true, ConvergenceSummary.MEAN_ABSOLUTE,
				"Spatial mean absolute difference (%)");
		plotMethodReference(outputDir, references, true, ConvergenceSummary.MAXIMUM_ABSOLUTE,
				"Maximum absolute difference (%)");
		plotRealizationPairs(outputDir, realizationPairs, ConvergenceSummary.MEAN_ABSOLUTE,
				"Spatial mean absolute difference (%)");
		plotRealizationPairs(outputDir, realizationPairs, ConvergenceSummary.MAXIMUM_ABSOLUTE,
				"Maximum absolute difference (%)");
	}

	private static void plotReference(File outputDir, String periodName, List<ReferenceSummary> rows,
			String reference, String referencePrefix, ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		List<ReferenceSummary> matching = rows.stream().filter(row -> row.method() == SOBOL
				&& row.reference().equals(reference)
				&& row.spatialSummary().equals(spatialSummary)).toList();
		Preconditions.checkState(!matching.isEmpty(), "No rows for %s, %s", reference, spatialSummary);
		int[] counts = matching.stream().mapToInt(ReferenceSummary::sampleCount).distinct().sorted().toArray();
//		String title = periodName+", "+(reference.equals(SOBOL_REFERENCE)
//				? "pooled Sobol consensus" : MCS_REFERENCE+" reference");
		String title = reference.equals(SOBOL_REFERENCE) ? "Pooled Sobol consensus" : MCS_REFERENCE+" reference";
		String prefix = "convergence_"+referencePrefix+"_"+summaryPrefix(spatialSummary);
		writePlot(outputDir, prefix, title, "Sample count", yLabel, counts,
				Arrays.stream(counts).mapToObj(Integer::toString).toArray(String[]::new), matching);
	}

	private static void plotMethodReference(File outputDir, List<ReferenceSummary> rows,
			boolean sobolConsensus, ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		List<MethodSummary> matching = new ArrayList<>();
		for (int m=0; m<COMPARISON_METHODS.size(); m++) {
			SamplingMethod method = COMPARISON_METHODS.get(m);
			String reference = sobolConsensus && method == SOBOL ? SOBOL_REFERENCE
					: sobolConsensus ? POOLED_SOBOL_REFERENCE : MCS_REFERENCE;
			for (ReferenceSummary row : rows) {
				if (row.method() == method && row.sampleCount() == 4096 && row.reference().equals(reference)
						&& row.spatialSummary() == spatialSummary) {
					matching.add(new MethodSummary(m, row.metric(), row.spatialSummary(), row.realizations(),
							row.minimum(), row.median(), row.maximum()));
				}
			}
		}
		Preconditions.checkState(!matching.isEmpty(), "No 4096-sample method comparisons");
		String referencePrefix = sobolConsensus ? "sobol_consensus" : "mcs_reference";
		String title = "4096 samples versus "+(sobolConsensus ? "pooled Sobol consensus" : MCS_REFERENCE);
		writePlot(outputDir, "method_comparison_"+referencePrefix+"_"+summaryPrefix(spatialSummary),
				title, "Sampling method", yLabel, new int[] { 0, 1 },
				COMPARISON_METHODS.stream().map(SamplingMethod::getShortName).toArray(String[]::new), matching);
	}

	private static void plotRealizationPairs(File outputDir, List<RealizationPairSummary> rows,
			ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		List<MethodSummary> matching = new ArrayList<>();
		for (int m=0; m<COMPARISON_METHODS.size(); m++) {
			SamplingMethod method = COMPARISON_METHODS.get(m);
			for (RealizationPairSummary row : rows) {
				if (row.method() == method && row.sampleCount() == 4096
						&& row.spatialSummary() == spatialSummary) {
					matching.add(new MethodSummary(m, row.metric(), row.spatialSummary(), row.realizationPairs(),
							row.minimum(), row.median(), row.maximum()));
				}
			}
		}
		Preconditions.checkState(!matching.isEmpty(), "No 4096-sample realization-pair comparisons");
		writePlot(outputDir, "method_comparison_realization_pairs_"+summaryPrefix(spatialSummary),
				"Differences between 4096-sample realizations", "Sampling method", yLabel,
				new int[] { 0, 1 },
				COMPARISON_METHODS.stream().map(SamplingMethod::getShortName).toArray(String[]::new), matching);
	}

	private static void plotDoubling(File outputDir, String periodName, List<DoublingSummary> rows,
			ConvergenceSummary spatialSummary, String yLabel) throws IOException {
		List<DoublingSummary> matching = rows.stream()
				.filter(row -> row.spatialSummary().equals(spatialSummary)).toList();
		Preconditions.checkState(!matching.isEmpty(), "No doubling rows for %s", spatialSummary);
		int[] upperCounts = matching.stream().mapToInt(DoublingSummary::upperCount).distinct().sorted().toArray();
		String[] labels = new String[upperCounts.length];
		for (int i=0; i<labels.length; i++)
			labels[i] = (upperCounts[i]/2)+"-"+upperCounts[i];
		List<SummaryRow> generic = new ArrayList<>(matching);
//		String title = periodName+", paired sample-count increases";
		String title = "Paired sample-count increases";
		writePlot(outputDir, "convergence_paired_doubling_"+summaryPrefix(spatialSummary),
				title, "Sample-count increase", yLabel,
				upperCounts, labels, generic);
	}

	private static void writePlot(File outputDir, String prefix, String title, String xLabel, String yLabel,
			int[] counts, String[] countLabels, List<? extends SummaryRow> rows) throws IOException {
		Map<ConvergenceMetric, List<? extends SummaryRow>> byMetric = new LinkedHashMap<>();
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			List<? extends SummaryRow> metricRows = rows.stream().filter(row -> row.metric().equals(metric)).toList();
			if (!metricRows.isEmpty())
				byMetric.put(metric, metricRows);
		}

		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<XY_DataSet> medianFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> medianChars = new ArrayList<>();
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
			UncertainArbDiscFunc uncertainty = new UncertainArbDiscFunc(median, lower, upper);
			funcs.add(uncertainty);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					new Color(color.getRed(), color.getGreen(), color.getBlue(), 70)));
			median.setName(entry.getKey().shortLabel);
			medianFuncs.add(median);
			PlotSymbol sym = METRIC_SYMBOLS.get(entry.getKey());
			medianChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sym, 5f, color));
		}
		// Put every envelope in the dataset first so all median lines render above all shading.
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
		Range xRange = new Range(-0.2d, counts.length-0.8d);
		gp.drawGraphPanel(plot, false, true, xRange, Y_RANGE);
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

	private static List<DoublingSummary> loadDoublingSummaries(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		Map<DoublingGroup, List<Double>> groups = new LinkedHashMap<>();
		for (int row=1; row<csv.getNumRows(); row++) {
			SamplingMethod method = SamplingMethod.valueOf(csv.get(row, 1));
			int lowerCount = Integer.parseInt(csv.get(row, 4));
			int upperCount = Integer.parseInt(csv.get(row, 5));
			ConvergenceMetric metric = ConvergenceMetric.fromLabel(csv.get(row, 6));
			addValue(groups, new DoublingGroup(method, lowerCount, upperCount, metric,
					ConvergenceSummary.MEAN_ABSOLUTE), Double.parseDouble(csv.get(row, 8)));
			addValue(groups, new DoublingGroup(method, lowerCount, upperCount, metric,
					ConvergenceSummary.MAXIMUM_ABSOLUTE), Double.parseDouble(csv.get(row, 10)));
		}
		List<DoublingSummary> rows = new ArrayList<>();
		for (Map.Entry<DoublingGroup, List<Double>> entry : groups.entrySet()) {
			DoublingGroup group = entry.getKey();
			double[] values = entry.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			rows.add(new DoublingSummary(group.method(), group.lowerCount(), group.upperCount(), group.metric(),
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
		return summary == ConvergenceSummary.MEAN_ABSOLUTE ? "mean_abs" : "max_abs";
	}

	private interface SummaryRow {
		int count();
		ConvergenceMetric metric();
		double minimum();
		double median();
		double maximum();
	}

	private record ReferenceSummary(SamplingMethod method, int sampleCount, String reference, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary,
			int realizations, double minimum, double median, double maximum) implements SummaryRow {
		@Override public int count() { return sampleCount; }
	}

	private record DoublingSummary(SamplingMethod method, int lowerCount, int upperCount, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary,
			int realizations, double minimum, double median, double maximum) implements SummaryRow {
		@Override public int count() { return upperCount; }
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

	private record DoublingGroup(SamplingMethod method, int lowerCount, int upperCount,
			ConvergenceMetric metric,
			ConvergenceSummary spatialSummary) {}

	private record RealizationPairGroup(SamplingMethod method, int sampleCount,
			ConvergenceMetric metric, ConvergenceSummary spatialSummary) {}
}
