package scratch.kevin.sampling;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.data.Range;
import org.opensha.commons.util.ColorUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
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

	private static final Range Y_RANGE = new Range(5e-3, 2e1);
	private static final Range SIGNED_Y_RANGE = new Range(-1d, 1d);
	private static final double[] CONVERGENCE_REFERENCE_ERRORS = {0.1d, 1d, 10d};
//	private static final double[] CONVERGENCE_REFERENCE_ERRORS = {0.01, 0.1d, 1d, 10d, 100};
	private static final String SOBOL_REFERENCE = HazardConvergenceCalcs.LOO_SOBOL_REFERENCE_NAME;
	private static final String POOLED_SOBOL_REFERENCE = HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME;
	private static final String MCS_REFERENCE = HazardConvergenceCalcs.MCS_REFERENCE_NAME;
	private static final SamplingMethod MCS = SamplingMethod.MONTE_CARLO;
	private static final SamplingMethod SOBOL = SamplingMethod.OWEN_SCRAMBLED_SOBOL;

	private static final boolean INCLUDE_SIGNED_VARIABILITY_UNCERTAINTIES = true;

//	private static final ConvergenceMetric[] PLOT_METRICS = ConvergenceMetric.values();
	static final ConvergenceMetric[] PLOT_METRICS = {
			ConvergenceMetric.MEAN_HAZARD,
			ConvergenceMetric.STANDARD_DEVIATION,
			ConvergenceMetric.CENTRAL_68_RANGE,
			ConvergenceMetric.CENTRAL_95_RANGE
	};

	private static final Map<ConvergenceMetric, Color> METRIC_COLORS = Map.of(
			ConvergenceMetric.MEAN_HAZARD, Colors.tab_blue,
			ConvergenceMetric.STANDARD_DEVIATION, Colors.tab_orange,
			ConvergenceMetric.CENTRAL_68_RANGE, Colors.tab_green,
			ConvergenceMetric.CENTRAL_95_RANGE, Colors.tab_red,
			ConvergenceMetric.IQR, Colors.tab_purple);

	private static final Map<ConvergenceMetric, PlotSymbol> METRIC_SYMBOLS = Map.of(
			ConvergenceMetric.MEAN_HAZARD, PlotSymbol.FILLED_CIRCLE,
			ConvergenceMetric.STANDARD_DEVIATION, PlotSymbol.FILLED_INV_TRIANGLE,
			ConvergenceMetric.CENTRAL_68_RANGE, PlotSymbol.FILLED_SQUARE,
			ConvergenceMetric.CENTRAL_95_RANGE, PlotSymbol.FILLED_DIAMOND,
			ConvergenceMetric.IQR, PlotSymbol.FILLED_TRIANGLE);

	static String getMethodName(SamplingMethod method) {
		if (method == SamplingMethod.OWEN_SCRAMBLED_SOBOL)
			return SamplingMethod.SOBOL.getShortName();
		return method.getShortName();
	}

	private static final Map<SamplingMethod, String> METHOD_FILE_PREFIXES;
	static {
		Map<SamplingMethod, String> prefixes = new HashMap<>();
		for (SamplingMethod method : SamplingMethod.values())
			prefixes.put(method, getMethodName(method).toLowerCase().replaceAll("-", "_").replace("'", ""));
		METHOD_FILE_PREFIXES = prefixes;
	}

	private static final boolean PLOT_INDV_MEANS = true;

	public static void main(String[] args) throws IOException {
		File convergenceDir = new File(PaperPaths.FIGURES_DIR, "hazard_convergence");
		plotPeriod(new File(convergenceDir, "pga_two_in_50"), "PGA");
		plotPeriod(new File(convergenceDir, "1s_sa_two_in_50"), "1s SA");
	}

	static void plotPeriod(File outputDir, String periodName) throws IOException {
		File referenceFile = new File(outputDir, "reference_comparisons.csv");
		List<ReferenceSummary> references = loadReferenceSummaries(referenceFile);
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
		for (boolean sobolPool : new boolean[] {true, false}) {
			plotMethodReference(outputDir, references, sobolPool, ConvergenceSummary.MEAN_ABSOLUTE,
					"Absolute difference (%)");
			plotMethodReference(outputDir, references, sobolPool, ConvergenceSummary.MEAN_SIGNED,
					"Signed bias (%)");
		}
		plotRealizationPairs(outputDir, realizationPairs, ConvergenceSummary.MEAN_ABSOLUTE,
				"Absolute difference (%)");
		plotMeanHazardDoublingImprovement(outputDir, loadReferenceValues(referenceFile));
	}

	/**
	 * Plots the reduction in spatially averaged mean-hazard error when the sample count is doubled. Sobol' values use
	 * nested prefixes from the same scramble. MCS values pair each disjoint span with the twice-as-long span having the
	 * same starting index. This keeps realization-to-realization variation out of each improvement ratio.
	 */
	private static void plotMeanHazardDoublingImprovement(File outputDir, List<ReferenceValue> rows)
			throws IOException {
		int maxSobolRunSize = rows.stream()
				.filter(row -> row.method() == SOBOL && row.reference().equals(MCS_REFERENCE))
				.mapToInt(ReferenceValue::maximumRunSize).max().orElse(0);
		List<PairedImprovement> improvements = new ArrayList<>();

		// Use only the longest Sobol' runs so every doubling is drawn from a consistent collection of scrambles.
		Map<String, Map<Integer, ReferenceValue>> sobolByRun = new LinkedHashMap<>();
		for (ReferenceValue row : rows) {
			if (row.method() == SOBOL && row.maximumRunSize() == maxSobolRunSize
					&& row.reference().equals(MCS_REFERENCE) && row.metric() == ConvergenceMetric.MEAN_HAZARD)
				sobolByRun.computeIfAbsent(row.run(), unused -> new HashMap<>()).put(row.sampleCount(), row);
		}
		Map<Integer, List<Double>> sobolRatios = new LinkedHashMap<>();
		for (Map<Integer, ReferenceValue> runRows : sobolByRun.values()) {
			for (ReferenceValue lower : runRows.values()) {
				ReferenceValue upper = runRows.get(2*lower.sampleCount());
				if (upper != null)
					addImprovement(sobolRatios, lower, upper);
			}
		}
		addImprovementSummaries(improvements, SOBOL, sobolRatios);

		// Primary MCS spans use a fixed-size leave-span-out reference with an independent reserve replacement.
		Map<Long, ReferenceValue> mcsByStartAndCount = new HashMap<>();
		for (ReferenceValue row : rows) {
			if (row.method() == MCS && row.reference().equals(HazardConvergenceCalcs.REPLACED_MCS_REFERENCE_NAME)
					&& row.metric() == ConvergenceMetric.MEAN_HAZARD)
				mcsByStartAndCount.put(spanKey(row.startIndex(), row.sampleCount()), row);
		}
		Map<Integer, List<Double>> mcsRatios = new LinkedHashMap<>();
		for (ReferenceValue lower : mcsByStartAndCount.values()) {
			// Match the first N-point half to its containing 2N-point span. Restricting to aligned starts avoids
			// counting two strongly related child spans against the same parent as independent observations.
			if (lower.startIndex() % (2*lower.sampleCount()) != 0)
				continue;
			ReferenceValue upper = mcsByStartAndCount.get(spanKey(lower.startIndex(), 2*lower.sampleCount()));
			if (upper != null)
				addImprovement(mcsRatios, lower, upper);
		}
		addImprovementSummaries(improvements, MCS, mcsRatios);

		int[] lowerCounts = improvements.stream().mapToInt(PairedImprovement::lowerCount).distinct().sorted().toArray();
		if (lowerCounts.length == 0)
			return;
		System.out.println("Building plot: paired mean-hazard doubling improvement");
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		ArbitrarilyDiscretizedFunc noImprovement = horizontalLine(lowerCounts.length, 1d);
		funcs.add(noImprovement);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.75f, Color.GRAY));
		addRateReference(funcs, chars, lowerCounts.length, Math.sqrt(2d), "N^-1/2 (1.41x)", PlotLineType.DASHED);
		addRateReference(funcs, chars, lowerCounts.length, 2d, "N^-1 (2x)", PlotLineType.DOTTED);
		addRateReference(funcs, chars, lowerCounts.length, Math.pow(2d, 1.5d), "N^-3/2 (2.83x)",
				PlotLineType.DOTTED_AND_DASHED);

		for (SamplingMethod method : new SamplingMethod[] { MCS, SOBOL }) {
			Color color = method == MCS ? Colors.tab_red : Colors.tab_blue;
			PlotSymbol symbol = method == MCS ? PlotSymbol.FILLED_CIRCLE : PlotSymbol.FILLED_SQUARE;
			ArbitrarilyDiscretizedFunc median = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lower = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upper = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<lowerCounts.length; i++) {
				int count = lowerCounts[i];
				PairedImprovement summary = improvements.stream()
						.filter(candidate -> candidate.method() == method && candidate.lowerCount() == count)
						.findFirst().orElse(null);
				if (summary == null)
					continue;
				median.set((double)i, summary.median());
				double factor = Double.isFinite(summary.logStandardDeviation())
						? Math.exp(summary.logStandardDeviation()) : 1d;
				lower.set((double)i, summary.median()/factor);
				upper.set((double)i, summary.median()*factor);
				System.out.println("\t"+getMethodName(method)+" "+count+" -> "+(2*count)+": median="
						+(float)summary.median()+", pairs="+summary.pairs());
			}
			if (median.size() == 0)
				continue;
			UncertainArbDiscFunc uncertainty = new UncertainArbDiscFunc(median, lower, upper);
			funcs.add(uncertainty);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
					ColorUtils.transparent(color, 55)));
			median.setName(getMethodName(method));
			funcs.add(median);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, symbol, 5f, color));
		}

		PlotSpec plot = new PlotSpec(funcs, chars, "Mean hazard convergence rate", "Initial sample count, N",
				"Error reduction, E(N) / E(2N)");
		plot.setLegendInset(RectangleAnchor.BOTTOM, 0.5, 0.025, 0.9, false);
		HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
		gp.getPlotPrefs().setPlotLabelFontSize(10);
		gp.getPlotPrefs().setLegendFontSize(8);
		gp.getPlotPrefs().setLegendLineLength(6d);
		gp.getPlotPrefs().setTickLabelFontSize(8);
		Range xRange = lowerCounts.length == 1 ? new Range(-0.5, 0.5) : new Range(-0.2d, lowerCounts.length-0.8d);
		gp.drawGraphPanel(plot, false, false, xRange, new Range(0.5d, 4d));
		PlotUtils.setXTick(gp, 1d);
		((NumberAxis)gp.getXAxis()).setNumberFormatOverride(categoryFormat(
				Arrays.stream(lowerCounts).mapToObj(Integer::toString).toArray(String[]::new)));
		PlotUtils.writePrintPlots(outputDir, "mean_hazard_doubling_improvement", gp,
				PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d, 3.4d, 300, true, true, false);
	}

	private static void addImprovement(Map<Integer, List<Double>> ratios,
			ReferenceValue lower, ReferenceValue upper) {
		if (lower.meanAbsolutePercentChange() > 0d && upper.meanAbsolutePercentChange() > 0d)
			ratios.computeIfAbsent(lower.sampleCount(), unused -> new ArrayList<>())
					.add(lower.meanAbsolutePercentChange()/upper.meanAbsolutePercentChange());
	}

	private static void addImprovementSummaries(List<PairedImprovement> summaries,
			SamplingMethod method, Map<Integer, List<Double>> ratios) {
		for (Map.Entry<Integer, List<Double>> entry : ratios.entrySet()) {
			double[] values = entry.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			summaries.add(new PairedImprovement(method, entry.getKey(), values.length,
					StatUtils.percentile(values, 50d), HazardConvergenceCalcs.logStandardDeviation(values)));
		}
	}

	private static long spanKey(int startIndex, int sampleCount) {
		return ((long)startIndex << 32) ^ (sampleCount & 0xffffffffL);
	}

	private static ArbitrarilyDiscretizedFunc horizontalLine(int count, double y) {
		ArbitrarilyDiscretizedFunc line = new ArbitrarilyDiscretizedFunc();
		line.set(-0.2d, y);
		line.set(count-0.8d, y);
		return line;
	}

	private static void addRateReference(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars,
			int count, double factor, String name, PlotLineType lineType) {
		ArbitrarilyDiscretizedFunc reference = horizontalLine(count, factor);
		reference.setName(name);
		funcs.add(reference);
		chars.add(new PlotCurveCharacterstics(lineType, 1f, Color.DARK_GRAY));
	}

	private static String referenceFor(SamplingMethod method, boolean sobolPool) {
		return sobolPool ? (method == SOBOL ? SOBOL_REFERENCE : POOLED_SOBOL_REFERENCE)
				: (method == MCS ? HazardConvergenceCalcs.REPLACED_MCS_REFERENCE_NAME : MCS_REFERENCE);
	}

	private static boolean matchesReference(ReferenceSummary row, SamplingMethod method, boolean sobolPool) {
		if (sobolPool && method == SOBOL)
			// A fixed-size consensus only leaves out Sobol runs of that size; other sizes use the full pool.
			return row.reference().equals(SOBOL_REFERENCE) || row.reference().equals(POOLED_SOBOL_REFERENCE);
		if (!sobolPool && method == MCS)
			// Primary-pool spans use reserve replacement; the independent reserve prefix uses the full pool.
			return row.reference().equals(HazardConvergenceCalcs.REPLACED_MCS_REFERENCE_NAME)
					|| row.reference().equals(MCS_REFERENCE);
		return row.reference().equals(referenceFor(method, sobolPool));
	}

	private static boolean includeSummary(ConvergenceSummary actual, ConvergenceSummary requested) {
		return actual == requested || requested == ConvergenceSummary.MEAN_ABSOLUTE
				&& actual == ConvergenceSummary.MAXIMUM_ABSOLUTE;
	}

	private static void plotReference(File outputDir, List<ReferenceSummary> rows,
			SamplingMethod method, boolean sobolPool, ConvergenceSummary summary) throws IOException {
		List<ReferenceSummary> matching = rows.stream().filter(row -> row.method() == method
				&& matchesReference(row, method, sobolPool)
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
		writePlot(outputDir, prefix, getMethodName(method)+" vs "+(sobolPool ? "Sobol' pool" : "MCS pool"),
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
			int index = labels.size();
			for (ReferenceSummary row : rows) {
				if (row.method() == method && row.sampleCount() == sampleCount
						&& matchesReference(row, method, sobolConsensus)
						&& includeSummary(row.spatialSummary(), spatialSummary)) {
					matching.add(new MethodSummary(index, row.metric(), row.spatialSummary(), row.realizations(),
							row.mean(), row.standardDeviation(), row.minimum(), row.median(), row.maximum(),
							row.logStandardDeviation(), row.individualValues()));
				}
			}
			if (matching.stream().anyMatch(row -> row.count() == index))
				labels.add(getMethodName(method));
		}
		if (matching.isEmpty())
			return;
		String referencePrefix = sobolConsensus ? "pooled_sobol" : "pooled_mcs";
//		String title = sampleCount+" samples versus "
//				+(sobolConsensus ? POOLED_SOBOL_REFERENCE : MCS_REFERENCE);
		String title = sampleCount+" samples vs "
				+(sobolConsensus ? "Sobol' pool" : "MCS pool");
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
							row.mean(), row.standardDeviation(), row.minimum(), row.median(), row.maximum(),
							row.logStandardDeviation(), row.individualValues()));
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
		System.out.println("Building plot: "+title);
		boolean signed = primary == ConvergenceSummary.MEAN_SIGNED;
		Map<ConvergenceMetric, List<? extends SummaryRow>> byMetric = new LinkedHashMap<>();
		for (ConvergenceMetric metric : PLOT_METRICS) {
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
		if (primary == ConvergenceSummary.MEAN_ABSOLUTE && xLabel.equals("Sample count")) {
			int referenceCount = counts[0];
			for (int r=0; r<CONVERGENCE_REFERENCE_ERRORS.length; r++) {
				double refError = CONVERGENCE_REFERENCE_ERRORS[r];
				ArbitrarilyDiscretizedFunc reference = new ArbitrarilyDiscretizedFunc();
				for (int i=-1; i<=counts.length; i++) {
					double count;
					if (i == -1)
						count = counts[0]/2;
					else if (i == counts.length)
						count = counts[i-1]*2;
					else
						count = counts[i];
					reference.set((double)i, refError*Math.sqrt((double)referenceCount/count));
				}
				if (r == 0)
					reference.setName("N^-½");
//					reference.setName("N^⁻½");
//					reference.setName("ₙ⁻½");
//					reference.setName("1/√N");
				funcs.add(reference);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.DARK_GRAY));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.GRAY));
			}
		}
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
			DefaultXY_DataSet indvMeans = entry.getKey() == ConvergenceMetric.MEAN_HAZARD && PLOT_INDV_MEANS ? new DefaultXY_DataSet() : null;
			for (int i=0; i<counts.length; i++) {
				int count = counts[i];
				// More than one accepted reference label can represent the same plotted point. In particular,
				// native Sobol runs use the full fixed-size Sobol pool while prefixes of runs in that pool use
				// leave-one-out references. Combine their realization values before calculating plot statistics.
				SummaryRow row = combineRows(entry.getValue().stream()
						.filter(candidate -> candidate.count() == count).toList());
				if (indvMeans != null) {
					double[] values = row.individualValues();
					if (values.length > 3000) {
						double[] subValues = new double[2000];
						Random r = new Random(values.length);
						for (int j=0; j<subValues.length; j++)
							subValues[j] = values[r.nextInt(values.length)];
						System.out.println("\tDownsampling "+values.length+" down to "+subValues.length);
						values = subValues;
					}
					for (double value : values)
						indvMeans.set((double)i, value);

					System.out.println("\t"+values.length+" mean values for "+title+", count="+count);
				}
				double center = signed ? row.mean() : row.median();
				median.set((double)i, center);
				if (signed) {
					double standardDeviation = row.standardDeviation();
					if (!Double.isFinite(standardDeviation))
						standardDeviation = 0d;
					lower.set((double)i, center-standardDeviation);
					upper.set((double)i, center+standardDeviation);
				} else {
					double logSD = row.logStandardDeviation();
					double factor = Double.isFinite(logSD) ? Math.exp(logSD) : 1d;
					lower.set((double)i, row.median()/factor);
					upper.set((double)i, row.median()*factor);
				}
			}
			Color color = METRIC_COLORS.get(entry.getKey());
			PlotSymbol sym = METRIC_SYMBOLS.get(entry.getKey());
			PlotSymbol outlineSym = PlotSymbol.getOutlineSymbol(sym);
			// Positive quantities use a shaded, median-centered multiplicative range defined by one standard
			// deviation of the log-transformed values. For signed biases, retain that shading for mean hazard but
			// draw the other mean +/- standard-deviation ranges as unobtrusive dotted bounds to avoid overlapping
			// shaded regions.
			if (!signed || entry.getKey() == ConvergenceMetric.MEAN_HAZARD) {
				UncertainArbDiscFunc uncertainty = new UncertainArbDiscFunc(median, lower, upper);
				funcs.add(uncertainty);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
						ColorUtils.transparent(color, 70)));
			} else if (INCLUDE_SIGNED_VARIABILITY_UNCERTAINTIES) {
				funcs.add(lower);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f,
						ColorUtils.transparent(color, 120)));
				funcs.add(upper);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f,
						ColorUtils.transparent(color, 120)));
			}
			if (indvMeans != null) {
				medianFuncs.add(indvMeans);
				medianChars.add(new PlotCurveCharacterstics(sym, 1.5f, color.darker().darker()));
			}
			if (outlineSym != null) {
				// add slightly darker outline overlay
				ArbitrarilyDiscretizedFunc clone = median.deepClone();
				clone.setName(null);
				medianFuncs.add(clone);
				medianChars.add(new PlotCurveCharacterstics(outlineSym, 5f, color.darker().darker()));
			}
			median.setName(entry.getKey().shortLabel);
			medianFuncs.add(median);
			if (primary == ConvergenceSummary.MEAN_ABSOLUTE) {
				ArbitrarilyDiscretizedFunc typicalWorst = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<counts.length; i++) {
					int count = counts[i];
					List<? extends SummaryRow> matchingMaximums = rows.stream()
							.filter(row -> row.metric() == entry.getKey() && row.count() == count
									&& row.spatialSummary() == ConvergenceSummary.MAXIMUM_ABSOLUTE)
							.toList();
					// Each value is the largest spatial error in one realization. Their median represents the
					// typical realization's worst site without increasing merely because more runs are available.
					if (!matchingMaximums.isEmpty())
						typicalWorst.set((double)i, combineRows(matchingMaximums).median());
				}
				if (typicalWorst.size() > 0) {
					maxFuncs.add(typicalWorst); // Unnamed, so it adds no legend entry.
					maxChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 0.5f,
							PlotSymbol.getOutlineSymbol(sym), 3f, color));
				}
			}
			medianChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sym, 5f, color));
		}
		// Put every envelope in the dataset first so all median lines render above all shading.
		funcs.addAll(maxFuncs);
		chars.addAll(maxChars);
		// add copies for names (then clear the names)
		for (int i=0; i<medianFuncs.size(); i++) {
			XY_DataSet func = medianFuncs.get(i);
			if (func.getName() != null && !func.getName().isBlank()) {
				XY_DataSet legendFunc = new DefaultXY_DataSet(-1000, -1000);
				legendFunc.setName(func.getName());
				funcs.add(legendFunc);
				chars.add(medianChars.get(i));
				func.setName(null);
			}
		}
		// reverse them so that mean is on top
		Collections.reverse(medianFuncs);
		Collections.reverse(medianChars);
		funcs.addAll(medianFuncs);
		chars.addAll(medianChars);

		PlotSpec plot = new PlotSpec(funcs, chars, title, xLabel, yLabel);
//		plot.setLegendVisible(true);
//		plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
//		plot.setLegendInset(RectangleAnchor.TOP);
//		plot.setLegendInset(RectangleAnchor.TOP, 0.5, 0.975, 0.9, false);
		plot.setLegendInset(RectangleAnchor.BOTTOM, 0.5, 0.025, 0.925, false);
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
					group.spatialSummary(), values.length, StatUtils.mean(values),
					HazardConvergenceCalcs.standardDeviation(values), StatUtils.min(values),
					StatUtils.percentile(values, 50d), StatUtils.max(values),
					HazardConvergenceCalcs.logStandardDeviation(values), values));
		}
		return rows;
	}

	private static List<ReferenceValue> loadReferenceValues(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		List<ReferenceValue> rows = new ArrayList<>(csv.getNumRows()-1);
		for (int row=1; row<csv.getNumRows(); row++) {
			rows.add(new ReferenceValue(csv.get(row, 0), SamplingMethod.valueOf(csv.get(row, 1)),
					Integer.parseInt(csv.get(row, 3)), Integer.parseInt(csv.get(row, 4)), csv.get(row, 5),
					ConvergenceMetric.fromLabel(csv.get(row, 7)), Double.parseDouble(csv.get(row, 9)),
					Integer.parseInt(csv.get(row, 16)), Integer.parseInt(csv.get(row, 17))));
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
					group.spatialSummary(), values.length, StatUtils.mean(values),
					HazardConvergenceCalcs.standardDeviation(values), StatUtils.min(values),
					StatUtils.percentile(values, 50d), StatUtils.max(values),
					HazardConvergenceCalcs.logStandardDeviation(values), values));
		}
		return rows;
	}

	private static <K> void addValue(Map<K, List<Double>> groups, K key, double value) {
		groups.computeIfAbsent(key, unused -> new ArrayList<>()).add(value);
	}

	static SummaryRow combineRows(List<? extends SummaryRow> rows) {
		if (rows.isEmpty())
			throw new IllegalArgumentException("Cannot combine an empty set of summary rows");
		SummaryRow first = rows.get(0);
		int size = rows.stream().mapToInt(row -> row.individualValues().length).sum();
		double[] values = new double[size];
		int offset = 0;
		for (SummaryRow row : rows) {
			if (row.count() != first.count() || row.metric() != first.metric()
					|| row.spatialSummary() != first.spatialSummary())
				throw new IllegalArgumentException("Cannot combine summary rows for different plotted quantities");
			System.arraycopy(row.individualValues(), 0, values, offset, row.individualValues().length);
			offset += row.individualValues().length;
		}
		return new CombinedSummary(first.count(), first.metric(), first.spatialSummary(),
				StatUtils.mean(values), HazardConvergenceCalcs.standardDeviation(values), StatUtils.min(values),
				StatUtils.percentile(values, 50d), StatUtils.max(values),
				HazardConvergenceCalcs.logStandardDeviation(values), values);
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
		double mean();
		double standardDeviation();
		double minimum();
		double median();
		double maximum();
		double logStandardDeviation();
		double[] individualValues();
	}

	record ReferenceSummary(SamplingMethod method, int sampleCount, String reference, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary,
			int realizations, double mean, double standardDeviation, double minimum, double median, double maximum,
			double logStandardDeviation, double[] individualValues) implements SummaryRow {
		@Override public int count() { return sampleCount; }
	}

	private record ReferenceValue(String run, SamplingMethod method, int maximumRunSize, int sampleCount,
			String reference, ConvergenceMetric metric, double meanAbsolutePercentChange,
			int startIndex, int endIndex) {}

	private record PairedImprovement(SamplingMethod method, int lowerCount, int pairs,
			double median, double logStandardDeviation) {}

	private record MethodSummary(int methodIndex, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary, int realizations,
			double mean, double standardDeviation, double minimum, double median, double maximum,
			double logStandardDeviation, double[] individualValues) implements SummaryRow {
		@Override public int count() { return methodIndex; }
	}

	private record CombinedSummary(int count, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary, double mean, double standardDeviation,
			double minimum, double median, double maximum, double logStandardDeviation,
			double[] individualValues) implements SummaryRow {}


	private record RealizationPairSummary(SamplingMethod method, int sampleCount, ConvergenceMetric metric,
			ConvergenceSummary spatialSummary, int realizationPairs,
			double mean, double standardDeviation, double minimum, double median, double maximum,
			double logStandardDeviation, double[] individualValues) {}

	private record ReferenceGroup(SamplingMethod method, int sampleCount, String reference,
			ConvergenceMetric metric,
			ConvergenceSummary spatialSummary) {}

	private record RealizationPairGroup(SamplingMethod method, int sampleCount,
			ConvergenceMetric metric, ConvergenceSummary spatialSummary) {}
}
