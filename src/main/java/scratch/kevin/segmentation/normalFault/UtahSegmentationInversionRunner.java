package scratch.kevin.segmentation.normalFault;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.time.LocalDate;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Properties;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SegmentationProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationRuptureExemptions;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;

import com.google.common.base.Preconditions;
import com.google.gson.GsonBuilder;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

/**
 * Utah comparison based on scratch.kevin.nshm23.HardcodedInversionFactoryRunner. No arguments runs MID and both
 * geometric fits sequentially, using the standard factory, cache, and report output parent. Other branch values,
 * iteration budget, and solver settings match across cases. Completed solutions are reused on subsequent launches.
 *
 * <p>Options: --prepare-only, --index-only, --case MID|MANUAL_60K|AUTOMATED_5KM|ALL, --threads N,
 * --iterations-per-rupture N, --output-dir DIR, --no-reports.
 * Preparation writes regional configurations for inspection; cluster-specific solving builds its own subset configs.
 */
public class UtahSegmentationInversionRunner {

	private static final File OUTPUT_PARENT = new File("/home/kevin/markdown/inversions");
	private static final File CACHE_DIR = new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache");

	private enum Calculation {
		MID(null),
		MANUAL_60K(NSHM23_NormalFaultSegmentationModels.MANUAL_60K),
		AUTOMATED_5KM(NSHM23_NormalFaultSegmentationModels.AUTOMATED_5KM);

		final NSHM23_NormalFaultSegmentationModels geometry;
		Calculation(NSHM23_NormalFaultSegmentationModels geometry) { this.geometry = geometry; }
		String directoryName() { return name().toLowerCase(Locale.ROOT); }
		String displayName() { return geometry == null ? NSHM23_SegmentationModels.MID.getName() : geometry.getName(); }
	}

	public static void main(String[] args) throws IOException {
		int threads = 16;
		long iterations = NSHM23_InvConfigFactory.NUM_ITERS_PER_RUP_DEFAULT;
		boolean prepareOnly = false;
		boolean indexOnly = false;
		boolean reports = true;
		List<Calculation> calculations = Arrays.asList(Calculation.values());
		String date = LocalDate.now(ZoneId.of("America/Los_Angeles")).format(DateTimeFormatter.ofPattern("yyyy_MM_dd"));
		File output = new File(OUTPUT_PARENT, date+"-nshm23-utah-normal_fault_segmentation-unconstrained_rup_exemption");
		for (int i=0; i<args.length; i++) {
			switch (args[i]) {
			case "--prepare-only": prepareOnly = true; break;
			case "--index-only": indexOnly = true; break;
			case "--no-reports": reports = false; break;
			case "--threads": threads = Integer.parseInt(value(args, ++i)); break;
			case "--iterations-per-rupture": iterations = Long.parseLong(value(args, ++i)); break;
			case "--output-dir": output = new File(value(args, ++i)); break;
			case "--case":
				String selected = value(args, ++i).toUpperCase(Locale.ROOT);
				calculations = selected.equals("ALL") ? Arrays.asList(Calculation.values())
						: List.of(Calculation.valueOf(selected));
				break;
			default: throw new IllegalArgumentException("Unknown option: "+args[i]);
			}
		}
		Preconditions.checkArgument(threads > 0 && iterations > 0, "Threads and iterations must be positive");
		Preconditions.checkState(output.isDirectory() || output.mkdirs(), "Cannot create %s", output);
		checkSettings(output, threads, iterations);
		System.out.println("Utah segmentation comparison: "+output.getAbsolutePath());
		writeIndex(output, threads, iterations);
		if (indexOnly) return;

		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(CACHE_DIR);
		factory.setNumItersPerRup(iterations);
		FaultSystemRupSet reference = null;
		for (Calculation calculation : calculations) {
			File directory = new File(output, calculation.directoryName());
			Preconditions.checkState(directory.isDirectory() || directory.mkdir(), "Cannot create %s", directory);
			File solutionFile = new File(directory, "solution.zip");
			if (solutionFile.exists()) {
				System.out.println("Reusing completed solution: "+solutionFile);
				FaultSystemRupSet completed = FaultSystemRupSet.load(new File(directory, "rupSet.zip"));
				if (reference == null) reference = completed;
				else checkSameRuptures(reference, completed);
				continue;
			}
			LogicTreeBranch<LogicTreeNode> branch = buildBranch(calculation);
			System.out.println("Preparing "+calculation+": "+branch);
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, threads);
			if (reference == null) reference = rupSet;
			else checkSameRuptures(reference, rupSet);
			InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, threads);
			Preconditions.checkNotNull(config, "No inversion configuration for %s", calculation);
			InversionConfiguration.writeJSON(config, new File(directory, "config.json"));
			writeExemptionAudit(directory, rupSet, config);
			Files.writeString(new File(directory, "branch.json").toPath(),
					new GsonBuilder().setPrettyPrinting().create().toJson(branch, LogicTreeBranch.class),
					StandardCharsets.UTF_8);
			rupSet.write(new File(directory, "rupSet.zip"));
			if (prepareOnly) {
				validateInputs(directory, rupSet, config);
				System.out.println("Prepared "+calculation+": "+rupSet.getNumSections()+" subsections, "
						+rupSet.getNumRuptures()+" ruptures; completion="+config.getCompletionCriteria());
				continue;
			}
			System.out.println("Running "+calculation);
			FaultSystemSolution solution = Inversions.run(rupSet, factory, branch, threads, null);
			Preconditions.checkNotNull(solution, "No solution for %s", calculation);
			solution.write(solutionFile);
			writeIndex(output, threads, iterations);
		}
		if (!prepareOnly && reports) writeReports(output, calculations);
		writeIndex(output, threads, iterations);
	}

	private static void writeExemptionAudit(File directory, FaultSystemRupSet rupSet, InversionConfiguration config)
			throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("from_id", "to_id", "from_name", "to_name", "minimum_subsection_count", "exempt_rupture_count");
		int sides = 0;
		for (InversionConstraint constraint : config.getConstraints()) {
			if (!(constraint instanceof SegmentationProbabilityConstraint)) continue;
			SegmentationRuptureExemptions exemptions =
					((SegmentationProbabilityConstraint)constraint).getRuptureExemptions();
			for (IDPairing side : new java.util.TreeSet<>(exemptions.getExemptSides())) {
				BitSet exempt = exemptions.getExemptRuptures(side.getID1(), side.getID2());
				int minimum = rupSet.getSectionsIndicesForRup(exempt.nextSetBit(0)).size();
				csv.addLine(""+side.getID1(), ""+side.getID2(),
						rupSet.getFaultSectionData(side.getID1()).getSectionName(),
						rupSet.getFaultSectionData(side.getID2()).getSectionName(), ""+minimum, ""+exempt.cardinality());
				sides++;
			}
		}
		csv.writeToFile(new File(directory, "minimum_rupture_exemptions.csv"));
		System.out.println("Minimum-rupture exemption applies to "+sides+" internal transition sides");
	}

	private static void validateInputs(File directory, FaultSystemRupSet rupSet, InversionConfiguration config)
			throws IOException {
		InversionInputGenerator inputs = new InversionInputGenerator(rupSet, config);
		inputs.generateInputs();
		List<String> summary = new ArrayList<>();
		DoubleMatrix2D[] matrices = {inputs.getA(), inputs.getA_ineq()};
		double[][] targets = {inputs.getD(), inputs.getD_ineq()};
		for (int i=0; i<matrices.length; i++) {
			DoubleMatrix2D matrix = matrices[i];
			if (matrix == null) continue;
			matrix.forEachNonZero((row, column, coefficient) -> {
				Preconditions.checkState(Double.isFinite(coefficient), "Nonfinite matrix coefficient");
				return coefficient;
			});
			for (double target : targets[i]) Preconditions.checkState(Double.isFinite(target), "Nonfinite target");
			summary.add((i == 0 ? "Equality" : "Inequality")+": "+matrix.rows()+" rows, "+matrix.columns()
					+" columns, "+matrix.cardinality()+" nonzero entries; all coefficients and targets finite.");
		}
		Files.write(new File(directory, "input_summary.txt").toPath(), summary, StandardCharsets.UTF_8);
	}

	private static String value(String[] args, int index) {
		Preconditions.checkArgument(index < args.length, "Missing option value");
		return args[index];
	}

	private static LogicTreeBranch<LogicTreeNode> buildBranch(Calculation calculation) {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>(NSHM23_LogicTreeBranch.levelsOnFault);
		if (calculation.geometry != null)
			levels.set(levels.indexOf(NSHM23_LogicTreeBranch.SEG), LogicTreeLevel.forEnum(
					NSHM23_NormalFaultSegmentationModels.class, "Normal Fault Segmentation", "NormalSeg"));
		levels.add(NSHM23_LogicTreeBranch.SINGLE_STATES);
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
			branch.setValue(node instanceof NSHM23_SegmentationModels && calculation.geometry != null
					? calculation.geometry : node);
		branch.setValue(NSHM23_SingleStates.UT);
		return branch;
	}

	private static void checkSameRuptures(FaultSystemRupSet reference, FaultSystemRupSet current) {
		Preconditions.checkState(reference.areSectionsEquivalentTo(current.getFaultSectionDataList()),
				"Comparison cases must use identical fault sections");
		Preconditions.checkState(reference.getNumRuptures() == current.getNumRuptures(), "Rupture counts differ");
		for (int r=0; r<reference.getNumRuptures(); r++)
			Preconditions.checkState(reference.getSectionsIndicesForRup(r).equals(current.getSectionsIndicesForRup(r))
					&& reference.getMagForRup(r) == current.getMagForRup(r), "Rupture %s differs between cases", r);
	}

	private static void checkSettings(File output, int threads, long iterations) throws IOException {
		Properties settings = new Properties();
		settings.setProperty("threads", ""+threads);
		settings.setProperty("iterations_per_rupture", ""+iterations);
		settings.setProperty("cluster_specific", ""+NSHM23_InvConfigFactory.SOLVE_CLUSTERS_INDIVIDUALLY_DEFAULT);
		settings.setProperty("reference_branch", NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.buildFileName());
		settings.setProperty("state", NSHM23_SingleStates.UT.name());
		settings.setProperty("geometric_exemption_policy",
				NSHM23_NormalFaultSegmentationModels.AUTOMATED_5KM.getRuptureExemptionPolicy().name());
		File file = new File(output, "run-settings.properties");
		if (file.exists()) {
			Properties previous = new Properties();
			try (Reader reader = Files.newBufferedReader(file.toPath(), StandardCharsets.UTF_8)) {
				previous.load(reader);
			}
			Preconditions.checkState(settings.equals(previous), "Settings changed; select a new --output-dir");
		} else {
			try (Writer writer = Files.newBufferedWriter(file.toPath(), StandardCharsets.UTF_8)) {
				settings.store(writer, "Common Utah inversion settings");
			}
		}
	}

	private static void writeReports(File output, List<Calculation> calculations) throws IOException {
		File baselineFile = new File(new File(output, Calculation.MID.directoryName()), "solution.zip");
		FaultSystemSolution baseline = baselineFile.exists() ? FaultSystemSolution.load(baselineFile) : null;
		for (Calculation calculation : calculations) {
			File directory = new File(output, calculation.directoryName());
			File solutionFile = new File(directory, "solution.zip");
			if (!solutionFile.exists()) continue;
			FaultSystemSolution solution = calculation == Calculation.MID
					? baseline : FaultSystemSolution.load(solutionFile);
			ReportMetadata metadata = new ReportMetadata(new RupSetMetadata(calculation.displayName(), solution),
					calculation != Calculation.MID && baseline != null
							? new RupSetMetadata("NSHM23 MID", baseline) : null);
			ReportPageGen report = new ReportPageGen(metadata, new File(directory, "report"),
					ReportPageGen.getDefaultSolutionPlots(PlotLevel.REVIEW));
			report.setReplot(true);
			report.generatePage();
		}
	}

	private static void writeIndex(File output, int threads, long iterations) throws IOException {
		List<String> lines = new ArrayList<>(List.of("# Utah normal-fault segmentation comparison", "",
				"Common branch: `"+NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.buildFileName()+"`, restricted with "
				+"`NSHM23_SingleStates.UT` (retains Utah-connected neighboring sections). "
				+"Only segmentation changes.", "",
				"Threads: "+threads+"; iterations per rupture: "+iterations+"; "
				+"standard cluster-specific NSHM23 solver. "
				+"The factory also applies its usual minimum iteration budget based on section count.", "",
				"Geometric choices use the 10-degree bend threshold, normal-to-normal geometry, and MID background "
				+"without hardcoded Wasatch rules. Original MID retains its production Wasatch rules.", "",
				"If a subsection has no eligible rupture that can carry rate by itself under the original segmentation "
				+"bounds, all of its minimum-size ruptures (by subsection count) are exempt from internal bend bounds "
				+"on that subsection's side. Their participation still counts toward the reference rate. Larger "
				+"crossings and between-parent jumps remain constrained. MFD adjustments use the same policy.", "",
				"Each `config.json` is the regional inspection configuration; "
				+"cluster solves build subset configurations. "
				+"Independent stochastic realizations can differ even with matched calculation settings.", "",
				"Launch `scratch.kevin.segmentation.normalFault.UtahSegmentationInversionRunner` with no arguments "
				+"to run all three cases. Add `--prepare-only` to write inputs, "
				+"or `--case MID|MANUAL_60K|AUTOMATED_5KM` "
				+"to select a case. Use `--output-dir` to continue in this directory on a later date, "
				+"or `--index-only` to refresh this index without preparing or running calculations.", ""));
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Case", "Status", "Inputs", "Solution", "Report");
		for (Calculation calculation : Calculation.values()) {
			String path = calculation.directoryName();
			File directory = new File(output, path);
			boolean solved = new File(directory, "solution.zip").exists();
			boolean prepared = new File(directory, "rupSet.zip").exists();
			boolean report = new File(directory, "report/README.md").exists();
			table.addLine(calculation.displayName(), solved ? "Completed" : prepared ? "Prepared" : "Pending",
					prepared ? "[Rupture set]("+path+"/rupSet.zip), [config]("+path+"/config.json), "
							+"[exemptions]("+path+"/minimum_rupture_exemptions.csv)" : "",
					solved ? "[Solution]("+path+"/solution.zip)" : "",
					report ? "[Report]("+path+"/report/README.md)" : "");
		}
		lines.addAll(table.build());
		appendMisfitSummary(lines, output);
		MarkdownUtils.writeReadmeAndHTML(lines, output);
	}

	private static void appendMisfitSummary(List<String> lines, File output) throws IOException {
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Case", "Constraint", "Rows", "Mean z", "Mean absolute z", "RMS z", "Largest absolute z");
		for (Calculation calculation : Calculation.values()) {
			File solutionFile = new File(new File(output, calculation.directoryName()), "solution.zip");
			if (!solutionFile.isFile()) continue;
			InversionMisfitStats misfits = FaultSystemSolution.load(solutionFile).getModule(InversionMisfitStats.class);
			if (misfits == null) continue;
			for (MisfitStats stats : misfits.getStats()) {
				if (stats.numRows == 0 || stats.range.inequality
						|| stats.range.weightingType != ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)
					continue;
				String constraint = switch (stats.range.shortName) {
				case "UncertSlipRate" -> "Slip rate";
				case PaleoRateInversionConstraint.SHORT_NAME -> "Paleoseismic event rate";
				case PaleoSlipInversionConstraint.SHORT_NAME -> "Paleoseismic slip";
				default -> null;
				};
				if (constraint != null)
					table.addLine(calculation.displayName(), constraint, Integer.toString(stats.numRows),
							formatMisfit(stats.mean), formatMisfit(stats.absMean), formatMisfit(stats.rmse),
							formatMisfit(Math.max(Math.abs(stats.min), Math.abs(stats.max))));
			}
		}
		if (table.getNumLines() > 1) {
			lines.addAll(List.of("", "## Slip-rate and paleoseismic misfits", "",
					"Each z-score is the solution prediction minus the inversion target, divided by the "
					+"standard deviation used in that constraint. Overall inversion weights are removed. "
					+"Rows within each constraint contribute equally. Negative mean z indicates underprediction.", ""));
			lines.addAll(table.build());
		}
	}

	private static String formatMisfit(double value) {
		return String.format(Locale.ROOT, "%.3f", value);
	}
}
