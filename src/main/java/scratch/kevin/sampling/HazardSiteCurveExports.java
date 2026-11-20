package scratch.kevin.sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.statistics.descriptive.Quantile;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.CSVWriter;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SiteLogicTreeHazardCurveCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;

import com.google.common.base.Preconditions;

/** Exports selected map-grid sites from the paper's sampled runs into simplified site-curve archives. */
public class HazardSiteCurveExports {

	private static final List<ExportSite> SITES = List.of(
			new ExportSite("Manu'a Islands", new Location(-14.2, 190.4)),
			new ExportSite("Samoa", new Location(-13.6, 187.6)),
			new ExportSite("Tonga Deformation Front", new Location(-15.0, 187.6)));
	private static final double[] PERIODS = { 0d, 1d };
	private static final double[] FRACTILES = { 0.025, 0.16, 0.5, 0.84, 0.975 };
	private static final Quantile QUANTILE = Quantile.withDefaults();

	/**
	 * Optional arguments: output directory, then reader thread count. Defaults to site_curves beside the source runs
	 * and up to eight readers. This reads source archives only; no hazard calculations are performed.
	 */
	public static void main(String[] args) throws IOException {
		Preconditions.checkArgument(args.length <= 2,
				"Usage: HazardSiteCurveExports [output-directory] [reader-threads]");
		File outputDir = args.length > 0 ? new File(args[0]) : new File(PaperPaths.INVS_SUB_DIR, "site_curves");
		int threads = args.length > 1 ? Integer.parseInt(args[1])
				: Math.min(8, Runtime.getRuntime().availableProcessors());
		List<File> mcsRuns = new ArrayList<>();
		for (int count : HazardConvergenceCalcs.runDirs.row(SamplingMethod.MONTE_CARLO)
				.keySet().stream().sorted().toList())
			mcsRuns.addAll(HazardConvergenceCalcs.runDirs.get(SamplingMethod.MONTE_CARLO, count));
		exportGroup(outputDir, "mcs_pool", mcsRuns, SITES, PERIODS, threads);
		exportGroup(outputDir, "mcs_reserve", List.of(HazardConvergenceCalcs.MCS_RESERVE_DIR),
				SITES, PERIODS, threads);
		File sobolRun = HazardConvergenceCalcs.runDirs.get(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 16384).getFirst();
		exportGroup(outputDir, sobolRun.getName(), List.of(sobolRun), SITES, PERIODS, threads);
	}

	static record ExportSite(String name, Location location) {
		String prefix() {
			return FileNameUtils.simplify(name);
		}
	}

	private record SourceRun(File directory, int size, int firstSampleIndex, SamplingMethod method, long seed) {}
	private record SelectedCurves(double[] x, double[][] y) {}

	static void exportGroup(File outputDir, String prefix, List<File> runDirs, List<ExportSite> sites,
			double[] periods, int threads) throws IOException {
		Preconditions.checkArgument(!runDirs.isEmpty() && !sites.isEmpty() && periods.length > 0);
		Preconditions.checkArgument(threads > 0, "Reader thread count must be positive");
		Files.createDirectories(outputDir.toPath());
		Preconditions.checkArgument(runDirs.stream().map(File::getAbsolutePath).distinct().count() == runDirs.size(),
				"Duplicate source runs");
		Preconditions.checkArgument(sites.stream().map(ExportSite::prefix).distinct().count() == sites.size(),
				"Site filename prefixes must be unique");
		GriddedRegion region = loadRegion(runDirs.getFirst());
		int[] indexes = resolveSites(region, sites);
		List<SourceRun> runs = new ArrayList<>();
		List<String> levelNames = null;
		int total = 0;
		for (File dir : runDirs) {
			resolveSites(loadRegion(dir), sites);
			LogicTree<?> tree = loadTree(dir);
			List<String> names = tree.getLevels().stream().map(LogicTreeLevel::getShortName).toList();
			if (levelNames == null)
				levelNames = names;
			Preconditions.checkState(levelNames.equals(names), "Analysis-tree levels differ in %s", dir);
			Preconditions.checkState(tree.size() > 0, "Empty analysis tree: %s", dir);
			for (int b = 0; b < tree.size(); b++)
				Preconditions.checkState(Math.abs(tree.getBranchWeight(b) * tree.size() - 1d) < 1e-6,
						"Expected equally weighted samples: %s branch %s", dir, b);
			runs.add(new SourceRun(dir, tree.size(), total, tree.getSamplingMethod(), tree.getSamplingRandomSeed()));
			total = Math.addExact(total, tree.size());
		}
		System.out.println("Exporting " + prefix + ": " + total + " samples from " + runs.size() + " runs");
		File fullFile = new File(outputDir, prefix + "_curves.zip");
		File summaryFile = new File(outputDir, prefix + "_summaries.zip");
		File fullTemp = File.createTempFile(prefix + "_curves_", ".tmp", outputDir);
		File summaryTemp = File.createTempFile(prefix + "_summaries_", ".tmp", outputDir);
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		try {
			try (ZipOutputStream full = new ZipOutputStream(Files.newOutputStream(fullTemp.toPath()));
					ZipOutputStream summary = new ZipOutputStream(Files.newOutputStream(summaryTemp.toPath()))) {
				for (ZipOutputStream output : List.of(full, summary))
					writeMetadata(output, runs, sites, indexes, levelNames, total);
				for (double period : periods) {
					double[] x = null;
					double[][][] values = null;
					for (SourceRun run : runs) {
						System.out.println("  Loading period " + period + ": " + run.directory().getName());
						LogicTree<?> tree = loadTree(run.directory());
						GriddedRegion runRegion = loadRegion(run.directory());
						int[] runIndexes = resolveSites(runRegion, sites);
						try (SelectedCurveReader reader = new SelectedCurveReader(run.directory(), runRegion, sites,
								runIndexes, period)) {
							// Bound outstanding tasks and preserve branch ordering regardless of completion order.
							for (int start = 0; start < run.size(); start += 128) {
								List<Future<SelectedCurves>> pending = new ArrayList<>();
								int end = Math.min(run.size(), start + 128);
								for (int b = start; b < end; b++) {
									LogicTreeBranch<?> branch = tree.getBranch(b);
									pending.add(executor.submit(() -> reader.load(branch)));
								}
								for (int b = start; b < end; b++) {
									SelectedCurves curves = await(pending.get(b - start));
									if (x == null) {
										x = curves.x();
										values = new double[sites.size()][x.length][total];
									}
									Preconditions.checkState(Arrays.equals(x, curves.x()),
											"Curve intensity grids differ in %s branch %s", run.directory(), b);
									for (int s = 0; s < sites.size(); s++)
										for (int i = 0; i < x.length; i++)
											values[s][i][run.firstSampleIndex() + b] = curves.y()[s][i];
								}
							}
						}
					}
					for (int s = 0; s < sites.size(); s++) {
						String entry = MPJ_SiteLogicTreeHazardCurveCalc.getSitePeriodPrefix(sites.get(s).prefix(), period)
								+ ".csv";
						writeCurves(full, entry, x, values[s], total);
						writeSummary(summary, entry, period, x, values[s]);
					}
				}
			}
			Files.move(fullTemp.toPath(), fullFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			Files.move(summaryTemp.toPath(), summaryFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			System.out.println("Wrote " + fullFile + " and " + summaryFile);
		} finally {
			executor.shutdownNow();
			Files.deleteIfExists(fullTemp.toPath());
			Files.deleteIfExists(summaryTemp.toPath());
		}
	}

	private static SelectedCurves await(Future<SelectedCurves> future) throws IOException {
		try {
			return future.get();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted reading site curves", e);
		} catch (ExecutionException e) {
			throw new IOException("Failed reading site curves", e.getCause());
		}
	}

	private static LogicTree<?> loadTree(File runDir) throws IOException {
		return LogicTree.read(new File(runDir, "logic_tree_analysis.json"));
	}

	private static GriddedRegion loadRegion(File runDir) throws IOException {
		return GriddedRegion.fromFeature(Feature.read(new File(runDir, "gridded_region.geojson")));
	}

	private static int[] resolveSites(GriddedRegion region, List<ExportSite> sites) {
		int[] indexes = new int[sites.size()];
		for (int s = 0; s < sites.size(); s++) {
			ExportSite site = sites.get(s);
			indexes[s] = region.indexForLocation(site.location());
			Preconditions.checkState(indexes[s] >= 0
					&& LocationUtils.areSimilar(site.location(), region.locationForIndex(indexes[s])),
					"Site %s (%s) is not an exact grid location", site.name(), site.location());
		}
		return indexes;
	}

	private static CSVWriter writer(ZipOutputStream zip, String entry, List<String> header) throws IOException {
		zip.putNextEntry(new ZipEntry(entry));
		CSVWriter writer = new CSVWriter(zip, true);
		writer.write(header);
		return writer;
	}

	private static void finishEntry(ZipOutputStream zip, CSVWriter writer) throws IOException {
		writer.flush();
		zip.closeEntry();
	}

	private static void writeMetadata(ZipOutputStream zip, List<SourceRun> runs, List<ExportSite> sites,
			int[] indexes, List<String> levels, int total) throws IOException {
		CSVWriter writer = writer(zip, "sites.csv",
				List.of("Site Name", "Latitude", "Longitude", "Site Filename Prefix"));
		for (ExportSite site : sites)
			writer.write(List.of(site.name(), Double.toString(site.location().lat),
					Double.toString(site.location().lon), site.prefix()));
		finishEntry(zip, writer);
		writer = writer(zip, "site_grid_indexes.csv", List.of("Site Filename Prefix", "First Run Grid Index"));
		for (int s = 0; s < sites.size(); s++)
			writer.write(List.of(sites.get(s).prefix(), Integer.toString(indexes[s])));
		finishEntry(zip, writer);
		writer = writer(zip, "runs.csv", List.of("Run ID", "Run Directory", "Sample Count",
				"First Sample Index", "Sampling Method", "Sampling Random Seed"));
		for (SourceRun run : runs)
			writer.write(List.of(run.directory().getName(), run.directory().getAbsolutePath(),
					Integer.toString(run.size()), Integer.toString(run.firstSampleIndex()),
					run.method() == null ? "" : run.method().name(), Long.toString(run.seed())));
		finishEntry(zip, writer);
		writer = writer(zip, "samples.csv", List.of("Branch Index", "Run ID", "Original Branch Index",
				"Original Branch Filename"));
		for (SourceRun run : runs) {
			LogicTree<?> tree = loadTree(run.directory());
			for (int b = 0; b < tree.size(); b++)
				writer.write(List.of(Integer.toString(run.firstSampleIndex() + b), run.directory().getName(),
						Integer.toString(b), tree.getBranch(b).buildFileName()));
		}
		finishEntry(zip, writer);
		List<String> header = new ArrayList<>(List.of("Branch Index", "Branch Weight"));
		header.addAll(levels);
		writer = writer(zip, "logic_tree.csv", header);
		for (SourceRun run : runs) {
			LogicTree<?> tree = loadTree(run.directory());
			for (int b = 0; b < tree.size(); b++) {
				List<String> row = new ArrayList<>();
				row.add(Integer.toString(run.firstSampleIndex() + b));
				row.add(Double.toString(1d / total));
				for (LogicTreeNode node : tree.getBranch(b))
					row.add(node.getShortName());
				writer.write(row);
			}
		}
		finishEntry(zip, writer);
		zip.putNextEntry(new ZipEntry("README.txt"));
		zip.write(("Samples are equally weighted; branch indexes are zero-based and follow runs.csv order.\n"
				+ "logic_tree.csv contains the source analysis-tree node short names. No merged JSON tree is constructed.\n"
				+ "samples.csv maps every output row to its source run and branch.\n"
				+ "Curves retain the native intensity grid and source exceedance-probability units and duration.\n"
				+ "Summaries are pointwise curve statistics, not statistics of inverted hazard-map values.\n"
				+ "Percentiles use Apache Commons Statistics Quantile HF8.\n"
				+ "For a single Sobol run, the first N curve rows preserve its N-point prefix; reweight to 1/N.\n")
				.getBytes(StandardCharsets.UTF_8));
		zip.closeEntry();
	}

	private static void writeCurves(ZipOutputStream zip, String entry, double[] x, double[][] values,
			int total) throws IOException {
		List<String> header = new ArrayList<>(List.of("Branch Index", "Branch Weight"));
		for (double intensity : x)
			header.add(Double.toString(intensity));
		CSVWriter writer = writer(zip, entry, header);
		String weight = Double.toString(1d / total);
		for (int b = 0; b < total; b++) {
			List<String> row = new ArrayList<>(header.size());
			row.add(Integer.toString(b));
			row.add(weight);
			for (double[] atIntensity : values)
				row.add(Double.toString(atIntensity[b]));
			writer.write(row);
		}
		finishEntry(zip, writer);
	}

	private static void writeSummary(ZipOutputStream zip, String entry, double period, double[] x,
			double[][] values) throws IOException {
		String label = period == 0d ? "PGA (g)" : (float)period + "s SA (g)";
		CSVWriter writer = writer(zip, entry,
				List.of(label, "Mean", "Median", "Min", "p2.5", "p16", "p84", "p97.5", "Max"));
		for (int i = 0; i < x.length; i++) {
			double[] atIntensity = values[i];
			double mean = 0d;
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;
			for (double value : atIntensity) {
				mean += value;
				min = Math.min(min, value);
				max = Math.max(max, value);
			}
			mean /= atIntensity.length;
			// Full curves have already been written; quantile selection can safely rearrange this array.
			double[] q = QUANTILE.evaluate(atIntensity, FRACTILES);
			writer.write(List.of(Double.toString(x[i]), Double.toString(mean), Double.toString(q[2]),
					Double.toString(min), Double.toString(q[0]), Double.toString(q[1]),
					Double.toString(q[3]), Double.toString(q[4]), Double.toString(max)));
		}
		finishEntry(zip, writer);
	}

	private static final class SelectedCurveReader implements AutoCloseable {
		private final File runDir;
		private final GriddedRegion region;
		private final List<ExportSite> sites;
		private final int[] indexes;
		private final double period;
		private final ZipFile archive;

		SelectedCurveReader(File runDir, GriddedRegion region, List<ExportSite> sites, int[] indexes,
				double period) throws IOException {
			this.runDir = runDir;
			this.region = region;
			this.sites = sites;
			this.indexes = indexes;
			this.period = period;
			File file = new File(runDir, "results_hazard_curves.zip");
			archive = file.isFile() ? new ZipFile(file) : null;
		}

		SelectedCurves load(LogicTreeBranch<?> branch) throws IOException {
			String filename = SolHazardMapCalc.getCSV_FileName("curves", period);
			InputStream raw;
			boolean gzip;
			if (archive != null) {
				String prefix = branch.buildFileName() + "/" + filename;
				ZipEntry entry = archive.getEntry(prefix + ".gz");
				if (entry == null)
					entry = archive.getEntry(prefix);
				Preconditions.checkNotNull(entry, "Missing %s[.gz] in %s", prefix, archive.getName());
				raw = archive.getInputStream(entry);
				gzip = entry.getName().endsWith(".gz");
			} else {
				File directory = new File(branch.getBranchDirectory(new File(runDir, "results"), false),
						"hazard_" + (float)region.getSpacing() + "deg_grid_seis_INCLUDE");
				File file = new File(directory, filename);
				if (!file.isFile())
					file = new File(directory, filename + ".gz");
				raw = new FileInputStream(file);
				gzip = file.getName().endsWith(".gz");
			}
			try (InputStream stream = raw;
					BufferedReader reader = new BufferedReader(new InputStreamReader(
							gzip ? new GZIPInputStream(stream) : stream, StandardCharsets.UTF_8))) {
				String line = reader.readLine();
				Preconditions.checkNotNull(line, "Empty curve CSV for %s", branch);
				List<String> header = CSVFile.loadLine(line);
				double[] x = new double[header.size() - 3];
				for (int i = 0; i < x.length; i++) {
					x[i] = Double.parseDouble(header.get(i + 3));
					Preconditions.checkState(Double.isFinite(x[i]) && x[i] > 0d
							&& (i == 0 || x[i] > x[i - 1]), "Invalid intensity grid for %s", branch);
				}
				Preconditions.checkState(x.length > 0, "No intensities for %s", branch);
				double[][] y = new double[sites.size()][];
				int found = 0;
				while ((line = reader.readLine()) != null) {
					// Only parse the index on unselected rows, avoiding their many numeric curve columns.
					int comma = line.indexOf(',');
					Preconditions.checkState(comma > 0, "Malformed curve row for %s", branch);
					int index = Integer.parseInt(line.substring(0, comma).trim());
					for (int s = 0; s < sites.size(); s++) {
						if (index != indexes[s])
							continue;
						Preconditions.checkState(y[s] == null, "Duplicate grid index %s for %s", index, branch);
						List<String> row = CSVFile.loadLine(line);
						Preconditions.checkState(row.size() == header.size(), "Curve row size mismatch for %s", branch);
						Location location = new Location(Double.parseDouble(row.get(1)), Double.parseDouble(row.get(2)));
						Preconditions.checkState(LocationUtils.areSimilar(location, sites.get(s).location()),
								"Site coordinate mismatch for %s in %s", sites.get(s).name(), branch);
						y[s] = new double[x.length];
						for (int i = 0; i < x.length; i++) {
							double value = Double.parseDouble(row.get(i + 3));
							Preconditions.checkState(Double.isFinite(value) && value >= 0d && value <= 1d,
									"Invalid probability %s for %s", value, branch);
							y[s][i] = value;
						}
						found++;
					}
					if (found == sites.size())
						break;
				}
				Preconditions.checkState(found == sites.size(), "Missing selected sites for %s", branch);
				return new SelectedCurves(x, y);
			}
		}

		@Override
		public void close() throws IOException {
			if (archive != null)
				archive.close();
		}
	}
}
