package scratch.kevin.sampling;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.sampling.ArrayPointSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeLevel.FileBackedLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.FileBackedNode;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SiteLogicTreeHazardCurveCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;

import scratch.kevin.sampling.HazardSiteCurveExports.ExportSite;

public class HazardSiteCurveExportsTest {

	@Rule
	public TemporaryFolder temporary = new TemporaryFolder();

	@Test
	public void pooledExportPreservesSampleOrderAndAnalysisLabels() throws IOException {
		GriddedRegion region = new GriddedRegion(new Location(-14.2, 190.4), new Location(-14.0, 190.6),
				0.2, GriddedRegion.ANCHOR_0_0);
		List<ExportSite> sites = List.of(new ExportSite("First Site", region.locationForIndex(0)),
				new ExportSite("Last Site", region.locationForIndex(region.getNodeCount() - 1)));
		File first = sourceRun("first", 0, region, true);
		File second = sourceRun("second", 4, region, false);
		File output = temporary.newFolder("output");
		HazardSiteCurveExports.exportGroup(output, "pooled", List.of(first, second), sites,
				new double[] { 0d, 1d }, 2);
		try (ZipFile full = new ZipFile(new File(output, "pooled_curves.zip"));
				ZipFile summary = new ZipFile(new File(output, "pooled_summaries.zip"))) {
			assertNull(full.getEntry("logic_tree.json"));
			CSVFile<String> logic = read(full, "logic_tree.csv");
			assertEquals(9, logic.getNumRows());
			assertEquals("Choice", logic.get(0, 2));
			assertEquals("Choice0", logic.get(1, 2));
			assertEquals("Choice3", logic.get(8, 2));
			CSVFile<String> samples = read(full, "samples.csv");
			assertEquals("second", samples.get(5, 1));
			assertEquals(0, samples.getInt(5, 2));
			for (double period : new double[] { 0d, 1d }) {
				for (ExportSite site : sites) {
					int index = region.indexForLocation(site.location());
					String entry = MPJ_SiteLogicTreeHazardCurveCalc.getSitePeriodPrefix(site.prefix(), period) + ".csv";
					CSVFile<String> curves = read(full, entry);
					CSVFile<String> stats = read(summary, entry);
					assertEquals(9, curves.getNumRows());
					for (int b = 0; b < 8; b++) {
						assertEquals(b, curves.getInt(b + 1, 0));
						assertEquals(0.125, curves.getDouble(b + 1, 1), 0d);
						assertEquals(probability(b, index, period), curves.getDouble(b + 1, 2), 0d);
					}
					double shift = 0.001 * index + 0.1 * period;
					assertEquals(0.045 + shift, stats.getDouble(1, 1), 1e-15);
					assertEquals(0.045 + shift, stats.getDouble(1, 2), 1e-15);
					assertEquals(0.01 + shift, stats.getDouble(1, 3), 1e-15);
					assertEquals(0.01 + shift, stats.getDouble(1, 4), 1e-15);
					// HF8 at p=0.16 for eight evenly spaced observations: rank 1.666... (one-based).
					assertEquals(1d / 60d + shift, stats.getDouble(1, 5), 1e-15);
					assertEquals(0.08 + shift, stats.getDouble(1, 8), 1e-15);
				}
			}
		}
	}

	@Test(expected = IllegalStateException.class)
	public void rejectsOffGridSite() throws IOException {
		GriddedRegion region = new GriddedRegion(new Location(-14.2, 190.4), new Location(-14.0, 190.6),
				0.2, GriddedRegion.ANCHOR_0_0);
		File run = sourceRun("run", 0, region, true);
		HazardSiteCurveExports.exportGroup(temporary.newFolder("bad"), "bad", List.of(run),
				List.of(new ExportSite("Off grid", new Location(-14.15, 190.4))), new double[] { 0d }, 1);
	}

	private File sourceRun(String name, int offset, GriddedRegion region, boolean archived) throws IOException {
		File run = temporary.newFolder(name);
		if (!archived)
			Files.createDirectory(new File(run, "results").toPath());
		Feature.write(region.toFeature(), new File(run, "gridded_region.geojson"));
		List<FileBackedNode> nodes = new ArrayList<>();
		for (int b = 0; b < 4; b++)
			nodes.add(new FileBackedNode("Choice " + b, "Choice" + b, 0.25, "choice" + b));
		FileBackedLevel level = new FileBackedLevel("Choice", "Choice", nodes);
		LogicTree<LogicTreeNode> tree = LogicTree.buildSampled(List.of(level),
				new ArrayPointSet(new double[][] { { 0.1 }, { 0.3 }, { 0.6 }, { 0.9 } }));
		tree.write(new File(run, "logic_tree_analysis.json"));
		try (ZipOutputStream zip = archived ? new ZipOutputStream(
				Files.newOutputStream(new File(run, "results_hazard_curves.zip").toPath())) : null) {
			for (int b = 0; b < 4; b++) {
				for (double period : new double[] { 0d, 1d }) {
					StringBuilder csv = new StringBuilder("Index,Latitude,Longitude,0.1,1.0\n");
					for (int n = 0; n < region.getNodeCount(); n++) {
						Location loc = region.locationForIndex(n);
						double y = probability(offset + b, n, period);
						csv.append(n).append(',').append(loc.lat).append(',').append(loc.lon)
								.append(',').append(y).append(',').append(y / 2d).append('\n');
					}
					boolean gzip = b % 2 == 0;
					byte[] bytes = csv.toString().getBytes(StandardCharsets.UTF_8);
					if (gzip) {
						ByteArrayOutputStream compressed = new ByteArrayOutputStream();
						try (GZIPOutputStream gz = new GZIPOutputStream(compressed)) {
							gz.write(bytes);
						}
						bytes = compressed.toByteArray();
					}
					String fileName = SolHazardMapCalc.getCSV_FileName("curves", period) + (gzip ? ".gz" : "");
					if (archived) {
						zip.putNextEntry(new ZipEntry(tree.getBranch(b).buildFileName() + "/" + fileName));
						zip.write(bytes);
						zip.closeEntry();
					} else {
						File dir = new File(tree.getBranch(b).getBranchDirectory(new File(run, "results"), true),
								"hazard_0.2deg_grid_seis_INCLUDE");
						Files.createDirectories(dir.toPath());
						Files.write(new File(dir, fileName).toPath(), bytes);
					}
				}
			}
		}
		return run;
	}

	private static double probability(int sample, int site, double period) {
		return 0.01 * (sample + 1) + 0.001 * site + 0.1 * period;
	}

	private static CSVFile<String> read(ZipFile zip, String entry) throws IOException {
		try (var stream = zip.getInputStream(zip.getEntry(entry))) {
			return CSVFile.readStream(stream, true);
		}
	}
}
