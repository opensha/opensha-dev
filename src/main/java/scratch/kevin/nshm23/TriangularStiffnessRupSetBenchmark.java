package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.CoulombRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.RupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.RupSetSubsectioningModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchGenerationStrategy;

import com.google.gson.GsonBuilder;

/**
 * Standalone, cold-cache NSHM23 rupture-set experiment. Run each arm in a fresh
 * JVM with the same resource limits. The extra rectangular/cap=0 arm separates
 * the effect of removing the rectangular self-stiffness cap from triangulation.
 * No computational caches are read or written; ordinary input data are retained.
 */
public class TriangularStiffnessRupSetBenchmark extends NSHM23_InvConfigFactory {

	private final Path output;
	private final PatchGenerationStrategy strategy;
	private final double cap;
	private final Map<String, Object> metrics = new LinkedHashMap<>();
	private SubSectStiffnessCalculator stiffness;

	private TriangularStiffnessRupSetBenchmark(Path output, PatchGenerationStrategy strategy, double cap) {
		this.output = output;
		this.strategy = strategy;
		this.cap = cap;
		setCacheDir(null);
		setAutoCache(false);
	}

	@Override
	protected synchronized FaultSystemRupSet buildGenericRupSet(LogicTreeBranch<?> branch, int threads) {
		RupSetFaultModel fm = branch.requireValue(RupSetFaultModel.class);
		try {
			List<? extends FaultSection> sections = fm.getDefaultDeformationModel().build(fm,
					branch.requireValue(RupSetSubsectioningModel.class), branch);
			GeoJSONFaultReader.writeFaultSections(output.resolve("input_subsections.geojson").toFile(), sections);
			metrics.put("subsections", sections.size());
			stiffness = new SubSectStiffnessCalculator(sections, 2d, 3e4, 3e4, 0.5d,
					strategy, cap);
			long start = System.nanoTime();
			int patches = 0;
			try (var writer = Files.newBufferedWriter(output.resolve("patch_counts.csv"))) {
				writer.write("section_id,patches\n");
				for (FaultSection section : sections) {
					int count = stiffness.getPatches(section).size();
					patches += count;
					writer.write(section.getSectionId() + "," + count + "\n");
				}
			}
			metrics.put("patches", patches);
			metrics.put("patch_build_seconds", secondsSince(start));
			System.out.println("Prepared " + sections.size() + " sections and " + patches + " patches");
			RupSetConfig config = new CoulombRupSetConfig(sections, fm.getFilePrefix().toLowerCase(),
					branch.requireValue(RupSetScalingRelationship.class)) {
				@Override public SubSectStiffnessCalculator getStiffnessCalc() { return stiffness; }
				// The base implementation supplies a global default for null, so override
				// the accessor as well as disabling automatic cache writes.
				@Override protected File getCacheDir() { return null; }
			};
			config.setAutoCache(false);
			// Reuse this fresh in-memory calculator in the factory's branch update.
			// This also avoids its fallback disk-cache lookup.
			distAzCache.put(fm, config.getDistAzCalc());
			start = System.nanoTime();
			FaultSystemRupSet result = config.build(threads);
			metrics.put("generic_build_seconds", secondsSince(start));
			metrics.put("cached_distances", config.getDistAzCalc().getNumCachedDistances());
			metrics.put("cached_azimuths", config.getDistAzCalc().getNumCachedAzimuths());
			return result;
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	private static double secondsSince(long start) {
		return (System.nanoTime() - start) * 1e-9;
	}

	public static void main(String[] args) throws IOException {
		if (args.length != 3)
			throw new IllegalArgumentException("Usage: <legacy|triangles|rectangle_cap0> <output-dir> <threads>");
		Path output = Path.of(args[1]).toAbsolutePath();
		Files.createDirectories(output);
		if (Files.exists(output.resolve("rupture_set.zip")))
			throw new IllegalStateException("Refusing to overwrite a completed experiment: " + output);
		PatchGenerationStrategy strategy = switch (args[0]) {
			case "legacy", "rectangle_cap0" -> PatchGenerationStrategy.AUTO;
			case "triangles" -> PatchGenerationStrategy.TRIANGULAR;
			default -> throw new IllegalArgumentException("Unknown arm: " + args[0]);
		};
		int threads = Integer.parseInt(args[2]);
		var experiment = new TriangularStiffnessRupSetBenchmark(output, strategy, args[0].equals("legacy") ? 1d : 0d);
		var metrics = experiment.metrics;
		metrics.put("arm", args[0]);
		metrics.put("start_utc", Instant.now().toString());
		metrics.put("java_version", System.getProperty("java.version"));
		metrics.put("jvm_arguments", ManagementFactory.getRuntimeMXBean().getInputArguments());
		metrics.put("threads", threads);
		metrics.put("strategy", strategy.name());
		metrics.put("self_stiffness_cap", experiment.cap);
		metrics.put("grid_spacing_km", 2d);
		metrics.put("disk_cache_reads", false);
		metrics.put("disk_cache_writes", false);
		long start = System.nanoTime();
		try {
			var branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
			metrics.put("branch", branch.toString());
			FaultSystemRupSet rupSet = experiment.buildRuptureSet(branch, threads);
			metrics.put("build_seconds", secondsSince(start));
			metrics.put("ruptures", rupSet.getNumRuptures());
			metrics.put("min_magnitude", rupSet.getMinMag());
			metrics.put("max_magnitude", rupSet.getMaxMag());
			rupSet.requireModule(PlausibilityConfiguration.class).writeJSON(output.resolve("plausibility.json").toFile());
			long writeStart = System.nanoTime();
			rupSet.write(output.resolve("rupture_set.zip").toFile());
			metrics.put("write_seconds", secondsSince(writeStart));
			metrics.put("zip_bytes", Files.size(output.resolve("rupture_set.zip")));
			metrics.put("success", true);
		} catch (Throwable failure) {
			metrics.put("success", false);
			metrics.put("failure", failure.toString());
			failure.printStackTrace();
		} finally {
			metrics.put("total_seconds", secondsSince(start));
			metrics.put("finish_utc", Instant.now().toString());
			if (ManagementFactory.getOperatingSystemMXBean() instanceof com.sun.management.OperatingSystemMXBean os)
				metrics.put("process_cpu_seconds", os.getProcessCpuTime() * 1e-9);
			Files.writeString(output.resolve("metrics.json"), new GsonBuilder().setPrettyPrinting().create().toJson(metrics));
		}
		if (!Boolean.TRUE.equals(metrics.get("success")))
			System.exit(1);
		Files.writeString(output.resolve("COMPLETE"), Instant.now() + "\n");
		System.out.println("EXPERIMENT COMPLETE: " + output);
		// Some library executors are persistent; all results are closed at this point.
		System.exit(0);
	}
}
