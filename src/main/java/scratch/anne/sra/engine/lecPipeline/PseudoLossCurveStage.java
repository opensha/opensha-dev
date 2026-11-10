package scratch.anne.sra.engine.lecPipeline;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.domain.risk.LossExceedanceCurve;
import scratch.anne.sra.engine.logicTree.StreamingBranchLecBuilder;
import scratch.anne.sra.io.readers.ParquetMetadataReader;
import scratch.anne.sra.io.readers.RuptureLossLookupReader;
import scratch.anne.sra.io.readers.RuptureLossLookupReader.RuptureAttributes;
import scratch.anne.sra.io.writers.PseudoLossCurveExporter;

/**
 * Pipeline stage 3: full discrete-rank pseudo loss curves, streamed from
 * stage 2's branch contributions.
 *
 * <p>Replicates {@code PseudoLossCurveExportExample}'s body. No existing
 * class is modified; this class only calls the same public readers/
 * writers the original script already used.</p>
 *
 * <p>Gated by existence: if BOTH {@link PipelineConfig#getPseudoLossCurveCsv}
 * (for {@code config.getCsvSampleBranches()}) and
 * {@link PipelineConfig#getPseudoLossCurveParquet} (covering every
 * branch) already exist, this stage is skipped -- no SolutionLogicTree
 * load or Parquet streaming needed just to check.</p>
 *
 * <p>{@code csvSampleBranches} only ever limits the CSV output; the
 * Parquet output always covers every branch in the contributions file,
 * regardless of that set -- see {@code PipelineConfig}'s class javadoc
 * on the {@code _subset-n-N} vs {@code _sample-k-M} naming distinction.</p>
 */
public final class PseudoLossCurveStage {

    private PseudoLossCurveStage() {}

    public static void run(PipelineConfig config) throws SQLException, IOException {

        Path outputCsv = config.getPseudoLossCurveCsv();
        Path outputParquet = config.getPseudoLossCurveParquet();

        if (Files.exists(outputCsv) && Files.exists(outputParquet)) {
            System.out.println("[Stage 3] Skipping -- already exist: " + outputCsv + ", " + outputParquet);
            return;
        }

        System.out.println("[Stage 3] Not found -- running pseudo loss curve export...");

        Path inputCsv = config.getPerRuptureEalCsv();
        Path contributionsParquet = config.getBranchContributionsParquet();

        Map<RuptureKey, Double> fixedLossLookup = RuptureLossLookupReader.read(inputCsv);
        Map<RuptureKey, RuptureAttributes> ruptureAttributes = RuptureLossLookupReader.readRuptureAttributes(inputCsv);

        // propagate the SAME provenance (numBranches, solutionLogicTreeZip,
        // fileTag, full HazardParameters) that BranchContributionsStage
        // already wrote into contributionsParquet's own file-level
        // metadata -- rather than independently reconstructing it (which
        // would require reloading the SolutionLogicTree, exactly what
        // this stage is otherwise careful to avoid)
        Map<String, String> contribMetadata = ParquetMetadataReader.readKeyValueMetadata(contributionsParquet);
        String metadataComment = contribMetadata.isEmpty() ? null : formatMetadataComment(contribMetadata);

        Set<Integer> targetBranches = config.getCsvSampleBranches();

        System.out.println("[Stage 3] Extracting full curves for branches: " + targetBranches
                + " (CSV only -- Parquet output always covers every branch)");

        // gather each requested branch's curve fully in memory -- fine
        // for a handful of branches (the wide CSV output needs the full
        // set at once to compute its shared union loss axis), unlike a
        // full-scale export. The Parquet output is written per branch
        // as curves stream past, below, and is never collected here, so
        // it stays flat regardless of total branch count.
        Map<Integer, LossExceedanceCurve> curvesByBranch = new HashMap<>();

        try (PseudoLossCurveExporter exporter = new PseudoLossCurveExporter(outputParquet, contribMetadata)) {
            StreamingBranchLecBuilder.streamDiscreteRankCurves(contributionsParquet, fixedLossLookup, (branchIndex, curve) -> {
                // always written, for every branch, regardless of targetBranches
                exporter.writeParquet(branchIndex, curve, ruptureAttributes);

                if (targetBranches.contains(branchIndex))
                    curvesByBranch.put(branchIndex, curve);
            });

            // wide/pivoted CSV, limited to targetBranches -- carries the
            // same cumulative-likelihood information as the Parquet
            // output for these branches, human-readable, with
            // sourceId/ruptureId/fltSysRupIndex/magnitude attributed per
            // row (see PseudoLossCurveExporter)
            exporter.writeCsv(curvesByBranch, ruptureAttributes, metadataComment, outputCsv);
        }

        System.out.println("[Stage 3] Wrote curves CSV (branches " + targetBranches + ") to " + outputCsv);
        System.out.println("[Stage 3] Wrote per-rupture curve keys Parquet (all branches) to " + outputParquet);
    }

    /**
     * Flattens a key-value map into one {@code key=val,key2=val2,...}
     * line, for the CSV comment line -- same shape as
     * {@code BranchContributionsStage}'s own inline metadata strings,
     * just built from an already-flattened map instead of individually
     * appended fields.
     */
    private static String formatMetadataComment(Map<String, String> metadata) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, String> e : metadata.entrySet()) {
            if (sb.length() > 0)
                sb.append(',');
            sb.append(e.getKey()).append('=').append(e.getValue());
        }
        return sb.toString();
    }
}
