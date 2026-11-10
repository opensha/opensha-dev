package scratch.anne.sra.io.writers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import scratch.anne.sra.engine.logicTree.BranchLecDistributionBuilder.LevelDistribution;
import scratch.anne.sra.util.StringUtil;

/**
 * Writes the ensemble-level percentile summary CSV: one row per metric
 * (EAL, plus one row per requested likelihood level), each with its own
 * weighted mean and weighted percentiles across all branches.
 *
 * <p>Pure writer -- no file reading, no Parquet/streaming knowledge, no
 * percentile math. Takes already-built LevelDistributions (from
 * BranchLecDistributionBuilder) and formats them. Same role as
 * BranchEalSummaryExporter / BranchContributionExporter: a small,
 * focused output class with no upstream responsibilities.</p>
 */
public final class LogicTreeBranchPercentileSummaryExporter {

    private LogicTreeBranchPercentileSummaryExporter() {}

    /**
     * @param ealDistribution        the EAL metric's own LevelDistribution, written as the first row
     * @param levelDistributions     one LevelDistribution per requested likelihood level, in the order to write
     * @param percentilesOfInterest  same array used to build the distributions -- controls column order
     * @param metadataComment        optional single-line {@code # key=val,...} provenance comment
     *                               written immediately before the header row -- see
     *                               {@link LogicTreeBranchEalSummaryExporter#write}'s javadoc for the
     *                               full rationale. Pass {@code null} or blank to omit it.
     */
    public static void write(LevelDistribution ealDistribution,
                              Map<Double, LevelDistribution> levelDistributions,
                              double[] percentilesOfInterest,
                              String metadataComment,
                              Path outputCsv) throws IOException {

        System.out.println("Writing ensemble percentile summary to: " + outputCsv);

        try (BufferedWriter w = Files.newBufferedWriter(outputCsv)) {
            if (metadataComment != null && !metadataComment.isBlank()) {
                w.write("# " + metadataComment);
                w.newLine();
            }

            StringBuilder header = new StringBuilder("likelihood,mean");
            for (double p : percentilesOfInterest)
                header.append(",p").append(StringUtil.formatPercentile(p));
            header.append(",missingCount,missingFraction");
            w.write(header.toString());
            w.newLine();

            writeRow("eal", ealDistribution, percentilesOfInterest, w);

            for (LevelDistribution dist : levelDistributions.values())
                writeRow(StringUtil.formatLikelihood(dist.requestedLevel()), dist, percentilesOfInterest, w);
        }
    }

    private static void writeRow(String label, LevelDistribution dist, double[] percentilesOfInterest,
                                  BufferedWriter w) throws IOException {
        StringBuilder line = new StringBuilder(label).append(',').append(dist.mean());
        for (double p : percentilesOfInterest)
            line.append(',').append(dist.percentiles().get(p));
        line.append(',').append(dist.missingCount()).append(',').append(dist.missingFraction());
        w.write(line.toString());
        w.newLine();
    }
}