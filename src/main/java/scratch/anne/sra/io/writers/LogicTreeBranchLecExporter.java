package scratch.anne.sra.io.writers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import scratch.anne.sra.util.StringUtil;

/**
 * Writes the per-branch LEC summary CSV: the existing branch-eals.csv
 * content plus one loss column per requested likelihood level.
 *
 * <p>Pure writer -- no Parquet/streaming knowledge, no percentile or
 * curve math. Takes the already-built per-branch loss lookup (from
 * streaming each branch's LossExceedanceCurve and calling
 * sampleAtLevels) and the original branch-eals.csv rows, and merges
 * them row by row. Same role as BranchEalSummaryExporter /
 * BranchLecPercentileSummaryExporter: a small, focused output class
 * with no upstream responsibilities.</p>
 */
public final class LogicTreeBranchLecExporter {

    private LogicTreeBranchLecExporter() {}

    /**
     * @param branchEALsCsv          existing per-branch summary CSV to extend (read, not modified)
     * @param likelihoods            likelihood levels, in the order matching each entry of lossByBranchAndLevel's arrays
     * @param lossByBranchAndLevel   branchIndex -> array of losses, one per likelihood level (NaN if that
     *                               branch never reached a given level -- see LossExceedanceCurve.sampleAtLevels)
     * @param outputCsv              destination file
     */
    public static void write(Path branchEALsCsv, double[] likelihoods,
                              Map<Integer, double[]> lossByBranchAndLevel, Path outputCsv) throws IOException {

        System.out.println("Writing per-branch summary to: " + outputCsv);

        try (var lines = Files.lines(branchEALsCsv);
             BufferedWriter w = Files.newBufferedWriter(outputCsv)) {

            List<String> allLines = lines.toList();

            // carry branchEALsCsv's own leading "# key=val,..." provenance
            // comment (if any -- see LogicTreeBranchEalSummaryExporter)
            // forward into this file unchanged, so a reader doesn't lose
            // that record just because this file extends the original
            int headerIndex = 0;
            if (!allLines.isEmpty() && allLines.get(0).startsWith("#")) {
                w.write(allLines.get(0));
                w.newLine();
                headerIndex = 1;
            }

            String header = allLines.get(headerIndex);

            StringBuilder newHeader = new StringBuilder(header);
            for (double level : likelihoods)
                newHeader.append(",loss@").append(StringUtil.formatLikelihood(level));
            w.write(newHeader.toString());
            w.newLine();

            int missingBranches = 0;
            for (int row = headerIndex + 1; row < allLines.size(); row++) {
                String line = allLines.get(row);
                int branchIndex = Integer.parseInt(line.split(",")[0].trim());

                double[] losses = lossByBranchAndLevel.get(branchIndex);

                StringBuilder newLine = new StringBuilder(line);
                if (losses == null) {
                    missingBranches++;
                    for (int i = 0; i < likelihoods.length; i++)
                        newLine.append(',');
                } else {
                    for (double loss : losses)
                        newLine.append(',').append(Double.isNaN(loss) ? "" : loss);
                }
                w.write(newLine.toString());
                w.newLine();
            }

            if (missingBranches > 0)
                System.out.println("WARNING: " + missingBranches + " branches from " + branchEALsCsv
                        + " had no rows in the streamed contributions");
        }
    }
}