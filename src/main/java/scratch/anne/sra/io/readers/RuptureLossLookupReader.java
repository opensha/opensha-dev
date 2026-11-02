package scratch.anne.sra.io.readers;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import scratch.anne.sra.domain.hazard.RuptureKey;

/**
 * Reads the per-rupture data block of a RuptureResultsExporter CSV
 * directly into a fixed {@code Map<RuptureKey, Double>} of ruptureLoss
 * values, WITHOUT reconstructing a full RuptureLossResultsCollection.
 *
 * <p>Per confirmed requirements: only ruptureLoss is needed from this
 * file. likelihood will be replaced by each branch's own rate in the
 * branch loop, and lossContribution/relativeContribution are derived
 * from likelihood*ruptureLoss -- neither is meaningful input once
 * likelihood is being swapped out. This reader intentionally ignores
 * those three columns rather than reconstructing a RuptureResult that
 * would need them.</p>
 *
 * <p>The output map plugs directly into {@link BranchEALCalculator}'s
 * constructor, in place of {@link BranchEALCalculator#buildFixedLossLookup}
 * (which builds the same shape from an in-memory, already-frozen
 * RuptureLossResultsCollection -- this is the equivalent for the
 * already-exported-to-disk case).</p>
 */
public final class RuptureLossLookupReader {

    private RuptureLossLookupReader() {}

    /**
     * Parses sourceId, ruptureId, ruptureLoss from the exported CSV's data
     * rows (skipping the "# HazardParameters" metadata block, the blank
     * lines, and the totalExpectedLoss line, exactly matching
     * RuptureResultsExporter.exportCsv's own layout).
     *
     * <p>Columns are resolved by name from the header row rather than by
     * fixed position, so this reader keeps working if
     * RuptureResultsExporter's column set or ordering changes (e.g. new
     * columns inserted before ruptureLoss).</p>
     *
     * @throws IllegalArgumentException if the expected data header row
     *         ("sourceId,ruptureId,...") is never found, or if a required
     *         column is missing from it
     */
    public static Map<RuptureKey, Double> read(Path exportedCsv) throws IOException {

        Map<RuptureKey, Double> lookup = new HashMap<>();

        try (BufferedReader br = Files.newBufferedReader(exportedCsv)) {

            String line;
            Map<String, Integer> colIndex = null;

            // Skip everything up through the data header row.
            while ((line = br.readLine()) != null) {
                if (line.trim().startsWith("sourceId")) {
                    colIndex = parseHeader(line);
                    break;
                }
            }

            if (colIndex == null) {
                throw new IllegalArgumentException(
                        "Could not find 'sourceId,ruptureId,...' data header row in: " + exportedCsv);
            }

            int sourceIdCol = requireColumn(colIndex, "sourceId", exportedCsv);
            int ruptureIdCol = requireColumn(colIndex, "ruptureId", exportedCsv);
            int ruptureLossCol = requireColumn(colIndex, "ruptureLoss", exportedCsv);
            int minColumns = Math.max(sourceIdCol, Math.max(ruptureIdCol, ruptureLossCol)) + 1;

            // Remaining lines are data rows, in whatever column order the
            // header row declared.
            while ((line = br.readLine()) != null) {
                if (line.isBlank())
                    continue;

                String[] tokens = line.split(",");
                if (tokens.length < minColumns) {
                    throw new IllegalArgumentException(
                            "Malformed data row (expected at least " + minColumns
                                    + " columns): '" + line + "' in " + exportedCsv);
                }

                int sourceId = Integer.parseInt(tokens[sourceIdCol].trim());
                int ruptureId = Integer.parseInt(tokens[ruptureIdCol].trim());
                double ruptureLoss = Double.parseDouble(tokens[ruptureLossCol].trim());

                lookup.put(new RuptureKey(sourceId, ruptureId), ruptureLoss);
            }
        }

        return lookup;
    }

    /**
     * Fault-system rupture index and magnitude for one rupture, as
     * recorded in a RuptureResultsExporter CSV.
     *
     * @param fltSysRupIndex nullable -- see {@link #readRuptureAttributes}
     */
    public record RuptureAttributes(Integer fltSysRupIndex, double magnitude) {}

    /**
     * Parses sourceId, ruptureId, fltSysRupIndex, magnitude from the
     * exported CSV's data rows, same scope/layout assumptions as
     * {@link #read}.
     *
     * <p>fltSysRupIndex is nullable per rupture (e.g. gridded/background
     * sources have no fault-system rupture index, and if the export
     * wasn't built from a fault-system solution at all, EVERY row's
     * value is empty) -- see
     * {@code RuptureLossResultsCollection.getFltSysRupIndex}. The column
     * header itself is always present (written unconditionally by
     * RuptureResultsExporter); only individual values may be empty.
     * magnitude, by contrast, is never empty -- every rupture has one
     * regardless of ERF type.</p>
     *
     * @throws IllegalArgumentException if the expected data header row
     *         ("sourceId,ruptureId,...") is never found, or if a required
     *         column is missing from it
     */
    public static Map<RuptureKey, RuptureAttributes> readRuptureAttributes(Path exportedCsv) throws IOException {

        Map<RuptureKey, RuptureAttributes> lookup = new HashMap<>();

        try (BufferedReader br = Files.newBufferedReader(exportedCsv)) {

            String line;
            Map<String, Integer> colIndex = null;

            while ((line = br.readLine()) != null) {
                if (line.trim().startsWith("sourceId")) {
                    colIndex = parseHeader(line);
                    break;
                }
            }

            if (colIndex == null) {
                throw new IllegalArgumentException(
                        "Could not find 'sourceId,ruptureId,...' data header row in: " + exportedCsv);
            }

            int sourceIdCol = requireColumn(colIndex, "sourceId", exportedCsv);
            int ruptureIdCol = requireColumn(colIndex, "ruptureId", exportedCsv);
            int fltSysRupIndexCol = requireColumn(colIndex, "fltSysRupIndex", exportedCsv);
            int magnitudeCol = requireColumn(colIndex, "magnitude", exportedCsv);
            int minColumns = Math.max(Math.max(sourceIdCol, ruptureIdCol), Math.max(fltSysRupIndexCol, magnitudeCol)) + 1;

            while ((line = br.readLine()) != null) {
                if (line.isBlank())
                    continue;

                String[] tokens = line.split(",");
                if (tokens.length < minColumns) {
                    throw new IllegalArgumentException(
                            "Malformed data row (expected at least " + minColumns
                                    + " columns): '" + line + "' in " + exportedCsv);
                }

                int sourceId = Integer.parseInt(tokens[sourceIdCol].trim());
                int ruptureId = Integer.parseInt(tokens[ruptureIdCol].trim());
                String fltSysRupIndexStr = tokens[fltSysRupIndexCol].trim();
                Integer fltSysRupIndex = fltSysRupIndexStr.isEmpty() ? null : Integer.parseInt(fltSysRupIndexStr);
                double magnitude = Double.parseDouble(tokens[magnitudeCol].trim());

                lookup.put(new RuptureKey(sourceId, ruptureId), new RuptureAttributes(fltSysRupIndex, magnitude));
            }
        }

        return lookup;
    }

    /**
     * Builds a column-name -> index map from the data header row.
     */
    private static Map<String, Integer> parseHeader(String headerLine) {
        String[] headers = headerLine.split(",");
        Map<String, Integer> colIndex = new HashMap<>();
        for (int i = 0; i < headers.length; i++) {
            colIndex.put(headers[i].trim(), i);
        }
        return colIndex;
    }

    /**
     * Looks up a required column's index by name, failing clearly if the
     * header row doesn't contain it.
     */
    private static int requireColumn(Map<String, Integer> colIndex, String name, Path exportedCsv) {
        Integer index = colIndex.get(name);
        if (index == null) {
            throw new IllegalArgumentException(
                    "Data header row is missing required column '" + name + "' in: " + exportedCsv);
        }
        return index;
    }
}