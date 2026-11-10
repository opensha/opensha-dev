package scratch.anne.sra.io.writers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;

/**
 * Shared helper for writing a single {@code key,value} metadata line,
 * used by the CSV exporters that write a leading {@code # ...} metadata
 * block (e.g. {@code RuptureResultsExporter}, {@code RuptureProbabilityCurvesExporter}).
 */
final class CsvParamWriter {

    private CsvParamWriter() {}

    /**
     * Writes {@code param,value} followed by a newline.
     *
     * @throws UncheckedIOException if the write fails
     */
    static void writeParam(BufferedWriter writer, String param, String value) {
        try {
            writer.write(param + "," + value);
            writer.newLine();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
