package scratch.anne.sra.io.writers;

/** Output format selector for exporters that support more than one, e.g. {@link RuptureResultsExporter}. */
public enum FileFormat {
    CSV,
    PARQUET
}