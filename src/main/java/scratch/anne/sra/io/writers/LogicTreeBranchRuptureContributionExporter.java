package scratch.anne.sra.io.writers;

import org.apache.avro.Schema;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.parquet.avro.AvroParquetWriter;
import org.apache.parquet.hadoop.ParquetFileWriter;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.io.LocalOutputFile;
import org.apache.parquet.io.OutputFile;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Map;

/**
 * Writes the raw per-branch, per-rupture (loss, rate) contribution file --
 * NOT an accumulated/ranked curve of any kind. One row per rupture that
 * had nonzero rate on a given branch, in whatever order the branch loop
 * naturally produces them (fault ruptures then gridded ruptures, per the
 * existing mapping loops).
 *
 * <p>This is the single frozen artifact that every downstream LEC
 * question -- discrete-rank curves, lognormal-COV smearing, fixed
 * exceedance-probability levels, or any future accumulation strategy --
 * replays against. None of that logic runs here; this class only writes
 * rows.</p>
 *
 * <p>Usage: open once per export run, call {@link #write} for every
 * nonzero-rate rupture contribution as the branch loop computes them,
 * close once at the very end (try-with-resources across the whole
 * branch loop, not per-branch).</p>
 *
 * <p>Since this is the one artifact every downstream consumer (
 * {@code StreamingBranchLecBuilder}, {@code PseudoLossCurveStage},
 * {@code LecSummaryStage}, any future accumulation strategy) reads, it
 * is also the natural place to record run-level provenance -- which
 * solution logic tree zip these branches came from (which may itself
 * already be a downsampled subset of a much larger original tree), how
 * many of its branches were actually run, the full
 * {@code HazardParameters} that produced the true-mean ERF these branch
 * ruptures are mapped against (see
 * {@code HazardParametersMetadataBuilder}, the same field set
 * {@code RuptureResultsExporter} already writes into
 * {@code _eal-per-rup.csv}'s own metadata block), AND that same CSV's
 * portfolio/response-model provenance (see
 * {@code PortfolioMetadataBuilder}/{@code ResponseModelMetadataBuilder}
 * -- which portfolio CSV and vulnerability library JSON these losses
 * were convolved from, plus the portfolio's own total asset value)
 * -- recording the IDENTICAL provenance here too, rather than only the
 * branch/zip fields, means a reader of this Parquet alone can tell
 * which ERF, portfolio, and response-model library the fixed
 * per-rupture losses it joins against actually came from -- via
 * Parquet's own file-level key-value metadata (see
 * {@link #LogicTreeBranchRuptureContributionExporter(Path, Map)}),
 * rather than needing a separate CSV comment line that every consumer
 * would have to know to read back from a DIFFERENT file. This is read
 * back the same way any Parquet reader exposes file metadata (e.g.
 * pyarrow's {@code ParquetFile(path).metadata.metadata}, or DuckDB's
 * {@code parquet_metadata()} / {@code parquet_kv_metadata()}).</p>
 */
public final class LogicTreeBranchRuptureContributionExporter implements AutoCloseable {

    private static final Schema SCHEMA = new Schema.Parser().parse("""
            {
              "type": "record",
              "name": "BranchRuptureContribution",
              "namespace": "scratch.anne.sra",
              "fields": [
                {"name": "branchIndex", "type": "int"},
                {"name": "sourceId",    "type": "int"},
                {"name": "ruptureId",   "type": "int"},
                {"name": "loss",        "type": "double"},
                {"name": "rate",        "type": "double"}
              ]
            }
            """);

    private final ParquetWriter<GenericRecord> writer;

    /** No run-level metadata. */
    public LogicTreeBranchRuptureContributionExporter(Path file) throws IOException {
        this(file, Collections.emptyMap());
    }

    /**
     * @param extraMetadata run-level provenance written into the
     *        Parquet file's own footer key-value metadata (e.g.
     *        {@code numBranches}, {@code solutionLogicTreeZip},
     *        {@code fileTag} -- see {@code PipelineConfig}). Every key
     *        and value must be a {@code String} (Parquet's own
     *        requirement); pass an empty map for none.
     */
    public LogicTreeBranchRuptureContributionExporter(Path file, Map<String, String> extraMetadata) throws IOException {
        OutputFile out = new LocalOutputFile(file.toAbsolutePath().normalize());
        this.writer = AvroParquetWriter
                .<GenericRecord>builder(out)
                .withSchema(SCHEMA)
                .withCompressionCodec(CompressionCodecName.SNAPPY)
                .withWriteMode(ParquetFileWriter.Mode.OVERWRITE)
                .withExtraMetaData(extraMetadata)
                .build();
    }

    /**
     * Writes one rupture's contribution for one branch. Skips silently if
     * rate is zero (nothing for downstream accumulation to do with a
     * zero-rate row -- same "skip the lookup entirely" convention as
     * BranchEALCalculator.accumulate).
     */
    public void write(int branchIndex, int sourceId, int ruptureId, double loss, double rate) {
        if (rate == 0d)
            return;
        GenericRecord record = new GenericData.Record(SCHEMA);
        record.put("branchIndex", branchIndex);
        record.put("sourceId", sourceId);
        record.put("ruptureId", ruptureId);
        record.put("loss", loss);
        record.put("rate", rate);
        try {
            writer.write(record);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}