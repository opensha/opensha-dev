package scratch.anne.sra.io.writers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.avro.Schema;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.parquet.avro.AvroParquetWriter;
import org.apache.parquet.hadoop.ParquetFileWriter;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.io.LocalOutputFile;
import org.apache.parquet.io.OutputFile;

import scratch.anne.sra.domain.hazard.LossExceedanceCurve;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.io.readers.RuptureLossLookupReader.RuptureAttributes;

/**
 * Writes the FULL discrete-rank pseudo loss curve for one or more
 * branches.
 *
 * <p>Two independent output formats, same underlying per-branch,
 * per-rank data ({@code loss}, {@code cumulativeLikelihood}, and the
 * specific rupture responsible):</p>
 * <ul>
 *   <li>{@link #writeCsv} -- human-readable CSV:
 *       "loss, sourceId, ruptureId, fltSysRupIndex, magnitude, branch_0,
 *       branch_1, ..." with one row per distinct loss value across the
 *       included branches and one column per branch. Each branch column
 *       is that branch's cumulative likelihood carried forward
 *       (step-function-wise, via
 *       {@link LossExceedanceCurve#cumulativeLikelihoodAtOrBelow}) to a
 *       loss value that may not be one of that branch's own rank
 *       points -- but since loss values come from the fixed true-mean
 *       per-rupture loss lookup, a given loss row corresponds to one
 *       specific true-mean rupture (modulo coincidental loss-value
 *       collisions between distinct ruptures), so sourceId/ruptureId/
 *       fltSysRupIndex/magnitude are attributed once per row, not per
 *       branch column. Carries exactly the same cumulative-likelihood
 *       information as the Parquet output below for the branches
 *       included (every branch's own (loss, cumulativeLikelihood) pair
 *       appears somewhere in this table). Only meant for a handful of
 *       branches at once -- the union loss axis and every included
 *       branch's curve are built and held in memory together to produce
 *       it, unlike the streamed Parquet output.</li>
 *   <li>{@link #writeParquet} -- long format, one row per
 *       (branchIndex, rank), with the specific rupture (sourceId,
 *       ruptureId, fltSysRupIndex, magnitude) responsible for that
 *       rank's loss threshold. Columnar, compressed, scales to the full
 *       logic tree -- same streaming/DuckDB pattern as
 *       {@link LogicTreeBranchRuptureContributionExporter}.</li>
 * </ul>
 *
 * <p>Neither method's {@link LossExceedanceCurve} input carries
 * fltSysRupIndex or magnitude directly (neither is part of a RuptureKey,
 * and nothing upstream -- {@code LecAccumulation.discreteRank}, the
 * contributions Parquet schema, {@code StreamingBranchLecBuilder} --
 * threads them through). Both methods instead take a separate
 * {@code Map<RuptureKey, RuptureAttributes>} lookup (see
 * {@code RuptureLossLookupReader#readRuptureAttributes}, built from the
 * same fixed true-mean per-rupture export the loss values themselves
 * come from) and resolve it per rupture key. A rupture absent from the
 * lookup (gridded/background sources, or an export not built from a
 * fault-system solution at all) writes a null/empty fltSysRupIndex and
 * magnitude (magnitude is otherwise never empty for a rupture present in
 * the lookup).</p>
 *
 * <p>Usage: open once per export run (try-with-resources across the
 * whole branch loop). Call {@link #writeParquet} once per branch as
 * each branch's {@link LossExceedanceCurve} becomes available --
 * writing per branch as curves stream in (rather than collecting every
 * curve first) keeps memory flat regardless of how many branches are
 * exported. {@link #writeCsv} instead takes every targeted branch's
 * curve at once, since it needs the full set to compute the shared
 * union loss axis. Close once at the end.</p>
 */
public final class PseudoLossCurveExporter implements AutoCloseable {

    private static final Schema SCHEMA = new Schema.Parser().parse("""
            {
              "type": "record",
              "name": "PseudoLossCurveRow",
              "namespace": "scratch.anne.risk_system_vb",
              "fields": [
                {"name": "branchIndex",         "type": "int"},
                {"name": "rank",                "type": "int"},
                {"name": "loss",                "type": "double"},
                {"name": "cumulativeLikelihood", "type": "double"},
                {"name": "sourceId",            "type": ["null", "int"], "default": null},
                {"name": "fltSysRupIndex",      "type": ["null", "int"], "default": null},
                {"name": "ruptureId",           "type": ["null", "int"], "default": null},
                {"name": "magnitude",           "type": ["null", "double"], "default": null}
              ]
            }
            """);

    private final ParquetWriter<GenericRecord> parquetWriter;

    /**
     * Opens a Parquet-only exporter, with no run-level metadata.
     * {@link #writeCsv} doesn't depend on anything opened by this
     * constructor and can still be called on an instance created this
     * way.
     */
    public PseudoLossCurveExporter(Path parquetFile) throws IOException {
        this(parquetFile, java.util.Collections.emptyMap());
    }

    /**
     * @param extraMetadata run-level provenance written into the
     *        Parquet file's own footer key-value metadata (see
     *        {@code LogicTreeBranchRuptureContributionExporter}'s
     *        javadoc for the full rationale -- typically propagated
     *        forward from the branch contributions Parquet this
     *        curve was streamed from, via
     *        {@code ParquetMetadataReader}, rather than independently
     *        reconstructed). Every key and value must be a
     *        {@code String} (Parquet's own requirement); pass an
     *        empty map for none.
     */
    public PseudoLossCurveExporter(Path parquetFile, Map<String, String> extraMetadata) throws IOException {
        if (parquetFile != null) {
            OutputFile out = new LocalOutputFile(parquetFile.toAbsolutePath().normalize());
            this.parquetWriter = AvroParquetWriter
                    .<GenericRecord>builder(out)
                    .withSchema(SCHEMA)
                    .withCompressionCodec(CompressionCodecName.SNAPPY)
                    .withWriteMode(ParquetFileWriter.Mode.OVERWRITE)
                    .withExtraMetaData(extraMetadata)
                    .build();
        } else {
            this.parquetWriter = null;
        }
    }

    /**
     * Writes every rank of one branch's curve to the Parquet output.
     *
     * @param attributesByKey fltSysRupIndex/magnitude per RuptureKey
     *        (see {@code RuptureLossLookupReader#readRuptureAttributes});
     *        a rupture absent from this map writes null fltSysRupIndex
     *        and magnitude
     */
    public void writeParquet(int branchIndex, LossExceedanceCurve curve, Map<RuptureKey, RuptureAttributes> attributesByKey) {
        if (parquetWriter == null)
            throw new IllegalStateException("No Parquet output file was configured for this exporter");
        try {
            for (int i = 0; i < curve.size(); i++) {
                RuptureKey key = curve.keyAt(i);
                RuptureAttributes attrs = key == null ? null : attributesByKey.get(key);
                GenericRecord record = new GenericData.Record(SCHEMA);
                record.put("branchIndex", branchIndex);
                record.put("rank", i);
                record.put("loss", curve.lossAt(i));
                record.put("cumulativeLikelihood", curve.cumulativeLikelihoodAt(i));
                record.put("sourceId", key == null ? null : key.sourceId());
                record.put("fltSysRupIndex", attrs == null ? null : attrs.fltSysRupIndex());
                record.put("ruptureId", key == null ? null : key.ruptureId());                
                record.put("magnitude", attrs == null ? null : attrs.magnitude());
                parquetWriter.write(record);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Writes a CSV for the given branches: "loss", "sourceId",
     * "fltSysRupIndex", "ruptureId", "magnitude" columns followed by one
     * "branch_&lt;index&gt;" column per branch (ascending branch-index
     * order), with one row per distinct loss value across every included
     * branch's curve.
     *
     * <p>A branch's cell is {@code NaN} when that row's own rupture
     * (sourceId/ruptureId) never occurs in that branch at all (zero rate
     * / not mapped), as opposed to a legitimate 0 cumulative likelihood
     * (every rupture that DOES occur in that branch has a larger loss
     * than this row's). Both would otherwise read as the same step-
     * function value via {@link LossExceedanceCurve#cumulativeLikelihoodAtOrBelow}
     * -- this check distinguishes them by looking for the row's exact
     * RuptureKey among that branch's own rank points, not just by the
     * resulting value.</p>
     *
     * @param curvesByBranch every branch to include, keyed by branch index
     * @param attributesByKey fltSysRupIndex/magnitude per RuptureKey
     *        (see {@code RuptureLossLookupReader#readRuptureAttributes});
     *        a rupture absent from this map writes empty fltSysRupIndex
     *        and magnitude cells
     * @param metadataComment optional single-line {@code # key=val,...}
     *        provenance comment written immediately before the header
     *        row (e.g. propagated from the branch contributions
     *        Parquet's own file-level metadata via
     *        {@code ParquetMetadataReader}) -- see
     *        {@code LogicTreeBranchEalSummaryExporter#write}'s javadoc
     *        for the full rationale. Pass {@code null} or blank to omit it.
     */
    public void writeCsv(Map<Integer, LossExceedanceCurve> curvesByBranch, Map<RuptureKey, RuptureAttributes> attributesByKey,
            String metadataComment, Path outputCsv) throws IOException {
        List<Integer> branchOrder = new ArrayList<>(curvesByBranch.keySet());
        branchOrder.sort(null);

        // one distinct loss value -> whichever rupture key produced it,
        // for the sourceId/fltSysRupIndex/ruptureId/magnitude columns.
        // Losses come from the fixed true-mean per-rupture lookup, so
        // (modulo coincidental loss-value collisions between distinct
        // ruptures) every branch sharing a loss value is referring to
        // the same rupture -- the first key encountered for a given loss
        // is used.
        Map<Double, RuptureKey> keyByLoss = new HashMap<>();
        TreeSet<Double> allLosses = new TreeSet<>();
        for (LossExceedanceCurve curve : curvesByBranch.values()) {
            for (int i = 0; i < curve.size(); i++) {
                double loss = curve.lossAt(i);
                allLosses.add(loss);
                keyByLoss.putIfAbsent(loss, curve.keyAt(i));
            }
        }
        double[] lossAxis = allLosses.stream().mapToDouble(Double::doubleValue).toArray();

        // per branch, the set of RuptureKeys that actually occur in that
        // branch's own curve -- used below to tell "this rupture is
        // absent from this branch" (NaN) apart from "this rupture is
        // present in this branch but every one of its own ruptures has a
        // loss <= this row's loss, carrying a possibly-0 cumulative
        // likelihood forward" (a legitimate value, not a missing one)
        Map<Integer, Set<RuptureKey>> keysPresentByBranch = new HashMap<>();
        for (Map.Entry<Integer, LossExceedanceCurve> entry : curvesByBranch.entrySet()) {
            LossExceedanceCurve curve = entry.getValue();
            Set<RuptureKey> keys = new HashSet<>();
            for (int i = 0; i < curve.size(); i++) {
                RuptureKey key = curve.keyAt(i);
                if (key != null)
                    keys.add(key);
            }
            keysPresentByBranch.put(entry.getKey(), keys);
        }

        try (BufferedWriter w = Files.newBufferedWriter(outputCsv)) {
            if (metadataComment != null && !metadataComment.isBlank()) {
                w.write("# " + metadataComment);
                w.newLine();
            }

            StringBuilder header = new StringBuilder("loss,sourceId,fltSysRupIndex,ruptureId,magnitude");
            for (int b : branchOrder)
                header.append(",branch_").append(b);
            w.write(header.toString());
            w.newLine();

            for (double loss : lossAxis) {
                RuptureKey key = keyByLoss.get(loss);
                RuptureAttributes attrs = key == null ? null : attributesByKey.get(key);

                StringBuilder row = new StringBuilder(String.valueOf(loss));
                row.append(',').append(key == null ? "" : key.sourceId());
                row.append(',').append(attrs == null || attrs.fltSysRupIndex() == null ? "" : attrs.fltSysRupIndex());
                row.append(',').append(key == null ? "" : key.ruptureId());
                row.append(',').append(attrs == null ? "" : attrs.magnitude());
                for (int b : branchOrder) {
                    boolean ruptureAbsentFromBranch = key != null
                            && !keysPresentByBranch.get(b).contains(key);
                    double value = ruptureAbsentFromBranch
                            ? Double.NaN
                            : curvesByBranch.get(b).cumulativeLikelihoodAtOrBelow(loss);
                    row.append(',').append(value);
                }
                w.write(row.toString());
                w.newLine();
            }
        }
    }

    @Override
    public void close() throws IOException {
        if (parquetWriter != null)
            parquetWriter.close();
    }
}
