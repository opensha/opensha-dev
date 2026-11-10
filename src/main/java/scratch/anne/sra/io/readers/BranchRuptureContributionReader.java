package scratch.anne.sra.io.readers;

import org.apache.avro.generic.GenericRecord;
import org.apache.parquet.avro.AvroParquetReader;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.io.LocalInputFile;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads back the per-branch, per-rupture (loss, rate) contribution
 * Parquet file written by {@code LogicTreeBranchRuptureContributionExporter}
 * -- the one frozen artifact every downstream LEC accumulation strategy
 * (discrete-rank curves, lognormal-COV smearing, etc.) replays against.
 *
 * <p>{@link #readAll} returns full rows with rupture identity intact;
 * {@link #readGroupedByBranch} drops identity and groups by branch,
 * the shape {@code LecAccumulationStrategy} actually needs.</p>
 */
public final class BranchRuptureContributionReader {

    private BranchRuptureContributionReader() {}

    /** Full rows, identity preserved -- for diagnostics, dedupe checks, tracing a specific rupture. */
    public static List<BranchRuptureContributionRow> readAll(Path parquetFile) throws IOException {
        List<BranchRuptureContributionRow> rows = new ArrayList<>();
        try (ParquetReader<GenericRecord> reader = AvroParquetReader
                .<GenericRecord>builder(new LocalInputFile(parquetFile.toAbsolutePath().normalize()))
                .build()) {

            GenericRecord record;
            while ((record = reader.read()) != null) {
                rows.add(new BranchRuptureContributionRow(
                        (int) record.get("branchIndex"),
                        (int) record.get("sourceId"),
                        (int) record.get("ruptureId"),
                        (double) record.get("loss"),
                        (double) record.get("rate")));
            }
        }
        return rows;
    }

    /** Grouped by branch, identity dropped -- the shape LecAccumulationStrategy actually needs. */
    public static Map<Integer, List<LossRateContribution>> readGroupedByBranch(Path parquetFile) throws IOException {
        Map<Integer, List<LossRateContribution>> byBranch = new HashMap<>();
        for (BranchRuptureContributionRow row : readAll(parquetFile)) {
            byBranch.computeIfAbsent(row.branchIndex(), k -> new ArrayList<>())
                    .add(new LossRateContribution(row.loss(), row.rate()));
        }
        return byBranch;
    }
    
    /** Full row as stored in the contributions Parquet -- includes rupture identity. */
    public record BranchRuptureContributionRow(int branchIndex, int sourceId, int ruptureId, double loss, double rate) {}
    
    /** One rupture's (loss, rate) contribution to a branch's LEC. */
    public record LossRateContribution(double loss, double rate) {}
}
