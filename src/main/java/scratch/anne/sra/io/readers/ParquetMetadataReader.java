package scratch.anne.sra.io.readers;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.parquet.hadoop.ParquetFileReader;
import org.apache.parquet.hadoop.metadata.ParquetMetadata;
import org.apache.parquet.io.LocalInputFile;

/**
 * Reads back the file-level key-value metadata written by
 * {@code LogicTreeBranchRuptureContributionExporter} (numBranches,
 * solutionLogicTreeZip, fileTag, plus the full HazardParameters
 * provenance fields -- see that class's and
 * {@code HazardParametersMetadataBuilder}'s javadoc), so downstream
 * pipeline stages (e.g. {@code PseudoLossCurveStage},
 * {@code LecSummaryStage}) can propagate the SAME provenance into their
 * own output files without needing to independently reconstruct it
 * (e.g. by reloading the SolutionLogicTree, which those stages are
 * otherwise careful to avoid).
 */
public final class ParquetMetadataReader {

    /**
     * Keys that {@code AvroParquetWriter} (via {@code AvroWriteSupport}/
     * {@code ParquetWriter} itself) injects into every file's footer
     * metadata automatically -- NOT user-level provenance. If one of
     * these is read back and then handed straight to a NEW
     * {@code AvroParquetWriter}'s own {@code withExtraMetaData}, that
     * writer tries to inject the same key itself, causing
     * {@code IllegalArgumentException: Duplicate metadata key}. Always
     * stripped by {@link #readKeyValueMetadata}, so callers never need
     * to know these keys exist at all before re-propagating a metadata
     * map into another Avro Parquet writer.
     *
     * <p>{@code parquet.avro.schema}: the serialized Avro schema itself
     * (written by {@code AvroWriteSupport.init}).</p>
     * <p>{@code writer.model.name}: which object model wrote the file
     * (always {@code "avro"} for {@code AvroParquetWriter}, written by
     * {@code ParquetWriter}'s own constructor via
     * {@code ParquetWriter.OBJECT_MODEL_NAME_PROP} -- confirmed from an
     * actual written file's footer, not assumed).</p>
     */
    private static final Set<String> RESERVED_AVRO_KEYS = Set.of(
            "parquet.avro.schema",
            "writer.model.name"
    );

    private ParquetMetadataReader() {}

    /**
     * @return the Parquet file's own footer key-value metadata, with
     *         any Avro-internal reserved keys (see
     *         {@link #RESERVED_AVRO_KEYS}) already stripped -- never
     *         {@code null}; empty if the file has no OTHER metadata
     */
    public static Map<String, String> readKeyValueMetadata(Path parquetFile) throws IOException {
        try (ParquetFileReader reader = ParquetFileReader.open(
                new LocalInputFile(parquetFile.toAbsolutePath().normalize()))) {
            ParquetMetadata footer = reader.getFooter();
            Map<String, String> raw = footer.getFileMetaData().getKeyValueMetaData();

            Map<String, String> filtered = new HashMap<>(raw);
            filtered.keySet().removeAll(RESERVED_AVRO_KEYS);
            return filtered;
        }
    }
}
