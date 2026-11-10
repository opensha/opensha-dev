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

import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.HazardParametersMetadataBuilder;
import scratch.anne.sra.engine.perRupture.RuptureLossResultsCollection;
import scratch.anne.sra.engine.perRupture.RuptureLossResultsCollection.RuptureResult;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

/**
 * Writes a {@link RuptureLossResultsCollection} to CSV or Parquet, one
 * row per rupture (sourceId, fltSysRupIndex, ruptureId, magnitude,
 * likelihood, ruptureLoss, lossContribution, relativeContribution).
 *
 * <p>Both formats carry the same file-level provenance metadata: the
 * standard {@code HazardParameters} fields (via
 * {@link HazardParametersMetadataBuilder}) plus {@code totalExpectedLoss},
 * merged with whatever {@code extraMetadata} the caller supplies (e.g.
 * portfolio/response-model provenance -- see
 * {@link scratch.anne.sra.domain.portfolio.PortfolioMetadataBuilder}/
 * {@link scratch.anne.sra.domain.structuralResponse.ResponseModelMetadataBuilder}).
 * For CSV this is a leading {@code # HazardParameters} key,value block;
 * for Parquet it's the file's own footer key-value metadata.</p>
 */
public final class RuptureResultsExporter {

    private static final Schema SCHEMA = new Schema.Parser().parse("""
            {
              "type": "record",
              "name": "RuptureResult",
              "namespace": "scratch.anne.sra",
              "fields": [
                {"name": "sourceId",             "type": "int"},
                {"name": "fltSysRupIndex",       "type": ["null", "int"], "default": null},
                {"name": "ruptureId",            "type": "int"},
                {"name": "magnitude",            "type": "double"},
                {"name": "likelihood",           "type": "double"},
                {"name": "ruptureLoss",          "type": "double"},
                {"name": "lossContribution",     "type": "double"},
                {"name": "relativeContribution", "type": "double"}
              ]
            }
            """);

    // ------------------------------------------------------------
    // PUBLIC ENTRY
    // ------------------------------------------------------------
    public void export(RuptureLossResultsCollection collection,
                       Path file,
                       FileFormat format) {
        export(collection, file, format, Map.of());
    }

    /**
     * @param extraMetadata additional {@code key -> value} reproducibility
     *                      metadata to merge in alongside the standard
     *                      HazardParameters fields -- e.g.
     *                      {@link scratch.anne.sra.domain.portfolio.PortfolioMetadataBuilder}/
     *                      {@link scratch.anne.sra.domain.structuralResponse.ResponseModelMetadataBuilder}'s
     *                      output. Pass {@code Map.of()} if this
     *                      exporter's caller has no portfolio/response
     *                      model file to attribute (e.g. {@code ExploreERF},
     *                      which has no portfolio at all).
     */
    public void export(RuptureLossResultsCollection collection,
                       Path file,
                       FileFormat format,
                       Map<String, String> extraMetadata) {
        switch (format) {
            case CSV     -> exportCsv(collection, file, extraMetadata);
            case PARQUET -> exportParquet(collection, file, extraMetadata);
            default      -> throw new IllegalArgumentException("Unknown format: " + format);
        }
    }

    // ------------------------------------------------------------
    // Shared: builds the full reproducibility metadata map
    // ------------------------------------------------------------
    private Map<String, String> buildMetadata(HazardParameters hp, double totalExpectedLoss,
                                                Map<String, String> extraMetadata) {
        // delegates the standard erfName/erfDuration/hazardMetric/gmmName/
        // solutionFile/erfConfig.*/filterConfig.* fields to the shared
        // builder (also used by LogicTreeBranchRuptureContributionExporter,
        // so both files record IDENTICAL HazardParameters provenance);
        // totalExpectedLoss is specific to this exporter, added on top,
        // followed by whatever portfolio/response-model provenance the
        // caller supplied (empty map if none -- e.g. ExploreERF, which
        // has no portfolio at all).
        Map<String, String> meta = new java.util.HashMap<>(HazardParametersMetadataBuilder.build(hp));
        meta.put("totalExpectedLoss", String.valueOf(totalExpectedLoss));
        meta.putAll(extraMetadata);
        return meta;
    }

    // ------------------------------------------------------------
    // PARQUET EXPORT
    // ------------------------------------------------------------
    private void exportParquet(RuptureLossResultsCollection collection, Path file,
                                Map<String, String> extraMetadata) {
        HazardParameters hp = collection.getHazardParameters();

        Map<String, String> meta = buildMetadata(hp, collection.getTotalExpectedLoss(), extraMetadata);

        OutputFile out = new LocalOutputFile(file.toAbsolutePath().normalize());

        try (ParquetWriter<GenericRecord> writer = AvroParquetWriter
                .<GenericRecord>builder(out)
                .withSchema(SCHEMA)
                .withCompressionCodec(CompressionCodecName.SNAPPY)
                .withWriteMode(ParquetFileWriter.Mode.OVERWRITE)
                .withExtraMetaData(meta)
                .build()) {

            for (RuptureResult r : collection.results()) {
                GenericRecord record = new GenericData.Record(SCHEMA);
                record.put("sourceId",            r.key().sourceId());
                record.put("fltSysRupIndex",      r.fltSysRupIndex());
                record.put("ruptureId",           r.key().ruptureId());
                record.put("magnitude",           r.magnitude());
                record.put("likelihood",          r.likelihood());
                record.put("ruptureLoss",         r.ruptureLoss());
                record.put("lossContribution",    r.lossContribution());
                record.put("relativeContribution", r.relativeContribution());
                writer.write(record);
            }

        } catch (IOException e) {
            throw new RuntimeException("Parquet export failed", e);
        }
        System.out.println("Parquet export complete: " + file.toAbsolutePath());
    }

    // ------------------------------------------------------------
    // CSV EXPORT
    // ------------------------------------------------------------
    private void exportCsv(RuptureLossResultsCollection collection, Path file,
                            Map<String, String> extraMetadata) {
        HazardParameters hazardParameters = collection.getHazardParameters();

        Map<String, String> meta = buildMetadata(hazardParameters, collection.getTotalExpectedLoss(), extraMetadata);
        // totalExpectedLoss is written separately below, in its own
        // blank-line-delimited block, matching the original file's
        // layout -- pull it out here so it isn't double-written above.
        String totalExpectedLossStr = meta.remove("totalExpectedLoss");

        try (BufferedWriter writer = Files.newBufferedWriter(file)) {
            writer.write("# HazardParameters");
            writer.newLine();

            // Fixed fields first, in a stable order.
            writeParam(writer, "erfName",      meta.remove("erfName"));
            writeParam(writer, "erfDuration",  meta.remove("erfDuration"));
            writeParam(writer, "hazardMetric", meta.remove("hazardMetric"));
            writeParam(writer, "gmmName",      meta.remove("gmmName"));
            writeParam(writer, "solutionFile", meta.remove("solutionFile"));

            // Remaining entries: erfConfig.*, filterConfig.*, and any
            // extraMetadata the caller supplied (e.g. portfolio.*/
            // responseModel.* -- see PortfolioMetadataBuilder/
            // ResponseModelMetadataBuilder); order not guaranteed, same
            // as the original filterConfig loop.
            for (Map.Entry<String, String> e : meta.entrySet()) {
                writeParam(writer, e.getKey(), e.getValue());
            }

            writer.newLine();
            writeParam(writer, "totalExpectedLoss", totalExpectedLossStr);
            writer.newLine();

            writer.write("sourceId,fltSysRupIndex,ruptureId,magnitude,likelihood,ruptureLoss,lossContribution,relativeContribution");
            writer.newLine();

            for (RuptureResult r : collection.results()) {
                writer.write(formatCsvRow(r));
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException("CSV export failed", e);
        }
    }

    private String formatCsvRow(RuptureResult r) {
        return r.key().sourceId()
        		+ "," + (r.fltSysRupIndex() == null ? "" : r.fltSysRupIndex())
                + "," + r.key().ruptureId()                
                + "," + r.magnitude()
                + "," + r.likelihood()
                + "," + r.ruptureLoss()
                + "," + r.lossContribution()
                + "," + r.relativeContribution();
    }

    private void writeParam(BufferedWriter writer, String param, String value) {
        CsvParamWriter.writeParam(writer, param, value);
    }
}