package scratch.anne.sra.io.writers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import scratch.anne.sra.domain.hazard.HazardCurve;
import scratch.anne.sra.domain.hazard.HazardParameters;
import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.engine.per_rupture.RuptureProbabilityCurveCollection;
import scratch.anne.sra.engine.per_rupture.RuptureProbabilityCurveCollection.RuptureExceedanceResult;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

/**
 * Writes a {@link RuptureProbabilityCurveCollection} to CSV.
 *
 * <p>Unlike {@code RuptureResultsExporter} (many assets, one loss value per
 * rupture) or {@code HazardCurvesExporter} (many sites/IMTs, one curve per
 * site+IMT), this collection is single-asset / single-site / single-IMT: one
 * {@link SiteKey} and one {@link ImKey} apply to every row, and both are
 * stored on the collection itself. Both are written once in the header,
 * alongside the {@link HazardParameters} metadata, rather than repeated per
 * row.</p>
 *
 * <p>Requires the collection to already be frozen (see
 * {@link RuptureProbabilityCurveCollection#freeze()}, now called
 * automatically at the end of {@code fromERF}), since the direct and
 * aggregate hazard curves are only available after freezing.</p>
 *
 * <p>Layout:
 * <pre>
 * # HazardParameters
 * erfName,...
 * erfDuration,...
 * hazardMetric,...
 * gmmName,...
 * [filter config entries]
 * &lt;blank&gt;
 * # Site
 * lat,...
 * lon,...
 * vs30,...
 * &lt;blank&gt;
 * # ImKey
 * imt,...
 * &lt;blank&gt;
 * ,,,,imls,v1,v2,v3,...
 * ,,,,direct hcurve,v1,v2,v3,...
 * &lt;blank&gt;
 * ,,,,aggregate hcurve,v1,v2,v3,...
 * &lt;blank&gt;
 * sourceId,ruptureId,likelihood,magnitude,distance
 * ...per-rupture rows, with curve values appended...
 * </pre>
 * The {@code imls}/{@code direct hcurve}/{@code aggregate hcurve} rows have
 * 4 leading empty fields so the label lands in the 5th column and values
 * line up under the curve columns of the per-rupture rows below (which have
 * 5 leading columns: sourceId, ruptureId, likelihood, magnitude, distance).</p>
 */
public final class RuptureProbabilityCurvesExporter {

    public void export(RuptureProbabilityCurveCollection collection, Path file) {
        HazardParameters hazardParameters = collection.getHazardParameters();
        Map<String, String> filterConfig = hazardParameters.getFilterConfig();
        SiteKey siteKey = collection.getSiteKey();
        ImKey imKey = collection.getImKey();

        double[] imls = imKey.getLinearValues();

        try (BufferedWriter writer = Files.newBufferedWriter(file)) {

            // ---------------- HEADER: hazard parameters ----------------
            writer.write("# HazardParameters");
            writer.newLine();
            writeParam(writer, "erfName", hazardParameters.getErfName());
            writeParam(writer, "erfDuration", String.valueOf(hazardParameters.getErfDuration()));
            writeParam(writer, "hazardMetric", hazardParameters.getHazardMetric().name());
            writeParam(writer, "gmmName", hazardParameters.getGmmName());
            if (filterConfig != null) {
                for (Map.Entry<String, String> e : filterConfig.entrySet()) {
                    writeParam(writer, e.getKey(), e.getValue());
                }
            }
            writer.newLine();

            // ---------------- HEADER: site ----------------
            writer.write("# Site");
            writer.newLine();
            writeParam(writer, "lat", String.valueOf(siteKey.getLat()));
            writeParam(writer, "lon", String.valueOf(siteKey.getLon()));
            writeParam(writer, "vs30", String.valueOf(siteKey.getVs30()));
            writer.newLine();

            // ---------------- HEADER: IM ----------------
            writer.write("# ImKey");
            writer.newLine();
            writeParam(writer, "imt", imKey.getImtString());
            writer.newLine();

            // ---------------- IMLS / DIRECT / AGGREGATE (indented) ----------------
            writer.write(formatIndentedRow("imls", imls));
            writer.newLine();
            writer.write(formatIndentedRow("direct hcurve", collection.getDirectHazardCurve().getHazard()));
            writer.newLine();
            writer.write(formatIndentedRow("aggregate hcurve", collection.getAggregateHazardCurve().getHazard()));
            writer.newLine();
            writer.newLine();

            // ---------------- DATA HEADER ROW ----------------
            writer.write("sourceId,ruptureId,likelihood,magnitude,distance");
            writer.newLine();

            // ---------------- DATA ROWS: PER-RUPTURE ----------------
            for (RuptureExceedanceResult r : collection.results()) {
                writer.write(formatRuptureRow(r));
                writer.newLine();
            }

        } catch (IOException e) {
            throw new RuntimeException("CSV export failed", e);
        }
    }

    /**
     * Formats a label + value array with 4 leading empty fields, so the
     * label lands in the 5th column and values line up under the curve
     * columns of the per-rupture rows (which have 5 leading columns before
     * curve values start).
     */
    private String formatIndentedRow(String label, double[] values) {
        StringBuilder row = new StringBuilder(",,,,");
        row.append(label);
        for (double v : values)
            row.append(",").append(v);
        return row.toString();
    }

    private String formatRuptureRow(RuptureExceedanceResult r) {
        RuptureKey key = r.key();
        HazardCurve curve = r.curve();

        StringBuilder row = new StringBuilder();
        row.append(key.sourceId())
                .append(",").append(key.ruptureId())
                .append(",").append(r.likelihood())
                .append(",").append(r.magnitude())
                .append(",").append(r.distance());

        for (double h : curve.getHazard())
            row.append(",").append(h);

        return row.toString();
    }

    private void writeParam(BufferedWriter writer, String param, String value) {
        CsvParamWriter.writeParam(writer, param, value);
    }
}