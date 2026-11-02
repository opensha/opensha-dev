package scratch.anne.sra.domain.hazard;

import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import scratch.anne.sra.domain.hazard.HazardParameters.HazardMetricType;
import scratch.anne.sra.io.writers.HazardCurvesExporter;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.SiteKey;

public final class HazardCurveCollection {

    // ------------------------------------------------------------
    // METADATA (single source of truth)
    // ------------------------------------------------------------

    private final HazardParameters hazardParameters;

    // ------------------------------------------------------------
    // DATA
    // ------------------------------------------------------------

    private final Map<SiteKey, Map<ImKey, HazardCurve>> data;

    private boolean frozen = false;

    // ------------------------------------------------------------
    // CONSTRUCTOR
    // ------------------------------------------------------------

    public HazardCurveCollection(HazardParameters hazardParameters) {
        this.hazardParameters = Objects.requireNonNull(hazardParameters);

        if (hazardParameters.getHazardMetricType() != HazardMetricType.EXCEEDANCE) {
            throw new IllegalArgumentException(
                "HazardCurveCollection requires EXCEEDANCE hazard metric type, but got: "
                + hazardParameters.getHazardMetricType()
            );
        }

        this.data = new ConcurrentHashMap<>();
    }

    // ------------------------------------------------------------
    // MUTATION PHASE
    // ------------------------------------------------------------

    public void put(SiteKey site, ImKey im, HazardCurve hazard) {
        if (frozen) {
            throw new IllegalStateException("HazardCurveCollection is frozen");
        }

        data.computeIfAbsent(site, s -> new ConcurrentHashMap<>())
        .put(im, hazard);
    }

    // ------------------------------------------------------------
    // FREEZE
    // ------------------------------------------------------------

    public void freeze() {
        if (frozen) return;

        for (Map.Entry<SiteKey, Map<ImKey, HazardCurve>> e : data.entrySet()) {
            e.setValue(Collections.unmodifiableMap(e.getValue()));
        }

        frozen = true;
    }

    public boolean isFrozen() {
        return frozen;
    }

    // ------------------------------------------------------------
    // ACCESSORS
    // ------------------------------------------------------------

    public HazardParameters getParameters() {
        return hazardParameters;
    }

    public Map<SiteKey, Map<ImKey, HazardCurve>> asMap() {
        return Collections.unmodifiableMap(data);
    }

    public Set<SiteKey> getSiteKeys() {
        return data.keySet();
    }

    public Map<ImKey, HazardCurve> getSiteView(SiteKey site) {
        return data.getOrDefault(site, Map.of());
    }

    public HazardCurve get(SiteKey site, ImKey im) {
        return data.getOrDefault(site, Map.of()).get(im);
    }

    // ------------------------------------------------------------
    // ITERATION (important for export)
    // ------------------------------------------------------------

    public void forEachCurve(HazardConsumer consumer) {
        data.forEach((site, imMap) ->
            imMap.forEach((im, curve) ->
                consumer.accept(site, im, curve)
            )
        );
    }

    @FunctionalInterface
    public interface HazardConsumer {
        void accept(SiteKey site, ImKey im, HazardCurve curve);
    }

    // ------------------------------------------------------------
    // EXPORT HOOK (key addition)
    // ------------------------------------------------------------

    public void export(HazardCurvesExporter exporter, Path file) {
        exporter.export(this, file);
    }
}