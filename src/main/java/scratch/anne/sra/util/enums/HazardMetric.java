package scratch.anne.sra.util.enums;

/**
 * The kind of hazard quantity used by a run: rate or probability, of
 * exceedance or occurrence.
 */
public enum HazardMetric {
    RATE_EXCEEDANCE(HazardMetricType.EXCEEDANCE),
    PROBABILITY_EXCEEDANCE(HazardMetricType.EXCEEDANCE),
    RATE_OCCURENCE(HazardMetricType.OCCURRENCE),
    PROBABILITY_OCCURENCE(HazardMetricType.OCCURRENCE);

    private final HazardMetricType metricType;

    HazardMetric(HazardMetricType metricType) {
        this.metricType = metricType;
    }

    public HazardMetricType getHazardMetricType() {
        return metricType;
    }
}
