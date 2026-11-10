package scratch.anne.sra.util.enums;

/**
 * The kind of hazard quantity used by a risk calculation: rate or probability, considering
 * exceedance or occurrence.
 */
public enum HazardMetric {
    RATE_EXCEEDANCE(HazardMetric.HazardMetricType.EXCEEDANCE),
    PROBABILITY_EXCEEDANCE(HazardMetric.HazardMetricType.EXCEEDANCE),
    RATE_OCCURENCE(HazardMetric.HazardMetricType.OCCURRENCE),
    PROBABILITY_OCCURENCE(HazardMetric.HazardMetricType.OCCURRENCE);

    private final HazardMetricType metricType;

    HazardMetric(HazardMetricType metricType) {
        this.metricType = metricType;
    }

    public HazardMetricType getHazardMetricType() {
        return metricType;
    }

    /**
     * Whether a {@link HazardMetric} describes an exceedance quantity
     * (rate/probability of exceeding a level) or an occurrence quantity
     * (rate/probability of occurrence).
     */
    public enum HazardMetricType {
        EXCEEDANCE,
        OCCURRENCE
    }
}
