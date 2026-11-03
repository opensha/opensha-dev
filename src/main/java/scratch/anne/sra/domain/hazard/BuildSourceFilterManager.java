package scratch.anne.sra.domain.hazard;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.sha.calc.sourceFilters.FixedDistanceCutoffFilter;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.calc.sourceFilters.params.SourceFiltersParam;

/**
 * Builds a {@link SourceFilterManager} from a flat {@code key -> value}
 * filter config (e.g. {@code max_distance}, {@code min_mag}), as parsed
 * from a run-parameters CSV by
 * {@link scratch.anne.sra.io.readers.ParseRiskRunParametersCSV#extractSourceFilterConfig}.
 */
public class BuildSourceFilterManager {

    private final List<SourceFilters> enabledFilters = new ArrayList<>();

    private Double maxDistance;

    public BuildSourceFilterManager enable(SourceFilters filter) {
        enabledFilters.add(filter);
        return this;
    }

    public BuildSourceFilterManager maxDistance(double distance) {
        this.maxDistance = distance;
        enable(SourceFilters.FIXED_DIST_CUTOFF);
        return this;
    }

    public BuildSourceFilterManager fromConfig(Map<String,String> config) {

        String maxDistance = config.get("max_distance");

        if (maxDistance != null && !maxDistance.isBlank()) {
            maxDistance(Double.parseDouble(maxDistance));
        }

        return this;
    }

    public SourceFilterManager build() {

        SourceFilterManager manager =
                SourceFiltersParam.getDefault();

        if (maxDistance != null) {

            manager.setEnabled(
                    SourceFilters.FIXED_DIST_CUTOFF,
                    true);

            FixedDistanceCutoffFilter filter =
                    (FixedDistanceCutoffFilter)
                            manager.getFilterInstance(
                                    SourceFilters.FIXED_DIST_CUTOFF);

            filter.setMaxDistance(maxDistance);
        }

        return manager;
    }
}