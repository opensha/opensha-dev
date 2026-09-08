package scratch.kevin.sampling;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceMetric;
import scratch.kevin.sampling.HazardConvergenceCalcs.ConvergenceSummary;
import scratch.kevin.sampling.HazardConvergencePlots.ReferenceSummary;

public class HazardConvergencePlotsTest {
	@Rule public TemporaryFolder output = new TemporaryFolder();

	@Test public void discoversCountsWithoutLHSAndWithSingleMethod() throws Exception {
		List<ReferenceSummary> rows = new ArrayList<>();
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			for (int count : new int[] {512, 1024})
				rows.add(new ReferenceSummary(SamplingMethod.OWEN_SCRAMBLED_SOBOL, count,
						HazardConvergenceCalcs.MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
						2, 0.1, 0.2, 0.3));
			rows.add(new ReferenceSummary(SamplingMethod.MONTE_CARLO, 512,
					HazardConvergenceCalcs.LOO_MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
					4, 0.2, 0.3, 0.4));
			rows.add(new ReferenceSummary(SamplingMethod.LATIN_HYPERCUBE, 512,
					HazardConvergenceCalcs.MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
					4, 0.1, 0.2, 0.3));
		}
		for (ReferenceSummary row : new ArrayList<>(rows))
			rows.add(new ReferenceSummary(row.method(), row.sampleCount(), row.reference(), row.metric(),
					ConvergenceSummary.MAXIMUM_ABSOLUTE, row.realizations(), 0.5, 1, 5));
		HazardConvergencePlots.plotMethodReference(output.getRoot(), rows, false,
				ConvergenceSummary.MEAN_ABSOLUTE, "Spatial mean absolute difference (%)");
		for (int count : new int[] {512, 1024}) {
			assertTrue(new java.io.File(output.getRoot(),
					"method_comparison_"+count+"_mcs_reference_mean_abs.png").isFile());
		}
	}
}
