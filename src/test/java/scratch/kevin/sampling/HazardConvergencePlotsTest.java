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
import scratch.kevin.sampling.HazardConvergencePlots.SummaryRow;

public class HazardConvergencePlotsTest {
	@Rule public TemporaryFolder output = new TemporaryFolder();

	@Test public void discoversCountsWithoutLHSAndWithSingleMethod() throws Exception {
		List<ReferenceSummary> rows = new ArrayList<>();
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			for (int count : new int[] {512, 1024})
				rows.add(new ReferenceSummary(SamplingMethod.OWEN_SCRAMBLED_SOBOL, count,
						HazardConvergenceCalcs.MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
						2, 0.2, 0.05, 0.1, 0.2, 0.3, 0.25, new double[] {0.15, 0.25}));
			rows.add(new ReferenceSummary(SamplingMethod.MONTE_CARLO, 512,
					HazardConvergenceCalcs.LOO_MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
					4, 0.3, 0.05, 0.2, 0.3, 0.4, 0.2, new double[] {0.2, 0.25, 0.35, 0.4}));
			rows.add(new ReferenceSummary(SamplingMethod.LATIN_HYPERCUBE, 512,
					HazardConvergenceCalcs.MCS_REFERENCE_NAME, metric, ConvergenceSummary.MEAN_ABSOLUTE,
					4, 0.2, 0.05, 0.1, 0.2, 0.3, 0.25, new double[] {0.1, 0.15, 0.25, 0.3}));
		}
		for (ReferenceSummary row : new ArrayList<>(rows))
			rows.add(new ReferenceSummary(row.method(), row.sampleCount(), row.reference(), row.metric(),
					ConvergenceSummary.MAXIMUM_ABSOLUTE, row.realizations(), 1d, 0.25, 0.5, 1d, 5d, 0.2,
					row.individualValues()));
		HazardConvergencePlots.plotMethodReference(output.getRoot(), rows, false,
				ConvergenceSummary.MEAN_ABSOLUTE, "Spatial mean absolute difference (%)");
		for (int count : new int[] {512, 1024}) {
			assertTrue(new java.io.File(output.getRoot(),
					"method_comparison_"+count+"_pooled_mcs_mean_abs.png").isFile());
		}
	}

	@Test public void combinesFullPoolAndLeaveOneOutRealizations() {
		ReferenceSummary fullPool = new ReferenceSummary(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 4096,
				HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME, ConvergenceMetric.MEAN_HAZARD,
				ConvergenceSummary.MEAN_ABSOLUTE, 4, 2.5, 0d, 1d, 2.5, 4d, 0d,
				new double[] {1d, 2d, 3d, 4d});
		ReferenceSummary leaveOneOut = new ReferenceSummary(SamplingMethod.OWEN_SCRAMBLED_SOBOL, 4096,
				HazardConvergenceCalcs.LOO_SOBOL_REFERENCE_NAME, ConvergenceMetric.MEAN_HAZARD,
				ConvergenceSummary.MEAN_ABSOLUTE, 5, 7d, 0d, 5d, 7d, 9d, 0d,
				new double[] {5d, 6d, 7d, 8d, 9d});

		SummaryRow combined = HazardConvergencePlots.combineRows(List.of(fullPool, leaveOneOut));
		assertArrayEquals(new double[] {1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d, 9d},
				combined.individualValues(), 0d);
		assertEquals(5d, combined.mean(), 0d);
		assertEquals(5d, combined.median(), 0d);
		assertEquals(1d, combined.minimum(), 0d);
		assertEquals(9d, combined.maximum(), 0d);
	}
}
